/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 *
 * \brief This file contains declarations of high-level functions used
 * by mdrun to compute energies and forces for listed interactions.
 *
 * Clients of libgromacs that want to evaluate listed interactions
 * should call functions declared here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_GPU_H
#define GMX_LISTED_FORCES_LISTED_FORCES_GPU_H

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"

class DeviceContext;
class DeviceStream;

struct gmx_enerdata_t;
struct gmx_ffparams_t;
struct gmx_mtop_t;
struct t_inputrec;
struct gmx_wallcycle;

namespace gmx
{

struct NBAtomDataGpu;
template<typename>
class ArrayRef;
class StepWorkload;

/*! \brief The number on bonded function types supported on GPUs */
static constexpr int numFTypesOnGpu = 8;

/*! \brief List of all bonded function types supported on GPUs
 *
 * \note This list should be in sync with the actual GPU code.
 * \note Perturbed interactions are not supported on GPUs.
 * \note The function types in the list are ordered on increasing value.
 * \note Currently bonded are only supported with CUDA and SYCL, not with OpenCL.
 */
constexpr std::array<int, numFTypesOnGpu> fTypesOnGpu = { F_BONDS,  F_ANGLES, F_UREY_BRADLEY,
                                                          F_PDIHS,  F_RBDIHS, F_IDIHS,
                                                          F_PIDIHS, F_LJ14 };

/*! \brief Checks whether the GROMACS build allows to compute bonded interactions on a GPU.
 *
 * \param[out] error  If non-null, the diagnostic message when bondeds cannot run on a GPU.
 *
 * \returns true when this build can run bonded interactions on a GPU, false otherwise.
 *
 * \throws std::bad_alloc when out of memory.
 */
bool buildSupportsListedForcesGpu(std::string* error);

/*! \brief Checks whether the input system allows to compute bonded interactions on a GPU.
 *
 * \param[in]  ir     Input system.
 * \param[in]  mtop   Complete system topology to search for supported interactions.
 * \param[out] error  If non-null, the error message if the input is not supported on GPU.
 *
 * \returns true if PME can run on GPU with this input, false otherwise.
 */
bool inputSupportsListedForcesGpu(const t_inputrec& ir, const gmx_mtop_t& mtop, std::string* error);

class ListedForcesGpu
{
public:
    /*! \brief Construct the manager with constant data and the stream to use.
     *
     * \param[in] ffparams                   Force-field parameters.
     * \param[in] electrostaticsScaleFactor  Scaling factor for the electrostatic potential
     * \param[in] numEnergyGroupsForListedForces  The number of energy groups used for listed forces
     *                                       (Coulomb constant, multiplied by the Fudge factor).
     * \param[in] deviceContext              GPU device context.
     * \param[in] deviceStream               GPU device stream.
     * \param[in] wcycle                     The wallclock counter.
     *
     * \note Only assigning all energies to energy group pair 0,0 is supported.
     *       Passing numEnergyGroupsForListedForces>1 will lead to an assertion failure.
     */
    ListedForcesGpu(const gmx_ffparams_t& ffparams,
                    float                 electrostaticsScaleFactor,
                    int                   numEnergyGroupsForListedForces,
                    const DeviceContext&  deviceContext,
                    const DeviceStream&   deviceStream,
                    gmx_wallcycle*        wcycle);
    //! Destructor
    ~ListedForcesGpu();

    /*! \brief Update flag whether there are bonded interactions suitable for the GPU.
     *
     * Intended to be called early during search steps so domainWork flags can be populated.
     */
    void updateHaveInteractions(const InteractionDefinitions& idef);

    /*! \brief Update lists of interactions from idef suitable for the GPU,
     * using the data structures prepared for PP work.
     *
     * Intended to be called after each neighbour search
     * stage. Copies the bonded interactions assigned to the GPU
     * to device data structures, and updates device buffers that
     * may have been updated after search.
     *
     * \param[in]     nbnxnAtomOrder   Mapping between rvec and NBNXM formats.
     * \param[in]     idef             List of interactions to compute.
     * \param[in,out] nbnxmAtomDataGpu Nbnxm GPU atom data (XQ and force buffers).
     */
    void updateInteractionListsAndDeviceBuffers(ArrayRef<const int>           nbnxnAtomOrder,
                                                const InteractionDefinitions& idef,
                                                NBAtomDataGpu*                nbnxmAtomDataGpu);
    /*! \brief
     * Update PBC data.
     *
     * Converts PBC data from t_pbc into the PbcAiuc format and stores the latter.
     *
     * \param[in] pbcType The type of the periodic boundary.
     * \param[in] box     The periodic boundary box matrix.
     * \param[in] canMoleculeSpanPbc  Whether one molecule can have atoms in different PBC cells.
     */
    void setPbc(PbcType pbcType, const matrix box, bool canMoleculeSpanPbc);

    /*! \brief Returns whether there are bonded interactions
     * assigned to the GPU
     *
     * \returns If the list of interaction has elements.
     */
    bool haveInteractions() const;

    /*! \brief Launches bonded kernel on a GPU
     *
     * \param[in]  stepWork  Simulation step work to determine if energy/virial are to be computed on this step.
     */
    void launchKernel(const gmx::StepWorkload& stepWork);

    /*! \brief Sets the PBC and launches bonded kernel on a GPU
     *
     * \param[in] pbcType The type of the periodic boundary.
     * \param[in] box     The periodic boundary box matrix.
     * \param[in] canMoleculeSpanPbc  Whether one molecule can have atoms in different PBC cells.
     * \param[in] stepWork  Simulation step work to determine if energy/virial are to be computed on this step.
     */
    void setPbcAndlaunchKernel(PbcType                  pbcType,
                               const matrix             box,
                               bool                     canMoleculeSpanPbc,
                               const gmx::StepWorkload& stepWork);

    /*! \brief Launches the transfer of computed bonded energies.
     */
    void launchEnergyTransfer();

    /*! \brief Waits on the energy transfer, and accumulates bonded energies to \c enerd.
     *
     * \param[in,out] enerd The energy data object to add energy terms to.
     */
    void waitAccumulateEnergyTerms(gmx_enerdata_t* enerd);

    /*! \brief Clears the device side energy buffer
     */
    void clearEnergies();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif // GMX_LISTED_FORCES_LISTED_FORCES_GPU_H
