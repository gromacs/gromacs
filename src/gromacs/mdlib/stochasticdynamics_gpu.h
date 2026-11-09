/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * \brief Declarations for GPU implementation of Langevin (SD) integrator.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_STOCHASTICDYNAMICS_GPU_H
#define GMX_MDLIB_STOCHASTICDYNAMICS_GPU_H

#include "config.h"

#if GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gputraits.cuh"
#elif GMX_GPU_HIP
#    include "gromacs/gpu_utils/gputraits_hip.h"
#elif GMX_GPU_SYCL
#    include "gromacs/gpu_utils/gputraits_sycl.h"
#endif

#include <memory>

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/mdlib/stochasticdynamics.h"

class DeviceContext;
class DeviceStream;

namespace gmx
{

//! The number of table bits used to generate the normal distribution table.
static constexpr int sc_normalDistributionTableBits = 14;

class StochasticDynamicsGpu
{

public:
    /*! \brief Constructor.
     *
     * \param[in] deviceContext       Device context (dummy in CUDA).
     * \param[in] deviceStream        Device stream to use.
     * \param[in] numTempCouplGroups  Number of temperature groups.
     * \param[in] delta_t             Timestep.
     * \param[in] ref_t               The reference temperature for each temperature group.
     * \param[in] tau_t               The time constant of the temperature coupling.
     */
    StochasticDynamicsGpu(const DeviceContext& deviceContext,
                          const DeviceStream&  deviceStream,
                          int                  numTempCouplGroups,
                          float                delta_t,
                          const float*         ref_t,
                          const float*         tau_t);
    ~StochasticDynamicsGpu();

    /*! \brief Integrate
     *
     * Integrates the equation of motion using Stochastic Dynamics algorithm.
     * Updates coordinates and velocities on the GPU. The current coordinates are saved for constraints.
     *
     * \param[in,out] d_x                      Coordinates to update
     * \param[out]    d_xp                     Place to save the values of initial coordinates to.
     * \param[in,out] d_v                      Velocities (will be updated).
     * \param[in]     d_f                      Forces.
     * \param[in]     dt                       Timestep.
     * \param[in]     seed                     Random seed.
     * \param[in]     step                     The step number in the simulation.
     * \param[in]     updateStage              Whether to do forces only or friction and noise only (combined is not used yet).
     */
    void integrate(DeviceBuffer<RVec> d_x,
                   DeviceBuffer<RVec> d_xp,
                   DeviceBuffer<RVec> d_v,
                   DeviceBuffer<RVec> d_f,
                   real               dt,
                   int64_t            seed,
                   int64_t            step,
                   SDUpdateStage      updateStage);

    /*! \brief Set the integrator
     *
     * Allocates memory for inverse masses, and, if needed for temperature scaling factor(s),
     * temperature coupling groups and global atom indices. Copies inverse masses, temperature
     * coupling groups and global atom indices to the GPU.
     *
     * The copies are asynchronous, so the device stream has to be synchronized between calls.
     *
     * \param[in] numAtoms          Number of atoms in the system.
     * \param[in] inverseMasses     Inverse masses of atoms.
     * \param[in] tempCouplGroups   Maps the atom index to temperature coupling group.
     * \param[in] globalAtomIndices Maps the local atom index to the global atom index, which
     *                              selects the random number stream of each atom. Should be
     *                              empty when the local and global atom orders are the same.
     */
    void set(int                            numAtoms,
             ArrayRef<const real>           inverseMasses,
             ArrayRef<const unsigned short> tempCouplGroups,
             ArrayRef<const int>            globalAtomIndices);

private:
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! GPU stream
    const DeviceStream& deviceStream_;
    //! Number of atoms
    int numAtoms_;

    //! 1/mass for all atoms (GPU)
    DeviceBuffer<float> d_inverseMasses_;
    //! Current size of the reciprocal masses array
    int numInverseMasses_ = -1;
    //! Maximum size of the reciprocal masses array
    int numInverseMassesAlloc_ = -1;

    //! Temperature group mapping
    DeviceBuffer<unsigned short> d_tempCouplGroups_;
    //! Current size of the temperature coupling groups array.
    int numTempCouplGroups_ = -1;
    //! Maximum size of the temperature coupling groups array.
    int numTempCouplGroupsAlloc_ = -1;

    //! Whether the local atom order differs from the global one, so d_globalAtomIndices_ is used.
    bool haveGlobalAtomIndices_ = false;
    //! Pinned host copy of the global atom indices, needed for the asynchronous copy to the GPU.
    HostVector<int> h_globalAtomIndices_;
    //! Global atom index for each local atom (GPU).
    DeviceBuffer<int> d_globalAtomIndices_ = nullptr;
    //! Current size of the global atom indices array.
    int numGlobalAtomIndices_ = -1;
    //! Maximum size of the global atom indices array.
    int numGlobalAtomIndicesAlloc_ = -1;

    //! The sigma of the stochastic dynamics noise.
    DeviceBuffer<float> d_sdSigmaV_;
    //! EM (alpha in the SD equation). We lose precision (compared to the CPU version) by using float.
    DeviceBuffer<float> d_sdConstEm_;
    //! A device copy of the normal distribution table.
    DeviceBuffer<float> d_distributionTable_;
};

} // namespace gmx

#endif
