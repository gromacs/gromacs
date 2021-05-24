/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief Declares class for GPU implementation of SETTLE
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_SETTLE_GPU_H
#define GMX_MDLIB_SETTLE_GPU_H

#include "gmxpre.h"

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/settle.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc_aiuc.h"
#include "gromacs/topology/topology.h"

class InteractionDefinitions;

namespace gmx
{

//! \internal \brief Indices of atoms in a water molecule
struct WaterMolecule
{
    //! Oxygen atom
    int ow1;
    //! First hydrogen atom
    int hw2;
    //! Second hydrogen atom
    int hw3;
};


/*! \internal \brief Class with interfaces and data for GPU version of SETTLE. */
class SettleGpu
{

public:
    /*! \brief Create SETTLE object
     *
     *  Extracts masses for oxygen and hydrogen as well as the O-H and H-H target distances
     *  from the topology data (mtop), check their values for consistency and calls the
     *  following constructor.
     *
     * \param[in] mtop           Topology of the system to gen the masses for O and H atoms and
     *                           target O-H and H-H distances. These values are also checked for
     *                           consistency.
     * \param[in] deviceContext  Device context (dummy in CUDA).
     * \param[in] deviceStream   Device stream to use.
     */
    SettleGpu(const gmx_mtop_t& mtop, const DeviceContext& deviceContext, const DeviceStream& deviceStream);

    ~SettleGpu();

    /*! \brief Apply SETTLE.
     *
     * Applies SETTLE to coordinates and velocities, stored on GPU. Data at pointers d_xp and
     * d_v change in the GPU memory. The results are not automatically copied back to the CPU
     * memory. Method uses this class data structures which should be updated when needed using
     * update method.
     *
     * \param[in]     d_x               Coordinates before timestep (in GPU memory)
     * \param[in,out] d_xp              Coordinates after timestep (in GPU memory). The
     *                                  resulting constrained coordinates will be saved here.
     * \param[in]     updateVelocities  If the velocities should be updated.
     * \param[in,out] d_v               Velocities to update (in GPU memory, can be nullptr
     *                                  if not updated)
     * \param[in]     invdt             Reciprocal timestep (to scale Lagrange
     *                                  multipliers when velocities are updated)
     * \param[in]     computeVirial     If virial should be updated.
     * \param[in,out] virialScaled      Scaled virial tensor to be updated.
     * \param[in]     pbcAiuc           PBC data.
     */
    void apply(const DeviceBuffer<Float3>& d_x,
               DeviceBuffer<Float3>        d_xp,
               bool                        updateVelocities,
               DeviceBuffer<Float3>        d_v,
               real                        invdt,
               bool                        computeVirial,
               tensor                      virialScaled,
               const PbcAiuc&              pbcAiuc);

    /*! \brief
     * Update data-structures (e.g. after NB search step).
     *
     * Updates the constraints data and copies it to the GPU. Should be
     * called if the particles were sorted, redistributed between domains, etc.
     * Does not recycle the data preparation routines from the CPU version.
     * All three atoms from single water molecule should be handled by the same GPU.
     *
     * SETTLEs atom ID's is taken from idef.il[F_SETTLE].iatoms.
     *
     * \param[in] idef    System topology
     */
    void set(const InteractionDefinitions& idef);

private:
    //! GPU context object
    const DeviceContext& deviceContext_;
    //! GPU stream
    const DeviceStream& deviceStream_;

    //! Scaled virial tensor (9 reals, GPU)
    std::vector<float> h_virialScaled_;
    //! Scaled virial tensor (9 reals, GPU)
    DeviceBuffer<float> d_virialScaled_;

    //! Number of settles
    int numSettles_ = 0;

    //! Indexes of atoms (.i for oxygen, .j and.k for hydrogens, CPU)
    std::vector<WaterMolecule> h_atomIds_;
    //! Indexes of atoms (.i for oxygen, .j and.k for hydrogens, GPU)
    DeviceBuffer<WaterMolecule> d_atomIds_;
    //! Current size of the array of atom IDs
    int numAtomIds_ = -1;
    //! Allocated size for the array of atom IDs
    int numAtomIdsAlloc_ = -1;

    //! Settle parameters
    SettleParameters settleParameters_;
};

} // namespace gmx

#endif // GMX_MDLIB_SETTLE_GPU_H
