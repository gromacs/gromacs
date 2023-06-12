/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 *
 * \brief Declares backend-specific functions for GPU implementation of SETTLE.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_SETTLE_GPU_INTERNAL_H
#define GMX_MDLIB_SETTLE_GPU_INTERNAL_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/settle_gpu.h"

namespace gmx
{

/*! \brief Apply SETTLE.
 *
 * Applies SETTLE to coordinates and velocities, stored on GPU. Data at pointers d_xp and
 * d_v change in the GPU memory. The results are not automatically copied back to the CPU
 * memory. Method uses this class data structures which should be updated when needed using
 * update method.
 *
 * \param[in]     numSettles        Number of SETTLE constraints.
 * \param[in]     d_atomIds         Device buffer with indices of atoms to be SETTLEd.
 * \param[in]     settleParameters  Parameters for SETTLE constraints.
 * \param[in]     d_x               Coordinates before timestep (in GPU memory)
 * \param[in,out] d_xp              Coordinates after timestep (in GPU memory). The
 *                                  resulting constrained coordinates will be saved here.
 * \param[in]     updateVelocities  If the velocities should be updated.
 * \param[in,out] d_v               Velocities to update (in GPU memory, can be nullptr
 *                                  if not updated)
 * \param[in]     invdt             Reciprocal timestep (to scale Lagrange
 *                                  multipliers when velocities are updated)
 * \param[in]     computeVirial     If virial should be updated.
 * \param[in,out] d_virialScaled      Scaled virial tensor to be updated.
 * \param[in]     pbcAiuc           PBC data.
 * \param[in]     deviceStream      Device stream to launch kernel in.
 */
void launchSettleGpuKernel(int                                numSettles,
                           const DeviceBuffer<WaterMolecule>& d_atomIds,
                           const SettleParameters&            settleParameters,
                           const DeviceBuffer<Float3>&        d_x,
                           DeviceBuffer<Float3>               d_xp,
                           bool                               updateVelocities,
                           DeviceBuffer<Float3>               d_v,
                           real                               invdt,
                           bool                               computeVirial,
                           DeviceBuffer<float>                d_virialScaled,
                           const PbcAiuc&                     pbcAiuc,
                           const DeviceStream&                deviceStream);

} // namespace gmx

#endif // GMX_MDLIB_SETTLE_GPU_INTERNAL_H
