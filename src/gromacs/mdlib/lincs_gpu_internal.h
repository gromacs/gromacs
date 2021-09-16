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
 * \brief Declare backend-specific LINCS GPU functions
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Alan Gray <alang@nvidia.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_LINCS_GPU_INTERNAL_H
#define GMX_MDLIB_LINCS_GPU_INTERNAL_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"

class DeviceStream;

namespace gmx
{

struct LincsGpuKernelParameters;

//! Number of threads in a GPU block
constexpr static int c_threadsPerBlock = 256;

/*! \brief Backend-specific function to launch LINCS kernel.
 *
 * \param kernelParams LINCS parameters.
 * \param d_x Initial coordinates before the integration.
 * \param d_xp Coordinates after the integration which will be updated.
 * \param updateVelocities Whether to also update velocities.
 * \param d_v Velocities to update (ignored if \p updateVelocities is \c false).
 * \param invdt Reciprocal of timestep.
 * \param computeVirial Whether to compute the virial.
 * \param deviceStream Device stream for kernel launch.
 */
void launchLincsGpuKernel(LincsGpuKernelParameters*   kernelParams,
                          const DeviceBuffer<Float3>& d_x,
                          DeviceBuffer<Float3>        d_xp,
                          bool                        updateVelocities,
                          DeviceBuffer<Float3>        d_v,
                          real                        invdt,
                          bool                        computeVirial,
                          const DeviceStream&         deviceStream);

} // namespace gmx

#endif // GMX_MDLIB_LINCS_GPU_INTERNAL_H
