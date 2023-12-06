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
/*! \libinternal \file
 *
 * \brief Declarations for backend specific GPU functions for Leap-Frog.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_LEAPFROG_GPU_INTERNAL_H
#define GMX_MDLIB_LEAPFROG_GPU_INTERNAL_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/mdlib/leapfrog_gpu.h"

namespace gmx
{

/*! \brief Backend-specific function to launch GPU Leap Frog kernel.
 *
 * \param numAtoms Total number of atoms.
 * \param[in,out] d_x Buffer containing initial coordinates, and where the updated ones will be written.
 * \param[out] d_x0 Buffer where a copy of the initial coordinates will be written.
 * \param[in,out] d_v Buffer containing initial velocities, and where the updated ones will be written.
 * \param[in]  d_f Buffer containing forces.
 * \param[in] d_inverseMasses Buffer containing atoms' reciprocal masses.
 * \param dt Timestep.
 * \param doTemperatureScaling Whether temperature scaling is needed.
 * \param numTempScaleValues Number of different T-couple values.
 * \param d_tempScaleGroups Mapping of atoms into temperature scaling groups.
 * \param d_lambdas Temperature scaling factors (one per group).
 * \param parrinelloRahmanVelocityScaling The properties of the Parrinello-Rahman velocity scaling matrix.
 * \param prVelocityScalingMatrixDiagonal Diagonal elements of Parrinello-Rahman velocity scaling matrix.
 * \param deviceStream Device stream for kernel launch.
 */
void launchLeapFrogKernel(int                             numAtoms,
                          DeviceBuffer<Float3>            d_x,
                          DeviceBuffer<Float3>            d_x0,
                          DeviceBuffer<Float3>            d_v,
                          DeviceBuffer<Float3>            d_f,
                          DeviceBuffer<float>             d_inverseMasses,
                          float                           dt,
                          bool                            doTemperatureScaling,
                          int                             numTempScaleValues,
                          DeviceBuffer<unsigned short>    d_tempScaleGroups,
                          DeviceBuffer<float>             d_lambdas,
                          ParrinelloRahmanVelocityScaling parrinelloRahmanVelocityScaling,
                          Float3                          prVelocityScalingMatrixDiagonal,
                          const DeviceStream&             deviceStream);


} // namespace gmx

#endif
