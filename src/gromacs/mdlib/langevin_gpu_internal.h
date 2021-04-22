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
/*! \internal \file
 *
 * \brief Implements Langevin (SD) integrator using CUDA
 *
 * This file contains implementation of the Langevin (SD) integrator
 * using CUDA, including class initialization, data-structures management
 * and GPU kernel.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 *
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_LANGEVIN_GPU_INTERNAL_H
#define GMX_MDLIB_LANGEVIN_GPU_INTERNAL_H

#include "gmxpre.h"

#include "langevin_gpu.h"

namespace gmx
{
/*! \brief Backend-specific function to launch GPU Langevin (SD) kernel
 *
 *  The coordinates and velocities are updated on the GPU.
 *
 *  Each GPU thread works with a single particle. Empty declaration is needed to
 *  avoid "no previous prototype for function" clang warning.
 *
 * \param[in]     numAtoms            Total number of atoms.
 * \param[in,out] d_x                 Buffer containing initial coordinates, and where the updated ones will be written.
 * \param[out]    d_xp                Buffer where a copy of the initial coordinates will be written.
 * \param[in,out] d_v                 Buffer containing initial velocities, and where the updated ones will be written.
 * \param[in]     d_f                 Buffer containing forces.
 * \param[in]     d_inverseMasses     Buffer containing atoms' reciprocal masses.
 * \param[in]     dt                  Timestep.
 * \param[in]     seed                Random seed for sd integrator.
 * \param[in]     step                The step number in the simulation.
 * \param[in]     d_tempCouplGroups   Mapping of atoms into temperate coupling groups.
 * \param[in]     d_sdSigmaV          The sigma of the stochastic dynamics noise.
 * \param[in]     d_sdConstEm         EM (alpha in the SD equation). We lose precision (compared to the CPU version) by using float.
 * \param[in]     d_distributionTable The normal distribution table.
 * \param[in]     updateType          The update type of the SD integrator.
 * \param         deviceStream        Device stream for kernel launch.
 */
void launchLangevinKernel(int                                numAtoms,
                          DeviceBuffer<Float3>               d_x,
                          DeviceBuffer<Float3>               d_xp,
                          DeviceBuffer<Float3>               d_v,
                          const DeviceBuffer<Float3>         d_f,
                          const DeviceBuffer<float>          d_inverseMasses,
                          float                              dt,
                          int                                seed,
                          int                                step,
                          const DeviceBuffer<unsigned short> d_tempCouplGroups,
                          const DeviceBuffer<float>          d_sdSigmaV,
                          const DeviceBuffer<float>          d_sdConstEm,
                          const DeviceBuffer<float>          d_distributionTable,
                          SDUpdate                           updateType,
                          const DeviceStream&                deviceStream);

} // namespace gmx

#endif
