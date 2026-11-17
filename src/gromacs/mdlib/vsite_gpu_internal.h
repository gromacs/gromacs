/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Declares GPU kernel launch functions for virtual site construction and force spreading.
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_VSITE_GPU_INTERNAL_H
#define GMX_MDLIB_VSITE_GPU_INTERNAL_H

#include "gromacs/gpu_utils/devicebuffer_datatype.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/mdlib/vsite_gpu.h"
#include "gromacs/topology/ifunc.h"

class DeviceStream;
struct PbcAiuc;

namespace gmx
{

/*!
 * \brief Launch the virtual site construction kernel.
 *
 * Dispatches to the appropriate kernel based on the vsite type. Each kernel reads the positions
 * of the constructing atoms and writes the computed vsite position to the \p d_x buffer.
 * Does nothing if there are no vsites of type \p vsiteType.
 *
 * \tparam        vsiteType      Type of virtual site to construct.
 * \param[in]     kernelParams   Per-type vsite counts and PBC info.
 * \param[in]     kernelBuffers  Per-type device buffers.
 * \param[in,out] d_x            Coordinate buffer.
 * \param[in]     deviceStream   CUDA stream to launch kernel in.
 */
template<InteractionFunction vsiteType>
void launchConstructVSiteKernel(const VsiteGpuKernelParameters& kernelParams,
                                const VsiteGpuKernelBuffers&    kernelBuffers,
                                DeviceBuffer<Float3>            d_x,
                                const DeviceStream&             deviceStream);

/*! \brief Launch the force spreading kernel.
 *
 * Dispatches to the appropriate kernel based on the vsite type. Each kernel reads the force on each
 * vsite and adds the appropriately weighted force to its constructing atoms.
 * Does nothing if there are no vsites of type \p vsiteType.
 *
 * \tparam        vsiteType      Type of virtual site to spread forces for.
 * \param[in]     kernelParams   Per-type vsite counts.
 * \param[in]     kernelBuffers  Per-type device buffers.
 * \param[in,out] d_f            Force buffer; vsite forces read and added to parent forces.
 * \param[in]     deviceStream   CUDA stream to launch kernel in.
 */
template<InteractionFunction vsiteType>
void launchSpreadVSiteForceKernel(const VsiteGpuKernelParameters& kernelParams,
                                  const VsiteGpuKernelBuffers&    kernelBuffers,
                                  DeviceBuffer<Float3>            d_f,
                                  const DeviceStream&             deviceStream);

/*! \brief Launch the fused virtual site construction kernel.
 *
 * Constructs the vsites of all supported types in a single kernel launch, each thread
 * constructing one vsite. Since all types are constructed concurrently, this may only be
 * called when no vsite is a constructing atom of another vsite; otherwise the per-type
 * launchers have to be used in the order of \c vsiteTypesOnGpu.
 * Does nothing if there are no vsites at all.
 *
 * \param[in]     kernelParams   Per-type vsite counts, thread ranges and PBC info.
 * \param[in]     kernelBuffers  Per-type device buffers.
 * \param[in,out] d_x            Coordinate buffer.
 * \param[in]     deviceStream   CUDA stream to launch kernel in.
 */
void launchConstructVSitesFusedKernel(const VsiteGpuKernelParameters& kernelParams,
                                      const VsiteGpuKernelBuffers&    kernelBuffers,
                                      DeviceBuffer<Float3>            d_x,
                                      const DeviceStream&             deviceStream);

/*! \brief Launch the fused force spreading kernel.
 *
 * Spreads the forces of the vsites of all supported types in a single kernel launch, each
 * thread spreading the force of one vsite. Since all types are spread concurrently, this may
 * only be called when no vsite is a constructing atom of another vsite; otherwise the per-type
 * launchers have to be used in the reverse order of \c vsiteTypesOnGpu.
 * Does nothing if there are no vsites at all.
 *
 * \param[in]     kernelParams   Per-type vsite counts and thread ranges.
 * \param[in]     kernelBuffers  Per-type device buffers.
 * \param[in,out] d_f            Force buffer; vsite forces read and added to parent forces.
 * \param[in]     deviceStream   CUDA stream to launch kernel in.
 */
void launchSpreadVSiteForcesFusedKernel(const VsiteGpuKernelParameters& kernelParams,
                                        const VsiteGpuKernelBuffers&    kernelBuffers,
                                        DeviceBuffer<Float3>            d_f,
                                        const DeviceStream&             deviceStream);

} // namespace gmx

#endif // GMX_MDLIB_VSITE_GPU_INTERNAL_H
