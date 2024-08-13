/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 *  \brief
 *  Common code for GPU buffer operations, namely the coordinate layout conversion
 *
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/gpu_utils/device_stream.h"
#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/nbnxm_gpu_buffer_ops_internal.h"
#if GMX_GPU_CUDA
#    include "gromacs/nbnxm/cuda/nbnxm_cuda_types.h"
#elif GMX_GPU_OPENCL
#    include "gromacs/nbnxm/opencl/nbnxm_ocl_types.h"
#elif GMX_GPU_SYCL
#    include "gromacs/nbnxm/sycl/nbnxm_sycl_types.h"
#endif
#include "gromacs/utility/exceptions.h"

namespace gmx
{

void nbnxn_gpu_x_to_nbat_x(const Grid&           grid,
                           NbnxmGpu*             nb,
                           DeviceBuffer<RVec>    d_x,
                           GpuEventSynchronizer* xReadyOnDevice,
                           const AtomLocality    locality,
                           int                   gridId,
                           int                   numColumnsMax,
                           bool                  mustInsertNonLocalDependency)
{
    GMX_ASSERT(bool(GMX_GPU_CUDA) || bool(GMX_GPU_SYCL),
               "NBNXM X buffer operations only supported in CUDA and SYCL");
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");
    InteractionLocality interactionLoc = atomToInteractionLocality(locality);

    const DeviceStream& deviceStream = *nb->deviceStreams[interactionLoc];

    const int numAtoms = grid.srcAtomEnd() - grid.srcAtomBegin();

    // Only insert wait on the first iteration of the loop.
    if (xReadyOnDevice != nullptr)
    {
        xReadyOnDevice->enqueueWaitEvent(deviceStream);
    }

    // avoid empty kernel launch, skip to inserting stream dependency
    if (numAtoms != 0)
    {
        GMX_ASSERT(d_x, "Need a valid device pointer");
        launchNbnxmKernelTransformXToXq(grid, nb, d_x, deviceStream, numColumnsMax, gridId);
    }

    if (mustInsertNonLocalDependency)
    {
        nbnxnInsertNonlocalGpuDependency(nb, interactionLoc);
    }
}

} // namespace gmx
