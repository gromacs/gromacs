/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
#ifndef GMX_EWALD_PME_PERSISTENT_DATA_H
#define GMX_EWALD_PME_PERSISTENT_DATA_H

#include "config.h"

#include <memory>

struct gmx_device_info_t;

/*! \brief \internal
 * PME persistent host data, which should be initialized once per device context for the whole execution.
 * Using the PmePersistentDataHandle typedef below is recommended for retaining the data.
 *
 * Primary purpose will be to not recompile GPU kernels for each OpenCL unit test.
 * In CUDA, this just assigns the kernel function pointers.
 * This also relies on the fact that most of the kernels are used
 * (handles to 9 kernels are prepared, but only 4 are used during typical simulation).
 * If PME GPU supported different interpolation orders, even smaller share of all possible kernels would be used.
 * This also doesn't manage cuFFT/clFFT kernels, which depend on the PME grid dimensions.
 *
 * So far there is no CPU-only data here.
 */
struct PmeGpuPersistentData //TODO rename?
{
#if GMX_GPU != GMX_GPU_NONE

    //! Conveniently all the PME kernels use the same single argument type
#if GMX_GPU == GMX_GPU_CUDA
    using PmeKernelHandle = void(*)(const struct PmeGpuCudaKernelParams);
#elif GMX_GPU == GMX_GPU_OPENCL
    using PmeKernelHandle = cl_kernel;
#endif
    // TODO: All these kernels are compiled during pme_gpu_init() only for the given PME order!
    // (and only order of 4 is supported now, anyway).
    // spreading kernels also have hardcoded X/Y indices wrapping parameters as a placeholder for implementing
    // 1/2D decomposition.
    PmeKernelHandle splineKernel;
    PmeKernelHandle spreadKernel;
    PmeKernelHandle splineAndSpreadKernel;
    // Same for gather: hardcoded X/Y unwrap parameters, order of 4,
    // + it can reduce with previous forces in the host buffer, or ignore it.
    PmeKernelHandle gatherReduceWithInputKernel;
    PmeKernelHandle gatherKernel;
    // solve kernel doesn't care about spline order, but can optionally compute energy and virial,
    // and supports XYZ and YZX grid orderings.
    PmeKernelHandle solveYZXKernel;
    PmeKernelHandle solveXYZKernel;
    PmeKernelHandle solveYZXEnergyKernel;
    PmeKernelHandle solveXYZEnergyKernel;

#endif // GMX_GPU

    PmeGpuPersistentData() = delete;
    explicit PmeGpuPersistentData(struct gmx_pme_t *pme);
    ~PmeGpuPersistentData();
    PmeGpuPersistentData &operator=(const PmeGpuPersistentData &) = delete;
    PmeGpuPersistentData(const PmeGpuPersistentData &)            = delete;
};

//! The handle for manipulating persistent data. Storing the shared pointer in unit tests environment guarantees the persistency.
using PmePersistentDataHandle = std::shared_ptr<PmeGpuPersistentData>;

//! Factory function used by unit tests to build persistent data for the device at once
PmePersistentDataHandle buildPersistentPmeData(gmx_device_info_t *deviceInfo);

#endif
