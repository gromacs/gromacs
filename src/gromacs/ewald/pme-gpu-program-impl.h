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
/*! \internal \file
 * \brief
 * Declares PmeGpuProgramImpl, which stores PME GPU (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_PME_GPU_PROGRAM_IMPL_H
#define GMX_EWALD_PME_PME_GPU_PROGRAM_IMPL_H

#include "config.h"

#include "gromacs/utility/classhelpers.h"

#if GMX_GPU == GMX_GPU_CUDA
#include "gromacs/gpu_utils/gputraits.cuh"
#elif GMX_GPU == GMX_GPU_OPENCL
#include "gromacs/gpu_utils/gputraits_ocl.h"
#elif GMX_GPU == GMX_GPU_NONE
// TODO place in gputraits_stub.h
using Context = void *;
#endif

struct gmx_device_info_t;

/*! \internal
 * \brief
 * PME GPU persistent host program/kernel data, which should be initialized once for the whole execution.
 *
 * Primary purpose of this is to not recompile GPU kernels for each OpenCL unit test,
 * while the relevant GPU context (e.g. cl_context) instance persists.
 * In CUDA, this just assigns the kernel function pointers.
 * This also implicitly relies on the fact that reasonable share of the kernels are always used.
 * If there were more template parameters, even smaller share of all possible kernels would be used.
 *
 * \todo In future if we would need to react to either user input or
 * auto-tuning to compile different kernels, then we might wish to
 * revisit the number of kernels we pre-compile, and/or the management
 * of their lifetime.
 *
 * This also doesn't manage cuFFT/clFFT kernels, which depend on the PME grid dimensions.
 *
 * TODO: pass cl_context to the constructor and not create it inside.
 * See also Redmine #2522.
 */
struct PmeGpuProgramImpl
{
    /*! \brief
     * This is a handle to the GPU context, which is just a dummy in CUDA,
     * but is created/destroyed by this class in OpenCL.
     * TODO: Later we want to be able to own the context at a higher level and not here,
     * but this class would still need the non-owning context handle to build the kernels.
     */
    Context context;

    //! Conveniently all the PME kernels use the same single argument type
#if GMX_GPU == GMX_GPU_CUDA
    using PmeKernelHandle = void(*)(const struct PmeGpuCudaKernelParams);
#elif GMX_GPU == GMX_GPU_OPENCL
    using PmeKernelHandle = cl_kernel;
#else
    using PmeKernelHandle = void *;
#endif

    /*! \brief
     * Maximum synchronous GPU thread group execution width.
     * "Warp" is a CUDA term which we end up reusing in OpenCL kernels as well.
     * For CUDA, this is a static value that comes from gromacs/gpu_utils/cuda_arch_utils.cuh;
     * for OpenCL, we have to query it dynamically.
     */
    size_t warpSize;

    //@{
    /**
     * Spread/spline kernels are compiled only for order of 4.
     * Spreading kernels also have hardcoded X/Y indices wrapping parameters,
     * as a placeholder for implementing 1/2D decomposition.
     */
    size_t          spreadWorkGroupSize;

    PmeKernelHandle splineKernel;
    PmeKernelHandle spreadKernel;
    PmeKernelHandle splineAndSpreadKernel;
    //@}

    //@{
    /** Same for gather: hardcoded X/Y unwrap parameters, order of 4, plus
     * it can either reduce with previous forces in the host buffer, or ignore them.
     */
    size_t          gatherWorkGroupSize;

    PmeKernelHandle gatherReduceWithInputKernel;
    PmeKernelHandle gatherKernel;
    //@}

    //@{
    /** Solve kernel doesn't care about the interpolation order, but can optionally
     * compute energy and virial, and supports XYZ and YZX grid orderings.
     */
    size_t          solveMaxWorkGroupSize;

    PmeKernelHandle solveYZXKernel;
    PmeKernelHandle solveXYZKernel;
    PmeKernelHandle solveYZXEnergyKernel;
    PmeKernelHandle solveXYZEnergyKernel;
    //@}

    PmeGpuProgramImpl() = delete;
    //! Constructor for the given device
    explicit PmeGpuProgramImpl(const gmx_device_info_t *deviceInfo);
    ~PmeGpuProgramImpl();
    GMX_DISALLOW_COPY_AND_ASSIGN(PmeGpuProgramImpl);

    private:
        // Compiles kernels, if supported. Called by the constructor.
        void compileKernels(const gmx_device_info_t *deviceInfo);
};

#endif
