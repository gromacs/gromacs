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
 * Implements PmeGpuProgramImpl, which stores permanent PME GPU context-derived data,
 * such as (compiled) kernel handles.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#include "gmxpre.h"

#include "pme-gpu-program-impl.h"

#include "config.h"

#include "pme-gpu-constants.h"
#include "pme-gpu-internal.h" // for GridOrdering enum
#include "pme-gpu-types-host.h"
#include "pme-grid.h"

#if GMX_GPU == GMX_GPU_OPENCL

#include "gromacs/gpu_utils/ocl_compiler.h"
#include "gromacs/utility/stringutil.h"

#elif GMX_GPU == GMX_GPU_CUDA

#include "gromacs/gpu_utils/cuda_arch_utils.cuh" // only for warp_size

//@{
/**
 * PME CUDA kernels forward declarations. Kernels are documented in their respective files.
 * TODO: just include kernel files instead?
 */
template <
    const int order,
    const bool computeSplines,
    const bool spreadCharges,
    const bool wrapX,
    const bool wrapY
    >
void pme_spline_and_spread_kernel(const PmeGpuCudaKernelParams kernelParams);

template<
    GridOrdering gridOrdering,
    bool computeEnergyAndVirial
    >
void pme_solve_kernel(const PmeGpuCudaKernelParams kernelParams);

template <
    const int order,
    const bool overwriteForces,
    const bool wrapX,
    const bool wrapY
    >
void pme_gather_kernel(const PmeGpuCudaKernelParams kernelParams);
//@}
#endif

PmeGpuProgramImpl::PmeGpuProgramImpl(const gmx_device_info_t *deviceInfo)
{
#if GMX_GPU == GMX_GPU_CUDA

    GMX_UNUSED_VALUE(deviceInfo);

    warpSize = warp_size;

    // PME interpolation order
    constexpr int  pmeOrder = 4;
    GMX_UNUSED_VALUE(pmeOrder);
    // These hardcoded spread/gather parameters refer to not-implemented PME GPU 2D decomposition in X/Y
    constexpr bool wrapX = true;
    constexpr bool wrapY = true;
    GMX_UNUSED_VALUE(wrapX);
    GMX_UNUSED_VALUE(wrapY);
    splineAndSpreadKernel       = pme_spline_and_spread_kernel<pmeOrder, true, true, wrapX, wrapY>;
    splineKernel                = pme_spline_and_spread_kernel<pmeOrder, true, false, wrapX, wrapY>;
    spreadKernel                = pme_spline_and_spread_kernel<pmeOrder, false, true, wrapX, wrapY>;
    gatherKernel                = pme_gather_kernel<pmeOrder, true, wrapX, wrapY>;
    gatherReduceWithInputKernel = pme_gather_kernel<pmeOrder, false, wrapX, wrapY>;
    solveXYZKernel              = pme_solve_kernel<GridOrdering::XYZ, false>;
    solveXYZEnergyKernel        = pme_solve_kernel<GridOrdering::XYZ, true>;
    solveYZXKernel              = pme_solve_kernel<GridOrdering::YZX, false>;
    solveYZXEnergyKernel        = pme_solve_kernel<GridOrdering::YZX, true>;

#elif GMX_GPU == GMX_GPU_OPENCL

    // Context creation (which should happen outside of this class, actually)
    cl_platform_id        platformId = deviceInfo->ocl_gpu_id.ocl_platform_id;
    cl_device_id          deviceId   = deviceInfo->ocl_gpu_id.ocl_device_id;
    cl_context_properties contextProperties[3];
    contextProperties[0] = CL_CONTEXT_PLATFORM;
    contextProperties[1] = (cl_context_properties) platformId;
    contextProperties[2] = 0; /* Terminates the list of properties */

    cl_int  clError;
    context = clCreateContext(contextProperties, 1, &deviceId, nullptr, nullptr, &clError);
    if (clError != CL_SUCCESS)
    {
        const std::string errorString = gmx::formatString("Failed to create context for PME on GPU #%s:\n OpenCL error %d: %s",
                                                          deviceInfo->device_name, clError, ocl_get_error_string(clError).c_str());
        GMX_THROW(gmx::InternalError(errorString));
    }

    warpSize = gmx::ocl::getWarpSize(context, deviceId);

    compileKernels(deviceInfo);
#endif
}

PmeGpuProgramImpl::~PmeGpuProgramImpl()
{
 #if GMX_GPU == GMX_GPU_OPENCL
    // TODO: log releasing errors
    cl_int gmx_used_in_debug stat = 0;
    stat |= clReleaseKernel(splineAndSpreadKernel);
    stat |= clReleaseKernel(splineKernel);
    stat |= clReleaseKernel(spreadKernel);
    stat |= clReleaseKernel(gatherKernel);
    stat |= clReleaseKernel(gatherReduceWithInputKernel);
    stat |= clReleaseKernel(solveXYZKernel);
    stat |= clReleaseKernel(solveXYZEnergyKernel);
    stat |= clReleaseKernel(solveYZXKernel);
    stat |= clReleaseKernel(solveYZXEnergyKernel);
    stat |= clReleaseProgram(program);
    stat |= clReleaseContext(context);
    GMX_ASSERT(stat == CL_SUCCESS, gmx::formatString("Failed to release PME OpenCL resources %d: %s",
                                                     stat, ocl_get_error_string(stat).c_str()).c_str());
#endif
}

void PmeGpuProgramImpl::compileKernels(const gmx_device_info_t *deviceInfo)
{
#if GMX_GPU == GMX_GPU_OPENCL
    program  = nullptr;
    /* Need to catch std::bad_alloc here and during compilation string handling. */
    try
    {
        /* Here we pass macros and static const int variables defined in include
         * files outside as macros, to avoid including those files
         * in the JIT compilation that happens at runtime.
         */

        constexpr int order = 4;

//FIXME the spline indcies shoudl go into an inlien function like they did in CUDA kernels

//FIXME thsi is duplicated in pme-gpu-internal.cpp and should be a method
        const int  c_gatherMaxThreadsPerBlock = std::min(c_gatherMaxWarpsPerBlock * warpSize,
                                                         deviceInfo->maxWorkGroupSize);
        const int  c_spreadMaxThreadsPerBlock = std::min(c_spreadMaxWarpsPerBlock * warpSize,
                                                         deviceInfo->maxWorkGroupSize);
        const int  c_solveMaxThreadsPerBlock = std::min(c_solveMaxWarpsPerBlock * warpSize,
                                                        deviceInfo->maxWorkGroupSize);

        const std::string commonDefines = gmx::formatString(
                    "-Dwarp_size=%d "
                    "-Dorder=%d "
                    "-DPME_SPREADGATHER_ATOMS_PER_WARP=%d "
                    "-DPME_SPREADGATHER_THREADS_PER_ATOM=%d "
                    "-Dc_usePadding=%d "
                    "-Dc_skipNeutralAtoms=%d "
                    "-Dc_pmeMaxUnitcellShift=%f "       // forwarding from pme-grid.h, used for spline computation table sizes only
                    "-DDIM=%d -DXX=%d -DYY=%d -DZZ=%d " // forwarding from vectypes.h
                    "-DwrapX=true -DwrapY=true "        // decomposition parameter placeholders
                    "-Dc_virialAndEnergyCount=%d "
                    "-Dc_spreadMaxThreadsPerBlock=%d -Dc_solveMaxThreadsPerBlock=%d -Dc_gatherMaxThreadsPerBlock=%d ",
                    warpSize,
                    order,
                    warpSize / PME_SPREADGATHER_THREADS_PER_ATOM,
                    PME_SPREADGATHER_THREADS_PER_ATOM,
                    c_usePadding,
                    c_skipNeutralAtoms,
                    static_cast<float>(c_pmeMaxUnitcellShift),
                    DIM, XX, YY, ZZ,
                    c_virialAndEnergyCount,
                    c_spreadMaxThreadsPerBlock, c_solveMaxThreadsPerBlock, c_gatherMaxThreadsPerBlock);
        try
        {
            /* TODO when we have a proper MPI-aware logging module,
               the log output here should be written there */
            program = gmx::ocl::compileProgram(stderr,
                                               "src/gromacs/ewald",
                                               "pme-program.cl",
                                               commonDefines,
                                               context,
                                               deviceInfo->ocl_gpu_id.ocl_device_id,
                                               deviceInfo->vendor_e);
        }
        catch (gmx::GromacsException &e)
        {
            e.prependContext(gmx::formatString("Failed to compile PME kernels for GPU #%s\n",
                                               deviceInfo->device_name));
            throw;
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    constexpr cl_uint      expectedKernelCount = 9;
    // Has to be equal or larger than the number of kernel instances.
    // If it is not, CL_INVALID_VALUE will be thrown.
    std::vector<cl_kernel> kernels(expectedKernelCount, nullptr);
    cl_uint                actualKernelCount = 0;
    cl_int                 clError           = clCreateKernelsInProgram(program, kernels.size(), kernels.data(), &actualKernelCount);
    if (clError != CL_SUCCESS)
    {
        const std::string errorString = gmx::formatString("Failed to create kernels for PME on GPU #%s:\n OpenCL error %d: %s",
                                                          deviceInfo->device_name, clError, ocl_get_error_string(clError).c_str());
        GMX_THROW(gmx::InternalError(errorString));
    }
    kernels.resize(actualKernelCount);

    std::array<char, 100> kernelNamesBuffer;
    for (const auto &kernel : kernels)
    {
        clError = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME,
                                  kernelNamesBuffer.size(), kernelNamesBuffer.data(), nullptr);
        if (clError != CL_SUCCESS)
        {
            const std::string errorString = gmx::formatString("Failed to parse kernels for PME on GPU #%s:\n OpenCL error %d: %s",
                                                              deviceInfo->device_name, clError, ocl_get_error_string(clError).c_str());
            GMX_THROW(gmx::InternalError(errorString));
        }

        // The names below must correspond to those defined in pme-program.cl
        // TODO use a map with string key instead?
        if (!strcmp(kernelNamesBuffer.data(), "pmeSplineKernel"))
        {
            splineKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSplineAndSpreadKernel"))
        {
            splineAndSpreadKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSpreadKernel"))
        {
            spreadKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeGatherKernel"))
        {
            gatherKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeGatherReduceWithInputKernel"))
        {
            gatherReduceWithInputKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSolveYZXKernel"))
        {
            solveYZXKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSolveYZXEnergyKernel"))
        {
            solveYZXEnergyKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSolveXYZKernel"))
        {
            solveXYZKernel = kernel;
        }
        else if (!strcmp(kernelNamesBuffer.data(), "pmeSolveXYZEnergyKernel"))
        {
            solveXYZEnergyKernel = kernel;
        }
    }
#else
    GMX_UNUSED_VALUE(deviceInfo);
#endif
}
