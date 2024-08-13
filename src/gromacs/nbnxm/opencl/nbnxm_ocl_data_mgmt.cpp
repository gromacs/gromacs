/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 *  \brief Define OpenCL implementation of nbnxm_gpu_data_mgmt.h
 *
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_nbnxm
 */
#include "gmxpre.h"

#include <cassert>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/hardware/device_information.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/gpu_jit_support.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_gpu.h"
#include "gromacs/nbnxm/nbnxm_gpu_data_mgmt.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "nbnxm_ocl_types.h"

namespace gmx
{

/*! \brief This parameter should be determined heuristically from the
 * kernel execution times
 *
 * This value is best for small systems on a single AMD Radeon R9 290X
 * (and about 5% faster than 40, which is the default for CUDA
 * devices). Larger simulation systems were quite insensitive to the
 * value of this parameter.
 */
static unsigned int gpu_min_ci_balanced_factor = 50;

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static cl_kernel nbnxn_gpu_create_kernel(NbnxmGpu* nb, const char* kernel_name)
{
    cl_kernel kernel;
    cl_int    cl_error;

    kernel = clCreateKernel(nb->dev_rundata->program, kernel_name, &cl_error);
    if (CL_SUCCESS != cl_error)
    {
        gmx_fatal(FARGS,
                  "Failed to create kernel '%s' for GPU #%s: OpenCL error %d",
                  kernel_name,
                  nb->deviceContext_->deviceInfo().device_name,
                  cl_error);
    }

    return kernel;
}

/*! \brief Initializes the OpenCL kernel pointers of the nbnxn_ocl_ptr_t input data structure. */
static void nbnxn_gpu_init_kernels(NbnxmGpu* nb)
{
    /* Init to 0 main kernel arrays */
    /* They will be later on initialized in select_nbnxn_kernel */
    // TODO: consider always creating all variants of the kernels here so that there is no
    // need for late call to clCreateKernel -- if that gives any advantage?
    std::memset(nb->kernel_ener_noprune_ptr, 0, sizeof(nb->kernel_ener_noprune_ptr));
    std::memset(nb->kernel_ener_prune_ptr, 0, sizeof(nb->kernel_ener_prune_ptr));
    std::memset(nb->kernel_noener_noprune_ptr, 0, sizeof(nb->kernel_noener_noprune_ptr));
    std::memset(nb->kernel_noener_prune_ptr, 0, sizeof(nb->kernel_noener_prune_ptr));

    /* Init pruning kernels
     *
     * TODO: we could avoid creating kernels if dynamic pruning is turned off,
     * but ATM that depends on force flags not passed into the initialization.
     */
    nb->kernel_pruneonly[epruneFirst] = nbnxn_gpu_create_kernel(nb, "nbnxn_kernel_prune_opencl");
    nb->kernel_pruneonly[epruneRolling] =
            nbnxn_gpu_create_kernel(nb, "nbnxn_kernel_prune_rolling_opencl");
}

void gpu_init_platform_specific(NbnxmGpu* nb)
{
    /* set device info, just point it to the right GPU among the detected ones */
    nb->dev_rundata = new gmx_device_runtime_data_t();

    /* Enable LJ param manual prefetch for AMD or Intel or if we request through env. var.
     * TODO: decide about NVIDIA
     */
    nb->bPrefetchLjParam = (std::getenv("GMX_OCL_DISABLE_I_PREFETCH") == nullptr)
                           && ((nb->deviceContext_->deviceInfo().deviceVendor == DeviceVendor::Amd)
                               || (nb->deviceContext_->deviceInfo().deviceVendor == DeviceVendor::Intel)
                               || (std::getenv("GMX_OCL_ENABLE_I_PREFETCH") != nullptr));

    /* NOTE: in CUDA we pick L1 cache configuration for the nbnxn kernels here,
     * but sadly this is not supported in OpenCL (yet?). Consider adding it if
     * it becomes supported.
     */
    nbnxn_gpu_compile_kernels(nb);
    nbnxn_gpu_init_kernels(nb);
}

/*! \brief Releases an OpenCL kernel pointer */
static void free_kernel(cl_kernel* kernel_ptr)
{
    GMX_ASSERT(kernel_ptr, "Need a valid kernel pointer");

    if (*kernel_ptr)
    {
        cl_int cl_error = clReleaseKernel(*kernel_ptr);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clReleaseKernel failed: " + ocl_get_error_string(cl_error)).c_str());

        *kernel_ptr = nullptr;
    }
}

/*! \brief Releases a list of OpenCL kernel pointers */
static void free_kernels(cl_kernel* kernels, int count)
{
    for (int i = 0; i < count; i++)
    {
        free_kernel(kernels + i);
    }
}

/*! \brief Free the OpenCL program.
 *
 *  The function releases the OpenCL program assuciated with the
 *  device that the calling PP rank is running on.
 *
 *  \param program [in]  OpenCL program to release.
 */
static void freeGpuProgram(cl_program program)
{
    if (program)
    {
        cl_int cl_error = clReleaseProgram(program);
        GMX_RELEASE_ASSERT(cl_error == CL_SUCCESS,
                           ("clReleaseProgram failed: " + ocl_get_error_string(cl_error)).c_str());
        program = nullptr;
    }
}

//! This function is documented in the header file
void gpu_free_platform_specific(NbnxmGpu* nb)
{
    /* Free kernels */
    // NOLINTNEXTLINE(bugprone-sizeof-expression)
    int kernel_count = sizeof(nb->kernel_ener_noprune_ptr) / sizeof(nb->kernel_ener_noprune_ptr[0][0]);
    free_kernels(nb->kernel_ener_noprune_ptr[0], kernel_count);

    // NOLINTNEXTLINE(bugprone-sizeof-expression)
    kernel_count = sizeof(nb->kernel_ener_prune_ptr) / sizeof(nb->kernel_ener_prune_ptr[0][0]);
    free_kernels(nb->kernel_ener_prune_ptr[0], kernel_count);

    // NOLINTNEXTLINE(bugprone-sizeof-expression)
    kernel_count = sizeof(nb->kernel_noener_noprune_ptr) / sizeof(nb->kernel_noener_noprune_ptr[0][0]);
    free_kernels(nb->kernel_noener_noprune_ptr[0], kernel_count);

    // NOLINTNEXTLINE(bugprone-sizeof-expression)
    kernel_count = sizeof(nb->kernel_noener_prune_ptr) / sizeof(nb->kernel_noener_prune_ptr[0][0]);
    free_kernels(nb->kernel_noener_prune_ptr[0], kernel_count);

    freeGpuProgram(nb->dev_rundata->program);
    delete nb->dev_rundata;
}

//! This function is documented in the header file
int gpu_min_ci_balanced(NbnxmGpu* nb)
{
    if (nb != nullptr)
    {
        return gpu_min_ci_balanced_factor * nb->deviceContext_->deviceInfo().compute_units
               / getDeviceComputeUnitFactor(nb->deviceContext_->deviceInfo());
    }
    else
    {
        return 0;
    }
}

} // namespace gmx
