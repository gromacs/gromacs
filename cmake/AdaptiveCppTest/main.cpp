/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#include <sycl/sycl.hpp>

// Always fail if we are not compiling with AdaptiveCpp/hipSYCL
#if !defined(__HIPSYCL__) && !defined(__ADAPTIVECPP__)
#    error "Neither __HIPSYCL__ nor __ADAPTIVECPP__ macro not defined. Please check your AdaptiveCpp/hipSYCL installation."
#endif

/* Next, we check four compile paths and issue a warning is each one is triggered.
 * These warnings can later be checked in the compiler logs */

int main()
{
    sycl::queue q;
    // Empty kernel to probe compiler definitions in multipass mode
    q.parallel_for<class K>(sycl::range<1>{ 16 },
                            [=](sycl::id<1> itemIdx) {
#if defined(__SYCL_DEVICE_ONLY__) && defined(__NVPTX__)
#    warning GMX_SYCL_TEST_HAVE_CUDA_TARGET
#endif
#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
#    warning GMX_SYCL_TEST_HAVE_HIP_TARGET
#    if __AMDGCN_WAVEFRONT_SIZE == 64
#        warning GMX_SYCL_TEST_HAVE_HIP_WAVE64_TARGET
#    elif __AMDGCN_WAVEFRONT_SIZE == 32
#        warning GMX_SYCL_TEST_HAVE_HIP_WAVE32_TARGET
#    endif
#endif
#if defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))
#    warning GMX_SYCL_TEST_HAVE_SPIRV_TARGET
#endif
                            });
    return 0;
}

// Test SSCP/generic/single-pass compiler by checking ACpp-specific macro
#if defined(__HIPSYCL_ENABLE_LLVM_SSCP_TARGET__) || defined(__ACPP_ENABLE_LLVM_SSCP_TARGET__)
#    warning GMX_SYCL_TEST_HAVE_GENERIC_TARGET
#endif
