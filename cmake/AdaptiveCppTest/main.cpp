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

// Always fail if we are not compiling with AdaptiveCpp/hipSYCL
#if !defined(__HIPSYCL__) && !defined(__ADAPTIVECPP__)
#    error "Neither __HIPSYCL__ nor __ADAPTIVECPP__ macro not defined. Please check your AdaptiveCpp/hipSYCL installation."
#endif

/* Next, we optionally check three backends.
 * If CHECK_CUDA_TARGET is defined:
 *  - If we are not compiling for CUDA (because we did not specify CUDA devices among targets),
 *  this test compilation will fail.
 *  - If we are compiling for CUDA, the compilation will proceed.
 * Same for HIP, LevelZero (SPIR-V) and SSCP (generic).
 *
 * This allows us to compile this test file with different -DCHECK_x_TARGET flags to see which
 * backends we are compiling for and report it to CMake.
 * */

// AdaptiveCpp 23.10.0 only defines __HIPSYCL_ENABLE_x_TARGET__
#if defined(CHECK_CUDA_TARGET) && !defined(__HIPSYCL_ENABLE_CUDA_TARGET__)
#    error "CUDA target not enabled";
#endif

#if defined(CHECK_HIP_TARGET) && !defined(__HIPSYCL_ENABLE_HIP_TARGET__)
#    error "HIP target not enabled";
#endif

#if defined(CHECK_LEVELZERO_TARGET) && !defined(__HIPSYCL_ENABLE_SPIRV_TARGET__)
#    error "LevelZero (SPIR-V) target not enabled"
#endif

#if defined(CHECK_GENERIC_TARGET) && !defined(__HIPSYCL_ENABLE_LLVM_SSCP_TARGET__)
#    error "Generic (SSCP) target not enabled"
#endif


int main() {}
