/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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

// Always fail if we are not compiling with hipSYCL
#if !defined(__HIPSYCL__)
#    error "__HIPSYCL__ macro not defined. Please check your hipSYCL installation."
#endif

/* Next, we optionally check three backends.
 * If CHECK_CUDA_TARGET is defined:
 *  - If we are not compiling for CUDA (because we did not specify CUDA devices among targets),
 *  this test compilation will fail.
 *  - If we are compiling for CUDA, the compilation will proceed.
 * Same for HIP and LevelZero.
 *
 * This allows us to compile this test file with different -DCHECK_x_TARGET flags to see which
 * backends we are compiling for and report it to CMake.
 * */

#if defined(CHECK_CUDA_TARGET) && !defined(__HIPSYCL_ENABLE_CUDA_TARGET__)
#    error "CUDA target not enabled";
#endif

#if defined(CHECK_HIP_TARGET) && !defined(__HIPSYCL_ENABLE_HIP_TARGET__)
#    error "HIP target not enabled";
#endif

#if defined(CHECK_LEVELZERO_TARGET) && !defined(__HIPSYCL_ENABLE_SPIRV_TARGET__)
#    error "LevelZero (SPIR-V) target not enabled"
#endif

int main() {}
