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
#ifndef GMX_GPU_UTILS_DEVICE_MACROS_H
#define GMX_GPU_UTILS_DEVICE_MACROS_H

#include "config.h"

/* The following lightweight attribute macros (GMX_HOST_ATTRIBUTE,
 * GMX_DEVICE_ATTRIBUTE, GMX_KERNEL_ATTRIBUTE, GMX_HOSTDEVICE_ATTRIBUTE,
 * GMX_FUNC_ATTRIBUTE, GMX_DEVICE_FUNC_ATTRIBUTE, GMX_ALWAYS_INLINE_ATTRIBUTE)
 * are defined in this header so they can be used without pulling in the
 * heavy GPU toolkit headers.
 */
#if (GMX_GPU_CUDA || GMX_GPU_HIP) && (defined(__CUDACC__) || defined(__HIPCC__))
#    define GMX_HOST_ATTRIBUTE __host__
#    define GMX_DEVICE_ATTRIBUTE __device__
#    define GMX_KERNEL_ATTRIBUTE __global__
#else
#    define GMX_HOST_ATTRIBUTE
#    define GMX_DEVICE_ATTRIBUTE
#    define GMX_KERNEL_ATTRIBUTE
#endif

#if defined(_MSC_VER) // MSVC does not support __attribute__ and always_inline
#    define GMX_ALWAYS_INLINE_ATTRIBUTE __forceinline
#else
#    define GMX_ALWAYS_INLINE_ATTRIBUTE __attribute__((always_inline))
#endif

#define GMX_HOSTDEVICE_ATTRIBUTE GMX_HOST_ATTRIBUTE GMX_DEVICE_ATTRIBUTE
#define GMX_FUNC_ATTRIBUTE GMX_HOSTDEVICE_ATTRIBUTE GMX_ALWAYS_INLINE_ATTRIBUTE
#define GMX_DEVICE_FUNC_ATTRIBUTE GMX_DEVICE_ATTRIBUTE GMX_ALWAYS_INLINE_ATTRIBUTE

/*!\brief Define a assert that can be used in both host and device code
 *
 * Use a plain device-side assert during the GPU device compilation pass
 * and regular (GROMACS) throw in host code.
 */
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
#    define GMX_HOST_DEVICE_THROW(e) assert(false)
#else
#    define GMX_HOST_DEVICE_THROW(e) GMX_THROW(InternalError(e))
#endif

#endif // GMX_GPU_UTILS_DEVICE_MACROS_H
