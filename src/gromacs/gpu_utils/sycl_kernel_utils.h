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

#ifndef GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H
#define GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H

#include "gmxsycl.h"

/*! \file
 *  \brief SYCL kernel helper functions.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

/*! \brief Access mode to use for atomic accessors.
 *
 * Intel DPCPP compiler has \c sycl::atomic_ref, but has no \c sycl::atomic_fetch_add for floats.
 * However, \c atomic_ref can not be constructed from \c sycl::atomic, so we can not use
 * atomic accessors. Thus, we use \c mode::read_write accessors and \c atomic_ref.
 *
 * hipSYCL does not have \c sycl::atomic_ref, but has \c sycl::atomic_fetch_add for floats, which
 * requires using atomic accessors. Thus, we use \c mode::atomic accessors.
 *
 * The \ref atomicFetchAdd function could be used for doing operations on such accessors.
 */
static constexpr auto mode_atomic = GMX_SYCL_DPCPP ? cl::sycl::access::mode::read_write :
                                                   /* GMX_SYCL_HIPSYCL */ cl::sycl::access::mode::atomic;

//! \brief Full warp active thread mask used in CUDA warp-level primitives.
static constexpr unsigned int c_cudaFullWarpMask = 0xffffffff;

/*! \brief Convenience wrapper to do atomic addition to a global buffer.
 *
 * The implementation differences between DPCPP and hipSYCL are explained in \ref mode_atomic.
 */
template<class IndexType>
static inline void atomicFetchAdd(DeviceAccessor<float, mode_atomic> acc, const IndexType idx, const float val)
{
#if GMX_SYCL_DPCPP
    sycl_2020::atomic_ref<float, sycl_2020::memory_order::relaxed, sycl_2020::memory_scope::device, cl::sycl::access::address_space::global_space>
            fout_atomic(acc[idx]);
    fout_atomic.fetch_add(val);
#elif GMX_SYCL_HIPSYCL
#    ifdef SYCL_DEVICE_ONLY
    /* While there is support for float atomics on device, the host implementation uses
     * Clang's __atomic_fetch_add intrinsic, that, at least in Clang 11, does not support
     * floats. Luckily, we don't want to run on host. */
    // The pragmas below can be removed once we switch to sycl::atomic
#        pragma clang diagnostic push
#        pragma clang diagnostic ignored "-Wdeprecated-declarations"
    acc[idx].fetch_add(val);
#        pragma clang diagnostic push
#    else
    GMX_ASSERT(false, "hipSYCL host codepath not supported");
    GMX_UNUSED_VALUE(val);
    GMX_UNUSED_VALUE(acc);
    GMX_UNUSED_VALUE(idx);
#    endif
#endif
}

/*! \brief Issue an intra sub-group barrier.
 *
 * Equivalent with CUDA's \c syncwarp(c_cudaFullWarpMask).
 *
 */
static inline void subGroupBarrier(const cl::sycl::nd_item<1> itemIdx)
{
#if GMX_SYCL_HIPSYCL
    cl::sycl::group_barrier(itemIdx.get_sub_group(), cl::sycl::memory_scope::sub_group);
#else
    itemIdx.get_sub_group().barrier();
#endif
}

namespace sycl_2020
{
#if GMX_SYCL_HIPSYCL
__device__ __host__ static inline float shift_left(sycl_2020::sub_group,
                                                   float                                var,
                                                   sycl_2020::sub_group::linear_id_type delta)
{
    // No sycl::sub_group::shift_left / shuffle_down in hipSYCL yet
#    ifdef SYCL_DEVICE_ONLY
#        if defined(HIPSYCL_PLATFORM_CUDA) && defined(__HIPSYCL_ENABLE_CUDA_TARGET__)
    return __shfl_down_sync(c_cudaFullWarpMask, var, delta);
#        elif defined(HIPSYCL_PLATFORM_ROCM) && defined(__HIPSYCL_ENABLE_HIP_TARGET__)
    // Do we need more ifdefs? https://github.com/ROCm-Developer-Tools/HIP/issues/1491
    return __shfl_down(var, delta);
#        else
#            error "Unsupported hipSYCL target"
#        endif
#    else
    // Should never be called
    GMX_UNUSED_VALUE(var);
    GMX_UNUSED_VALUE(delta);
    assert(false);
    return NAN;
#    endif
}
#elif GMX_SYCL_DPCPP
static inline float shift_left(sycl_2020::sub_group sg, float var, sycl_2020::sub_group::linear_id_type delta)
{
    return sg.shuffle_down(var, delta);
}
#endif

#if GMX_SYCL_HIPSYCL
__device__ __host__ static inline float shift_right(sycl_2020::sub_group,
                                                    float                                var,
                                                    sycl_2020::sub_group::linear_id_type delta)
{
    // No sycl::sub_group::shift_right / shuffle_up in hipSYCL yet
#    ifdef SYCL_DEVICE_ONLY
#        if defined(HIPSYCL_PLATFORM_CUDA) && defined(__HIPSYCL_ENABLE_CUDA_TARGET__)
    return __shfl_up_sync(c_cudaFullWarpMask, var, delta);
#        elif defined(HIPSYCL_PLATFORM_ROCM) && defined(__HIPSYCL_ENABLE_HIP_TARGET__)
    // Do we need more ifdefs? https://github.com/ROCm-Developer-Tools/HIP/issues/1491
    return __shfl_up(var, delta);
#        else
#            error "Unsupported hipSYCL target"
#        endif
#    else
    // Should never be called
    assert(false);
    GMX_UNUSED_VALUE(var);
    GMX_UNUSED_VALUE(delta);
    return NAN;
#    endif
}
#elif GMX_SYCL_DPCPP
static inline float shift_right(sycl_2020::sub_group sg, float var, sycl_2020::sub_group::linear_id_type delta)
{
    return sg.shuffle_up(var, delta);
}
#endif
} // namespace sycl_2020

#endif /* GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H */
