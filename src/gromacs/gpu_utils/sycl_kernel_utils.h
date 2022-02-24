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

#ifndef GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H
#define GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H

#include "gmxsycl.h"

/*! \file
 *  \brief SYCL kernel helper functions.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

//! \brief Full warp active thread mask used in CUDA warp-level primitives.
static constexpr unsigned int c_cudaFullWarpMask = 0xffffffff;

/*! \brief Convenience wrapper to do atomic addition to a global buffer.
 */
template<typename T, sycl_2020::memory_scope MemoryScope = sycl_2020::memory_scope::device>
static inline void atomicFetchAdd(T& val, const T delta)
{
    sycl_2020::atomic_ref<T, sycl_2020::memory_order::relaxed, MemoryScope, sycl::access::address_space::global_space> ref(
            val);
    ref.fetch_add(delta);
}

/*! \brief Convenience wrapper to do atomic loads from a global buffer.
 */
template<typename T, sycl_2020::memory_scope MemoryScope = sycl_2020::memory_scope::device>
static inline T atomicLoad(T& val)
{
    sycl_2020::atomic_ref<T, sycl_2020::memory_order::relaxed, MemoryScope, sycl::access::address_space::global_space> ref(
            val);
    return ref.load();
}


/*! \brief Issue an intra sub-group barrier.
 *
 * Equivalent with CUDA's \c syncwarp(c_cudaFullWarpMask).
 *
 */
template<int Dim>
static inline void subGroupBarrier(const sycl::nd_item<Dim> itemIdx)
{
#if GMX_SYCL_HIPSYCL
    sycl::group_barrier(itemIdx.get_sub_group(), sycl::memory_scope::sub_group);
#else
    itemIdx.get_sub_group().barrier();
#endif
}

namespace sycl_2020
{
#if GMX_SYCL_HIPSYCL
__device__ __host__ static inline float shift_left(sycl::sub_group, float var, sycl::sub_group::linear_id_type delta)
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
static inline float shift_left(sycl::sub_group sg, float var, sycl::sub_group::linear_id_type delta)
{
    return sg.shuffle_down(var, delta);
}
#endif

#if GMX_SYCL_HIPSYCL
__device__ __host__ static inline float shift_right(sycl::sub_group, float var, sycl::sub_group::linear_id_type delta)
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
static inline float shift_right(sycl::sub_group sg, float var, sycl::sub_group::linear_id_type delta)
{
    return sg.shuffle_up(var, delta);
}
#endif

#if GMX_SYCL_HIPSYCL
/*! \brief Polyfill for sycl::isfinite missing from hipSYCL
 *
 * Does not follow GROMACS style because it should follow the name for
 * which it is a polyfill. */
template<typename Real>
__device__ __host__ static inline bool isfinite(Real value)
{
    // This is not yet implemented in hipSYCL pending
    // https://github.com/illuhad/hipSYCL/issues/636
#    ifdef SYCL_DEVICE_ONLY
#        if defined(HIPSYCL_PLATFORM_CUDA) && defined(__HIPSYCL_ENABLE_CUDA_TARGET__)
    return ::isfinite(value);
#        elif defined(HIPSYCL_PLATFORM_ROCM) && defined(__HIPSYCL_ENABLE_HIP_TARGET__)
    return ::isfinite(value);
#        else
#            error "Unsupported hipSYCL target"
#        endif
#    else
    // Should never be called
    assert(false);
    GMX_UNUSED_VALUE(value);
    return false;
#    endif
}
#elif GMX_SYCL_DPCPP
template<typename Real>
static inline bool isfinite(Real value)
{
    return sycl::isfinite(value);
}

#endif

#if GMX_SYCL_HIPSYCL

/*! \brief Polyfill for sycl::vec::load buggy in hipSYCL
 *
 * Loads from the address \c ptr offset in elements of type T by
 * NumElements * offset, into the components of \c v.
 *
 * Can probably be removed when
 * https://github.com/illuhad/hipSYCL/issues/647 is resolved. */
template<sycl::access::address_space AddressSpace, typename T, int NumElements>
static inline void loadToVec(size_t                                 offset,
                             sycl::multi_ptr<const T, AddressSpace> ptr,
                             sycl::vec<T, NumElements>*             v)
{
    for (int i = 0; i < NumElements; ++i)
    {
        (*v)[i] = ptr.get()[offset * NumElements + i];
    }
}

/*! \brief Polyfill for sycl::vec::store buggy in hipSYCL
 *
 * Loads from the address \c ptr offset in elements of type T by
 * NumElements * offset, into the components of \c v.
 *
 * Can probably be removed when
 * https://github.com/illuhad/hipSYCL/issues/647 is resolved. */
template<sycl::access::address_space AddressSpace, typename T, int NumElements>
static inline void storeFromVec(const sycl::vec<T, NumElements>& v,
                                size_t                           offset,
                                sycl::multi_ptr<T, AddressSpace> ptr)
{
    for (int i = 0; i < NumElements; ++i)
    {
        ptr.get()[offset * NumElements + i] = v[i];
    }
}

#elif GMX_SYCL_DPCPP

/*! \brief Polyfill for sycl::vec::load buggy in hipSYCL
 *
 * Loads from the address \c ptr offset in elements of type T by
 * NumElements * offset, into the components of \c v.
 *
 * Can probably be removed when
 * https://github.com/illuhad/hipSYCL/issues/647 is resolved. */
template<sycl::access::address_space AddressSpace, typename T, int NumElements>
static inline void loadToVec(size_t offset,
                             sycl::multi_ptr<const T, AddressSpace> ptr,
                             sycl::vec<T, NumElements>* v)
{
    v->load(offset, ptr);
}

/*! \brief Polyfill for sycl::vec::store buggy in hipSYCL
 *
 * Loads from the address \c ptr offset in elements of type T by
 * NumElements * offset, into the components of \c v.
 *
 * Can probably be removed when
 * https://github.com/illuhad/hipSYCL/issues/647 is resolved. */
template<sycl::access::address_space AddressSpace, typename T, int NumElements>
static inline void storeFromVec(const sycl::vec<T, NumElements>& v,
                                size_t offset,
                                sycl::multi_ptr<T, AddressSpace> ptr)
{
    v.store(offset, ptr);
}

#endif

} // namespace sycl_2020

#endif /* GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H */
