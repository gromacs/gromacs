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

// Prior to AdaptiveCpp 24.06, only HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP is defined
// TODO: Remove once we have AdaptiveCpp 24.06 as the minimum supported version
#if GMX_SYCL_ACPP && !defined(ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP) \
        && defined(HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP)
#    define ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP
#endif

#if defined(SYCL_EXT_ONEAPI_ASSERT) && SYCL_EXT_ONEAPI_ASSERT && !defined(__AMDGCN__)
#    define SYCL_ASSERT(condition) assert(condition)
#else
/* Assertions are not defined in SYCL standard, but they are available as oneAPI extension sycl_ext_oneapi_assert.
 * Technically, asserts should work just fine with hipSYCL, since they are supported by CUDA since sm_20.
 * But with some settings (Clang 14, hipSYCL 0.9.2, RelWithAssert), CUDA build fails at link time:
 * ptxas fatal   : Unresolved extern function '__assert_fail'
 * So, we just disable kernel asserts unless they are promised to be available.
 * Also disable asserts on AMD GPUs, since they cause `hipErrorIllegalState` (oneAPI 2024.1-2024.2, ROCm 5.4-6.1) */
#    define SYCL_ASSERT(condition)
#endif

#if GMX_SYCL_ACPP && ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP
HIPSYCL_UNIVERSAL_TARGET
static inline void atomicAddOptimizedAmd(float gmx_unused* ptr, const float gmx_unused delta)
{
#    if defined(__gfx908__) // Special function for AMD MI100
    CLANG_DIAGNOSTIC_IGNORE("-Wdeprecated-declarations")
    atomicAddNoRet(ptr, delta);
    CLANG_DIAGNOSTIC_RESET
#    elif defined(__gfx90a__) // Special function for AMD MI200
    unsafeAtomicAdd(ptr, delta);
#    else
    atomicAdd(ptr, delta);
#    endif
}
#endif


constexpr bool compilingForHost()
{
    // Skip compiling for CPU. Makes compiling this file ~10% faster for oneAPI/CUDA or
    // hipSYCL/CUDA. For DPC++, any non-CPU targets must be explicitly allowed in the #if below.
#if GMX_SYCL_ACPP
    __hipsycl_if_target_host(return true;);
#endif
    // We need to list all valid device targets here
#if GMX_SYCL_DPCPP \
        && !(defined(__NVPTX__) || defined(__AMDGCN__) || defined(__SPIR__) || defined(__SPIRV__))
    return true;
#endif
    return false;
}

template<int expectedSubGroupSize>
constexpr bool compilingForSubGroupSize()
{
#if defined(__SYCL_DEVICE_ONLY__) && defined(__NVPTX__)
    return expectedSubGroupSize == 32;
#elif defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
    return expectedSubGroupSize == __AMDGCN_WAVEFRONT_SIZE;
#elif defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))
    return true; // Assume that we have set reqd_sub_group_size attribute for the kernel
#else
    return false; // Unknown architecture
#endif
}

template<int expectedSubGroupSize>
constexpr bool skipKernelCompilation()
{
    if constexpr (compilingForHost())
    {
        return true;
    }
    if constexpr (!compilingForSubGroupSize<expectedSubGroupSize>())
    {
        return true;
    }
    /* Currently, the only SPIR-V target is Intel; for this we don't need 64-wide kernels.
     * This will require changing if we ever have other SPIR-V targets. */
#if defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))
    if constexpr (expectedSubGroupSize > 32)
    {
        return true;
    }
#endif
    return false;
}

template<typename T, sycl::memory_scope MemoryScope, sycl::access::address_space AddressSpace>
static inline void atomicAddDefault(T& val, const T delta)
{
    using sycl::memory_order;
    sycl::atomic_ref<T, memory_order::relaxed, MemoryScope, AddressSpace> ref(val);
    ref.fetch_add(delta);
}

/*! \brief Convenience wrapper to do atomic addition to a global buffer.
 */
template<typename T,
         sycl::memory_scope          MemoryScope  = sycl::memory_scope::device,
         sycl::access::address_space AddressSpace = sycl::access::address_space::global_space>
static inline void atomicFetchAdd(T& val, const T delta)
{
    using sycl::access::address_space;
    // Check if we need/can call the optimized atomicAdd for AMD devices, see #4465.
    if constexpr (GMX_SYCL_ACPP && std::is_same_v<T, float> && AddressSpace == address_space::global_space)
    {
#if GMX_SYCL_ACPP && ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP
        atomicAddOptimizedAmd(&val, delta);
#else
        atomicAddDefault<T, MemoryScope, AddressSpace>(val, delta);
#endif
    }
    else if constexpr (std::is_same_v<T, Float3>)
    {
        atomicFetchAdd<float, MemoryScope, AddressSpace>(val[XX], delta[XX]);
        atomicFetchAdd<float, MemoryScope, AddressSpace>(val[YY], delta[YY]);
        atomicFetchAdd<float, MemoryScope, AddressSpace>(val[ZZ], delta[ZZ]);
    }
    else
    {
        atomicAddDefault<T, MemoryScope, AddressSpace>(val, delta);
    }
}

template<typename T>
static inline void atomicFetchAddLocal(T& val, const T delta)
{
    atomicFetchAdd<T, sycl::memory_scope::work_group, sycl::access::address_space::local_space>(val, delta);
}

/*! \brief Convenience wrapper to do atomic loads from a global buffer.
 */
template<typename T, sycl::memory_scope MemoryScope = sycl::memory_scope::device>
static inline T atomicLoad(T& val)
{
#if GMX_SYCL_ACPP && GMX_ACPP_HAVE_CUDA_TARGET
    /* Some versions of Clang do not support atomicLoad in NVPTX backend, and die with ICE,
     * e.g. Clang 14, hipSYCL 0.9.2, CUDA 11.5, Debug:
     * fatal error: error in backend: Cannot select: 0xd870450: i32,ch = AtomicLoad<(load seq_cst (s32) from %ir.27)>.
     *
     * Since we use relaxed memory order, normal loads should be safe for small datatypes.
     * Allegedly, datatypes up to 128 bytes should be fine, but we only use floats, so we have
     * a very concervative 4-byte limit here.
     * As a demonstration of correctness, we don't use atomic loads in CUDA and it's doing fine.
     * See https://github.com/AdaptiveCpp/AdaptiveCpp/issues/752 */
    static_assert(sizeof(T) <= 4);
    return val;
#else
    using sycl::access::address_space;
    sycl::atomic_ref<T, sycl::memory_order::relaxed, MemoryScope, address_space::global_space> ref(val);
    return ref.load();
#endif
}

/*! \brief Issue an intra sub-group barrier.
 *
 * Equivalent with CUDA's \c syncwarp(c_cudaFullWarpMask).
 *
 */
template<int Dim>
static inline void subGroupBarrier(const sycl::nd_item<Dim> itemIdx)
{
#if GMX_SYCL_ACPP
    sycl::group_barrier(itemIdx.get_sub_group(), sycl::memory_scope::sub_group);
#else
    itemIdx.get_sub_group().barrier();
#endif
}

#if defined(__SYCL_DEVICE_ONLY__) && defined(__AMDGCN__)
/* !\brief Cross-lane move operation using AMD DPP (Data-Parallel Primitives).
 *
 * Uses the __builtin_amdgcn_update_dpp intrinsic which expressed the data
 * movement but the compiler will combine these with subsequent instructions
 * if possible.
 *
 * Note that this is a generic implementation for any type T (for current use
 * it could be more simple).
 *
 * Ref: https://gpuopen.com/learn/amd-gcn-assembly-cross-lane-operations
 */
template<class T, int dppCtrl, int rowMask = 0xf, int bankMask = 0xf, bool boundCtrl = true>
#    if GMX_SYCL_ACPP
__device__ __host__
#    endif
        __attribute__((always_inline)) T
        amdDppUpdateShfl(const T& input)
{
    constexpr int c_wordCount = (sizeof(T) + sizeof(int) - 1) / sizeof(int);

    struct V
    {
        int word[c_wordCount];
    };
    V wordList = __builtin_bit_cast(V, input);
#    pragma unroll
    for (int i = 0; i < c_wordCount; i++)
    {
        wordList.word[i] =
                __builtin_amdgcn_update_dpp(0, wordList.word[i], dppCtrl, rowMask, bankMask, boundCtrl);
    }

    return __builtin_bit_cast(T, wordList);
}
#endif

#if GMX_SYCL_ACPP && !(defined(ACPP_VERSION_MAJOR) && ACPP_VERSION_MAJOR >= 24)
namespace sycl
{
/*! \brief Popcount instruction for SYCL
 *
 * sycl::popcount is missing in AdaptiveCpp prior to 24.02.0.
 * This is very basic version using Clang built-in.
 */
template<typename T>
T popcount(T x)
{
    if constexpr (sizeof(T) == 4)
    {
        return __builtin_popcount(x);
    }
    else
    {
        static_assert(sizeof(T) == 8);
        return __builtin_popcountll(x);
    }
}
} // namespace sycl
#endif

#endif /* GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H */
