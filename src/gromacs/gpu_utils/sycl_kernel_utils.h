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

#if defined(SYCL_EXT_ONEAPI_ASSERT) && SYCL_EXT_ONEAPI_ASSERT
#    define SYCL_ASSERT(condition) assert(condition)
#else
/* Assertions are not defined in SYCL standard, but they are available as oneAPI extension sycl_ext_oneapi_assert.
 * Technically, asserts should work just fine with hipSYCL, since they are supported by CUDA since sm_20.
 * But with some settings (Clang 14, hipSYCL 0.9.2, RelWithAssert), CUDA build fails at link time:
 * ptxas fatal   : Unresolved extern function '__assert_fail'
 * So, we just disable kernel asserts unless they are promised to be available. */
#    define SYCL_ASSERT(condition)
#endif

#if GMX_SYCL_HIPSYCL && HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP
HIPSYCL_UNIVERSAL_TARGET
static inline void atomicAddOptimizedAmd(float gmx_unused* ptr, const float gmx_unused delta)
{
#    if defined(__gfx908__) // Special function for AMD MI100
#        pragma clang diagnostic push
#        pragma clang diagnostic ignored "-Wdeprecated-declarations"
    atomicAddNoRet(ptr, delta);
#        pragma clang diagnostic pop
#    elif defined(__gfx90a__) // Special function for AMD MI200
    unsafeAtomicAdd(ptr, delta); // Not checked on real hardware, see #4465
#    else
    atomicAdd(ptr, delta);
#    endif
}
#endif


constexpr bool compilingForHost()
{
    // Skip compiling for CPU. Makes compiling this file ~10% faster for oneAPI/CUDA or
    // hipSYCL/CUDA. For DPC++, any non-CPU targets must be explicitly allowed in the #if below.
#if GMX_SYCL_HIPSYCL
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
    if constexpr (GMX_SYCL_HIPSYCL && std::is_same_v<T, float> && AddressSpace == address_space::global_space)
    {
#if GMX_SYCL_HIPSYCL && HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP
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
#if GMX_SYCL_HIPSYCL && GMX_HIPSYCL_HAVE_CUDA_TARGET
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

/*! \brief Special packed Float3 flavor to help compiler optimizations on AMD CDNA2 devices.
 *
 * Full FP32 performance of AMD CDNA2 devices, like MI200-series, can only be achieved
 * when operating on float2, in a SIMD2-fashion. Compiler (at least up to ROCm 5.6)
 * can use packed math automatically for normal Float3, but generates a lot of
 * data movement between normal and packed registers. Using this class helps avoid
 * this problem.
 *
 * The approach is based on the similar solution used by AMD and StreamHPC in their port.
 *
 * Currently only used in NBNXM kernels if GMX_NBNXM_ENABLE_PACKED_FLOAT3 is enabled
 *
 * \todo This class shall be removed as soon as the compiler is improved.
 *
 * See issue #4854 for more details.
 */
struct AmdPackedFloat3
{
    typedef float __attribute__((ext_vector_type(2))) Native_float2_;

    /* According to C++ standard, we should give names to all
     * the types and fields declared below. This, however, makes
     * this code very verbose, and harms readability in a major
     * way while this code is aimed to be used in a pretty niche
     * case with relatively small selection of compilers
     * (flavors of Clang 14-18, maybe later). So, we prefer
     * to disable the warnings. */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wnested-anon-types"
#pragma clang diagnostic ignored "-Wgnu-anonymous-struct"
    struct __attribute__((packed))
    {
        union
        {
            Native_float2_ xy_;
            struct
            {
                float x_, y_;
            };
        };
        float z_;
    };
#pragma clang diagnostic pop
    template<typename Index>
    __attribute__((always_inline)) float operator[](Index i) const
    {
        switch (i)
        {
            case 0: return xy_.x;
            case 1: return xy_.y;
            default: SYCL_ASSERT(i == 2); return z_;
        }
    }
    template<typename Index>
    __attribute__((always_inline)) float& operator[](Index i)
    {
        switch (i)
        {
            case 0: return x_;
            case 1: return y_;
            default: SYCL_ASSERT(i == 2); return z_;
        }
    }

    __attribute__((always_inline)) float          x() const { return xy_.x; }
    __attribute__((always_inline)) float          y() const { return xy_.y; }
    __attribute__((always_inline)) Native_float2_ xy() const { return xy_; }
    __attribute__((always_inline)) float          z() const { return z_; }

    AmdPackedFloat3() = default;

    AmdPackedFloat3(float x, float y, float z) : xy_{ x, y }, z_{ z } {}

    AmdPackedFloat3(Native_float2_ xy, float z) : xy_{ xy }, z_{ z } {}

    AmdPackedFloat3(Float3 r) : xy_{ r[0], r[1] }, z_{ r[2] } {}

    explicit operator Float3() const { return Float3{ xy_.x, xy_.y, z_ }; }

    __attribute__((always_inline)) AmdPackedFloat3& operator=(const AmdPackedFloat3& x)
    {
        xy_ = x.xy_;
        z_  = x.z_;
        return *this;
    }

    //! Allow inplace addition for AmdPackedFloat3
    __attribute__((always_inline)) AmdPackedFloat3& operator+=(const AmdPackedFloat3& right)
    {
        return *this = *this + right;
    }
    //! Allow inplace subtraction for AmdPackedFloat3
    __attribute__((always_inline)) AmdPackedFloat3& operator-=(const AmdPackedFloat3& right)
    {
        return *this = *this - right;
    }
    //! Allow vector addition
    __attribute__((always_inline)) AmdPackedFloat3 operator+(const AmdPackedFloat3& right) const
    {
        return { xy_ + right.xy(), z_ + right.z() };
    }
    //! Allow vector subtraction
    __attribute__((always_inline)) AmdPackedFloat3 operator-(const AmdPackedFloat3& right) const
    {
        return { xy_ - right.xy(), z_ - right.z() };
    }
    //! Scale vector by a scalar
    __attribute__((always_inline)) AmdPackedFloat3& operator*=(const float& right)
    {
        xy_ *= right;
        z_ *= right;
        return *this;
    }

    //! Length^2 of vector
    __attribute__((always_inline)) float norm2() const { return dot(*this); }

    //! Return dot product
    __attribute__((always_inline)) float dot(const AmdPackedFloat3& right) const
    {
        return x() * right.x() + y() * right.y() + z() * right.z();
    }
};
static_assert(sizeof(AmdPackedFloat3) == 12);

__attribute__((always_inline)) static AmdPackedFloat3 operator*(const AmdPackedFloat3& v, const float& s)
{
    return { v.xy() * s, v.z() * s };
}
__attribute__((always_inline)) static AmdPackedFloat3 operator*(const float& s, const AmdPackedFloat3& v)
{
    return { v.xy() * s, v.z() * s };
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
#    if GMX_SYCL_HIPSYCL
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

#endif /* GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H */
