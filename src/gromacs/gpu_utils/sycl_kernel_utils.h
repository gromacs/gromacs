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

#include <config.h>

#include <type_traits>

#include "gmxsycl.h"
#include "gputraits_sycl.h"

#if defined(__AMDGCN__)
#    include "hip_sycl_kernel_utils.h"
#endif

#if defined(SYCL_EXT_ONEAPI_WORK_GROUP_STATIC)
#    include <sycl/ext/oneapi/work_group_static.hpp>
#endif

/*! \file
 *  \brief SYCL kernel helper functions.
 *
 *  \author Andrey Alekseenko <al42and@gmail.com>
 */

//! \brief Full warp active thread mask used in CUDA warp-level primitives.
static constexpr unsigned int c_cudaFullWarpMask = 0xffffffff;

// Prior to AdaptiveCpp 24.06, only HIPSYCL_* are defined
// TODO: Remove once we have AdaptiveCpp 24.06 as the minimum supported version
#if GMX_SYCL_ACPP && !defined(ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP) \
        && defined(HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP)
#    define ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP HIPSYCL_LIBKERNEL_IS_DEVICE_PASS_HIP
#endif
#if GMX_SYCL_ACPP && !defined(ACPP_UNIVERSAL_TARGET) && defined(HIPSYCL_UNIVERSAL_TARGET)
#    define ACPP_UNIVERSAL_TARGET HIPSYCL_UNIVERSAL_TARGET
#endif

#if defined(SYCL_EXT_ONEAPI_ASSERT) && SYCL_EXT_ONEAPI_ASSERT && !defined(__AMDGCN__)
#    define SYCL_ASSERT(condition) assert(condition)
#else
/* Assertions are not defined in SYCL standard, but they are available as oneAPI extension sycl_ext_oneapi_assert.
 * Technically, asserts should work just fine with AdaptiveCpp, since they are supported by CUDA since sm_20.
 * But with some settings (Clang 14, hipSYCL 0.9.2, RelWithAssert), CUDA build fails at link time:
 * ptxas fatal   : Unresolved extern function '__assert_fail'
 * So, we just disable kernel asserts unless they are promised to be available.
 * Also disable asserts on AMD GPUs, since they cause `hipErrorIllegalState` (oneAPI 2024.1-2024.2, ROCm 5.4-6.1) */
#    define SYCL_ASSERT(condition)
#endif

#if GMX_SYCL_ACPP && ACPP_LIBKERNEL_IS_DEVICE_PASS_HIP
ACPP_UNIVERSAL_TARGET
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
    // AdaptiveCpp/CUDA. For DPC++, any non-CPU targets must be explicitly allowed in the #if below.
#if GMX_SYCL_ACPP
#    if !GMX_ACPP_HAVE_GENERIC_TARGET
    __hipsycl_if_target_host(return true;);
#    else
    return false;
#    endif
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
    return expectedSubGroupSize == deviceWavefrontSize();
#elif defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))
    return true; // Assume that we have set reqd_sub_group_size attribute for the kernel
#elif GMX_ACPP_HAVE_GENERIC_TARGET
    return true;
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
#if (defined(__SYCL_DEVICE_ONLY__) && (defined(__SPIR__) || defined(__SPIRV__))) && !GMX_ACPP_HAVE_GENERIC_TARGET
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

template<typename T>
static inline void atomicFetchAddLocal(T* val, const T delta)
{
    atomicFetchAddLocal<T>(*val, delta);
}

// \brief Staggered atomic force component accumulation into global memory to reduce clashes
//
// Reduce the number of atomic clashes by a theoretical max 3x by having consecutive threads
// accumulate different force components at the same time.
static inline void staggeredAtomicAddForce(sycl::global_ptr<Float3> gm_f, Float3 f, const int localId)
{
    __builtin_assume(localId >= 0);

    using Int3  = sycl::int3;
    Int3 offset = { 0, 1, 2 };

    // Shift force components x (0), y (1), and z (2) left by 2, 1, and 0, respectively
    // to end up with zxy, yzx, xyz on consecutive threads.
    f      = (localId % 3 == 0) ? Float3(f[1], f[2], f[0]) : f;
    offset = (localId % 3 == 0) ? Int3(offset[1], offset[2], offset[0]) : offset;
    f      = (localId % 3 <= 1) ? Float3(f[1], f[2], f[0]) : f;
    offset = (localId % 3 <= 1) ? Int3(offset[1], offset[2], offset[0]) : offset;

    atomicFetchAdd(gm_f[0][offset[0]], f[0]);
    atomicFetchAdd(gm_f[0][offset[1]], f[1]);
    atomicFetchAdd(gm_f[0][offset[2]], f[2]);
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
     * a very conservative 4-byte limit here.
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

/*! \brief Shim to smooth over implementation differences
 *
 * This lets the classic oneAPI implementation and the handler-free
 * extension for workgroup-static (a.k.a. local accessor) storage work
 * smoothly in the same templated kernel when the memory is not
 * needed for the current template instantiation. */
template<typename T>
struct NullptrWrapper
{
};

#if GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION
static constexpr bool sc_useCommandGroupHandler = false;
#else // This is the default (and classic) implementation
static constexpr bool sc_useCommandGroupHandler = true;
#endif

namespace detail
{

/*! \brief Type trait for implementing StaticLocalStorage
 *
 * Specializations of this type allow non-trivial vector-style types
 * to work with work_group_static (which requires that the type of its
 * array is trivial). */
template<typename T>
struct UnderlyingType
{
    static constexpr size_t sc_valueTypeMultiple = 1;
    using Type                                   = T;
};

// Template specialization for sycl::vec
template<typename T, size_t size>
struct UnderlyingType<sycl::vec<T, size>>
{
    static constexpr size_t sc_valueTypeMultiple = size;
    using Type                                   = T;
};

// Template specialization for gmx::BasicVector
template<typename T>
struct UnderlyingType<gmx::BasicVector<T>>
{
    static constexpr size_t sc_valueTypeMultiple = DIM;
    using Type                                   = T;
};

} // namespace detail

/*! \brief Shim for the two implementations of local storage
 *
 * A sycl::local_accessor works as a facet of a command-group handler,
 * so when using handler-free submission, we must also use the
 * extension for work-group static storage, which does not depend on a
 * command-group handler.
 *
 * Because we have the above two implementation flavours for SYCL
 * device kernels, this class permits the code to specify the base
 * type, buffer size, and condition under which it exists in one
 * place, and then use that correctly for both implementations of
 * local storage.
 *
 * Note that the declaration of a sycl::local_accessor must be in
 * the host code, whereas the declaration of a work_group_static
 * must be in the kernel code.
 **/
template<typename T, size_t size, bool condition = true>
class StaticLocalStorage
{
public:
    //! Helper typedef
    using Underlying = detail::UnderlyingType<T>;
#if defined(SYCL_EXT_ONEAPI_WORK_GROUP_STATIC)
    //! Helper typedef for the static array that defines the storage of a work_group_static
    using StaticArray = typename Underlying::Type[size * Underlying::sc_valueTypeMultiple];
    //! Convenience typedef for a work_group_static
    using WorkGroupStatic = sycl::ext::oneapi::experimental::work_group_static<StaticArray>;
#else
    using WorkGroupStatic = std::nullptr_t;
#endif
    /*! \brief The type to declare local storage in host code, given
     * the build configuration
     *
     * With sc_useCommandGroupHandler, this type
     * (i.e. sycl::local_accessor) can only be successfully declared
     * in the host code. */
    using ActiveHostStorage =
            std::conditional_t<sc_useCommandGroupHandler, sycl::local_accessor<T, 1>, std::nullptr_t>;
    /*! \brief The type to declare local storage in device code, given
     * the build configuration
     *
     * With !sc_useCommandGroupHandler, this type (ie. work_group_static)
     * can only be successfully declared in the device kernel. */
    using ActiveDeviceStorage =
            std::conditional_t<sc_useCommandGroupHandler, std::nullptr_t, WorkGroupStatic>;
    //! Type of local storage in the host code, given \c condition
    using HostStorage = std::conditional_t<condition, ActiveHostStorage, std::nullptr_t>;
    //! Type of local storage in the device code, given \c condition
    using DeviceStorage = std::conditional_t<condition, ActiveDeviceStorage, std::nullptr_t>;
    /*! \brief Helper function to make storage for the local_accessor case
     *
     * May only be called on the host.
     *
     * Note that the equivalent approach for work_group_static
     * compiles and fails at runtime, because the address of the local
     * storage is invalid. */
    static HostStorage makeHostStorage(sycl::handler& cgh)
    {
        if constexpr (condition)
        {
            if constexpr (sc_useCommandGroupHandler)
            {
                return ActiveHostStorage{ sycl::range<1>(size), cgh };
            }
            else
            {
                return nullptr;
            }
        }
        else
        {
            return nullptr;
        }
    }
    /*! \brief Helper function to make storage for the work_group_static case
     *
     * Should only be called on the host, to match its other implementation. */
    static HostStorage makeHostStorage(std::nullptr_t /* cgh */) { return nullptr; }
    //! Helper typedef for get_pointer() implementation
    using LocalPtr = sycl::multi_ptr<T, sycl::access::address_space::local_space>;
    /*! \brief Provide a consistent way to get a sycl::local_ptr when using sycl::local_accessor
     *
     * Note that the sycl::local_accessor is created on the host and
     * is const in the device kernel, but the sycl::local_ptr obtained
     * from it can be either const or non-const according to need. */
    static LocalPtr get_pointer(const sycl::local_accessor<T, 1>& hostStorage,
                                std::nullptr_t& /* deviceStorage */) noexcept
    {
        return hostStorage.get_pointer();
    }
    /*! \brief Provide a consistent way to get a sycl::local_ptr when using a work_group_static
     *
     * Note that the work_group_static is created on the device and
     * can be treated as either const or non-const according to
     * need. But since the device kernels generally need read-write access, we
     * only need the non-const overload. */
#if defined(SYCL_EXT_ONEAPI_WORK_GROUP_STATIC)
    static LocalPtr get_pointer(const std::nullptr_t& /* hostStorage */, WorkGroupStatic& deviceStorage)
    {
        // Extract the pointer to raw value underlying T
        typename detail::UnderlyingType<T>::Type* rawData = deviceStorage;
        // Potentially reinterpret it as a pointer-to-vector, ie. T*
        T* maybeVectorData = reinterpret_cast<T*>(rawData);
        // Construct the sycl::local_ptr<T>
        return LocalPtr{ maybeVectorData };
    }
#endif

    /*! \brief Provide a consistent way to get a sycl::local_ptr
     * containing nullptr when the storage is empty. */
    static LocalPtr get_pointer(const std::nullptr_t& /* hostStorage */,
                                std::nullptr_t& /* deviceStorage */) noexcept
    {
        return nullptr;
    }
};

#endif /* GMX_GPU_UTILS_SYCL_KERNEL_UTILS_H */
