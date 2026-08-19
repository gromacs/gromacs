/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH
#define GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH

/*! \file
 *  \brief CUDA device util functions (callable from GPU kernels).
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"

/*! Load directly or using __ldg() when supported. */
template<typename T>
__device__ __forceinline__ T LDG(const T* ptr)
{
    /* CC >=3.5 supports constant loads through texture or L1 */
    return __ldg(ptr);
}

/*! \brief Fetch the value by \p index from the texture object.
 *
 * \tparam T            Raw data type
 * \param[in] texObj    Table texture object
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template<typename T>
static __forceinline__ __device__ T fetchFromTexture(const cudaTextureObject_t texObj, int index)
{
    assert(index >= 0);
    // NOLINTNEXTLINE(misc-static-assert)
    assert(!c_disableCudaTextures);
    return tex1Dfetch<T>(texObj, index);
}

/*! \brief Fetch the value by \p index from the parameter lookup table.
 *
 *  Depending on what is supported, it fetches parameters either
 *  using direct load or texture objects.
 *
 * \tparam T            Raw data type
 * \param[in] d_ptr     Device pointer to the raw table memory
 * \param[in] texObj    Table texture object
 * \param[in] index     Non-negative element index
 * \returns             The value from the table at \p index
 */
template<typename T>
static __forceinline__ __device__ T fetchFromParamLookupTable(const T*                  d_ptr,
                                                              const cudaTextureObject_t texObj,
                                                              int                       index)
{
    assert(index >= 0);
    T result;
#if DISABLE_CUDA_TEXTURES
    GMX_UNUSED_VALUE(texObj);
    result = LDG(d_ptr + index);
#else
    GMX_UNUSED_VALUE(d_ptr);
    result = fetchFromTexture<T>(texObj, index);
#endif
    return result;
}

/*!
 * \brief GPU-scoped acquire load of a 32-bit value from global memory.
 *
 * PTX: ld.acquire.gpu.global.u32
 *
 * Semantics:
 * - Acquire ordering at GPU scope: subsequent global reads/writes by this thread
 *   cannot be reordered before this load. Ensures visibility of prior writes by
 *   other thread blocks on the same GPU (e.g., grid-synchronization counters)
 *   before consuming dependent data.
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 *
 * \param[in] ptr  Address of the 32-bit value in global memory
 * \return The loaded 32-bit value
 */
inline __device__ uint32_t loadAcquireGpuAsm(const uint32_t* ptr)
{
    uint32_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.acquire.gpu.global.u32 %0, [%1];" : "=r"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped acquire load of a 64-bit value from global memory.
 *  we need system scoped ld.acquire as the signal reads though local to
 *  the gpu but are written by remote gpu.
 *
 * PTX: ld.acquire.sys.global.u64
 *
 * Semantics:
 * - Acquire ordering at system scope: subsequent global reads/writes by this thread
 *   cannot be reordered before this load. Ensures visibility of remote GPU writes
 *   (e.g., NVSHMEM signals) before consuming dependent data.
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 */
inline __device__ uint64_t loadAcquireSysAsm(const uint64_t* ptr)
{
    uint64_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.acquire.sys.global.u64 %0, [%1];" : "=l"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped relaxed load of a 64-bit value from global memory.
 *
 * PTX: ld.relaxed.sys.global.u64
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 */
inline __device__ uint64_t loadRelaxedSysAsm(const uint64_t* ptr)
{
    uint64_t retval = 0;
#if __CUDA_ARCH__ >= 700
    asm("ld.relaxed.sys.global.u64 %0, [%1];" : "=l"(retval) : "l"(ptr) : "memory");
#endif
    return retval;
}

/*!
 * \brief System-scoped release store of a 64-bit value to global memory.
 *
 * PTX: st.release.sys.global.u64
 *
 * Semantics:
 * - Release ordering at system scope: all prior global writes by this thread
 *   become visible to system peers before the store is observable. Used to
 *   publish completion signals after making packed data globally visible.
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 */
inline __device__ void storeReleaseSysAsm(uint64_t* ptr, const uint64_t val)
{
#if __CUDA_ARCH__ >= 700
    asm("st.release.sys.global.u64 [%0], %1;" : : "l"(ptr), "l"(val) : "memory");
#endif
}

/*!
 * \brief System-scoped relaxed store of a 64-bit value to global memory.
 *
 * PTX: st.relaxed.sys.global.u64
 *
 * Semantics:
 * - Relaxed ordering at system scope: does not impose ordering on prior global writes
 *   from this thread; only this store is guaranteed visible at system scope.
 * - Appropriate when no earlier data from this thread needs to be ordered with the signal.
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 */
inline __device__ void storeRelaxedSysAsm(uint64_t* ptr, const uint64_t val)
{
#if __CUDA_ARCH__ >= 700
    asm("st.relaxed.sys.global.u64 [%0], %1;" : : "l"(ptr), "l"(val) : "memory");
#endif
}


/*!
 * \brief Atomically increment counter with release semantics.
 *
 * PTX: atom.inc.release.gpu.global.u32 old, [addr], modMinusOne
 *
 * - The counter at addr is incremented modulo (modMinusOne + 1) and the OLD value is returned.
 * - The .release qualifier ensures that all prior global writes by this thread are made visible
 *   at GPU scope before the atomic is observed by others. This acts as a cache-flushing
 *   fence for data the threadBlock produced before signalling its arrival.
 * - We use modulo (numBlocks - 1) + 1 so that the last arriving block can be detected by
 *   comparing the returned old value against (numBlocks - 1).
 *
 * \warning Only supported with sm_70 and newer; no-op on older architectures!
 */
inline __device__ uint32_t atomicIncReleaseGpu(uint32_t* addr, int32_t modMinusOne)
{
    uint32_t old = 0;
#if __CUDA_ARCH__ >= 700
    asm("atom.inc.release.gpu.global.u32 %0,[%1],%2;"
        : "=r"(old)
        : "l"(addr), "r"(modMinusOne)
        : "memory");
#endif
    return old;
}

#endif /* GMX_GPU_UTILS_CUDA_KERNEL_UTILS_CUH */
