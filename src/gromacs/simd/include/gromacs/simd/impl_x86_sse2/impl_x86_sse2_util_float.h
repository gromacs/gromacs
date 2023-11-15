/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#ifndef GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H

#include "config.h"

#include <emmintrin.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_simd_float.h"
#include "gromacs/utility/basedefinitions.h"


namespace gmx
{

template<int align>
static inline void gmx_simdcall gatherLoadTranspose(const float*       base,
                                                    const std::int32_t offset[],
                                                    SimdFloat*         v0,
                                                    SimdFloat*         v1,
                                                    SimdFloat*         v2,
                                                    SimdFloat*         v3)
{
    assert(std::size_t(base + align * offset[0]) % 16 == 0);
    assert(std::size_t(base + align * offset[1]) % 16 == 0);
    assert(std::size_t(base + align * offset[2]) % 16 == 0);
    assert(std::size_t(base + align * offset[3]) % 16 == 0);

    v0->simdInternal_ = _mm_load_ps(base + align * offset[0]);
    v1->simdInternal_ = _mm_load_ps(base + align * offset[1]);
    v2->simdInternal_ = _mm_load_ps(base + align * offset[2]);
    v3->simdInternal_ = _mm_load_ps(base + align * offset[3]);

    _MM_TRANSPOSE4_PS(v0->simdInternal_, v1->simdInternal_, v2->simdInternal_, v3->simdInternal_);
}

template<int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float* base, const std::int32_t offset[], SimdFloat* v0, SimdFloat* v1)
{
    __m128 t1, t2;

    v0->simdInternal_ =
            _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(base + align * offset[0])));
    v1->simdInternal_ =
            _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(base + align * offset[1])));
    t1 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(base + align * offset[2])));
    t2 = _mm_castpd_ps(_mm_load_sd(reinterpret_cast<const double*>(base + align * offset[3])));
    t1 = _mm_unpacklo_ps(v0->simdInternal_, t1);
    t2 = _mm_unpacklo_ps(v1->simdInternal_, t2);
    v0->simdInternal_ = _mm_unpacklo_ps(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_ps(t1, t2);
}

static const int c_simdBestPairAlignmentFloat = 2;

template<int align>
static inline void gmx_simdcall gatherLoadUTranspose(const float*       base,
                                                     const std::int32_t offset[],
                                                     SimdFloat*         v0,
                                                     SimdFloat*         v1,
                                                     SimdFloat*         v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8;

    if (align % 4 != 0)
    {
        // general case, not aligned to 4-byte boundary
        t1 = _mm_loadu_ps(base + align * offset[0]);
        t2 = _mm_loadu_ps(base + align * offset[1]);
        t3 = _mm_loadu_ps(base + align * offset[2]);
        t4 = _mm_loadu_ps(base + align * offset[3]);
    }
    else
    {
        // aligned to 4-byte boundary or more
        t1 = _mm_load_ps(base + align * offset[0]);
        t2 = _mm_load_ps(base + align * offset[1]);
        t3 = _mm_load_ps(base + align * offset[2]);
        t4 = _mm_load_ps(base + align * offset[3]);
    }
    t5  = _mm_unpacklo_ps(t1, t2);
    t6  = _mm_unpacklo_ps(t3, t4);
    t7  = _mm_unpackhi_ps(t1, t2);
    t8  = _mm_unpackhi_ps(t3, t4);
    *v0 = _mm_movelh_ps(t5, t6);
    *v1 = _mm_movehl_ps(t6, t5);
    *v2 = _mm_movelh_ps(t7, t8);
}


template<int align>
static inline void gmx_simdcall
transposeScatterStoreU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    __m128 t1, t2;

    // general case, not aligned to 4-byte boundary
    t1 = _mm_unpacklo_ps(v0.simdInternal_, v1.simdInternal_);
    t2 = _mm_unpackhi_ps(v0.simdInternal_, v1.simdInternal_);
    _mm_storel_pi(reinterpret_cast<__m64*>(base + align * offset[0]), t1);
    _mm_store_ss(base + align * offset[0] + 2, v2.simdInternal_);
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[1]), t1);
    _mm_store_ss(base + align * offset[1] + 2,
                 _mm_shuffle_ps(v2.simdInternal_, v2.simdInternal_, _MM_SHUFFLE(1, 1, 1, 1)));
    _mm_storel_pi(reinterpret_cast<__m64*>(base + align * offset[2]), t2);
    _mm_store_ss(base + align * offset[2] + 2,
                 _mm_shuffle_ps(v2.simdInternal_, v2.simdInternal_, _MM_SHUFFLE(2, 2, 2, 2)));
    _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[3]), t2);
    _mm_store_ss(base + align * offset[3] + 2,
                 _mm_shuffle_ps(v2.simdInternal_, v2.simdInternal_, _MM_SHUFFLE(3, 3, 3, 3)));
}


template<int align>
static inline void gmx_simdcall
transposeScatterIncrU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

    if (align < 4)
    {
        t5  = _mm_unpacklo_ps(v1.simdInternal_, v2.simdInternal_);
        t6  = _mm_unpackhi_ps(v1.simdInternal_, v2.simdInternal_);
        t7  = _mm_shuffle_ps(v0.simdInternal_, t5, _MM_SHUFFLE(1, 0, 0, 0));
        t8  = _mm_shuffle_ps(v0.simdInternal_, t5, _MM_SHUFFLE(3, 2, 0, 1));
        t9  = _mm_shuffle_ps(v0.simdInternal_, t6, _MM_SHUFFLE(1, 0, 0, 2));
        t10 = _mm_shuffle_ps(v0.simdInternal_, t6, _MM_SHUFFLE(3, 2, 0, 3));

        t1 = _mm_load_ss(base + align * offset[0]);
        t1 = _mm_loadh_pi(t1, reinterpret_cast<__m64*>(base + align * offset[0] + 1));
        t1 = _mm_add_ps(t1, t7);
        _mm_store_ss(base + align * offset[0], t1);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[0] + 1), t1);

        t2 = _mm_load_ss(base + align * offset[1]);
        t2 = _mm_loadh_pi(t2, reinterpret_cast<__m64*>(base + align * offset[1] + 1));
        t2 = _mm_add_ps(t2, t8);
        _mm_store_ss(base + align * offset[1], t2);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[1] + 1), t2);

        t3 = _mm_load_ss(base + align * offset[2]);
        t3 = _mm_loadh_pi(t3, reinterpret_cast<__m64*>(base + align * offset[2] + 1));
        t3 = _mm_add_ps(t3, t9);
        _mm_store_ss(base + align * offset[2], t3);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[2] + 1), t3);

        t4 = _mm_load_ss(base + align * offset[3]);
        t4 = _mm_loadh_pi(t4, reinterpret_cast<__m64*>(base + align * offset[3] + 1));
        t4 = _mm_add_ps(t4, t10);
        _mm_store_ss(base + align * offset[3], t4);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[3] + 1), t4);
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations

        t1 = _mm_unpacklo_ps(v0.simdInternal_, v2.simdInternal_); // x0 z0 x1 z1
        t2 = _mm_unpackhi_ps(v0.simdInternal_, v2.simdInternal_); // x2 z2 x3 z3
        t3 = _mm_unpacklo_ps(v1.simdInternal_, _mm_setzero_ps()); // y0  0 y1  0
        t4 = _mm_unpackhi_ps(v1.simdInternal_, _mm_setzero_ps()); // y2  0 y3  0
        t5 = _mm_unpacklo_ps(t1, t3);                             // x0 y0 z0  0
        t6 = _mm_unpackhi_ps(t1, t3);                             // x1 y1 z1  0
        t7 = _mm_unpacklo_ps(t2, t4);                             // x2 y2 z2  0
        t8 = _mm_unpackhi_ps(t2, t4);                             // x3 y3 z3  0

        if (align % 4 == 0)
        {
            // alignment is a multiple of 4
            _mm_store_ps(base + align * offset[0], _mm_add_ps(_mm_load_ps(base + align * offset[0]), t5));
            _mm_store_ps(base + align * offset[1], _mm_add_ps(_mm_load_ps(base + align * offset[1]), t6));
            _mm_store_ps(base + align * offset[2], _mm_add_ps(_mm_load_ps(base + align * offset[2]), t7));
            _mm_store_ps(base + align * offset[3], _mm_add_ps(_mm_load_ps(base + align * offset[3]), t8));
        }
        else
        {
            // alignment >=5, but not a multiple of 4
            _mm_storeu_ps(base + align * offset[0],
                          _mm_add_ps(_mm_loadu_ps(base + align * offset[0]), t5));
            _mm_storeu_ps(base + align * offset[1],
                          _mm_add_ps(_mm_loadu_ps(base + align * offset[1]), t6));
            _mm_storeu_ps(base + align * offset[2],
                          _mm_add_ps(_mm_loadu_ps(base + align * offset[2]), t7));
            _mm_storeu_ps(base + align * offset[3],
                          _mm_add_ps(_mm_loadu_ps(base + align * offset[3]), t8));
        }
    }
}

template<int align>
static inline void gmx_simdcall
transposeScatterDecrU(float* base, const std::int32_t offset[], SimdFloat v0, SimdFloat v1, SimdFloat v2)
{
    // This implementation is identical to the increment version, apart from using subtraction instead
    __m128 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

    if (align < 4)
    {
        t5  = _mm_unpacklo_ps(v1.simdInternal_, v2.simdInternal_);
        t6  = _mm_unpackhi_ps(v1.simdInternal_, v2.simdInternal_);
        t7  = _mm_shuffle_ps(v0.simdInternal_, t5, _MM_SHUFFLE(1, 0, 0, 0));
        t8  = _mm_shuffle_ps(v0.simdInternal_, t5, _MM_SHUFFLE(3, 2, 0, 1));
        t9  = _mm_shuffle_ps(v0.simdInternal_, t6, _MM_SHUFFLE(1, 0, 0, 2));
        t10 = _mm_shuffle_ps(v0.simdInternal_, t6, _MM_SHUFFLE(3, 2, 0, 3));

        t1 = _mm_load_ss(base + align * offset[0]);
        t1 = _mm_loadh_pi(t1, reinterpret_cast<__m64*>(base + align * offset[0] + 1));
        t1 = _mm_sub_ps(t1, t7);
        _mm_store_ss(base + align * offset[0], t1);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[0] + 1), t1);

        t2 = _mm_load_ss(base + align * offset[1]);
        t2 = _mm_loadh_pi(t2, reinterpret_cast<__m64*>(base + align * offset[1] + 1));
        t2 = _mm_sub_ps(t2, t8);
        _mm_store_ss(base + align * offset[1], t2);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[1] + 1), t2);

        t3 = _mm_load_ss(base + align * offset[2]);
        t3 = _mm_loadh_pi(t3, reinterpret_cast<__m64*>(base + align * offset[2] + 1));
        t3 = _mm_sub_ps(t3, t9);
        _mm_store_ss(base + align * offset[2], t3);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[2] + 1), t3);

        t4 = _mm_load_ss(base + align * offset[3]);
        t4 = _mm_loadh_pi(t4, reinterpret_cast<__m64*>(base + align * offset[3] + 1));
        t4 = _mm_sub_ps(t4, t10);
        _mm_store_ss(base + align * offset[3], t4);
        _mm_storeh_pi(reinterpret_cast<__m64*>(base + align * offset[3] + 1), t4);
    }
    else
    {
        // Extra elements means we can use full width-4 load/store operations

        t1 = _mm_unpacklo_ps(v0.simdInternal_, v2.simdInternal_); // x0 z0 x1 z1
        t2 = _mm_unpackhi_ps(v0.simdInternal_, v2.simdInternal_); // x2 z2 x3 z3
        t3 = _mm_unpacklo_ps(v1.simdInternal_, _mm_setzero_ps()); // y0  0 y1  0
        t4 = _mm_unpackhi_ps(v1.simdInternal_, _mm_setzero_ps()); // y2  0 y3  0
        t5 = _mm_unpacklo_ps(t1, t3);                             // x0 y0 z0  0
        t6 = _mm_unpackhi_ps(t1, t3);                             // x1 y1 z1  0
        t7 = _mm_unpacklo_ps(t2, t4);                             // x2 y2 z2  0
        t8 = _mm_unpackhi_ps(t2, t4);                             // x3 y3 z3  0

        if (align % 4 == 0)
        {
            // alignment is a multiple of 4
            _mm_store_ps(base + align * offset[0], _mm_sub_ps(_mm_load_ps(base + align * offset[0]), t5));
            _mm_store_ps(base + align * offset[1], _mm_sub_ps(_mm_load_ps(base + align * offset[1]), t6));
            _mm_store_ps(base + align * offset[2], _mm_sub_ps(_mm_load_ps(base + align * offset[2]), t7));
            _mm_store_ps(base + align * offset[3], _mm_sub_ps(_mm_load_ps(base + align * offset[3]), t8));
        }
        else
        {
            // alignment >=5, but not a multiple of 4
            _mm_storeu_ps(base + align * offset[0],
                          _mm_sub_ps(_mm_loadu_ps(base + align * offset[0]), t5));
            _mm_storeu_ps(base + align * offset[1],
                          _mm_sub_ps(_mm_loadu_ps(base + align * offset[1]), t6));
            _mm_storeu_ps(base + align * offset[2],
                          _mm_sub_ps(_mm_loadu_ps(base + align * offset[2]), t7));
            _mm_storeu_ps(base + align * offset[3],
                          _mm_sub_ps(_mm_loadu_ps(base + align * offset[3]), t8));
        }
    }
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline void gmx_simdcall expandScalarsToTriplets(SimdFloat  scalar,
                                                        SimdFloat* triplets0,
                                                        SimdFloat* triplets1,
                                                        SimdFloat* triplets2)
{
    triplets0->simdInternal_ =
            _mm_shuffle_ps(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE(1, 0, 0, 0));
    triplets1->simdInternal_ =
            _mm_shuffle_ps(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE(2, 2, 1, 1));
    triplets2->simdInternal_ =
            _mm_shuffle_ps(scalar.simdInternal_, scalar.simdInternal_, _MM_SHUFFLE(3, 3, 3, 2));
}
#endif


template<int align>
static inline void gmx_simdcall gatherLoadBySimdIntTranspose(const float* base,
                                                             SimdFInt32   offset,
                                                             SimdFloat*   v0,
                                                             SimdFloat*   v1,
                                                             SimdFloat*   v2,
                                                             SimdFloat*   v3)
{
    // For present-generation x86 CPUs it appears to be faster to simply
    // store the SIMD integer to memory and then use the normal load operations.
    // This is likely because (a) the extract function is expensive, and (b)
    // the alignment scaling can often be done as part of the load instruction
    // (which is even cheaper than doing it in SIMD registers).
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_FINT32_WIDTH];
    _mm_store_si128(reinterpret_cast<__m128i*>(ioffset), offset.simdInternal_);
    gatherLoadTranspose<align>(base, ioffset, v0, v1, v2, v3);
}

template<int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float* base, SimdFInt32 offset, SimdFloat* v0, SimdFloat* v1)
{
    // For present-generation x86 CPUs it appears to be faster to simply
    // store the SIMD integer to memory and then use the normal load operations.
    // This is likely because (a) the extract function is expensive, and (b)
    // the alignment scaling can often be done as part of the load instruction
    // (which is even cheaper than doing it in SIMD registers).
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_FINT32_WIDTH];
    _mm_store_si128(reinterpret_cast<__m128i*>(ioffset), offset.simdInternal_);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}


template<int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float* base, SimdFInt32 offset, SimdFloat* v0, SimdFloat* v1)
{
    // For present-generation x86 CPUs it appears to be faster to simply
    // store the SIMD integer to memory and then use the normal load operations.
    // This is likely because (a) the extract function is expensive, and (b)
    // the alignment scaling can often be done as part of the load instruction
    // (which is even cheaper than doing it in SIMD registers).
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t ioffset[GMX_SIMD_FINT32_WIDTH];
    _mm_store_si128(reinterpret_cast<__m128i*>(ioffset), offset.simdInternal_);
    gatherLoadTranspose<align>(base, ioffset, v0, v1);
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline float gmx_simdcall reduceIncr4ReturnSum(float* m, SimdFloat v0, SimdFloat v1, SimdFloat v2, SimdFloat v3)
{
    _MM_TRANSPOSE4_PS(v0.simdInternal_, v1.simdInternal_, v2.simdInternal_, v3.simdInternal_);
    v0.simdInternal_ = _mm_add_ps(v0.simdInternal_, v1.simdInternal_);
    v2.simdInternal_ = _mm_add_ps(v2.simdInternal_, v3.simdInternal_);
    v0.simdInternal_ = _mm_add_ps(v0.simdInternal_, v2.simdInternal_);
    v2.simdInternal_ = _mm_add_ps(v0.simdInternal_, _mm_load_ps(m));

    assert(std::size_t(m) % 16 == 0);
    _mm_store_ps(m, v2.simdInternal_);

    __m128 b = _mm_add_ps(v0.simdInternal_,
                          _mm_shuffle_ps(v0.simdInternal_, v0.simdInternal_, _MM_SHUFFLE(1, 0, 3, 2)));
    b        = _mm_add_ss(b, _mm_shuffle_ps(b, b, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float*>(&b);
}
#endif

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_UTIL_FLOAT_H
