/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_avx_512_general.h"
#include "impl_x86_avx_512_simd_float.h"

namespace gmx
{

static const int c_simdBestPairAlignmentFloat = 2;

namespace
{
// Multiply function optimized for powers of 2, for which it is done by
// shifting. Currently up to 8 is accelerated. Could be accelerated for any
// number with a constexpr log2 function.
template<int n>
SimdFInt32 fastMultiply(SimdFInt32 x)
{
    if (n == 2)
    {
        return _mm256_slli_epi32(x.simdInternal_, 1);
    }
    else if (n == 4)
    {
        return _mm256_slli_epi32(x.simdInternal_, 2);
    }
    else if (n == 8)
    {
        return _mm256_slli_epi32(x.simdInternal_, 3);
    }
    else
    {
        return x * n;
    }
}

template<int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *, SimdFInt32)
{
    //Nothing to do. Termination of recursion.
}
}

template <int align, typename ... Targs>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *base, SimdFInt32 offset, SimdFloat *v, Targs... Fargs)
{
    // For align 1 or 2: No multiplication of offset is needed
    if (align > 2)
    {
        offset = fastMultiply<align>(offset);
    }
    // For align 2: Scale of 2*sizeof(float) is used (maximum supported scale)
    constexpr int align_ = (align > 2) ? 1 : align;
    v->simdInternal_ = _mm256_i32gather_ps(base, offset.simdInternal_, sizeof(float)*align_);
    // Gather remaining elements. Avoid extra multiplication (new align is 1 or 2).
    gatherLoadBySimdIntTranspose<align_>(base+1, offset, Fargs ...);
}

template <int align, typename ... Targs>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float *base, SimdFInt32 offset, Targs... Fargs)
{
    gatherLoadBySimdIntTranspose<align>(base, offset, Fargs ...);
}

template <int align, typename ... Targs>
static inline void gmx_simdcall
gatherLoadTranspose(const float *base, const std::int32_t offset[], Targs... Fargs)
{
    gatherLoadBySimdIntTranspose<align>(base, simdLoad(offset, SimdFInt32Tag()), Fargs ...);
}

template <int align, typename ... Targs>
static inline void gmx_simdcall
gatherLoadUTranspose(const float *base, const std::int32_t offset[], Targs... Fargs)
{
    gatherLoadTranspose<align>(base, offset, Fargs ...);
}

template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float *              base,
                       const std::int32_t   offset[],
                       SimdFloat            v0,
                       SimdFloat            v1,
                       SimdFloat            v2)
{
    SimdFInt32 simdoffset = simdLoad(offset, SimdFInt32Tag());
    if (align > 2)
    {
        simdoffset = fastMultiply<align>(simdoffset);
    }
    constexpr size_t scale = (align > 2) ? sizeof(float) : sizeof(float) * align;

    _mm256_i32scatter_ps(base,       simdoffset.simdInternal_, v0.simdInternal_, scale);
    _mm256_i32scatter_ps(&(base[1]), simdoffset.simdInternal_, v1.simdInternal_, scale);
    _mm256_i32scatter_ps(&(base[2]), simdoffset.simdInternal_, v2.simdInternal_, scale);
}

template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    __m256 t[4], t5, t6, t7, t8;
    int    i;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    o[8];
    store(o, fastMultiply<align>(simdLoad(offset, SimdFInt32Tag())));
    if (align < 4)
    {
        t5   = _mm256_unpacklo_ps(v0.simdInternal_, v1.simdInternal_);
        t6   = _mm256_unpackhi_ps(v0.simdInternal_, v1.simdInternal_);
        t[0] = _mm256_shuffle_ps(t5, v2.simdInternal_, _MM_SHUFFLE(0, 0, 1, 0));
        t[1] = _mm256_shuffle_ps(t5, v2.simdInternal_, _MM_SHUFFLE(1, 1, 3, 2));
        t[2] = _mm256_shuffle_ps(t6, v2.simdInternal_, _MM_SHUFFLE(2, 2, 1, 0));
        t[3] = _mm256_shuffle_ps(t6, v2.simdInternal_, _MM_SHUFFLE(3, 3, 3, 2));
        for (i = 0; i < 4; i++)
        {
            _mm_mask_storeu_ps(base + o[i], avx256Int2Mask(7),
                               _mm_add_ps(_mm_loadu_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
            _mm_mask_storeu_ps(base + o[4 + i], avx256Int2Mask(7), 
                               _mm_add_ps(_mm_loadu_ps(base + o[4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
        }
    }
    else
    {
        //One could use shuffle here too if it is OK to overwrite the padded elements for alignment
        t5    = _mm256_unpacklo_ps(v0.simdInternal_, v2.simdInternal_);
        t6    = _mm256_unpackhi_ps(v0.simdInternal_, v2.simdInternal_);
        t7    = _mm256_unpacklo_ps(v1.simdInternal_, _mm256_setzero_ps());
        t8    = _mm256_unpackhi_ps(v1.simdInternal_, _mm256_setzero_ps());
        t[0]  = _mm256_unpacklo_ps(t5, t7);                             // x0 y0 z0  0 | x4 y4 z4 0
        t[1]  = _mm256_unpackhi_ps(t5, t7);                             // x1 y1 z1  0 | x5 y5 z5 0
        t[2]  = _mm256_unpacklo_ps(t6, t8);                             // x2 y2 z2  0 | x6 y6 z6 0
        t[3]  = _mm256_unpackhi_ps(t6, t8);                             // x3 y3 z3  0 | x7 y7 z7 0
        if (align % 4 == 0)
        {
            for (i = 0; i < 4; i++)
            {
                _mm_store_ps(base + o[i], _mm_add_ps(_mm_load_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
                _mm_store_ps(base + o[ 4 + i],
                             _mm_add_ps(_mm_load_ps(base + o[ 4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
            }
        }
        else
        {
            for (i = 0; i < 4; i++)
            {
                _mm_storeu_ps(base + o[i], _mm_add_ps(_mm_loadu_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
                _mm_storeu_ps(base + o[ 4 + i],
                              _mm_add_ps(_mm_loadu_ps(base + o[ 4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
            }
        }
    }
}

template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    __m256 t[4], t5, t6, t7, t8;
    int    i;
    alignas(GMX_SIMD_ALIGNMENT) std::int32_t    o[8];
    store(o, fastMultiply<align>(simdLoad(offset, SimdFInt32Tag())));
    if (align < 4)
    {
        t5   = _mm256_unpacklo_ps(v0.simdInternal_, v1.simdInternal_);
        t6   = _mm256_unpackhi_ps(v0.simdInternal_, v1.simdInternal_);
        t[0] = _mm256_shuffle_ps(t5, v2.simdInternal_, _MM_SHUFFLE(0, 0, 1, 0));
        t[1] = _mm256_shuffle_ps(t5, v2.simdInternal_, _MM_SHUFFLE(1, 1, 3, 2));
        t[2] = _mm256_shuffle_ps(t6, v2.simdInternal_, _MM_SHUFFLE(2, 2, 1, 0));
        t[3] = _mm256_shuffle_ps(t6, v2.simdInternal_, _MM_SHUFFLE(3, 3, 3, 2));
        for (i = 0; i < 4; i++)
        {
            _mm_mask_storeu_ps(base + o[i], avx256Int2Mask(7), 
                                          _mm_sub_ps(_mm_loadu_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
            _mm_mask_storeu_ps(base + o[ 4 + i], avx256Int2Mask(7), 
                                          _mm_sub_ps(_mm_loadu_ps(base + o[ 4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
        }
    }
    else
    {
        //One could use shuffle here too if it is OK to overwrite the padded elements for alignment
        t5    = _mm256_unpacklo_ps(v0.simdInternal_, v2.simdInternal_);
        t6    = _mm256_unpackhi_ps(v0.simdInternal_, v2.simdInternal_);
        t7    = _mm256_unpacklo_ps(v1.simdInternal_, _mm256_setzero_ps());
        t8    = _mm256_unpackhi_ps(v1.simdInternal_, _mm256_setzero_ps());
        t[0]  = _mm256_unpacklo_ps(t5, t7);                             // x0 y0 z0  0 | x4 y4 z4 0
        t[1]  = _mm256_unpackhi_ps(t5, t7);                             // x1 y1 z1  0 | x5 y5 z5 0
        t[2]  = _mm256_unpacklo_ps(t6, t8);                             // x2 y2 z2  0 | x6 y6 z6 0
        t[3]  = _mm256_unpackhi_ps(t6, t8);                             // x3 y3 z3  0 | x7 y7 z7 0
        if (align % 4 == 0)
        {
            for (i = 0; i < 4; i++)
            {
                _mm_store_ps(base + o[i], _mm_sub_ps(_mm_load_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
                _mm_store_ps(base + o[ 4 + i],
                             _mm_sub_ps(_mm_load_ps(base + o[ 4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
            }
        }
        else
        {
            for (i = 0; i < 4; i++)
            {
                _mm_storeu_ps(base + o[i], _mm_sub_ps(_mm_loadu_ps(base + o[i]), _mm256_castps256_ps128(t[i])));
                _mm_storeu_ps(base + o[ 4 + i],
                              _mm_sub_ps(_mm_loadu_ps(base + o[ 4 + i]), _mm256_extractf32x4_ps(t[i], 1)));
            }
        }
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    triplets0->simdInternal_ = _mm256_permutexvar_ps(_mm256_set_epi32(2, 2, 1, 1, 1, 0, 0, 0),
                                                     scalar.simdInternal_);
    triplets1->simdInternal_ = _mm256_permutexvar_ps(_mm256_set_epi32(5, 4, 4, 4, 3, 3, 3, 2),
                                                     scalar.simdInternal_);
    triplets2->simdInternal_ = _mm256_permutexvar_ps(_mm256_set_epi32(7, 7, 7, 6, 6, 6, 5, 5),
                                                     scalar.simdInternal_);
}


static inline float gmx_simdcall
reduceIncr4ReturnSum(float *    m,
                     SimdFloat  v0,
                     SimdFloat  v1,
                     SimdFloat  v2,
                     SimdFloat  v3)
{
        __m128 t0, t2;

    assert(std::size_t(m) % 16 == 0);

    v0.simdInternal_ = _mm256_hadd_ps(v0.simdInternal_, v1.simdInternal_);
    v2.simdInternal_ = _mm256_hadd_ps(v2.simdInternal_, v3.simdInternal_);
    v0.simdInternal_ = _mm256_hadd_ps(v0.simdInternal_, v2.simdInternal_);
    t0               = _mm_add_ps(_mm256_castps256_ps128(v0.simdInternal_), _mm256_extractf128_ps(v0.simdInternal_, 0x1));

    t2 = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t2);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

static inline SimdFloat gmx_simdcall
loadDualHsimd(const float * m0,
              const float * m1)
{
    assert(std::size_t(m0) % 16 == 0);
    assert(std::size_t(m1) % 16 == 0);

    return {
        _mm256_insertf32x4(_mm256_castps128_ps256(_mm_load_ps(m0)), _mm_load_ps(m1), 1)
    };
}

static inline SimdFloat gmx_simdcall
loadDuplicateHsimd(const float * m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm256_broadcast_f32x4(_mm_load_ps(m))
    };
}

static inline SimdFloat gmx_simdcall
loadU1DualHsimd(const float * m)
{
    __m128 t0, t1;
    t0 = _mm_broadcast_ss(m);
    t1 = _mm_broadcast_ss(m+1);
    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t1, 0x1)
    };
}


static inline void gmx_simdcall
storeDualHsimd(float *     m0,
               float *     m1,
               SimdFloat   a)
{
    assert(std::size_t(m0) % 16 == 0);
    assert(std::size_t(m1) % 16 == 0);

    _mm_store_ps(m0, _mm256_castps256_ps128(a.simdInternal_));
    _mm_store_ps(m1, _mm256_extractf32x4_ps(a.simdInternal_, 1));
}

static inline void gmx_simdcall
incrDualHsimd(float *     m0,
              float *     m1,
              SimdFloat   a)
{
    assert(std::size_t(m0) % 16 == 0);
    assert(std::size_t(m1) % 16 == 0);

    __m128 x;

    // Lower half
    x = _mm_load_ps(m0);
    x = _mm_add_ps(x, _mm256_castps256_ps128(a.simdInternal_));
    _mm_store_ps(m0, x);

    // Upper half
    x = _mm_load_ps(m1);
    x = _mm_add_ps(x, _mm256_extractf32x4_ps(a.simdInternal_, 1));
    _mm_store_ps(m1, x);
}

static inline void gmx_simdcall
decrHsimd(float *    m,
          SimdFloat  a)
{
    assert(std::size_t(m) % 16 == 0);
    __m128 asum = _mm_add_ps(_mm256_castps256_ps128(a.simdInternal_), _mm256_extractf128_ps(a.simdInternal_, 0x1));
    _mm_store_ps(m, _mm_sub_ps(_mm_load_ps(m), asum));
}


template <int align>
static inline void gmx_simdcall
gatherLoadTransposeHsimd(const float *        base0,
                         const float *        base1,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1)
{
    __m128 t0, t1, t2, t3, t4, t5, t6, t7;
    __m256 tA, tB, tC, tD;

    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base0) % 8 == 0);
    assert(std::size_t(base1) % 8 == 0);
    assert(align % 2 == 0);

    t0  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[0]));
    t1  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[1]));
    t2  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[2]));
    t3  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base0 + align * offset[3]));
    t4  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[0]));
    t5  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[1]));
    t6  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[2]));
    t7  = _mm_loadl_pi(_mm_setzero_ps(), reinterpret_cast<const __m64 *>(base1 + align * offset[3]));

    tA  = _mm256_insertf128_ps(_mm256_castps128_ps256(t0), t4, 0x1);
    tB  = _mm256_insertf128_ps(_mm256_castps128_ps256(t1), t5, 0x1);
    tC  = _mm256_insertf128_ps(_mm256_castps128_ps256(t2), t6, 0x1);
    tD  = _mm256_insertf128_ps(_mm256_castps128_ps256(t3), t7, 0x1);

    tA                = _mm256_unpacklo_ps(tA, tC);
    tB                = _mm256_unpacklo_ps(tB, tD);
    v0->simdInternal_ = _mm256_unpacklo_ps(tA, tB);
    v1->simdInternal_ = _mm256_unpackhi_ps(tA, tB);
}

static inline float gmx_simdcall
reduceIncr4ReturnSumHsimd(float *     m,
                          SimdFloat   v0,
                          SimdFloat   v1)
{
    __m128 t0, t1;

    v0.simdInternal_ = _mm256_hadd_ps(v0.simdInternal_, v1.simdInternal_);
    t0               = _mm256_extractf128_ps(v0.simdInternal_, 0x1);
    t0               = _mm_hadd_ps(_mm256_castps256_ps128(v0.simdInternal_), t0);
    t0               = _mm_permute_ps(t0, _MM_SHUFFLE(3, 1, 2, 0));

    assert(std::size_t(m) % 16 == 0);

    t1   = _mm_add_ps(t0, _mm_load_ps(m));
    _mm_store_ps(m, t1);

    t0 = _mm_add_ps(t0, _mm_permute_ps(t0, _MM_SHUFFLE(1, 0, 3, 2)));
    t0 = _mm_add_ss(t0, _mm_permute_ps(t0, _MM_SHUFFLE(0, 3, 2, 1)));
    return *reinterpret_cast<float *>(&t0);
}

static inline SimdFloat gmx_simdcall
loadU4NOffset(const float* f, int offset)
{
    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_loadu_ps(f)), _mm_loadu_ps(f+offset), 0x1)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_UTIL_FLOAT_H
