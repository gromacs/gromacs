/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_AVX_512_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_512_UTIL_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_avx_512_general.h"
#include "impl_x86_avx_512_simd_double.h"

namespace gmx
{

// On MIC it is better to use scatter operations, so we define the load routines
// that use a SIMD offset variable first.

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *  base,
                             SimdDInt32      simdoffset,
                             SimdDouble *    v0,
                             SimdDouble *    v1,
                             SimdDouble *    v2,
                             SimdDouble *    v3)
{
    assert(std::size_t(base) % 32 == 0);
    assert(align % 4 == 0);

    // All instructions might be latency ~4 on MIC, so we use shifts where we
    // only need a single instruction (since the shift parameter is an immediate),
    // but multiplication otherwise.
    if (align == 4)
    {
        simdoffset = simdoffset << 2;
    }
    else if (align == 8)
    {
        simdoffset = simdoffset << 3;
    }
    else
    {
        simdoffset = simdoffset * SimdDInt32(align);
    }

    // The 8 corresponds to sizeof(double), but debug builds with gcc require
    // an immediate since sizeof() is taken as a function.
    v0->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base,   8);
    v1->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+1, 8);
    v2->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+2, 8);
    v3->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+3, 8);
}

template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const double *  base,
                              SimdDInt32      simdoffset,
                              SimdDouble *    v0,
                              SimdDouble *    v1)
{
    // All instructions might be latency ~4 on MIC, so we use shifts where we
    // only need a single instruction (since the shift parameter is an immediate),
    // but multiplication otherwise.
    if (align == 2)
    {
        simdoffset = simdoffset << 1;
    }
    else if (align == 4)
    {
        simdoffset = simdoffset << 2;
    }
    else if (align == 8)
    {
        simdoffset = simdoffset << 3;
    }
    else
    {
        simdoffset = simdoffset * SimdDInt32(align);
    }

    v0->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base,   8);
    v1->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+1, 8);
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const double *    base,
                             SimdDInt32        simdoffset,
                             SimdDouble *      v0,
                             SimdDouble *      v1)
{
    assert(std::size_t(base) % 16 == 0);
    assert(align % 2 == 0);
    gatherLoadUBySimdIntTranspose<align>(base, simdoffset, v0, v1);
}




template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1,
                    SimdDouble *          v2,
                    SimdDouble *          v3)
{
    assert(std::size_t(offset) % 64 == 0);
    gatherLoadBySimdIntTranspose<align>(base, simdLoadDI(offset), v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1)
{
    assert(std::size_t(offset) % 64 == 0);
    gatherLoadBySimdIntTranspose<align>(base, simdLoadDI(offset), v0, v1);
}

static const int c_simdBestPairAlignmentDouble = 2;

template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const double *        base,
                     const std::int32_t    offset[],
                     SimdDouble *          v0,
                     SimdDouble *          v1,
                     SimdDouble *          v2)
{
    SimdDInt32 simdoffset;

    assert(std::size_t(offset) % 32 == 0);

    simdoffset = simdLoadDI(offset);

    // All instructions might be latency ~4 on MIC, so we use shifts where we
    // only need a single instruction (since the shift parameter is an immediate),
    // but multiplication otherwise.
    if (align == 4)
    {
        simdoffset = simdoffset << 2;
    }
    else if (align == 8)
    {
        simdoffset = simdoffset << 3;
    }
    else
    {
        simdoffset = simdoffset * SimdDInt32(align);
    }

    v0->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base,   8);
    v1->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+1, 8);
    v2->simdInternal_ = _mm512_i32gather_pd(simdoffset.simdInternal_, base+2, 8);
}

template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(double *            base,
                       const std::int32_t  offset[],
                       SimdDouble          v0,
                       SimdDouble          v1,
                       SimdDouble          v2)
{
    SimdDInt32 simdoffset;

    assert(std::size_t(offset) % 32 == 0);

    simdoffset = simdLoadDI(offset);

    // All instructions might be latency ~4 on MIC, so we use shifts where we
    // only need a single instruction (since the shift parameter is an immediate),
    // but multiplication otherwise.
    if (align == 4)
    {
        simdoffset = simdoffset << 2;
    }
    else if (align == 8)
    {
        simdoffset = simdoffset << 3;
    }
    else
    {
        simdoffset = simdoffset * SimdDInt32(align);
    }

    _mm512_i32scatter_pd(base,   simdoffset.simdInternal_, v0.simdInternal_, 8);
    _mm512_i32scatter_pd(base+1, simdoffset.simdInternal_, v1.simdInternal_, 8);
    _mm512_i32scatter_pd(base+2, simdoffset.simdInternal_, v2.simdInternal_, 8);
}

template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(double *            base,
                      const std::int32_t  offset[],
                      SimdDouble          v0,
                      SimdDouble          v1,
                      SimdDouble          v2)
{
    __m512d t[4], t5, t6, t7, t8;
    GMX_ALIGNED(std::int64_t, 8)    o[8];
    _mm512_store_epi64(o, _mm512_cvtepi32_epi64(_mm256_mullo_epi32(_mm256_load_si256((const __m256i*)(offset  )), _mm256_set1_epi32(align))));
    t5   = _mm512_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t6   = _mm512_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t7   = _mm512_unpacklo_pd(v2.simdInternal_, _mm512_setzero_pd());
    t8   = _mm512_unpackhi_pd(v2.simdInternal_, _mm512_setzero_pd());
    t[0] = _mm512_mask_permutex_pd(t5, avx512Int2Mask(0xCC), t7, 0x4E);
    t[1] = _mm512_mask_permutex_pd(t6, avx512Int2Mask(0xCC), t8, 0x4E);
    t[2] = _mm512_mask_permutex_pd(t7, avx512Int2Mask(0x33), t5, 0x4E);
    t[3] = _mm512_mask_permutex_pd(t8, avx512Int2Mask(0x33), t6, 0x4E);
    if (align < 4)
    {
        for (int i = 0; i < 4; i++)
        {
            _mm512_mask_storeu_pd(base + o[0 + i], avx512Int2Mask(7), _mm512_castpd256_pd512(
                                          _mm256_add_pd(_mm256_loadu_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i]))));
            _mm512_mask_storeu_pd(base + o[4 + i], avx512Int2Mask(7), _mm512_castpd256_pd512(
                                          _mm256_add_pd(_mm256_loadu_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1))));
        }
    }
    else
    {
        if (align % 4 == 0)
        {
            for (int i = 0; i < 4; i++)
            {
                _mm256_store_pd(base + o[0 + i],
                                _mm256_add_pd(_mm256_load_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i])));
                _mm256_store_pd(base + o[4 + i],
                                _mm256_add_pd(_mm256_load_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1)));
            }
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                _mm256_storeu_pd(base + o[0 + i],
                                 _mm256_add_pd(_mm256_loadu_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i])));
                _mm256_storeu_pd(base + o[4 + i],
                                 _mm256_add_pd(_mm256_loadu_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1)));
            }
        }
    }
}

template <int align>
static inline void gmx_simdcall
transposeScatterDecrU(double *            base,
                      const std::int32_t  offset[],
                      SimdDouble          v0,
                      SimdDouble          v1,
                      SimdDouble          v2)
{
    __m512d t[4], t5, t6, t7, t8;
    GMX_ALIGNED(std::int64_t, 8)    o[8];
    _mm512_store_epi64(o, _mm512_cvtepi32_epi64(_mm256_mullo_epi32(_mm256_load_si256((const __m256i*)(offset  )), _mm256_set1_epi32(align))));
    t5   = _mm512_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t6   = _mm512_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t7   = _mm512_unpacklo_pd(v2.simdInternal_, _mm512_setzero_pd());
    t8   = _mm512_unpackhi_pd(v2.simdInternal_, _mm512_setzero_pd());
    t[0] = _mm512_mask_permutex_pd(t5, avx512Int2Mask(0xCC), t7, 0x4E);
    t[2] = _mm512_mask_permutex_pd(t7, avx512Int2Mask(0x33), t5, 0x4E);
    t[1] = _mm512_mask_permutex_pd(t6, avx512Int2Mask(0xCC), t8, 0x4E);
    t[3] = _mm512_mask_permutex_pd(t8, avx512Int2Mask(0x33), t6, 0x4E);
    if (align < 4)
    {
        for (int i = 0; i < 4; i++)
        {
            _mm512_mask_storeu_pd(base + o[0 + i], avx512Int2Mask(7), _mm512_castpd256_pd512(
                                          _mm256_sub_pd(_mm256_loadu_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i]))));
            _mm512_mask_storeu_pd(base + o[4 + i], avx512Int2Mask(7), _mm512_castpd256_pd512(
                                          _mm256_sub_pd(_mm256_loadu_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1))));
        }
    }
    else
    {
        if (align % 4 == 0)
        {
            for (int i = 0; i < 4; i++)
            {
                _mm256_store_pd(base + o[0 + i],
                                _mm256_sub_pd(_mm256_load_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i])));
                _mm256_store_pd(base + o[4 + i],
                                _mm256_sub_pd(_mm256_load_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1)));
            }
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                _mm256_storeu_pd(base + o[0 + i],
                                 _mm256_sub_pd(_mm256_loadu_pd(base + o[0 + i]), _mm512_castpd512_pd256(t[i])));
                _mm256_storeu_pd(base + o[4 + i],
                                 _mm256_sub_pd(_mm256_loadu_pd(base + o[4 + i]), _mm512_extractf64x4_pd(t[i], 1)));
            }
        }
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdDouble    scalar,
                        SimdDouble *  triplets0,
                        SimdDouble *  triplets1,
                        SimdDouble *  triplets2)
{
    triplets0->simdInternal_ = _mm512_castsi512_pd(_mm512_permutexvar_epi32(_mm512_set_epi32(5, 4, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0, 1, 0, 1, 0),
                                                                            _mm512_castpd_si512(scalar.simdInternal_)));
    triplets1->simdInternal_ = _mm512_castsi512_pd(_mm512_permutexvar_epi32(_mm512_set_epi32(11, 10, 9, 8, 9, 8, 9, 8, 7, 6, 7, 6, 7, 6, 5, 4),
                                                                            _mm512_castpd_si512(scalar.simdInternal_)));
    triplets2->simdInternal_ = _mm512_castsi512_pd(_mm512_permutexvar_epi32(_mm512_set_epi32(15, 14, 15, 14, 15, 14, 13, 12, 13, 12, 13, 12, 11, 10, 11, 10),
                                                                            _mm512_castpd_si512(scalar.simdInternal_)));
}


static inline double gmx_simdcall
reduceIncr4ReturnSum(double *    m,
                     SimdDouble  v0,
                     SimdDouble  v1,
                     SimdDouble  v2,
                     SimdDouble  v3)
{
    __m512d t0, t2;
    __m256d t3, t4;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm512_add_pd(v0.simdInternal_, _mm512_permute_pd(v0.simdInternal_, 0x55));
    t2 = _mm512_add_pd(v2.simdInternal_, _mm512_permute_pd(v2.simdInternal_, 0x55));
    t0 = _mm512_mask_add_pd(t0, avx512Int2Mask(0xAA), v1.simdInternal_, _mm512_permute_pd(v1.simdInternal_, 0x55));
    t2 = _mm512_mask_add_pd(t2, avx512Int2Mask(0xAA), v3.simdInternal_, _mm512_permute_pd(v3.simdInternal_, 0x55));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0x4E));
    t0 = _mm512_mask_add_pd(t0, avx512Int2Mask(0xF0), t2, _mm512_shuffle_f64x2(t2, t2, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_shuffle_f64x2(t0, t0, 0xB1));
    t0 = _mm512_mask_shuffle_f64x2(t0, avx512Int2Mask(0x0C), t0, t0, 0xEE);

    t3 = _mm512_castpd512_pd256(t0);
    t4 = _mm256_load_pd(m);
    t4 = _mm256_add_pd(t4, t3);
    _mm256_store_pd(m, t4);

    t3 = _mm256_add_pd(t3, _mm256_permutex_pd(t3, 0x4E));
    t3 = _mm256_add_pd(t3, _mm256_permutex_pd(t3, 0xB1));

    return _mm_cvtsd_f64(_mm256_castpd256_pd128(t3));
}

static inline SimdDouble gmx_simdcall
loadDualHsimd(const double * m0,
              const double * m1)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    return {
               _mm512_insertf64x4(_mm512_castpd256_pd512(_mm256_load_pd(m0)),
                                  _mm256_load_pd(m1), 1)
    };
}

static inline SimdDouble gmx_simdcall
loadDuplicateHsimd(const double * m)
{
    assert(std::size_t(m) % 32 == 0);

    return {
               _mm512_broadcast_f64x4(_mm256_load_pd(m))
    };
}

static inline SimdDouble gmx_simdcall
load1DualHsimd(const double * m)
{
    return {
               _mm512_insertf64x4(_mm512_broadcastsd_pd(_mm_loadu_pd(m)),      //note: this does not load 128bit each
                                  _mm256_broadcastsd_pd(_mm_loadu_pd(m+1)), 1) //because both intriniscs get combined
    };
}


static inline void gmx_simdcall
storeDualHsimd(double *     m0,
               double *     m1,
               SimdDouble   a)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    _mm256_store_pd(m0, _mm512_castpd512_pd256(a.simdInternal_));
    _mm256_store_pd(m1, _mm512_extractf64x4_pd(a.simdInternal_, 1));
}

static inline void gmx_simdcall
incrDualHsimd(double *     m0,
              double *     m1,
              SimdDouble   a)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    __m256d x;

    // Lower half
    x = _mm256_load_pd(m0);
    x = _mm256_add_pd(x, _mm512_castpd512_pd256(a.simdInternal_));
    _mm256_store_pd(m0, x);

    // Upper half
    x = _mm256_load_pd(m1);
    x = _mm256_add_pd(x, _mm512_extractf64x4_pd(a.simdInternal_, 1));
    _mm256_store_pd(m1, x);
}

static inline void gmx_simdcall
decrHsimd(double *    m,
          SimdDouble  a)
{
    __m256d t;

    assert(std::size_t(m) % 32 == 0);

    a.simdInternal_ = _mm512_add_pd(a.simdInternal_, _mm512_shuffle_f64x2(a.simdInternal_, a.simdInternal_, 0xEE));
    t               = _mm256_load_pd(m);
    t               = _mm256_sub_pd(t, _mm512_castpd512_pd256(a.simdInternal_));
    _mm256_store_pd(m, t);
}


template <int align>
static inline void gmx_simdcall
gatherLoadTransposeHsimd(const double *       base0,
                         const double *       base1,
                         const std::int32_t   offset[],
                         SimdDouble *         v0,
                         SimdDouble *         v1)
{
    __m128i  idx0, idx1;
    __m256i  idx;
    __m512d  tmp1, tmp2;

    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base0) % 16 == 0);
    assert(std::size_t(base1) % 16 == 0);
    assert(std::size_t(align) % 2  == 0);

    idx0 = _mm_load_si128(reinterpret_cast<const __m128i*>(offset));

    idx0 = _mm_mullo_epi32(idx0, _mm_set1_epi32(align));
    idx1 = _mm_add_epi32(idx0, _mm_set1_epi32(1));

    idx = _mm256_inserti128_si256(_mm256_castsi128_si256(idx0), idx1, 1);

    tmp1 = _mm512_i32gather_pd(idx, base0, 8);
    tmp2 = _mm512_i32gather_pd(idx, base1, 8);

    v0->simdInternal_ = _mm512_shuffle_f64x2(tmp1, tmp2, 0x44 );
    v1->simdInternal_ = _mm512_shuffle_f64x2(tmp1, tmp2, 0xEE );
}

static inline double gmx_simdcall
reduceIncr4ReturnSumHsimd(double *     m,
                          SimdDouble   v0,
                          SimdDouble   v1)
{
    __m512d  t0;
    __m256d  t2, t3;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm512_add_pd(v0.simdInternal_, _mm512_permutex_pd(v0.simdInternal_, 0x4E));
    t0 = _mm512_mask_add_pd(t0, avx512Int2Mask(0xCC), v1.simdInternal_, _mm512_permutex_pd(v1.simdInternal_, 0x4E));
    t0 = _mm512_add_pd(t0, _mm512_permutex_pd(t0, 0xB1));
    t0 = _mm512_mask_shuffle_f64x2(t0, avx512Int2Mask(0xAA), t0, t0, 0xEE);

    t2 = _mm512_castpd512_pd256(t0);
    t3 = _mm256_load_pd(m);
    t3 = _mm256_add_pd(t3, t2);
    _mm256_store_pd(m, t3);

    t2 = _mm256_add_pd(t2, _mm256_permutex_pd(t2, 0x4E));
    t2 = _mm256_add_pd(t2, _mm256_permutex_pd(t2, 0xB1));

    return _mm_cvtsd_f64(_mm256_castpd256_pd128(t2));
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_512_UTIL_DOUBLE_H
