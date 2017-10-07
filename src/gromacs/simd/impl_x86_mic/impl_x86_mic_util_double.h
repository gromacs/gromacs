/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_MIC_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_MIC_UTIL_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_mic_simd_double.h"

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
    assert((size_t)base % 32 == 0);
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

    v0->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base,   sizeof(double));
    v1->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+1, sizeof(double));
    v2->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+2, sizeof(double));
    v3->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+3, sizeof(double));
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

    v0->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base,   sizeof(double));
    v1->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+1, sizeof(double));
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
    gatherLoadBySimdIntTranspose<align>(base, simdLoad(offset, SimdDInt32Tag()), v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const double *        base,
                    const std::int32_t    offset[],
                    SimdDouble *          v0,
                    SimdDouble *          v1)
{
    gatherLoadBySimdIntTranspose<align>(base, simdLoad(offset, SimdDInt32Tag()), v0, v1);
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

    simdoffset = simdLoad(offset, SimdDInt32Tag());

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

    v0->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base,   sizeof(double));
    v1->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+1, sizeof(double));
    v2->simdInternal_ = _mm512_i32logather_pd(simdoffset.simdInternal_, base+2, sizeof(double));
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

    simdoffset = simdLoad(offset, SimdDInt32Tag());

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

    _mm512_i32loscatter_pd(base,   simdoffset.simdInternal_, v0.simdInternal_, sizeof(double));
    _mm512_i32loscatter_pd(base+1, simdoffset.simdInternal_, v1.simdInternal_, sizeof(double));
    _mm512_i32loscatter_pd(base+2, simdoffset.simdInternal_, v2.simdInternal_, sizeof(double));
}

template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(double *            base,
                      const std::int32_t  offset[],
                      SimdDouble          v0,
                      SimdDouble          v1,
                      SimdDouble          v2)
{
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata0[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata1[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata2[GMX_SIMD_DOUBLE_WIDTH];

    store(rdata0, v0);
    store(rdata1, v1);
    store(rdata2, v2);

    for (int i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[ align * offset[i] + 0] += rdata0[i];
        base[ align * offset[i] + 1] += rdata1[i];
        base[ align * offset[i] + 2] += rdata2[i];
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
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata0[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata1[GMX_SIMD_DOUBLE_WIDTH];
    GMX_ALIGNED(double, GMX_SIMD_DOUBLE_WIDTH)  rdata2[GMX_SIMD_DOUBLE_WIDTH];

    store(rdata0, v0);
    store(rdata1, v1);
    store(rdata2, v2);

    for (int i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[ align * offset[i] + 0] -= rdata0[i];
        base[ align * offset[i] + 1] -= rdata1[i];
        base[ align * offset[i] + 2] -= rdata2[i];
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdDouble    scalar,
                        SimdDouble *  triplets0,
                        SimdDouble *  triplets1,
                        SimdDouble *  triplets2)
{
    triplets0->simdInternal_ = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0, 1, 0, 1, 0),
                                                                           _mm512_castpd_si512(scalar.simdInternal_)));
    triplets1->simdInternal_ = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(11, 10, 9, 8, 9, 8, 9, 8, 7, 6, 7, 6, 7, 6, 5, 4),
                                                                           _mm512_castpd_si512(scalar.simdInternal_)));
    triplets2->simdInternal_ = _mm512_castsi512_pd(_mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 15, 14, 15, 14, 13, 12, 13, 12, 13, 12, 11, 10, 11, 10),
                                                                           _mm512_castpd_si512(scalar.simdInternal_)));
}


static inline double gmx_simdcall
reduceIncr4ReturnSum(double *    m,
                     SimdDouble  v0,
                     SimdDouble  v1,
                     SimdDouble  v2,
                     SimdDouble  v3)
{
    double  d;
    __m512d t0, t1, t2, t3;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0x33), v0.simdInternal_, v2.simdInternal_), _MM_SWIZ_REG_BADC);
    t2 = _mm512_mask_blend_pd(_mm512_int2mask(0x33), v2.simdInternal_, v0.simdInternal_);
    t1 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0x33), v1.simdInternal_, v3.simdInternal_), _MM_SWIZ_REG_BADC);
    t3 = _mm512_mask_blend_pd(_mm512_int2mask(0x33), v3.simdInternal_, v1.simdInternal_);
    t0 = _mm512_add_pd(t0, t2);
    t1 = _mm512_add_pd(t1, t3);

    t2 = _mm512_swizzle_pd(_mm512_mask_blend_pd(_mm512_int2mask(0b01010101), t0, t1), _MM_SWIZ_REG_CDAB);
    t3 = _mm512_mask_blend_pd(_mm512_int2mask(0b01010101), t1, t0);
    t2 = _mm512_add_pd(t2, t3);

    t2 = _mm512_add_pd(t2, _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(t2), _MM_PERM_BADC)));

    t0 = _mm512_mask_extload_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    t0 = _mm512_add_pd(t0, t2);
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0xF), t0);

    t2 = _mm512_add_pd(t2, _mm512_swizzle_pd(t2, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_pd(t2, _mm512_swizzle_pd(t2, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_pd(&d, _mm512_mask2int(0x01), t2);
    return d;
}

static inline SimdDouble gmx_simdcall
loadDualHsimd(const double * m0,
              const double * m1)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    return _mm512_mask_extload_pd(_mm512_extload_pd(m0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE), _mm512_int2mask(0xF0),
                                  m1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
}

static inline SimdDouble gmx_simdcall
loadDuplicateHsimd(const double * m)
{
    assert(std::size_t(m) % 32 == 0);

    return _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
}

static inline SimdDouble gmx_simdcall
load1DualHsimd(const double * m)
{
    return _mm512_mask_extload_pd(_mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE), _mm512_int2mask(0xF0),
                                  m+1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
}


static inline void gmx_simdcall
storeDualHsimd(double *     m0,
               double *     m1,
               SimdDouble   a)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    _mm512_mask_packstorelo_pd(m0, _mm512_int2mask(0x0F), a.simdInternal_);
    _mm512_mask_packstorelo_pd(m1, _mm512_int2mask(0xF0), a.simdInternal_);
}

static inline void gmx_simdcall
incrDualHsimd(double *     m0,
              double *     m1,
              SimdDouble   a)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    __m512d x;

    // Update lower half
    x = _mm512_extload_pd(m0, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    x = _mm512_add_pd(x, a.simdInternal_);
    _mm512_mask_packstorelo_pd(m0, _mm512_int2mask(0x0F), x);

    // Update upper half
    x = _mm512_extload_pd(m1, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    x = _mm512_add_pd(x, a.simdInternal_);
    _mm512_mask_packstorelo_pd(m1, _mm512_int2mask(0xF0), x);
}

static inline void gmx_simdcall
decrHsimd(double *    m,
          SimdDouble  a)
{
    __m512d t;

    assert(std::size_t(m) % 32 == 0);

    t               = _mm512_extload_pd(m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    a.simdInternal_ = _mm512_add_pd(a.simdInternal_, _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(a.simdInternal_), _MM_PERM_BADC)));
    t               = _mm512_sub_pd(t, a.simdInternal_);
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0x0F), t);
}


template <int align>
static inline void gmx_simdcall
gatherLoadTransposeHsimd(const double *       base0,
                         const double *       base1,
                         const std::int32_t   offset[],
                         SimdDouble *         v0,
                         SimdDouble *         v1)
{
    __m512i  idx0, idx1, idx;
    __m512d  tmp1, tmp2;

    assert(std::size_t(offset) % 16 == 0);
    assert(std::size_t(base0) % 16 == 0);
    assert(std::size_t(base1) % 16 == 0);
    assert(std::size_t(align) % 2  == 0);

    idx0 = _mm512_extload_epi32(offset, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_permute4f128_epi32(idx0, _mm512_int2mask(0x00F0), idx1, _MM_PERM_AAAA);

    tmp1 = _mm512_i32logather_pd(idx, base0, sizeof(double));
    tmp2 = _mm512_i32logather_pd(idx, base1, sizeof(double));

    v0->simdInternal_ = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(tmp1), _mm512_int2mask(0xFF00), _mm512_castpd_ps(tmp2), _MM_PERM_BABA));
    v1->simdInternal_ = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(tmp2), _mm512_int2mask(0x00FF), _mm512_castpd_ps(tmp1), _MM_PERM_DCDC));
}

static inline double gmx_simdcall
reduceIncr4ReturnSumHsimd(double *     m,
                          SimdDouble   v0,
                          SimdDouble   v1)
{
    double   d;
    __m512d  t0, t1;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm512_add_pd(v0.simdInternal_, _mm512_swizzle_pd(v0.simdInternal_, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_pd(t0, _mm512_int2mask(0xCC), v1.simdInternal_, _mm512_swizzle_pd(v1.simdInternal_, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_CDAB));
    t0 = _mm512_castps_pd(_mm512_mask_permute4f128_ps(_mm512_castpd_ps(t0), _mm512_int2mask(0xCCCC),
                                                      _mm512_castpd_ps(t0), _MM_PERM_DCDC));

    t1 = _mm512_mask_extload_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE);
    t1 = _mm512_add_pd(t1, t0);
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_pd(t0, _mm512_swizzle_pd(t0, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_pd(&d, _mm512_mask2int(0x03), t0);
    return d;
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_MIC_UTIL_DOUBLE_H
