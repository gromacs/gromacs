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

#ifndef GMX_SIMD_IMPL_X86_MIC_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_MIC_UTIL_FLOAT_H

#include "config.h"

#include <cassert>
#include <cstdint>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_mic_simd_float.h"

namespace gmx
{

// On MIC it is better to use scatter operations, so we define the load routines
// that use a SIMD offset variable first.

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *  base,
                             SimdFInt32     simdoffset,
                             SimdFloat *    v0,
                             SimdFloat *    v1,
                             SimdFloat *    v2,
                             SimdFloat *    v3)
{
    assert(std::size_t(base) % 16 == 0);
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
        simdoffset = simdoffset * SimdFInt32(align);
    }

    v0->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base,   sizeof(float));
    v1->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+1, sizeof(float));
    v2->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+2, sizeof(float));
    v3->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+3, sizeof(float));
}

template <int align>
static inline void gmx_simdcall
gatherLoadUBySimdIntTranspose(const float *  base,
                              SimdFInt32     simdoffset,
                              SimdFloat *    v0,
                              SimdFloat *    v1)
{
    // All instructions might be latency ~4 on MIC, so we use shifts where we
    // only need a single instruction (since the shift parameter is an immediate),
    // but multiplication otherwise.
    // For align == 2 we can merge the constant into the scale parameter,
    // which can take constants up to 8 in total.
    if (align == 2)
    {
        v0->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base,   align * sizeof(float));
        v1->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+1, align * sizeof(float));
    }
    else
    {
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
            simdoffset = simdoffset * SimdFInt32(align);
        }
        v0->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base,   sizeof(float));
        v1->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+1, sizeof(float));
    }
}

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      simdoffset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    assert(std::size_t(base) % 8 == 0);
    assert(align % 2 == 0);
    gatherLoadUBySimdIntTranspose<align>(base, simdoffset, v0, v1);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1,
                    SimdFloat *          v2,
                    SimdFloat *          v3)
{
    gatherLoadBySimdIntTranspose<align>(base, simdLoad(offset, SimdFInt32Tag()), v0, v1, v2, v3);
}

template <int align>
static inline void gmx_simdcall
gatherLoadTranspose(const float *        base,
                    const std::int32_t   offset[],
                    SimdFloat *          v0,
                    SimdFloat *          v1)
{
    gatherLoadBySimdIntTranspose<align>(base, simdLoad(offset, SimdFInt32Tag()), v0, v1);
}

static const int c_simdBestPairAlignmentFloat = 2;

template <int align>
static inline void gmx_simdcall
gatherLoadUTranspose(const float *        base,
                     const std::int32_t   offset[],
                     SimdFloat *          v0,
                     SimdFloat *          v1,
                     SimdFloat *          v2)
{
    SimdFInt32 simdoffset;

    assert(std::size_t(offset) % 64 == 0);

    simdoffset = simdLoad(offset, SimdFInt32Tag());

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
        simdoffset = simdoffset * SimdFInt32(align);
    }

    v0->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base,   sizeof(float));
    v1->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+1, sizeof(float));
    v2->simdInternal_ = _mm512_i32gather_ps(simdoffset.simdInternal_, base+2, sizeof(float));
}


template <int align>
static inline void gmx_simdcall
transposeScatterStoreU(float *              base,
                       const std::int32_t   offset[],
                       SimdFloat            v0,
                       SimdFloat            v1,
                       SimdFloat            v2)
{
    SimdFInt32 simdoffset;

    assert(std::size_t(offset) % 64 == 0);

    simdoffset = simdLoad(offset, SimdFInt32Tag());

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
        simdoffset = simdoffset * SimdFInt32(align);
    }

    _mm512_i32scatter_ps(base,   simdoffset.simdInternal_, v0.simdInternal_, sizeof(float));
    _mm512_i32scatter_ps(base+1, simdoffset.simdInternal_, v1.simdInternal_, sizeof(float));
    _mm512_i32scatter_ps(base+2, simdoffset.simdInternal_, v2.simdInternal_, sizeof(float));
}


template <int align>
static inline void gmx_simdcall
transposeScatterIncrU(float *              base,
                      const std::int32_t   offset[],
                      SimdFloat            v0,
                      SimdFloat            v1,
                      SimdFloat            v2)
{
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata0[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata1[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata2[GMX_SIMD_FLOAT_WIDTH];

    store(rdata0, v0);
    store(rdata1, v1);
    store(rdata2, v2);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[ align * offset[i] + 0] += rdata0[i];
        base[ align * offset[i] + 1] += rdata1[i];
        base[ align * offset[i] + 2] += rdata2[i];
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
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata0[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata1[GMX_SIMD_FLOAT_WIDTH];
    GMX_ALIGNED(float, GMX_SIMD_FLOAT_WIDTH)  rdata2[GMX_SIMD_FLOAT_WIDTH];

    store(rdata0, v0);
    store(rdata1, v1);
    store(rdata2, v2);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[ align * offset[i] + 0] -= rdata0[i];
        base[ align * offset[i] + 1] -= rdata1[i];
        base[ align * offset[i] + 2] -= rdata2[i];
    }
}

static inline void gmx_simdcall
expandScalarsToTriplets(SimdFloat    scalar,
                        SimdFloat *  triplets0,
                        SimdFloat *  triplets1,
                        SimdFloat *  triplets2)
{
    triplets0->simdInternal_ = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0),
                                                                           _mm512_castps_si512(scalar.simdInternal_)));
    triplets1->simdInternal_ = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 6, 6, 6, 5, 5),
                                                                           _mm512_castps_si512(scalar.simdInternal_)));
    triplets2->simdInternal_ = _mm512_castsi512_ps(_mm512_permutevar_epi32(_mm512_set_epi32(15, 15, 15, 14, 14, 14, 13, 13, 13, 12, 12, 12, 11, 11, 11, 10),
                                                                           _mm512_castps_si512(scalar.simdInternal_)));
}


static inline float gmx_simdcall
reduceIncr4ReturnSum(float *    m,
                     SimdFloat  v0,
                     SimdFloat  v1,
                     SimdFloat  v2,
                     SimdFloat  v3)
{
    float  f;
    __m512 t0, t1, t2, t3;

    assert(std::size_t(m) % 16 == 0);

    t0 = _mm512_add_ps(v0.simdInternal_, _mm512_swizzle_ps(v0.simdInternal_, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v2.simdInternal_, _mm512_swizzle_ps(v2.simdInternal_, _MM_SWIZ_REG_BADC));
    t1 = _mm512_add_ps(v1.simdInternal_, _mm512_swizzle_ps(v1.simdInternal_, _MM_SWIZ_REG_BADC));
    t1 = _mm512_mask_add_ps(t1, _mm512_int2mask(0xCCCC), v3.simdInternal_, _mm512_swizzle_ps(v3.simdInternal_, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t2 = _mm512_mask_add_ps(t2, _mm512_int2mask(0xAAAA), t1, _mm512_swizzle_ps(t1, _MM_SWIZ_REG_CDAB));

    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_BADC));
    t2 = _mm512_add_ps(t2, _mm512_permute4f128_ps(t2, _MM_PERM_CDAB));

    t0 = _mm512_mask_extload_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE);
    t0 = _mm512_add_ps(t0, t2);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0xF), t0);

    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_BADC));
    t2 = _mm512_add_ps(t2, _mm512_swizzle_ps(t2, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t2);
    return f;
}

static inline SimdFloat gmx_simdcall
loadDualHsimd(const float * m0,
              const float * m1)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    return _mm512_castpd_ps(_mm512_mask_extload_pd(_mm512_extload_pd(reinterpret_cast<const double *>(m0), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE),
                                                   _mm512_int2mask(0xF0), reinterpret_cast<const double *>(m1), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static inline SimdFloat gmx_simdcall
loadDuplicateHsimd(const float * m)
{
    assert(std::size_t(m) % 32 == 0);

    return _mm512_castpd_ps(_mm512_extload_pd(reinterpret_cast<const double *>(m), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
}

static inline SimdFloat gmx_simdcall
load1DualHsimd(const float * m)
{
    return _mm512_mask_extload_ps(_mm512_extload_ps(m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE), _mm512_int2mask(0xFF00),
                                  m+1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
}


static inline void gmx_simdcall
storeDualHsimd(float *     m0,
               float *     m1,
               SimdFloat   a)
{
    __m512 t0;

    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    _mm512_mask_packstorelo_ps(m0, _mm512_int2mask(0x00FF), a.simdInternal_);
    _mm512_mask_packstorelo_ps(m1, _mm512_int2mask(0xFF00), a.simdInternal_);
}

static inline void gmx_simdcall
incrDualHsimd(float *     m0,
              float *     m1,
              SimdFloat   a)
{
    assert(std::size_t(m0) % 32 == 0);
    assert(std::size_t(m1) % 32 == 0);

    __m512 x;

    // Update lower half
    x = _mm512_castpd_ps(_mm512_extload_pd(reinterpret_cast<const double *>(m0), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
    x = _mm512_add_ps(x, a.simdInternal_);
    _mm512_mask_packstorelo_ps(m0, _mm512_int2mask(0x00FF), x);

    // Update upper half
    x = _mm512_castpd_ps(_mm512_extload_pd(reinterpret_cast<const double *>(m1), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
    x = _mm512_add_ps(x, a.simdInternal_);
    _mm512_mask_packstorelo_ps(m1, _mm512_int2mask(0xFF00), x);
}

static inline void gmx_simdcall
decrHsimd(float *    m,
          SimdFloat  a)
{
    __m512 t;

    assert(std::size_t(m) % 32 == 0);

    t = _mm512_castpd_ps(_mm512_extload_pd(reinterpret_cast<const double *>(m), _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE));
    a = _mm512_add_ps(a.simdInternal_, _mm512_permute4f128_ps(a.simdInternal_, _MM_PERM_BADC));
    t = _mm512_sub_ps(t, a.simdInternal_);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0x00FF), t);
}


template <int align>
static inline void gmx_simdcall
gatherLoadTransposeHsimd(const float *        base0,
                         const float *        base1,
                         const std::int32_t   offset[],
                         SimdFloat *          v0,
                         SimdFloat *          v1)
{
    __m512i idx0, idx1, idx;
    __m512  tmp1, tmp2;

    assert(std::size_t(offset) % 32 == 0);
    assert(std::size_t(base0) % 8 == 0);
    assert(std::size_t(base1) % 8 == 0);
    assert(std::size_t(align) % 2 == 0);

    idx0 = _mm512_loadunpacklo_epi32(_mm512_undefined_epi32(), offset);

    idx0 = _mm512_mullo_epi32(idx0, _mm512_set1_epi32(align));
    idx1 = _mm512_add_epi32(idx0, _mm512_set1_epi32(1));

    idx = _mm512_mask_permute4f128_epi32(idx0, _mm512_int2mask(0xFF00), idx1, _MM_PERM_BABA);

    tmp1 = _mm512_i32gather_ps(idx, base0, sizeof(float));
    tmp2 = _mm512_i32gather_ps(idx, base1, sizeof(float));

    v0->simdInternal_ = _mm512_mask_permute4f128_ps(tmp1, _mm512_int2mask(0xFF00), tmp2, _MM_PERM_BABA);
    v1->simdInternal_ = _mm512_mask_permute4f128_ps(tmp2, _mm512_int2mask(0x00FF), tmp1, _MM_PERM_DCDC);
}

static inline float gmx_simdcall
reduceIncr4ReturnSumHsimd(float *     m,
                          SimdFloat   v0,
                          SimdFloat   v1)
{
    float  f;
    __m512 t0, t1;

    assert(std::size_t(m) % 32 == 0);

    t0 = _mm512_add_ps(v0.simdInternal_, _mm512_swizzle_ps(v0.simdInternal_, _MM_SWIZ_REG_BADC));
    t0 = _mm512_mask_add_ps(t0, _mm512_int2mask(0xCCCC), v1.simdInternal_, _mm512_swizzle_ps(v1.simdInternal_, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));
    t0 = _mm512_add_ps(t0, _mm512_castpd_ps(_mm512_swizzle_pd(_mm512_castps_pd(t0), _MM_SWIZ_REG_BADC)));
    t0 = _mm512_mask_permute4f128_ps(t0, _mm512_int2mask(0xAAAA), t0, _MM_PERM_BADC);
    t1 = _mm512_mask_extload_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), m, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE);
    t1 = _mm512_add_ps(t1, t0);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0xF), t1);

    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_BADC));
    t0 = _mm512_add_ps(t0, _mm512_swizzle_ps(t0, _MM_SWIZ_REG_CDAB));

    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), t0);
    return f;
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_MIC_UTIL_FLOAT_H
