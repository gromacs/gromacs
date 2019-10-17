/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2019, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_MIC_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_X86_MIC_SIMD4_FLOAT_H

#include "config.h"

#include <cassert>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_mic_simd_float.h"

namespace gmx
{

class Simd4Float
{
public:
    Simd4Float() {}

    Simd4Float(float f) : simdInternal_(_mm512_set1_ps(f)) {}

    // Internal utility constructor to simplify return statements
    Simd4Float(__m512 simd) : simdInternal_(simd) {}

    __m512 simdInternal_;
};

class Simd4FBool
{
public:
    Simd4FBool() {}

    // Internal utility constructor to simplify return statements
    Simd4FBool(__mmask16 simd) : simdInternal_(simd) {}

    __mmask16 simdInternal_;
};

static inline Simd4Float gmx_simdcall load4(const float* m)
{
    assert(size_t(m) % 16 == 0);
    return { _mm512_mask_extload_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), m,
                                    _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE) };
}

static inline void gmx_simdcall store4(float* m, Simd4Float a)
{
    assert(size_t(m) % 16 == 0);
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0xF), a.simdInternal_);
}

static inline Simd4Float gmx_simdcall load4U(const float* m)
{
    return { _mm512_mask_loadunpackhi_ps(
            _mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), m),
            _mm512_int2mask(0xF), m + 16) };
}

static inline void gmx_simdcall store4U(float* m, Simd4Float a)
{
    _mm512_mask_packstorelo_ps(m, _mm512_int2mask(0xF), a.simdInternal_);
    _mm512_mask_packstorehi_ps(m + 16, _mm512_int2mask(0xF), a.simdInternal_);
}

static inline Simd4Float gmx_simdcall simd4SetZeroF()
{
    return { _mm512_setzero_ps() };
}

static inline Simd4Float gmx_simdcall operator&(Simd4Float a, Simd4Float b)
{
    return { _mm512_castsi512_ps(_mm512_mask_and_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0xF),
                                                       _mm512_castps_si512(a.simdInternal_),
                                                       _mm512_castps_si512(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall andNot(Simd4Float a, Simd4Float b)
{
    return { _mm512_castsi512_ps(_mm512_mask_andnot_epi32(
            _mm512_undefined_epi32(), _mm512_int2mask(0xF), _mm512_castps_si512(a.simdInternal_),
            _mm512_castps_si512(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator|(Simd4Float a, Simd4Float b)
{
    return { _mm512_castsi512_ps(_mm512_mask_or_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0xF),
                                                      _mm512_castps_si512(a.simdInternal_),
                                                      _mm512_castps_si512(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator^(Simd4Float a, Simd4Float b)
{
    return { _mm512_castsi512_ps(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0xF),
                                                       _mm512_castps_si512(a.simdInternal_),
                                                       _mm512_castps_si512(b.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall operator+(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_add_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), a.simdInternal_,
                                b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_sub_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), a.simdInternal_,
                                b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall operator-(Simd4Float x)
{
    return { _mm512_mask_addn_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), x.simdInternal_,
                                 _mm512_setzero_ps()) };
}

static inline Simd4Float gmx_simdcall operator*(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_mul_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), a.simdInternal_,
                                b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { _mm512_mask_fmadd_ps(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { _mm512_mask_fmsub_ps(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnma(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { _mm512_mask_fnmadd_ps(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall fnms(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return { _mm512_mask_fnmsub_ps(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_) };
}

static inline Simd4Float gmx_simdcall rsqrt(Simd4Float x)
{
    return { _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), x.simdInternal_) };
}

static inline Simd4Float gmx_simdcall abs(Simd4Float x)
{
    return { _mm512_castsi512_ps(_mm512_mask_andnot_epi32(
            _mm512_undefined_epi32(), _mm512_int2mask(0xF),
            _mm512_castps_si512(_mm512_set1_ps(GMX_FLOAT_NEGZERO)), _mm512_castps_si512(x.simdInternal_))) };
}

static inline Simd4Float gmx_simdcall max(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_gmax_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), a.simdInternal_,
                                 b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall min(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_gmin_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), a.simdInternal_,
                                 b.simdInternal_) };
}

static inline Simd4Float gmx_simdcall round(Simd4Float x)
{
    return { _mm512_mask_round_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), x.simdInternal_,
                                  _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE) };
}

static inline Simd4Float gmx_simdcall trunc(Simd4Float x)
{
    return { _mm512_mask_round_ps(_mm512_undefined_ps(), _mm512_int2mask(0xF), x.simdInternal_,
                                  _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE) };
}

static inline float gmx_simdcall dotProduct(Simd4Float a, Simd4Float b)
{
    __m512 x = _mm512_mask_mul_ps(_mm512_setzero_ps(), _mm512_int2mask(0x7), a.simdInternal_,
                                  b.simdInternal_);
    x        = _mm512_add_ps(x, _mm512_swizzle_ps(x, _MM_SWIZ_REG_BADC));
    x        = _mm512_add_ps(x, _mm512_swizzle_ps(x, _MM_SWIZ_REG_CDAB));
    float f;
    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), x);
    return f;
}

static inline void gmx_simdcall transpose(Simd4Float* v0, Simd4Float* v1, Simd4Float* v2, Simd4Float* v3)
{
    v0->simdInternal_ = _mm512_mask_permute4f128_ps(v0->simdInternal_, _mm512_int2mask(0x00F0),
                                                    v1->simdInternal_, _MM_PERM_AAAA);
    v2->simdInternal_ = _mm512_mask_permute4f128_ps(v2->simdInternal_, _mm512_int2mask(0x00F0),
                                                    v3->simdInternal_, _MM_PERM_AAAA);
    v0->simdInternal_ = _mm512_mask_permute4f128_ps(v0->simdInternal_, _mm512_int2mask(0xFF00),
                                                    v2->simdInternal_, _MM_PERM_BABA);
    v0->simdInternal_ = _mm512_castsi512_ps(_mm512_permutevar_epi32(
            _mm512_set_epi32(15, 11, 7, 3, 14, 10, 6, 2, 13, 9, 5, 1, 12, 8, 4, 0),
            _mm512_castps_si512(v0->simdInternal_)));
    v1->simdInternal_ = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F),
                                                    v0->simdInternal_, _MM_PERM_BBBB);
    v2->simdInternal_ = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F),
                                                    v0->simdInternal_, _MM_PERM_CCCC);
    v3->simdInternal_ = _mm512_mask_permute4f128_ps(_mm512_setzero_ps(), _mm512_int2mask(0x000F),
                                                    v0->simdInternal_, _MM_PERM_DDDD);
}

// Picky, picky, picky:
// icc-16 complains about "Illegal value of immediate argument to intrinsic"
// unless we use
// 1) Ordered-quiet for ==
// 2) Unordered-quiet for !=
// 3) Ordered-signaling for < and <=

static inline Simd4FBool gmx_simdcall operator==(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ) };
}

static inline Simd4FBool gmx_simdcall operator!=(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_NEQ_UQ) };
}

static inline Simd4FBool gmx_simdcall operator<(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_LT_OS) };
}

static inline Simd4FBool gmx_simdcall operator<=(Simd4Float a, Simd4Float b)
{
    return { _mm512_mask_cmp_ps_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_LE_OS) };
}

static inline Simd4FBool gmx_simdcall operator&&(Simd4FBool a, Simd4FBool b)
{
    return { _mm512_kand(a.simdInternal_, b.simdInternal_) };
}

static inline Simd4FBool gmx_simdcall operator||(Simd4FBool a, Simd4FBool b)
{
    return { _mm512_kor(a.simdInternal_, b.simdInternal_) };
}

static inline bool gmx_simdcall anyTrue(Simd4FBool a)
{
    return (_mm512_mask2int(a.simdInternal_) & 0xF) != 0;
}

static inline Simd4Float gmx_simdcall selectByMask(Simd4Float a, Simd4FBool m)
{
    return { _mm512_mask_mov_ps(_mm512_setzero_ps(), m.simdInternal_, a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall selectByNotMask(Simd4Float a, Simd4FBool m)
{
    return { _mm512_mask_mov_ps(_mm512_setzero_ps(), _mm512_knot(m.simdInternal_), a.simdInternal_) };
}

static inline Simd4Float gmx_simdcall blend(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    return { _mm512_mask_blend_ps(sel.simdInternal_, a.simdInternal_, b.simdInternal_) };
}

static inline float gmx_simdcall reduce(Simd4Float a)
{
    __m512 x = a.simdInternal_;
    x        = _mm512_add_ps(x, _mm512_swizzle_ps(x, _MM_SWIZ_REG_BADC));
    x        = _mm512_add_ps(x, _mm512_swizzle_ps(x, _MM_SWIZ_REG_CDAB));
    float f;
    _mm512_mask_packstorelo_ps(&f, _mm512_mask2int(0x1), x);
    return f;
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_MIC_SIMD4_FLOAT_H
