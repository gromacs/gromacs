/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_X86_MIC_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_MIC_SIMD4_DOUBLE_H

#include "config.h"

#include <cassert>

#include <immintrin.h>

#include "gromacs/utility/basedefinitions.h"

#include "impl_x86_mic_simd_double.h"

namespace gmx
{

class Simd4Double
{
    public:
        Simd4Double() {}

        Simd4Double(double d) : simdInternal_(_mm512_set1_pd(d)) {}

        // Internal utility constructor to simplify return statements
        Simd4Double(__m512d simd) : simdInternal_(simd) {}

        __m512d  simdInternal_;
};

class Simd4DBool
{
    public:
        Simd4DBool() {}

        // Internal utility constructor to simplify return statements
        Simd4DBool(__mmask16 simd) : simdInternal_(simd) {}

        __mmask16  simdInternal_;
};

static inline Simd4Double gmx_simdcall
load4(const double *m)
{
    assert(size_t(m) % 32 == 0);
    return {
               _mm512_mask_extload_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), m, _MM_UPCONV_PD_NONE, _MM_BROADCAST_4X8, _MM_HINT_NONE)
    };
}

static inline void gmx_simdcall
store4(double *m, Simd4Double a)
{
    assert(size_t(m) % 32 == 0);
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0xF), a.simdInternal_);
}

static inline Simd4Double gmx_simdcall
load4U(const double *m)
{
    return {
               _mm512_mask_loadunpackhi_pd(_mm512_mask_loadunpacklo_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), m), _mm512_int2mask(0xF), m+8)
    };
}

static inline void gmx_simdcall
store4U(double *m, Simd4Double a)
{
    _mm512_mask_packstorelo_pd(m, _mm512_int2mask(0xF), a.simdInternal_);
    _mm512_mask_packstorehi_pd(m+8, _mm512_int2mask(0xF), a.simdInternal_);
}

static inline Simd4Double gmx_simdcall
simd4SetZeroD()
{
    return {
               _mm512_setzero_pd()
    };
}

static inline Simd4Double gmx_simdcall
operator&(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_castsi512_pd(_mm512_mask_and_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF), _mm512_castpd_si512(a.simdInternal_),
                                                         _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
andNot(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_castsi512_pd(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF), _mm512_castpd_si512(a.simdInternal_),
                                                            _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
operator|(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_castsi512_pd(_mm512_mask_or_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF), _mm512_castpd_si512(a.simdInternal_),
                                                        _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
operator^(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_castsi512_pd(_mm512_mask_xor_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF), _mm512_castpd_si512(a.simdInternal_),
                                                         _mm512_castpd_si512(b.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
operator+(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_add_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator-(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_sub_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator-(Simd4Double x)
{
    return {
               _mm512_mask_addn_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), x.simdInternal_, _mm512_setzero_pd())
    };
}

static inline Simd4Double gmx_simdcall
operator*(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_mul_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm512_mask_fmadd_pd(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm512_mask_fmsub_pd(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fnma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm512_mask_fnmadd_pd(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fnms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm512_mask_fnmsub_pd(a.simdInternal_, _mm512_int2mask(0xF), b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
rsqrt(Simd4Double x)
{
    return {
               _mm512_mask_cvtpslo_pd(_mm512_undefined_pd(),
                                      _mm512_int2mask(0xF),
                                      _mm512_mask_rsqrt23_ps(_mm512_undefined_ps(),
                                                             _mm512_int2mask(0xF),
                                                             _mm512_mask_cvtpd_pslo(_mm512_undefined_ps(),
                                                                                    _mm512_int2mask(0xF), x.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
abs(Simd4Double x)
{
    return {
               _mm512_castsi512_pd(_mm512_mask_andnot_epi32(_mm512_undefined_epi32(), _mm512_int2mask(0x00FF),
                                                            _mm512_castpd_si512(_mm512_set1_pd(GMX_DOUBLE_NEGZERO)),
                                                            _mm512_castpd_si512(x.simdInternal_)))

    };
}

static inline Simd4Double gmx_simdcall
max(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_gmax_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
min(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_gmin_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
round(Simd4Double x)
{
    return {
               _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), x.simdInternal_, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE)
    };
}

static inline Simd4Double gmx_simdcall
trunc(Simd4Double x)
{
    return {
               _mm512_mask_roundfxpnt_adjust_pd(_mm512_undefined_pd(), _mm512_int2mask(0xF), x.simdInternal_, _MM_FROUND_TO_ZERO, _MM_EXPADJ_NONE)
    };
}

static inline float gmx_simdcall
dotProduct(Simd4Double a, Simd4Double b)
{
    return _mm512_mask_reduce_add_pd(_mm512_int2mask(7),
                                     _mm512_mask_mul_pd(_mm512_undefined_pd(), _mm512_int2mask(7),
                                                        a.simdInternal_, b.simdInternal_));
}

static inline void gmx_simdcall
transpose(Simd4Double * v0, Simd4Double * v1,
          Simd4Double * v2, Simd4Double * v3)
{
    __m512i t0 = _mm512_mask_permute4f128_epi32(_mm512_castpd_si512(v0->simdInternal_), 0xFF00,
                                                _mm512_castpd_si512(v1->simdInternal_), _MM_PERM_BABA);
    __m512i t1 = _mm512_mask_permute4f128_epi32(_mm512_castpd_si512(v2->simdInternal_), 0xFF00,
                                                _mm512_castpd_si512(v3->simdInternal_), _MM_PERM_BABA);

    t0 = _mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 7, 6, 13, 12, 5, 4, 11, 10, 3, 2, 9, 8, 1, 0), t0);
    t1 = _mm512_permutevar_epi32(_mm512_set_epi32(15, 14, 7, 6, 13, 12, 5, 4, 11, 10, 3, 2, 9, 8, 1, 0), t1);

    v0->simdInternal_ = _mm512_mask_swizzle_pd(_mm512_castsi512_pd(t0), _mm512_int2mask(0xCC),
                                               _mm512_castsi512_pd(t1), _MM_SWIZ_REG_BADC);
    v1->simdInternal_ = _mm512_mask_swizzle_pd(_mm512_castsi512_pd(t1), _mm512_int2mask(0x33),
                                               _mm512_castsi512_pd(t0), _MM_SWIZ_REG_BADC);

    v2->simdInternal_ = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(v0->simdInternal_), _MM_PERM_DCDC));
    v3->simdInternal_ = _mm512_castps_pd(_mm512_permute4f128_ps(_mm512_castpd_ps(v1->simdInternal_), _MM_PERM_DCDC));
}

// Picky, picky, picky:
// icc-16 complains about "Illegal value of immediate argument to intrinsic"
// unless we use
// 1) Ordered-quiet for ==
// 2) Unordered-quiet for !=
// 3) Ordered-signaling for < and <=

static inline Simd4DBool gmx_simdcall
operator==(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator!=(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_NEQ_UQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator<(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_LT_OS)
    };
}

static inline Simd4DBool gmx_simdcall
operator<=(Simd4Double a, Simd4Double b)
{
    return {
               _mm512_mask_cmp_pd_mask(_mm512_int2mask(0xF), a.simdInternal_, b.simdInternal_, _CMP_LE_OS)
    };
}

static inline Simd4DBool gmx_simdcall
operator&&(Simd4DBool a, Simd4DBool b)
{
    return {
               _mm512_kand(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4DBool gmx_simdcall
operator||(Simd4DBool a, Simd4DBool b)
{
    return {
               _mm512_kor(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(Simd4DBool a)
{
    return (_mm512_mask2int(a.simdInternal_) & 0xF) != 0;
}

static inline Simd4Double gmx_simdcall
selectByMask(Simd4Double a, Simd4DBool m)
{
    return {
               _mm512_mask_mov_pd(_mm512_setzero_pd(), m.simdInternal_, a.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
selectByNotMask(Simd4Double a, Simd4DBool m)
{
    return {
               _mm512_mask_mov_pd(_mm512_setzero_pd(), _mm512_knot(m.simdInternal_), a.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
blend(Simd4Double a, Simd4Double b, Simd4DBool sel)
{
    return {
               _mm512_mask_blend_pd(sel.simdInternal_, a.simdInternal_, b.simdInternal_)
    };
}

static inline float gmx_simdcall
reduce(Simd4Double a)
{
    return _mm512_mask_reduce_add_pd(_mm512_int2mask(0xF), a.simdInternal_);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_MIC_SIMD4_DOUBLE_H
