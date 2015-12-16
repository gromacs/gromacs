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

#ifndef GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H

#include "config.h"

#include <emmintrin.h>

#include "gromacs/utility/real.h"

#include "impl_x86_sse2_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble          __m128d
#define simdLoadD            _mm_load_pd
#define simdLoad1D           _mm_load1_pd
#define simdSet1D            _mm_set1_pd
#define simdStoreD           _mm_store_pd
#define simdLoadUD           _mm_loadu_pd
#define simdStoreUD          _mm_storeu_pd
#define simdSetZeroD         _mm_setzero_pd
#define simdAddD             _mm_add_pd
#define simdSubD             _mm_sub_pd
#define simdMulD             _mm_mul_pd
#define simdFmaddD(a, b, c)    _mm_add_pd(_mm_mul_pd(a, b), c)
#define simdFmsubD(a, b, c)    _mm_sub_pd(_mm_mul_pd(a, b), c)
#define simdFnmaddD(a, b, c)   _mm_sub_pd(c, _mm_mul_pd(a, b))
#define simdFnmsubD(a, b, c)   _mm_sub_pd(_mm_setzero_pd(), simdFmaddD(a, b, c))
#define simdAndD             _mm_and_pd
#define simdAndNotD          _mm_andnot_pd
#define simdOrD              _mm_or_pd
#define simdXorD             _mm_xor_pd
#define simdRsqrtD(x)        _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x)))
/* Don't use FMA for sqrt N-R iterations - this saves 1 instruction without FMA hardware */
#define simdRcpD(x)          _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x)))
#define simdAbsD(x)         _mm_andnot_pd(_mm_set1_pd(GMX_DOUBLE_NEGZERO), x)
#define simdNegD(x)         _mm_xor_pd(x, _mm_set1_pd(GMX_DOUBLE_NEGZERO))
#define simdMaxD             _mm_max_pd
#define simdMinD             _mm_min_pd
#define simdRoundD(x)        _mm_cvtepi32_pd(_mm_cvtpd_epi32(x))
#define simdTruncD(x)        _mm_cvtepi32_pd(_mm_cvttpd_epi32(x))
#define simdFractionD(x)     _mm_sub_pd(x, simdTruncD(x))
#define simdGetExponentD    simdGetExponentD_sse2
#define simdGetMantissaD    simdGetMantissaD_sse2
#define simdSetExponentD    simdSetExponentD_sse2
/* integer datatype corresponding to double: SimdDInt32 */
#define SimdDInt32          __m128i
#define simdLoadDI(m)        _mm_loadl_epi64((const __m128i *)(m))
#define simdSet1DI           _mm_set1_epi32
#define simdStoreDI(m, x)     _mm_storel_epi64((__m128i *)(m), x)
#define simdLoadUDI(m)       _mm_loadl_epi64((const __m128i *)(m))
#define simdStoreUDI(m, x)    _mm_storel_epi64((__m128i *)(m), x)
#define simdSetZeroDI        _mm_setzero_si128
#define simdCvtD2I           _mm_cvtpd_epi32
#define simdCvttD2I          _mm_cvttpd_epi32
#define simdCvtI2D           _mm_cvtepi32_pd
#define simdExtractDI(x, i)   _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (i)))
/* Integer logical ops on SimdDInt32 */
#define simdSlliDI           _mm_slli_epi32
#define simdSrliDI           _mm_srli_epi32
#define simdAndDI            _mm_and_si128
#define simdAndNotDI         _mm_andnot_si128
#define simdOrDI             _mm_or_si128
#define simdXorDI            _mm_xor_si128
/* Integer arithmetic ops on integer datatype corresponding to double */
#define simdAddDI            _mm_add_epi32
#define simdSubDI            _mm_sub_epi32
#define simdMulDI            simdMulDI_sse2
/* Boolean & comparison operations on SimdDouble */
#define SimdDBool            __m128d
#define simdCmpEqD            _mm_cmpeq_pd
#define simdCmpLtD            _mm_cmplt_pd
#define simdCmpLeD            _mm_cmple_pd
#define simdAndDB             _mm_and_pd
#define simdOrDB              _mm_or_pd
#define simdAnyTrueDB         _mm_movemask_pd
#define simdMaskD        _mm_and_pd
#define simdMaskNotD(a, sel) _mm_andnot_pd(sel, a)
#define simdBlendD(a, b, sel)  _mm_or_pd(_mm_andnot_pd(sel, a), _mm_and_pd(sel, b))
#define simdReduceD(a)        simdReduceD_sse2(a)

/* Boolean & comparison operations on SimdDInt32 */
#define SimdDIBool           __m128i
#define simdCmpEqDI           _mm_cmpeq_epi32
#define simdCmpLtDI           _mm_cmplt_epi32
#define simdAndDIB            _mm_and_si128
#define simdOrDIB             _mm_or_si128
#define simdAnyTrueDIB(x)     _mm_movemask_epi8(_mm_shuffle_epi32(x, _MM_SHUFFLE(1, 0, 1, 0)))
#define simdMaskDI       _mm_and_si128
#define simdMaskNotDI(a, sel)  _mm_andnot_si128(sel, a)
#define simdBlendDI(a, b, sel) _mm_or_si128(_mm_andnot_si128(sel, a), _mm_and_si128(sel, b))
#define simdCvtDB2DIB(x)      _mm_shuffle_epi32(_mm_castpd_si128(x), _MM_SHUFFLE(2, 0, 2, 0))
#define simdCvtDIB2DB(x)      _mm_castsi128_pd(_mm_shuffle_epi32(x, _MM_SHUFFLE(1, 1, 0, 0)))
/* Float/double conversion */
#define simdCvtF2DD(f, d0, d1)  { *d0 = _mm_cvtps_pd(f); *d1 = _mm_cvtps_pd(_mm_movehl_ps(f, f)); }
#define simdCvtDD2F(d0, d1)    _mm_movelh_ps(_mm_cvtpd_ps(d0), _mm_cvtpd_ps(d1))


/****************************************************
 * DOUBLE PRECISION IMPLEMENTATION HELPER FUNCTIONS *
 ****************************************************/
static inline __m128d gmx_simdcall
simdGetExponentD_sse2(__m128d x)
{
    /* Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds */
    const __m128d expmask      = _mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );
    const __m128i expbias      = _mm_set1_epi32(1023);
    __m128i       iexp;

    iexp   = _mm_castpd_si128(_mm_and_pd(x, expmask));
    iexp   = _mm_sub_epi32(_mm_srli_epi64(iexp, 52), expbias);
    iexp   = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(3, 1, 2, 0) );
    return _mm_cvtepi32_pd(iexp);
}

static inline __m128d gmx_simdcall
simdGetMantissaD_sse2(__m128d x)
{
    /* Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds */
    const __m128d mantmask = _mm_castsi128_pd( _mm_set_epi32(0x000FFFFF, 0xFFFFFFFF, 0x000FFFFF, 0xFFFFFFFF) );
    const __m128d one      = _mm_set1_pd(1.0);

    x = _mm_and_pd(x, mantmask);
    return _mm_or_pd(x, one);
}

static inline __m128d gmx_simdcall
simdSetExponentD_sse2(__m128d x)
{
    const __m128i  expbias      = _mm_set1_epi32(1023);
    __m128i        iexp         = _mm_cvtpd_epi32(x);

    /* After conversion integers will be in slot 0,1. Move them to 0,2 so
     * we can do a 64-bit shift and get them to the dp exponents. */
    iexp = _mm_shuffle_epi32(iexp, _MM_SHUFFLE(3, 1, 2, 0));
    iexp = _mm_slli_epi64(_mm_add_epi32(iexp, expbias), 52);
    return _mm_castsi128_pd(iexp);
}

static inline __m128i gmx_simdcall
simdMulDI_sse2(__m128i a, __m128i b)
{
    __m128i c;

    a = _mm_unpacklo_epi32(a, _mm_setzero_si128());       /* 0 a[1] 0 a[0] */
    b = _mm_unpacklo_epi32(b, _mm_setzero_si128());       /* 0 b[1] 0 b[0] */

    c  = _mm_mul_epu32(a, b);                             /* 0 a[1]*b[1] 0 a[0]*b[0] */
    return _mm_shuffle_epi32(c, _MM_SHUFFLE(3, 1, 2, 0)); /* 0 0 a[1]*b[1] a[0]*b[0] */
}

static inline double gmx_simdcall
simdReduceD_sse2(__m128d a)
{
    __m128d b;
    double  f;

    b = _mm_add_sd(a, _mm_shuffle_pd(a, a, _MM_SHUFFLE2(1, 1)));
    _mm_store_sd(&f, b);
    return f;
}

#endif /* GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H */
