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

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <emmintrin.h>

#include "impl_x86_sse2_simd_float.h"

namespace gmx
{

struct SimdDouble
{
    __m128d r;
};

struct SimdDInt32
{
    __m128i i;
};

struct SimdDBool
{
    __m128d b;
};

struct SimdDIBool
{
    __m128i b;
};

static inline SimdDouble gmx_simdcall
simdLoadD(const double *m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_pd(m)
    };
}

static inline SimdDouble gmx_simdcall
simdLoad1D(const double *m)
{
    return {
               _mm_load1_pd(m)
    };
}

static inline SimdDouble gmx_simdcall
simdSet1D(double r)
{
    return {
               _mm_set1_pd(r)
    };
}

static inline SimdDouble gmx_simdcall
simdSetZeroD()
{
    return {
               _mm_setzero_pd()
    };
}

static inline void gmx_simdcall
simdStoreD(double *m, SimdDouble a)
{
    assert(std::size_t(m) % 16 == 0);
    _mm_store_pd(m, a.r);
}

static inline SimdDouble gmx_simdcall
simdLoadUD(const double *m)
{
    return {
               _mm_loadu_pd(m)
    };
}

static inline void gmx_simdcall
simdStoreUD(double *m, SimdDouble a) { _mm_storeu_pd(m, a.r); }

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
{
    assert(std::size_t(m) % 8 == 0);
    return {
               _mm_loadl_epi64(reinterpret_cast<const __m128i *>(m))
    };
}

static inline SimdDInt32 gmx_simdcall
simdSet1DI(std::int32_t b)
{
    return {
               _mm_set1_epi32(b)
    };
}

static inline SimdDInt32 gmx_simdcall
simdSetZeroDI()
{
    return {
               _mm_setzero_si128()
    };
}

static inline void gmx_simdcall
simdStoreDI(std::int32_t * m, SimdDInt32 a)
{
    assert(std::size_t(m) % 8 == 0);
    _mm_storel_epi64(reinterpret_cast<__m128i *>(m), a.i);
}

static inline SimdDInt32 gmx_simdcall
simdLoadUDI(const std::int32_t *m) { return simdLoadDI(m); }

static inline void gmx_simdcall
simdStoreUDI(std::int32_t * m, SimdDInt32 a) { simdStoreDI(m, a); }

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
template<int index> gmx_simdcall
static inline std::int32_t
simdExtractDI(SimdDInt32 a)
{
    return _mm_cvtsi128_si32( _mm_srli_si128(a.i, 4 * index) );
}
#endif

static inline SimdDouble gmx_simdcall
simdAndD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_and_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdAndNotD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_andnot_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdOrD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_or_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdXorD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_xor_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdAddD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_add_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdSubD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_sub_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdMulD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_mul_pd(a.r, b.r)
    };
}

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline SimdDouble gmx_simdcall
simdFmaddD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_add_pd(_mm_mul_pd(a.r, b.r), c.r)
    };
}

static inline SimdDouble gmx_simdcall
simdFmsubD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(_mm_mul_pd(a.r, b.r), c.r)
    };
}

static inline SimdDouble gmx_simdcall
simdFnmaddD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(c.r, _mm_mul_pd(a.r, b.r))
    };
}

static inline SimdDouble gmx_simdcall
simdFnmsubD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm_sub_pd(_mm_setzero_pd(), _mm_add_pd(_mm_mul_pd(a.r, b.r), c.r))
    };
}
#endif

static inline SimdDouble gmx_simdcall
simdRsqrtD(SimdDouble x)
{
    return {
               _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x.r)))
    };
}

static inline SimdDouble gmx_simdcall
simdRcpD(SimdDouble x)
{
    return {
               _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x.r)))
    };
}

static inline SimdDouble gmx_simdcall
simdMulMaskD(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm_and_pd(_mm_mul_pd(a.r, b.r), m.b)
    };
}

static inline SimdDouble gmx_simdcall
simdFmaddMaskD(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm_and_pd(_mm_add_pd(_mm_mul_pd(a.r, b.r), c.r), m.b)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
simdRsqrtMaskD(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.r = _mm_or_pd(_mm_andnot_pd(m.b, _mm_set1_pd(1.0)), _mm_and_pd(m.b, x.r));
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x.r))), m.b)
    };
}

static inline SimdDouble gmx_simdcall
simdRcpMaskD(SimdDouble x, SimdDBool m)
{
    // The result will always be correct since we mask the result with m, but
    // for debug builds we also want to make sure not to generate FP exceptions
#ifndef NDEBUG
    x.r = _mm_or_pd(_mm_andnot_pd(m.b, _mm_set1_pd(1.0)), _mm_and_pd(m.b, x.r));
#endif
    return {
               _mm_and_pd(_mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(x.r))), m.b)
    };
}
#endif

static inline SimdDouble gmx_simdcall
simdAbsD(SimdDouble x)
{
    return {
               _mm_andnot_pd( _mm_set1_pd(GMX_DOUBLE_NEGZERO), x.r )
    };
}

static inline SimdDouble gmx_simdcall
simdNegD(SimdDouble x)
{
    return {
               _mm_xor_pd(x.r, _mm_set1_pd(GMX_DOUBLE_NEGZERO))
    };
}

static inline SimdDouble gmx_simdcall
simdMaxD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_max_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdMinD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_min_pd(a.r, b.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
simdRoundD(SimdDouble x)
{
    return {
               _mm_cvtepi32_pd( _mm_cvtpd_epi32(x.r) )
    };
}

static inline SimdDouble gmx_simdcall
simdTruncD(SimdDouble x)
{
    return {
               _mm_cvtepi32_pd( _mm_cvttpd_epi32(x.r) )
    };
}

static inline SimdDouble gmx_simdcall
simdFractionD(SimdDouble x)
{
    return {
               _mm_sub_pd(x.r, _mm_cvtepi32_pd( _mm_cvttpd_epi32(x.r) ))
    };
}
#endif

static inline SimdDouble gmx_simdcall
simdGetExponentD(SimdDouble x)
{
    // Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds
    const __m128d exponentMask = _mm_castsi128_pd( _mm_set_epi32(0x7FF00000, 0x00000000, 0x7FF00000, 0x00000000) );
    const __m128i exponentBias = _mm_set1_epi32(1023);
    __m128i       iExponent;

    iExponent   = _mm_castpd_si128(_mm_and_pd(x.r, exponentMask));
    iExponent   = _mm_sub_epi32(_mm_srli_epi64(iExponent, 52), exponentBias);
    iExponent   = _mm_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0) );

    return {
               _mm_cvtepi32_pd(iExponent)
    };
}

static inline SimdDouble gmx_simdcall
simdGetMantissaD(SimdDouble x)
{
    // Don't use _mm_set1_epi64x() - on MSVC it is only supported for 64-bit builds
    const __m128d mantissaMask = _mm_castsi128_pd( _mm_set_epi32(0x000FFFFF, 0xFFFFFFFF, 0x000FFFFF, 0xFFFFFFFF) );
    const __m128d one          = _mm_set1_pd(1.0);

    return {
               _mm_or_pd(_mm_and_pd(x.r, mantissaMask), one)
    };
}

static inline SimdDouble gmx_simdcall
simdSetExponentD(SimdDouble x)
{
    const __m128i  exponentBias = _mm_set1_epi32(1023);
    __m128i        iExponent    = _mm_cvtpd_epi32(x.r);

    // After conversion integers will be in slot 0,1. Move them to 0,2 so
    // we can do a 64-bit shift and get them to the dp exponents.
    iExponent = _mm_shuffle_epi32(iExponent, _MM_SHUFFLE(3, 1, 2, 0));
    iExponent = _mm_slli_epi64(_mm_add_epi32(iExponent, exponentBias), 52);

    return {
               _mm_castsi128_pd(iExponent)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpEqD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmpeq_pd(a.r, b.r)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpNzD(SimdDouble a)
{
    return {
               _mm_cmpneq_pd(a.r, _mm_setzero_pd())
    };
}

static inline SimdDBool gmx_simdcall
simdCmpLtD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmplt_pd(a.r, b.r)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpLeD(SimdDouble a, SimdDouble b)
{
    return {
               _mm_cmple_pd(a.r, b.r)
    };
}

static inline SimdDBool gmx_simdcall
simdAndDB(SimdDBool a, SimdDBool b)
{
    return {
               _mm_and_pd(a.b, b.b)
    };
}

static inline SimdDBool gmx_simdcall
simdOrDB(SimdDBool a, SimdDBool b)
{
    return {
               _mm_or_pd(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueDB(SimdDBool a) { return _mm_movemask_pd(a.b) != 0; }

static inline SimdDouble gmx_simdcall
simdMaskD(SimdDouble a, SimdDBool mask)
{
    return {
               _mm_and_pd(a.r, mask.b)
    };
}

static inline SimdDouble gmx_simdcall
simdMaskNotD(SimdDouble a, SimdDBool mask)
{
    return {
               _mm_andnot_pd(mask.b, a.r)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDouble gmx_simdcall
simdBlendD(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm_or_pd(_mm_andnot_pd(sel.b, a.r), _mm_and_pd(sel.b, b.r))
    };
}
#endif

// Override for AVX-128-FMA and higher
#if GMX_SIMD_X86_SSE2 || GMX_SIMD_X86_SSE4_1
static inline double gmx_simdcall
simdReduceD(SimdDouble a)
{
    __m128d b = _mm_add_sd(a.r, _mm_shuffle_pd(a.r, a.r, _MM_SHUFFLE2(1, 1)));
    return *reinterpret_cast<double *>(&b);
}
#endif

static inline SimdDInt32 gmx_simdcall
simdSlliDI(SimdDInt32 a, int n)
{
    return {
               _mm_slli_epi32(a.i, n)
    };
}

static inline SimdDInt32 gmx_simdcall
simdSrliDI(SimdDInt32 a, int n)
{
    return {
               _mm_srli_epi32(a.i, n)
    };
}

static inline SimdDInt32 gmx_simdcall
simdAndDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_and_si128(a.i, b.i)
    };
}

static inline SimdDInt32 gmx_simdcall
simdAndNotDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_andnot_si128(a.i, b.i)
    };
}

static inline SimdDInt32 gmx_simdcall
simdOrDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_or_si128(a.i, b.i)
    };
}

static inline SimdDInt32 gmx_simdcall
simdXorDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_xor_si128(a.i, b.i)
    };
}

static inline SimdDInt32 gmx_simdcall
simdAddDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_add_epi32(a.i, b.i)
    };
}

static inline SimdDInt32 gmx_simdcall
simdSubDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_sub_epi32(a.i, b.i)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDInt32 gmx_simdcall
simdMulDI(SimdDInt32 a, SimdDInt32 b)
{

    __m128i tmpA = _mm_unpacklo_epi32(a.i, _mm_setzero_si128()); // 0 a[1] 0 a[0]
    __m128i tmpB = _mm_unpacklo_epi32(b.i, _mm_setzero_si128()); // 0 b[1] 0 b[0]

    __m128i tmpC  = _mm_mul_epu32(tmpA, tmpB);                   // 0 a[1]*b[1] 0 a[0]*b[0]

    return {
               _mm_shuffle_epi32(tmpC, _MM_SHUFFLE(3, 1, 2, 0))
    };
}
#endif

static inline SimdDIBool gmx_simdcall
simdCmpEqDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_cmpeq_epi32(a.i, b.i)
    };
}

static inline SimdDIBool gmx_simdcall
simdCmpLtDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_cmplt_epi32(a.i, b.i)
    };
}

static inline SimdDIBool gmx_simdcall
simdAndDIB(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm_and_si128(a.b, b.b)
    };
}

static inline SimdDIBool gmx_simdcall
simdOrDIB(SimdDIBool a, SimdDIBool b)
{
    return {
               _mm_or_si128(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueDIB(SimdDIBool a)
{
    return _mm_movemask_epi8(_mm_shuffle_epi32(a.b, _MM_SHUFFLE(1, 0, 1, 0))) != 0;
}

static inline SimdDInt32 gmx_simdcall
simdMaskDI(SimdDInt32 a, SimdDIBool mask)
{
    return {
               _mm_and_si128(a.i, mask.b)
    };
}

static inline SimdDInt32 gmx_simdcall
simdMaskNotDI(SimdDInt32 a, SimdDIBool mask)
{
    return {
               _mm_andnot_si128(mask.b, a.i)
    };
}

// Override for SSE4.1 and higher
#if GMX_SIMD_X86_SSE2
static inline SimdDInt32 gmx_simdcall
simdBlendDI(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm_or_si128(_mm_andnot_si128(sel.b, a.i), _mm_and_si128(sel.b, b.i))
    };
}
#endif

static inline SimdDInt32 gmx_simdcall
simdCvtD2I(SimdDouble a)
{
    return {
               _mm_cvtpd_epi32(a.r)
    };
}

static inline SimdDInt32 gmx_simdcall
simdCvttD2I(SimdDouble a)
{
    return {
               _mm_cvttpd_epi32(a.r)
    };
}

static inline SimdDouble gmx_simdcall
simdCvtI2D(SimdDInt32 a)
{
    return {
               _mm_cvtepi32_pd(a.i)
    };
}

static inline SimdDIBool gmx_simdcall
simdCvtDB2DIB(SimdDBool a)
{
    return {
               _mm_shuffle_epi32(_mm_castpd_si128(a.b), _MM_SHUFFLE(2, 0, 2, 0))
    };
}

static inline SimdDBool gmx_simdcall
simdCvtDIB2DB(SimdDIBool a)
{
    return {
               _mm_castsi128_pd(_mm_shuffle_epi32(a.b, _MM_SHUFFLE(1, 1, 0, 0)))
    };
}

static inline void gmx_simdcall
simdCvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
    d0->r = _mm_cvtps_pd(f.r);
    d1->r = _mm_cvtps_pd(_mm_movehl_ps(f.r, f.r));
}

static inline SimdFloat gmx_simdcall
simdCvtDD2F(SimdDouble d0, SimdDouble d1)
{
    return {
               _mm_movelh_ps(_mm_cvtpd_ps(d0.r), _mm_cvtpd_ps(d1.r))
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE2_SIMD_DOUBLE_H
