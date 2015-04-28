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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

#include "impl_x86_avx_256_simd_float.h"




namespace gmx
{

struct SimdDouble
{
    __m256d r;
};

struct SimdDInt32
{
    __m128i i;
};

struct SimdDBool
{
    __m256d b;
};

struct SimdDIBool
{
    __m128i b;
};

static inline SimdDouble gmx_simdcall
simdLoadD(const double *m)
{
    assert(std::size_t(m) % 32 == 0);
    return {
               _mm256_load_pd(m)
    };
}

static inline SimdDouble gmx_simdcall
simdLoad1D(const double *m)
{
    return {
               _mm256_broadcast_sd(m)
    };
}

static inline SimdDouble gmx_simdcall
simdSet1D(double r)
{
    return {
               _mm256_set1_pd(r)
    };
}

static inline SimdDouble gmx_simdcall
simdSetZeroD()
{
    return {
               _mm256_setzero_pd()
    };
}

static inline void gmx_simdcall
simdStoreD(double *m, SimdDouble a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_pd(m, a.r);
}

static inline SimdDouble gmx_simdcall
simdLoadUD(const double *m)
{
    return {
               _mm256_loadu_pd(m)
    };
}

static inline void gmx_simdcall
simdStoreUD(double *m, SimdDouble a) { _mm256_storeu_pd(m, a.r); }

static inline SimdDInt32 gmx_simdcall
simdLoadDI(const std::int32_t * m)
{
    assert(std::size_t(m) % 16 == 0);
    return {
               _mm_load_si128(reinterpret_cast<const __m128i *>(m))
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
    assert(std::size_t(m) % 16 == 0);
    _mm_store_si128(reinterpret_cast<__m128i *>(m), a.i);
}

static inline SimdDInt32 gmx_simdcall
simdLoadUDI(const std::int32_t *m)
{
    return {
               _mm_loadu_si128(reinterpret_cast<const __m128i *>(m))
    };
}

static inline void gmx_simdcall
simdStoreUDI(std::int32_t * m, SimdDInt32 a)
{
    _mm_storeu_si128(reinterpret_cast<__m128i *>(m), a.i);
}

template<int index> gmx_simdcall
static inline std::int32_t
simdExtractDI(SimdDInt32 a)
{
    return _mm_extract_epi32(a.i, index);
}

static inline SimdDouble gmx_simdcall
simdAndD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_and_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdAndNotD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_andnot_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdOrD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_or_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdXorD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_xor_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdAddD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_add_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdSubD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_sub_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdMulD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_mul_pd(a.r, b.r)
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdDouble gmx_simdcall
simdFmaddD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_add_pd(_mm256_mul_pd(a.r, b.r), c.r)
    };
}

static inline SimdDouble gmx_simdcall
simdFmsubD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(_mm256_mul_pd(a.r, b.r), c.r)
    };
}

static inline SimdDouble gmx_simdcall
simdFnmaddD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(c.r, _mm256_mul_pd(a.r, b.r))
    };
}

static inline SimdDouble gmx_simdcall
simdFnmsubD(SimdDouble a, SimdDouble b, SimdDouble c)
{
    return {
               _mm256_sub_pd(_mm256_setzero_pd(), _mm256_add_pd(_mm256_mul_pd(a.r, b.r), c.r))
    };
}
#endif

static inline SimdDouble gmx_simdcall
simdRsqrtD(SimdDouble x)
{
    return {
               _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.r)))
    };
}

static inline SimdDouble gmx_simdcall
simdRcpD(SimdDouble x)
{
    return {
               _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x.r)))
    };
}

static inline SimdDouble gmx_simdcall
simdMulMaskD(SimdDouble a, SimdDouble b, SimdDBool m)
{
    return {
               _mm256_and_pd(_mm256_mul_pd(a.r, b.r), m.b)
    };
}

static inline SimdDouble
simdFmaddMaskD(SimdDouble a, SimdDouble b, SimdDouble c, SimdDBool m)
{
    return {
               _mm256_and_pd(_mm256_add_pd(_mm256_mul_pd(a.r, b.r), c.r), m.b)
    };
}

static inline SimdDouble
simdRsqrtMaskD(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.r = _mm256_blendv_pd(_mm256_set1_pd(1.0), x.r, m.b);
#endif
    return {
               _mm256_and_pd(_mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.r))), m.b)
    };
}

static inline SimdDouble
simdRcpMaskD(SimdDouble x, SimdDBool m)
{
#ifndef NDEBUG
    x.r = _mm256_blendv_pd(_mm256_set1_pd(1.0), x.r, m.b);
#endif
    return {
               _mm256_and_pd(_mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x.r))), m.b)
    };
}

static inline SimdDouble gmx_simdcall
simdAbsD(SimdDouble x)
{
    return {
               _mm256_andnot_pd( _mm256_set1_pd(GMX_DOUBLE_NEGZERO), x.r )
    };
}

static inline SimdDouble gmx_simdcall
simdNegD(SimdDouble x)
{
    return {
               _mm256_xor_pd(x.r, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
    };
}

static inline SimdDouble gmx_simdcall
simdMaxD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_max_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdMinD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_min_pd(a.r, b.r)
    };
}

static inline SimdDouble gmx_simdcall
simdRoundD(SimdDouble x)
{
    return {
               _mm256_round_pd(x.r, _MM_FROUND_NINT)
    };
}

static inline SimdDouble gmx_simdcall
simdTruncD(SimdDouble x)
{
    return {
               _mm256_round_pd(x.r, _MM_FROUND_TRUNC)
    };
}

static inline SimdDouble gmx_simdcall
simdFractionD(SimdDouble x)
{
    return {
               _mm256_sub_pd(x.r, _mm256_round_pd(x.r, _MM_FROUND_TRUNC))
    };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline SimdDouble gmx_simdcall
simdGetExponentD(SimdDouble x)
{
    const __m256d exponentMask      = _mm256_castsi256_pd( _mm256_set1_epi64x(0x7FF0000000000000LL));
    const __m128i exponentBias      = _mm_set1_epi32(1023);
    __m256i       iExponent;
    __m128i       iExponentLow, iExponentHigh;

    iExponent     = _mm256_castpd_si256(_mm256_and_pd(x.r, exponentMask));
    iExponentHigh = _mm256_extractf128_si256(iExponent, 0x1);
    iExponentLow  = _mm256_castsi256_si128(iExponent);
    iExponentLow  = _mm_srli_epi64(iExponentLow, 52);
    iExponentHigh = _mm_srli_epi64(iExponentHigh, 52);
    iExponentLow  = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(1, 1, 2, 0));
    iExponentHigh = _mm_shuffle_epi32(iExponentHigh, _MM_SHUFFLE(2, 0, 1, 1));
    iExponentLow  = _mm_or_si128(iExponentLow, iExponentHigh);
    iExponentLow  = _mm_sub_epi32(iExponentLow, exponentBias);
    return {
               _mm256_cvtepi32_pd(iExponentLow)
    };
}

static inline SimdDouble gmx_simdcall
simdSetExponentD(SimdDouble x)
{
    const __m128i exponentBias = _mm_set1_epi32(1023);
    __m128i       iExponentLow, iExponentHigh;

    iExponentLow  = _mm256_cvtpd_epi32(x.r);
    iExponentLow  = _mm_add_epi32(iExponentLow, exponentBias);
    iExponentHigh = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(3, 3, 2, 2));
    iExponentLow  = _mm_shuffle_epi32(iExponentLow, _MM_SHUFFLE(1, 1, 0, 0));
    iExponentHigh = _mm_slli_epi64(iExponentHigh, 52);
    iExponentLow  = _mm_slli_epi64(iExponentLow, 52);
    return {
               _mm256_castsi256_pd(_mm256_insertf128_si256(_mm256_castsi128_si256(iExponentLow), iExponentHigh, 0x1))
    };
}
#endif

static inline SimdDouble gmx_simdcall
simdGetMantissaD(SimdDouble x)
{
    const __m256d mantissaMask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL));
    const __m256d one          = _mm256_set1_pd(1.0);

    return {
               _mm256_or_pd(_mm256_and_pd(x.r, mantissaMask), one)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpEqD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.r, b.r, _CMP_EQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpNzD(SimdDouble a)
{
    return {
               _mm256_cmp_pd(a.r, _mm256_setzero_pd(), _CMP_NEQ_OQ)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpLtD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.r, b.r, _CMP_LT_OQ)
    };
}

static inline SimdDBool gmx_simdcall
simdCmpLeD(SimdDouble a, SimdDouble b)
{
    return {
               _mm256_cmp_pd(a.r, b.r, _CMP_LE_OQ)
    };
}

static inline SimdDBool gmx_simdcall
simdAndDB(SimdDBool a, SimdDBool b)
{
    return {
               _mm256_and_pd(a.b, b.b)
    };
}

static inline SimdDBool gmx_simdcall
simdOrDB(SimdDBool a, SimdDBool b)
{
    return {
               _mm256_or_pd(a.b, b.b)
    };
}

static inline bool gmx_simdcall
simdAnyTrueDB(SimdDBool a) { return _mm256_movemask_pd(a.b) != 0; }

static inline SimdDouble gmx_simdcall
simdMaskD(SimdDouble a, SimdDBool mask)
{
    return {
               _mm256_and_pd(a.r, mask.b)
    };
}

static inline SimdDouble gmx_simdcall
simdMaskNotD(SimdDouble a, SimdDBool mask)
{
    return {
               _mm256_andnot_pd(mask.b, a.r)
    };
}

static inline SimdDouble gmx_simdcall
simdBlendD(SimdDouble a, SimdDouble b, SimdDBool sel)
{
    return {
               _mm256_blendv_pd(a.r, b.r, sel.b)
    };
}

static inline double gmx_simdcall
simdReduceD(SimdDouble a)
{
    __m128d a0, a1;
    a.r = _mm256_add_pd(a.r, _mm256_permute_pd(a.r, 0b0101 ));
    a0  = _mm256_castpd256_pd128(a.r);
    a1  = _mm256_extractf128_pd(a.r, 0x1);
    a0  = _mm_add_sd(a0, a1);

    return *reinterpret_cast<double *>(&a0);
}

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

static inline SimdDInt32 gmx_simdcall
simdMulDI(SimdDInt32 a, SimdDInt32 b)
{
    return {
               _mm_mullo_epi32(a.i, b.i)
    };
}

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
simdAnyTrueDIB(SimdDIBool a) { return _mm_movemask_epi8(_mm_shuffle_epi32(a.b, _MM_SHUFFLE(1, 0, 1, 0))) != 0; }

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

static inline SimdDInt32 gmx_simdcall
simdBlendDI(SimdDInt32 a, SimdDInt32 b, SimdDIBool sel)
{
    return {
               _mm_blendv_epi8(a.i, b.i, sel.b)
    };
}

static inline SimdDInt32 gmx_simdcall
simdCvtD2I(SimdDouble a)
{
    return {
               _mm256_cvtpd_epi32(a.r)
    };
}

static inline SimdDInt32 gmx_simdcall
simdCvttD2I(SimdDouble a)
{
    return {
               _mm256_cvttpd_epi32(a.r)
    };
}

static inline SimdDouble gmx_simdcall
simdCvtI2D(SimdDInt32 a)
{
    return {
               _mm256_cvtepi32_pd(a.i)
    };
}

static inline SimdDIBool gmx_simdcall
simdCvtDB2DIB(SimdDBool a)
{
    __m128i a1 = _mm256_extractf128_si256(_mm256_castpd_si256(a.b), 0x1);
    __m128i a0 = _mm256_castsi256_si128(_mm256_castpd_si256(a.b));
    a0 = _mm_shuffle_epi32(a0, _MM_SHUFFLE(2, 0, 2, 0));
    a1 = _mm_shuffle_epi32(a1, _MM_SHUFFLE(2, 0, 2, 0));

    return {
               _mm_blend_epi16(a0, a1, 0xF0)
    };
}

static inline SimdDBool gmx_simdcall
simdCvtDIB2DB(SimdDIBool a)
{
    __m128d lo = _mm_castsi128_pd(_mm_unpacklo_epi32(a.b, a.b));
    __m128d hi = _mm_castsi128_pd(_mm_unpackhi_epi32(a.b, a.b));

    return {
               _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 0x1)
    };
}

static inline void gmx_simdcall
simdCvtF2DD(SimdFloat f, SimdDouble *d0, SimdDouble *d1)
{
    d0->r = _mm256_cvtps_pd(_mm256_castps256_ps128(f.r));
    d1->r = _mm256_cvtps_pd(_mm256_extractf128_ps(f.r, 0x1));
}

static inline SimdFloat gmx_simdcall
simdCvtDD2F(SimdDouble d0, SimdDouble d1)
{
    __m128 f0 = _mm256_cvtpd_ps(d0.r);
    __m128 f1 = _mm256_cvtpd_ps(d1.r);
    return {
               _mm256_insertf128_ps(_mm256_castps128_ps256(f0), f1, 0x1)
    };
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD_DOUBLE_H
