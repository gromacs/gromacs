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

#ifndef GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <immintrin.h>

namespace gmx
{
    
struct Simd4Double
{
    __m256d r;
};

struct Simd4DBool
{
    __m256d b;
};

static inline Simd4Double gmx_simdcall
simd4LoadD(const double *m)
{
    assert(std::size_t(m) % 32 == 0);
    return { _mm256_load_pd(m) };
}

static inline Simd4Double gmx_simdcall
simd4Load1D(const double *m)
{
    return { _mm256_broadcast_sd(m) };
}

static inline Simd4Double gmx_simdcall
simd4Set1D(double r)
{
    return { _mm256_set1_pd(r) };
}

static inline Simd4Double gmx_simdcall
simd4SetZeroD()
{
    return { _mm256_setzero_pd() };
}

static inline void gmx_simdcall
simd4StoreD(double *m, Simd4Double a)
{
    assert(std::size_t(m) % 32 == 0);
    _mm256_store_pd(m, a.r);
}

static inline Simd4Double gmx_simdcall
simd4LoadUD(const double *m)
{
    return {_mm256_loadu_pd(m) };
}

static inline void gmx_simdcall
simd4StoreUD(double *m, Simd4Double a) { _mm256_storeu_pd(m, a.r); }

static inline Simd4Double gmx_simdcall
simd4AndD(Simd4Double a, Simd4Double b)
{
    return { _mm256_and_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4AndNotD(Simd4Double a, Simd4Double b)
{
    return { _mm256_andnot_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4OrD(Simd4Double a, Simd4Double b)
{
    return { _mm256_or_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4XorD(Simd4Double a, Simd4Double b)
{
    return { _mm256_xor_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4AddD(Simd4Double a, Simd4Double b)
{
    return { _mm256_add_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4SubD(Simd4Double a, Simd4Double b)
{
    return { _mm256_sub_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4MulD(Simd4Double a, Simd4Double b)
{
    return { _mm256_mul_pd(a.r, b.r) };
}

// Override for AVX2 and higher
#if GMX_SIMD_X86_AVX_256
static inline Simd4Double gmx_simdcall
simd4FmaddD(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return { _mm256_add_pd(_mm256_mul_pd(a.r, b.r), c.r) };
}

static inline Simd4Double gmx_simdcall
simd4FmsubD(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return { _mm256_sub_pd(_mm256_mul_pd(a.r, b.r), c.r) };
}

static inline Simd4Double gmx_simdcall
simd4FnmaddD(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return { _mm256_sub_pd(c.r, _mm256_mul_pd(a.r, b.r)) };
}

static inline Simd4Double gmx_simdcall
simd4FnmsubD(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return { _mm256_sub_pd(_mm256_setzero_pd(), _mm256_add_pd(_mm256_mul_pd(a.r, b.r), c.r)) };
}
#endif

static inline Simd4Double gmx_simdcall
simd4RsqrtD(Simd4Double x)
{
    return { _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.r))) };
}

static inline Simd4Double gmx_simdcall
simd4AbsD(Simd4Double x)
{
    return { _mm256_andnot_pd( _mm256_set1_pd(GMX_DOUBLE_NEGZERO), x.r ) };
}

static inline Simd4Double gmx_simdcall
simd4NegD(Simd4Double x)
{
    return { _mm256_xor_pd(x.r, _mm256_set1_pd(GMX_DOUBLE_NEGZERO)) };
}

static inline Simd4Double gmx_simdcall
simd4MaxD(Simd4Double a, Simd4Double b)
{
    return { _mm256_max_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4MinD(Simd4Double a, Simd4Double b)
{
    return { _mm256_min_pd(a.r, b.r) };
}

static inline Simd4Double gmx_simdcall
simd4RoundD(Simd4Double x)
{
    return { _mm256_round_pd(x.r, _MM_FROUND_NINT) };
}

static inline Simd4Double gmx_simdcall
simd4TruncD(Simd4Double x)
{
    return { _mm256_round_pd(x.r, _MM_FROUND_TRUNC) };
}

static inline float gmx_simdcall
simd4DotProductD(Simd4Double a, Simd4Double b)
{
    __m128d tmp1, tmp2;
    a.r  = _mm256_mul_pd(a.r, b.r);
    tmp1 = _mm256_castpd256_pd128(a.r);
    tmp2 = _mm256_extractf128_pd(a.r, 0x1);

    tmp1 = _mm_add_pd(tmp1, _mm_permute_pd(tmp1, _MM_SHUFFLE2(0, 1)));
    tmp1 = _mm_add_pd(tmp1, tmp2);
    return *reinterpret_cast<double *>(&tmp1);
}

static inline void gmx_simdcall
simd4Transpose(Simd4Double * v0, Simd4Double * v1,
               Simd4Double * v2, Simd4Double * v3)
{
    __m256d t1, t2, t3, t4;
    t1    = _mm256_unpacklo_pd(v0->r, v1->r);
    t2    = _mm256_unpackhi_pd(v0->r, v1->r);
    t3    = _mm256_unpacklo_pd(v2->r, v3->r);
    t4    = _mm256_unpackhi_pd(v2->r, v3->r);
    v0->r = _mm256_permute2f128_pd(t1, t3, 0x20);
    v1->r = _mm256_permute2f128_pd(t2, t4, 0x20);
    v2->r = _mm256_permute2f128_pd(t1, t3, 0x31);
    v3->r = _mm256_permute2f128_pd(t2, t4, 0x31);
}

static inline Simd4DBool gmx_simdcall
simd4CmpEqD(Simd4Double a, Simd4Double b)
{
    return { _mm256_cmp_pd(a.r, b.r, _CMP_EQ_OQ) };
}

static inline Simd4DBool gmx_simdcall
simd4CmpLtD(Simd4Double a, Simd4Double b)
{
    return { _mm256_cmp_pd(a.r, b.r, _CMP_LT_OQ) };
}

static inline Simd4DBool gmx_simdcall
simd4CmpLeD(Simd4Double a, Simd4Double b)
{
    return { _mm256_cmp_pd(a.r, b.r, _CMP_LE_OQ) };
}

static inline Simd4DBool gmx_simdcall
simd4AndDB(Simd4DBool a, Simd4DBool b)
{
    return { _mm256_and_pd(a.b, b.b) };
}

static inline Simd4DBool gmx_simdcall
simd4OrDB(Simd4DBool a, Simd4DBool b)
{
    return { _mm256_or_pd(a.b, b.b) };
}

static inline bool gmx_simdcall
simd4AnyTrueDB(Simd4DBool a) { return _mm256_movemask_pd(a.b) != 0; }

static inline Simd4Double gmx_simdcall
simd4MaskD(Simd4Double a, Simd4DBool mask)
{
    return { _mm256_and_pd(a.r, mask.b) };
}

static inline Simd4Double gmx_simdcall
simd4MaskNotD(Simd4Double a, Simd4DBool mask)
{
    return { _mm256_andnot_pd(mask.b, a.r) };
}

static inline Simd4Double gmx_simdcall
simd4BlendD(Simd4Double a, Simd4Double b, Simd4DBool sel)
{
    return { _mm256_blendv_pd(a.r, b.r, sel.b) };
}

static inline float gmx_simdcall
simd4ReduceD(Simd4Double a)
{
    __m128d a0, a1;
    // test with shuffle & add as an alternative to hadd later
    a.r = _mm256_hadd_pd(a.r, a.r);
    a0  = _mm256_castpd256_pd128(a.r);
    a1  = _mm256_extractf128_pd(a.r, 0x1);
    a0  = _mm_add_sd(a0, a1);
    return *reinterpret_cast<double *>(&a0);
}

} // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_256_SIMD4_DOUBLE_H
