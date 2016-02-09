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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <immintrin.h>
#include <x86intrin.h>

namespace gmx
{

class Simd4Double
{
    public:
        Simd4Double() {}

        Simd4Double(double d) : simdInternal_(_mm256_set1_pd(d)) {}

        // Internal utility constructor to simplify return statements
        Simd4Double(__m256d simd) : simdInternal_(simd) {}

        __m256d  simdInternal_;
};

class Simd4DBool
{
    public:
        Simd4DBool() {}

        //! \brief Construct from scalar bool
        Simd4DBool(bool b) : simdInternal_(_mm256_castsi256_pd(_mm256_set1_epi32( b ? 0xFFFFFFFF : 0))) {}

        // Internal utility constructor to simplify return statements
        Simd4DBool(__m256d simd) : simdInternal_(simd) {}

        __m256d  simdInternal_;
};

static inline Simd4Double gmx_simdcall
load4(const double *m)
{
    assert(size_t(m) % 32 == 0);
    return {
               _mm256_load_pd(m)
    };
}

static inline void gmx_simdcall
store4(double *m, Simd4Double a)
{
    assert(size_t(m) % 32 == 0);
    _mm256_store_pd(m, a.simdInternal_);
}

static inline Simd4Double gmx_simdcall
load4U(const double *m)
{
    return {
               _mm256_loadu_pd(m)
    };
}

static inline void gmx_simdcall
store4U(double *m, Simd4Double a)
{
    _mm256_storeu_pd(m, a.simdInternal_);
}

static inline Simd4Double gmx_simdcall
simd4SetZeroD()
{
    return {
               _mm256_setzero_pd()
    };
}

static inline Simd4Double gmx_simdcall
operator&(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
andNot(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_andnot_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator|(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator^(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_xor_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator+(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_add_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator-(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_sub_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
operator-(Simd4Double x)
{
    return {
               _mm256_xor_pd(x.simdInternal_, _mm256_set1_pd(GMX_DOUBLE_NEGZERO))
    };
}

static inline Simd4Double gmx_simdcall
operator*(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_mul_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm256_macc_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm256_msub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fnma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm256_nmacc_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
fnms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return {
               _mm256_nmsub_pd(a.simdInternal_, b.simdInternal_, c.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
rsqrt(Simd4Double x)
{
    return {
               _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x.simdInternal_)))
    };
}

static inline Simd4Double gmx_simdcall
abs(Simd4Double x)
{
    return {
               _mm256_andnot_pd( _mm256_set1_pd(GMX_DOUBLE_NEGZERO), x.simdInternal_ )
    };
}

static inline Simd4Double gmx_simdcall
max(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_max_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
min(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_min_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
round(Simd4Double x)
{
    return {
               _mm256_round_pd(x.simdInternal_, _MM_FROUND_NINT)
    };
}

static inline Simd4Double gmx_simdcall
trunc(Simd4Double x)
{
    return {
               _mm256_round_pd(x.simdInternal_, _MM_FROUND_TRUNC)
    };
}

static inline float gmx_simdcall
dotProduct(Simd4Double a, Simd4Double b)
{
    __m128d tmp1, tmp2;
    a.simdInternal_  = _mm256_mul_pd(a.simdInternal_, b.simdInternal_);
    tmp1             = _mm256_castpd256_pd128(a.simdInternal_);
    tmp2             = _mm256_extractf128_pd(a.simdInternal_, 0x1);

    tmp1 = _mm_add_pd(tmp1, _mm_permute_pd(tmp1, _MM_SHUFFLE2(0, 1)));
    tmp1 = _mm_add_pd(tmp1, tmp2);
    return *reinterpret_cast<double *>(&tmp1);
}

static inline void gmx_simdcall
transpose(Simd4Double * v0, Simd4Double * v1,
          Simd4Double * v2, Simd4Double * v3)
{
    __m256d t1, t2, t3, t4;
    t1                = _mm256_unpacklo_pd(v0->simdInternal_, v1->simdInternal_);
    t2                = _mm256_unpackhi_pd(v0->simdInternal_, v1->simdInternal_);
    t3                = _mm256_unpacklo_pd(v2->simdInternal_, v3->simdInternal_);
    t4                = _mm256_unpackhi_pd(v2->simdInternal_, v3->simdInternal_);
    v0->simdInternal_ = _mm256_permute2f128_pd(t1, t3, 0x20);
    v1->simdInternal_ = _mm256_permute2f128_pd(t2, t4, 0x20);
    v2->simdInternal_ = _mm256_permute2f128_pd(t1, t3, 0x31);
    v3->simdInternal_ = _mm256_permute2f128_pd(t2, t4, 0x31);
}

static inline Simd4DBool gmx_simdcall
operator==(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_EQ_OQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator!=(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_NEQ_OQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator<(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_LT_OQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator<=(Simd4Double a, Simd4Double b)
{
    return {
               _mm256_cmp_pd(a.simdInternal_, b.simdInternal_, _CMP_LE_OQ)
    };
}

static inline Simd4DBool gmx_simdcall
operator&&(Simd4DBool a, Simd4DBool b)
{
    return {
               _mm256_and_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline Simd4DBool gmx_simdcall
operator||(Simd4DBool a, Simd4DBool b)
{
    return {
               _mm256_or_pd(a.simdInternal_, b.simdInternal_)
    };
}

static inline bool gmx_simdcall
anyTrue(Simd4DBool a) { return _mm256_movemask_pd(a.simdInternal_) != 0; }

static inline Simd4Double gmx_simdcall
selectByMask(Simd4Double a, Simd4DBool mask)
{
    return {
               _mm256_and_pd(a.simdInternal_, mask.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
selectByNotMask(Simd4Double a, Simd4DBool mask)
{
    return {
               _mm256_andnot_pd(mask.simdInternal_, a.simdInternal_)
    };
}

static inline Simd4Double gmx_simdcall
blend(Simd4Double a, Simd4Double b, Simd4DBool sel)
{
    return {
               _mm256_blendv_pd(a.simdInternal_, b.simdInternal_, sel.simdInternal_)
    };
}

static inline float gmx_simdcall
reduce(Simd4Double a)
{
    __m128d a0, a1;
    /* test with shuffle & add as an alternative to hadd later */
    a.simdInternal_ = _mm256_hadd_pd(a.simdInternal_, a.simdInternal_);
    a0              = _mm256_castpd256_pd128(a.simdInternal_);
    a1              = _mm256_extractf128_pd(a.simdInternal_, 0x1);
    a0              = _mm_add_sd(a0, a1);
    return *reinterpret_cast<double *>(&a0);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_128_FMA_SIMD4_DOUBLE_H
