/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifndef _gmx_x86_avx_256_h_
#define _gmx_x86_avx_256_h_


#include <immintrin.h>
#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h> /* FMA */
#endif


#include <stdio.h>

#include "types/simple.h"


#define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

#define _GMX_MM_BLEND256D(b3, b2, b1, b0) (((b3) << 3) | ((b2) << 2) | ((b1) << 1) | ((b0)))
#define _GMX_MM_PERMUTE(fp3, fp2, fp1, fp0) (((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0)))
#define _GMX_MM_PERMUTE256D(fp3, fp2, fp1, fp0) (((fp3) << 3) | ((fp2) << 2) | ((fp1) << 1) | ((fp0)))
#define _GMX_MM_PERMUTE128D(fp1, fp0)         (((fp1) << 1) | ((fp0)))


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {           \
        __m128d __gmx_t1 = row0;                         \
        row0           = _mm_unpacklo_pd(row0, row1);     \
        row1           = _mm_unpackhi_pd(__gmx_t1, row1); \
}

#define GMX_MM256_FULLTRANSPOSE4_PD(row0, row1, row2, row3) \
    {                                                        \
        __m256d _t0, _t1, _t2, _t3;                          \
        _t0  = _mm256_unpacklo_pd((row0), (row1));           \
        _t1  = _mm256_unpackhi_pd((row0), (row1));           \
        _t2  = _mm256_unpacklo_pd((row2), (row3));           \
        _t3  = _mm256_unpackhi_pd((row2), (row3));           \
        row0 = _mm256_permute2f128_pd(_t0, _t2, 0x20);       \
        row1 = _mm256_permute2f128_pd(_t1, _t3, 0x20);       \
        row2 = _mm256_permute2f128_pd(_t0, _t2, 0x31);       \
        row3 = _mm256_permute2f128_pd(_t1, _t3, 0x31);       \
    }

#if (defined (_MSC_VER) || defined(__INTEL_COMPILER))
#  define gmx_mm_castsi128_ps(a) _mm_castsi128_ps(a)
#  define gmx_mm_castps_si128(a) _mm_castps_si128(a)
#  define gmx_mm_castps_ps128(a) (a)
#  define gmx_mm_castsi128_pd(a) _mm_castsi128_pd(a)
#  define gmx_mm_castpd_si128(a) _mm_castpd_si128(a)
#elif defined(__GNUC__)
#  define gmx_mm_castsi128_ps(a) ((__m128)(a))
#  define gmx_mm_castps_si128(a) ((__m128i)(a))
#  define gmx_mm_castps_ps128(a) ((__m128)(a))
#  define gmx_mm_castsi128_pd(a) ((__m128d)(a))
#  define gmx_mm_castpd_si128(a) ((__m128i)(a))
#else
static __m128  gmx_mm_castsi128_ps(__m128i a)
{
    return *(__m128 *) &a;
}
static __m128i gmx_mm_castps_si128(__m128 a)
{
    return *(__m128i *) &a;
}
static __m128  gmx_mm_castps_ps128(__m128 a)
{
    return *(__m128 *) &a;
}
static __m128d gmx_mm_castsi128_pd(__m128i a)
{
    return *(__m128d *) &a;
}
static __m128i gmx_mm_castpd_si128(__m128d a)
{
    return *(__m128i *) &a;
}
#endif

static gmx_inline __m256
gmx_mm256_unpack128lo_ps(__m256 xmm1, __m256 xmm2)
{
    return _mm256_permute2f128_ps(xmm1, xmm2, 0x20);
}

static gmx_inline __m256
gmx_mm256_unpack128hi_ps(__m256 xmm1, __m256 xmm2)
{
    return _mm256_permute2f128_ps(xmm1, xmm2, 0x31);
}

static gmx_inline __m256
gmx_mm256_set_m128(__m128 hi, __m128 lo)
{
    return _mm256_insertf128_ps(_mm256_castps128_ps256(lo), hi, 0x1);
}


static gmx_inline __m256
gmx_mm256_load4_ps(float const * p)
{
    __m128 a;

    a = _mm_load_ps(p);
    return _mm256_insertf128_ps(_mm256_castps128_ps256(a), a, 0x1);
}


static __m256d
gmx_mm256_unpack128lo_pd(__m256d xmm1, __m256d xmm2)
{
    return _mm256_permute2f128_pd(xmm1, xmm2, 0x20);
}

static __m256d
gmx_mm256_unpack128hi_pd(__m256d xmm1, __m256d xmm2)
{
    return _mm256_permute2f128_pd(xmm1, xmm2, 0x31);
}

static __m256d
gmx_mm256_set_m128d(__m128d hi, __m128d lo)
{
    return _mm256_insertf128_pd(_mm256_castpd128_pd256(lo), hi, 0x1);
}


static __m128 gmx_mm256_sum4h_m128(__m256 x, __m256 y)
{
    __m256 sum;

    sum = _mm256_add_ps(x, y);
    return _mm_add_ps(_mm256_castps256_ps128(sum), _mm256_extractf128_ps(sum, 0x1));
}


static void
gmx_mm_printxmm_ps(const char *s, __m128 xmm)
{
    float f[4];

    _mm_storeu_ps(f, xmm);
    printf("%s: %15.10e %15.10e %15.10e %15.10e\n", s, f[0], f[1], f[2], f[3]);
}


static void
gmx_mm_printxmmsum_ps(const char *s, __m128 xmm)
{
    float f[4];

    _mm_storeu_ps(f, xmm);
    printf("%s (sum): %15.10g\n", s, f[0]+f[1]+f[2]+f[3]);
}


static void
gmx_mm_printxmm_pd(const char *s, __m128d xmm)
{
    double f[2];

    _mm_storeu_pd(f, xmm);
    printf("%s: %30.20e %30.20e\n", s, f[0], f[1]);
}

static void
gmx_mm_printxmmsum_pd(const char *s, __m128d xmm)
{
    double f[2];

    _mm_storeu_pd(f, xmm);
    printf("%s (sum): %15.10g\n", s, f[0]+f[1]);
}


static void
gmx_mm_printxmm_epi32(const char *s, __m128i xmmi)
{
    int i[4];

    _mm_storeu_si128((__m128i *)i, xmmi);
    printf("%10s: %2d %2d %2d %2d\n", s, i[0], i[1], i[2], i[3]);
}

static void
gmx_mm256_printymm_ps(const char *s, __m256 ymm)
{
    float f[8];

    _mm256_storeu_ps(f, ymm);
    printf("%s: %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n", s, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);
}

static void
gmx_mm256_printymmsum_ps(const char *s, __m256 ymm)
{
    float f[8];

    _mm256_storeu_ps(f, ymm);
    printf("%s (sum): %15.10g\n", s, f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]);
}


static void
gmx_mm256_printymm_pd(const char *s, __m256d ymm)
{
    double f[4];

    _mm256_storeu_pd(f, ymm);
    printf("%s: %16.12f %16.12f %16.12f %16.12f\n", s, f[0], f[1], f[2], f[3]);
}

static void
gmx_mm256_printymmsum_pd(const char *s, __m256d ymm)
{
    double f[4];

    _mm256_storeu_pd(f, ymm);
    printf("%s (sum): %15.10g\n", s, f[0]+f[1]+f[2]+f[3]);
}



static void
gmx_mm256_printymm_epi32(const char *s, __m256i ymmi)
{
    int i[8];

    _mm256_storeu_si256((__m256i *)i, ymmi);
    printf("%10s: %2d %2d %2d %2d %2d %2d %2d %2d\n", s, i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]);
}



static int gmx_mm_check_and_reset_overflow(void)
{
    int MXCSR;
    int sse_overflow;

    MXCSR = _mm_getcsr();
    /* The overflow flag is bit 3 in the register */
    if (MXCSR & 0x0008)
    {
        sse_overflow = 1;
        /* Set the overflow flag to zero */
        MXCSR = MXCSR & 0xFFF7;
        _mm_setcsr(MXCSR);
    }
    else
    {
        sse_overflow = 0;
    }

    return sse_overflow;
}

/* Work around gcc bug with wrong type for mask formal parameter to maskload/maskstore */
#ifdef GMX_X86_AVX_GCC_MASKLOAD_BUG
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), _mm_castsi128_ps(mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), _mm_castsi128_ps(mask), (x))
#    define gmx_mm256_maskload_ps(mem, mask)    _mm256_maskload_ps((mem), _mm256_castsi256_ps(mask))
#    define gmx_mm256_maskstore_ps(mem, mask, x) _mm256_maskstore_ps((mem), _mm256_castsi256_ps(mask), (x))
#else
#    define gmx_mm_maskload_ps(mem, mask)       _mm_maskload_ps((mem), (mask))
#    define gmx_mm_maskstore_ps(mem, mask, x)    _mm_maskstore_ps((mem), (mask), (x))
#    define gmx_mm256_maskload_ps(mem, mask)    _mm256_maskload_ps((mem), (mask))
#    define gmx_mm256_maskstore_ps(mem, mask, x) _mm256_maskstore_ps((mem), (mask), (x))
#endif


#endif /* _gmx_x86_avx_256_h_ */
