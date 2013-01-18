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
#ifndef _gmx_x86_avx_128_fma_h_
#define _gmx_x86_avx_128_fma_h_


#include <immintrin.h>
#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h> /* FMA */
#endif
#ifdef HAVE_INTRIN_H
#include <intrin.h> /* FMA MSVC */
#endif


#include <stdio.h>

#include "types/simple.h"


#define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))

#define _GMX_MM_BLEND(b3, b2, b1, b0) (((b3) << 3) | ((b2) << 2) | ((b1) << 1) | ((b0)))

#define _GMX_MM_PERMUTE128D(fp1, fp0)         (((fp1) << 1) | ((fp0)))


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {           \
        __m128d __gmx_t1 = row0;                         \
        row0           = _mm_unpacklo_pd(row0, row1);     \
        row1           = _mm_unpackhi_pd(__gmx_t1, row1); \
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

#if GMX_EMULATE_AMD_FMA
/* Wrapper routines so we can do test builds on non-FMA or non-AMD hardware */
static __m128
_mm_macc_ps(__m128 a, __m128 b, __m128 c)
{
    return _mm_add_ps(c, _mm_mul_ps(a, b));
}

static __m128
_mm_nmacc_ps(__m128 a, __m128 b, __m128 c)
{
    return _mm_sub_ps(c, _mm_mul_ps(a, b));
}

static __m128
_mm_msub_ps(__m128 a, __m128 b, __m128 c)
{
    return _mm_sub_ps(_mm_mul_ps(a, b), c);
}

static __m128d
_mm_macc_pd(__m128d a, __m128d b, __m128d c)
{
    return _mm_add_pd(c, _mm_mul_pd(a, b));
}

static __m128d
_mm_nmacc_pd(__m128d a, __m128d b, __m128d c)
{
    return _mm_sub_pd(c, _mm_mul_pd(a, b));
}

static __m128d
_mm_msub_pd(__m128d a, __m128d b, __m128d c)
{
    return _mm_sub_pd(_mm_mul_pd(a, b), c);
}
#endif /* AMD FMA emulation support */

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



#endif /* _gmx_x86_avx_128_fma_h_ */
