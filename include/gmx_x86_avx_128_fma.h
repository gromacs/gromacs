/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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

#define _GMX_MM_BLEND(b3,b2,b1,b0) (((b3) << 3) | ((b2) << 2) | ((b1) << 1) | ((b0)))

#define _GMX_MM_PERMUTE128D(fp1,fp0)         (((fp1) << 1) | ((fp0)))


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {           \
    __m128d __gmx_t1 = row0;                         \
    row0           = _mm_unpacklo_pd(row0,row1);     \
    row1           = _mm_unpackhi_pd(__gmx_t1,row1); \
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

/* Gromacs wrapper defines to enable FMA instructions using either
 * - AMD FMA4
 * - Intel FMA3 (different intrinsics names)
 * - Emulated FMA for debugging
 */
#ifdef __FMA4__
#    define gmx_mm_fmadd_ps(_a,_b,_c)   _mm_macc_ps(_a,_b,_c)
#    define gmx_mm_fmsub_ps(_a,_b,_c)   _mm_msub_ps(_a,_b,_c)
#    define gmx_mm_fnmadd_ps(_a,_b,_c)  _mm_nmacc_ps(_a,_b,_c)
#    define gmx_mm_fnmsub_ps(_a,_b,_c)  _mm_nmsub_ps(_a,_b,_c)
#    define gmx_mm_fmadd_ss(_a,_b,_c)   _mm_macc_ss(_a,_b,_c)
#    define gmx_mm_fmsub_ss(_a,_b,_c)   _mm_msub_ss(_a,_b,_c)
#    define gmx_mm_fnmadd_ss(_a,_b,_c)  _mm_nmacc_ss(_a,_b,_c)
#    define gmx_mm_fnmsub_ss(_a,_b,_c)  _mm_nmsub_ss(_a,_b,_c)
#    define gmx_mm_fmadd_pd(_a,_b,_c)   _mm_macc_pd(_a,_b,_c)
#    define gmx_mm_fmsub_pd(_a,_b,_c)   _mm_msub_pd(_a,_b,_c)
#    define gmx_mm_fnmadd_pd(_a,_b,_c)  _mm_nmacc_pd(_a,_b,_c)
#    define gmx_mm_fnmsub_pd(_a,_b,_c)  _mm_nmsub_pd(_a,_b,_c)
#    define gmx_mm_fmadd_sd(_a,_b,_c)   _mm_macc_sd(_a,_b,_c)
#    define gmx_mm_fmsub_sd(_a,_b,_c)   _mm_msub_sd(_a,_b,_c)
#    define gmx_mm_fnmadd_sd(_a,_b,_c)  _mm_nmacc_sd(_a,_b,_c)
#    define gmx_mm_fnmsub_sd(_a,_b,_c)  _mm_nmsub_sd(_a,_b,_c)
#elif defined(__INTEL_COMPILER)
#    define gmx_mm_fmadd_ps(_a,_b,_c)   _mm_fmadd_ps(_a,_b,_c)
#    define gmx_mm_fmsub_ps(_a,_b,_c)   _mm_fmsub_ps(_a,_b,_c)
#    define gmx_mm_fnmadd_ps(_a,_b,_c)  _mm_fnmadd_ps(_a,_b,_c)
#    define gmx_mm_fnmsub_ps(_a,_b,_c)  _mm_fnmsub_ps(_a,_b,_c)
#    define gmx_mm_fmadd_ss(_a,_b,_c)   _mm_fmadd_ss(_a,_b,_c)
#    define gmx_mm_fmsub_ss(_a,_b,_c)   _mm_fmsub_ss(_a,_b,_c)
#    define gmx_mm_fnmadd_ss(_a,_b,_c)  _mm_fnmadd_ss(_a,_b,_c)
#    define gmx_mm_fnmsub_ss(_a,_b,_c)  _mm_fnmsub_ss(_a,_b,_c)
#    define gmx_mm_fmadd_pd(_a,_b,_c)   _mm_fmadd_pd(_a,_b,_c)
#    define gmx_mm_fmsub_pd(_a,_b,_c)   _mm_fmsub_pd(_a,_b,_c)
#    define gmx_mm_fnmadd_pd(_a,_b,_c)  _mm_fnmadd_pd(_a,_b,_c)
#    define gmx_mm_fnmsub_pd(_a,_b,_c)  _mm_fnmsub_pd(_a,_b,_c)
#    define gmx_mm_fmadd_sd(_a,_b,_c)   _mm_fmadd_sd(_a,_b,_c)
#    define gmx_mm_fmsub_sd(_a,_b,_c)   _mm_fmsub_sd(_a,_b,_c)
#    define gmx_mm_fnmadd_sd(_a,_b,_c)  _mm_fnmadd_sd(_a,_b,_c)
#    define gmx_mm_fnmsub_sd(_a,_b,_c)  _mm_fnmsub_sd(_a,_b,_c)
#elif defined GMX_EMULATE_FMA
#    define gmx_mm_fmadd_ps(_a,_b,_c)   _mm_add_ps(_c,_mm_mul_ps(_a,_b))
#    define gmx_mm_fmsub_ps(_a,_b,_c)   _mm_sub_ps(_mm_mul_ps(_a,_b),_c)
#    define gmx_mm_fnmadd_ps(_a,_b,_c)  _mm_sub_ps(_c,_mm_mul_ps(_a,_b))
#    define gmx_mm_fmadd_ss(_a,_b,_c)   _mm_add_ss(_c,_mm_mul_ss(_a,_b))
#    define gmx_mm_fmsub_ss(_a,_b,_c)   _mm_sub_ss(_mm_mul_ss(_a,_b),_c)
#    define gmx_mm_fnmadd_ss(_a,_b,_c)  _mm_sub_ss(_c,_mm_mul_ss(_a,_b))
#    define gmx_mm_fmadd_pd(_a,_b,_c)   _mm_add_pd(_c,_mm_mul_pd(_a,_b))
#    define gmx_mm_fmsub_pd(_a,_b,_c)   _mm_sub_pd(_mm_mul_pd(_a,_b),_c)
#    define gmx_mm_fnmadd_pd(_a,_b,_c)  _mm_sub_pd(_c,_mm_mul_pd(_a,_b))
#    define gmx_mm_fmadd_sd(_a,_b,_c)   _mm_add_sd(_c,_mm_mul_sd(_a,_b))
#    define gmx_mm_fmsub_sd(_a,_b,_c)   _mm_sub_sd(_mm_mul_sd(_a,_b),_c)
#    define gmx_mm_fnmadd_sd(_a,_b,_c)  _mm_sub_sd(_c,_mm_mul_sd(_a,_b))
static gmx_inline __m128
gmx_mm_fnmsub_ps( __m128 a, __m128 b, __m128 c)
{
    __m128  signbit  = gmx_mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_sub_ps(_mm_xor_ps(_mm_mul_ps(a,b),signbit),c);
}
static gmx_inline __m128
gmx_mm_fnmsub_ss( __m128 a, __m128 b, __m128 c)
{
    __m128  signbit  = gmx_mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    return _mm_sub_ss(_mm_xor_ss(_mm_mul_ss(a,b),signbit),c);
}
static gmx_inline __m128d
gmx_mm_fnmsub_pd( __m128d a, __m128d b, __m128d c)
{
    __m128d signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );
    return _mm_sub_pd(_mm_xor_pd(_mm_mul_pd(a,b),signbit),c);
}
static gmx_inline __m128d
gmx_mm_fnmsub_sd( __m128d a, __m128d b, __m128d c)
{
    __m128d signbit   = gmx_mm_castsi128_pd( _mm_set_epi32(0x80000000,0x00000000,0x80000000,0x00000000) );
    return _mm_sub_sd(_mm_xor_sd(_mm_mul_sd(a,b),signbit),c);
}
#endif


static void
gmx_mm_printxmm_ps(const char *s,__m128 xmm)
{
    float f[4];

    _mm_storeu_ps(f,xmm);
    printf("%s: %15.10e %15.10e %15.10e %15.10e\n",s,f[0],f[1],f[2],f[3]);
}


static void
gmx_mm_printxmmsum_ps(const char *s,__m128 xmm)
{
    float f[4];

    _mm_storeu_ps(f,xmm);
    printf("%s (sum): %15.10g\n",s,f[0]+f[1]+f[2]+f[3]);
}


static void
gmx_mm_printxmm_pd(const char *s,__m128d xmm)
{
    double f[2];

    _mm_storeu_pd(f,xmm);
    printf("%s: %30.20e %30.20e\n",s,f[0],f[1]);
}

static void
gmx_mm_printxmmsum_pd(const char *s,__m128d xmm)
{
    double f[2];

    _mm_storeu_pd(f,xmm);
    printf("%s (sum): %15.10g\n",s,f[0]+f[1]);
}


static void
gmx_mm_printxmm_epi32(const char *s,__m128i xmmi)
{
    int i[4];

    _mm_storeu_si128((__m128i *)i,xmmi);
    printf("%10s: %2d %2d %2d %2d\n",s,i[0],i[1],i[2],i[3]);
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



#endif /* _gmx_x86_avx_128_fma_h_ */
