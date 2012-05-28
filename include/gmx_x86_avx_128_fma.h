/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifndef _gmx_x86_avx_128_fma_h_
#define _gmx_x86_avx_128_fma_h_


#include <immintrin.h>
#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h> /* FMA */
#endif


#include <stdio.h>

#ifndef gmx_inline
#  ifdef __cplusplus
#    define gmx_inline inline
#  elif defined(__GNUC__)
/* GCC */
#    define gmx_inline   __inline__
#  elif (defined(__INTEL_COMPILER) || defined(__ECC)) && defined(__ia64__)
/* ICC */
#    define gmx_inline __inline__
#  elif defined(__PATHSCALE__)
/* Pathscale */
#    define gmx_inline __inline__
#  elif defined(__PGIC__)
/* Portland */
#    define gmx_inline __inline
#  elif defined _MSC_VER
/* MSVC */
#    define gmx_inline __inline
#  elif defined(__xlC__)
/* IBM */
#    define gmx_inline __inline
#  else
#    define gmx_inline
#  endif
#endif

#define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))


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
static __m128  gmx_mm_castsi128_ps(__m128i a) { return *(__m128 *) &a;  } 
static __m128i gmx_mm_castps_si128(__m128 a)  { return *(__m128i *) &a; } 
static __m128  gmx_mm_castps_ps128(__m128 a) { return *(__m128 *) &a;  } 
static __m128d gmx_mm_castsi128_pd(__m128i a) { return *(__m128d *) &a;  } 
static __m128i gmx_mm_castpd_si128(__m128d a)  { return *(__m128i *) &a; } 
#endif

#ifndef _MSC_VER
/* The warning directive is not supported by MSVC, and that compiler
 * does not support overriding built-in functions anyway...
 */
#if !defined(HAVE_x86INTRIN_H) || !defined(__FMA4__)
#warning Emulating FMA instructions - this is probably not what you want!
/* Wrapper routines so we can do test builds on non-FMA hardware */
static __m128
_mm_macc_ps(__m128 a, __m128 b, __m128 c)
{
    
    return _mm_add_ps(c,_mm_mul_ps(a,b));
}

static __m128
_mm_nmacc_ps(__m128 a, __m128 b, __m128 c)
{
    
    return _mm_sub_ps(c,_mm_mul_ps(a,b));
}

static __m128d
_mm_macc_pd(__m128d a, __m128d b, __m128d c)
{
    
    return _mm_add_pd(c,_mm_mul_pd(a,b));
}

static __m128d
_mm_nmacc_pd(__m128d a, __m128d b, __m128d c)
{
    
    return _mm_sub_pd(c,_mm_mul_pd(a,b));
}
#endif /* FMA4 support */

#endif /* _MSC_VER */

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
