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
#ifndef _gmx_x86_avx_256_h_
#define _gmx_x86_avx_256_h_


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

#define _GMX_MM_PERMUTE(fp3,fp2,fp1,fp0) (((fp3) << 6) | ((fp2) << 4) | ((fp1) << 2) | ((fp0)))
#define _GMX_MM_PERMUTE256D(fp3,fp2,fp1,fp0) (((fp3) << 3) | ((fp2) << 2) | ((fp1) << 1) | ((fp0)))
#define _GMX_MM_PERMUTE128D(fp1,fp0)         (((fp1) << 1) | ((fp0)))


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {           \
    __m128d __gmx_t1 = row0;                         \
    row0           = _mm_unpacklo_pd(row0,row1);     \
    row1           = _mm_unpackhi_pd(__gmx_t1,row1); \
}



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

static void
gmx_mm256_printymm_ps(const char *s,__m256 ymm)
{
	float f[8];
	
	_mm256_storeu_ps(f,ymm);
	printf("%s: %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n",s,f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7]);	
}

static void
gmx_mm256_printymmsum_ps(const char *s,__m256 ymm)
{
	float f[8];
	
	_mm256_storeu_ps(f,ymm);
	printf("%s (sum): %15.10g\n",s,f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]);	
}


static void
gmx_mm256_printymm_pd(const char *s,__m256d ymm)
{
	double f[4];
	
	_mm256_storeu_pd(f,ymm);
	printf("%s: %16.12f %16.12f %16.12f %16.12f\n",s,f[0],f[1],f[2],f[3]);	
}

static void
gmx_mm256_printymmsum_pd(const char *s,__m256d ymm)
{
	double f[4];
	
	_mm256_storeu_pd(f,ymm);
	printf("%s (sum): %15.10g\n",s,f[0]+f[1]+f[2]+f[3]);	
}



static void
gmx_mm256_printymm_epi32(const char *s,__m256i ymmi)
{
    int i[8];
    
    _mm256_storeu_si256((__m256i *)i,ymmi);
    printf("%10s: %2d %2d %2d %2d %2d %2d %2d %2d\n",s,i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);      
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



#endif /* _gmx_x86_avx_256_h_ */
