/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of GROMACS.
 * Copyright (c) 2012-  
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _gmx_x86_simd256_single_h_
#define _gmx_x86_simd256_single_h_

/* This file includes the highest possible level of x86 (math) acceleration */

#ifdef GMX_X86_AVX_256
#include "gmx_x86_avx_256.h"
#include "gmx_math_x86_avx_256_single.h"
#else
#ifdef GMX_X86_AVX_128_FMA
#include "gmx_x86_avx_128_fma.h"
#include "gmx_math_x86_avx_128_fma_single.h"
#else
#ifdef GMX_X86_SSE4_1
#include "gmx_x86_sse4_1.h"
#include "gmx_math_x86_sse4_1_single.h"
#else
#ifdef GMX_X86_SSE2
#include "gmx_x86_sse2.h"
#include "gmx_math_x86_sse2_single.h"
#else
#error No x86 acceleration defined
#endif
#endif
#endif
#endif


static inline __m128
gmx_mm_calc_rsq_ps(__m128 dx, __m128 dy, __m128 dz)
{
    return _mm_add_ps( _mm_add_ps( _mm_mul_ps(dx,dx), _mm_mul_ps(dy,dy) ), _mm_mul_ps(dz,dz) );
}

/* Normal sum of four __m128 registers */
#define gmx_mm_sum4_ps(t0,t1,t2,t3)  _mm_add_ps(_mm_add_ps(t0,t1),_mm_add_ps(t2,t3))

#ifdef GMX_X86_AVX_256

static inline __m256
gmx_mm256_calc_rsq_ps(__m256 dx, __m256 dy, __m256 dz)
{
    return _mm256_add_ps( _mm256_add_ps( _mm256_mul_ps(dx,dx), _mm256_mul_ps(dy,dy) ), _mm256_mul_ps(dz,dz) );
}

/* Normal sum of four __m256 registers */
#define gmx_mm256_sum4_ps(t0,t1,t2,t3)  _mm256_add_ps(_mm256_add_ps(t0,t1),_mm256_add_ps(t2,t3))

#endif

#define GMX_MM_IPROD_PS(ax,ay,az,bx,by,bz)                 \
    _mm_add_ps(_mm_add_ps(_mm_mul_ps(ax,bx),_mm_mul_ps(ay,by)),_mm_mul_ps(az,bz))

#define GMX_MM_NORM2_PS(ax,ay,az) GMX_MM_IPROD_PS(ax,ay,az,ax,ay,az)

#define GMX_MM_CPROD_PS(ax,ay,az,bx,by,bz,cx,cy,cz)        \
    {                                                          \
    cx = _mm_sub_ps(_mm_mul_ps(ay,bz),_mm_mul_ps(az,by));  \
    cy = _mm_sub_ps(_mm_mul_ps(az,bx),_mm_mul_ps(ax,bz));  \
    cz = _mm_sub_ps(_mm_mul_ps(ax,by),_mm_mul_ps(ay,bx));  \
    }

#endif /* _gmx_x86_simd256_single_h_ */
