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
#ifndef _gmx_x86_simd_double_h_
#define _gmx_x86_simd_double_h_

/* This file includes the highest possible level of x86 (math) acceleration */

#include "../utility/gmx_header_config.h"

#ifdef GMX_X86_AVX_256
#include "gmx_x86_avx_256.h"
#include "gmx_math_x86_avx_256_double.h"
#else
#ifdef GMX_X86_AVX_128_FMA
#include "gmx_x86_avx_128_fma.h"
#include "gmx_math_x86_avx_128_fma_double.h"
#else
#ifdef GMX_X86_SSE4_1
#include "gmx_x86_sse4_1.h"
#include "gmx_math_x86_sse4_1_double.h"
#else
#ifdef GMX_X86_SSE2
#include "gmx_x86_sse2.h"
#include "gmx_math_x86_sse2_double.h"
#else
#error No x86 acceleration defined
#endif
#endif
#endif
#endif

static inline __m128d
gmx_mm_calc_rsq_pd(__m128d dx, __m128d dy, __m128d dz)
{
    return _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx), _mm_mul_pd(dy,dy) ), _mm_mul_pd(dz,dz) );
}

/* Normal sum of four __m128d registers */
#define gmx_mm_sum4_pd(t0,t1,t2,t3)  _mm_add_pd(_mm_add_pd(t0,t1),_mm_add_pd(t2,t3))

#ifdef GMX_X86_AVX_256

static inline __m256d
gmx_mm256_calc_rsq_pd(__m256d dx, __m256d dy, __m256d dz)
{
    return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(dx,dx), _mm256_mul_pd(dy,dy) ), _mm256_mul_pd(dz,dz) );
}

/* Normal sum of four xmm registers */
#define gmx_mm256_sum4_pd(t0,t1,t2,t3)  _mm256_add_pd(_mm256_add_pd(t0,t1),_mm256_add_pd(t2,t3))

#endif

#endif /* _gmx_x86_simd_double_h_ */
