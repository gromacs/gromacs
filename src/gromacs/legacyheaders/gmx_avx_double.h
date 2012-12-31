/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
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
#ifndef _gmx_avx_double_h_
#define _gmx_avx_double_h_

/* We require AVX now! */

#include <immintrin.h> /* AVX */

static inline __m256d
gmx_mm256_invsqrt_pd(__m256d x)
{
    /* There is no double precision AVX rsqrt instruction.
     * But using a single precision rsqrt still gives the full precision.
     */
    const __m256d half    = _mm256_set_pd(0.5,0.5,0.5,0.5);
    const __m256d three   = _mm256_set_pd(3.0,3.0,3.0,3.0);

    __m256d lu = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)));

    lu = _mm256_mul_pd(half,_mm256_mul_pd(_mm256_sub_pd(three,_mm256_mul_pd(_mm256_mul_pd(lu,lu),x)),lu));
    return _mm256_mul_pd(half,_mm256_mul_pd(_mm256_sub_pd(three,_mm256_mul_pd(_mm256_mul_pd(lu,lu),x)),lu));
}

static inline __m256d
gmx_mm256_calc_rsq_pd(__m256d dx, __m256d dy, __m256d dz)
{
    return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(dx,dx), _mm256_mul_pd(dy,dy) ), _mm256_mul_pd(dz,dz) );
}

/* Normal sum of four xmm registers */
#define gmx_mm256_sum4_pd(t0,t1,t2,t3)  _mm256_add_pd(_mm256_add_pd(t0,t1),_mm256_add_pd(t2,t3))

#endif /* gmx_avx_double_h_ */
