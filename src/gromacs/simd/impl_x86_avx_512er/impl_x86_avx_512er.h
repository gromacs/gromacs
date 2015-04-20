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

#ifndef GMX_SIMD_IMPL_X86_AVX_512ER_H
#define GMX_SIMD_IMPL_X86_AVX_512ER_H

#include <math.h>

#include <immintrin.h>

/* Intel AVX-512ER */

/* This implementation inherits 99% from AVX-512F, but adds extended-precision
 * lookups for 1/sqrt(x) and 1x, as well as single-precision versions of
 * exp(x) and log(x).
 */

/* Inherit most stuff from AVX-512F */
#include "gromacs/simd/impl_x86_avx_512f/impl_x86_avx_512f.h"

/* Override some AVX-512F settings */
/* Implementation details */
#undef  GMX_SIMD_RSQRT_BITS
#define GMX_SIMD_RSQRT_BITS         28
#undef  GMX_SIMD_RCP_BITS
#define GMX_SIMD_RCP_BITS           28

#undef  gmx_simd_rsqrt_f
#define gmx_simd_rsqrt_f           _mm512_rsqrt28_ps
#undef  gmx_simd_rcp_f
#define gmx_simd_rcp_f             _mm512_rcp28_ps


#undef  gmx_simd_rsqrt_d
#define gmx_simd_rsqrt_d           _mm512_rsqrt28_pd
#undef  gmx_simd_rcp_d
#define gmx_simd_rcp_d             _mm512_rcp28_pd


#undef  gmx_simd4_rsqrt_f
#define gmx_simd4_rsqrt_f(x)       _mm512_castps512_ps128(_mm512_rsqrt28_ps(_mm512_castps128_ps512(x)))

#undef  gmx_simd4_rsqrt_d
#define gmx_simd4_rsqrt_d(x)       _mm512_castpd512_pd256(_mm512_rsqrt28_pd(_mm512_castpd256_pd512(x)))

static gmx_inline __m512
gmx_simd_exp_f_x86_avx_512er(__m512 x)
{
    const gmx_simd_float_t  argscale    = gmx_simd_set1_f(1.44269504088896341f);
    const gmx_simd_float_t  invargscale = gmx_simd_set1_f(-0.69314718055994528623f);

    __m512                  xscaled = _mm512_mul_ps(x, argscale);
    __m512                  r       = _mm512_exp2a23_ps(xscaled);

    /* exp2a23_ps provides 23 bits of accuracy, but we ruin some of that with our argument
     * scaling. To correct this, we find the difference between the scaled argument and
     * the true one (extended precision arithmetics does not appear to be necessary to
     * fulfill our accuracy requirements) and then multiply by the exponent of this
     * correction since exp(a+b)=exp(a)*exp(b).
     * Note that this only adds two instructions (and maybe some constant loads).
     */
    x         = gmx_simd_fmadd_f(invargscale, xscaled, x);
    /* x will now be a _very_ small number, so approximate exp(x)=1+x.
     * We should thus apply the correction as r'=r*(1+x)=r+r*x
     */
    r         = gmx_simd_fmadd_f(r, x, r);
    return r;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512ER_H */
