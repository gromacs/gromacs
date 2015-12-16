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

#ifndef GMX_SIMD_IMPL_X86_AVX_512ER_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_X86_AVX_512ER_SIMD_FLOAT_H

#include "config.h"

#include <immintrin.h>

#include "impl_x86_avx_512er_common.h"

#undef  simdRsqrtF
#define simdRsqrtF           _mm512_rsqrt28_ps

#undef  simdRcpF
#define simdRcpF             _mm512_rcp28_ps

#undef  simdExpF
#define simdExpF(x)          simdExpF_x86_avx_512er(x)

/* Implementation helper */

static inline __m512 gmx_simdcall
simdExpF_x86_avx_512er(__m512 x)
{
    const SimdFloat         argscale    = simdSet1F(1.44269504088896341f);
    const SimdFloat         invargscale = simdSet1F(-0.69314718055994528623f);

    __m512                  xscaled = _mm512_mul_ps(x, argscale);
    __m512                  r       = _mm512_exp2a23_ps(xscaled);

    /* exp2a23_ps provides 23 bits of accuracy, but we ruin some of that with our argument
     * scaling. To correct this, we find the difference between the scaled argument and
     * the true one (extended precision arithmetics does not appear to be necessary to
     * fulfill our accuracy requirements) and then multiply by the exponent of this
     * correction since exp(a+b)=exp(a)*exp(b).
     * Note that this only adds two instructions (and maybe some constant loads).
     */
    x         = simdFmaddF(invargscale, xscaled, x);
    /* x will now be a _very_ small number, so approximate exp(x)=1+x.
     * We should thus apply the correction as r'=r*(1+x)=r+r*x
     */
    r         = simdFmaddF(r, x, r);
    return r;
}

#endif /* GMX_SIMD_IMPL_X86_AVX_512ER_SIMD_FLOAT_H */
