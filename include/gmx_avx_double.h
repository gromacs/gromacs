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
    const __m256d half    = _mm256_set_pd(0.5, 0.5, 0.5, 0.5);
    const __m256d three   = _mm256_set_pd(3.0, 3.0, 3.0, 3.0);

    __m256d       lu = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)));

    lu = _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu, lu), x)), lu));
    return _mm256_mul_pd(half, _mm256_mul_pd(_mm256_sub_pd(three, _mm256_mul_pd(_mm256_mul_pd(lu, lu), x)), lu));
}

static inline __m256d
gmx_mm256_calc_rsq_pd(__m256d dx, __m256d dy, __m256d dz)
{
    return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy) ), _mm256_mul_pd(dz, dz) );
}

/* Normal sum of four xmm registers */
#define gmx_mm256_sum4_pd(t0, t1, t2, t3)  _mm256_add_pd(_mm256_add_pd(t0, t1), _mm256_add_pd(t2, t3))

#endif /* gmx_avx_double_h_ */
