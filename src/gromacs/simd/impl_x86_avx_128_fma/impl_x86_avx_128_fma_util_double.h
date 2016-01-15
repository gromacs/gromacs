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

#ifndef GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_DOUBLE_H

#include "config.h"

#include <cassert>
#include <cstddef>

#include <immintrin.h>
#include <x86intrin.h>

#include "gromacs/simd/impl_x86_sse4_1/impl_x86_sse4_1_util_double.h"

namespace gmx
{

static inline void gmx_simdcall
expandScalarsToTriplets(SimdDouble    scalar,
                        SimdDouble *  triplets0,
                        SimdDouble *  triplets1,
                        SimdDouble *  triplets2)
{
    triplets0->simdInternal_ = _mm_permute_pd(scalar.simdInternal_, _MM_SHUFFLE2(0, 0));
    triplets1->simdInternal_ = _mm_permute_pd(scalar.simdInternal_, _MM_SHUFFLE2(1, 0));
    triplets2->simdInternal_ = _mm_permute_pd(scalar.simdInternal_, _MM_SHUFFLE2(1, 1));
}

static inline double
reduceIncr4ReturnSum(double *    m,
                     SimdDouble  v0,
                     SimdDouble  v1,
                     SimdDouble  v2,
                     SimdDouble  v3)
{
    __m128d t1, t2, t3, t4;

    t1 = _mm_unpacklo_pd(v0.simdInternal_, v1.simdInternal_);
    t2 = _mm_unpackhi_pd(v0.simdInternal_, v1.simdInternal_);
    t3 = _mm_unpacklo_pd(v2.simdInternal_, v3.simdInternal_);
    t4 = _mm_unpackhi_pd(v2.simdInternal_, v3.simdInternal_);

    t1 = _mm_add_pd(t1, t2);
    t3 = _mm_add_pd(t3, t4);

    assert(std::size_t(m) % 16 == 0);

    t2 = _mm_add_pd(t1, _mm_load_pd(m));
    t4 = _mm_add_pd(t3, _mm_load_pd(m + 2));
    _mm_store_pd(m, t2);
    _mm_store_pd(m + 2, t4);

    t1 = _mm_add_pd(t1, t3);

    t2 = _mm_add_sd(t1, _mm_permute_pd(t1, _MM_SHUFFLE2(1, 1)));
    return *reinterpret_cast<double *>(&t2);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_AVX_128_FMA_UTIL_DOUBLE_H
