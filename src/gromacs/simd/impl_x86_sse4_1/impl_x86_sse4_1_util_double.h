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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H
#define GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_util_double.h"

#include "impl_x86_sse4_1_simd_double.h"

namespace gmx
{

template <int align>
static inline void
gatherLoadBySimdIntTranspose(const double *  base,
                             SimdDInt32      offset,
                             SimdDouble *    v0,
                             SimdDouble *    v1,
                             SimdDouble *    v2,
                             SimdDouble *    v3)
{
    __m128d t1, t2, t3, t4;
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 4)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 8)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset.simdInternal_ = _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 16)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 4);
    }

    if (align == 4 || align == 8 || align == 12 || align == 16)
    {
        assert(std::size_t(base + extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + extract<1>(offset)) % 16 == 0);

        t1  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 1));
        t3  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 0) + 2);
        t4  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 1) + 2);
    }
    else
    {
        assert(std::size_t(base + align * extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<1>(offset)) % 16 == 0);

        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 1));
        t3  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 0) + 2);
        t4  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 1) + 2);
    }
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
    v2->simdInternal_ = _mm_unpacklo_pd(t3, t4);
    v3->simdInternal_ = _mm_unpackhi_pd(t3, t4);
}

template <int align>
static inline void
gatherLoadBySimdIntTranspose(const double *    base,
                             SimdDInt32        offset,
                             SimdDouble *      v0,
                             SimdDouble *      v1)
{
    __m128d t1, t2;
    /* Use optimized bit-shift multiply for the most common alignments */
    if (align == 2)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 1);
    }
    else if (align == 4)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 6)
    {
        /* multiply by 3, then by 2 */
        offset.simdInternal_ = _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 1);
    }
    else if (align == 8)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 3);
    }
    else if (align == 12)
    {
        /* multiply by 3, then by 4 */
        offset.simdInternal_ = _mm_add_epi32(offset.simdInternal_, _mm_slli_epi32(offset.simdInternal_, 1));
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 2);
    }
    else if (align == 16)
    {
        offset.simdInternal_ = _mm_slli_epi32(offset.simdInternal_, 4);
    }

    if (align == 2 || align == 4 || align == 6 ||
        align == 8 || align == 12 || align == 16)
    {
        assert(std::size_t(base + extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + extract<1>(offset)) % 16 == 0);

        t1  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 0));
        t2  = _mm_load_pd(base + _mm_extract_epi32(offset.simdInternal_, 1));
    }
    else
    {
        assert(std::size_t(base + align * extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<1>(offset)) % 16 == 0);

        t1  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 0));
        t2  = _mm_load_pd(base + align * _mm_extract_epi32(offset.simdInternal_, 1));
    }
    v0->simdInternal_ = _mm_unpacklo_pd(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_pd(t1, t2);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_UTIL_DOUBLE_H
