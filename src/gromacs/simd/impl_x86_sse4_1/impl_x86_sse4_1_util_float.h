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

#ifndef GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H
#define GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H

#include "config.h"

#include <smmintrin.h>

#include "gromacs/simd/impl_x86_sse2/impl_x86_sse2_util_float.h"

#include "impl_x86_sse4_1_simd_float.h"

namespace gmx
{

template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *  base,
                             SimdFInt32     offset,
                             SimdFloat *    v0,
                             SimdFloat *    v1,
                             SimdFloat *    v2,
                             SimdFloat *    v3)
{
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
        assert(std::size_t(base + extract<2>(offset)) % 16 == 0);
        assert(std::size_t(base + extract<3>(offset)) % 16 == 0);

        v0->simdInternal_ = _mm_load_ps(base + _mm_extract_epi32(offset.simdInternal_, 0) );
        v1->simdInternal_ = _mm_load_ps(base + _mm_extract_epi32(offset.simdInternal_, 1) );
        v2->simdInternal_ = _mm_load_ps(base + _mm_extract_epi32(offset.simdInternal_, 2) );
        v3->simdInternal_ = _mm_load_ps(base + _mm_extract_epi32(offset.simdInternal_, 3) );
    }
    else
    {
        assert(std::size_t(base + align * extract<0>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<1>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<2>(offset)) % 16 == 0);
        assert(std::size_t(base + align * extract<3>(offset)) % 16 == 0);

        v0->simdInternal_ = _mm_load_ps(base + align * _mm_extract_epi32(offset.simdInternal_, 0) );
        v1->simdInternal_ = _mm_load_ps(base + align * _mm_extract_epi32(offset.simdInternal_, 1) );
        v2->simdInternal_ = _mm_load_ps(base + align * _mm_extract_epi32(offset.simdInternal_, 2) );
        v3->simdInternal_ = _mm_load_ps(base + align * _mm_extract_epi32(offset.simdInternal_, 3) );
    }
    _MM_TRANSPOSE4_PS(v0->simdInternal_, v1->simdInternal_, v2->simdInternal_, v3->simdInternal_);
}


template <int align>
static inline void gmx_simdcall
gatherLoadBySimdIntTranspose(const float *   base,
                             SimdFInt32      offset,
                             SimdFloat *     v0,
                             SimdFloat *     v1)
{
    __m128 t1, t2;

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
        v0->simdInternal_ = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset.simdInternal_, 0))));
        v1->simdInternal_ = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset.simdInternal_, 1))));
        t1                = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset.simdInternal_, 2))));
        t2                = _mm_castpd_ps(_mm_load_sd((const double *)(base + _mm_extract_epi32(offset.simdInternal_, 3))));
    }
    else
    {
        v0->simdInternal_ = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset.simdInternal_, 0))));
        v1->simdInternal_ = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset.simdInternal_, 1))));
        t1                = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset.simdInternal_, 2))));
        t2                = _mm_castpd_ps(_mm_load_sd((const double *)(base + align * _mm_extract_epi32(offset.simdInternal_, 3))));
    }
    t1                = _mm_unpacklo_ps(v0->simdInternal_, t1);
    t2                = _mm_unpacklo_ps(v1->simdInternal_, t2);
    v0->simdInternal_ = _mm_unpacklo_ps(t1, t2);
    v1->simdInternal_ = _mm_unpackhi_ps(t1, t2);
}

}      // namespace gmx

#endif // GMX_SIMD_IMPL_X86_SSE4_1_UTIL_FLOAT_H
