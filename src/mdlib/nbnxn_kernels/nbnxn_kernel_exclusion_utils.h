/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
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
#ifndef _nbnxn_kernel_exclusion_utils_h_
#define _nbnxn_kernel_exclusion_utils_h_

/* This files contains all typedefs and static inline functions for
   handling diagonal, Newton and topology exclusions in the SIMD
   kernels. The details of the implementations depend on the feature sets
   of the SIMD hardware. */

#ifdef GMX_SIMD_REFERENCE_PLAIN_C

/* Code for handling loading exclusions and converting them into
   interactions. The x86 code might use either integer- or real-type
   SIMD, but the reference code does not need to know. */

typedef gmx_simd_ref_epi32            gmx_simd_ref_exclfilter;
#define gmx_exclfilter                gmx_simd_ref_exclfilter
static const unsigned FILTER_STRIDE = GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE;

#define gmx_load1_exclfilter(e)       gmx_simd_ref_load1_exclfilter(e)
#define gmx_load_exclusion_filter(e)  gmx_simd_ref_load_exclusion_filter((int *) e)
#define gmx_checkbitmask_pb(m0, m1)   gmx_simd_ref_checkbitmask_pb(m0, m1)

static gmx_inline gmx_simd_ref_exclfilter
gmx_simd_ref_load1_exclfilter(int src)
{
    gmx_simd_ref_exclfilter a;
    int                     i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = src;
    }

    return a;
}

static gmx_inline gmx_simd_ref_exclfilter
gmx_simd_ref_load_exclusion_filter(const unsigned *src)
{
    gmx_simd_ref_exclfilter a;
    int                     i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        a.r[i] = src[i];
    }

    return a;
}

/* For topology exclusion-pair checking we need: ((a & b) ? True :
 * False). The x86 implementations use hardware-suitable integer-
 * and/or real-valued SIMD operations and a bit-wise "and" to do
 * this. The reference implementation normally uses logical operations
 * for logic, but in this case the i- and j-atom exclusion masks
 * computed during searching expect to be combined with bit-wise
 * "and".
 *
 * If the same bit is set in both input masks, return TRUE, else
 * FALSE. This function is only called with a single bit set in b.
 */
static gmx_inline gmx_simd_ref_pb
gmx_simd_ref_checkbitmask_pb(gmx_simd_ref_exclfilter a, gmx_simd_ref_exclfilter b)
{
    gmx_simd_ref_pb c;
    int             i;

    for (i = 0; i < GMX_SIMD_REF_WIDTH; i++)
    {
        c.r[i] = ((a.r[i] & b.r[i]) != 0);
    }

    return c;
}

#else /* GMX_SIMD_REFERENCE_PLAIN_C */

#ifdef GMX_X86_SSE2

#if !defined GMX_X86_AVX_256 || defined GMX_USE_HALF_WIDTH_SIMD_HERE
#ifndef GMX_DOUBLE

#define gmx_exclfilter gmx_epi32
static const unsigned FILTER_STRIDE = GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE;

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm_set1_epi32(e);
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    return _mm_load_si128((__m128i *) i);
}

static gmx_inline gmx_mm_pb
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return gmx_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()));
}

#else /* GMX_DOUBLE */

#define gmx_exclfilter gmx_epi32
static const unsigned FILTER_STRIDE = GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE;

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm_set1_epi32(e);
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    return _mm_load_si128((__m128i *) i);
}

static gmx_inline gmx_mm_pb
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return gmx_mm_castsi128_pd(_mm_cmpeq_epi32(_mm_andnot_si128(m0, m1), _mm_setzero_si128()));
}

#endif /* GMX_DOUBLE */

#else
/* We have GMX_X86_AVX_256 and not GMX_USE_HALF_WIDTH_SIMD_HERE,
 * so we use 256-bit SIMD.
 */

#ifndef GMX_DOUBLE

#define gmx_exclfilter gmx_mm_pr
static const unsigned FILTER_STRIDE = 1;

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm256_castsi256_ps(_mm256_set1_epi32(e));
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    return gmx_load_pr((real *) (i));
}

static gmx_inline gmx_mm_pb
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    return _mm256_cmp_ps(_mm256_cvtepi32_ps(_mm256_castps_si256(_mm256_and_ps(m0, m1))), _mm256_setzero_ps(), 0x0c);
}

#else /* GMX_DOUBLE */

#define gmx_exclfilter gmx_mm_pr
static const unsigned FILTER_STRIDE = 2;

static gmx_inline gmx_exclfilter
gmx_load1_exclfilter(int e)
{
    return _mm256_castsi256_pd(_mm256_set1_epi32(e));
}

static gmx_inline gmx_exclfilter
gmx_load_exclusion_filter(const unsigned *i)
{
    return gmx_load_pr((real *) (i));
}

static gmx_inline gmx_mm_pb
gmx_checkbitmask_pb(gmx_exclfilter m0, gmx_exclfilter m1)
{
    /* With <= 16 bits used the cast and conversion should not be
     * required, since only mantissa bits are set and that would give
     * a non-zero float, but with the Intel compiler this does not
     * work correctly. Because AVX does not have int->double
     * conversion, we convert via float. */
    return _mm256_cmp_pd(_mm256_castps_pd(_mm256_cvtepi32_ps(_mm256_castpd_si256(_mm256_and_pd(m0, m1)))), _mm256_setzero_pd(), 0x0c);
}

#endif /* GMX_DOUBLE */
#endif
#endif /* GMX_X86_SSE2 */
#endif /* GMX_SIMD_REFERENCE_PLAIN_C */


#endif /* _nbnxn_kernel_exclusion_utils_h_ */
