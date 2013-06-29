/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS Development Team
 * Copyright (c) 2012,2013 by the GROMACS development team, led by
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

/*! \file
 * \brief
 * The types in this file are intended to be used for supporting
 * writing new architecture-independent SIMD intrinsics code. To
 * support a new architecture, defining types equivalent to those here
 * and and implementing functions that use them should be (nearly) all
 * that is needed.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_simd
 */

#ifndef _gmx_simd_types_h_
#define _gmx_simd_types_h_

#include "typedefs.h"


/* For testing the plain-C reference "SIMD" code, uncomment the next
 * line and modify src/gromacs/legacyheaders/types/nb_verlet.h. Make
 * sure there is no other SIMD active (i.e. cmake
 * -DGMX_CPU_ACCELERATION=None). */
/* #define GMX_SIMD_REFERENCE_PLAIN_C */
#ifdef GMX_SIMD_REFERENCE_PLAIN_C
/* float/double SIMD register type */
#define gmx_mm_pr  gmx_simd_ref_pr

/* boolean SIMD register type */
#define gmx_mm_pb  gmx_simd_ref_pb

/* integer SIMD register type, only for table indexing and exclusion masks */
#define gmx_epi32  gmx_simd_ref_epi32
#define GMX_SIMD_EPI32_WIDTH  GMX_SIMD_REF_EPI32_WIDTH

/* Type for handling loading and applying exclusion masks. */
#define gmx_exclmask               gmx_simd_ref_exclmask
#define EXCL_MASK_STRIDE           (GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE)

#endif /* GMX_SIMD_REFERENCE_PLAIN_C */


#ifdef GMX_X86_SSE2
/* This is for general x86 SIMD instruction sets that also support SSE2 */
#define GMX_HAVE_SIMD_MACROS

/* Include the highest supported x86 SIMD intrisics + math functions */
#ifdef GMX_X86_AVX_256
#include "gmx_x86_avx_256.h"
#else
#ifdef GMX_X86_AVX_128_FMA
#include "gmx_x86_avx_128_fma.h"
#else
#ifdef GMX_X86_SSE4_1
#include "gmx_x86_sse4_1.h"
#else
#ifdef GMX_X86_SSE2
#include "gmx_x86_sse2.h"
#else
#error No x86 acceleration defined
#endif
#endif
#endif
#endif
/* exp and trigonometric functions are included above */
#define GMX_SIMD_HAVE_EXP
#define GMX_SIMD_HAVE_TRIGONOMETRIC

#if !defined GMX_X86_AVX_256 || defined GMX_USE_HALF_WIDTH_SIMD_HERE

#ifndef GMX_DOUBLE

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m128

#define gmx_mm_pb  __m128

#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

#define gmx_exclmask                gmx_epi32
#define EXCL_MASK_STRIDE            (GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE)

#else /* ifndef GMX_DOUBLE */

#define GMX_SIMD_WIDTH_HERE  2

#define gmx_mm_pr  __m128d

#define gmx_mm_pb  __m128d

#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

#define gmx_exclmask                gmx_epi32
#define EXCL_MASK_STRIDE            (GMX_SIMD_EPI32_WIDTH/GMX_SIMD_WIDTH_HERE)

#endif /* ifndef GMX_DOUBLE */

#else
/* We have GMX_X86_AVX_256 and not GMX_USE_HALF_WIDTH_SIMD_HERE,
 * so we use 256-bit SIMD.
 */

#ifndef GMX_DOUBLE

#define GMX_SIMD_WIDTH_HERE  8

#define gmx_mm_pr  __m256

#define gmx_mm_pb  __m256

#define gmx_epi32  __m256i
#define GMX_SIMD_EPI32_WIDTH  8

#define gmx_exclmask                gmx_mm_pr
#ifdef GMX_DOUBLE
#define EXCL_MASK_STRIDE            2
#else
#define EXCL_MASK_STRIDE            1
#endif

#else

#define GMX_SIMD_WIDTH_HERE  4

#define gmx_mm_pr  __m256d

#define gmx_mm_pb  __m256d

/* We use 128-bit integer registers because of missing 256-bit operations */
#define gmx_epi32  __m128i
#define GMX_SIMD_EPI32_WIDTH  4

#define gmx_exclmask                gmx_mm_pr
#ifdef GMX_DOUBLE
#define EXCL_MASK_STRIDE            2
#else
#define EXCL_MASK_STRIDE            1
#endif

#endif /* GMX_DOUBLE */

#endif /* 128- or 256-bit x86 SIMD */

#endif /* GMX_X86_SSE2 */

#endif /* _gmx_simd_types_h_ */
