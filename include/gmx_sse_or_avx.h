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

/* Undefine all defines used below so we can include this file multiple times
 * with different settings from the same source file.
 */

#undef SSE_OR_AVX_WIDTH

#undef gmx_epi32

#undef gmx_mm_pr

#undef gmx_load_pr
#undef gmx_load1_pr
#undef gmx_set1_pr
#undef gmx_setzero_pr
#undef gmx_store_pr
#undef gmx_storeu_pr

#undef gmx_add_pr
#undef gmx_sub_pr
#undef gmx_mul_pr
#undef gmx_max_pr
#undef gmx_cmplt_pr
#undef gmx_and_pr
#undef gmx_or_pr
#undef gmx_andnot_pr

#undef gmx_floor_pr
#undef gmx_blendv_pr

#undef gmx_movemask_pr

#undef gmx_mm_castsi128_pr

#undef gmx_cvttpr_epi32
#undef gmx_cvtepi32_pr

#undef gmx_invsqrt_pr
#undef gmx_calc_rsq_pr
#undef gmx_sum4_pr


/* By defining GMX_MM128_HERE or GMX_MM256_HERE before including this file
 * the same intrinsics, with defines, can be compiled for either 128 or 256
 * bit wide SSE or AVX instructions.
 * The gmx_ prefix is replaced by _mm_ or _mm256_ (SSE or AVX).
 * The _pr suffix is replaced by _ps or _pd (single or double precision).
 * Note that compiler settings will decide if 128-bit intrinsics will
 * be translated into SSE or AVX instructions.
 */

#if !defined GMX_MM128_HERE && !defined GMX_MM256_HERE
"You should define GMX_MM128_HERE or GMX_MM256_HERE"
#endif

#if defined GMX_MM128_HERE && defined GMX_MM256_HERE
"You should not define both GMX_MM128_HERE and GMX_MM256_HERE"
#endif

#ifdef GMX_MM128_HERE

#define gmx_epi32  __m128i

#ifndef GMX_DOUBLE

#include "gmx_sse2_single.h"

#define SSE_OR_AVX_WIDTH  4

#define gmx_mm_pr  __m128

#define gmx_load_pr       _mm_load_ps
#define gmx_load1_pr      _mm_load1_ps
#define gmx_set1_pr       _mm_set1_ps
#define gmx_setzero_pr    _mm_setzero_ps
#define gmx_store_pr      _mm_store_ps
#define gmx_storeu_pr     _mm_storeu_ps

#define gmx_add_pr        _mm_add_ps
#define gmx_sub_pr        _mm_sub_ps
#define gmx_mul_pr        _mm_mul_ps
#define gmx_max_pr        _mm_max_ps
#define gmx_cmplt_pr      _mm_cmplt_ps
#define gmx_and_pr        _mm_and_ps
#define gmx_or_pr         _mm_or_ps
#define gmx_andnot_pr     _mm_andnot_ps
/* SSE 4.1 */
#define gmx_floor_pr      _mm_floor_ps
#define gmx_blendv_pr     _mm_blendv_ps

#define gmx_movemask_pr   _mm_movemask_ps

#define gmx_mm_castsi128_pr gmx_mm_castsi128_ps

#define gmx_cvttpr_epi32  _mm_cvttps_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_ps

#define gmx_invsqrt_pr    gmx_mm_invsqrt_ps
#define gmx_calc_rsq_pr   gmx_mm_calc_rsq_ps
#define gmx_sum4_pr       gmx_mm_sum4_ps

#else /* ifndef GMX_DOUBLE */

#include "gmx_sse2_double.h"

#define SSE_OR_AVX_WIDTH  2

#define gmx_mm_pr  __m128d

#define gmx_load_pr       _mm_load_pd
#define gmx_load1_pr      _mm_load1_pd
#define gmx_set1_pr       _mm_set1_pd
#define gmx_setzero_pr    _mm_setzero_pd
#define gmx_store_pr      _mm_store_pd
#define gmx_storeu_pr     _mm_storeu_pd

#define gmx_add_pr        _mm_add_pd
#define gmx_sub_pr        _mm_sub_pd
#define gmx_mul_pr        _mm_mul_pd
#define gmx_max_pr        _mm_max_pd
#define gmx_cmplt_pr      _mm_cmplt_pd
#define gmx_and_pr        _mm_and_pd
#define gmx_or_pr         _mm_or_pd
#define gmx_andnot_pr     _mm_andnot_pd
/* SSE 4.1 */
#define gmx_floor_pr      _mm_floor_pd
#define gmx_blendv_pr     _mm_blendv_pd

#define gmx_movemask_pr   _mm_movemask_pd

#define gmx_mm_castsi128_pr gmx_mm_castsi128_pd

#define gmx_cvttpr_epi32  _mm_cvttpd_epi32
#define gmx_cvtepi32_pr   _mm_cvtepi32_pd

#define gmx_invsqrt_pr    gmx_mm_invsqrt_pd
#define gmx_calc_rsq_pr   gmx_mm_calc_rsq_pd
#define gmx_sum4_pr       gmx_mm_sum4_pd

#endif /* ifndef GMX_DOUBLE */

#endif /* GMX_MM128_HERE */

#ifdef GMX_MM256_HERE

#define gmx_epi32 __m256i

#ifndef GMX_DOUBLE

#include "gmx_avx_single.h"

#define SSE_OR_AVX_WIDTH  8

#define gmx_mm_pr  __m256

#define gmx_load_pr       _mm256_load_ps
#define gmx_load1_pr(x)   _mm256_set1_ps((x)[0])
#define gmx_set1_pr       _mm256_set1_ps
#define gmx_setzero_pr    _mm256_setzero_ps
#define gmx_store_pr      _mm256_store_ps
#define gmx_storeu_pr     _mm256_storeu_ps

#define gmx_add_pr        _mm256_add_ps
#define gmx_sub_pr        _mm256_sub_ps
#define gmx_mul_pr        _mm256_mul_ps
#define gmx_max_pr        _mm256_max_ps
/* Not-equal (ordered, non-signaling)  */
#define gmx_cmpneq_pr(x,y)  _mm256_cmp_ps(x,y,0x0c)
/* Less-than (ordered, non-signaling)  */
#define gmx_cmplt_pr(x,y) _mm256_cmp_ps(x,y,0x11)
#define gmx_and_pr        _mm256_and_ps
#define gmx_or_pr         _mm256_or_ps
#define gmx_andnot_pr     _mm256_andnot_ps

#define gmx_floor_pr      _mm256_floor_ps
#define gmx_blendv_pr     _mm256_blendv_ps

#define gmx_movemask_pr   _mm256_movemask_ps

#define gmx_mm_castsi256_pr _mm256_castsi256_ps

#define gmx_cvttpr_epi32  _mm256_cvttps_epi32

#define gmx_invsqrt_pr    gmx_mm256_invsqrt_ps
#define gmx_calc_rsq_pr   gmx_mm256_calc_rsq_ps
#define gmx_sum4_pr       gmx_mm256_sum4_ps

#else

#include "gmx_avx_double.h"

#define SSE_OR_AVX_WIDTH  4

#define gmx_mm_pr  __m256d

#define gmx_load_pr       _mm256_load_pd
#define gmx_load1_pr(x)   _mm256_set1_pd((x)[0])
#define gmx_set1_pr       _mm256_set1_pd
#define gmx_setzero_pr    _mm256_setzero_pd
#define gmx_store_pr      _mm256_store_pd
#define gmx_storeu_pr     _mm256_storeu_pd

#define gmx_add_pr        _mm256_add_pd
#define gmx_sub_pr        _mm256_sub_pd
#define gmx_mul_pr        _mm256_mul_pd
#define gmx_max_pr        _mm256_max_pd
/* Not-equal (ordered, non-signaling)  */
#define gmx_cmpneq_pr(x,y)  _mm256_cmp_pd(x,y,0x0c)
/* Less-than (ordered, non-signaling)  */
#define gmx_cmplt_pr(x,y) _mm256_cmp_pd(x,y,0x11)
#define gmx_and_pr        _mm256_and_pd
#define gmx_or_pr         _mm256_or_pd
#define gmx_andnot_pr     _mm256_andnot_pd

#define gmx_floor_pr      _mm256_floor_pd
#define gmx_blendv_pr     _mm256_blendv_pd

#define gmx_movemask_pr   _mm256_movemask_pd

#define gmx_mm_castsi256_pr _mm256_castsi256_pd

#define gmx_cvttpr_epi32  _mm256_cvttpd_epi32

#define gmx_invsqrt_pr    gmx_mm256_invsqrt_pd
#define gmx_calc_rsq_pr   gmx_mm256_calc_rsq_pd
#define gmx_sum4_pr       gmx_mm256_sum4_pd

#endif

#endif /* GMX_MM256_HERE */
