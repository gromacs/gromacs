/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* We require SSE2 now! */

#include <math.h>
#include <xmmintrin.h> /* SSE */
#include <emmintrin.h> /* SSE2 */
#ifdef GMX_SSE3
#  include <pmmintrin.h> /* SSE3 */
#endif
#ifdef GMX_SSE4
#  include <smmintrin.h> /* SSE4.1 */
#endif


#include <stdio.h>
                                      


static inline void
gmx_printxmm_pd(const char *s,__m128d xmm)
{
	double f[2];
	
	_mm_storeu_pd(f,xmm);
	printf("%s: %15.10g %15.10g\n",s,f[0],f[1]);	
}


#define GMX_MM_TRANSPOSE2_PD(row0, row1) {            \
     __m128d __gmx_t1 = row0;                         \
     row0           = _mm_unpacklo_pd(row0,row1);     \
     row1           = _mm_unpackhi_pd(__gmx_t1,row1); \
}


#define GMX_MM_LOAD_1JCOORD_1ATOM_PD(ptr1,jx1,jy1,jz1) {               \
	 jx1            = _mm_loadu_pd(ptr1);                              \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                        \
     jz1            = _mm_load_sd(ptr1+2);                             \
}


#define GMX_MM_LOAD_1JCOORD_2ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
}


#define GMX_MM_LOAD_1JCOORD_3ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
	 jx3            = _mm_loadu_pd(ptr1+6);                           \
     jy3            = _mm_unpackhi_pd(jx3,jx3);                       \
     jz3            = _mm_load_sd(ptr1+8);                            \
}


#define GMX_MM_LOAD_1JCOORD_4ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_pd(ptr1);                             \
     jy1            = _mm_unpackhi_pd(jx1,jx1);                       \
     jz1            = _mm_loadu_pd(ptr1+2);                           \
     jx2            = _mm_unpackhi_pd(jz1,jz1);                       \
     jy2            = _mm_loadu_pd(ptr1+4);                           \
     jz2            = _mm_unpackhi_pd(jy2,jy2);                       \
	 jx3            = _mm_loadu_pd(ptr1+6);                           \
     jy3            = _mm_unpackhi_pd(jx3,jx3);                       \
     jz3            = _mm_loadu_pd(ptr1+8);                           \
     jx4            = _mm_unpackhi_pd(jz3,jz3);                       \
     jy4            = _mm_loadu_pd(ptr1+10);                          \
     jz4            = _mm_unpackhi_pd(jy4,jy4);                       \
}

#define GMX_MM_LOAD_2JCOORD_1ATOM_PD(ptr1,ptr2,jx1,jy1,jz1) {   \
     __m128d _tmp1;                                             \
	 _tmp1           = _mm_loadu_pd(ptr1);                      \
     jy1             = _mm_loadu_pd(ptr2);                      \
     jz1             = _mm_load_sd(ptr1+2);                     \
     jx1             = _mm_unpacklo_pd(_tmp1,jy1);              \
     jy1             = _mm_unpackhi_pd(_tmp1,jy1);              \
     jz1             = _mm_loadh_pd(jz1,ptr2+2);                \
}


#define GMX_MM_LOAD_2JCOORD_2ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128d _tmp1, _tmp2,_tmp3;                                            \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
}


#define GMX_MM_LOAD_2JCOORD_3ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128d _tmp1, _tmp2, _tmp3, _tmp4, _tmp5;                             \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
	 _tmp4          = _mm_loadu_pd(ptr1+6);                                 \
	 jy3            = _mm_loadu_pd(ptr2+6);                                 \
	 jz3            = _mm_load_sd(ptr1+8);                                  \
	 _tmp5          = _mm_load_sd(ptr2+8);                                  \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
     jx3            = _mm_unpacklo_pd(_tmp4,jy3);                           \
     jy3            = _mm_unpackhi_pd(_tmp4,jy3);                           \
     jz3            = _mm_unpacklo_pd(jz2,_tmp5);                           \
}


#define GMX_MM_LOAD_2JCOORD_4ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128d _tmp1, _tmp2,_tmp3, _tmp4, _tmp5, _tmp6;                       \
	 _tmp1          = _mm_loadu_pd(ptr1);                                   \
	 jy1            = _mm_loadu_pd(ptr2);                                   \
	 _tmp2          = _mm_loadu_pd(ptr1+2);                                 \
	 jx2            = _mm_loadu_pd(ptr2+2);                                 \
	 _tmp3          = _mm_loadu_pd(ptr1+4);                                 \
	 jz2            = _mm_loadu_pd(ptr2+4);                                 \
	 _tmp4          = _mm_loadu_pd(ptr1+6);                                 \
	 jy3            = _mm_loadu_pd(ptr2+6);                                 \
	 _tmp5          = _mm_loadu_pd(ptr1+8);                                 \
	 jx4            = _mm_loadu_pd(ptr2+8);                                 \
	 _tmp6          = _mm_loadu_pd(ptr1+10);                                \
	 jz4            = _mm_loadu_pd(ptr2+10);                                \
     jx1            = _mm_unpacklo_pd(_tmp1,jy1);                           \
     jy1            = _mm_unpackhi_pd(_tmp1,jy1);                           \
     jz1            = _mm_unpacklo_pd(_tmp2,jx2);                           \
     jx2            = _mm_unpackhi_pd(_tmp2,jx2);                           \
     jy2            = _mm_unpacklo_pd(_tmp3,jz2);                           \
     jz2            = _mm_unpackhi_pd(_tmp3,jz2);                           \
     jx3            = _mm_unpacklo_pd(_tmp4,jy3);                           \
     jy3            = _mm_unpackhi_pd(_tmp4,jy3);                           \
     jz3            = _mm_unpacklo_pd(_tmp5,jx4);                           \
     jx4            = _mm_unpackhi_pd(_tmp5,jx4);                           \
     jy4            = _mm_unpacklo_pd(_tmp6,jz4);                           \
     jz4            = _mm_unpackhi_pd(_tmp6,jz4);                           \
}


#define GMX_MM_UPDATE_1JCOORD_1ATOM_PD(ptr1,jx1,jy1,jz1) {            \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                       \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));      \
     _mm_store_sd(ptr1+2, _mm_add_sd( _mm_load_sd(ptr1+2), jz1));     \
}


#define GMX_MM_UPDATE_1JCOORD_2ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));        \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));    \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));    \
}


#define GMX_MM_UPDATE_1JCOORD_3ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     jx3            = _mm_unpacklo_pd(jx3,jy3);                         \
     _mm_storeu_pd(ptr1, _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));        \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));    \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));    \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), jx3 ));    \
     _mm_store_sd(ptr1+8, _mm_add_pd( _mm_load_sd(ptr1+8), jz3 ));      \
}


#define GMX_MM_UPDATE_1JCOORD_4ATOMS_PD(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     jx1            = _mm_unpacklo_pd(jx1,jy1);                         \
     jz1            = _mm_unpacklo_pd(jz1,jx2);                         \
     jy2            = _mm_unpacklo_pd(jy2,jz2);                         \
     jx3            = _mm_unpacklo_pd(jx3,jy3);                         \
     jz3            = _mm_unpacklo_pd(jz3,jx4);                         \
     jy4            = _mm_unpacklo_pd(jy4,jz4);                         \
     _mm_storeu_pd(ptr1,    _mm_add_pd( _mm_loadu_pd(ptr1), jx1 ));     \
     _mm_storeu_pd(ptr1+2,  _mm_add_pd( _mm_loadu_pd(ptr1+2), jz1 ));   \
     _mm_storeu_pd(ptr1+4,  _mm_add_pd( _mm_loadu_pd(ptr1+4), jy2 ));   \
     _mm_storeu_pd(ptr1+6,  _mm_add_pd( _mm_loadu_pd(ptr1+6), jx3 ));   \
     _mm_storeu_pd(ptr1+8,  _mm_add_pd( _mm_loadu_pd(ptr1+8), jz3 ));   \
     _mm_storeu_pd(ptr1+10, _mm_add_pd( _mm_loadu_pd(ptr1+10), jy4 ));  \
}


#define GMX_MM_UPDATE_2JCOORD_1ATOM_PD(ptr1,ptr2,jx1,jy1,jz1) {         \
     __m128d _tmp1,_tmp2,_tmp3;                                         \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp3          = _mm_unpackhi_pd(jz1,jz1);                         \
     _mm_storeu_pd(ptr1,  _mm_add_pd( _mm_loadu_pd(ptr1), _tmp1 ));     \
     _mm_storeu_pd(ptr2,  _mm_add_pd( _mm_loadu_pd(ptr2), _tmp2 ));     \
     _mm_store_sd(ptr1+2, _mm_add_pd( _mm_load_sd(ptr1+2), jz1 ));      \
     _mm_store_sd(ptr2+2, _mm_add_sd( _mm_load_sd(ptr2+2), _tmp3 ));    \
}


#define GMX_MM_UPDATE_2JCOORD_2ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128d _tmp1,_tmp2,_tmp3;                                         \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
}


#define GMX_MM_UPDATE_2JCOORD_3ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128d _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                             \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _tmp4          = _mm_unpacklo_pd(jx3,jy3);                         \
     jy3            = _mm_unpackhi_pd(jx3,jy3);                         \
     _tmp5          = _mm_unpackhi_pd(jz3,jz3);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), _tmp4 ));  \
     _mm_storeu_pd(ptr2+6, _mm_add_pd( _mm_loadu_pd(ptr2+6),   jy3 ));  \
     _mm_store_sd(ptr1+8, _mm_add_sd( _mm_load_sd(ptr1+8),   jz3 ));    \
     _mm_store_sd(ptr2+8, _mm_add_sd( _mm_load_sd(ptr2+8), _tmp5 ));    \
}


#define GMX_MM_UPDATE_2JCOORD_4ATOMS_PD(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128d _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6;                       \
     _tmp1          = _mm_unpacklo_pd(jx1,jy1);                         \
     jy1            = _mm_unpackhi_pd(jx1,jy1);                         \
     _tmp2          = _mm_unpacklo_pd(jz1,jx2);                         \
     jx2            = _mm_unpackhi_pd(jz1,jx2);                         \
     _tmp3          = _mm_unpacklo_pd(jy2,jz2);                         \
     jz2            = _mm_unpackhi_pd(jy2,jz2);                         \
     _tmp4          = _mm_unpacklo_pd(jx3,jy3);                         \
     jy3            = _mm_unpackhi_pd(jx3,jy3);                         \
     _tmp5          = _mm_unpacklo_pd(jz3,jx4);                         \
     jx4            = _mm_unpackhi_pd(jz3,jx4);                         \
     _tmp6          = _mm_unpacklo_pd(jy4,jz4);                         \
     jz4            = _mm_unpackhi_pd(jy4,jz4);                         \
     _mm_storeu_pd(ptr1,   _mm_add_pd( _mm_loadu_pd(ptr1),   _tmp1 ));  \
     _mm_storeu_pd(ptr2,   _mm_add_pd( _mm_loadu_pd(ptr2),     jy1 ));  \
     _mm_storeu_pd(ptr1+2, _mm_add_pd( _mm_loadu_pd(ptr1+2), _tmp2 ));  \
     _mm_storeu_pd(ptr2+2, _mm_add_pd( _mm_loadu_pd(ptr2+2),   jx2 ));  \
     _mm_storeu_pd(ptr1+4, _mm_add_pd( _mm_loadu_pd(ptr1+4), _tmp3 ));  \
     _mm_storeu_pd(ptr2+4, _mm_add_pd( _mm_loadu_pd(ptr2+4),   jz2 ));  \
     _mm_storeu_pd(ptr1+6, _mm_add_pd( _mm_loadu_pd(ptr1+6), _tmp4 ));  \
     _mm_storeu_pd(ptr2+6, _mm_add_pd( _mm_loadu_pd(ptr2+6),   jy3 ));  \
     _mm_storeu_pd(ptr1+8, _mm_add_pd( _mm_loadu_pd(ptr1+8), _tmp5 ));  \
     _mm_storeu_pd(ptr2+8, _mm_add_pd( _mm_loadu_pd(ptr2+8),   jx4 ));  \
     _mm_storeu_pd(ptr1+10,_mm_add_pd( _mm_loadu_pd(ptr1+10),_tmp6 ));  \
     _mm_storeu_pd(ptr2+10,_mm_add_pd( _mm_loadu_pd(ptr2+10),  jz4 ));  \
}


static inline __m128d
gmx_mm_scalarprod_pd(__m128d x, __m128d y, __m128d z)
{
	return (__m128d) _mm_add_pd(_mm_add_pd(_mm_mul_pd(x,x),_mm_mul_pd(y,y)),_mm_mul_pd(z,z));
}


static inline __m128d
gmx_mm_recip_pd(__m128d x)
{
	const __m128d two  = (const __m128d) {2.0,2.0};
	
	/* Lookup instruction only exists in single precision, convert back and forth... */
	__m128d lu = _mm_cvtps_pd(_mm_rcp_ps( _mm_cvtpd_ps(x)));
	
	/* Perform two N-R steps for double precision */
	lu         = _mm_mul_pd(lu,_mm_sub_pd(two,_mm_mul_pd(x,lu)));
	return (__m128d) _mm_mul_pd(lu,_mm_sub_pd(two,_mm_mul_pd(x,lu)));
}



static inline __m128d
gmx_mm_invsqrt_pd(__m128d x)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	const __m128d three = (const __m128d) {3.0,3.0};

	/* Lookup instruction only exists in single precision, convert back and forth... */
	__m128d lu = _mm_cvtps_pd(_mm_rsqrt_ps( _mm_cvtpd_ps(x)));
	
	lu = _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu,lu),x)),lu));
	return (__m128d) _mm_mul_pd(half,_mm_mul_pd(_mm_sub_pd(three,_mm_mul_pd(_mm_mul_pd(lu,lu),x)),lu));
}

static inline __m128d
gmx_mm_calc_rsq_pd(__m128d dx, __m128d dy, __m128d dz)
{
    return _mm_add_pd( _mm_add_pd( _mm_mul_pd(dx,dx), _mm_mul_pd(dy,dy) ), _mm_mul_pd(dz,dz) );
}



/* Routine to be called with rswitch/rcut at the beginning of a kernel
 * to set up the 7 constants used for analytic 5th order switch calculations.
 */
static inline void
gmx_mm_setup_switch5_constants_pd(__m128d rswitch, __m128d rcut,
								  __m128d *switch_C3, __m128d *switch_C4, __m128d *switch_C5,
								  __m128d *switch_D2, __m128d *switch_D3, __m128d *switch_D4)
{
	const __m128d  cm6  = (const __m128d) { -6.0, -6.0};
	const __m128d cm10  = (const __m128d) {-10.0,-10.0};
	const __m128d  c15  = (const __m128d) { 15.0, 15.0};
	const __m128d cm30  = (const __m128d) {-30.0,-30.0};
	const __m128d  c60  = (const __m128d) { 60.0, 60.0};

	__m128d d,dinv,dinv2,dinv3,dinv4,dinv5;
	
	d       = _mm_sub_pd(rcut,rswitch);
	dinv    = gmx_mm_recip_pd(d);
	dinv2   = _mm_mul_pd(dinv,dinv);
	dinv3   = _mm_mul_pd(dinv2,dinv);
	dinv4   = _mm_mul_pd(dinv2,dinv2);
	dinv5   = _mm_mul_pd(dinv3,dinv2);
	
	*switch_C3 = _mm_mul_pd(cm10,dinv3);
	*switch_C4 = _mm_mul_pd(c15,dinv4);
	*switch_C5 = _mm_mul_pd(cm6,dinv5);
	*switch_D2 = _mm_mul_pd(cm30,dinv3);
	*switch_D3 = _mm_mul_pd(c60,dinv4);
	*switch_D4 = _mm_mul_pd(cm30,dinv5);
}


#define GMX_MM_SET_SWITCH5_PD(r,rswitch,rcut,sw,dsw,sw_C3,sw_C4,sw_C5,sw_D2,sw_D3,sw_D4) {   \
    const __m128d  _sw_one  = (const __m128d) {1.0,1.0};                                     \
    __m128d d,d2;                                                                            \
    d     = _mm_max_pd(r,rswitch);                                                           \
    d     = _mm_min_pd(d,rcut);                                                              \
    d     = _mm_sub_pd(d,rswitch);                                                           \
    d2    = _mm_mul_pd(d,d);                                                                 \
    sw    = _mm_mul_pd(d,sw_C5);                                                             \
    dsw   = _mm_mul_pd(d,sw_D4);                                                             \
    sw    = _mm_add_pd(sw,sw_C4);                                                            \
    dsw   = _mm_add_pd(dsw,sw_D3);                                                           \
    sw    = _mm_mul_pd(sw,d);                                                                \
    dsw   = _mm_mul_pd(dsw,d);                                                               \
    sw    = _mm_add_pd(sw,sw_C3);                                                            \
    dsw   = _mm_add_pd(dsw,sw_D2);                                                           \
    sw    = _mm_mul_pd(sw,_mm_mul_pd(d,d2));                                                 \
    dsw   = _mm_mul_pd(dsw,d2);                                                              \
    sw    = _mm_add_pd(sw,_sw_one);                                                          \
}


/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_coulomb_pd(__m128d rinv, __m128d qq,__m128d *vctot)
{
	__m128d vcoul = _mm_mul_pd(qq,rinv);
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return (__m128d) vcoul;
}


static inline void
gmx_mm_int_coulomb_noforce_pd(__m128d rinv, __m128d qq,__m128d *vctot)
{
	__m128d vcoul = _mm_mul_pd(qq,rinv);
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return;
}

/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_coulombrf_pd(__m128d rinv, __m128d rsq, __m128d krf, __m128d crf, __m128d qq,__m128d *vctot)
{
	const __m128d two  = (const __m128d) {2.0,2.0};
	__m128d vcoul,krsq;
	
	krsq   = _mm_mul_pd(krf,rsq);
	vcoul  = _mm_mul_pd(qq, _mm_sub_pd(_mm_add_pd(rinv,krsq),crf));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	return (__m128d) _mm_mul_pd(qq, _mm_sub_pd(rinv, _mm_mul_pd(two,krsq)));
}


static inline void
gmx_mm_int_coulombrf_noforce_pd(__m128d rinv, __m128d rsq, __m128d krf, __m128d crf, __m128d qq,__m128d *vctot)
{
	__m128d vcoul,krsq;
	
	krsq   = _mm_mul_pd(krf,rsq);
	vcoul  = _mm_mul_pd(qq, _mm_sub_pd(_mm_add_pd(rinv,krsq),crf));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	return;
}



/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128d
gmx_mm_int_lj_pd(__m128d rinvsq, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
	const __m128d six    = (const __m128d) {6.0,6.0};
	const __m128d twelve = (const __m128d) {12.0,12.0};
	
	__m128d rinvsix,vvdw6,vvdw12;
		
	rinvsix  = _mm_mul_pd(_mm_mul_pd(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_pd(c6,rinvsix);  
	vvdw12   = _mm_mul_pd(c12, _mm_mul_pd(rinvsix,rinvsix));
	*vvdwtot = _mm_add_pd(*vvdwtot , _mm_sub_pd(vvdw12,vvdw6));
	
	return (__m128d) _mm_sub_pd( _mm_mul_pd(twelve,vvdw12),_mm_mul_pd(six,vvdw6));
}
		   

static inline void
gmx_mm_int_lj_potonly_pd(__m128d rinvsq, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
	__m128d rinvsix,vvdw6,vvdw12;
	
	rinvsix  = _mm_mul_pd(_mm_mul_pd(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_pd(c6,rinvsix);  
	vvdw12   = _mm_mul_pd(c12, _mm_mul_pd(rinvsix,rinvsix));
	*vvdwtot = _mm_add_pd(*vvdwtot , _mm_sub_pd(vvdw12,vvdw6));
	
	return;
}

#ifdef GMX_SSE4
#  define gmx_mm_extract_epi32(x, imm) _mm_extract_epi32(x,imm)
#else
#  define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#endif



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_coulomb_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d *vctot)
{
    __m128d  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i  n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_pd(VFtab + 4* n_a);
	F        = _mm_load_pd(VFtab + 4* n_b);
	G        = _mm_load_pd(VFtab + 4* n_a + 2);
	H        = _mm_load_pd(VFtab + 4* n_b + 2);
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	F        = _mm_add_pd(F, _mm_add_pd(G,H));  /* Fp    */
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Y, _mm_mul_pd(eps,F)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	F        = _mm_mul_pd(qq, _mm_add_pd(F, _mm_add_pd(G, _mm_add_pd(H,H))));
	
	return (__m128d) _mm_mul_pd(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_lj_pd(__m128d r, __m128d tabscale, double * VFtab, int offset, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a,n_b;
	
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset);
	Fd       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset);
	Gd       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 2);
	Hd       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 2);
	Yr       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 4);
	Fr       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 4);
	Gr       = _mm_load_pd(VFtab + 4*(offset+2)* n_a + 4*offset + 6);
	Hr       = _mm_load_pd(VFtab + 4*(offset+2)* n_b + 4*offset + 6);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fd        = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr        = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_2_table_coulomb_and_lj_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d c6, __m128d c12, 
								  __m128d *vctot, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	Yc       = _mm_load_pd(VFtab + 12* n_a);
	Fc       = _mm_load_pd(VFtab + 12* n_b);
	Gc       = _mm_load_pd(VFtab + 12* n_a + 2);
	Hc       = _mm_load_pd(VFtab + 12* n_b + 2);
	Yd       = _mm_load_pd(VFtab + 12* n_a + 4);
	Fd       = _mm_load_pd(VFtab + 12* n_b + 4);
	Gd       = _mm_load_pd(VFtab + 12* n_a + 6);
	Hd       = _mm_load_pd(VFtab + 12* n_b + 6);
	Yr       = _mm_load_pd(VFtab + 12* n_a + 8);
	Fr       = _mm_load_pd(VFtab + 12* n_b + 8);
	Gr       = _mm_load_pd(VFtab + 12* n_a + 10);
	Hr       = _mm_load_pd(VFtab + 12* n_b + 10);
	GMX_MM_TRANSPOSE2_PD(Yc,Fc);
	GMX_MM_TRANSPOSE2_PD(Gc,Hc);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	Hc       = _mm_mul_pd(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_pd(Gc,eps);               /* Geps  */
	Fc       = _mm_add_pd(Fc, _mm_add_pd(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Yc, _mm_mul_pd(eps,Fc)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fc       = _mm_mul_pd(qq, _mm_add_pd(Fc, _mm_add_pd(Gc, _mm_add_pd(Hc,Hc))));
	Fd       = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr       = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fc,_mm_add_pd(Fd,Fr)),tabscale);
}




/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_coulomb_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d *vctot)
{
    __m128d  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Y        = _mm_load_pd(VFtab + 4* n_a);
	F        = _mm_setzero_pd();
	G        = _mm_load_pd(VFtab + 4* n_a + 2);
	H        = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);

	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	F        = _mm_add_pd(F, _mm_add_pd(G,H));  /* Fp    */
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Y, _mm_mul_pd(eps,F)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	F        = _mm_mul_pd(qq, _mm_add_pd(F, _mm_add_pd(G, _mm_add_pd(H,H))));
	
	return (__m128d) _mm_mul_pd(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_lj_pd(__m128d r, __m128d tabscale, double * VFtab, int offset, __m128d c6, __m128d c12, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset);
	Fd       = _mm_setzero_pd();
	Gd       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 2);
	Hd       = _mm_setzero_pd();
	Yr       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 4);
	Fr       = _mm_setzero_pd();
	Gr       = _mm_load_pd(VFtab + 4*(offset+2)*n_a + 4*offset + 6);
	Hr       = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);
	
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fd        = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr        = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128d
gmx_mm_int_1_table_coulomb_and_lj_pd(__m128d r, __m128d tabscale, double * VFtab, __m128d qq, __m128d c6, __m128d c12, 
									 __m128d *vctot, __m128d *vvdwtot)
{
    __m128d  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i  n0;
	int     n_a;
	
    rt       = _mm_mul_pd(r,tabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Yc       = _mm_load_pd(VFtab + 12* n_a);
	Fc       = _mm_setzero_pd();
	Gc       = _mm_load_pd(VFtab + 12* n_a + 2);
	Hc       = _mm_setzero_pd();
	Yd       = _mm_load_pd(VFtab + 12* n_a + 4);
	Fd       = _mm_setzero_pd();
	Gd       = _mm_load_pd(VFtab + 12* n_a + 6);
	Hd       = _mm_setzero_pd();
	Yr       = _mm_load_pd(VFtab + 12* n_a + 8);
	Fr       = _mm_setzero_pd();
	Gr       = _mm_load_pd(VFtab + 12* n_a + 10);
	Hr       = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Yc,Fc);
	GMX_MM_TRANSPOSE2_PD(Gc,Hc);
	GMX_MM_TRANSPOSE2_PD(Yd,Fd);
	GMX_MM_TRANSPOSE2_PD(Gd,Hd);
	GMX_MM_TRANSPOSE2_PD(Yr,Fr);
	GMX_MM_TRANSPOSE2_PD(Gr,Hr);	
	Hc       = _mm_mul_pd(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_pd(Gc,eps);               /* Geps  */
	Fc       = _mm_add_pd(Fc, _mm_add_pd(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_pd(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_pd(Gd,eps);               /* Geps  */
	Fd       = _mm_add_pd(Fd, _mm_add_pd(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_pd(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_pd(Gr,eps);               /* Geps  */
	Fr       = _mm_add_pd(Fr, _mm_add_pd(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_pd(qq, _mm_add_pd(Yc, _mm_mul_pd(eps,Fc)));
	*vctot   = _mm_add_pd(*vctot,vcoul);
	
	vvdw6    = _mm_mul_pd(c6,  _mm_add_pd(Yd, _mm_mul_pd(eps,Fd)));
	vvdw12   = _mm_mul_pd(c12, _mm_add_pd(Yr, _mm_mul_pd(eps,Fr)));
	*vvdwtot = _mm_add_pd(*vvdwtot, _mm_add_pd(vvdw6,vvdw12));
	
	Fc       = _mm_mul_pd(qq, _mm_add_pd(Fc, _mm_add_pd(Gc, _mm_add_pd(Hc,Hc))));
	Fd       = _mm_mul_pd(c6,  _mm_add_pd(Fd, _mm_add_pd(Gd, _mm_add_pd(Hd,Hd))));
	Fr       = _mm_mul_pd(c12, _mm_add_pd(Fr, _mm_add_pd(Gr, _mm_add_pd(Hr,Hr))));
	
	return (__m128d) _mm_mul_pd( _mm_add_pd(Fc,_mm_add_pd(Fd,Fr)),tabscale);
}



/* Return force should be multiplied by +rinv to get fscal */
static inline __m128d
gmx_mm_int_2_genborn_pd(__m128d r, __m128d isai, 
						double * isaj1, double *isaj2, 
						__m128d gbtabscale, double * GBtab, __m128d qq, __m128d *dvdasum, 
						double *dvdaj1, double *dvdaj2,
						__m128d *vgbtot)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	
    __m128d  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i  n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_sd(isaj1);
	t2       = _mm_load_sd(isaj2);
	isaj     = _mm_unpacklo_pd(isaj,t2);  /* t2 t1 */
	
	isaprod     = _mm_mul_pd(isai,isaj);
	qq          = _mm_mul_pd(qq,isaprod);
	gbtabscale  = _mm_mul_pd( isaprod, gbtabscale );
	
	rt       = _mm_mul_pd(r,gbtabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_pd(GBtab + 4* n_a);
	F        = _mm_load_pd(GBtab + 4* n_b);
	G        = _mm_load_pd(GBtab + 4* n_a + 2);
	H        = _mm_load_pd(GBtab + 4* n_b + 2);
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);
	G        = _mm_mul_pd(G,eps);               /* Geps  */
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	F        = _mm_add_pd(_mm_add_pd(F,G),H);  /* Fp    */
	
	VV       = _mm_add_pd(Y, _mm_mul_pd(eps,F));
	FF       = _mm_add_pd(_mm_add_pd(F,G), _mm_add_pd(H,H));
	
	vgb      = _mm_mul_pd(qq, VV);
	*vgbtot  = _mm_sub_pd(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_pd(_mm_mul_pd(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_pd(half,(vgb+ftmp*r));
	
	*dvdasum = _mm_add_pd(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_pd(_mm_mul_pd(dvdatmp,isaj), isaj);
	
	/* Update 2 dada[j] values */
	Y        = _mm_load_sd(dvdaj1);
	F        = _mm_load_sd(dvdaj2);
	t2       = _mm_unpackhi_pd(dvdatmp,dvdatmp);
	
	_mm_store_sd( dvdaj1 , _mm_add_sd( Y, dvdatmp ) );
	_mm_store_sd( dvdaj2 , _mm_add_sd( F, t2 ) );
	
	return (__m128d) ftmp;
}

/* Return force should be multiplied by +rinv to get fscal */
static inline __m128d
gmx_mm_int_1_genborn_pd(__m128d r, __m128d isai, 
						double * isaj1, 
						__m128d gbtabscale, double * GBtab, __m128d qq, __m128d *dvdasum, 
						double *dvdaj1, 
						__m128d *vgbtot)
{
	const __m128d half  = (const __m128d) {0.5,0.5};
	
    __m128d  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i  n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_sd(isaj1);
	
	isaprod     = _mm_mul_pd(isai,isaj);
	qq          = _mm_mul_pd(qq,isaprod);
	gbtabscale  = _mm_mul_pd( isaprod, gbtabscale );
	
	rt       = _mm_mul_pd(r,gbtabscale); 
	n0       = _mm_cvttpd_epi32(rt);
	eps      = _mm_sub_pd(rt, _mm_cvtepi32_pd(n0));
	eps2     = _mm_mul_pd(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Y        = _mm_load_pd(GBtab + 4* n_a);
	F        = _mm_setzero_pd();
	G        = _mm_load_pd(GBtab + 4* n_a + 2);
	H        = _mm_setzero_pd();
	GMX_MM_TRANSPOSE2_PD(Y,F);
	GMX_MM_TRANSPOSE2_PD(G,H);

	G        = _mm_mul_pd(G,eps);               /* Geps  */
	H        = _mm_mul_pd(H,eps2);              /* Heps2 */
	F        = _mm_add_pd(_mm_add_pd(F,G),H);  /* Fp    */
	
	VV       = _mm_add_pd(Y, _mm_mul_pd(eps,F));
	FF       = _mm_add_pd(_mm_add_pd(F,G), _mm_add_pd(H,H));
	
	vgb      = _mm_mul_pd(qq, VV);
	*vgbtot  = _mm_sub_pd(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_pd(_mm_mul_pd(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_pd(half,(vgb+ftmp*r));
	
	*dvdasum = _mm_add_pd(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_pd(_mm_mul_pd(dvdatmp,isaj), isaj);
	
	/* Update 1 dada[j] values */
	Y        = _mm_load_sd(dvdaj1);
	
	_mm_store_sd( dvdaj1 , _mm_add_sd( Y, dvdatmp ) );
	
	return (__m128d) ftmp;
}




static inline void
gmx_mm_update_iforce_1atom_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							  double *fptr,
							  double *fshiftptr)
{
	__m128d t1,t2,t3;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	/* fiz1 is fine as it is */
#else
	/* SSE2 */
	/* transpose data */
	t1 = fix1;
	fix1 = _mm_unpacklo_pd(fix1,fiy1); /* y0 x0 */
	fiy1 = _mm_unpackhi_pd(t1,fiy1);   /* y1 x1 */
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_sd( fiz1, _mm_unpackhi_pd(fiz1,fiz1 ));
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_store_sd( fptr+2, _mm_add_sd( _mm_load_sd(fptr+2), fiz1 ));
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}

static inline void
gmx_mm_update_iforce_2atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	
	t1 = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); /* x and y sums */	
	fiz1 = _mm_add_sd(fiz1, _mm_unpackhi_pd(fiy2,fiy2)); /* z sum */

	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}



static inline void
gmx_mm_update_iforce_3atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   __m128d fix3, __m128d fiy3, __m128d fiz3,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1,t2;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
	fix3 = _mm_hadd_pd(fix3,fiy3);   
	/* fiz3 is fine as it is */
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	t1 = fix3;
	fix3 = _mm_unpacklo_pd(fix3,fiy3); /* y0 x0 */
	fiy3 = _mm_unpackhi_pd(t1,fiy3);   /* y1 x1 */
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
	
	fix3 = _mm_add_pd(fix3,fiy3);
	fiz3 = _mm_add_sd( fiz3, _mm_unpackhi_pd(fiz3,fiz3));
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	_mm_storeu_pd( fptr+6, _mm_add_pd( _mm_loadu_pd(fptr+6), fix3 ));
	_mm_store_sd( fptr+8, _mm_add_sd( _mm_load_sd(fptr+8), fiz3 ));
	
	fix1 = _mm_add_pd(fix1,fix3);
	t1   = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); /* x and y sums */	

	t2   = _mm_shuffle_pd(fiy2,fiy2,_MM_SHUFFLE2(1,1));
	fiz1 = _mm_add_sd(fiz1,fiz3);
	fiz1 = _mm_add_sd(fiz1,t2); /* z sum */
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}



static inline void
gmx_mm_update_iforce_4atoms_pd(__m128d fix1, __m128d fiy1, __m128d fiz1,
							   __m128d fix2, __m128d fiy2, __m128d fiz2,
							   __m128d fix3, __m128d fiy3, __m128d fiz3,
							   __m128d fix4, __m128d fiy4, __m128d fiz4,
							   double *fptr,
							   double *fshiftptr)
{
	__m128d t1,t2;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_pd(fix1,fiy1);   
	fiz1 = _mm_hadd_pd(fiz1,fix2);   
	fiy2 = _mm_hadd_pd(fiy2,fiz2);   
	fix3 = _mm_hadd_pd(fix3,fiy3);   
	fiz3 = _mm_hadd_pd(fiz3,fix4);   
	fiy4 = _mm_hadd_pd(fiy4,fiz4);   
#else
	/* SSE2 */
	/* transpose data */
	GMX_MM_TRANSPOSE2_PD(fix1,fiy1);
	GMX_MM_TRANSPOSE2_PD(fiz1,fix2);
	GMX_MM_TRANSPOSE2_PD(fiy2,fiz2);
	GMX_MM_TRANSPOSE2_PD(fix3,fiy3);
	GMX_MM_TRANSPOSE2_PD(fiz3,fix4);
	GMX_MM_TRANSPOSE2_PD(fiy4,fiz4);
	
	fix1 = _mm_add_pd(fix1,fiy1);
	fiz1 = _mm_add_pd(fiz1,fix2);
	fiy2 = _mm_add_pd(fiy2,fiz2);
	fix3 = _mm_add_pd(fix3,fiy3);
	fiz3 = _mm_add_pd(fiz3,fix4);
	fiy4 = _mm_add_pd(fiy4,fiz4);
#endif
	_mm_storeu_pd( fptr, _mm_add_pd( _mm_loadu_pd(fptr), fix1 ));
	_mm_storeu_pd( fptr+2, _mm_add_pd( _mm_loadu_pd(fptr+2), fiz1 ));
	_mm_storeu_pd( fptr+4, _mm_add_pd( _mm_loadu_pd(fptr+4), fiy2 ));
	_mm_storeu_pd( fptr+6, _mm_add_pd( _mm_loadu_pd(fptr+6), fix3 ));
	_mm_storeu_pd( fptr+8, _mm_add_pd( _mm_loadu_pd(fptr+8), fiz3 ));
	_mm_storeu_pd( fptr+10, _mm_add_pd( _mm_loadu_pd(fptr+10), fiy4 ));
	
	t1 = _mm_shuffle_pd(fiz1,fiy2,_MM_SHUFFLE2(0,1)); 
	fix1 = _mm_add_pd(fix1,t1); 
	t2 = _mm_shuffle_pd(fiz3,fiy4,_MM_SHUFFLE2(0,1)); 
	fix3 = _mm_add_pd(fix3,t2); 
	fix1 = _mm_add_pd(fix1,fix3); /* x and y sums */
	
	
	fiz1 = _mm_add_sd(fiz1, _mm_unpackhi_pd(fiy2,fiy2)); 
	fiz3 = _mm_add_sd(fiz3, _mm_unpackhi_pd(fiy4,fiy4)); 
	fiz1 = _mm_add_sd(fiz1,fiz3); /* z sum */
	
	
	_mm_storeu_pd( fshiftptr, _mm_add_pd( _mm_loadu_pd(fshiftptr), fix1 ));
	_mm_store_sd( fshiftptr+2, _mm_add_sd( _mm_load_sd(fshiftptr+2), fiz1 ));
}


static inline void
gmx_mm_update_1pot_pd(__m128d pot1, double *ptr1)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_pd(pot1,pot1);
#else
	/* SSE2 */
	pot1 = _mm_add_pd(pot1, _mm_unpackhi_pd(pot1,pot1));
#endif
	_mm_store_sd(ptr1,_mm_add_sd(pot1,_mm_load_sd(ptr1)));
}
				   
static inline void
gmx_mm_update_2pot_pd(__m128d pot1, double *ptr1, __m128d pot2, double *ptr2)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_pd(pot1,pot2); 
#else
	/* SSE2 */
	GMX_MM_TRANSPOSE2_PD(pot1,pot2);
	pot1 = _mm_add_pd(pot1,pot2);
#endif
	pot2 = _mm_unpackhi_pd(pot1,pot1);
	
	_mm_store_sd(ptr1,_mm_add_sd(pot1,_mm_load_sd(ptr1)));
	_mm_store_sd(ptr2,_mm_add_sd(pot2,_mm_load_sd(ptr2)));
}


