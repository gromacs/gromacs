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

#include <xmmintrin.h> /* SSE */
#ifdef GMX_SSE2
#  include <emmintrin.h> /* SSE2 */
#endif
#ifdef GMX_SSE3
#  include <pmmintrin.h> /* SSE3 */
#endif
#ifdef GMX_SSE4
#  include <smmintrin.h> /* SSE4.1 */
#endif


#include <stdio.h>
  


static inline void
gmx_printxmm_ps(const char *s,__m128 xmm)
{
	float f[4];
	
	_mm_storeu_ps(f,xmm);
	printf("%s: %15.10g %15.10g %15.10g %15.10g\n",s,f[0],f[1],f[2],f[3]);	
}



#define GMX_MM_LOAD_1JCOORD_1ATOM_PS(ptr1,jx1,jy1,jz1) {              \
	 jx1            = _mm_load_ss(ptr1);                               \
     jy1            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+1));\
     jz1            = _mm_shuffle_ps(jy1,jy1,_MM_SHUFFLE(0,0,0,1));    \
}


#define GMX_MM_LOAD_1JCOORD_2ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
	 jx1            = _mm_loadu_ps(ptr1);                              \
     jy1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,1));    \
     jz1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));    \
     jx2            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,3));    \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));\
     jz2            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,1));    \
}


#define GMX_MM_LOAD_1JCOORD_3ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
	 jx1            = _mm_loadu_ps(ptr1);                                          \
     jy1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,1));                \
     jz1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));                \
     jx2            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,3));                \
     jy2            = _mm_loadu_ps(ptr1+4);                                        \
     jz2            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,1));                \
     jx3            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,2));                \
     jy3            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,3));                \
     jz3            = _mm_load_ss(ptr1+8);                                         \
}


#define GMX_MM_LOAD_1JCOORD_4ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_ps(ptr1);                                                      \
     jy1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,1));                            \
     jz1            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));                            \
     jx2            = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,3));                            \
     jy2            = _mm_loadu_ps(ptr1+4);                                                    \
     jz2            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,1));                            \
     jx3            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,2));                            \
     jy3            = _mm_shuffle_ps(jy2,jy2,_MM_SHUFFLE(0,0,0,3));                            \
     jz3            = _mm_loadu_ps(ptr1+8);                                                    \
     jx4            = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));                            \
     jy4            = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,2));                            \
     jz4            = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,3));                            \
}


#define GMX_MM_LOAD_2JCOORD_1ATOM_PS(ptr1,ptr2,jx1,jy1,jz1) {  \
     __m128 _tmp1,_tmp2;                                        \
	 _tmp1           = _mm_load_ss(ptr1);                       \
     _tmp2           = _mm_load_ss(ptr2);                       \
	 _tmp1           = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));   \
     _tmp2           = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));   \
     jx1             = _mm_unpacklo_ps(_tmp1,_tmp2);            \
     jy1             = _mm_unpackhi_ps(_tmp1,_tmp2);            \
	 jx1             = _mm_unpacklo_ps(_tmp1,_tmp2);            \
     jz1             = _mm_movehl_ps(jz1,jy1);                  \
}


#define GMX_MM_LOAD_2JCOORD_2ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                   \
	 _tmp1          = _mm_loadu_ps(ptr1);                                   \
     jy1            = _mm_loadu_ps(ptr2);                                   \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));     \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));     \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                           \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                           \
     jy2            = _mm_unpacklo_ps(jy2,_tmp2);                           \
     jy1            = _mm_movehl_ps(jx1,jx1);                               \
     jx2            = _mm_movehl_ps(jz1,jz1);                               \
     jz2            = _mm_movehl_ps(jy2,jy2);                               \
}


#define GMX_MM_LOAD_2JCOORD_3ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                        \
	 _tmp1          = _mm_loadu_ps(ptr1);                                               \
     jy1            = _mm_loadu_ps(ptr2);                                               \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                             \
     jz2            = _mm_loadu_ps(ptr2+4);                                             \
     jz3            = _mm_load_ss(ptr1+8);                                              \
     _tmp3          = _mm_load_ss(ptr2+8);                                              \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                                       \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                                       \
     jy2            = _mm_unpacklo_ps(_tmp2,jz2);                                       \
     jx3            = _mm_unpackhi_ps(_tmp2,jz2);                                       \
     jy1            = _mm_movehl_ps(jx1,jx1);                                           \
     jx2            = _mm_movehl_ps(jz1,jz1);                                           \
     jz2            = _mm_movehl_ps(jy2,jy2);                                           \
     jy3            = _mm_movehl_ps(jx3,jx3);                                           \
     jz3            = _mm_unpacklo_ps(jz3,_tmp3);                                       \
}


#define GMX_MM_LOAD_2JCOORD_4ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1, _tmp2, _tmp3,_tmp4;                                                  \
	 _tmp1          = _mm_loadu_ps(ptr1);                                               \
     jy1            = _mm_loadu_ps(ptr2);                                               \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                             \
     jz2            = _mm_loadu_ps(ptr2+4);                                             \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                             \
     _tmp4          = _mm_loadu_ps(ptr2+8);                                             \
     jx1            = _mm_unpacklo_ps(_tmp1,jy1);                                       \
     jz1            = _mm_unpackhi_ps(_tmp1,jy1);                                       \
     jy2            = _mm_unpacklo_ps(_tmp2,jz2);                                       \
     jx3            = _mm_unpackhi_ps(_tmp2,jz2);                                       \
     jz3            = _mm_unpacklo_ps(_tmp3,_tmp4);                                     \
     jy4            = _mm_unpackhi_ps(_tmp3,_tmp4);                                     \
     jy1            = _mm_movehl_ps(jx1,jx1);                                           \
     jx2            = _mm_movehl_ps(jz1,jz1);                                           \
     jz2            = _mm_movehl_ps(jy2,jy2);                                           \
     jy3            = _mm_movehl_ps(jx3,jx3);                                           \
     jx4            = _mm_movehl_ps(jz3,jz3);                                           \
     jz4            = _mm_movehl_ps(jy4,jy4);                                           \
}


#define GMX_MM_LOAD_3JCOORD_1ATOM_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp3,_tmp4;                                      \
	 jx1            = _mm_load_ss(ptr1);                            \
     jy1            = _mm_load_ss(ptr2);                            \
     jz1            = _mm_load_ss(ptr3);                            \
	 jx1            = _mm_loadh_pi(jx1,(__m64 *)(ptr1+1));          \
     jy1            = _mm_loadh_pi(jy1,(__m64 *)(ptr2+1));          \
     jz1            = _mm_loadh_pi(jz1,(__m64 *)(ptr3+1));          \
     _tmp1          = _mm_unpacklo_ps(jx1,jy1);                     \
     _tmp3          = _mm_unpackhi_ps(jx1,jy1);                     \
     _tmp4          = _mm_unpackhi_ps(jz1,jz1);                     \
     jx1            = _mm_movelh_ps(_tmp1,jz1);                     \
     jy1            = _mm_movelh_ps(_tmp3,_tmp4);                   \
     jz1            = _mm_movehl_ps(_tmp4,_tmp3);                   \
}


#define GMX_MM_LOAD_3JCOORD_2ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                        \
	 jx1            = _mm_loadu_ps(ptr1);                                        \
     jy1            = _mm_loadu_ps(ptr2);                                        \
     jz1            = _mm_loadu_ps(ptr3);                                        \
     jx2            = _mm_setzero_ps();                                          \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                         \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));          \
     jz2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));          \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));          \
     _tmp2          = _mm_unpacklo_ps(jy2,jz2);                                  \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_tmp1);                              \
     jy2            = _mm_movelh_ps(_tmp2,_tmp1);                                \
     jz2            = _mm_movehl_ps(_tmp1,_tmp2);                                \
}


#define GMX_MM_LOAD_3JCOORD_3ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) {  \
     __m128 _tmp1, _tmp2;                                                                     \
	 jx1            = _mm_loadu_ps(ptr1);                                                     \
     jy1            = _mm_loadu_ps(ptr2);                                                     \
     jz1            = _mm_loadu_ps(ptr3);                                                     \
     jx2            = _mm_setzero_ps();                                                       \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                      \
     jy2            = _mm_loadu_ps(ptr1+4);                                                   \
     jz2            = _mm_loadu_ps(ptr2+4);                                                   \
     jx3            = _mm_loadu_ps(ptr3+4);                                                   \
     jy3            = _mm_setzero_ps();                                                       \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                      \
     jz3            = _mm_load_ss(ptr1+8);                                                    \
     _tmp1          = _mm_load_ss(ptr2+8);                                                    \
     _tmp2          = _mm_load_ss(ptr3+8);                                                    \
     jz3            = _mm_unpacklo_ps(jz3,_tmp2);                                             \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_tmp1);                                           \
     jz3            = _mm_unpacklo_ps(jz3,_tmp1);                                             \
}


#define GMX_MM_LOAD_3JCOORD_4ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_ps(ptr1);                                                     \
     jy1            = _mm_loadu_ps(ptr2);                                                     \
     jz1            = _mm_loadu_ps(ptr3);                                                     \
     jx2            = _mm_setzero_ps();                                                       \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                      \
     jy2            = _mm_loadu_ps(ptr1+4);                                                   \
     jz2            = _mm_loadu_ps(ptr2+4);                                                   \
     jx3            = _mm_loadu_ps(ptr3+4);                                                   \
     jy3            = _mm_setzero_ps();                                                       \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                      \
     jz3            = _mm_loadu_ps(ptr1+8);                                                   \
     jx4            = _mm_loadu_ps(ptr2+8);                                                   \
     jy4            = _mm_loadu_ps(ptr3+8);                                                   \
     jz4            = _mm_setzero_ps();                                                       \
     _MM_TRANSPOSE4_PS(jz3,jx4,jy4,jz4);                                                      \
}



#define GMX_MM_LOAD_4JCOORD_1ATOM_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                               \
	 jx1            = _mm_load_ss(ptr1);                                 \
     _tmp1          = _mm_load_ss(ptr2);                                 \
     jy1            = _mm_load_ss(ptr3);                                 \
     jz1            = _mm_load_ss(ptr4);                                 \
	 jx1            = _mm_loadh_pi(jx1,(__m64 *)(ptr1+1));               \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr2+1));             \
     jy1            = _mm_loadh_pi(jy1,(__m64 *)(ptr3+1));               \
     jz1            = _mm_loadh_pi(jz1,(__m64 *)(ptr4+1));               \
     _tmp2          = _mm_unpacklo_ps(jx1,_tmp1);                        \
     _tmp3          = _mm_unpacklo_ps(jy1,jz1);                          \
     _tmp4          = _mm_unpackhi_ps(jx1,_tmp1);                        \
     _tmp5          = _mm_unpackhi_ps(jy1,jz1);                          \
     jx1            = _mm_movelh_ps(_tmp2,_tmp3);                        \
     jy1            = _mm_movelh_ps(_tmp4,_tmp5);                        \
     jz1            = _mm_movehl_ps(_tmp5,_tmp4);                        \
}


#define GMX_MM_LOAD_4JCOORD_2ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                                             \
	 jx1            = _mm_loadu_ps(ptr1);                                             \
     jy1            = _mm_loadu_ps(ptr2);                                             \
     jz1            = _mm_loadu_ps(ptr3);                                             \
     jx2            = _mm_loadu_ps(ptr4);                                             \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                              \
     jy2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));               \
     jz2            = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr2+4));               \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));               \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr4+4));               \
     _tmp1          = _mm_unpacklo_ps(jy2,_tmp1);                                     \
     _tmp2          = _mm_unpacklo_ps(jz2,_tmp2);                                     \
     jy2            = _mm_unpacklo_ps(_tmp1,_tmp2);                                   \
     jz2            = _mm_unpackhi_ps(_tmp1,_tmp2);                                   \
}


#define GMX_MM_LOAD_4JCOORD_3ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                                              \
	 jx1            = _mm_loadu_ps(ptr1);                                                     \
     jy1            = _mm_loadu_ps(ptr2);                                                     \
     jz1            = _mm_loadu_ps(ptr3);                                                     \
     jx2            = _mm_loadu_ps(ptr4);                                                     \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                      \
     jy2            = _mm_loadu_ps(ptr1+4);                                                   \
     jz2            = _mm_loadu_ps(ptr2+4);                                                   \
     jx3            = _mm_loadu_ps(ptr3+4);                                                   \
     jy3            = _mm_loadu_ps(ptr4+4);                                                   \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                      \
     jz3            = _mm_load_ss(ptr1+8);                                                    \
     _tmp1          = _mm_load_ss(ptr2+8);                                                    \
     _tmp2          = _mm_load_ss(ptr3+8);                                                    \
     _tmp3          = _mm_load_ss(ptr4+8);                                                    \
     jz3            = _mm_unpacklo_ps(jz3,_tmp2);                                             \
     _tmp1          = _mm_unpacklo_ps(_tmp1,_tmp3);                                           \
     jz3            = _mm_unpacklo_ps(jz3,_tmp1);                                             \
}


#define GMX_MM_LOAD_4JCOORD_4ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
	 jx1            = _mm_loadu_ps(ptr1);                                                     \
     jy1            = _mm_loadu_ps(ptr2);                                                     \
     jz1            = _mm_loadu_ps(ptr3);                                                     \
     jx2            = _mm_loadu_ps(ptr4);                                                     \
     _MM_TRANSPOSE4_PS(jx1,jy1,jz1,jx2);                                                      \
     jy2            = _mm_loadu_ps(ptr1+4);                                                   \
     jz2            = _mm_loadu_ps(ptr2+4);                                                   \
     jx3            = _mm_loadu_ps(ptr3+4);                                                   \
     jy3            = _mm_loadu_ps(ptr4+4);                                                   \
     _MM_TRANSPOSE4_PS(jy2,jz2,jx3,jy3);                                                      \
     jz3            = _mm_loadu_ps(ptr1+8);                                                   \
     jx4            = _mm_loadu_ps(ptr2+8);                                                   \
     jy4            = _mm_loadu_ps(ptr3+8);                                                   \
     jz4            = _mm_loadu_ps(ptr4+8);                                                   \
     _MM_TRANSPOSE4_PS(jz3,jx4,jy4,jz4);                                                      \
}


#define GMX_MM_UPDATE_1JCOORD_1ATOM_PS(ptr1,jx1,jy1,jz1) {           \
     __m128 _tmp1;                                                    \
     jy1            = _mm_unpacklo_ps(jy1,jz1);                       \
     jx1            = _mm_movelh_ps(jx1,jy1);                         \
     _tmp1          = _mm_load_ss(ptr1);                              \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                          \
     _mm_store_ss(ptr1,_tmp1);                                        \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
}


#define GMX_MM_UPDATE_1JCOORD_2ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1, _tmp2;                                               \
     _tmp1          = _mm_loadu_ps(ptr1);                               \
     _tmp2          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4)); \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                         \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                         \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                         \
     jx1            = _mm_movelh_ps(jx1,jz1);                           \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                            \
     _tmp2          = _mm_add_ps(_tmp2,jy2);                            \
     _mm_storeu_ps(ptr1,_tmp1);                                         \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp2);                            \
}


#define GMX_MM_UPDATE_1JCOORD_3ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1, _tmp2, _tmp3;                                        \
     _tmp1          = _mm_loadu_ps(ptr1);                               \
     _tmp2          = _mm_loadu_ps(ptr1+4);                             \
     _tmp3          = _mm_load_ss(ptr1+8);                              \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                         \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                         \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                         \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                         \
     jx1            = _mm_movelh_ps(jx1,jz1);                           \
     jy2            = _mm_movelh_ps(jy2,jx3);                           \
     _tmp1           = _mm_add_ps(_tmp1,jx1);                           \
     _tmp2           = _mm_add_ps(_tmp2,jy2);                           \
     _tmp3           = _mm_add_ss(_tmp3,jz3);                           \
     _mm_storeu_ps(ptr1,_tmp1);                                         \
     _mm_storeu_ps(ptr1+4,_tmp2);                                       \
     _mm_store_ss(ptr1+8,_tmp3);                                        \
}


#define GMX_MM_UPDATE_1JCOORD_4ATOMS_PS(ptr1,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1, _tmp2, _tmp3;                                        \
     _tmp1          = _mm_loadu_ps(ptr1);                               \
     _tmp2          = _mm_loadu_ps(ptr1+4);                             \
     _tmp3          = _mm_loadu_ps(ptr1+8);                             \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                         \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                         \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                         \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                         \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                         \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                         \
     jx1            = _mm_movelh_ps(jx1,jz1);                           \
     jy2            = _mm_movelh_ps(jy2,jx3);                           \
     jz3            = _mm_movelh_ps(jz3,jy4);                           \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                            \
     _tmp2          = _mm_add_ps(_tmp2,jy2);                            \
     _tmp3          = _mm_add_ps(_tmp3,jz3);                            \
     _mm_storeu_ps(ptr1,_tmp1);                                         \
     _mm_storeu_ps(ptr1+4,_tmp2);                                       \
     _mm_storeu_ps(ptr1+8,_tmp3);                                       \
}


#define GMX_MM_UPDATE_2JCOORD_1ATOM_PS(ptr1,ptr2,jx1,jy1,jz1) {        \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4;                                    \
     _tmp1          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1));   \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr2));              \
     _tmp2          = _mm_load_ss(ptr1+2);                              \
     _tmp3          = _mm_load_ss(ptr2+2);                              \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                         \
     _tmp4          = _mm_shuffle_ps(jz1,jz1,_MM_SHUFFLE(0,0,0,1));     \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                            \
     _mm_storel_pi((__m64 *)(ptr1),_tmp1);                              \
     _mm_storeh_pi((__m64 *)(ptr2),_tmp1);                              \
     _mm_store_ss(ptr1+2,_mm_add_ss(_tmp2,jz1));                        \
	 _mm_store_ss(ptr2+2,_mm_add_ss(_tmp3,_tmp4));                      \
}


#define GMX_MM_UPDATE_2JCOORD_2ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2) {  \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5;                                     \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr2);                                      \
     _tmp3          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));        \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr2+4));                   \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     _tmp4          = _mm_movelh_ps(jx1,jz1);                                  \
     _tmp5          = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp1          = _mm_add_ps(_tmp1,_tmp4);                                 \
     _tmp2          = _mm_add_ps(_tmp2,_tmp5);                                 \
     _tmp3          = _mm_add_ps(_tmp3,jy2);                                   \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr2,_tmp2);                                                \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp3);                                   \
	 _mm_storeh_pi((__m64 *)(ptr2+4),_tmp3);                                   \
}


#define GMX_MM_UPDATE_2JCOORD_3ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;         \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
     _tmp3          = _mm_load_ss(ptr1+8);                                     \
     _tmp4          = _mm_loadu_ps(ptr2);                                      \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
     _tmp6          = _mm_load_ss(ptr2+8);                                     \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
     _tmp7          = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));            \
     _tmp8          = _mm_movelh_ps(jx1,jz1);                                  \
     _tmp9          = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp10         = _mm_movelh_ps(jy2,jx3);                                  \
     _tmp11         = _mm_movehl_ps(jx3,jy2);                                  \
     _tmp1          = _mm_add_ps(_tmp1,_tmp8);                                 \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                   \
     _tmp4          = _mm_add_ps(_tmp4,_tmp9);                                 \
     _tmp5          = _mm_add_ps(_tmp5,_tmp11);                                \
     _tmp6          = _mm_add_ss(_tmp6,_tmp7);                                 \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr1+4,_tmp2);                                              \
     _mm_store_ss(ptr1+8,_tmp3);                                               \
     _mm_storeu_ps(ptr2,_tmp4);                                                \
     _mm_storeu_ps(ptr2+4,_tmp5);                                              \
     _mm_store_ss(ptr2+8,_tmp6);                                               \
}


#define GMX_MM_UPDATE_2JCOORD_4ATOMS_PS(ptr1,ptr2,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13; \
     _tmp1          = _mm_loadu_ps(ptr1);                                       \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                     \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                     \
     _tmp4          = _mm_loadu_ps(ptr2);                                       \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                     \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                     \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                 \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                 \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                 \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                 \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                 \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                 \
     _tmp8          = _mm_movelh_ps(jx1,jz1);                                   \
     _tmp9          = _mm_movehl_ps(jz1,jx1);                                   \
     _tmp10         = _mm_movelh_ps(jy2,jx3);                                   \
     _tmp11         = _mm_movehl_ps(jx3,jy2);                                   \
     _tmp12         = _mm_movelh_ps(jz3,jy4);                                   \
     _tmp13         = _mm_movehl_ps(jy4,jz3);                                   \
     _tmp1          = _mm_add_ps(_tmp1,_tmp8);                                  \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                 \
     _tmp3          = _mm_add_ps(_tmp3,_tmp12);                                 \
     _tmp4          = _mm_add_ps(_tmp4,_tmp9);                                  \
     _tmp5          = _mm_add_ps(_tmp5,_tmp11);                                 \
     _tmp6          = _mm_add_ps(_tmp6,_tmp13);                                 \
     _mm_storeu_ps(ptr1,_tmp1);                                                 \
     _mm_storeu_ps(ptr1+4,_tmp2);                                               \
     _mm_storeu_ps(ptr1+8,_tmp3);                                               \
     _mm_storeu_ps(ptr2,_tmp4);                                                 \
     _mm_storeu_ps(ptr2+4,_tmp5);                                               \
     _mm_storeu_ps(ptr2+8,_tmp6);                                               \
}


#define GMX_MM_UPDATE_3JCOORD_1ATOM_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7;                \
     _tmp1          = _mm_load_ss(ptr1);                              \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
     _tmp2          = _mm_load_ss(ptr2);                              \
     _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));          \
     _tmp3          = _mm_load_ss(ptr3);                              \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));          \
     _tmp4          = _mm_unpacklo_ps(jy1,jz1);                       \
     _tmp5          = _mm_unpackhi_ps(jy1,jz1);                       \
     _tmp6          = _mm_shuffle_ps(jx1,_tmp4,_MM_SHUFFLE(3,2,0,1)); \
     _tmp7          = _mm_shuffle_ps(jx1,jx1,_MM_SHUFFLE(0,0,0,2));   \
     jx1            = _mm_movelh_ps(jx1,_tmp4);                       \
     _tmp7          = _mm_movelh_ps(_tmp7,_tmp5);                     \
     _tmp1          = _mm_add_ps(_tmp1,jx1);                          \
     _tmp2          = _mm_add_ps(_tmp2,_tmp6);                        \
     _tmp3          = _mm_add_ps(_tmp3,_tmp7);                        \
     _mm_store_ss(ptr1,_tmp1);                                        \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
     _mm_store_ss(ptr2,_tmp2);                                        \
     _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                          \
     _mm_store_ss(ptr3,_tmp3);                                        \
     _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                          \
}


#define GMX_MM_UPDATE_3JCOORD_2ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;      \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr2);                                      \
     _tmp3          = _mm_loadu_ps(ptr3);                                      \
     _tmp4          = _mm_loadl_pi(_tmp4,(__m64 *)(ptr1+4));                   \
     _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr2+4));                   \
     _tmp5          = _mm_loadl_pi(_tmp5,(__m64 *)(ptr3+4));                   \
     _tmp6          = _mm_unpackhi_ps(jx1,jy1);                                \
	 jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     _tmp7          = _mm_unpackhi_ps(jz1,jx2);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     _tmp8          = _mm_unpackhi_ps(jy2,jz2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     _tmp9          = _mm_movelh_ps(jx1,jz1);                                  \
     _tmp10         = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp6          = _mm_movelh_ps(_tmp6,_tmp7);                              \
     _tmp1          = _mm_add_ps(_tmp1,_tmp9);                                 \
     _tmp2          = _mm_add_ps(_tmp2,_tmp10);                                \
     _tmp3          = _mm_add_ps(_tmp3,_tmp6);                                 \
     _tmp4          = _mm_add_ps(_tmp4,jy2);                                   \
     _tmp5          = _mm_add_ps(_tmp5,_tmp8);                                 \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr2,_tmp2);                                                \
     _mm_storeu_ps(ptr3,_tmp3);                                                \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp4);                                   \
     _mm_storeh_pi((__m64 *)(ptr2+4),_tmp4);                                   \
	 _mm_storel_pi((__m64 *)(ptr3+4),_tmp5);                                   \
}


#define GMX_MM_UPDATE_3JCOORD_3ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;      \
     __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;    \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
     _tmp3          = _mm_load_ss(ptr1+8);                                     \
     _tmp4          = _mm_loadu_ps(ptr2);                                      \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
     _tmp6          = _mm_load_ss(ptr2+8);                                     \
     _tmp7          = _mm_loadu_ps(ptr3);                                      \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
     _tmp9          = _mm_load_ss(ptr3+8);                                     \
     _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
     _tmp14         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));            \
     _tmp15         = _mm_movehl_ps(jz3,jz3);                                  \
     _tmp16         = _mm_movelh_ps(jx1,jz1);                                  \
     _tmp17         = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                            \
     _tmp18         = _mm_movelh_ps(jy2,jx3);                                  \
     _tmp19         = _mm_movehl_ps(jx3,jy2);                                  \
     _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                            \
     _tmp1          = _mm_add_ps(_tmp1,_tmp16);                                \
     _tmp2          = _mm_add_ps(_tmp2,_tmp18);                                \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                   \
     _tmp4          = _mm_add_ps(_tmp4,_tmp17);                                \
     _tmp5          = _mm_add_ps(_tmp5,_tmp19);                                \
     _tmp6          = _mm_add_ss(_tmp6,_tmp14);                                \
     _tmp7          = _mm_add_ps(_tmp7,_tmp10);                                \
     _tmp8          = _mm_add_ps(_tmp8,_tmp12);                                \
     _tmp9          = _mm_add_ss(_tmp9,_tmp15);                                \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr1+4,_tmp2);                                              \
     _mm_store_ss(ptr1+8,_tmp3);                                               \
     _mm_storeu_ps(ptr2,_tmp4);                                                \
     _mm_storeu_ps(ptr2+4,_tmp5);                                              \
     _mm_store_ss(ptr2+8,_tmp6);                                               \
     _mm_storeu_ps(ptr3,_tmp7);                                                \
     _mm_storeu_ps(ptr3+4,_tmp8);                                              \
     _mm_store_ss(ptr3+8,_tmp9);                                               \
}


#define GMX_MM_UPDATE_3JCOORD_4ATOMS_PS(ptr1,ptr2,ptr3,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;         \
     __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21;      \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                    \
     _tmp4          = _mm_loadu_ps(ptr2);                                      \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                    \
     _tmp7          = _mm_loadu_ps(ptr3);                                      \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
     _tmp9          = _mm_loadu_ps(ptr3+8);                                    \
     _tmp10         = _mm_unpackhi_ps(jx1,jy1);                                \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     _tmp11         = _mm_unpackhi_ps(jz1,jx2);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     _tmp12         = _mm_unpackhi_ps(jy2,jz2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     _tmp13         = _mm_unpackhi_ps(jx3,jy3);                                \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
     _tmp14         = _mm_unpackhi_ps(jz3,jx4);                                \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                \
     _tmp15         = _mm_unpackhi_ps(jy4,jz4);                                \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                \
     _tmp16         = _mm_movelh_ps(jx1,jz1);                                  \
     _tmp17         = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp10         = _mm_movelh_ps(_tmp10,_tmp11);                            \
     _tmp18         = _mm_movelh_ps(jy2,jx3);                                  \
     _tmp19         = _mm_movehl_ps(jx3,jy2);                                  \
     _tmp12         = _mm_movelh_ps(_tmp12,_tmp13);                            \
     _tmp20         = _mm_movelh_ps(jz3,jy4);                                  \
     _tmp21         = _mm_movehl_ps(jy4,jz3);                                  \
     _tmp14         = _mm_movelh_ps(_tmp14,_tmp15);                            \
     _tmp1          = _mm_add_ps(_tmp1,_tmp16);                                \
     _tmp2          = _mm_add_ps(_tmp2,_tmp18);                                \
     _tmp3          = _mm_add_ps(_tmp3,_tmp20);                                \
     _tmp4          = _mm_add_ps(_tmp4,_tmp17);                                \
     _tmp5          = _mm_add_ps(_tmp5,_tmp19);                                \
     _tmp6          = _mm_add_ps(_tmp6,_tmp21);                                \
     _tmp7          = _mm_add_ps(_tmp7,_tmp10);                                \
     _tmp8          = _mm_add_ps(_tmp8,_tmp12);                                \
     _tmp9          = _mm_add_ps(_tmp9,_tmp14);                                \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr1+4,_tmp2);                                              \
     _mm_storeu_ps(ptr1+8,_tmp3);                                              \
     _mm_storeu_ps(ptr2,_tmp4);                                                \
     _mm_storeu_ps(ptr2+4,_tmp5);                                              \
     _mm_storeu_ps(ptr2+8,_tmp6);                                              \
     _mm_storeu_ps(ptr3,_tmp7);                                                \
     _mm_storeu_ps(ptr3+4,_tmp8);                                              \
     _mm_storeu_ps(ptr3+8,_tmp9);                                              \
}



#define GMX_MM_UPDATE_4JCOORD_1ATOM_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;      \
     _tmp1          = _mm_load_ss(ptr1);                              \
     _tmp1          = _mm_loadh_pi(_tmp1,(__m64 *)(ptr1+1));          \
     _tmp2          = _mm_load_ss(ptr2);                              \
     _tmp2          = _mm_loadh_pi(_tmp2,(__m64 *)(ptr2+1));          \
     _tmp3          = _mm_load_ss(ptr3);                              \
     _tmp3          = _mm_loadh_pi(_tmp3,(__m64 *)(ptr3+1));          \
     _tmp4          = _mm_load_ss(ptr4);                              \
     _tmp4          = _mm_loadh_pi(_tmp4,(__m64 *)(ptr4+1));          \
     _tmp5          = _mm_unpacklo_ps(jy1,jz1);                       \
     _tmp6          = _mm_unpackhi_ps(jy1,jz1);                       \
     _tmp7          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(1,0,0,0)); \
     _tmp8          = _mm_shuffle_ps(jx1,_tmp5,_MM_SHUFFLE(3,2,0,1)); \
     _tmp9          = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(1,0,0,2)); \
     _tmp10         = _mm_shuffle_ps(jx1,_tmp6,_MM_SHUFFLE(3,2,0,3)); \
     _tmp1          = _mm_add_ps(_tmp1,_tmp7);                        \
     _tmp2          = _mm_add_ps(_tmp2,_tmp8);                        \
     _tmp3          = _mm_add_ps(_tmp3,_tmp9);                        \
     _tmp4          = _mm_add_ps(_tmp4,_tmp10);                       \
     _mm_store_ss(ptr1,_tmp1);                                        \
     _mm_storeh_pi((__m64 *)(ptr1+1),_tmp1);                          \
     _mm_store_ss(ptr2,_tmp2);                                        \
     _mm_storeh_pi((__m64 *)(ptr2+1),_tmp2);                          \
     _mm_store_ss(ptr3,_tmp3);                                        \
     _mm_storeh_pi((__m64 *)(ptr3+1),_tmp3);                          \
     _mm_store_ss(ptr4,_tmp4);                                        \
     _mm_storeh_pi((__m64 *)(ptr4+1),_tmp4);                          \
}


#define GMX_MM_UPDATE_4JCOORD_2ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11,_tmp12,_tmp13;   \
     _tmp1          = _mm_loadu_ps(ptr1);                                       \
     _tmp2          = _mm_loadu_ps(ptr2);                                       \
     _tmp3          = _mm_loadu_ps(ptr3);                                       \
     _tmp4          = _mm_loadu_ps(ptr4);                                       \
     _tmp5          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr1+4));         \
     _tmp5          = _mm_loadh_pi(_tmp5,(__m64 *)(ptr2+4));                    \
     _tmp6          = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)(ptr3+4));         \
     _tmp6          = _mm_loadh_pi(_tmp6,(__m64 *)(ptr4+4));                    \
     _tmp7          = _mm_unpackhi_ps(jx1,jy1);                                 \
	 jx1            = _mm_unpacklo_ps(jx1,jy1);                                 \
     _tmp8          = _mm_unpackhi_ps(jz1,jx2);                                 \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                 \
     _tmp9          = _mm_unpackhi_ps(jy2,jz2);                                 \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                 \
     _tmp10         = _mm_movelh_ps(jx1,jz1);                                   \
     _tmp11         = _mm_movehl_ps(jz1,jx1);                                   \
     _tmp12         = _mm_movelh_ps(_tmp7,_tmp8);                               \
     _tmp13         = _mm_movehl_ps(_tmp8,_tmp7);                               \
     _tmp1          = _mm_add_ps(_tmp1,_tmp10);                                 \
     _tmp2          = _mm_add_ps(_tmp2,_tmp11);                                 \
     _tmp3          = _mm_add_ps(_tmp3,_tmp12);                                 \
     _tmp4          = _mm_add_ps(_tmp4,_tmp13);                                 \
     _tmp5          = _mm_add_ps(_tmp5,jy2);                                    \
     _tmp6          = _mm_add_ps(_tmp6,_tmp9);                                  \
     _mm_storeu_ps(ptr1,_tmp1);                                                 \
     _mm_storeu_ps(ptr2,_tmp2);                                                 \
     _mm_storeu_ps(ptr3,_tmp3);                                                 \
     _mm_storeu_ps(ptr4,_tmp4);                                                 \
     _mm_storel_pi((__m64 *)(ptr1+4),_tmp5);                                    \
     _mm_storeh_pi((__m64 *)(ptr2+4),_tmp5);                                    \
	 _mm_storel_pi((__m64 *)(ptr3+4),_tmp6);                                    \
	 _mm_storeh_pi((__m64 *)(ptr4+4),_tmp6);                                    \
}


#define GMX_MM_UPDATE_4JCOORD_3ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10;       \
     __m128 _tmp11,_tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19;     \
     __m128 _tmp20,_tmp21,_tmp22,_tmp23,_tmp24,_tmp25;                          \
     _tmp1          = _mm_loadu_ps(ptr1);                                       \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                     \
     _tmp3          = _mm_load_ss(ptr1+8);                                      \
     _tmp4          = _mm_loadu_ps(ptr2);                                       \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                     \
     _tmp6          = _mm_load_ss(ptr2+8);                                      \
     _tmp7          = _mm_loadu_ps(ptr3);                                       \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                     \
     _tmp9          = _mm_load_ss(ptr3+8);                                      \
     _tmp10         = _mm_loadu_ps(ptr4);                                       \
     _tmp11         = _mm_loadu_ps(ptr4+4);                                     \
     _tmp12         = _mm_load_ss(ptr4+8);                                      \
     _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                 \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                 \
     _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                 \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                 \
     _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                 \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                 \
     _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                 \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                 \
     _tmp17         = _mm_shuffle_ps(jz3,jz3,_MM_SHUFFLE(0,0,0,1));             \
     _tmp18         = _mm_movehl_ps(jz3,jz3);                                   \
     _tmp19         = _mm_shuffle_ps(_tmp18,_tmp18,_MM_SHUFFLE(0,0,0,1));       \
     _tmp20         = _mm_movelh_ps(jx1,jz1);                                   \
     _tmp21         = _mm_movehl_ps(jz1,jx1);                                   \
     _tmp22         = _mm_movelh_ps(_tmp13,_tmp14);                             \
     _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                             \
     _tmp23         = _mm_movelh_ps(jy2,jx3);                                   \
     _tmp24         = _mm_movehl_ps(jx3,jy2);                                   \
     _tmp25         = _mm_movelh_ps(_tmp15,_tmp16);                             \
     _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                             \
     _tmp1          = _mm_add_ps(_tmp1,_tmp20);                                 \
     _tmp2          = _mm_add_ps(_tmp2,_tmp23);                                 \
     _tmp3          = _mm_add_ss(_tmp3,jz3);                                    \
     _tmp4          = _mm_add_ps(_tmp4,_tmp21);                                 \
     _tmp5          = _mm_add_ps(_tmp5,_tmp24);                                 \
     _tmp6          = _mm_add_ss(_tmp6,_tmp17);                                 \
     _tmp7          = _mm_add_ps(_tmp7,_tmp22);                                 \
     _tmp8          = _mm_add_ps(_tmp8,_tmp25);                                 \
     _tmp9          = _mm_add_ss(_tmp9,_tmp18);                                 \
     _tmp10         = _mm_add_ps(_tmp10,_tmp14);                                \
     _tmp11         = _mm_add_ps(_tmp11,_tmp16);                                \
     _tmp12         = _mm_add_ss(_tmp12,_tmp19);                                \
     _mm_storeu_ps(ptr1,_tmp1);                                                 \
     _mm_storeu_ps(ptr1+4,_tmp2);                                               \
     _mm_store_ss(ptr1+8,_tmp3);                                                \
     _mm_storeu_ps(ptr2,_tmp4);                                                 \
     _mm_storeu_ps(ptr2+4,_tmp5);                                               \
     _mm_store_ss(ptr2+8,_tmp6);                                                \
     _mm_storeu_ps(ptr3,_tmp7);                                                 \
     _mm_storeu_ps(ptr3+4,_tmp8);                                               \
     _mm_store_ss(ptr3+8,_tmp9);                                                \
     _mm_storeu_ps(ptr4,_tmp10);                                                \
     _mm_storeu_ps(ptr4+4,_tmp11);                                              \
     _mm_store_ss(ptr4+8,_tmp12);                                               \
}


#define GMX_MM_UPDATE_4JCOORD_4ATOMS_PS(ptr1,ptr2,ptr3,ptr4,jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4) { \
     __m128 _tmp1,_tmp2,_tmp3,_tmp4,_tmp5,_tmp6,_tmp7,_tmp8,_tmp9,_tmp10,_tmp11;         \
     __m128 _tmp12,_tmp13,_tmp14,_tmp15,_tmp16,_tmp17,_tmp18,_tmp19,_tmp20,_tmp21,_tmp22;\
     __m128 _tmp23,_tmp24;                                                     \
     _tmp1          = _mm_loadu_ps(ptr1);                                      \
     _tmp2          = _mm_loadu_ps(ptr1+4);                                    \
     _tmp3          = _mm_loadu_ps(ptr1+8);                                    \
     _tmp4          = _mm_loadu_ps(ptr2);                                      \
     _tmp5          = _mm_loadu_ps(ptr2+4);                                    \
     _tmp6          = _mm_loadu_ps(ptr2+8);                                    \
     _tmp7          = _mm_loadu_ps(ptr3);                                      \
     _tmp8          = _mm_loadu_ps(ptr3+4);                                    \
     _tmp9          = _mm_loadu_ps(ptr3+8);                                    \
     _tmp10         = _mm_loadu_ps(ptr4);                                      \
     _tmp11         = _mm_loadu_ps(ptr4+4);                                    \
     _tmp12         = _mm_loadu_ps(ptr4+8);                                    \
     _tmp13         = _mm_unpackhi_ps(jx1,jy1);                                \
     jx1            = _mm_unpacklo_ps(jx1,jy1);                                \
     _tmp14         = _mm_unpackhi_ps(jz1,jx2);                                \
     jz1            = _mm_unpacklo_ps(jz1,jx2);                                \
     _tmp15         = _mm_unpackhi_ps(jy2,jz2);                                \
     jy2            = _mm_unpacklo_ps(jy2,jz2);                                \
     _tmp16         = _mm_unpackhi_ps(jx3,jy3);                                \
     jx3            = _mm_unpacklo_ps(jx3,jy3);                                \
     _tmp17         = _mm_unpackhi_ps(jz3,jx4);                                \
     jz3            = _mm_unpacklo_ps(jz3,jx4);                                \
     _tmp18         = _mm_unpackhi_ps(jy4,jz4);                                \
     jy4            = _mm_unpacklo_ps(jy4,jz4);                                \
     _tmp19         = _mm_movelh_ps(jx1,jz1);                                  \
     jz1            = _mm_movehl_ps(jz1,jx1);                                  \
     _tmp20         = _mm_movelh_ps(_tmp13,_tmp14);                            \
     _tmp14         = _mm_movehl_ps(_tmp14,_tmp13);                            \
     _tmp21         = _mm_movelh_ps(jy2,jx3);                                  \
     jx3            = _mm_movehl_ps(jx3,jy2);                                  \
     _tmp22         = _mm_movelh_ps(_tmp15,_tmp16);                            \
     _tmp16         = _mm_movehl_ps(_tmp16,_tmp15);                            \
     _tmp23         = _mm_movelh_ps(jz3,jy4);                                  \
     jy4            = _mm_movehl_ps(jy4,jz3);                                  \
     _tmp24         = _mm_movelh_ps(_tmp17,_tmp18);                            \
     _tmp18         = _mm_movehl_ps(_tmp18,_tmp17);                            \
     _tmp1          = _mm_add_ps(_tmp1,_tmp19);                                \
     _tmp2          = _mm_add_ps(_tmp2,_tmp21);                                \
     _tmp3          = _mm_add_ps(_tmp3,_tmp23);                                \
     _tmp4          = _mm_add_ps(_tmp4,jz1);                                   \
     _tmp5          = _mm_add_ps(_tmp5,jx3);                                   \
     _tmp6          = _mm_add_ps(_tmp6,jy4);                                   \
     _tmp7          = _mm_add_ps(_tmp7,_tmp20);                                \
     _tmp8          = _mm_add_ps(_tmp8,_tmp22);                                \
     _tmp9          = _mm_add_ps(_tmp9,_tmp24);                                \
     _tmp10         = _mm_add_ps(_tmp10,_tmp14);                               \
     _tmp11         = _mm_add_ps(_tmp11,_tmp16);                               \
     _tmp12         = _mm_add_ps(_tmp12,_tmp18);                               \
     _mm_storeu_ps(ptr1,_tmp1);                                                \
     _mm_storeu_ps(ptr1+4,_tmp2);                                              \
     _mm_storeu_ps(ptr1+8,_tmp3);                                              \
     _mm_storeu_ps(ptr2,_tmp4);                                                \
     _mm_storeu_ps(ptr2+4,_tmp5);                                              \
     _mm_storeu_ps(ptr2+8,_tmp6);                                              \
     _mm_storeu_ps(ptr3,_tmp7);                                                \
     _mm_storeu_ps(ptr3+4,_tmp8);                                              \
     _mm_storeu_ps(ptr3+8,_tmp9);                                              \
     _mm_storeu_ps(ptr4,_tmp10);                                               \
     _mm_storeu_ps(ptr4+4,_tmp11);                                             \
     _mm_storeu_ps(ptr4+8,_tmp12);                                             \
}



static inline __m128
gmx_mm_scalarprod_ps(__m128 x, __m128 y, __m128 z)
{
	return _mm_add_ps(_mm_add_ps(_mm_mul_ps(x,x),_mm_mul_ps(y,y)),_mm_mul_ps(z,z));
}


static inline __m128
gmx_mm_recip_ps(__m128 x)
{
	const __m128 two  = {2.0,2.0,2.0,2.0};
	
	__m128 lu = _mm_rcp_ps(x);
	return _mm_mul_ps(lu,_mm_sub_ps(two,_mm_mul_ps(x,lu)));
}



static inline __m128
gmx_mm_invsqrt_ps(__m128 x)
{
	const __m128 half  = {0.5,0.5,0.5,0.5};
	const __m128 three = {3.0,3.0,3.0,3.0};

	__m128 lu = _mm_rsqrt_ps(x);
	
	return _mm_mul_ps(half,_mm_mul_ps(_mm_sub_ps(three,_mm_mul_ps(_mm_mul_ps(lu,lu),x)),lu));
}


/* Routine to be called with rswitch/rcut at the beginning of a kernel
 * to set up the 7 constants used for analytic 5th order switch calculations.
 */
static inline void
gmx_mm_setup_switch5_constants_ps(__m128 rswitch, __m128 rcut,
								  __m128 *switch_C3, __m128 *switch_C4, __m128 *switch_C5,
								  __m128 *switch_D2, __m128 *switch_D3, __m128 *switch_D4)
{
	const __m128  cm6  = { -6.0, -6.0, -6.0, -6.0};
	const __m128 cm10  = {-10.0,-10.0,-10.0,-10.0};
	const __m128  c15  = { 15.0, 15.0, 15.0, 15.0};
	const __m128 cm30  = {-30.0,-30.0,-30.0,-30.0};
	const __m128  c60  = { 60.0, 60.0, 60.0, 60.0};

	__m128 d,dinv,dinv2,dinv3,dinv4,dinv5;
	
	d       = _mm_sub_ps(rcut,rswitch);
	dinv    = gmx_mm_recip_ps(d);
	dinv2   = _mm_mul_ps(dinv,dinv);
	dinv3   = _mm_mul_ps(dinv2,dinv);
	dinv4   = _mm_mul_ps(dinv2,dinv2);
	dinv5   = _mm_mul_ps(dinv3,dinv2);
	
	*switch_C3 = _mm_mul_ps(cm10,dinv3);
	*switch_C4 = _mm_mul_ps(c15,dinv4);
	*switch_C5 = _mm_mul_ps(cm6,dinv5);
	*switch_D2 = _mm_mul_ps(cm30,dinv3);
	*switch_D3 = _mm_mul_ps(c60,dinv4);
	*switch_D4 = _mm_mul_ps(cm30,dinv5);
}


#define GMX_MM_SET_SWITCH5_PS(r,rswitch,rcut,sw,dsw,sw_C3,sw_C4,sw_C5,sw_D2,sw_D3,sw_D4) {   \
    const __m128  _sw_one  = {  1.0,  1.0,  1.0,  1.0};                                      \
    __m128 d,d2;                                                                             \
    d     = _mm_max_ps(r,rswitch);                                                           \
    d     = _mm_min_ps(d,rcut);                                                              \
    d     = _mm_sub_ps(d,rswitch);                                                           \
    d2    = _mm_mul_ps(d,d);                                                                 \
    sw    = _mm_mul_ps(d,sw_C5);                                                             \
    dsw   = _mm_mul_ps(d,sw_D4);                                                             \
    sw    = _mm_add_ps(sw,sw_C4);                                                            \
    dsw   = _mm_add_ps(dsw,sw_D3);                                                           \
    sw    = _mm_mul_ps(sw,d);                                                                \
    dsw   = _mm_mul_ps(dsw,d);                                                               \
    sw    = _mm_add_ps(sw,sw_C3);                                                            \
    dsw   = _mm_add_ps(dsw,sw_D2);                                                           \
    sw    = _mm_mul_ps(sw,_mm_mul_ps(d,d2));                                                 \
    dsw   = _mm_mul_ps(dsw,d2);                                                              \
    sw    = _mm_add_ps(sw,_sw_one);                                                          \
}


/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128
gmx_mm_int_coulomb_ps(__m128 rinv, __m128 qq,__m128 *vctot)
{
	__m128 vcoul = _mm_mul_ps(qq,rinv);
	*vctot   = _mm_add_ps(*vctot,vcoul);
	return vcoul;
}


static inline void
gmx_mm_int_coulomb_noforce_ps(__m128 rinv, __m128 qq,__m128 *vctot)
{
	__m128 vcoul = _mm_mul_ps(qq,rinv);
	*vctot   = _mm_add_ps(*vctot,vcoul);
	return;
}

/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128
gmx_mm_int_coulombrf_ps(const __m128 rinv, const __m128 rsq, const __m128 krf, const __m128 crf, const __m128 qq,__m128 *vctot)
{
	const __m128 two  = {2.0,2.0,2.0,2.0};
	__m128 vcoul,krsq;
	
	krsq   = _mm_mul_ps(krf,rsq);
	vcoul  = _mm_mul_ps(qq, _mm_sub_ps(_mm_add_ps(rinv,krsq),crf));
	*vctot = _mm_add_ps(*vctot,vcoul);
	
	return _mm_mul_ps(qq, _mm_sub_ps(rinv, _mm_mul_ps(two,krsq)));
}


static inline void
gmx_mm_int_coulombrf_noforce_ps(__m128 rinv, __m128 rsq, __m128 krf, __m128 crf, __m128 qq,__m128 *vctot)
{
	__m128 vcoul,krsq;
	
	krsq   = _mm_mul_ps(krf,rsq);
	vcoul  = _mm_mul_ps(qq, _mm_sub_ps(_mm_add_ps(rinv,krsq),crf));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	return;
}


/* GB */




/* GB + RF */


/* Returns fscaltmp, multiply with rinvsq to get fscal! */
static inline __m128
gmx_mm_int_lj_ps(__m128 rinvsq, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
	const __m128 six    = {6.0,6.0,6.0,6.0};
	const __m128 twelve = {12.0,12.0,12.0,12.0};
	
	__m128 rinvsix,vvdw6,vvdw12;
		
	rinvsix  = _mm_mul_ps(_mm_mul_ps(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_ps(c6,rinvsix);  
	vvdw12   = _mm_mul_ps(c12, _mm_mul_ps(rinvsix,rinvsix));
	*vvdwtot = _mm_add_ps(*vvdwtot , _mm_sub_ps(vvdw12,vvdw6));
	
	return _mm_sub_ps( _mm_mul_ps(twelve,vvdw12),_mm_mul_ps(six,vvdw6));
}
		   

static inline void
gmx_mm_int_lj_potonly_ps(__m128 rinvsq, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
	__m128 rinvsix,vvdw6,vvdw12;
	
	rinvsix  = _mm_mul_ps(_mm_mul_ps(rinvsq,rinvsq),rinvsq);
	vvdw6    = _mm_mul_ps(c6,rinvsix);  
	vvdw12   = _mm_mul_ps(c12, _mm_mul_ps(rinvsix,rinvsix));
	*vvdwtot = _mm_add_ps(*vvdwtot , _mm_sub_ps(vvdw12,vvdw6));
	
	return;
}

#ifdef GMX_SSE4
#  define gmx_mm_extract_epi32(x, imm) _mm_extract_epi32(x,imm)
#else
#  define gmx_mm_extract_epi32(x, imm) _mm_cvtsi128_si32(_mm_srli_si128((x), 4 * (imm)))
#endif


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_4_table_coulomb_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 *vctot)
{
    __m128  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);

	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	n_d      = gmx_mm_extract_epi32(n0,3);
	Y        = _mm_load_ps(VFtab + 4* n_a);
	F        = _mm_load_ps(VFtab + 4* n_b);
	G        = _mm_load_ps(VFtab + 4* n_c);
	H        = _mm_load_ps(VFtab + 4* n_d);
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	F        = _mm_add_ps(F, _mm_add_ps(G,H));  /* Fp    */
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Y, _mm_mul_ps(eps,F)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	F        = _mm_mul_ps(qq, _mm_add_ps(F, _mm_add_ps(G, _mm_add_ps(H,H))));
	
	return _mm_mul_ps(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_4_table_lj_ps(__m128 r, __m128 tabscale, float * VFtab, int offset, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	n_d      = gmx_mm_extract_epi32(n0,3);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + 4*offset);
	Fd       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + 4*offset);
	Gd       = _mm_load_ps(VFtab + 4*(offset+2)* n_c + 4*offset);
	Hd       = _mm_load_ps(VFtab + 4*(offset+2)* n_d + 4*offset);
	Yr       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + 4*offset + 4);
	Fr       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + 4*offset + 4);
	Gr       = _mm_load_ps(VFtab + 4*(offset+2)* n_c + 4*offset + 4);
	Hr       = _mm_load_ps(VFtab + 4*(offset+2)* n_d + 4*offset + 4);
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fd        = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr        = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_4_table_coulomb_and_lj_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 c6, __m128 c12, 
								  __m128 *vctot, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	n_d      = gmx_mm_extract_epi32(n0,3);
	
	
	Yc       = _mm_load_ps(VFtab + 12* n_a);
	Fc       = _mm_load_ps(VFtab + 12* n_b);
	Gc       = _mm_load_ps(VFtab + 12* n_c);
	Hc       = _mm_load_ps(VFtab + 12* n_d);
	Yd       = _mm_load_ps(VFtab + 12* n_a + 4);
	Fd       = _mm_load_ps(VFtab + 12* n_b + 4);
	Gd       = _mm_load_ps(VFtab + 12* n_c + 4);
	Hd       = _mm_load_ps(VFtab + 12* n_d + 4);
	Yr       = _mm_load_ps(VFtab + 12* n_a + 8);
	Fr       = _mm_load_ps(VFtab + 12* n_b + 8);
	Gr       = _mm_load_ps(VFtab + 12* n_c + 8);
	Hr       = _mm_load_ps(VFtab + 12* n_d + 8);
	_MM_TRANSPOSE4_PS(Yc,Fc,Gc,Hc);
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hc       = _mm_mul_ps(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_ps(Gc,eps);               /* Geps  */
	Fc       = _mm_add_ps(Fc, _mm_add_ps(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Yc, _mm_mul_ps(eps,Fc)));
	*vctot   = _mm_add_ps(*vctot,vcoul);

	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fc       = _mm_mul_ps(qq, _mm_add_ps(Fc, _mm_add_ps(Gc, _mm_add_ps(Hc,Hc))));
	Fd       = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr       = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fc,_mm_add_ps(Fd,Fr)),tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_3_table_coulomb_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 *vctot)
{
    __m128  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i n0;
	int     n_a,n_b,n_c;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	Y        = _mm_load_ps(VFtab + 4* n_a);
	F        = _mm_load_ps(VFtab + 4* n_b);
	G        = _mm_load_ps(VFtab + 4* n_c);
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	F        = _mm_add_ps(F, _mm_add_ps(G,H));  /* Fp    */
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Y, _mm_mul_ps(eps,F)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	F        = _mm_mul_ps(qq, _mm_add_ps(F, _mm_add_ps(G, _mm_add_ps(H,H))));
	
	return _mm_mul_ps(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_3_table_lj_ps(__m128 r, __m128 tabscale, float * VFtab, int offset, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b,n_c;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset);
	Fd       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + offset);
	Gd       = _mm_load_ps(VFtab + 4*(offset+2)* n_c + offset);
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset + 4);
	Fr       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + offset + 4);
	Gr       = _mm_load_ps(VFtab + 4*(offset+2)* n_c + offset + 4);
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fd        = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr        = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_3_table_coulomb_and_lj_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 c6, __m128 c12, 
								  __m128 *vctot, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b,n_c;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	
	
	Yc       = _mm_load_ps(VFtab + 12* n_a);
	Fc       = _mm_load_ps(VFtab + 12* n_b);
	Gc       = _mm_load_ps(VFtab + 12* n_c);
	Hc       = _mm_setzero_ps();
	Yd       = _mm_load_ps(VFtab + 12* n_a + 4);
	Fd       = _mm_load_ps(VFtab + 12* n_b + 4);
	Gd       = _mm_load_ps(VFtab + 12* n_c + 4);
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 12* n_a + 8);
	Fr       = _mm_load_ps(VFtab + 12* n_b + 8);
	Gr       = _mm_load_ps(VFtab + 12* n_c + 8);
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yc,Fc,Gc,Hc);
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hc       = _mm_mul_ps(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_ps(Gc,eps);               /* Geps  */
	Fc       = _mm_add_ps(Fc, _mm_add_ps(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Yc, _mm_mul_ps(eps,Fc)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fc       = _mm_mul_ps(qq, _mm_add_ps(Fc, _mm_add_ps(Gc, _mm_add_ps(Hc,Hc))));
	Fd       = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr       = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fc,_mm_add_ps(Fd,Fr)),tabscale);
}





/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_2_table_coulomb_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 *vctot)
{
    __m128  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_ps(VFtab + 4* n_a);
	F        = _mm_load_ps(VFtab + 4* n_b);
	G        = _mm_setzero_ps();
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	F        = _mm_add_ps(F, _mm_add_ps(G,H));  /* Fp    */
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Y, _mm_mul_ps(eps,F)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	F        = _mm_mul_ps(qq, _mm_add_ps(F, _mm_add_ps(G, _mm_add_ps(H,H))));
	
	return _mm_mul_ps(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_2_table_lj_ps(__m128 r, __m128 tabscale, float * VFtab, int offset, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset);
	Fd       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + offset);
	Gd       = _mm_setzero_ps();
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset + 4);
	Fr       = _mm_load_ps(VFtab + 4*(offset+2)* n_b + offset + 4);
	Gr       = _mm_setzero_ps();
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fd        = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr        = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_2_table_coulomb_and_lj_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 c6, __m128 c12, 
								  __m128 *vctot, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a,n_b;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	
	Yc       = _mm_load_ps(VFtab + 12* n_a);
	Fc       = _mm_load_ps(VFtab + 12* n_b);
	Gc       = _mm_setzero_ps();
	Hc       = _mm_setzero_ps();
	Yd       = _mm_load_ps(VFtab + 12* n_a + 4);
	Fd       = _mm_load_ps(VFtab + 12* n_b + 4);
	Gd       = _mm_setzero_ps();
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 12* n_a + 8);
	Fr       = _mm_load_ps(VFtab + 12* n_b + 8);
	Gr       = _mm_setzero_ps();
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yc,Fc,Gc,Hc);
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hc       = _mm_mul_ps(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_ps(Gc,eps);               /* Geps  */
	Fc       = _mm_add_ps(Fc, _mm_add_ps(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Yc, _mm_mul_ps(eps,Fc)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fc       = _mm_mul_ps(qq, _mm_add_ps(Fc, _mm_add_ps(Gc, _mm_add_ps(Hc,Hc))));
	Fd       = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr       = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fc,_mm_add_ps(Fd,Fr)),tabscale);
}




/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_1_table_coulomb_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 *vctot)
{
    __m128  rt,eps,eps2,Y,F,G,H,vcoul;
	__m128i n0;
	int     n_a;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	Y        = _mm_load_ps(VFtab + 4* n_a);
	F        = _mm_setzero_ps();
	G        = _mm_setzero_ps();
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	F        = _mm_add_ps(F, _mm_add_ps(G,H));  /* Fp    */
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Y, _mm_mul_ps(eps,F)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	F        = _mm_mul_ps(qq, _mm_add_ps(F, _mm_add_ps(G, _mm_add_ps(H,H))));
	
	return _mm_mul_ps(F,tabscale);
}



/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_1_table_lj_ps(__m128 r, __m128 tabscale, float * VFtab, int offset, __m128 c6, __m128 c12, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	/* For a few cases, like TIP4p waters, there are particles with LJ-only interactions in a loop where
	 * the table data might contain both coulomb and LJ. To handle this case, we use an offset value of 0
	 * if the data is an LJ-only table, and 1 if it is actually a mixed coul+lj table.
	 */
	Yd       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset);
	Fd       = _mm_setzero_ps();
	Gd       = _mm_setzero_ps();
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 4*(offset+2)* n_a + offset + 4);
	Fr       = _mm_setzero_ps();
	Gr       = _mm_setzero_ps();
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fd        = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr        = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fd,Fr),tabscale);
}


/* Return force should be multiplied by -rinv to get fscal */
static inline __m128
gmx_mm_int_1_table_coulomb_and_lj_ps(__m128 r, __m128 tabscale, float * VFtab, __m128 qq, __m128 c6, __m128 c12, 
									 __m128 *vctot, __m128 *vvdwtot)
{
    __m128  rt,eps,eps2,vcoul,Yc,Fc,Gc,Hc,Yd,Fd,Gd,Hd,Yr,Fr,Gr,Hr,vvdw6,vvdw12;
	__m128i n0;
	int     n_a;
	
    rt       = _mm_mul_ps(r,tabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	
	Yc       = _mm_load_ps(VFtab + 12* n_a);
	Fc       = _mm_setzero_ps();
	Gc       = _mm_setzero_ps();
	Hc       = _mm_setzero_ps();
	Yd       = _mm_load_ps(VFtab + 12* n_a + 4);
	Fd       = _mm_setzero_ps();
	Gd       = _mm_setzero_ps();
	Hd       = _mm_setzero_ps();
	Yr       = _mm_load_ps(VFtab + 12* n_a + 8);
	Fr       = _mm_setzero_ps();
	Gr       = _mm_setzero_ps();
	Hr       = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Yc,Fc,Gc,Hc);
	_MM_TRANSPOSE4_PS(Yd,Fd,Gd,Hd);
	_MM_TRANSPOSE4_PS(Yr,Fr,Gr,Hr);
	Hc       = _mm_mul_ps(Hc,eps2);              /* Heps2 */
	Gc       = _mm_mul_ps(Gc,eps);               /* Geps  */
	Fc       = _mm_add_ps(Fc, _mm_add_ps(Gc,Hc));  /* Fp    */
	Hd       = _mm_mul_ps(Hd,eps2);              /* Heps2 */
	Gd       = _mm_mul_ps(Gd,eps);               /* Geps  */
	Fd       = _mm_add_ps(Fd, _mm_add_ps(Gd,Hd));  /* Fp    */
	Hr       = _mm_mul_ps(Hr,eps2);              /* Heps2 */
	Gr       = _mm_mul_ps(Gr,eps);               /* Geps  */
	Fr       = _mm_add_ps(Fr, _mm_add_ps(Gr,Hr));  /* Fp    */
	
	vcoul    = _mm_mul_ps(qq, _mm_add_ps(Yc, _mm_mul_ps(eps,Fc)));
	*vctot   = _mm_add_ps(*vctot,vcoul);
	
	vvdw6    = _mm_mul_ps(c6,  _mm_add_ps(Yd, _mm_mul_ps(eps,Fd)));
	vvdw12   = _mm_mul_ps(c12, _mm_add_ps(Yr, _mm_mul_ps(eps,Fr)));
	*vvdwtot = _mm_add_ps(*vvdwtot, _mm_add_ps(vvdw6,vvdw12));
	
	Fc       = _mm_mul_ps(qq, _mm_add_ps(Fc, _mm_add_ps(Gc, _mm_add_ps(Hc,Hc))));
	Fd       = _mm_mul_ps(c6,  _mm_add_ps(Fd, _mm_add_ps(Gd, _mm_add_ps(Hd,Hd))));
	Fr       = _mm_mul_ps(c12, _mm_add_ps(Fr, _mm_add_ps(Gr, _mm_add_ps(Hr,Hr))));
	
	return _mm_mul_ps( _mm_add_ps(Fc,_mm_add_ps(Fd,Fr)),tabscale);
}





/* Return force should be multiplied by +rinv to get fscal */
static inline __m128
gmx_mm_int_4_genborn_ps(__m128 r, __m128 isai, 
						float * isaj1, float *isaj2, float *isaj3, float *isaj4,
						__m128 gbtabscale, float * GBtab, __m128 qq, __m128 *dvdasum, 
						float *dvdaj1, float *dvdaj2, float *dvdaj3, float *dvdaj4, 
						__m128 *vgbtot)
{
	const __m128 half  = {0.5,0.5,0.5,0.5};
	
    __m128  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_ss(isaj1);
	t2       = _mm_load_ss(isaj2);
	t3       = _mm_load_ss(isaj3);
	t4       = _mm_load_ss(isaj4);
	isaj     = _mm_unpacklo_ps(isaj,t2);  /* - - t2 t1 */
	t3       = _mm_unpacklo_ps(t3,t4);  /* - - t4 t3 */
	isaj     = _mm_movelh_ps(isaj,t3); /* t4 t3 t2 t1 */
	
	isaprod     = _mm_mul_ps(isai,isaj);
	qq          = _mm_mul_ps(qq,isaprod);
	gbtabscale  = _mm_mul_ps( isaprod, gbtabscale );
	
	rt       = _mm_mul_ps(r,gbtabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
		
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	n_d      = gmx_mm_extract_epi32(n0,3);
	Y        = _mm_load_ps(GBtab + 4* n_a);
	F        = _mm_load_ps(GBtab + 4* n_b);
	G        = _mm_load_ps(GBtab + 4* n_c);
	H        = _mm_load_ps(GBtab + 4* n_d);
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	F        = _mm_add_ps(_mm_add_ps(F,G),H);  /* Fp    */
		
	VV       = _mm_add_ps(Y, _mm_mul_ps(eps,F));
	FF       = _mm_add_ps(_mm_add_ps(F,G), _mm_add_ps(H,H));

	vgb      = _mm_mul_ps(qq, VV);
	*vgbtot  = _mm_sub_ps(*vgbtot,vgb); /* Yes, the sign is correct */

	ftmp     = _mm_mul_ps(_mm_mul_ps(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_ps(half, _mm_add_ps(vgb,_mm_mul_ps(ftmp,r)));
	
	*dvdasum = _mm_add_ps(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_ps(_mm_mul_ps(dvdatmp,isaj), isaj);

	/* Update 4 dada[j] values */
	Y        = _mm_load_ss(dvdaj1);
	F        = _mm_load_ss(dvdaj2);
	G        = _mm_load_ss(dvdaj3);
	H        = _mm_load_ss(dvdaj4);
	t3       = _mm_movehl_ps(_mm_setzero_ps(),dvdatmp);
	t2       = _mm_shuffle_ps(dvdatmp,dvdatmp,_MM_SHUFFLE(0,0,0,1));
	t4       = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(0,0,0,1));
		
	_mm_store_ss( dvdaj1 , _mm_add_ss( Y, dvdatmp ) );
	_mm_store_ss( dvdaj2 , _mm_add_ss( F, t2 ) );
	_mm_store_ss( dvdaj3 , _mm_add_ss( G, t3 ) );
	_mm_store_ss( dvdaj4 , _mm_add_ss( H, t4 ) );
	
	return ftmp;
}



/* Return force should be multiplied by +rinv to get fscal */
static inline __m128
gmx_mm_int_3_genborn_ps(__m128 r, __m128 isai, 
						float * isaj1, float *isaj2, float *isaj3, 
						__m128 gbtabscale, float * GBtab, __m128 qq, __m128 *dvdasum, 
						float *dvdaj1, float *dvdaj2, float *dvdaj3,
						__m128 *vgbtot)
{
	const __m128 half  = {0.5,0.5,0.5,0.5};
	
    __m128  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_ss(isaj1);
	t2       = _mm_load_ss(isaj2);
	t3       = _mm_load_ss(isaj3);
	isaj     = _mm_unpacklo_ps(isaj,t2);  /* - - t2 t1 */
	t3       = _mm_unpacklo_ps(t3,t3);  /* - - t3 t3 */
	isaj     = _mm_movelh_ps(isaj,t3); /* t3 t3 t2 t1 */
	
	isaprod     = _mm_mul_ps(isai,isaj);
	qq          = _mm_mul_ps(qq,isaprod);
	gbtabscale  = _mm_mul_ps( isaprod, gbtabscale );
	
	rt       = _mm_mul_ps(r,gbtabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	n_c      = gmx_mm_extract_epi32(n0,2);
	Y        = _mm_load_ps(GBtab + 4* n_a);
	F        = _mm_load_ps(GBtab + 4* n_b);
	G        = _mm_load_ps(GBtab + 4* n_c);
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	F        = _mm_add_ps(_mm_add_ps(F,G),H);  /* Fp    */
	
	VV       = _mm_add_ps(Y, _mm_mul_ps(eps,F));
	FF       = _mm_add_ps(_mm_add_ps(F,G), _mm_add_ps(H,H));
	
	vgb      = _mm_mul_ps(qq, VV);
	*vgbtot  = _mm_sub_ps(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_ps(_mm_mul_ps(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_ps(half, _mm_add_ps(vgb,_mm_mul_ps(ftmp,r)));
	
	*dvdasum = _mm_add_ps(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_ps(_mm_mul_ps(dvdatmp,isaj), isaj);
	
	/* Update 3 dada[j] values */
	Y        = _mm_load_ss(dvdaj1);
	F        = _mm_load_ss(dvdaj2);
	G        = _mm_load_ss(dvdaj3);
	t3       = _mm_movehl_ps(_mm_setzero_ps(),dvdatmp);
	t2       = _mm_shuffle_ps(dvdatmp,dvdatmp,_MM_SHUFFLE(0,0,0,1));
	
	_mm_store_ss( dvdaj1 , _mm_add_ss( Y, dvdatmp ) );
	_mm_store_ss( dvdaj2 , _mm_add_ss( F, t2 ) );
	_mm_store_ss( dvdaj3 , _mm_add_ss( G, t3 ) );
	
	return ftmp;
}




/* Return force should be multiplied by +rinv to get fscal */
static inline __m128
gmx_mm_int_2_genborn_ps(__m128 r, __m128 isai, 
						float * isaj1, float *isaj2, 
						__m128 gbtabscale, float * GBtab, __m128 qq, __m128 *dvdasum, 
						float *dvdaj1, float *dvdaj2,
						__m128 *vgbtot)
{
	const __m128 half  = {0.5,0.5,0.5,0.5};
	
    __m128  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_ss(isaj1);
	t2       = _mm_load_ss(isaj2);
	isaj     = _mm_unpacklo_ps(isaj,t2);  /* - - t2 t1 */
	
	isaprod     = _mm_mul_ps(isai,isaj);
	qq          = _mm_mul_ps(qq,isaprod);
	gbtabscale  = _mm_mul_ps( isaprod, gbtabscale );
	
	rt       = _mm_mul_ps(r,gbtabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	n_b      = gmx_mm_extract_epi32(n0,1);
	Y        = _mm_load_ps(GBtab + 4* n_a);
	F        = _mm_load_ps(GBtab + 4* n_b);
	G        = _mm_setzero_ps();
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	F        = _mm_add_ps(_mm_add_ps(F,G),H);  /* Fp    */
	
	VV       = _mm_add_ps(Y, _mm_mul_ps(eps,F));
	FF       = _mm_add_ps(_mm_add_ps(F,G), _mm_add_ps(H,H));
	
	vgb      = _mm_mul_ps(qq, VV);
	*vgbtot  = _mm_sub_ps(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_ps(_mm_mul_ps(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_ps(half, _mm_add_ps(vgb,_mm_mul_ps(ftmp,r)));
	
	*dvdasum = _mm_add_ps(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_ps(_mm_mul_ps(dvdatmp,isaj), isaj);
	
	/* Update 2 dada[j] values */
	Y        = _mm_load_ss(dvdaj1);
	F        = _mm_load_ss(dvdaj2);
	t2       = _mm_shuffle_ps(dvdatmp,dvdatmp,_MM_SHUFFLE(0,0,0,1));
	
	_mm_store_ss( dvdaj1 , _mm_add_ss( Y, dvdatmp ) );
	_mm_store_ss( dvdaj2 , _mm_add_ss( F, t2 ) );
	
	return ftmp;
}

/* Return force should be multiplied by +rinv to get fscal */
static inline __m128
gmx_mm_int_1_genborn_ps(__m128 r, __m128 isai, 
						float * isaj1, 
						__m128 gbtabscale, float * GBtab, __m128 qq, __m128 *dvdasum, 
						float *dvdaj1, 
						__m128 *vgbtot)
{
	const __m128 half  = {0.5,0.5,0.5,0.5};
	
    __m128  rt,eps,eps2,Y,F,G,H,VV,FF,ftmp,isaprod,t2,t3,t4,isaj,vgb,dvdatmp;
	__m128i n0;
	int     n_a,n_b,n_c,n_d;
	
	/* Assemble isaj */
	isaj     = _mm_load_ss(isaj1);
	
	isaprod     = _mm_mul_ps(isai,isaj);
	qq          = _mm_mul_ps(qq,isaprod);
	gbtabscale  = _mm_mul_ps( isaprod, gbtabscale );
	
	rt       = _mm_mul_ps(r,gbtabscale); 
	n0       = _mm_cvttps_epi32(rt);
	eps      = _mm_sub_ps(rt, _mm_cvtepi32_ps(n0));
	eps2     = _mm_mul_ps(eps,eps);
	
	/* Extract indices from n0 */
	n_a      = gmx_mm_extract_epi32(n0,0);
	Y        = _mm_load_ps(GBtab + 4* n_a);
	F        = _mm_setzero_ps();
	G        = _mm_setzero_ps();
	H        = _mm_setzero_ps();
	_MM_TRANSPOSE4_PS(Y,F,G,H);
	G        = _mm_mul_ps(G,eps);               /* Geps  */
	H        = _mm_mul_ps(H,eps2);              /* Heps2 */
	F        = _mm_add_ps(_mm_add_ps(F,G),H);  /* Fp    */
	
	VV       = _mm_add_ps(Y, _mm_mul_ps(eps,F));
	FF       = _mm_add_ps(_mm_add_ps(F,G), _mm_add_ps(H,H));
	
	vgb      = _mm_mul_ps(qq, VV);
	*vgbtot  = _mm_sub_ps(*vgbtot,vgb); /* Yes, the sign is correct */
	
	ftmp     = _mm_mul_ps(_mm_mul_ps(qq, FF), gbtabscale);
	
	dvdatmp  = _mm_mul_ps(half, _mm_add_ps(vgb,_mm_mul_ps(ftmp,r)));
	
	*dvdasum = _mm_add_ps(*dvdasum,dvdatmp);
	
	dvdatmp  = _mm_mul_ps(_mm_mul_ps(dvdatmp,isaj), isaj);
	
	/* Update 1 dada[j] values */
	Y        = _mm_load_ss(dvdaj1);
	
	_mm_store_ss( dvdaj1 , _mm_add_ss( Y, dvdatmp ) );
	
	return ftmp;
}





static inline void
gmx_mm_update_iforce_1atom_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
							  float *fptr,
							  float *fshiftptr)
{
	__m128 t1,t2,t3;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_ps(fix1,fix1);   
	fiy1 = _mm_hadd_ps(fiy1,fiz1);
	
	fix1 = _mm_hadd_ps(fix1,fiy1); /* fiz1 fiy1 fix1 fix1 */ 
#else
	/* SSE2 */
	/* transpose data */
	t1 = fix1;
	_MM_TRANSPOSE4_PS(fix1,t1,fiy1,fiz1);  
	fix1 = _mm_add_ps(_mm_add_ps(fix1,t1), _mm_add_ps(fiy1,fiz1));
#endif
	t2 = _mm_load_ss(fptr);
	t2 = _mm_loadh_pi(t2,(__m64 *)(fptr+1));
	t3 = _mm_load_ss(fshiftptr);
	t3 = _mm_loadh_pi(t3,(__m64 *)(fshiftptr+1));
	
	t2 = _mm_add_ps(t2,fix1);
	t3 = _mm_add_ps(t3,fix1);
	
	_mm_store_ss(fptr,t2);
	_mm_storeh_pi((__m64 *)(fptr+1),t2);
	_mm_store_ss(fshiftptr,t3);
	_mm_storeh_pi((__m64 *)(fshiftptr+1),t3);
}

static inline void
gmx_mm_update_iforce_2atoms_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
							   __m128 fix2, __m128 fiy2, __m128 fiz2,
							   float *fptr,
							   float *fshiftptr)
{
	__m128 t1,t2,t3,t4;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_ps(fix1,fiy1);   
	fiz1 = _mm_hadd_ps(fiz1,fix2);
	fiy2 = _mm_hadd_ps(fiy2,fiz2);
	
	fix1 = _mm_hadd_ps(fix1,fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	fiy2 = _mm_hadd_ps(fiy2,fiy2); /*  -    -   fiz2 fiy2 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);  
	t1 = _mm_unpacklo_ps(fiy2,fiz2);
	t2 = _mm_unpackhi_ps(fiy2,fiz2);
		
	fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));
	t1   = _mm_add_ps(t1,t2);
	t2   = _mm_movehl_ps(t2,t1);
	fiy2 = _mm_add_ps(t1,t2);
#endif
	_mm_storeu_ps(fptr,   _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
	t1 = _mm_loadl_pi(t1,(__m64 *)(fptr+4));
	_mm_storel_pi((__m64 *)(fptr+4), _mm_add_ps(fiy2,t1));
	
	t4 = _mm_load_ss(fshiftptr+2);
	t4 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(0,0,3,2));   /* fiy2  -   fix2 fiz1 */
	t1 = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,1,0,0));       /* fiy2 fix2  -   fiz1 */
	t2 = _mm_shuffle_ps(fiy2,fix1,_MM_SHUFFLE(1,0,0,1));   /* fiy1 fix1  -   fiz2 */

	t1 = _mm_add_ps(t1,t2);
	t3 = _mm_add_ps(t3,t4);
	t1 = _mm_add_ps(t1,t3); /* y x - z */
	
	_mm_store_ss(fshiftptr+2,t1);
	_mm_storeh_pi((__m64 *)(fshiftptr),t1);
}



static inline void
gmx_mm_update_iforce_3atoms_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
							   __m128 fix2, __m128 fiy2, __m128 fiz2,
							   __m128 fix3, __m128 fiy3, __m128 fiz3,
							   float *fptr,
							   float *fshiftptr)
{
	__m128 t1,t2,t3,t4;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_ps(fix1,fiy1);   
	fiz1 = _mm_hadd_ps(fiz1,fix2);
	fiy2 = _mm_hadd_ps(fiy2,fiz2);
	fix3 = _mm_hadd_ps(fix3,fiy3);
	fiz3 = _mm_hadd_ps(fiz3,fiz3);
	
	fix1 = _mm_hadd_ps(fix1,fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	fiy2 = _mm_hadd_ps(fiy2,fix3); /* fiy3 fix3 fiz2 fiy2 */
	fiz3 = _mm_hadd_ps(fiz3,fiz3); /*  -    -    -   fiz3 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);  
	_MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);
	t2   = _mm_movehl_ps(t2,fiz3);
	t1   = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(0,0,0,1));
	t3   = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(0,0,0,1));
	
	fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));
	fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));
	fiz3 = _mm_add_ss(_mm_add_ps(fiz3,t1)  , _mm_add_ps(t2,t3));
#endif
	_mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
	_mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));
	_mm_store_ss (fptr+8,_mm_add_ss(fiz3,_mm_load_ss(fptr+8) ));
	
	t4 = _mm_load_ss(fshiftptr+2);
	t4 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(fiz3,fix1,_MM_SHUFFLE(1,0,0,0));   /* fiy1 fix1  -   fiz3 */
	t2 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(3,2,2,2));   /* fiy3 fix3  -   fiz1 */
	t3 = _mm_shuffle_ps(fiy2,fix1,_MM_SHUFFLE(3,3,0,1));   /* fix2 fix2 fiy2 fiz2 */
	t3 = _mm_shuffle_ps(t3  ,t3  ,_MM_SHUFFLE(1,2,0,0));   /* fiy2 fix2  -   fiz2 */

	t1 = _mm_add_ps(t1,t2);
	t3 = _mm_add_ps(t3,t4);
	t1 = _mm_add_ps(t1,t3); /* y x - z */
	
	_mm_store_ss(fshiftptr+2,t1);
	_mm_storeh_pi((__m64 *)(fshiftptr),t1);
}

static inline void
gmx_mm_update_iforce_4atoms_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
							   __m128 fix2, __m128 fiy2, __m128 fiz2,
							   __m128 fix3, __m128 fiy3, __m128 fiz3,
							   __m128 fix4, __m128 fiy4, __m128 fiz4,
							   float *fptr,
							   float *fshiftptr)
{
	__m128 t1,t2,t3,t4,t5;
	
#ifdef GMX_SSE3
	fix1 = _mm_hadd_ps(fix1,fiy1);   
	fiz1 = _mm_hadd_ps(fiz1,fix2);
	fiy2 = _mm_hadd_ps(fiy2,fiz2);
	fix3 = _mm_hadd_ps(fix3,fiy3);
	fiz3 = _mm_hadd_ps(fiz3,fix4);
	fiy4 = _mm_hadd_ps(fiy4,fiz4);
	
	fix1 = _mm_hadd_ps(fix1,fiz1); /* fix2 fiz1 fiy1 fix1 */ 
	fiy2 = _mm_hadd_ps(fiy2,fix3); /* fiy3 fix3 fiz2 fiy2 */
	fiz3 = _mm_hadd_ps(fiz3,fiy4); /* fiz4 fiy4 fix4 fiz3 */
#else
	/* SSE2 */
	/* transpose data */
	_MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);  
	_MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);
	_MM_TRANSPOSE4_PS(fiz3,fix4,fiy4,fiz4);
	
	fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));
	fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));
	fiz3 = _mm_add_ps(_mm_add_ps(fiz3,fix4), _mm_add_ps(fiy4,fiz4));
#endif
	_mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
	_mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));
	_mm_storeu_ps(fptr+8,_mm_add_ps(fiz3,_mm_loadu_ps(fptr+8)));
	
	t5 = _mm_load_ss(fshiftptr+2);
	t5 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));
	
	t1 = _mm_shuffle_ps(fix1,fix1,_MM_SHUFFLE(1,0,2,2));   /* fiy1 fix1  -   fiz1 */
	t2 = _mm_shuffle_ps(fiy2,fiy2,_MM_SHUFFLE(3,2,1,1));   /* fiy3 fix3  -   fiz2 */
	t3 = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(2,1,0,0));   /* fiy4 fix4  -   fiz3 */
	t4 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(0,0,3,3));   /* fiy2 fiy2 fix2 fix2 */
	t4 = _mm_shuffle_ps(fiz3,t4  ,_MM_SHUFFLE(2,0,3,3));   /* fiy2 fix2  -   fiz4 */
	
	t1 = _mm_add_ps(t1,t2);
	t3 = _mm_add_ps(t3,t4);
	t1 = _mm_add_ps(t1,t3); /* y x - z */
	t5 = _mm_add_ps(t5,t1);
	
	_mm_store_ss(fshiftptr+2,t5);
	_mm_storeh_pi((__m64 *)(fshiftptr),t5);
}


static inline void
gmx_mm_update_1pot_ps(__m128 pot1, float *ptr1)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_ps(pot1,pot1);
	pot1 = _mm_hadd_ps(pot1,pot1);
#else
	/* SSE2 */
	pot1 = _mm_add_ps(pot1,_mm_movehl_ps(pot1,pot1));
	pot1 = _mm_add_ps(pot1,_mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(0,0,0,1)));
#endif
	_mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));
}
				   
static inline void
gmx_mm_update_2pot_ps(__m128 pot1, float *ptr1, __m128 pot2, float *ptr2)
{
#ifdef GMX_SSE3
	pot1 = _mm_hadd_ps(pot1,pot2); 
	pot1 = _mm_hadd_ps(pot1,pot1); 
	pot2 = _mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(0,0,0,1));
#else
	/* SSE2 */
	__m128 t1,t2;
	t1   = _mm_movehl_ps(pot2,pot1); /* 2d 2c 1d 1c */
	t2   = _mm_movelh_ps(pot1,pot2); /* 2b 2a 1b 1a */
	t1   = _mm_add_ps(t1,t2);       /* 2  2  1  1  */
	t2   = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,3,1,1));
	pot1 = _mm_add_ps(t1,t2);       /* -  2  -  1  */
	pot2 = _mm_movehl_ps(t2,t1);    /* -  -  -  2  */
#endif
	_mm_store_ss(ptr1,_mm_add_ss(pot1,_mm_load_ss(ptr1)));
	_mm_store_ss(ptr2,_mm_add_ss(pot2,_mm_load_ss(ptr2)));
}


