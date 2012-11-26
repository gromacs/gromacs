/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 2011-2012, The GROMACS Development Team
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
#ifndef _kernelutil_x86_avx_128_fma_single_h_
#define _kernelutil_x86_avx_128_fma_single_h_


#include <math.h>

#include "gmx_x86_avx_128_fma.h"

/* Normal sum of four xmm registers */
#define gmx_mm_sum4_ps(t0,t1,t2,t3)  _mm_add_ps(_mm_add_ps(t0,t1),_mm_add_ps(t2,t3))

static gmx_inline int
gmx_mm_any_lt(__m128 a, __m128 b)
{
    return _mm_movemask_ps(_mm_cmplt_ps(a,b));
}

static gmx_inline __m128
gmx_mm_calc_rsq_ps(__m128 dx, __m128 dy, __m128 dz)
{
    return _mm_macc_ps(dx,dx,_mm_macc_ps(dy,dy,_mm_mul_ps(dz,dz)));
}

/* Load a single value from 1-4 places, merge into xmm register */

static gmx_inline __m128
gmx_mm_load_4real_swizzle_ps(const float * gmx_restrict ptrA,
                             const float * gmx_restrict ptrB,
                             const float * gmx_restrict ptrC,
                             const float * gmx_restrict ptrD)
{
    __m128 t1,t2;
    
    t1 = _mm_unpacklo_ps(_mm_load_ss(ptrA),_mm_load_ss(ptrC));
    t2 = _mm_unpacklo_ps(_mm_load_ss(ptrB),_mm_load_ss(ptrD));
    return _mm_unpacklo_ps(t1,t2);
}


static gmx_inline void
gmx_mm_store_4real_swizzle_ps(float * gmx_restrict ptrA,
                              float * gmx_restrict ptrB,
                              float * gmx_restrict ptrC,
                              float * gmx_restrict ptrD, __m128 xmm1)
{
    __m128 t2,t3,t4;

    t2       = _mm_permute_ps(xmm1,_MM_SHUFFLE(1,1,1,1));
    t3       = _mm_permute_ps(xmm1,_MM_SHUFFLE(2,2,2,2));
    t4       = _mm_permute_ps(xmm1,_MM_SHUFFLE(3,3,3,3));
    _mm_store_ss(ptrA,xmm1);
    _mm_store_ss(ptrB,t2);
    _mm_store_ss(ptrC,t3);
    _mm_store_ss(ptrD,t4);
}


static gmx_inline void
gmx_mm_increment_4real_swizzle_ps(float * gmx_restrict ptrA,
                                  float * gmx_restrict ptrB,
                                  float * gmx_restrict ptrC,
                                  float * gmx_restrict ptrD, __m128 xmm1)
{
    __m128 tmp;

    tmp = gmx_mm_load_4real_swizzle_ps(ptrA,ptrB,ptrC,ptrD);
    tmp = _mm_add_ps(tmp,xmm1);
    gmx_mm_store_4real_swizzle_ps(ptrA,ptrB,ptrC,ptrD,tmp);
}


static gmx_inline void
gmx_mm_load_4pair_swizzle_ps(const float * gmx_restrict p1,
                             const float * gmx_restrict p2,
                             const float * gmx_restrict p3,
                             const float * gmx_restrict p4,
                             __m128 * gmx_restrict c6, __m128 * gmx_restrict c12)
{
    __m128 t1,t2,t3,t4;
    t1   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p1);
    t2   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p2);
    t3   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p3);
    t4   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p4);
    t1   = _mm_unpacklo_ps(t1,t3);
    t2   = _mm_unpacklo_ps(t2,t4);
    *c6  = _mm_unpacklo_ps(t1,t2);
    *c12 = _mm_unpackhi_ps(t1,t2);
}




static gmx_inline void
gmx_mm_load_shift_and_1rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                         const float * gmx_restrict xyz,
                                         __m128 * gmx_restrict x1,
                                         __m128 * gmx_restrict y1,
                                         __m128 * gmx_restrict z1)
{
    __m128 t1,t2,t3,t4;
    
    t1   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)xyz_shift);
    t2   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)xyz);
    t3   = _mm_load_ss(xyz_shift+2);
    t4   = _mm_load_ss(xyz+2);
    t1   = _mm_add_ps(t1,t2);
    t3   = _mm_add_ss(t3,t4);
    
    *x1  = _mm_permute_ps(t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_permute_ps(t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_permute_ps(t3,_MM_SHUFFLE(0,0,0,0));
}


static gmx_inline void
gmx_mm_load_shift_and_3rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                         const float * gmx_restrict xyz,
                                         __m128 * gmx_restrict x1, __m128 * gmx_restrict y1, __m128 * gmx_restrict z1,
                                         __m128 * gmx_restrict x2, __m128 * gmx_restrict y2, __m128 * gmx_restrict z2,
                                         __m128 * gmx_restrict x3, __m128 * gmx_restrict y3, __m128 * gmx_restrict z3)
{
    __m128 tA,tB;
    __m128 t1,t2,t3,t4,t5,t6;
    
    tA   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)xyz_shift);
    tB   = _mm_load_ss(xyz_shift+2);
    
    t1   = _mm_loadu_ps(xyz);
    t2   = _mm_loadu_ps(xyz+4);
    t3   = _mm_load_ss(xyz+8);
    
    tA   = _mm_movelh_ps(tA,tB);
    t4   = _mm_permute_ps(tA,_MM_SHUFFLE(0,2,1,0));
    t5   = _mm_permute_ps(tA,_MM_SHUFFLE(1,0,2,1));
    t6   = _mm_permute_ps(tA,_MM_SHUFFLE(2,1,0,2));
    
    t1   = _mm_add_ps(t1,t4);
    t2   = _mm_add_ps(t2,t5);
    t3   = _mm_add_ss(t3,t6);
    
    *x1  = _mm_permute_ps(t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_permute_ps(t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_permute_ps(t1,_MM_SHUFFLE(2,2,2,2));
    *x2  = _mm_permute_ps(t1,_MM_SHUFFLE(3,3,3,3));
    *y2  = _mm_permute_ps(t2,_MM_SHUFFLE(0,0,0,0));
    *z2  = _mm_permute_ps(t2,_MM_SHUFFLE(1,1,1,1));
    *x3  = _mm_permute_ps(t2,_MM_SHUFFLE(2,2,2,2));
    *y3  = _mm_permute_ps(t2,_MM_SHUFFLE(3,3,3,3));
    *z3  = _mm_permute_ps(t3,_MM_SHUFFLE(0,0,0,0));
}


static gmx_inline void
gmx_mm_load_shift_and_4rvec_broadcast_ps(const float * gmx_restrict xyz_shift,
                                         const float * gmx_restrict xyz,
                                         __m128 * gmx_restrict x1, __m128 * gmx_restrict y1, __m128 * gmx_restrict z1,
                                         __m128 * gmx_restrict x2, __m128 * gmx_restrict y2, __m128 * gmx_restrict z2,
                                         __m128 * gmx_restrict x3, __m128 * gmx_restrict y3, __m128 * gmx_restrict z3,
                                         __m128 * gmx_restrict x4, __m128 * gmx_restrict y4, __m128 * gmx_restrict z4)
{
    __m128 tA,tB;
    __m128 t1,t2,t3,t4,t5,t6;
    
    tA   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)xyz_shift);
    tB   = _mm_load_ss(xyz_shift+2);
    
    t1   = _mm_loadu_ps(xyz);
    t2   = _mm_loadu_ps(xyz+4);
    t3   = _mm_loadu_ps(xyz+8);
    
    tA   = _mm_movelh_ps(tA,tB);
    t4   = _mm_permute_ps(tA,_MM_SHUFFLE(0,2,1,0));
    t5   = _mm_permute_ps(tA,_MM_SHUFFLE(1,0,2,1));
    t6   = _mm_permute_ps(tA,_MM_SHUFFLE(2,1,0,2));
    
    t1   = _mm_add_ps(t1,t4);
    t2   = _mm_add_ps(t2,t5);
    t3   = _mm_add_ps(t3,t6);
    
    *x1  = _mm_permute_ps(t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_permute_ps(t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_permute_ps(t1,_MM_SHUFFLE(2,2,2,2));
    *x2  = _mm_permute_ps(t1,_MM_SHUFFLE(3,3,3,3));
    *y2  = _mm_permute_ps(t2,_MM_SHUFFLE(0,0,0,0));
    *z2  = _mm_permute_ps(t2,_MM_SHUFFLE(1,1,1,1));
    *x3  = _mm_permute_ps(t2,_MM_SHUFFLE(2,2,2,2));
    *y3  = _mm_permute_ps(t2,_MM_SHUFFLE(3,3,3,3));
    *z3  = _mm_permute_ps(t3,_MM_SHUFFLE(0,0,0,0));
    *x4  = _mm_permute_ps(t3,_MM_SHUFFLE(1,1,1,1));
    *y4  = _mm_permute_ps(t3,_MM_SHUFFLE(2,2,2,2));
    *z4  = _mm_permute_ps(t3,_MM_SHUFFLE(3,3,3,3));
}


static gmx_inline void
gmx_mm_load_1rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                  __m128 * gmx_restrict x1, __m128 * gmx_restrict y1, __m128 * gmx_restrict z1)
{
    __m128 t1,t2,t3,t4;
    __m128i mask = _mm_set_epi32(0,-1,-1,-1);
    t1             = _mm_maskload_ps(ptrA,mask);
    t2             = _mm_maskload_ps(ptrB,mask);
    t3             = _mm_maskload_ps(ptrC,mask);
    t4             = _mm_maskload_ps(ptrD,mask);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *x1           = t1;
    *y1           = t2;
    *z1           = t3;
}


static gmx_inline void
gmx_mm_load_3rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                  __m128 * gmx_restrict x1, __m128 * gmx_restrict y1, __m128 * gmx_restrict z1,
                                  __m128 * gmx_restrict x2, __m128 * gmx_restrict y2, __m128 * gmx_restrict z2,
                                  __m128 * gmx_restrict x3, __m128 * gmx_restrict y3, __m128 * gmx_restrict z3)
{
    __m128 t1,t2,t3,t4;
    t1            = _mm_loadu_ps(ptrA);
    t2            = _mm_loadu_ps(ptrB);
    t3            = _mm_loadu_ps(ptrC);
    t4            = _mm_loadu_ps(ptrD);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *x1           = t1;
    *y1           = t2;
    *z1           = t3;
    *x2           = t4;
    t1            = _mm_loadu_ps(ptrA+4);
    t2            = _mm_loadu_ps(ptrB+4);
    t3            = _mm_loadu_ps(ptrC+4);
    t4            = _mm_loadu_ps(ptrD+4);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *y2           = t1;
    *z2           = t2;
    *x3           = t3;
    *y3           = t4;
    t1            = _mm_load_ss(ptrA+8);
    t2            = _mm_load_ss(ptrB+8);
    t3            = _mm_load_ss(ptrC+8);
    t4            = _mm_load_ss(ptrD+8);
    t1            = _mm_unpacklo_ps(t1,t3);
    t3            = _mm_unpacklo_ps(t2,t4);
    *z3           = _mm_unpacklo_ps(t1,t3);
}


static gmx_inline void
gmx_mm_load_4rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA, const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC, const float * gmx_restrict ptrD,
                                  __m128 * gmx_restrict x1, __m128 * gmx_restrict y1, __m128 * gmx_restrict z1,
                                  __m128 * gmx_restrict x2, __m128 * gmx_restrict y2, __m128 * gmx_restrict z2,
                                  __m128 * gmx_restrict x3, __m128 * gmx_restrict y3, __m128 * gmx_restrict z3,
                                  __m128 * gmx_restrict x4, __m128 * gmx_restrict y4, __m128 * gmx_restrict z4)
{
    __m128 t1,t2,t3,t4;
    t1            = _mm_loadu_ps(ptrA);
    t2            = _mm_loadu_ps(ptrB);
    t3            = _mm_loadu_ps(ptrC);
    t4            = _mm_loadu_ps(ptrD);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *x1           = t1;
    *y1           = t2;
    *z1           = t3;
    *x2           = t4;
    t1            = _mm_loadu_ps(ptrA+4);
    t2            = _mm_loadu_ps(ptrB+4);
    t3            = _mm_loadu_ps(ptrC+4);
    t4            = _mm_loadu_ps(ptrD+4);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *y2           = t1;
    *z2           = t2;
    *x3           = t3;
    *y3           = t4;
    t1            = _mm_loadu_ps(ptrA+8);
    t2            = _mm_loadu_ps(ptrB+8);
    t3            = _mm_loadu_ps(ptrC+8);
    t4            = _mm_loadu_ps(ptrD+8);
    _MM_TRANSPOSE4_PS(t1,t2,t3,t4);
    *z3           = t1;
    *x4           = t2;
    *y4           = t3;
    *z4           = t4;
}


static gmx_inline void
gmx_mm_decrement_1rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                       float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                       __m128 x1, __m128 y1, __m128 z1)
{
    __m128 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;
    t5          = _mm_unpacklo_ps(y1,z1);
    t6          = _mm_unpackhi_ps(y1,z1);
    t7          = _mm_shuffle_ps(x1,t5,_MM_SHUFFLE(1,0,0,0));
    t8          = _mm_shuffle_ps(x1,t5,_MM_SHUFFLE(3,2,0,1));
    t9          = _mm_shuffle_ps(x1,t6,_MM_SHUFFLE(1,0,0,2));
    t10         = _mm_shuffle_ps(x1,t6,_MM_SHUFFLE(3,2,0,3));
    t1          = _mm_load_ss(ptrA);
    t1          = _mm_loadh_pi(t1,(__m64 *)(ptrA+1));
    t1          = _mm_sub_ps(t1,t7);
    _mm_store_ss(ptrA,t1);
    _mm_storeh_pi((__m64 *)(ptrA+1),t1);
    t2          = _mm_load_ss(ptrB);
    t2          = _mm_loadh_pi(t2,(__m64 *)(ptrB+1));
    t2          = _mm_sub_ps(t2,t8);
    _mm_store_ss(ptrB,t2);
    _mm_storeh_pi((__m64 *)(ptrB+1),t2);
    t3          = _mm_load_ss(ptrC);
    t3          = _mm_loadh_pi(t3,(__m64 *)(ptrC+1));
    t3          = _mm_sub_ps(t3,t9);
    _mm_store_ss(ptrC,t3);
    _mm_storeh_pi((__m64 *)(ptrC+1),t3);
    t4          = _mm_load_ss(ptrD);
    t4          = _mm_loadh_pi(t4,(__m64 *)(ptrD+1));
    t4          = _mm_sub_ps(t4,t10);
    _mm_store_ss(ptrD,t4);
    _mm_storeh_pi((__m64 *)(ptrD+1),t4);
}


static gmx_inline void
gmx_mm_decrement_3rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                       float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                       __m128 x1, __m128 y1, __m128 z1,
                                       __m128 x2, __m128 y2, __m128 z2,
                                       __m128 x3, __m128 y3, __m128 z3)
{
    __m128 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    __m128 t11,t12,t13,t14,t15,t16,t17,t18,t19;
    __m128 t20,t21,t22,t23,t24,t25;
    t13         = _mm_unpackhi_ps(x1,y1);
    x1          = _mm_unpacklo_ps(x1,y1);
    t14         = _mm_unpackhi_ps(z1,x2);
    z1          = _mm_unpacklo_ps(z1,x2);
    t15         = _mm_unpackhi_ps(y2,z2);
    y2          = _mm_unpacklo_ps(y2,z2);
    t16         = _mm_unpackhi_ps(x3,y3);
    x3          = _mm_unpacklo_ps(x3,y3);
    t17         = _mm_permute_ps(z3,_MM_SHUFFLE(0,0,0,1));
    t18         = _mm_movehl_ps(z3,z3);
    t19         = _mm_permute_ps(t18,_MM_SHUFFLE(0,0,0,1));
    t20         = _mm_movelh_ps(x1,z1);
    t21         = _mm_movehl_ps(z1,x1);
    t22         = _mm_movelh_ps(t13,t14);
    t14         = _mm_movehl_ps(t14,t13);
    t23         = _mm_movelh_ps(y2,x3);
    t24         = _mm_movehl_ps(x3,y2);
    t25         = _mm_movelh_ps(t15,t16);
    t16         = _mm_movehl_ps(t16,t15);
    t1          = _mm_loadu_ps(ptrA);
    t2          = _mm_loadu_ps(ptrA+4);
    t3          = _mm_load_ss(ptrA+8);
    t1          = _mm_sub_ps(t1,t20);
    t2          = _mm_sub_ps(t2,t23);
    t3          = _mm_sub_ss(t3,z3);
    _mm_storeu_ps(ptrA,t1);
    _mm_storeu_ps(ptrA+4,t2);
    _mm_store_ss(ptrA+8,t3);
    t4          = _mm_loadu_ps(ptrB);
    t5          = _mm_loadu_ps(ptrB+4);
    t6          = _mm_load_ss(ptrB+8);
    t4          = _mm_sub_ps(t4,t21);
    t5          = _mm_sub_ps(t5,t24);
    t6          = _mm_sub_ss(t6,t17);
    _mm_storeu_ps(ptrB,t4);
    _mm_storeu_ps(ptrB+4,t5);
    _mm_store_ss(ptrB+8,t6);
    t7          = _mm_loadu_ps(ptrC);
    t8          = _mm_loadu_ps(ptrC+4);
    t9          = _mm_load_ss(ptrC+8);
    t7          = _mm_sub_ps(t7,t22);
    t8          = _mm_sub_ps(t8,t25);
    t9          = _mm_sub_ss(t9,t18);
    _mm_storeu_ps(ptrC,t7);
    _mm_storeu_ps(ptrC+4,t8);
    _mm_store_ss(ptrC+8,t9);
    t10         = _mm_loadu_ps(ptrD);
    t11         = _mm_loadu_ps(ptrD+4);
    t12         = _mm_load_ss(ptrD+8);
    t10         = _mm_sub_ps(t10,t14);
    t11         = _mm_sub_ps(t11,t16);
    t12         = _mm_sub_ss(t12,t19);
    _mm_storeu_ps(ptrD,t10);
    _mm_storeu_ps(ptrD+4,t11);
    _mm_store_ss(ptrD+8,t12);
}


static gmx_inline void
gmx_mm_decrement_4rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA, float * gmx_restrict ptrB,
                                       float * gmx_restrict ptrC, float * gmx_restrict ptrD,
                                       __m128 x1, __m128 y1, __m128 z1,
                                       __m128 x2, __m128 y2, __m128 z2,
                                       __m128 x3, __m128 y3, __m128 z3,
                                       __m128 x4, __m128 y4, __m128 z4)
{
    __m128 t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11;
    __m128 t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22;
    __m128 t23,t24;
    t13         = _mm_unpackhi_ps(x1,y1);
    x1          = _mm_unpacklo_ps(x1,y1);
    t14         = _mm_unpackhi_ps(z1,x2);
    z1          = _mm_unpacklo_ps(z1,x2);
    t15         = _mm_unpackhi_ps(y2,z2);
    y2          = _mm_unpacklo_ps(y2,z2);
    t16         = _mm_unpackhi_ps(x3,y3);
    x3          = _mm_unpacklo_ps(x3,y3);
    t17         = _mm_unpackhi_ps(z3,x4);
    z3          = _mm_unpacklo_ps(z3,x4);
    t18         = _mm_unpackhi_ps(y4,z4);
    y4          = _mm_unpacklo_ps(y4,z4);
    t19         = _mm_movelh_ps(x1,z1);
    z1          = _mm_movehl_ps(z1,x1);
    t20         = _mm_movelh_ps(t13,t14);
    t14         = _mm_movehl_ps(t14,t13);
    t21         = _mm_movelh_ps(y2,x3);
    x3          = _mm_movehl_ps(x3,y2);
    t22         = _mm_movelh_ps(t15,t16);
    t16         = _mm_movehl_ps(t16,t15);
    t23         = _mm_movelh_ps(z3,y4);
    y4          = _mm_movehl_ps(y4,z3);
    t24         = _mm_movelh_ps(t17,t18);
    t18         = _mm_movehl_ps(t18,t17);
    t1          = _mm_loadu_ps(ptrA);
    t2          = _mm_loadu_ps(ptrA+4);
    t3          = _mm_loadu_ps(ptrA+8);
    t1          = _mm_sub_ps(t1,t19);
    t2          = _mm_sub_ps(t2,t21);
    t3          = _mm_sub_ps(t3,t23);
    _mm_storeu_ps(ptrA,t1);
    _mm_storeu_ps(ptrA+4,t2);
    _mm_storeu_ps(ptrA+8,t3);
    t4          = _mm_loadu_ps(ptrB);
    t5          = _mm_loadu_ps(ptrB+4);
    t6          = _mm_loadu_ps(ptrB+8);
    t4          = _mm_sub_ps(t4,z1);
    t5          = _mm_sub_ps(t5,x3);
    t6          = _mm_sub_ps(t6,y4);
    _mm_storeu_ps(ptrB,t4);
    _mm_storeu_ps(ptrB+4,t5);
    _mm_storeu_ps(ptrB+8,t6);
    t7          = _mm_loadu_ps(ptrC);
    t8          = _mm_loadu_ps(ptrC+4);
    t9          = _mm_loadu_ps(ptrC+8);
    t7          = _mm_sub_ps(t7,t20);
    t8          = _mm_sub_ps(t8,t22);
    t9          = _mm_sub_ps(t9,t24);
    _mm_storeu_ps(ptrC,t7);
    _mm_storeu_ps(ptrC+4,t8);
    _mm_storeu_ps(ptrC+8,t9);
    t10         = _mm_loadu_ps(ptrD);
    t11         = _mm_loadu_ps(ptrD+4);
    t12         = _mm_loadu_ps(ptrD+8);
    t10         = _mm_sub_ps(t10,t14);
    t11         = _mm_sub_ps(t11,t16);
    t12         = _mm_sub_ps(t12,t18);
    _mm_storeu_ps(ptrD,t10);
    _mm_storeu_ps(ptrD+4,t11);
    _mm_storeu_ps(ptrD+8,t12);
}


static gmx_inline void
gmx_mm_update_iforce_1atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t2,t3;

    fix1 = _mm_hadd_ps(fix1,fix1);
    fiy1 = _mm_hadd_ps(fiy1,fiz1);

    fix1 = _mm_hadd_ps(fix1,fiy1); /* fiz1 fiy1 fix1 fix1 */

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

static gmx_inline void
gmx_mm_update_iforce_3atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      __m128 fix2, __m128 fiy2, __m128 fiz2,
                                      __m128 fix3, __m128 fiy3, __m128 fiz3,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t1,t2,t3,t4;

    fix1 = _mm_hadd_ps(fix1,fiy1);
    fiz1 = _mm_hadd_ps(fiz1,fix2);
    fiy2 = _mm_hadd_ps(fiy2,fiz2);
    fix3 = _mm_hadd_ps(fix3,fiy3);
    fiz3 = _mm_hadd_ps(fiz3,fiz3);

    fix1 = _mm_hadd_ps(fix1,fiz1); /* fix2 fiz1 fiy1 fix1 */
    fiy2 = _mm_hadd_ps(fiy2,fix3); /* fiy3 fix3 fiz2 fiy2 */
    fiz3 = _mm_hadd_ps(fiz3,fiz3); /*  -    -    -   fiz3 */

    _mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
    _mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));
    _mm_store_ss (fptr+8,_mm_add_ss(fiz3,_mm_load_ss(fptr+8) ));

    t4 = _mm_load_ss(fshiftptr+2);
    t4 = _mm_loadh_pi(t4,(__m64 *)(fshiftptr));

    t1 = _mm_shuffle_ps(fiz3,fix1,_MM_SHUFFLE(1,0,0,0));   /* fiy1 fix1  -   fiz3 */
    t2 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(3,2,2,2));   /* fiy3 fix3  -   fiz1 */
    t3 = _mm_shuffle_ps(fiy2,fix1,_MM_SHUFFLE(3,3,0,1));   /* fix2 fix2 fiy2 fiz2 */
    t3 = _mm_permute_ps(t3  ,_MM_SHUFFLE(1,2,0,0));        /* fiy2 fix2  -   fiz2 */

    t1 = _mm_add_ps(t1,t2);
    t3 = _mm_add_ps(t3,t4);
    t1 = _mm_add_ps(t1,t3); /* y x - z */

    _mm_store_ss(fshiftptr+2,t1);
    _mm_storeh_pi((__m64 *)(fshiftptr),t1);
}


static gmx_inline void
gmx_mm_update_iforce_4atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      __m128 fix2, __m128 fiy2, __m128 fiz2,
                                      __m128 fix3, __m128 fiy3, __m128 fiz3,
                                      __m128 fix4, __m128 fiy4, __m128 fiz4,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t1,t2,t3,t4,t5;

    fix1 = _mm_hadd_ps(fix1,fiy1);
    fiz1 = _mm_hadd_ps(fiz1,fix2);
    fiy2 = _mm_hadd_ps(fiy2,fiz2);
    fix3 = _mm_hadd_ps(fix3,fiy3);
    fiz3 = _mm_hadd_ps(fiz3,fix4);
    fiy4 = _mm_hadd_ps(fiy4,fiz4);

    fix1 = _mm_hadd_ps(fix1,fiz1); /* fix2 fiz1 fiy1 fix1 */
    fiy2 = _mm_hadd_ps(fiy2,fix3); /* fiy3 fix3 fiz2 fiy2 */
    fiz3 = _mm_hadd_ps(fiz3,fiy4); /* fiz4 fiy4 fix4 fiz3 */

    _mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
    _mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));
    _mm_storeu_ps(fptr+8,_mm_add_ps(fiz3,_mm_loadu_ps(fptr+8)));

    t5 = _mm_load_ss(fshiftptr+2);
    t5 = _mm_loadh_pi(t5,(__m64 *)(fshiftptr));

    t1 = _mm_permute_ps(fix1,_MM_SHUFFLE(1,0,2,2));
    t2 = _mm_permute_ps(fiy2,_MM_SHUFFLE(3,2,1,1));
    t3 = _mm_permute_ps(fiz3,_MM_SHUFFLE(2,1,0,0));
    t4 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(0,0,3,3));
    t4 = _mm_shuffle_ps(fiz3,t4  ,_MM_SHUFFLE(2,0,3,3));

    t1 = _mm_add_ps(t1,t2);
    t3 = _mm_add_ps(t3,t4);
    t1 = _mm_add_ps(t1,t3);
    t5 = _mm_add_ps(t5,t1);

    _mm_store_ss(fshiftptr+2,t5);
    _mm_storeh_pi((__m64 *)(fshiftptr),t5);
}



static gmx_inline void
gmx_mm_update_1pot_ps(__m128 pot1, float * gmx_restrict ptrA)
{
    pot1 = _mm_hadd_ps(pot1,pot1);
    pot1 = _mm_hadd_ps(pot1,pot1);
    _mm_store_ss(ptrA,_mm_add_ss(pot1,_mm_load_ss(ptrA)));
}

static gmx_inline void
gmx_mm_update_2pot_ps(__m128 pot1, float * gmx_restrict ptrA,
                      __m128 pot2, float * gmx_restrict ptrB)
{
    pot1 = _mm_hadd_ps(pot1,pot2);
    pot1 = _mm_hadd_ps(pot1,pot1);
    pot2 = _mm_permute_ps(pot1,_MM_SHUFFLE(0,0,0,1));
    _mm_store_ss(ptrA,_mm_add_ss(pot1,_mm_load_ss(ptrA)));
    _mm_store_ss(ptrB,_mm_add_ss(pot2,_mm_load_ss(ptrB)));
}


static gmx_inline void
gmx_mm_update_4pot_ps(__m128 pot1, float * gmx_restrict ptrA,
                      __m128 pot2, float * gmx_restrict ptrB,
                      __m128 pot3, float * gmx_restrict ptrC,
                      __m128 pot4, float * gmx_restrict ptrD)
{
    _MM_TRANSPOSE4_PS(pot1,pot2,pot3,pot4);
    pot1 = _mm_add_ps(_mm_add_ps(pot1,pot2),_mm_add_ps(pot3,pot4));
    pot2 = _mm_permute_ps(pot1,_MM_SHUFFLE(1,1,1,1));
    pot3 = _mm_permute_ps(pot1,_MM_SHUFFLE(2,2,2,2));
    pot4 = _mm_permute_ps(pot1,_MM_SHUFFLE(3,3,3,3));
    _mm_store_ss(ptrA,_mm_add_ss(pot1,_mm_load_ss(ptrA)));
    _mm_store_ss(ptrB,_mm_add_ss(pot2,_mm_load_ss(ptrB)));
    _mm_store_ss(ptrC,_mm_add_ss(pot3,_mm_load_ss(ptrC)));
    _mm_store_ss(ptrD,_mm_add_ss(pot4,_mm_load_ss(ptrD)));
}


#endif /* _kernelutil_x86_avx_128_fma_single_h_ */
