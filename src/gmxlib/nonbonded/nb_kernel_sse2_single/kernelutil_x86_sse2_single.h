/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011-2012, The GROMACS Development Team
 * Copyright (c) 2012, by the GROMACS development team, led by
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
#ifndef _kernelutil_x86_sse2_single_h_
#define _kernelutil_x86_sse2_single_h_

/* We require SSE2 now! */

#include <math.h>

#include "gmx_x86_sse2.h"


/* Normal sum of four xmm registers */
#define gmx_mm_sum4_ps(t0,t1,t2,t3)  _mm_add_ps(_mm_add_ps(t0,t1),_mm_add_ps(t2,t3))

static gmx_inline __m128
gmx_mm_calc_rsq_ps(__m128 dx, __m128 dy, __m128 dz)
{
    return _mm_add_ps( _mm_add_ps( _mm_mul_ps(dx,dx), _mm_mul_ps(dy,dy) ), _mm_mul_ps(dz,dz) );
}

static int
gmx_mm_any_lt(__m128 a, __m128 b)
{
    return _mm_movemask_ps(_mm_cmplt_ps(a,b));
}

/* Load a single value from 1-4 places, merge into xmm register */

static __m128
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

static void
gmx_mm_store_4real_swizzle_ps(float * gmx_restrict ptrA,
                              float * gmx_restrict ptrB,
                              float * gmx_restrict ptrC,
                              float * gmx_restrict ptrD,
                              __m128 xmm1)
{
    __m128 t2,t3,t4;

    t3       = _mm_movehl_ps(_mm_setzero_ps(),xmm1);
    t2       = _mm_shuffle_ps(xmm1,xmm1,_MM_SHUFFLE(1,1,1,1));
    t4       = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(1,1,1,1));
    _mm_store_ss(ptrA,xmm1);
    _mm_store_ss(ptrB,t2);
    _mm_store_ss(ptrC,t3);
    _mm_store_ss(ptrD,t4);
}

/* Similar to store, but increments value in memory */
static void
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


static void
gmx_mm_load_4pair_swizzle_ps(const float * gmx_restrict p1,
                             const float * gmx_restrict p2,
                             const float * gmx_restrict p3,
                             const float * gmx_restrict p4,
                             __m128 * gmx_restrict c6,
                             __m128 * gmx_restrict c12)
{
    __m128 t1,t2,t3,t4;

    t1   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p1);   /* - - c12a  c6a */
    t2   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p2);   /* - - c12b  c6b */
    t3   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p3);   /* - - c12c  c6c */
    t4   = _mm_loadl_pi(_mm_setzero_ps(),(__m64 *)p4);   /* - - c12d  c6d */
    t1   = _mm_unpacklo_ps(t1,t2);
    t2   = _mm_unpacklo_ps(t3,t4);
    *c6  = _mm_movelh_ps(t1,t2);
    *c12 = _mm_movehl_ps(t2,t1);
}

/* Routines to load 1-4 rvec from 4 places.
 * We mainly use these to load coordinates. The extra routines
 * are very efficient for the water-water loops, since we e.g.
 * know that a TIP4p water has 4 atoms, so we should load 12 floats+shuffle.
 */


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

    *x1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(0,0,0,0));
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
    t4   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(0,2,1,0));
    t5   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(1,0,2,1));
    t6   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(2,1,0,2));

    t1   = _mm_add_ps(t1,t4);
    t2   = _mm_add_ps(t2,t5);
    t3   = _mm_add_ss(t3,t6);

    *x1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(2,2,2,2));
    *x2  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,3,3,3));
    *y2  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(0,0,0,0));
    *z2  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(1,1,1,1));
    *x3  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(2,2,2,2));
    *y3  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(3,3,3,3));
    *z3  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(0,0,0,0));
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

    tA   = _mm_castpd_ps(_mm_load_sd((const double *)xyz_shift));
    tB   = _mm_load_ss(xyz_shift+2);

    t1   = _mm_loadu_ps(xyz);
    t2   = _mm_loadu_ps(xyz+4);
    t3   = _mm_loadu_ps(xyz+8);

    tA   = _mm_movelh_ps(tA,tB);
    t4   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(0,2,1,0));
    t5   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(1,0,2,1));
    t6   = _mm_shuffle_ps(tA,tA,_MM_SHUFFLE(2,1,0,2));

    t1   = _mm_add_ps(t1,t4);
    t2   = _mm_add_ps(t2,t5);
    t3   = _mm_add_ps(t3,t6);

    *x1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(0,0,0,0));
    *y1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(1,1,1,1));
    *z1  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(2,2,2,2));
    *x2  = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,3,3,3));
    *y2  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(0,0,0,0));
    *z2  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(1,1,1,1));
    *x3  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(2,2,2,2));
    *y3  = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(3,3,3,3));
    *z3  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(0,0,0,0));
    *x4  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(1,1,1,1));
    *y4  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(2,2,2,2));
    *z4  = _mm_shuffle_ps(t3,t3,_MM_SHUFFLE(3,3,3,3));
}


static void
gmx_mm_load_1rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA,
                                  const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC,
                                  const float * gmx_restrict ptrD,
                                  __m128 *      gmx_restrict x1,
                                  __m128 *      gmx_restrict y1,
                                  __m128 *      gmx_restrict z1)
{
    __m128 t1,t2,t3,t4,t5,t6,t7,t8;
    t1   = _mm_castpd_ps(_mm_load_sd((const double *)ptrA));
    t2   = _mm_castpd_ps(_mm_load_sd((const double *)ptrB));
    t3   = _mm_castpd_ps(_mm_load_sd((const double *)ptrC));
    t4   = _mm_castpd_ps(_mm_load_sd((const double *)ptrD));
    t5 = _mm_load_ss(ptrA+2);
    t6 = _mm_load_ss(ptrB+2);
    t7 = _mm_load_ss(ptrC+2);
    t8 = _mm_load_ss(ptrD+2);
    t1 = _mm_unpacklo_ps(t1,t2);
    t3 = _mm_unpacklo_ps(t3,t4);
    *x1 = _mm_movelh_ps(t1,t3);
    *y1 = _mm_movehl_ps(t3,t1);
    t5  = _mm_unpacklo_ps(t5,t6);
    t7  = _mm_unpacklo_ps(t7,t8);
    *z1 = _mm_movelh_ps(t5,t7);
}


static void
gmx_mm_load_3rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA,
                                  const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC,
                                  const float * gmx_restrict ptrD,
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


static void
gmx_mm_load_4rvec_4ptr_swizzle_ps(const float * gmx_restrict ptrA,
                                  const float * gmx_restrict ptrB,
                                  const float * gmx_restrict ptrC,
                                  const float * gmx_restrict ptrD,
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


static void
gmx_mm_decrement_1rvec_4ptr_swizzle_ps(float * gmx_restrict ptrA,
                                       float * gmx_restrict ptrB,
                                       float * gmx_restrict ptrC,
                                       float * gmx_restrict ptrD,
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



#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm_decrement_3rvec_4ptr_swizzle_ps(ptrA,ptrB,ptrC,ptrD, \
_x1,_y1,_z1,_x2,_y2,_z2,_x3,_y3,_z3) \
{\
__m128 _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10;\
__m128 _t11,_t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19;\
__m128 _t20,_t21,_t22,_t23,_t24,_t25;\
_t13         = _mm_unpackhi_ps(_x1,_y1);\
_x1          = _mm_unpacklo_ps(_x1,_y1);\
_t14         = _mm_unpackhi_ps(_z1,_x2);\
_z1          = _mm_unpacklo_ps(_z1,_x2);\
_t15         = _mm_unpackhi_ps(_y2,_z2);\
_y2          = _mm_unpacklo_ps(_y2,_z2);\
_t16         = _mm_unpackhi_ps(_x3,_y3);\
_x3          = _mm_unpacklo_ps(_x3,_y3);\
_t17         = _mm_shuffle_ps(_z3,_z3,_MM_SHUFFLE(0,0,0,1));\
_t18         = _mm_movehl_ps(_z3,_z3);\
_t19         = _mm_shuffle_ps(_t18,_t18,_MM_SHUFFLE(0,0,0,1));\
_t20         = _mm_movelh_ps(_x1,_z1);\
_t21         = _mm_movehl_ps(_z1,_x1);\
_t22         = _mm_movelh_ps(_t13,_t14);\
_t14         = _mm_movehl_ps(_t14,_t13);\
_t23         = _mm_movelh_ps(_y2,_x3);\
_t24         = _mm_movehl_ps(_x3,_y2);\
_t25         = _mm_movelh_ps(_t15,_t16);\
_t16         = _mm_movehl_ps(_t16,_t15);\
_t1          = _mm_loadu_ps(ptrA);\
_t2          = _mm_loadu_ps(ptrA+4);\
_t3          = _mm_load_ss(ptrA+8);\
_t1          = _mm_sub_ps(_t1,_t20);\
_t2          = _mm_sub_ps(_t2,_t23);\
_t3          = _mm_sub_ss(_t3,_z3);\
_mm_storeu_ps(ptrA,_t1);\
_mm_storeu_ps(ptrA+4,_t2);\
_mm_store_ss(ptrA+8,_t3);\
_t4          = _mm_loadu_ps(ptrB);\
_t5          = _mm_loadu_ps(ptrB+4);\
_t6          = _mm_load_ss(ptrB+8);\
_t4          = _mm_sub_ps(_t4,_t21);\
_t5          = _mm_sub_ps(_t5,_t24);\
_t6          = _mm_sub_ss(_t6,_t17);\
_mm_storeu_ps(ptrB,_t4);\
_mm_storeu_ps(ptrB+4,_t5);\
_mm_store_ss(ptrB+8,_t6);\
_t7          = _mm_loadu_ps(ptrC);\
_t8          = _mm_loadu_ps(ptrC+4);\
_t9          = _mm_load_ss(ptrC+8);\
_t7          = _mm_sub_ps(_t7,_t22);\
_t8          = _mm_sub_ps(_t8,_t25);\
_t9          = _mm_sub_ss(_t9,_t18);\
_mm_storeu_ps(ptrC,_t7);\
_mm_storeu_ps(ptrC+4,_t8);\
_mm_store_ss(ptrC+8,_t9);\
_t10         = _mm_loadu_ps(ptrD);\
_t11         = _mm_loadu_ps(ptrD+4);\
_t12         = _mm_load_ss(ptrD+8);\
_t10         = _mm_sub_ps(_t10,_t14);\
_t11         = _mm_sub_ps(_t11,_t16);\
_t12         = _mm_sub_ss(_t12,_t19);\
_mm_storeu_ps(ptrD,_t10);\
_mm_storeu_ps(ptrD+4,_t11);\
_mm_store_ss(ptrD+8,_t12);\
}
#else
/* Real function for sane compilers */
static void
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
    t17         = _mm_shuffle_ps(z3,z3,_MM_SHUFFLE(0,0,0,1));
    t18         = _mm_movehl_ps(z3,z3);
    t19         = _mm_shuffle_ps(t18,t18,_MM_SHUFFLE(0,0,0,1));
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
#endif


#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm_decrement_4rvec_4ptr_swizzle_ps(ptrA,ptrB,ptrC,ptrD, \
_x1,_y1,_z1,_x2,_y2,_z2,_x3,_y3,_z3,_x4,_y4,_z4) \
{\
__m128 _t1,_t2,_t3,_t4,_t5,_t6,_t7,_t8,_t9,_t10,_t11;\
__m128 _t12,_t13,_t14,_t15,_t16,_t17,_t18,_t19,_t20,_t21,_t22;\
__m128 _t23,_t24;\
_t13         = _mm_unpackhi_ps(_x1,_y1);\
_x1          = _mm_unpacklo_ps(_x1,_y1);\
_t14         = _mm_unpackhi_ps(_z1,_x2);\
_z1          = _mm_unpacklo_ps(_z1,_x2);\
_t15         = _mm_unpackhi_ps(_y2,_z2);\
_y2          = _mm_unpacklo_ps(_y2,_z2);\
_t16         = _mm_unpackhi_ps(_x3,_y3);\
_x3          = _mm_unpacklo_ps(_x3,_y3);\
_t17         = _mm_unpackhi_ps(_z3,_x4);\
_z3          = _mm_unpacklo_ps(_z3,_x4);\
_t18         = _mm_unpackhi_ps(_y4,_z4);\
_y4          = _mm_unpacklo_ps(_y4,_z4);\
_t19         = _mm_movelh_ps(_x1,_z1);\
_z1          = _mm_movehl_ps(_z1,_x1);\
_t20         = _mm_movelh_ps(_t13,_t14);\
_t14         = _mm_movehl_ps(_t14,_t13);\
_t21         = _mm_movelh_ps(_y2,_x3);\
_x3          = _mm_movehl_ps(_x3,_y2);\
_t22         = _mm_movelh_ps(_t15,_t16);\
_t16         = _mm_movehl_ps(_t16,_t15);\
_t23         = _mm_movelh_ps(_z3,_y4);\
_y4          = _mm_movehl_ps(_y4,_z3);\
_t24         = _mm_movelh_ps(_t17,_t18);\
_t18         = _mm_movehl_ps(_t18,_t17);\
_t1          = _mm_loadu_ps(ptrA);\
_t2          = _mm_loadu_ps(ptrA+4);\
_t3          = _mm_loadu_ps(ptrA+8);\
_t1          = _mm_sub_ps(_t1,_t19);\
_t2          = _mm_sub_ps(_t2,_t21);\
_t3          = _mm_sub_ps(_t3,_t23);\
_mm_storeu_ps(ptrA,_t1);\
_mm_storeu_ps(ptrA+4,_t2);\
_mm_storeu_ps(ptrA+8,_t3);\
_t4          = _mm_loadu_ps(ptrB);\
_t5          = _mm_loadu_ps(ptrB+4);\
_t6          = _mm_loadu_ps(ptrB+8);\
_t4          = _mm_sub_ps(_t4,_z1);\
_t5          = _mm_sub_ps(_t5,_x3);\
_t6          = _mm_sub_ps(_t6,_y4);\
_mm_storeu_ps(ptrB,_t4);\
_mm_storeu_ps(ptrB+4,_t5);\
_mm_storeu_ps(ptrB+8,_t6);\
_t7          = _mm_loadu_ps(ptrC);\
_t8          = _mm_loadu_ps(ptrC+4);\
_t9          = _mm_loadu_ps(ptrC+8);\
_t7          = _mm_sub_ps(_t7,_t20);\
_t8          = _mm_sub_ps(_t8,_t22);\
_t9          = _mm_sub_ps(_t9,_t24);\
_mm_storeu_ps(ptrC,_t7);\
_mm_storeu_ps(ptrC+4,_t8);\
_mm_storeu_ps(ptrC+8,_t9);\
_t10         = _mm_loadu_ps(ptrD);\
_t11         = _mm_loadu_ps(ptrD+4);\
_t12         = _mm_loadu_ps(ptrD+8);\
_t10         = _mm_sub_ps(_t10,_t14);\
_t11         = _mm_sub_ps(_t11,_t16);\
_t12         = _mm_sub_ps(_t12,_t18);\
_mm_storeu_ps(ptrD,_t10);\
_mm_storeu_ps(ptrD+4,_t11);\
_mm_storeu_ps(ptrD+8,_t12);\
}
#else
/* Real function for sane compilers */
static void
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
#endif


static gmx_inline void
gmx_mm_update_iforce_1atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t1,t2,t3;

    /* transpose data */
    t1 = fix1;
    _MM_TRANSPOSE4_PS(fix1,t1,fiy1,fiz1);
    fix1 = _mm_add_ps(_mm_add_ps(fix1,t1), _mm_add_ps(fiy1,fiz1));

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

#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm_update_iforce_3atom_swizzle_ps(fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3, \
                                              fptr,fshiftptr) \
{\
    __m128 _t1,_t2,_t3,_t4;\
\
    _MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);\
    _MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);\
    _t2   = _mm_movehl_ps(_mm_setzero_ps(),fiz3);\
    _t1   = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(0,0,0,1));\
    _t3   = _mm_shuffle_ps(_t2,_t2,_MM_SHUFFLE(0,0,0,1));\
    fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));\
    fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));\
    fiz3 = _mm_add_ss(_mm_add_ps(fiz3,_t1)  , _mm_add_ps(_t2,_t3));\
    _mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));\
    _mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));\
    _mm_store_ss (fptr+8,_mm_add_ss(fiz3,_mm_load_ss(fptr+8) ));\
    _t4 = _mm_load_ss(fshiftptr+2);\
    _t4 = _mm_loadh_pi(_t4,(__m64 *)(fshiftptr));\
    _t1 = _mm_shuffle_ps(fiz3,fix1,_MM_SHUFFLE(1,0,0,0));\
    _t2 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(3,2,2,2));\
    _t3 = _mm_shuffle_ps(fiy2,fix1,_MM_SHUFFLE(3,3,0,1));\
    _t3 = _mm_shuffle_ps(_t3  ,_t3  ,_MM_SHUFFLE(1,2,0,0));\
    _t1 = _mm_add_ps(_t1,_t2);\
    _t3 = _mm_add_ps(_t3,_t4);\
    _t1 = _mm_add_ps(_t1,_t3);\
    _mm_store_ss(fshiftptr+2,_t1);\
    _mm_storeh_pi((__m64 *)(fshiftptr),_t1);\
}
#else
/* Real function for sane compilers */
static gmx_inline void
gmx_mm_update_iforce_3atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      __m128 fix2, __m128 fiy2, __m128 fiz2,
                                      __m128 fix3, __m128 fiy3, __m128 fiz3,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t1,t2,t3,t4;

    /* transpose data */
    _MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);
    _MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);
    t2   = _mm_movehl_ps(_mm_setzero_ps(),fiz3);
    t1   = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(0,0,0,1));
    t3   = _mm_shuffle_ps(t2,t2,_MM_SHUFFLE(0,0,0,1));

    fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));
    fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));
    fiz3 = _mm_add_ss(_mm_add_ps(fiz3,t1)  , _mm_add_ps(t2,t3));

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
#endif

#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm_update_iforce_4atom_swizzle_ps(fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3,fix4,fiy4,fiz4, \
                                              fptr,fshiftptr) \
{\
    __m128 _t1,_t2,_t3,_t4,_t5;\
    _MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);\
    _MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);\
    _MM_TRANSPOSE4_PS(fiz3,fix4,fiy4,fiz4);\
    fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));\
    fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));\
    fiz3 = _mm_add_ps(_mm_add_ps(fiz3,fix4), _mm_add_ps(fiy4,fiz4));\
    _mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));\
    _mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));\
    _mm_storeu_ps(fptr+8,_mm_add_ps(fiz3,_mm_loadu_ps(fptr+8)));\
    _t5 = _mm_load_ss(fshiftptr+2);\
    _t5 = _mm_loadh_pi(_t5,(__m64 *)(fshiftptr));\
    _t1 = _mm_shuffle_ps(fix1,fix1,_MM_SHUFFLE(1,0,2,2));\
    _t2 = _mm_shuffle_ps(fiy2,fiy2,_MM_SHUFFLE(3,2,1,1));\
    _t3 = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(2,1,0,0));\
    _t4 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(0,0,3,3));\
    _t4 = _mm_shuffle_ps(fiz3,_t4  ,_MM_SHUFFLE(2,0,3,3));\
    _t1 = _mm_add_ps(_t1,_t2);\
    _t3 = _mm_add_ps(_t3,_t4);\
    _t1 = _mm_add_ps(_t1,_t3);\
    _t5 = _mm_add_ps(_t5,_t1);\
    _mm_store_ss(fshiftptr+2,_t5);\
    _mm_storeh_pi((__m64 *)(fshiftptr),_t5);\
}
#else
/* Real function for sane compilers */
static gmx_inline void
gmx_mm_update_iforce_4atom_swizzle_ps(__m128 fix1, __m128 fiy1, __m128 fiz1,
                                      __m128 fix2, __m128 fiy2, __m128 fiz2,
                                      __m128 fix3, __m128 fiy3, __m128 fiz3,
                                      __m128 fix4, __m128 fiy4, __m128 fiz4,
                                      float * gmx_restrict fptr,
                                      float * gmx_restrict fshiftptr)
{
    __m128 t1,t2,t3,t4,t5;

    /* transpose data */
    _MM_TRANSPOSE4_PS(fix1,fiy1,fiz1,fix2);
    _MM_TRANSPOSE4_PS(fiy2,fiz2,fix3,fiy3);
    _MM_TRANSPOSE4_PS(fiz3,fix4,fiy4,fiz4);

    fix1 = _mm_add_ps(_mm_add_ps(fix1,fiy1), _mm_add_ps(fiz1,fix2));
    fiy2 = _mm_add_ps(_mm_add_ps(fiy2,fiz2), _mm_add_ps(fix3,fiy3));
    fiz3 = _mm_add_ps(_mm_add_ps(fiz3,fix4), _mm_add_ps(fiy4,fiz4));

    _mm_storeu_ps(fptr,  _mm_add_ps(fix1,_mm_loadu_ps(fptr)  ));
    _mm_storeu_ps(fptr+4,_mm_add_ps(fiy2,_mm_loadu_ps(fptr+4)));
    _mm_storeu_ps(fptr+8,_mm_add_ps(fiz3,_mm_loadu_ps(fptr+8)));

    t5 = _mm_load_ss(fshiftptr+2);
    t5 = _mm_loadh_pi(t5,(__m64 *)(fshiftptr));

    t1 = _mm_shuffle_ps(fix1,fix1,_MM_SHUFFLE(1,0,2,2));
    t2 = _mm_shuffle_ps(fiy2,fiy2,_MM_SHUFFLE(3,2,1,1));
    t3 = _mm_shuffle_ps(fiz3,fiz3,_MM_SHUFFLE(2,1,0,0));
    t4 = _mm_shuffle_ps(fix1,fiy2,_MM_SHUFFLE(0,0,3,3));
    t4 = _mm_shuffle_ps(fiz3,t4  ,_MM_SHUFFLE(2,0,3,3));

    t1 = _mm_add_ps(t1,t2);
    t3 = _mm_add_ps(t3,t4);
    t1 = _mm_add_ps(t1,t3);
    t5 = _mm_add_ps(t5,t1);

    _mm_store_ss(fshiftptr+2,t5);
    _mm_storeh_pi((__m64 *)(fshiftptr),t5);
}
#endif


static void
gmx_mm_update_1pot_ps(__m128 pot1, float * gmx_restrict ptrA)
{
    pot1 = _mm_add_ps(pot1,_mm_movehl_ps(_mm_setzero_ps(),pot1));
    pot1 = _mm_add_ps(pot1,_mm_shuffle_ps(pot1,pot1,_MM_SHUFFLE(0,0,0,1)));
    _mm_store_ss(ptrA,_mm_add_ss(pot1,_mm_load_ss(ptrA)));
}

static void
gmx_mm_update_2pot_ps(__m128 pot1, float * gmx_restrict ptrA,
                      __m128 pot2, float * gmx_restrict ptrB)
{
    __m128 t1,t2;
    t1   = _mm_movehl_ps(pot2,pot1);
    t2   = _mm_movelh_ps(pot1,pot2);
    t1   = _mm_add_ps(t1,t2);
    t2   = _mm_shuffle_ps(t1,t1,_MM_SHUFFLE(3,3,1,1));
    pot1 = _mm_add_ps(t1,t2);
    pot2 = _mm_movehl_ps(t2,pot1);
    _mm_store_ss(ptrA,_mm_add_ss(pot1,_mm_load_ss(ptrA)));
    _mm_store_ss(ptrB,_mm_add_ss(pot2,_mm_load_ss(ptrB)));
}


#endif /* _kernelutil_x86_sse2_single_h_ */
