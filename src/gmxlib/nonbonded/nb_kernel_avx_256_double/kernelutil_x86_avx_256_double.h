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
#ifndef _kernelutil_x86_avx_256_double_h_
#define _kernelutil_x86_avx_256_double_h_


#include "gmx_x86_avx_256.h"


static int
gmx_mm256_any_lt(__m256d a, __m256d b)
{
    return _mm256_movemask_pd(_mm256_cmp_pd(a, b, _CMP_LT_OQ));
}

static gmx_inline __m256d
gmx_mm256_calc_rsq_pd(__m256d dx, __m256d dy, __m256d dz)
{
    return _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy) ), _mm256_mul_pd(dz, dz) );
}

/* Normal sum of four ymm registers */
#define gmx_mm256_sum4_pd(t0, t1, t2, t3)  _mm256_add_pd(_mm256_add_pd(t0, t1), _mm256_add_pd(t2, t3))


/* Load a single value from 1-4 places, merge into xmm register */
static __m256d
gmx_mm256_load_1real_pd(const double * gmx_restrict ptrA)
{
    return _mm256_castpd128_pd256(_mm_load_sd(ptrA));
}

static __m256d
gmx_mm256_load_2real_swizzle_pd(const double * gmx_restrict ptrA,
                                const double * gmx_restrict ptrB)
{
    __m128d tA, tB;

    tA = _mm_load_sd(ptrA);
    tB = _mm_load_sd(ptrB);

    return _mm256_castpd128_pd256(_mm_unpacklo_pd(tA, tB));
}


static __m256d
gmx_mm256_load_4real_swizzle_pd(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                const double * gmx_restrict ptrC, const double * gmx_restrict ptrD)
{
    __m128d t1, t2;

    t1 = _mm_unpacklo_pd(_mm_load_sd(ptrA), _mm_load_sd(ptrB));
    t2 = _mm_unpacklo_pd(_mm_load_sd(ptrC), _mm_load_sd(ptrD));
    return gmx_mm256_set_m128d(t2, t1);
}



static void
gmx_mm256_store_1real_pd(double * gmx_restrict ptrA, __m256d xmm1)
{
    _mm_store_sd(ptrA, _mm256_castpd256_pd128(xmm1));
}


static void
gmx_mm256_store_2real_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB, __m256d xmm1)
{
    __m256d t2;

    t2       = _mm256_permute_pd(xmm1, _GMX_MM_PERMUTE256D(1, 1, 1, 1));
    _mm_store_sd(ptrA, _mm256_castpd256_pd128(xmm1));
    _mm_store_sd(ptrB, _mm256_castpd256_pd128(t2));
}




static void
gmx_mm256_store_4real_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                 double * gmx_restrict ptrC, double * gmx_restrict ptrD, __m256d xmm1)
{
    __m256d t2;
    __m128d t3, t4;

    t2       = _mm256_permute_pd(xmm1, _GMX_MM_PERMUTE256D(1, 1, 1, 1));
    t3       = _mm256_extractf128_pd(xmm1, 0x1);
    t4       = _mm_permute_pd(t3, _GMX_MM_PERMUTE128D(1, 1));
    _mm_store_sd(ptrA, _mm256_castpd256_pd128(xmm1));
    _mm_store_sd(ptrB, _mm256_castpd256_pd128(t2));
    _mm_store_sd(ptrC, t3);
    _mm_store_sd(ptrD, t4);
}




static void
gmx_mm256_increment_1real_pd(double * gmx_restrict ptrA, __m256d xmm1)
{
    __m128d t1;

    t1   = _mm256_castpd256_pd128(xmm1);
    t1   = _mm_add_sd(t1, _mm_load_sd(ptrA));

    _mm_store_sd(ptrA, t1);
}


static void
gmx_mm256_increment_2real_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB, __m256d xmm1)
{
    __m128d t1, t2;

    t1   = _mm256_castpd256_pd128(xmm1);
    t2   = _mm_permute_pd(t1, _GMX_MM_PERMUTE128D(1, 1));

    t1   = _mm_add_sd(t1, _mm_load_sd(ptrA));
    t2   = _mm_add_sd(t2, _mm_load_sd(ptrB));

    _mm_store_sd(ptrA, t1);
    _mm_store_sd(ptrB, t2);
}


static void
gmx_mm256_increment_4real_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                     double * gmx_restrict ptrC, double * gmx_restrict ptrD, __m256d xmm1)
{
    __m128d t1, t2, t3, t4;

    t1   = _mm256_castpd256_pd128(xmm1);
    t2   = _mm_permute_pd(t1, _GMX_MM_PERMUTE128D(1, 1));
    t3   = _mm256_extractf128_pd(xmm1, 0x1);
    t4   = _mm_permute_pd(t3, _GMX_MM_PERMUTE128D(1, 1));

    t1   = _mm_add_sd(t1, _mm_load_sd(ptrA));
    t2   = _mm_add_sd(t2, _mm_load_sd(ptrB));
    t3   = _mm_add_sd(t3, _mm_load_sd(ptrC));
    t4   = _mm_add_sd(t4, _mm_load_sd(ptrD));

    _mm_store_sd(ptrA, t1);
    _mm_store_sd(ptrB, t2);
    _mm_store_sd(ptrC, t3);
    _mm_store_sd(ptrD, t4);
}



static void
gmx_mm256_load_1pair_swizzle_pd(const double * gmx_restrict p1, __m256d *c6, __m256d *c12)
{
    *c6     = _mm256_castpd128_pd256(_mm_load_sd(p1));
    *c12    = _mm256_castpd128_pd256(_mm_load_sd(p1+1));
}


static void
gmx_mm256_load_2pair_swizzle_pd(const double * gmx_restrict p1, const double * gmx_restrict p2, __m256d *c6, __m256d *c12)
{
    __m128d t1, t2, t3;

    t1   = _mm_loadu_pd(p1);
    t2   = _mm_loadu_pd(p2);
    *c6  = _mm256_castpd128_pd256(_mm_unpacklo_pd(t1, t2));
    *c12 = _mm256_castpd128_pd256(_mm_unpackhi_pd(t1, t2));
}



static void
gmx_mm256_load_4pair_swizzle_pd(const double * gmx_restrict p1, const double * gmx_restrict p2,
                                const double * gmx_restrict p3, const double * gmx_restrict p4,
                                __m256d * gmx_restrict c6, __m256d * gmx_restrict c12)
{
    __m256d t1, t2;

    t1   = gmx_mm256_set_m128d(_mm_loadu_pd(p3), _mm_loadu_pd(p1)); /* c12c  c6c | c12a  c6a */
    t2   = gmx_mm256_set_m128d(_mm_loadu_pd(p4), _mm_loadu_pd(p2)); /* c12d  c6d | c12b  c6b */

    *c6  = _mm256_unpacklo_pd(t1, t2);                              /* c6d c6c | c6b c6a */
    *c12 = _mm256_unpackhi_pd(t1, t2);                              /* c12d c12c | c12b c12a */
}


static gmx_inline void
gmx_mm256_load_shift_and_1rvec_broadcast_pd(const double * gmx_restrict xyz_shift,
                                            const double * gmx_restrict xyz,
                                            __m256d * gmx_restrict      x1,
                                            __m256d * gmx_restrict      y1,
                                            __m256d * gmx_restrict      z1)
{
    __m128d mem_xy, mem_z, mem_sxy, mem_sz, tx, ty, tz;

    mem_xy  = _mm_loadu_pd(xyz);
    mem_z   = _mm_load_sd(xyz+2);
    mem_sxy = _mm_loadu_pd(xyz_shift);
    mem_sz  = _mm_load_sd(xyz_shift+2);

    mem_xy  = _mm_add_pd(mem_xy, mem_sxy);
    mem_z   = _mm_add_pd(mem_z, mem_sz);

    tx  = _mm_shuffle_pd(mem_xy, mem_xy, _MM_SHUFFLE2(0, 0));
    ty  = _mm_shuffle_pd(mem_xy, mem_xy, _MM_SHUFFLE2(1, 1));
    tz  = _mm_shuffle_pd(mem_z, mem_z, _MM_SHUFFLE2(0, 0));

    *x1 = gmx_mm256_set_m128d(tx, tx);
    *y1 = gmx_mm256_set_m128d(ty, ty);
    *z1 = gmx_mm256_set_m128d(tz, tz);
}


static gmx_inline void
gmx_mm256_load_shift_and_3rvec_broadcast_pd(const double * gmx_restrict xyz_shift,
                                            const double * gmx_restrict xyz,
                                            __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                            __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                            __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3)
{
    __m128d t1, t2, t3, t4, t5, sxy, sz, szx, syz, tx, ty, tz;

    t1  = _mm_loadu_pd(xyz);
    t2  = _mm_loadu_pd(xyz+2);
    t3  = _mm_loadu_pd(xyz+4);
    t4  = _mm_loadu_pd(xyz+6);
    t5  = _mm_load_sd(xyz+8);

    sxy = _mm_loadu_pd(xyz_shift);
    sz  = _mm_load_sd(xyz_shift+2);
    szx = _mm_shuffle_pd(sz, sxy, _MM_SHUFFLE2(0, 0));
    syz = _mm_shuffle_pd(sxy, sz, _MM_SHUFFLE2(0, 1));

    t1  = _mm_add_pd(t1, sxy);
    t2  = _mm_add_pd(t2, szx);
    t3  = _mm_add_pd(t3, syz);
    t4  = _mm_add_pd(t4, sxy);
    t5  = _mm_add_sd(t5, sz);

    tx   = _mm_shuffle_pd(t1, t1, _MM_SHUFFLE2(0, 0));
    ty   = _mm_shuffle_pd(t1, t1, _MM_SHUFFLE2(1, 1));
    tz   = _mm_shuffle_pd(t2, t2, _MM_SHUFFLE2(0, 0));
    *x1  = gmx_mm256_set_m128d(tx, tx);
    *y1  = gmx_mm256_set_m128d(ty, ty);
    *z1  = gmx_mm256_set_m128d(tz, tz);
    tx   = _mm_shuffle_pd(t2, t2, _MM_SHUFFLE2(1, 1));
    ty   = _mm_shuffle_pd(t3, t3, _MM_SHUFFLE2(0, 0));
    tz   = _mm_shuffle_pd(t3, t3, _MM_SHUFFLE2(1, 1));
    *x2  = gmx_mm256_set_m128d(tx, tx);
    *y2  = gmx_mm256_set_m128d(ty, ty);
    *z2  = gmx_mm256_set_m128d(tz, tz);
    tx   = _mm_shuffle_pd(t4, t4, _MM_SHUFFLE2(0, 0));
    ty   = _mm_shuffle_pd(t4, t4, _MM_SHUFFLE2(1, 1));
    tz   = _mm_shuffle_pd(t5, t5, _MM_SHUFFLE2(0, 0));
    *x3  = gmx_mm256_set_m128d(tx, tx);
    *y3  = gmx_mm256_set_m128d(ty, ty);
    *z3  = gmx_mm256_set_m128d(tz, tz);
}


static gmx_inline void
gmx_mm256_load_shift_and_4rvec_broadcast_pd(const double * gmx_restrict xyz_shift,
                                            const double * gmx_restrict xyz,
                                            __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                            __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                            __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3,
                                            __m256d * gmx_restrict x4, __m256d * gmx_restrict y4, __m256d * gmx_restrict z4)
{
    __m128d t1, t2, t3, t4, t5, t6, sxy, sz, szx, syz, tx, ty, tz;

    t1  = _mm_loadu_pd(xyz);
    t2  = _mm_loadu_pd(xyz+2);
    t3  = _mm_loadu_pd(xyz+4);
    t4  = _mm_loadu_pd(xyz+6);
    t5  = _mm_loadu_pd(xyz+8);
    t6  = _mm_loadu_pd(xyz+10);

    sxy = _mm_loadu_pd(xyz_shift);
    sz  = _mm_load_sd(xyz_shift+2);
    szx = _mm_shuffle_pd(sz, sxy, _MM_SHUFFLE2(0, 0));
    syz = _mm_shuffle_pd(sxy, sz, _MM_SHUFFLE2(0, 1));

    t1  = _mm_add_pd(t1, sxy);
    t2  = _mm_add_pd(t2, szx);
    t3  = _mm_add_pd(t3, syz);
    t4  = _mm_add_pd(t4, sxy);
    t5  = _mm_add_pd(t5, szx);
    t6  = _mm_add_pd(t6, syz);

    tx   = _mm_shuffle_pd(t1, t1, _MM_SHUFFLE2(0, 0));
    ty   = _mm_shuffle_pd(t1, t1, _MM_SHUFFLE2(1, 1));
    tz   = _mm_shuffle_pd(t2, t2, _MM_SHUFFLE2(0, 0));
    *x1  = gmx_mm256_set_m128d(tx, tx);
    *y1  = gmx_mm256_set_m128d(ty, ty);
    *z1  = gmx_mm256_set_m128d(tz, tz);
    tx   = _mm_shuffle_pd(t2, t2, _MM_SHUFFLE2(1, 1));
    ty   = _mm_shuffle_pd(t3, t3, _MM_SHUFFLE2(0, 0));
    tz   = _mm_shuffle_pd(t3, t3, _MM_SHUFFLE2(1, 1));
    *x2  = gmx_mm256_set_m128d(tx, tx);
    *y2  = gmx_mm256_set_m128d(ty, ty);
    *z2  = gmx_mm256_set_m128d(tz, tz);
    tx   = _mm_shuffle_pd(t4, t4, _MM_SHUFFLE2(0, 0));
    ty   = _mm_shuffle_pd(t4, t4, _MM_SHUFFLE2(1, 1));
    tz   = _mm_shuffle_pd(t5, t5, _MM_SHUFFLE2(0, 0));
    *x3  = gmx_mm256_set_m128d(tx, tx);
    *y3  = gmx_mm256_set_m128d(ty, ty);
    *z3  = gmx_mm256_set_m128d(tz, tz);
    tx   = _mm_shuffle_pd(t5, t5, _MM_SHUFFLE2(1, 1));
    ty   = _mm_shuffle_pd(t6, t6, _MM_SHUFFLE2(0, 0));
    tz   = _mm_shuffle_pd(t6, t6, _MM_SHUFFLE2(1, 1));
    *x4  = gmx_mm256_set_m128d(tx, tx);
    *y4  = gmx_mm256_set_m128d(ty, ty);
    *z4  = gmx_mm256_set_m128d(tz, tz);
}


static void
gmx_mm256_load_1rvec_1ptr_swizzle_pd(const double * gmx_restrict p1,
                                     __m256d * gmx_restrict x, __m256d * gmx_restrict y, __m256d * gmx_restrict z)
{
    __m256d t1;

    t1            = _mm256_loadu_pd(p1);
    *x            = t1;
    *y            = _mm256_permute_pd(t1, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *z            = _mm256_castpd128_pd256(_mm256_extractf128_pd(t1, 0x1));
}


static void
gmx_mm256_load_3rvec_1ptr_swizzle_pd(const double * gmx_restrict p1,
                                     __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                     __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                     __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3)
{
    __m256d t1, t2, t3, t4;

    t1            = _mm256_loadu_pd(p1);
    t3            = _mm256_loadu_pd(p1+4);
    *x1           = t1;
    *y2           = t3;
    t2            = gmx_mm256_unpack128hi_pd(t1, t1);
    t4            = gmx_mm256_unpack128hi_pd(t3, t3);
    *z1           = t2;
    *x3           = t4;
    *y1           = _mm256_permute_pd(t1, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *z2           = _mm256_permute_pd(t3, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *x2           = _mm256_permute_pd(t2, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *y3           = _mm256_permute_pd(t4, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *z3           = _mm256_castpd128_pd256(_mm_load_sd(p1+8));
}

static void
gmx_mm256_load_4rvec_1ptr_swizzle_pd(const double * gmx_restrict p1,
                                     __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                     __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                     __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3,
                                     __m256d * gmx_restrict x4, __m256d * gmx_restrict y4, __m256d * gmx_restrict z4)
{
    __m256d t1, t2, t3, t4, t5, t6;

    t1            = _mm256_loadu_pd(p1);
    t2            = _mm256_loadu_pd(p1+4);
    t3            = _mm256_loadu_pd(p1+8);

    t4            = _mm256_castpd128_pd256(_mm256_extractf128_pd(t1, 0x1));
    t5            = _mm256_castpd128_pd256(_mm256_extractf128_pd(t2, 0x1));
    t6            = _mm256_castpd128_pd256(_mm256_extractf128_pd(t3, 0x1));

    *x1           = t1;
    *y2           = t2;
    *z3           = t3;
    *z1           = t4;
    *x3           = t5;
    *y4           = t6;

    *y1           = _mm256_permute_pd(t1, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *z2           = _mm256_permute_pd(t2, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *x4           = _mm256_permute_pd(t3, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *x2           = _mm256_permute_pd(t4, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *y3           = _mm256_permute_pd(t5, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    *z4           = _mm256_permute_pd(t6, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
}


static void
gmx_mm256_load_1rvec_4ptr_swizzle_pd(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                     const double * gmx_restrict ptrC, const double * gmx_restrict ptrD,
                                     __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1)
{
    __m256d t1, t2, t3, t4, t5, t6;

    t1           = _mm256_loadu_pd(ptrA);        /*   -  z1a | y1a x1a */
    t2           = _mm256_loadu_pd(ptrB);        /*   -  z1b | y1b x1b */
    t3           = _mm256_loadu_pd(ptrC);        /*   -  z1c | y1c x1c */
    t4           = _mm256_loadu_pd(ptrD);        /*   -  z1d | y1d x1d */

    t5           = _mm256_unpacklo_pd(t1, t2);   /*  z1b z1a | x1b x1a */
    t6           = _mm256_unpackhi_pd(t1, t2);   /*   -   -  | y1b y1a */
    t1           = _mm256_unpacklo_pd(t3, t4);   /*  z1c z1c | x1d x1c */
    t2           = _mm256_unpackhi_pd(t3, t4);   /*   -   -  | y1d y1c */

    *x1          = gmx_mm256_unpack128lo_pd(t5, t1);
    *y1          = gmx_mm256_unpack128lo_pd(t6, t2);
    *z1          = gmx_mm256_unpack128hi_pd(t5, t1);
}



static void
gmx_mm256_load_3rvec_4ptr_swizzle_pd(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                     const double * gmx_restrict ptrC, const double * gmx_restrict ptrD,
                                     __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                     __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                     __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3)
{
    __m256d t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;

    t1           = _mm256_loadu_pd(ptrA);                       /*  x2a z1a | y1a x1a */
    t2           = _mm256_loadu_pd(ptrB);                       /*  x2b z1b | y1b x1b */
    t3           = _mm256_loadu_pd(ptrC);                       /*  x2c z1c | y1c x1c */
    t4           = _mm256_loadu_pd(ptrD);                       /*  x2d z1d | y1d x1d */
    t5           = _mm256_loadu_pd(ptrA+4);                     /*  y3a x3a | z2a y2a */
    t6           = _mm256_loadu_pd(ptrB+4);                     /*  y3b x3b | z2b y2b */
    t7           = _mm256_loadu_pd(ptrC+4);                     /*  y3c x3c | z2c y2c */
    t8           = _mm256_loadu_pd(ptrD+4);                     /*  y3d x3d | z2d y2d */
    t9           = _mm256_castpd128_pd256(_mm_load_sd(ptrA+8)); /*   -   -  |  -  z3a */
    t10          = _mm256_castpd128_pd256(_mm_load_sd(ptrB+8)); /*   -   -  |  -  z3b */
    t11          = _mm256_castpd128_pd256(_mm_load_sd(ptrC+8)); /*   -   -  |  -  z3c */
    t12          = _mm256_castpd128_pd256(_mm_load_sd(ptrD+8)); /*   -   -  |  -  z3d */

    t13          = _mm256_unpacklo_pd(t1, t2);                  /*  z1b z1a | x1b x1a */
    t14          = _mm256_unpackhi_pd(t1, t2);                  /*  x2b x2a | y1b y1a */
    t1           = _mm256_unpacklo_pd(t3, t4);                  /*  z1d z1c | x1d x1c */
    t2           = _mm256_unpackhi_pd(t3, t4);                  /*  x2d x2c | y1d y1c */

    t3           = _mm256_unpacklo_pd(t5, t6);                  /*  x3b x3a | y2b y2a */
    t4           = _mm256_unpackhi_pd(t5, t6);                  /*  y3b y3a | z2b z2a */
    t5           = _mm256_unpacklo_pd(t7, t8);                  /*  x3d x3c | y2d y2c */
    t6           = _mm256_unpackhi_pd(t7, t8);                  /*  y3d y3c | z2d z2c */

    t9           = _mm256_unpacklo_pd(t9, t10);                 /*   -   -  | z3b z3a */
    t11          = _mm256_unpacklo_pd(t11, t12);                /*   -   -  | z3d z3c */

    *x1          = gmx_mm256_unpack128lo_pd(t13, t1);
    *y1          = gmx_mm256_unpack128lo_pd(t14, t2);
    *z1          = gmx_mm256_unpack128hi_pd(t13, t1);
    *x2          = gmx_mm256_unpack128hi_pd(t14, t2);
    *y2          = gmx_mm256_unpack128lo_pd(t3, t5);
    *z2          = gmx_mm256_unpack128lo_pd(t4, t6);
    *x3          = gmx_mm256_unpack128hi_pd(t3, t5);
    *y3          = gmx_mm256_unpack128hi_pd(t4, t6);
    *z3          = gmx_mm256_unpack128lo_pd(t9, t11);
}



static void
gmx_mm256_load_4rvec_4ptr_swizzle_pd(const double * gmx_restrict ptrA, const double * gmx_restrict ptrB,
                                     const double * gmx_restrict ptrC, const double * gmx_restrict ptrD,
                                     __m256d * gmx_restrict x1, __m256d * gmx_restrict y1, __m256d * gmx_restrict z1,
                                     __m256d * gmx_restrict x2, __m256d * gmx_restrict y2, __m256d * gmx_restrict z2,
                                     __m256d * gmx_restrict x3, __m256d * gmx_restrict y3, __m256d * gmx_restrict z3,
                                     __m256d * gmx_restrict x4, __m256d * gmx_restrict y4, __m256d * gmx_restrict z4)
{
    __m256d t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;

    t1           = _mm256_loadu_pd(ptrA);        /*  x2a z1a | y1a x1a */
    t2           = _mm256_loadu_pd(ptrB);        /*  x2b z1b | y1b x1b */
    t3           = _mm256_loadu_pd(ptrC);        /*  x2c z1c | y1c x1c */
    t4           = _mm256_loadu_pd(ptrD);        /*  x2d z1d | y1d x1d */
    t5           = _mm256_loadu_pd(ptrA+4);      /*  y3a x3a | z2a y2a */
    t6           = _mm256_loadu_pd(ptrB+4);      /*  y3b x3b | z2b y2b */
    t7           = _mm256_loadu_pd(ptrC+4);      /*  y3c x3c | z2c y2c */
    t8           = _mm256_loadu_pd(ptrD+4);      /*  y3d x3d | z2d y2d */
    t9           = _mm256_loadu_pd(ptrA+8);      /*  z4a y4a | x4a z3a */
    t10          = _mm256_loadu_pd(ptrB+8);      /*  z4b y4b | x4b z3b */
    t11          = _mm256_loadu_pd(ptrC+8);      /*  z4c y4c | x4c z3c */
    t12          = _mm256_loadu_pd(ptrD+8);      /*  z4d y4d | x4d z3d */

    t13          = _mm256_unpacklo_pd(t1, t2);   /*  z1b z1a | x1b x1a */
    t14          = _mm256_unpackhi_pd(t1, t2);   /*  x2b x2a | y1b y1a */
    t1           = _mm256_unpacklo_pd(t3, t4);   /*  z1d z1c | x1d x1c */
    t2           = _mm256_unpackhi_pd(t3, t4);   /*  x2d x2c | y1d y1c */

    t3           = _mm256_unpacklo_pd(t5, t6);   /*  x3b x3a | y2b y2a */
    t4           = _mm256_unpackhi_pd(t5, t6);   /*  y3b y3a | z2b z2a */
    t5           = _mm256_unpacklo_pd(t7, t8);   /*  x3d x3c | y2d y2c */
    t6           = _mm256_unpackhi_pd(t7, t8);   /*  y3d y3c | z2d z2c */

    t7           = _mm256_unpacklo_pd(t9, t10);  /*  y4b y4a | z3b z3a */
    t8           = _mm256_unpackhi_pd(t9, t10);  /*  z4b z4a | x4b x4a */
    t9           = _mm256_unpacklo_pd(t11, t12); /*  y4d y4c | z3d z3c */
    t10          = _mm256_unpackhi_pd(t11, t12); /*  z4d z4c | x4d x4c */

    *x1          = gmx_mm256_unpack128lo_pd(t13, t1);
    *y1          = gmx_mm256_unpack128lo_pd(t14, t2);
    *z1          = gmx_mm256_unpack128hi_pd(t13, t1);
    *x2          = gmx_mm256_unpack128hi_pd(t14, t2);
    *y2          = gmx_mm256_unpack128lo_pd(t3, t5);
    *z2          = gmx_mm256_unpack128lo_pd(t4, t6);
    *x3          = gmx_mm256_unpack128hi_pd(t3, t5);
    *y3          = gmx_mm256_unpack128hi_pd(t4, t6);
    *z3          = gmx_mm256_unpack128lo_pd(t7, t9);
    *x4          = gmx_mm256_unpack128lo_pd(t8, t10);
    *y4          = gmx_mm256_unpack128hi_pd(t7, t9);
    *z4          = gmx_mm256_unpack128hi_pd(t8, t10);
}



static void
gmx_mm256_decrement_1rvec_4ptr_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                          double * gmx_restrict ptrC, double * gmx_restrict ptrD,
                                          __m256d x1, __m256d y1, __m256d z1)
{
    __m256d t1, t2, tA, tB, tC, tD;
    __m256i mask;

    t1          = _mm256_unpacklo_pd(x1, y1);       /*  y1c x1c | y1a x1a */
    t2          = _mm256_unpackhi_pd(x1, y1);       /*  y1d x1d | y1b x1b */
    x1          = gmx_mm256_unpack128lo_pd(t1, z1); /*  -  z1a | y1a x1a */
    y1          = gmx_mm256_unpack128hi_pd(t1, z1); /*  -  z1c | y1c x1c */
    z1          = _mm256_permute_pd(z1, _GMX_MM_PERMUTE256D(0, 1, 0, 1));
    t1          = gmx_mm256_unpack128lo_pd(t2, z1); /*  -  z1b | y1b x1b */
    z1          = gmx_mm256_unpack128hi_pd(t2, z1); /*  -  z1d | y1d x1d */

    /* Construct a mask without executing any data loads */
    mask        = _mm256_castpd_si256(_mm256_blend_pd(_mm256_setzero_pd(),
                                                      _mm256_cmp_pd(_mm256_setzero_pd(), _mm256_setzero_pd(), _CMP_EQ_OQ), 0x7));

    tA          = _mm256_loadu_pd(ptrA);
    tB          = _mm256_loadu_pd(ptrB);
    tC          = _mm256_loadu_pd(ptrC);
    tD          = _mm256_loadu_pd(ptrD);

    tA          = _mm256_sub_pd(tA, x1);
    tB          = _mm256_sub_pd(tB, t1);
    tC          = _mm256_sub_pd(tC, y1);
    tD          = _mm256_sub_pd(tD, z1);

    _mm256_maskstore_pd(ptrA, mask, tA);
    _mm256_maskstore_pd(ptrB, mask, tB);
    _mm256_maskstore_pd(ptrC, mask, tC);
    _mm256_maskstore_pd(ptrD, mask, tD);
}



#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm256_decrement_3rvec_4ptr_swizzle_pd(ptrA, ptrB, ptrC, ptrD, \
                                                  _x1, _y1, _z1, _x2, _y2, _z2, _x3, _y3, _z3) \
    { \
        __m256d _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8, _t9, _t10; \
        __m128d _tA, _tB, _tC, _tD, _tE; \
        _t1          = _mm256_loadu_pd(ptrA); \
        _t2          = _mm256_loadu_pd(ptrB); \
        _t3          = _mm256_loadu_pd(ptrC); \
        _t4          = _mm256_loadu_pd(ptrD); \
        _t5          = _mm256_loadu_pd(ptrA+4); \
        _t6          = _mm256_loadu_pd(ptrB+4); \
        _t7          = _mm256_loadu_pd(ptrC+4); \
        _t8          = _mm256_loadu_pd(ptrD+4); \
        _tA          = _mm_load_sd(ptrA+8); \
        _tB          = _mm_load_sd(ptrB+8); \
        _tC          = _mm_load_sd(ptrC+8); \
        _tD          = _mm_load_sd(ptrD+8); \
        _t9          = _mm256_unpacklo_pd(_x1, _y1); \
        _x1          = _mm256_unpackhi_pd(_x1, _y1); \
        _y1          = _mm256_unpacklo_pd(_z1, _x2); \
        _z1          = _mm256_unpackhi_pd(_z1, _x2); \
        _x2          = _mm256_unpacklo_pd(_y2, _z2); \
        _y2          = _mm256_unpackhi_pd(_y2, _z2); \
        _z2          = _mm256_unpacklo_pd(_x3, _y3); \
        _x3          = _mm256_unpackhi_pd(_x3, _y3); \
        _t10         = gmx_mm256_unpack128lo_pd(_t9, _y1); \
        _y3          = gmx_mm256_unpack128hi_pd(_t9, _y1); \
        _t9          = gmx_mm256_unpack128lo_pd(_x1, _z1); \
        _y1          = gmx_mm256_unpack128hi_pd(_x1, _z1); \
        _x1          = gmx_mm256_unpack128lo_pd(_x2, _z2); \
        _z1          = gmx_mm256_unpack128hi_pd(_x2, _z2); \
        _x2          = gmx_mm256_unpack128lo_pd(_y2, _x3); \
        _z2          = gmx_mm256_unpack128hi_pd(_y2, _x3); \
        _t1          = _mm256_sub_pd(_t1, _t10); \
        _t2          = _mm256_sub_pd(_t2, _t9); \
        _t3          = _mm256_sub_pd(_t3, _y3); \
        _t4          = _mm256_sub_pd(_t4, _y1); \
        _t5          = _mm256_sub_pd(_t5, _x1); \
        _t6          = _mm256_sub_pd(_t6, _x2); \
        _t7          = _mm256_sub_pd(_t7, _z1); \
        _t8          = _mm256_sub_pd(_t8, _z2); \
        _tA          = _mm_sub_sd(_tA, _mm256_castpd256_pd128(_z3)); \
        _tB          = _mm_sub_sd(_tB, _mm_permute_pd(_mm256_castpd256_pd128(_z3), _GMX_MM_PERMUTE128D(1, 1))); \
        _tE          = _mm256_extractf128_pd(_z3, 0x1); \
        _tC          = _mm_sub_sd(_tC, _tE); \
        _tD          = _mm_sub_sd(_tD, _mm_permute_pd(_tE, _GMX_MM_PERMUTE128D(1, 1))); \
        _mm256_storeu_pd(ptrA, _t1); \
        _mm256_storeu_pd(ptrB, _t2); \
        _mm256_storeu_pd(ptrC, _t3); \
        _mm256_storeu_pd(ptrD, _t4); \
        _mm256_storeu_pd(ptrA+4, _t5); \
        _mm256_storeu_pd(ptrB+4, _t6); \
        _mm256_storeu_pd(ptrC+4, _t7); \
        _mm256_storeu_pd(ptrD+4, _t8); \
        _mm_store_sd(ptrA+8, _tA); \
        _mm_store_sd(ptrB+8, _tB); \
        _mm_store_sd(ptrC+8, _tC); \
        _mm_store_sd(ptrD+8, _tD); \
    }
#else
/* Real function for sane compilers */
static void
gmx_mm256_decrement_3rvec_4ptr_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                          double * gmx_restrict ptrC, double * gmx_restrict ptrD,
                                          __m256d x1, __m256d y1, __m256d z1,
                                          __m256d x2, __m256d y2, __m256d z2,
                                          __m256d x3, __m256d y3, __m256d z3)
{
    __m256d t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    __m128d tA, tB, tC, tD, tE;

    t1          = _mm256_loadu_pd(ptrA);
    t2          = _mm256_loadu_pd(ptrB);
    t3          = _mm256_loadu_pd(ptrC);
    t4          = _mm256_loadu_pd(ptrD);
    t5          = _mm256_loadu_pd(ptrA+4);
    t6          = _mm256_loadu_pd(ptrB+4);
    t7          = _mm256_loadu_pd(ptrC+4);
    t8          = _mm256_loadu_pd(ptrD+4);
    tA          = _mm_load_sd(ptrA+8);
    tB          = _mm_load_sd(ptrB+8);
    tC          = _mm_load_sd(ptrC+8);
    tD          = _mm_load_sd(ptrD+8);

    t9          = _mm256_unpacklo_pd(x1, y1);       /* y1c x1c | y1a x1a */
    x1          = _mm256_unpackhi_pd(x1, y1);       /* y1d x1d | y1b x1b */

    y1          = _mm256_unpacklo_pd(z1, x2);       /* x2c z1c | x2a z1a */
    z1          = _mm256_unpackhi_pd(z1, x2);       /* x2d z1d | x2b z1b */

    x2          = _mm256_unpacklo_pd(y2, z2);       /* z2c y2c | z2a y2a */
    y2          = _mm256_unpackhi_pd(y2, z2);       /* z2d y2d | z2b y2b */

    z2          = _mm256_unpacklo_pd(x3, y3);       /* y3c x3c | y3a x3a */
    x3          = _mm256_unpackhi_pd(x3, y3);       /* y3d x3d | y3b x3b */

    t10         = gmx_mm256_unpack128lo_pd(t9, y1); /* x2a z1a | y1a x1a */
    y3          = gmx_mm256_unpack128hi_pd(t9, y1); /* x2c z1c | y1c x1c */

    t9          = gmx_mm256_unpack128lo_pd(x1, z1); /* x2b z1b | y1b x1b */
    y1          = gmx_mm256_unpack128hi_pd(x1, z1); /* x2d z1d | y1d x1d */

    x1          = gmx_mm256_unpack128lo_pd(x2, z2); /* y3a x3a | z2a y2a */
    z1          = gmx_mm256_unpack128hi_pd(x2, z2); /* y3c x3c | z2c y2c */

    x2          = gmx_mm256_unpack128lo_pd(y2, x3); /* y3b x3b | z2b y2b */
    z2          = gmx_mm256_unpack128hi_pd(y2, x3); /* y3d x3d | z2d y2d */

    t1          = _mm256_sub_pd(t1, t10);
    t2          = _mm256_sub_pd(t2, t9);
    t3          = _mm256_sub_pd(t3, y3);
    t4          = _mm256_sub_pd(t4, y1);
    t5          = _mm256_sub_pd(t5, x1);
    t6          = _mm256_sub_pd(t6, x2);
    t7          = _mm256_sub_pd(t7, z1);
    t8          = _mm256_sub_pd(t8, z2);

    tA          = _mm_sub_sd(tA, _mm256_castpd256_pd128(z3));
    tB          = _mm_sub_sd(tB, _mm_permute_pd(_mm256_castpd256_pd128(z3), _GMX_MM_PERMUTE128D(1, 1)));
    tE          = _mm256_extractf128_pd(z3, 0x1);
    tC          = _mm_sub_sd(tC, tE);
    tD          = _mm_sub_sd(tD, _mm_permute_pd(tE, _GMX_MM_PERMUTE128D(1, 1)));

    /* Here we store a full 256-bit value and a separate 64-bit one; no overlap can happen */
    _mm256_storeu_pd(ptrA, t1);
    _mm256_storeu_pd(ptrB, t2);
    _mm256_storeu_pd(ptrC, t3);
    _mm256_storeu_pd(ptrD, t4);
    _mm256_storeu_pd(ptrA+4, t5);
    _mm256_storeu_pd(ptrB+4, t6);
    _mm256_storeu_pd(ptrC+4, t7);
    _mm256_storeu_pd(ptrD+4, t8);
    _mm_store_sd(ptrA+8, tA);
    _mm_store_sd(ptrB+8, tB);
    _mm_store_sd(ptrC+8, tC);
    _mm_store_sd(ptrD+8, tD);
}
#endif

#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm256_decrement_4rvec_4ptr_swizzle_pd(ptrA, ptrB, ptrC, ptrD, \
                                                  _x1, _y1, _z1, _x2, _y2, _z2, _x3, _y3, _z3, _x4, _y4, _z4) \
    { \
        __m256d _t1, _t2, _t3, _t4, _t5, _t6, _t7, _t8, _t9, _t10, _t11, _t12, _t13, _t14; \
        __m128d _tA, _tB, _tC, _tD, _tE; \
        _t1          = _mm256_loadu_pd(ptrA); \
        _t2          = _mm256_loadu_pd(ptrB); \
        _t3          = _mm256_loadu_pd(ptrC); \
        _t4          = _mm256_loadu_pd(ptrD); \
        _t5          = _mm256_loadu_pd(ptrA+4); \
        _t6          = _mm256_loadu_pd(ptrB+4); \
        _t7          = _mm256_loadu_pd(ptrC+4); \
        _t8          = _mm256_loadu_pd(ptrD+4); \
        _t9          = _mm256_loadu_pd(ptrA+8); \
        _t10         = _mm256_loadu_pd(ptrB+8); \
        _t11         = _mm256_loadu_pd(ptrC+8); \
        _t12         = _mm256_loadu_pd(ptrD+8); \
        _t13         = _mm256_unpacklo_pd(_x1, _y1); \
        _x1          = _mm256_unpackhi_pd(_x1, _y1); \
        _y1          = _mm256_unpacklo_pd(_z1, _x2); \
        _z1          = _mm256_unpackhi_pd(_z1, _x2); \
        _x2          = _mm256_unpacklo_pd(_y2, _z2); \
        _y2          = _mm256_unpackhi_pd(_y2, _z2); \
        _z2          = _mm256_unpacklo_pd(_x3, _y3); \
        _x3          = _mm256_unpackhi_pd(_x3, _y3); \
        _y3          = _mm256_unpacklo_pd(_z3, _x4); \
        _z3          = _mm256_unpackhi_pd(_z3, _x4); \
        _x4          = _mm256_unpacklo_pd(_y4, _z4); \
        _y4          = _mm256_unpackhi_pd(_y4, _z4); \
        _z4          = gmx_mm256_unpack128lo_pd(_t13, _y1); \
        _t13         = gmx_mm256_unpack128hi_pd(_t13, _y1); \
        _y1          = gmx_mm256_unpack128lo_pd(_x1, _z1); \
        _x1          = gmx_mm256_unpack128hi_pd(_x1, _z1); \
        _z1          = gmx_mm256_unpack128lo_pd(_x2, _z2); \
        _x2          = gmx_mm256_unpack128hi_pd(_x2, _z2); \
        _z2          = gmx_mm256_unpack128lo_pd(_y2, _x3); \
        _y2          = gmx_mm256_unpack128hi_pd(_y2, _x3); \
        _x3          = gmx_mm256_unpack128lo_pd(_y3, _x4); \
        _y3          = gmx_mm256_unpack128hi_pd(_y3, _x4); \
        _x4          = gmx_mm256_unpack128lo_pd(_z3, _y4); \
        _z3          = gmx_mm256_unpack128hi_pd(_z3, _y4); \
        _t1          = _mm256_sub_pd(_t1, _z4); \
        _t2          = _mm256_sub_pd(_t2, _y1); \
        _t3          = _mm256_sub_pd(_t3, _t13); \
        _t4          = _mm256_sub_pd(_t4, _x1); \
        _t5          = _mm256_sub_pd(_t5, _z1); \
        _t6          = _mm256_sub_pd(_t6, _z2); \
        _t7          = _mm256_sub_pd(_t7, _x2); \
        _t8          = _mm256_sub_pd(_t8, _y2); \
        _t9          = _mm256_sub_pd(_t9, _x3); \
        _t10         = _mm256_sub_pd(_t10, _x4); \
        _t11         = _mm256_sub_pd(_t11, _y3); \
        _t12         = _mm256_sub_pd(_t12, _z3); \
        _mm256_storeu_pd(ptrA, _t1); \
        _mm256_storeu_pd(ptrB, _t2); \
        _mm256_storeu_pd(ptrC, _t3); \
        _mm256_storeu_pd(ptrD, _t4); \
        _mm256_storeu_pd(ptrA+4, _t5); \
        _mm256_storeu_pd(ptrB+4, _t6); \
        _mm256_storeu_pd(ptrC+4, _t7); \
        _mm256_storeu_pd(ptrD+4, _t8); \
        _mm256_storeu_pd(ptrA+8, _t9); \
        _mm256_storeu_pd(ptrB+8, _t10); \
        _mm256_storeu_pd(ptrC+8, _t11); \
        _mm256_storeu_pd(ptrD+8, _t12); \
    }
#else
/* Real function for sane compilers */
static void
gmx_mm256_decrement_4rvec_4ptr_swizzle_pd(double * gmx_restrict ptrA, double * gmx_restrict ptrB,
                                          double * gmx_restrict ptrC, double * gmx_restrict ptrD,
                                          __m256d x1, __m256d y1, __m256d z1,
                                          __m256d x2, __m256d y2, __m256d z2,
                                          __m256d x3, __m256d y3, __m256d z3,
                                          __m256d x4, __m256d y4, __m256d z4)
{
    __m256d t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
    __m128d tA, tB, tC, tD, tE;

    t1          = _mm256_loadu_pd(ptrA);
    t2          = _mm256_loadu_pd(ptrB);
    t3          = _mm256_loadu_pd(ptrC);
    t4          = _mm256_loadu_pd(ptrD);
    t5          = _mm256_loadu_pd(ptrA+4);
    t6          = _mm256_loadu_pd(ptrB+4);
    t7          = _mm256_loadu_pd(ptrC+4);
    t8          = _mm256_loadu_pd(ptrD+4);
    t9          = _mm256_loadu_pd(ptrA+8);
    t10         = _mm256_loadu_pd(ptrB+8);
    t11         = _mm256_loadu_pd(ptrC+8);
    t12         = _mm256_loadu_pd(ptrD+8);

    t13         = _mm256_unpacklo_pd(x1, y1);        /* y1c x1c | y1a x1a */
    x1          = _mm256_unpackhi_pd(x1, y1);        /* y1d x1d | y1b x1b */
    y1          = _mm256_unpacklo_pd(z1, x2);        /* x2c z1c | x2a z1a */
    z1          = _mm256_unpackhi_pd(z1, x2);        /* x2d z1d | x2b z1b */
    x2          = _mm256_unpacklo_pd(y2, z2);        /* z2c y2c | z2a y2a */
    y2          = _mm256_unpackhi_pd(y2, z2);        /* z2d y2d | z2b y2b */
    z2          = _mm256_unpacklo_pd(x3, y3);        /* y3c x3c | y3a x3a */
    x3          = _mm256_unpackhi_pd(x3, y3);        /* y3d x3d | y3b x3b */
    y3          = _mm256_unpacklo_pd(z3, x4);        /* x4c z3c | x4a z3a */
    z3          = _mm256_unpackhi_pd(z3, x4);        /* x4d z3d | x4b z3b */
    x4          = _mm256_unpacklo_pd(y4, z4);        /* z4c y4c | z4a y4a */
    y4          = _mm256_unpackhi_pd(y4, z4);        /* z4d y4d | z4b y4b */

    z4          = gmx_mm256_unpack128lo_pd(t13, y1); /* x2a z1a | y1a x1a */
    t13         = gmx_mm256_unpack128hi_pd(t13, y1); /* x2c z1c | y1c x1c */
    y1          = gmx_mm256_unpack128lo_pd(x1, z1);  /* x2b z1b | y1b x1b */
    x1          = gmx_mm256_unpack128hi_pd(x1, z1);  /* x2d z1d | y1d x1d */
    z1          = gmx_mm256_unpack128lo_pd(x2, z2);  /* y3a x3a | z2a y2a */
    x2          = gmx_mm256_unpack128hi_pd(x2, z2);  /* y3c x3c | z2c y2c */
    z2          = gmx_mm256_unpack128lo_pd(y2, x3);  /* y3b x3b | z2b y2b */
    y2          = gmx_mm256_unpack128hi_pd(y2, x3);  /* y3d x3d | z2d y2d */
    x3          = gmx_mm256_unpack128lo_pd(y3, x4);  /* z4a y4a | x4a z3a */
    y3          = gmx_mm256_unpack128hi_pd(y3, x4);  /* z4c y4c | x4c z3c */
    x4          = gmx_mm256_unpack128lo_pd(z3, y4);  /* z4b y4b | x4b z3b */
    z3          = gmx_mm256_unpack128hi_pd(z3, y4);  /* z4d y4d | x4d z3d */

    t1          = _mm256_sub_pd(t1, z4);
    t2          = _mm256_sub_pd(t2, y1);
    t3          = _mm256_sub_pd(t3, t13);
    t4          = _mm256_sub_pd(t4, x1);
    t5          = _mm256_sub_pd(t5, z1);
    t6          = _mm256_sub_pd(t6, z2);
    t7          = _mm256_sub_pd(t7, x2);
    t8          = _mm256_sub_pd(t8, y2);
    t9          = _mm256_sub_pd(t9, x3);
    t10         = _mm256_sub_pd(t10, x4);
    t11         = _mm256_sub_pd(t11, y3);
    t12         = _mm256_sub_pd(t12, z3);

    /* Here we store a full 256-bit value and a separate 128-bit one; no overlap can happen */
    _mm256_storeu_pd(ptrA, t1);
    _mm256_storeu_pd(ptrB, t2);
    _mm256_storeu_pd(ptrC, t3);
    _mm256_storeu_pd(ptrD, t4);
    _mm256_storeu_pd(ptrA+4, t5);
    _mm256_storeu_pd(ptrB+4, t6);
    _mm256_storeu_pd(ptrC+4, t7);
    _mm256_storeu_pd(ptrD+4, t8);
    _mm256_storeu_pd(ptrA+8, t9);
    _mm256_storeu_pd(ptrB+8, t10);
    _mm256_storeu_pd(ptrC+8, t11);
    _mm256_storeu_pd(ptrD+8, t12);
}
#endif





static gmx_inline void
gmx_mm256_update_iforce_1atom_swizzle_pd(__m256d fix1, __m256d fiy1, __m256d fiz1,
                                         double * gmx_restrict fptr,
                                         double * gmx_restrict fshiftptr)
{
    __m256d t1, t2;
    __m128d tA, tB;
    fix1 = _mm256_hadd_pd(fix1, fiy1);
    fiz1 = _mm256_hadd_pd(fiz1, _mm256_setzero_pd());

    /* Add across the two lanes */
    tA   = _mm_add_pd(_mm256_castpd256_pd128(fix1), _mm256_extractf128_pd(fix1, 0x1));
    tB   = _mm_add_pd(_mm256_castpd256_pd128(fiz1), _mm256_extractf128_pd(fiz1, 0x1));

    fix1 = gmx_mm256_set_m128d(tB, tA); /* 0 fiz fiy fix */

    t1   = _mm256_loadu_pd(fptr);
    t2   = _mm256_loadu_pd(fshiftptr);

    t1   = _mm256_add_pd(t1, fix1);
    t2   = _mm256_add_pd(t2, fix1);

    _mm256_storeu_pd(fptr, t1);
    _mm256_storeu_pd(fshiftptr, t2);
}



#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm256_update_iforce_3atom_swizzle_pd(fix1, fiy1, fiz1, fix2, fiy2, fiz2, fix3, fiy3, fiz3, \
                                                 fptr, fshiftptr) \
    { \
        __m256d _t1, _t2, _t3, _t4; \
        __m128d _tz3, _tA, _tB, _tC, _tD; \
        fix1  = _mm256_hadd_pd(fix1, fiy1); \
        fiz1  = _mm256_hadd_pd(fiz1, fix2); \
        fiy2  = _mm256_hadd_pd(fiy2, fiz2); \
        fix3  = _mm256_hadd_pd(fix3, fiy3); \
        fiz3  = _mm256_hadd_pd(fiz3, _mm256_setzero_pd()); \
        _t1   = gmx_mm256_unpack128lo_pd(fix1, fiz1); \
        _t2   = gmx_mm256_unpack128hi_pd(fix1, fiz1); \
        _t1   = _mm256_add_pd(_t1, _t2); \
        _t3   = gmx_mm256_unpack128lo_pd(fiy2, fix3); \
        _t4   = gmx_mm256_unpack128hi_pd(fiy2, fix3); \
        _t3   = _mm256_add_pd(_t3, _t4); \
        _tz3  = _mm_add_pd(_mm256_castpd256_pd128(fiz3), _mm256_extractf128_pd(fiz3, 0x1)); \
        _t2   = _mm256_loadu_pd(fptr); \
        _t4   = _mm256_loadu_pd(fptr+4); \
        _tA   = _mm_load_sd(fptr+8); \
        _t2   = _mm256_add_pd(_t2, _t1); \
        _t4   = _mm256_add_pd(_t4, _t3); \
        _tA   = _mm_add_sd(_tA, _tz3); \
        _mm256_storeu_pd(fptr, _t2); \
        _mm256_storeu_pd(fptr+4, _t4); \
        _mm_store_sd(fptr+8, _tA); \
        _tB   = _mm256_extractf128_pd(_t1, 0x1); \
        _tC   = _mm256_extractf128_pd(_t3, 0x1); \
        _tz3  = _mm_add_sd(_tz3, _tB); \
        _tD   = _mm_permute_pd(_mm256_castpd256_pd128(_t3), _GMX_MM_PERMUTE128D(1, 1)); \
        _tz3  = _mm_add_sd(_tz3, _tD); \
        _tC   = _mm_add_pd(_tC, _mm256_castpd256_pd128(_t1)); \
        _tD   = _mm_shuffle_pd(_tB, _mm256_castpd256_pd128(_t3), _MM_SHUFFLE2(0, 1)); \
        _tC   = _mm_add_pd(_tC, _tD); \
        _tA   = _mm_loadu_pd(fshiftptr); \
        _tB   = _mm_load_sd(fshiftptr+2); \
        _tA   = _mm_add_pd(_tA, _tC); \
        _tB   = _mm_add_sd(_tB, _tz3); \
        _mm_storeu_pd(fshiftptr, _tA); \
        _mm_store_sd(fshiftptr+2, _tB); \
    }
#else
/* Real function for sane compilers */
static gmx_inline void
gmx_mm256_update_iforce_3atom_swizzle_pd(__m256d fix1, __m256d fiy1, __m256d fiz1,
                                         __m256d fix2, __m256d fiy2, __m256d fiz2,
                                         __m256d fix3, __m256d fiy3, __m256d fiz3,
                                         double * gmx_restrict fptr,
                                         double * gmx_restrict fshiftptr)
{
    __m256d t1, t2, t3, t4;
    __m128d tz3, tA, tB, tC, tD;

    fix1 = _mm256_hadd_pd(fix1, fiy1);                /*  Y1c-d X1c-d | Y1a-b X1a-b */
    fiz1 = _mm256_hadd_pd(fiz1, fix2);                /*  X2c-d Z1c-d | X2a-b Z1a-b */
    fiy2 = _mm256_hadd_pd(fiy2, fiz2);                /*  Z2c-d Y2c-d | Z2a-b Y2a-b */
    fix3 = _mm256_hadd_pd(fix3, fiy3);                /*  Y3c-d X3c-d | Y3a-b X3a-b */
    fiz3 = _mm256_hadd_pd(fiz3, _mm256_setzero_pd()); /*  0     Z3c-d | 0     Z3a-b */

    /* Add across the two lanes by swapping and adding back */
    t1   = gmx_mm256_unpack128lo_pd(fix1, fiz1);                                       /* X2a-b Z1a-b | Y1a-b X1a-b */
    t2   = gmx_mm256_unpack128hi_pd(fix1, fiz1);                                       /* X2c-d Z1c-d | Y1c-d X1c-d */
    t1   = _mm256_add_pd(t1, t2);                                                      /* x2 z1 | y1 x1 */

    t3   = gmx_mm256_unpack128lo_pd(fiy2, fix3);                                       /* Y3a-b X3a-b | Z2a-b Y2a-b */
    t4   = gmx_mm256_unpack128hi_pd(fiy2, fix3);                                       /* Y3c-d X3c-d | Z2c-d Y2c-d */
    t3   = _mm256_add_pd(t3, t4);                                                      /* y3 x3 | z2 y2 */

    tz3  = _mm_add_pd(_mm256_castpd256_pd128(fiz3), _mm256_extractf128_pd(fiz3, 0x1)); /* 0 z3 */

    t2   = _mm256_loadu_pd(fptr);
    t4   = _mm256_loadu_pd(fptr+4);
    tA   = _mm_load_sd(fptr+8);

    t2   = _mm256_add_pd(t2, t1);
    t4   = _mm256_add_pd(t4, t3);
    tA   = _mm_add_sd(tA, tz3);

    _mm256_storeu_pd(fptr, t2);
    _mm256_storeu_pd(fptr+4, t4);
    _mm_store_sd(fptr+8, tA);

    /* Add up shift force */
    /* t1:   x2 z1 | y1 x1 */
    /* t3:   y3 x3 | z2 y2 */
    /* tz3:           0 z3 */

    /* z component */
    tB   = _mm256_extractf128_pd(t1, 0x1);                                     /* x2 z1 */
    tC   = _mm256_extractf128_pd(t3, 0x1);                                     /* y3 x3 */
    tz3  = _mm_add_sd(tz3, tB);                                                /* 0  z1+z3 */
    tD   = _mm_permute_pd(_mm256_castpd256_pd128(t3), _GMX_MM_PERMUTE128D(1, 1));
    tz3  = _mm_add_sd(tz3, tD);                                                /* - z */

    tC   = _mm_add_pd(tC, _mm256_castpd256_pd128(t1));                         /* y1+y3 x1+x3 */

    tD   = _mm_shuffle_pd(tB, _mm256_castpd256_pd128(t3), _MM_SHUFFLE2(0, 1)); /* y2 x2 */
    tC   = _mm_add_pd(tC, tD);                                                 /* y x */

    tA   = _mm_loadu_pd(fshiftptr);
    tB   = _mm_load_sd(fshiftptr+2);
    tA   = _mm_add_pd(tA, tC);
    tB   = _mm_add_sd(tB, tz3);
    _mm_storeu_pd(fshiftptr, tA);
    _mm_store_sd(fshiftptr+2, tB);
}
#endif


#if defined (_MSC_VER) && defined(_M_IX86)
/* Macro work-around since 32-bit MSVC cannot handle >3 xmm/ymm parameters */
#define gmx_mm256_update_iforce_4atom_swizzle_pd(fix1, fiy1, fiz1, fix2, fiy2, fiz2, fix3, fiy3, fiz3, fix4, fiy4, fiz4, \
                                                 fptr, fshiftptr) \
    { \
        __m256d _t1, _t2, _t3, _t4, _t5, _t6; \
        __m128d _tA, _tB, _tC, _tD; \
        fix1  = _mm256_hadd_pd(fix1, fiy1); \
        fiz1  = _mm256_hadd_pd(fiz1, fix2); \
        fiy2  = _mm256_hadd_pd(fiy2, fiz2); \
        fix3  = _mm256_hadd_pd(fix3, fiy3); \
        fiz3  = _mm256_hadd_pd(fiz3, fix4); \
        fiy4  = _mm256_hadd_pd(fiy4, fiz4); \
        _t1   = gmx_mm256_unpack128lo_pd(fix1, fiz1); \
        _t2   = gmx_mm256_unpack128hi_pd(fix1, fiz1); \
        _t1   = _mm256_add_pd(_t1, _t2); \
        _t3   = gmx_mm256_unpack128lo_pd(fiy2, fix3); \
        _t4   = gmx_mm256_unpack128hi_pd(fiy2, fix3); \
        _t3   = _mm256_add_pd(_t3, _t4); \
        _t5   = gmx_mm256_unpack128lo_pd(fiz3, fiy4); \
        _t6   = gmx_mm256_unpack128hi_pd(fiz3, fiy4); \
        _t5   = _mm256_add_pd(_t5, _t6); \
        _t2   = _mm256_loadu_pd(fptr); \
        _t4   = _mm256_loadu_pd(fptr+4); \
        _t6   = _mm256_loadu_pd(fptr+8); \
        _t2   = _mm256_add_pd(_t2, _t1); \
        _t4   = _mm256_add_pd(_t4, _t3); \
        _t6   = _mm256_add_pd(_t6, _t5); \
        _mm256_storeu_pd(fptr, _t2); \
        _mm256_storeu_pd(fptr+4, _t4); \
        _mm256_storeu_pd(fptr+8, _t6); \
        _tA   = _mm256_extractf128_pd(_t1, 0x1); \
        _tB   = _mm256_extractf128_pd(_t3, 0x1); \
        _tC   = _mm256_extractf128_pd(_t5, 0x1); \
        _tB   = _mm_add_pd(_tB, _mm256_castpd256_pd128(_t1)); \
        _tA   = _mm_add_pd(_tA, _mm256_castpd256_pd128(_t5)); \
        _tC   = _mm_add_pd(_tC, _mm256_castpd256_pd128(_t3)); \
        _tD   = _mm_shuffle_pd(_tA, _tC, _MM_SHUFFLE2(0, 1)); \
        _tB   = _mm_add_pd(_tB, _tD); \
        _tC   = _mm_permute_pd(_tC, _GMX_MM_PERMUTE128D(1, 1)); \
        _tC   = _mm_add_sd(_tC, _tA); \
        _tA   = _mm_loadu_pd(fshiftptr); \
        _tD   = _mm_load_sd(fshiftptr+2); \
        _tA   = _mm_add_pd(_tA, _tB); \
        _tD   = _mm_add_sd(_tD, _tC); \
        _mm_storeu_pd(fshiftptr, _tA); \
        _mm_store_sd(fshiftptr+2, _tD); \
    }
#else
/* Real function for sane compilers */
static gmx_inline void
gmx_mm256_update_iforce_4atom_swizzle_pd(__m256d fix1, __m256d fiy1, __m256d fiz1,
                                         __m256d fix2, __m256d fiy2, __m256d fiz2,
                                         __m256d fix3, __m256d fiy3, __m256d fiz3,
                                         __m256d fix4, __m256d fiy4, __m256d fiz4,
                                         double * gmx_restrict fptr,
                                         double * gmx_restrict fshiftptr)
{
    __m256d t1, t2, t3, t4, t5, t6;
    __m128d tA, tB, tC, tD;

    fix1 = _mm256_hadd_pd(fix1, fiy1);                /*  Y1c-d X1c-d | Y1a-b X1a-b */
    fiz1 = _mm256_hadd_pd(fiz1, fix2);                /*  X2c-d Z1c-d | X2a-b Z1a-b */
    fiy2 = _mm256_hadd_pd(fiy2, fiz2);                /*  Z2c-d Y2c-d | Z2a-b Y2a-b */
    fix3 = _mm256_hadd_pd(fix3, fiy3);                /*  Y3c-d X3c-d | Y3a-b X3a-b */
    fiz3 = _mm256_hadd_pd(fiz3, fix4);                /*  X4c-d Z3c-d | X4a-b Z3a-b */
    fiy4 = _mm256_hadd_pd(fiy4, fiz4);                /*  Z4c-d Y4c-d | Z4a-b Y4a-b */

    /* Add across the two lanes by swapping and adding back */
    t1   = gmx_mm256_unpack128lo_pd(fix1, fiz1); /* X2a-b Z1a-b | Y1a-b X1a-b */
    t2   = gmx_mm256_unpack128hi_pd(fix1, fiz1); /* X2c-d Z1c-d | Y1c-d X1c-d */
    t1   = _mm256_add_pd(t1, t2);                /* x2 z1 | y1 x1 */

    t3   = gmx_mm256_unpack128lo_pd(fiy2, fix3); /* Y3a-b X3a-b | Z2a-b Y2a-b */
    t4   = gmx_mm256_unpack128hi_pd(fiy2, fix3); /* Y3c-d X3c-d | Z2c-d Y2c-d */
    t3   = _mm256_add_pd(t3, t4);                /* y3 x3 | z2 y2 */

    t5   = gmx_mm256_unpack128lo_pd(fiz3, fiy4); /* Z4a-b Y4a-b | X4a-b Z3a-b */
    t6   = gmx_mm256_unpack128hi_pd(fiz3, fiy4); /* Z4c-d Y4c-d | X4c-d Z3c-d */
    t5   = _mm256_add_pd(t5, t6);                /* z4 y4 | x4 z3 */

    t2   = _mm256_loadu_pd(fptr);
    t4   = _mm256_loadu_pd(fptr+4);
    t6   = _mm256_loadu_pd(fptr+8);

    t2   = _mm256_add_pd(t2, t1);
    t4   = _mm256_add_pd(t4, t3);
    t6   = _mm256_add_pd(t6, t5);

    _mm256_storeu_pd(fptr, t2);
    _mm256_storeu_pd(fptr+4, t4);
    _mm256_storeu_pd(fptr+8, t6);

    /* Add up shift force  */
    /* t1:   x2. z1. | y1. x1. */
    /* t3:   y3. x3. | z2 y2 */
    /* t5:   z4 y4 | x4. z3. */

    /* z component */
    tA   = _mm256_extractf128_pd(t1, 0x1);                /* x2 z1 */
    tB   = _mm256_extractf128_pd(t3, 0x1);                /* y3 x3 */
    tC   = _mm256_extractf128_pd(t5, 0x1);                /* z4 y4 */

    tB   = _mm_add_pd(tB, _mm256_castpd256_pd128(t1));    /*  y1+y3  x1+x3 */
    tA   = _mm_add_pd(tA, _mm256_castpd256_pd128(t5));    /*  x2+x4  z1+z3 */
    tC   = _mm_add_pd(tC, _mm256_castpd256_pd128(t3));    /*  z4+z2  y4+y2 */

    tD   = _mm_shuffle_pd(tA, tC, _MM_SHUFFLE2(0, 1));    /* y4+y2 x2+x4 */
    tB   = _mm_add_pd(tB, tD);                            /* y x */
    tC   = _mm_permute_pd(tC, _GMX_MM_PERMUTE128D(1, 1)); /*    - z4+z2 */
    tC   = _mm_add_sd(tC, tA);                            /* - z */

    tA   = _mm_loadu_pd(fshiftptr);
    tD   = _mm_load_sd(fshiftptr+2);
    tA   = _mm_add_pd(tA, tB);
    tD   = _mm_add_sd(tD, tC);
    _mm_storeu_pd(fshiftptr, tA);
    _mm_store_sd(fshiftptr+2, tD);
}
#endif



static void
gmx_mm256_update_1pot_pd(__m256d pot1, double * gmx_restrict ptrA)
{
    __m128d t1;

    pot1 = _mm256_hadd_pd(pot1, pot1);

    t1   = _mm_add_pd(_mm256_castpd256_pd128(pot1), _mm256_extractf128_pd(pot1, 0x1));

    _mm_store_sd(ptrA, _mm_add_sd(_mm_load_sd(ptrA), t1));
}

static void
gmx_mm256_update_2pot_pd(__m256d pot1, double * gmx_restrict ptrA,
                         __m256d pot2, double * gmx_restrict ptrB)
{
    __m128d t1, t2;

    pot1 = _mm256_hadd_pd(pot1, pot2);

    t1   = _mm_add_pd(_mm256_castpd256_pd128(pot1), _mm256_extractf128_pd(pot1, 0x1));

    t2   = _mm_permute_pd(t1, _GMX_MM_PERMUTE128D(1, 1));
    _mm_store_sd(ptrA, _mm_add_sd(_mm_load_sd(ptrA), t1));
    _mm_store_sd(ptrB, _mm_add_sd(_mm_load_sd(ptrB), t2));
}


#endif /* _kernelutil_x86_avx_256_double_h_ */
