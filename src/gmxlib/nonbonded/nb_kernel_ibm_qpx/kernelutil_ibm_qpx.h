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
#ifndef _kernelutil_ibm_qpx_h_
#define _kernelutil_ibm_qpx_h_


#include <math.h>

/* Normal sum of four qpx registers */
#define gmx_ibm_qpx_sum4(t0,t1,t2,t3)  vec_add(vec_add(t0,t1),vec_add(t2,t3))

#define GMX_IBM_QPX_TRANSPOSE4(_arg0,_arg1,_arg2,_arg3)   \
{                                                         \
    vector4double _t0,_t1,_t2,_t3,_unpacklo,_unpackhi;    \
    _unpacklo = vec_gpci(05140);                          \
    _unpackhi = vec_gpci(07362);                          \
    _t0       = vec_perm((_arg0),(_arg2),_unpacklo);      \
    _t1       = vec_perm((_arg0),(_arg2),_unpackhi);      \
    _t2       = vec_perm((_arg1),(_arg3),_unpacklo);      \
    _t3       = vec_perm((_arg1),(_arg3),_unpackhi);      \
    _arg0     = vec_perm(_t0,_t2,_unpacklo);              \
    _arg1     = vec_perm(_t0,_t2,_unpackhi);              \
    _arg2     = vec_perm(_t1,_t3,_unpacklo);              \
    _arg3     = vec_perm(_t1,_t3,_unpackhi);              \
}

static gmx_inline int
gmx_ibm_qpx_any_lt(vector4double a, vector4double b)
{
    vector4double t1;
    
    t1 = vec_cmplt(a,b);
    /* If comparison is true, sign will be positive */
    t1 = vec_and(t1,vec_sldw(t1,t1,2));
    t1 = vec_and(t1,vec_sldw(t1,t1,1));
    /* sign of first element is positive if it the comparison was true for any of the four */
    return (vec_extract(t1,0) > 0.0);
}

static gmx_inline vector4double
gmx_ibm_qpx_calc_rsq(vector4double dx, vector4double dy, vector4double dz)
{
    return vec_madd(dx,dx,vec_madd(dy,dy,vec_mul(dz,dz)));
}

/* Load a single value from 1-4 places, merge into qpx register */

static gmx_inline vector4double
gmx_ibm_qpx_load_4real_swizzle(const real * gmx_restrict ptrA,
                               const real * gmx_restrict ptrB,
                               const real * gmx_restrict ptrC,
                               const real * gmx_restrict ptrD)
{
    vector4double t1,t2,perm;

    perm = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    /* Erik, vec_lds requires that the EA be aligned */
    t1 = vec_perm(vec_lds(0,ptrA),vec_lds(0,ptrC),perm);
    t2 = vec_perm(vec_lds(0,ptrB),vec_lds(0,ptrD),perm);
    return vec_perm(t1,t2,perm);
}


static gmx_inline void
gmx_ibm_qpx_store_4real_swizzle(real * gmx_restrict ptrA,
                                real * gmx_restrict ptrB,
                                real * gmx_restrict ptrC,
                                real * gmx_restrict ptrD, vector4double qpx1)
{
    vector4double t2,t3,t4;

    t2       = vec_splat(qpx1,1);
    t3       = vec_splat(qpx1,2);
    t4       = vec_splat(qpx1,3);
    /* same - sts needs an aligned EA */
    vec_sts(qpx1,0,ptrA);
    vec_sts(t2,0,ptrB);
    vec_sts(t3,0,ptrC);
    vec_sts(t4,0,ptrD);
}


static gmx_inline void
gmx_ibm_qpx_increment_4real_swizzle(real * gmx_restrict ptrA,
                                    real * gmx_restrict ptrB,
                                    real * gmx_restrict ptrC,
                                    real * gmx_restrict ptrD, vector4double qpx1)
{
    vector4double tmp;

    tmp = gmx_ibm_qpx_load_4real_swizzle(ptrA,ptrB,ptrC,ptrD);
    tmp = vec_add(tmp,qpx1);
    gmx_ibm_qpx_store_4real_swizzle(ptrA,ptrB,ptrC,ptrD,tmp);
}


static gmx_inline void
gmx_ibm_qpx_load_4pair_swizzle(const real * gmx_restrict p1,
                               const real * gmx_restrict p2,
                               const real * gmx_restrict p3,
                               const real * gmx_restrict p4,
                               vector4double * gmx_restrict c6, vector4double * gmx_restrict c12)
{
    vector4double t1,t2,t3,t4,perm;
    perm = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */

    /* same, need aligned EA */
    t1   = vec_ld2(0,p1);
    t2   = vec_ld2(0,p2);
    t3   = vec_ld2(0,p3);
    t4   = vec_ld2(0,p4);
    t1   = vec_perm(t1,t3,perm);
    t2   = vec_perm(t2,t4,perm);
    *c6  = vec_perm(t1,t2,perm);
    *c12 = vec_perm(t1,t2,vec_gpci(07362)); /* Corresponds to _mm_unpackhi_ps */
}




static gmx_inline void
gmx_ibm_qpx_load_shift_and_1rvec_broadcast(const real * gmx_restrict xyz_shift,
                                           const real * gmx_restrict xyz,
                                           vector4double * gmx_restrict x1,
                                           vector4double * gmx_restrict y1,
                                           vector4double * gmx_restrict z1)
{
    vector4double t1,t2,t3,t4;

    t1   = vec_ld2(0,xyz_shift);
    t2   = vec_ld2(0,xyz);
    t3   = vec_lds(2,xyz_shift);
    t4   = vec_lds(2,xyz);
    t1   = vec_add(t1,t2);
    t3   = vec_add(t3,t4);

    *x1  = vec_splat(t1,0);
    *y1  = vec_splat(t1,1);
    *z1  = vec_splat(t3,0);
}


static gmx_inline void
gmx_ibm_qpx_load_shift_and_3rvec_broadcast(const real * gmx_restrict xyz_shift,
                                           const real * gmx_restrict xyz,
                                           vector4double * gmx_restrict x1, vector4double * gmx_restrict y1, vector4double * gmx_restrict z1,
                                           vector4double * gmx_restrict x2, vector4double * gmx_restrict y2, vector4double * gmx_restrict z2,
                                           vector4double * gmx_restrict x3, vector4double * gmx_restrict y3, vector4double * gmx_restrict z3)
{
    vector4double tA,tB,tC;
    vector4double t1,t2,t3,t4,t5,t6;

    /* same */
    tA   = vec_ld2(0,xyz_shift);    /* sX sY sX sY */
    tB   = vec_lds(2,xyz_shift);  /* sZ  -  -  - */

    t1   = vec_ld(0,xyz);
    t2   = vec_ld(4,xyz);
    t3   = vec_lds(8,xyz);

    tC   = vec_sldw(tA,tB,1); /* sY sX sY sZ */
    
    t4   = vec_sldw(tC,tA,1); /* sX sY sZ sX */
    t5   = vec_sldw(tC,tA,2); /* sY sZ sX sY */
    t6   = tB;                /* sZ  -  -  - */

    t1   = vec_add(t1,t4);
    t2   = vec_add(t2,t5);
    t3   = vec_add(t3,t6);

    *x1  = vec_splat(t1,0);
    *y1  = vec_splat(t1,1);
    *z1  = vec_splat(t1,2);
    *x2  = vec_splat(t1,3);
    *y2  = vec_splat(t2,0);
    *z2  = vec_splat(t2,1);
    *x3  = vec_splat(t2,2);
    *y3  = vec_splat(t2,3);
    *z3  = vec_splat(t3,0);
}


static gmx_inline void
gmx_ibm_qpx_load_shift_and_4rvec_broadcast(const real * gmx_restrict xyz_shift,
                                           const real * gmx_restrict xyz,
                                           vector4double * gmx_restrict x1, vector4double * gmx_restrict y1, vector4double * gmx_restrict z1,
                                           vector4double * gmx_restrict x2, vector4double * gmx_restrict y2, vector4double * gmx_restrict z2,
                                           vector4double * gmx_restrict x3, vector4double * gmx_restrict y3, vector4double * gmx_restrict z3,
                                           vector4double * gmx_restrict x4, vector4double * gmx_restrict y4, vector4double * gmx_restrict z4)
{
    vector4double tA,tB,tC;
    vector4double t1,t2,t3,t4,t5,t6;

    /* same */
    tA   = vec_ld2(0,xyz_shift);    /* sX sY sX sY */
    tB   = vec_lds(2,xyz_shift);  /* sZ  -  -  - */

    t1   = vec_ld(0,xyz);
    t2   = vec_ld(4,xyz);
    t3   = vec_ld(8,xyz);

    tC   = vec_sldw(tA,tB,1); /* sY sX sY sZ */
    
    t4   = vec_sldw(tC,tA,1); /* sX sY sZ sX */
    t5   = vec_sldw(tC,tA,2); /* sY sZ sX sY */
    t6   = vec_sldw(tC,t4,3); /* sZ sX sY sZ */

    t1   = vec_add(t1,t4);
    t2   = vec_add(t2,t5);
    t3   = vec_add(t3,t6);

    *x1  = vec_splat(t1,0);
    *y1  = vec_splat(t1,1);
    *z1  = vec_splat(t1,2);
    *x2  = vec_splat(t1,3);
    *y2  = vec_splat(t2,0);
    *z2  = vec_splat(t2,1);
    *x3  = vec_splat(t2,2);
    *y3  = vec_splat(t2,3);
    *z3  = vec_splat(t3,0);
    *x4  = vec_splat(t3,1);
    *y4  = vec_splat(t3,2);
    *z4  = vec_splat(t3,3);
}


static gmx_inline void
gmx_ibm_qpx_load_1rvec_4ptr_swizzle(const real * gmx_restrict ptrA, const real * gmx_restrict ptrB,
                                    const real * gmx_restrict ptrC, const real * gmx_restrict ptrD,
                                    vector4double * gmx_restrict x1, vector4double * gmx_restrict y1, vector4double * gmx_restrict z1)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8,unpacklo,unpackhi;
    /* same */
    t1             = vec_ld2(0,ptrA);    /* xA yA xA yA */
    t2             = vec_ld2(0,ptrB);    /* xB yB xB yB */
    t3             = vec_ld2(0,ptrC);    /* xC yC xC yC */
    t4             = vec_ld2(0,ptrD);    /* xD yD xD yD */

    t5             = vec_lds(2,ptrA);  /* zA  -  -  - */
    t6             = vec_lds(2,ptrB);  /* zB  -  -  - */
    t7             = vec_lds(2,ptrC);  /* zC  -  -  - */
    t8             = vec_lds(2,ptrD);  /* zD  -  -  - */

    t1             = vec_sldw(t1,t5,2); /* xA yA zA - */
    t2             = vec_sldw(t2,t6,2); /* xB yB zB - */
    t3             = vec_sldw(t3,t7,2); /* xC yC zC - */
    t4             = vec_sldw(t4,t8,2); /* xD yD zD - */
    
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */

    t5             = vec_perm(t1,t3,unpacklo); /* xA xC yA yC */
    t6             = vec_perm(t2,t4,unpacklo); /* xB xD yB yD */
    t7             = vec_perm(t1,t3,unpackhi); /* zA zC  -  - */
    t8             = vec_perm(t2,t4,unpackhi); /* zB zD  -   - */
    
    *x1            = vec_perm(t5,t6,unpacklo); /* xA xB xC xD */
    *y1            = vec_perm(t5,t6,unpackhi); /* yA yB yC yD */
    *z1            = vec_perm(t7,t8,unpacklo); /* zA zB zC zD */
}


static gmx_inline void
gmx_ibm_qpx_load_3rvec_4ptr_swizzle(const real * gmx_restrict ptrA, const real * gmx_restrict ptrB,
                                    const real * gmx_restrict ptrC, const real * gmx_restrict ptrD,
                                    vector4double * gmx_restrict x1, vector4double * gmx_restrict y1, vector4double * gmx_restrict z1,
                                    vector4double * gmx_restrict x2, vector4double * gmx_restrict y2, vector4double * gmx_restrict z2,
                                    vector4double * gmx_restrict x3, vector4double * gmx_restrict y3, vector4double * gmx_restrict z3)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8,unpacklo,unpackhi;
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */

    t1            = vec_ld(0,ptrA);
    t2            = vec_ld(0,ptrB);
    t3            = vec_ld(0,ptrC);
    t4            = vec_ld(0,ptrD);
    t5            = vec_perm(t1,t3,unpacklo);
    t6            = vec_perm(t2,t4,unpacklo);
    t7            = vec_perm(t1,t3,unpackhi);
    t8            = vec_perm(t2,t4,unpackhi);
    *x1           = vec_perm(t5,t6,unpacklo);
    *y1           = vec_perm(t5,t6,unpackhi);
    *z1           = vec_perm(t7,t8,unpacklo);
    *x2           = vec_perm(t7,t8,unpackhi);
    
    t1            = vec_ld(4,ptrA);
    t2            = vec_ld(4,ptrB);
    t3            = vec_ld(4,ptrC);
    t4            = vec_ld(4,ptrD);
    t5            = vec_perm(t1,t3,unpacklo);
    t6            = vec_perm(t2,t4,unpacklo);
    t7            = vec_perm(t1,t3,unpackhi);
    t8            = vec_perm(t2,t4,unpackhi);
    *y2           = vec_perm(t5,t6,unpacklo);
    *z2           = vec_perm(t5,t6,unpackhi);
    *x3           = vec_perm(t7,t8,unpacklo);
    *y3           = vec_perm(t7,t8,unpackhi);
    
    t1            = vec_lds(8,ptrA);
    t2            = vec_lds(8,ptrB);
    t3            = vec_lds(8,ptrC);
    t4            = vec_lds(8,ptrD);
    t1            = vec_perm(t1,t3,unpacklo);
    t3            = vec_perm(t2,t4,unpacklo);
    *z3           = vec_perm(t1,t3,unpacklo);
}


static gmx_inline void
gmx_ibm_qpx_load_4rvec_4ptr_swizzle(const real * gmx_restrict ptrA, const real * gmx_restrict ptrB,
                                    const real * gmx_restrict ptrC, const real * gmx_restrict ptrD,
                                    vector4double * gmx_restrict x1, vector4double * gmx_restrict y1, vector4double * gmx_restrict z1,
                                    vector4double * gmx_restrict x2, vector4double * gmx_restrict y2, vector4double * gmx_restrict z2,
                                    vector4double * gmx_restrict x3, vector4double * gmx_restrict y3, vector4double * gmx_restrict z3,
                                    vector4double * gmx_restrict x4, vector4double * gmx_restrict y4, vector4double * gmx_restrict z4)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8,unpacklo,unpackhi;
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */

    t1            = vec_ld(0,ptrA);
    t2            = vec_ld(0,ptrB);
    t3            = vec_ld(0,ptrC);
    t4            = vec_ld(0,ptrD);
    t5            = vec_perm(t1,t3,unpacklo);
    t6            = vec_perm(t2,t4,unpacklo);
    t7            = vec_perm(t1,t3,unpackhi);
    t8            = vec_perm(t2,t4,unpackhi);
    *x1           = vec_perm(t5,t6,unpacklo);
    *y1           = vec_perm(t5,t6,unpackhi);
    *z1           = vec_perm(t7,t8,unpacklo);
    *x2           = vec_perm(t7,t8,unpackhi);
    
    t1            = vec_ld(4,ptrA);
    t2            = vec_ld(4,ptrB);
    t3            = vec_ld(4,ptrC);
    t4            = vec_ld(4,ptrD);
    t5            = vec_perm(t1,t3,unpacklo);
    t6            = vec_perm(t2,t4,unpacklo);
    t7            = vec_perm(t1,t3,unpackhi);
    t8            = vec_perm(t2,t4,unpackhi);
    *y2           = vec_perm(t5,t6,unpacklo);
    *z2           = vec_perm(t5,t6,unpackhi);
    *x3           = vec_perm(t7,t8,unpacklo);
    *y3           = vec_perm(t7,t8,unpackhi);

    t1            = vec_ld(8,ptrA);
    t2            = vec_ld(8,ptrB);
    t3            = vec_ld(8,ptrC);
    t4            = vec_ld(8,ptrD);
    t5            = vec_perm(t1,t3,unpacklo);
    t6            = vec_perm(t2,t4,unpacklo);
    t7            = vec_perm(t1,t3,unpackhi);
    t8            = vec_perm(t2,t4,unpackhi);
    *z3           = vec_perm(t5,t6,unpacklo);
    *x4           = vec_perm(t5,t6,unpackhi);
    *y4           = vec_perm(t7,t8,unpacklo);
    *z4           = vec_perm(t7,t8,unpackhi);
}



static gmx_inline void
gmx_ibm_qpx_decrement_1rvec_4ptr_swizzle(real * gmx_restrict ptrA, real * gmx_restrict ptrB,
                                         real * gmx_restrict ptrC, real * gmx_restrict ptrD,
                                         vector4double x1, vector4double y1, vector4double z1)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8,t9;
    vector4double unpacklo,unpackhi;
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */

    t5          = vec_perm(x1,y1,unpacklo); /* xA yA xB yB */
    t6          = vec_perm(x1,y1,unpackhi); /* xC yC xD yD */
    t7          = vec_splat(z1,1);
    t8          = vec_splat(z1,2);
    t9          = vec_splat(z1,3);
    
    t1          = vec_sldw(vec_ld2(0,ptrA),vec_ld2(0,ptrB),2);
    t2          = vec_sldw(vec_ld2(0,ptrC),vec_ld2(0,ptrD),2);
    
    t1          = vec_sub(t1,t5);
    t2          = vec_sub(t2,t6);
    vec_st2(t1,0,ptrA);
    vec_st2(t2,0,ptrC);
    vec_sldw(t1,t1,2);
    vec_sldw(t2,t2,2);
    vec_st2(t1,0,ptrB);
    vec_st2(t2,0,ptrD);
    
    vec_sts( vec_sub(vec_lds(2,ptrA),z1), 2, ptrA);
    vec_sts( vec_sub(vec_lds(2,ptrB),t7), 2, ptrB);
    vec_sts( vec_sub(vec_lds(2,ptrC),t8), 2, ptrC);
    vec_sts( vec_sub(vec_lds(2,ptrD),t9), 2, ptrD);
}


static gmx_inline void
gmx_ibm_qpx_decrement_3rvec_4ptr_swizzle(real * gmx_restrict ptrA, real * gmx_restrict ptrB,
                                         real * gmx_restrict ptrC, real * gmx_restrict ptrD,
                                         vector4double x1, vector4double y1, vector4double z1,
                                         vector4double x2, vector4double y2, vector4double z2,
                                         vector4double x3, vector4double y3, vector4double z3)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8;
    vector4double tA,tB,tC,tD;
    vector4double unpacklo,unpackhi;
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */
    
    t1            = vec_perm(x1,z1,unpacklo);
    t2            = vec_perm(y1,x2,unpacklo);
    t3            = vec_perm(x1,z1,unpackhi);
    t4            = vec_perm(y1,x2,unpackhi);
    t5            = vec_perm(t1,t2,unpacklo);
    t6            = vec_perm(t1,t2,unpackhi);
    t7            = vec_perm(t3,t4,unpacklo);
    t8            = vec_perm(t3,t4,unpackhi);
    
    tA            = vec_ld(0,ptrA);
    tB            = vec_ld(0,ptrB);
    tC            = vec_ld(0,ptrC);
    tD            = vec_ld(0,ptrD);
    tA            = vec_sub(tA,t5);
    tB            = vec_sub(tB,t6);
    tC            = vec_sub(tC,t7);
    tD            = vec_sub(tD,t8);
    vec_st(tA,0,ptrA);
    vec_st(tB,0,ptrB);
    vec_st(tC,0,ptrC);
    vec_st(tD,0,ptrD);
    
    t1            = vec_perm(y2,x3,unpacklo);
    t2            = vec_perm(z2,y3,unpacklo);
    t3            = vec_perm(y2,x3,unpackhi);
    t4            = vec_perm(z2,y3,unpackhi);
    t5            = vec_perm(t1,t2,unpacklo);
    t6            = vec_perm(t1,t2,unpackhi);
    t7            = vec_perm(t3,t4,unpacklo);
    t8            = vec_perm(t3,t4,unpackhi);
    
    tA            = vec_ld(4,ptrA);
    tB            = vec_ld(4,ptrB);
    tC            = vec_ld(4,ptrC);
    tD            = vec_ld(4,ptrD);
    tA            = vec_sub(tA,t5);
    tB            = vec_sub(tB,t6);
    tC            = vec_sub(tC,t7);
    tD            = vec_sub(tD,t8);
    vec_st(tA,4,ptrA);
    vec_st(tB,4,ptrB);
    vec_st(tC,4,ptrC);
    vec_st(tD,4,ptrD);
    
    t2          = vec_splat(z3,1);
    t3          = vec_splat(z3,2);
    t4          = vec_splat(z3,3);
    vec_sts( vec_sub(vec_lds(8,ptrA),z3), 8, ptrA);
    vec_sts( vec_sub(vec_lds(8,ptrB),t2), 8, ptrB);
    vec_sts( vec_sub(vec_lds(8,ptrC),t3), 8, ptrC);
    vec_sts( vec_sub(vec_lds(8,ptrD),t4), 8, ptrD);

}


static gmx_inline void
gmx_ibm_qpx_decrement_4rvec_4ptr_swizzle(real * gmx_restrict ptrA, real * gmx_restrict ptrB,
                                         real * gmx_restrict ptrC, real * gmx_restrict ptrD,
                                         vector4double x1, vector4double y1, vector4double z1,
                                         vector4double x2, vector4double y2, vector4double z2,
                                         vector4double x3, vector4double y3, vector4double z3,
                                         vector4double x4, vector4double y4, vector4double z4)
{
    vector4double t1,t2,t3,t4,t5,t6,t7,t8;
    vector4double tA,tB,tC,tD;
    vector4double unpacklo,unpackhi;
    unpacklo = vec_gpci(05140); /* Corresponds to _mm_unpacklo_ps */
    unpackhi = vec_gpci(07362); /* Corresponds to _mm_unpackhi_ps */
    
    t1            = vec_perm(x1,z1,unpacklo);
    t2            = vec_perm(y1,x2,unpacklo);
    t3            = vec_perm(x1,z1,unpackhi);
    t4            = vec_perm(y1,x2,unpackhi);
    t5            = vec_perm(t1,t2,unpacklo);
    t6            = vec_perm(t1,t2,unpackhi);
    t7            = vec_perm(t3,t4,unpacklo);
    t8            = vec_perm(t3,t4,unpackhi);
    
    tA            = vec_ld(0,ptrA);
    tB            = vec_ld(0,ptrB);
    tC            = vec_ld(0,ptrC);
    tD            = vec_ld(0,ptrD);
    tA            = vec_sub(tA,t5);
    tB            = vec_sub(tB,t6);
    tC            = vec_sub(tC,t7);
    tD            = vec_sub(tD,t8);
    vec_st(tA,0,ptrA);
    vec_st(tB,0,ptrB);
    vec_st(tC,0,ptrC);
    vec_st(tD,0,ptrD);
    
    t1            = vec_perm(y2,x3,unpacklo);
    t2            = vec_perm(z2,y3,unpacklo);
    t3            = vec_perm(y2,x3,unpackhi);
    t4            = vec_perm(z2,y3,unpackhi);
    t5            = vec_perm(t1,t2,unpacklo);
    t6            = vec_perm(t1,t2,unpackhi);
    t7            = vec_perm(t3,t4,unpacklo);
    t8            = vec_perm(t3,t4,unpackhi);
    
    tA            = vec_ld(4,ptrA);
    tB            = vec_ld(4,ptrB);
    tC            = vec_ld(4,ptrC);
    tD            = vec_ld(4,ptrD);
    tA            = vec_sub(tA,t5);
    tB            = vec_sub(tB,t6);
    tC            = vec_sub(tC,t7);
    tD            = vec_sub(tD,t8);
    vec_st(tA,4,ptrA);
    vec_st(tB,4,ptrB);
    vec_st(tC,4,ptrC);
    vec_st(tD,4,ptrD);
    
    t1            = vec_perm(z3,y4,unpacklo);
    t2            = vec_perm(x4,z4,unpacklo);
    t3            = vec_perm(z3,y4,unpackhi);
    t4            = vec_perm(x4,z4,unpackhi);
    t5            = vec_perm(t1,t2,unpacklo);
    t6            = vec_perm(t1,t2,unpackhi);
    t7            = vec_perm(t3,t4,unpacklo);
    t8            = vec_perm(t3,t4,unpackhi);
    
    tA            = vec_ld(8,ptrA);
    tB            = vec_ld(8,ptrB);
    tC            = vec_ld(8,ptrC);
    tD            = vec_ld(8,ptrD);
    tA            = vec_sub(tA,t5);
    tB            = vec_sub(tB,t6);
    tC            = vec_sub(tC,t7);
    tD            = vec_sub(tD,t8);
    vec_st(tA,8,ptrA);
    vec_st(tB,8,ptrB);
    vec_st(tC,8,ptrC);
    vec_st(tD,8,ptrD);
}


static gmx_inline void
gmx_ibm_qpx_update_iforce_1atom_swizzle(vector4double fix1, vector4double fiy1, vector4double fiz1,
                                        real * gmx_restrict fptr,
                                        real * gmx_restrict fshiftptr)
{
    vector4double t1;
    vector4double tA,tB,tC,tD;
    
    /* transpose data */
    t1   = vec_xor(fix1,fix1); /* 0.0 */
    GMX_IBM_QPX_TRANSPOSE4(fix1,fiy1,fiz1,t1);
    fix1 = vec_add(vec_add(fix1,fiy1), vec_add(fiz1,t1));

    tA   = vec_sldw(vec_ld2(0,fptr),vec_lds(2,fptr),2);
    tB   = vec_sldw(vec_ld2(0,fshiftptr),vec_lds(2,fshiftptr),2);

    tA   = vec_add(tA,fix1);
    tB   = vec_add(tB,fix1);
    
    vec_st2(tA,0,fptr);
    vec_sts(vec_sldw(tA,tA,2),2,fptr);
    vec_st2(tB,0,fshiftptr);
    vec_sts(vec_sldw(tB,tB,2),2,fshiftptr);
}


static gmx_inline void
gmx_ibm_qpx_update_iforce_3atom_swizzle(vector4double fix1, vector4double fiy1, vector4double fiz1,
                                        vector4double fix2, vector4double fiy2, vector4double fiz2,
                                        vector4double fix3, vector4double fiy3, vector4double fiz3,
                                        real * gmx_restrict fptr,
                                        real * gmx_restrict fshiftptr)
{
    vector4double t1,t2,t3,t4;
    
    /* transpose data */
    GMX_IBM_QPX_TRANSPOSE4(fix1,fiy1,fiz1,fix2);
    GMX_IBM_QPX_TRANSPOSE4(fiy2,fiz2,fix3,fiy3);

    t1   = vec_splat(fiz3,1);
    t2   = vec_splat(fiz3,2);
    t3   = vec_splat(fiz3,3);
    fix1 = vec_add(vec_add(fix1,fiy1), vec_add(fiz1,fix2));  /* x1 y1 z1 x2 */
    fiy2 = vec_add(vec_add(fiy2,fiz2), vec_add(fix3,fiy3));  /* y2 z2 x3 y3 */
    fiz3 = vec_add(vec_add(fiz3,t1)  , vec_add(t2,t3));      /* z3  -  -  - */

    vec_st( vec_add(vec_ld(0,fptr),fix1),0,fptr);
    vec_st( vec_add(vec_ld(4,fptr),fiy2),4,fptr);
    vec_sts(vec_add(vec_lds(8,fptr),fiz3),8,fptr);

    t4   = vec_sldw(vec_ld2(0,fshiftptr),vec_lds(2,fshiftptr),2);

    /* fix1:   x1 y1 z1 - */
    t1   = vec_sldw(fix1,fiy2,3); /* x2 y2 z2 - */
    t2   = vec_sldw(fiy2,fiz3,2); /* x3 y3 z3 - */
    
    t4   = vec_add(vec_add(t4,fix1),vec_add(t1,t2)); /* Xsum Ysum Zsum - */
    vec_st2(t4,0,fshiftptr);
    vec_sts(vec_sldw(t4,t4,2),2,fshiftptr);
}


static gmx_inline void
gmx_ibm_qpx_update_iforce_4atom_swizzle(vector4double fix1, vector4double fiy1, vector4double fiz1,
                                        vector4double fix2, vector4double fiy2, vector4double fiz2,
                                        vector4double fix3, vector4double fiy3, vector4double fiz3,
                                        vector4double fix4, vector4double fiy4, vector4double fiz4,
                                        real * gmx_restrict fptr,
                                        real * gmx_restrict fshiftptr)
{
    vector4double t1,t2,t3,t4;
    
    /* transpose data */
    GMX_IBM_QPX_TRANSPOSE4(fix1,fiy1,fiz1,fix2);
    GMX_IBM_QPX_TRANSPOSE4(fiy2,fiz2,fix3,fiy3);
    GMX_IBM_QPX_TRANSPOSE4(fiz3,fix4,fiy4,fiz4);
    
    fix1 = vec_add(vec_add(fix1,fiy1), vec_add(fiz1,fix2));  /* x1 y1 z1 x2 */
    fiy2 = vec_add(vec_add(fiy2,fiz2), vec_add(fix3,fiy3));  /* y2 z2 x3 y3 */
    fiz3 = vec_add(vec_add(fiz3,fix4), vec_add(fiy4,fiz4));  /* z3 x4 y4 z4 */
    
    vec_st( vec_add(vec_ld(0,fptr),fix1),0,fptr);
    vec_st( vec_add(vec_ld(4,fptr),fiy2),4,fptr);
    vec_st( vec_add(vec_ld(8,fptr),fiz3),8,fptr);

    t4   = vec_sldw(vec_ld2(0,fshiftptr),vec_lds(2,fshiftptr),2);
    /* fix1:   x1 y1 z1 - */
    t1   = vec_sldw(fix1,fiy2,3); /* x2 y2 z2 - */
    t2   = vec_sldw(fiy2,fiz3,2); /* x3 y3 z3 - */
    t3   = vec_sldw(fiz3,fiz3,1); /* x4 y4 z4 - */
    
    t4   = vec_add(t4,vec_add(vec_add(fix1,t1),vec_add(t2,t3))); /* Xsum Ysum Zsum - */
    vec_st2(t4,0,fshiftptr);
    vec_sts(vec_sldw(t4,t4,2),2,fshiftptr);
}



static gmx_inline void
gmx_ibm_qpx_update_1pot(vector4double pot1, real * gmx_restrict ptrA)
{
    pot1 = vec_add(pot1,vec_sldw(pot1,pot1,2));
    pot1 = vec_add(pot1,vec_sldw(pot1,pot1,1));
    vec_sts( vec_add(vec_lds(0,ptrA),pot1), 0, ptrA);
}


static gmx_inline void
gmx_ibm_qpx_update_2pot(vector4double pot1, real * gmx_restrict ptrA,
                        vector4double pot2, real * gmx_restrict ptrB)
{
    vector4double t1,t2,t3,t4;

    t1   = vec_sldw(pot1,pot2,2); /* p1C p1D p2A p2B */
    t2   = vec_sldw(pot2,pot1,2); /* p2C p2D p1A p1B */
    t2   = vec_sldw(t2,t2,2);     /* p1A p1B p2C p2D */
    
    t1   = vec_add(t1,t2);  /* p1A+C p1B+D p2A+C p2B+D */
    t1   = vec_add(t1,vec_sldw(t1,t1,1));  /* p1A+B+C+D  -  p2A+B+C+D - */
    t2   = vec_splat(t1,2); /* p2A+B+C+D - - - */
                            
    vec_sts( vec_add(vec_lds(0,ptrA),t1), 0, ptrA);
    vec_sts( vec_add(vec_lds(0,ptrB),t2), 0, ptrB);
}

#endif /* _kernelutil_ibm_qpx_h_ */
