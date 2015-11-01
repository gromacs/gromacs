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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_DOUBLE_H

#include <math.h>

#include <altivec.h>

#include "impl_ibm_vsx_common.h"

/* IBM VSX SIMD instruction wrappers. Power7 and later.
 *
 * While this instruction set is similar to VMX, there are quite a few differences
 * that make it easier to understand if we start from scratch rather than by
 * including the VMX implementation and changing lots of things.
 */


/* Make sure we do not screw up c++ - undefine vector/bool, and rely on __vector,
 * which is present both on gcc and xlc.
 */
#undef vector

/* g++ is also unhappy with the clash of vector bool and the C++ reserved 'bool',
 * which is solved by undefining bool and reyling on __bool. However, that does
 * not work with xlc, which requires us to use bool. Solve the conflict by
 * defining a new vsx_bool.
 */
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
#    define vsx_bool __bool
#    undef  bool
#else
#    define vsx_bool bool
#endif


/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define SimdDouble          __vector double
#define simdLoadD(m)         (*(const SimdDouble *)(m))
#define simdStoreD(m, x)      { *(SimdDouble *)(m) = (x); }
#define simdLoad1D(m)        vec_splats((double)(*m))
#define simdSet1D(x)         vec_splats((double)(x))
#if defined(__ibmxl__) || defined(__xlC__)
#    define simdLoadUD(m)    vec_xld2(0, (double *)(m))
#    define simdStoreUD(m, x) vec_xstd2(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define simdLoadUD       simdLoadD
#    define simdStoreUD      simdStoreD
#endif
#define simdSetZeroD()       vec_splats(0.0)
#define simdAddD(a, b)       vec_add(a, b)
#define simdSubD(a, b)       vec_sub(a, b)
#define simdMulD(a, b)       vec_mul(a, b)
#define simdFmaddD(a, b, c)  vec_madd(a, b, c)
#define simdFmsubD(a, b, c)  vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define simdFnmaddD(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define simdFnmsubD(a, b, c) vec_nmadd(a, b, c)
#define simdAndD(a, b)       vec_and(a, b)
#define simdAndNotD(a, b)    vec_andc(b, a)
#define simdOrD(a, b)        vec_or(a, b)
#define simdXorD(a, b)       vec_xor(a, b)
#define simdRsqrtD(a)        vec_rsqrte(a)
#define simdRcpD(a)          vec_re(a)
#define simdAbsD(a)         vec_abs(a)
#define simdNegD(a)         (-(a))
#define simdMaxD(a, b)       vec_max(a, b)
#define simdMinD(a, b)       vec_min(a, b)
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
/* gcc up to at least version 4.9 does not support vec_round() in double precision. */
#    define simdRoundD(a)    ({ __vector double res; __asm__ ("xvrdpi %0,%1" : "=ww" (res) : "ww" ((__vector double) (a))); res; })
#else
/* IBM xlC */
#    define simdRoundD(a)    vec_round(a)
#endif
#define simdTruncD(a)        vec_trunc(a)
#define simdFractionD(x)     vec_sub(x, vec_trunc(x))
#define simdGetExponentD(a) simdGetExponentD_ibm_vsx(a)
#define simdGetMantissaD(a) simdGetMantissaD_ibm_vsx(a)
#define simdSetExponentD(a) simdSetExponentD_ibm_vsx(a)
/* integer datatype corresponding to double: SimdDInt32 */
#define SimdDInt32          __vector signed int
#define simdLoadDI(m)        simdLoadDI_ibm_vsx(m)
#define simdStoreDI(m, x)    simdStoreDI_ibm_vsx(m, x)
#define simdSet1DI(i)        vec_splats((int)(i))
#define simdLoadUDI          simdLoadDI
#define simdStoreUDI         simdStoreDI
#define simdSetZeroDI()      vec_splats((int)0)
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
/* gcc up to at least version 4.9 is missing intrinsics for double precision
 * to integer conversion, use inline asm instead.
 */
#    define simdCvttD2I(a)   simdCvttD2I_ibm_vsx(a)
#    define simdCvtI2D(a)    simdCvtI2D_ibm_vsx(a)
#else
/* IBM xlC */
#    define simdCvttD2I(a)       vec_cts(a, 0)
#    define simdCvtI2D(a)        vec_ctd(a, 0)
#endif
#define simdCvtD2I(a)         simdCvttD2I(simdRoundD(a))
#define simdExtractDI(a, i)   simdExtractFI_ibm_vsx(a, (i)*2)
/* Integer logical ops on SimdDInt32 */
#define simdSlliDI(a, i)      vec_sl(a, vec_splats((unsigned int)(i)))
#define simdSrliDI(a, i)      vec_sr(a, vec_splats((unsigned int)(i)))
#define simdAndDI(a, b)       vec_and(a, b)
#define simdAndNotDI(a, b)    vec_andc(b, a)
#define simdOrDI(a, b)        vec_or(a, b)
#define simdXorDI(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on SimdDInt32 */
#define simdAddDI(a, b)       vec_add(a, b)
#define simdSubDI(a, b)       vec_sub(a, b)
#define simdMulDI(a, b)       ((a)*(b))
/* Boolean & comparison operations on SimdDouble */
#define SimdDBool           __vector vsx_bool long long
#define simdCmpEqD(a, b)     vec_cmpeq(a, b)
#define simdCmpLtD(a, b)     vec_cmplt(a, b)
#define simdCmpLeD(a, b)     vec_cmple(a, b)
#define simdAndDB(a, b)      (__vector vsx_bool long long)vec_and((__vector signed int)a, (__vector signed int)b)
#define simdOrDB(a, b)       (__vector vsx_bool long long)vec_or((__vector signed int)a, (__vector signed int)b)
#define simdAnyTrueDB(a)     vec_any_ne((__vector vsx_bool int)a, (__vector vsx_bool int)vec_splats(0))
#define simdMaskD(a, sel)    vec_and(a, (__vector double)sel)
#define simdMaskNotD(a, sel) vec_andc(a, (__vector double)sel)
#define simdBlendD(a, b, sel)    vec_sel(a, b, sel)
#define simdReduceD(a)       simdReduceD_ibm_vsx(a)
/* Boolean & comparison operations on SimdFInt32 */
#define SimdDIBool          __vector vsx_bool int
#define simdCmpEqDI(a, b)    vec_cmpeq(a, b)
#define simdCmpLtDI(a, b)    vec_cmplt(a, b)
#define simdAndDIB(a, b)     vec_and(a, b)
#define simdOrDIB(a, b)      vec_or(a, b)
/* Since we have applied all operations to pairs of elements we can work on all elements here */
#define simdAnyTrueDIB(a)          vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define simdMaskDI(a, sel)    vec_and(a, (__vector signed int)sel)
#define simdMaskNotDI(a, sel) vec_andc(a, (__vector signed int)sel)
#define simdBlendDI(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define simdCvtDB2DIB(x)     (__vector vsx_bool int)(x)
#define simdCvtDIB2DB(x)     (__vector vsx_bool long long)(x)
/* Float/double conversion */
#define simdCvtF2DD(f, d0, d1)  simdCvtF2DD_ibm_vsx(f, d0, d1)
#define simdCvtDD2F(d0, d1)     simdCvtDD2F_ibm_vsx(d0, d1)

#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
/* gcc-4.9 is missing double-to-float/float-to-double conversions. */
#    define gmx_vsx_f2d(x) ({ __vector double res; __asm__ ("xvcvspdp %0,%1" : "=ww" (res) : "ww" ((__vector float) (x))); res; })
#    define gmx_vsx_d2f(x) ({ __vector float res; __asm__ ("xvcvdpsp %0,%1" : "=ww" (res) : "ww" ((__vector double) (x))); res; })
#else
/* f2d and d2f are indeed identical on xlC; it is selected by the argument and result type. */
#    define gmx_vsx_f2d(x)       vec_cvf(x)
#    define gmx_vsx_d2f(x)       vec_cvf(x)
#endif



/****************************************************
 * DOUBLE PREC. IMPLEMENTATION HELPER FUNCTIONS     *
 ****************************************************/
static inline SimdDInt32
simdLoadDI_ibm_vsx(const int *m)
{
#ifdef __xlC__
    /* old xlc version 12 does not understand long long VSX instructions */
    __vector signed int          t0, t1;
    const __vector unsigned char perm = { 0, 1, 2, 3, 0, 1, 2, 3, 16, 17, 18, 19, 16, 17, 18, 19 };
    t0 = vec_splats(m[0]);
    t1 = vec_splats(m[1]);
    return vec_perm(t0, t1, perm);
#else
    __vector long long int t0;
    t0 = vec_splats(*(long long int *)m);
    return vec_mergeh((__vector signed int)t0, (__vector signed int)t0);
#endif
}

static inline void
simdStoreDI_ibm_vsx(int *m, SimdDInt32 x)
{
#ifdef __xlC__
    /* old xlc version 12 does not understand long long VSX instructions */
    m[0] = vec_extract(x, 0);
    m[1] = vec_extract(x, 2);
#else
    __vector unsigned char perm = { 0, 1, 2, 3, 8, 9, 10, 11, 0, 1, 2, 3, 8, 9, 10, 11 };
    x                   = vec_perm(x, x, perm);
    *(long long int *)m = vec_extract((__vector long long int)x, 0);
#endif
}

#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
static inline SimdDInt32
simdCvttD2I_ibm_vsx(SimdDouble x)
{
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
    SimdDInt32                   ix;

    __asm__ ("xvcvdpsxws %0,%1" : "=ww" (ix) : "ww" ((__vector double) (x)));

    return vec_perm(ix, ix, perm);
}

static inline SimdDouble
simdCvtI2D_ibm_vsx(SimdDInt32 ix)
{
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
    SimdDouble                   x;
    ix = vec_perm(ix, ix, perm);
    __asm__ ("xvcvsxwdp %0,%1" : "=ww" (x) : "ww" ((__vector signed int) (ix)));
    return x;
}
#endif


static inline SimdDouble
simdGetExponentD_ibm_vsx(SimdDouble x)
{
#ifndef __BIG_ENDIAN__
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
#endif
    SimdDouble                   expmask = (__vector double)vec_splats(0x7FF0000000000000ULL);
    SimdDInt32                   i1023   = vec_splats(1023);
    SimdDInt32                   iexp;

    iexp = (__vector signed int)simdAndD(x, expmask);
    /* The data is in the upper half of each double (corresponding to elements 1/3).
     * First shift 52-32=20bits, and then permute to swap 0 with 1 and 2 with 3
     * For big endian they are in opposite order, so we avoid the swap.
     */
    iexp = simdSrliFI(iexp, 20);
#ifndef __BIG_ENDIAN__
    iexp = vec_perm(iexp, iexp, perm);
#endif
    iexp = vec_sub(iexp, i1023);
    /* Now we have the correct integer in elements 0 & 2. Never mind about elements 1,3 */
    return simdCvtI2D(iexp);
}

static inline SimdDouble
simdGetMantissaD_ibm_vsx(SimdDouble x)
{
    SimdDouble  expmask = (__vector double)vec_splats(0x7FF0000000000000ULL);

    x = simdAndNotD(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return simdOrD(x, vec_splats(1.0));
}

static inline SimdDouble
simdSetExponentD_ibm_vsx(SimdDouble x)
{
    SimdDInt32                   iexp  = simdCvtD2I(x);
    SimdDInt32                   i1023 = vec_splats(1023);
#ifdef __BIG_ENDIAN__
    const __vector unsigned char perm = {0, 1, 2, 3, 16, 17, 18, 19, 8, 9, 10, 11, 16, 17, 18, 19};
#else
    const __vector unsigned char perm = {16, 17, 18, 19, 0, 1, 2, 3, 16, 17, 18, 19, 8, 9, 10, 11};
#endif

    iexp = vec_add(iexp, i1023);
    /* exponent is now present in pairs of integers; 0011.
     * Elements 0/2 already correspond to the upper half of each double,
     * so we only need to shift by another 52-32=20 bits.
     * The remaining elements are set to zero.
     */
    iexp = vec_sl(iexp, (__vector unsigned int)vec_splats((int)20));
    iexp = vec_perm(iexp, vec_splats(0), perm);
    return (__vector double)iexp;
}

static inline double
simdReduceD_ibm_vsx(SimdDouble x)
{
    const __vector unsigned char perm = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
#ifdef __xlC__
    /* old xlc version 12 does not understand vec_perm() with double arguments */
    x = vec_add(x, (__vector double)vec_perm((__vector signed int)x, (__vector signed int)x, (__vector unsigned char)perm));
#else
    x = vec_add(x, vec_perm(x, x, (__vector unsigned char)perm));
#endif
    return vec_extract(x, 0);
}


/****************************************************
 * CONVERSION IMPLEMENTATION HELPER FUNCTIONS       *
 ****************************************************/
static inline void
simdCvtF2DD_ibm_vsx(SimdFloat f0,
                    SimdDouble * d0, SimdDouble * d1)
{
    __vector float fA, fB;
    fA  = vec_mergel(f0, f0); /* 0011 */
    fB  = vec_mergeh(f0, f0); /* 2233 */
    *d0 = gmx_vsx_f2d(fA);    /* 01 */
    *d1 = gmx_vsx_f2d(fB);    /* 23 */
}


static inline SimdFloat
simdCvtDD2F_ibm_vsx(SimdDouble d0, SimdDouble d1)
{
    __vector float fA, fB, fC, fD, fE;
    fA = gmx_vsx_d2f(d0);    /* 0x1x */
    fB = gmx_vsx_d2f(d1);    /* 2x3x */
    fC = vec_mergel(fA, fB); /* 02xx */
    fD = vec_mergeh(fA, fB); /* 13xx */
    fE = vec_mergel(fD, fC); /* 0123 */
    return fE;
}


/* Undefine our temporary work-arounds so they are not used by mistake */
#undef gmx_vsx_f2d
#undef gmx_vsx_d2f

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_DOUBLE_H */
