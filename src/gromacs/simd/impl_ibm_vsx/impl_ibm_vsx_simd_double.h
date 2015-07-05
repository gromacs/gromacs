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
#define gmx_simd_double_t          __vector double
#define gmx_simd_load_d(m)         (*(const gmx_simd_double_t *)(m))
#define gmx_simd_store_d(m, x)      { *(gmx_simd_double_t *)(m) = (x); }
#define gmx_simd_load1_d(m)        vec_splats((double)(*m))
#define gmx_simd_set1_d(x)         vec_splats((double)(x))
#if defined(__ibmxl__) || defined(__xlC__)
#    define gmx_simd_loadu_d(m)    vec_xld2(0, (double *)(m))
#    define gmx_simd_storeu_d(m, x) vec_xstd2(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define gmx_simd_loadu_d       gmx_simd_load_d
#    define gmx_simd_storeu_d      gmx_simd_store_d
#endif
#define gmx_simd_setzero_d()       vec_splats(0.0)
#define gmx_simd_add_d(a, b)       vec_add(a, b)
#define gmx_simd_sub_d(a, b)       vec_sub(a, b)
#define gmx_simd_mul_d(a, b)       vec_mul(a, b)
#define gmx_simd_fmadd_d(a, b, c)  vec_madd(a, b, c)
#define gmx_simd_fmsub_d(a, b, c)  vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_d(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_d(a, b, c) vec_nmadd(a, b, c)
#define gmx_simd_and_d(a, b)       vec_and(a, b)
#define gmx_simd_andnot_d(a, b)    vec_andc(b, a)
#define gmx_simd_or_d(a, b)        vec_or(a, b)
#define gmx_simd_xor_d(a, b)       vec_xor(a, b)
#define gmx_simd_rsqrt_d(a)        vec_rsqrte(a)
#define gmx_simd_rcp_d(a)          vec_re(a)
#define gmx_simd_fabs_d(a)         vec_abs(a)
#define gmx_simd_fneg_d(a)         (-(a))
#define gmx_simd_max_d(a, b)       vec_max(a, b)
#define gmx_simd_min_d(a, b)       vec_min(a, b)
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
/* gcc up to at least version 4.9 does not support vec_round() in double precision. */
#    define gmx_simd_round_d(a)    ({ __vector double res; __asm__ ("xvrdpi %0,%1" : "=ww" (res) : "ww" ((__vector double) (a))); res; })
#else
/* IBM xlC */
#    define gmx_simd_round_d(a)    vec_round(a)
#endif
#define gmx_simd_trunc_d(a)        vec_trunc(a)
#define gmx_simd_fraction_d(x)     vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_d(a) gmx_simd_get_exponent_d_ibm_vsx(a)
#define gmx_simd_get_mantissa_d(a) gmx_simd_get_mantissa_d_ibm_vsx(a)
#define gmx_simd_set_exponent_d(a) gmx_simd_set_exponent_d_ibm_vsx(a)
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t          __vector signed int
#define gmx_simd_load_di(m)        gmx_simd_load_di_ibm_vsx(m)
#define gmx_simd_store_di(m, x)    gmx_simd_store_di_ibm_vsx(m, x)
#define gmx_simd_set1_di(i)        vec_splats((int)(i))
#define gmx_simd_loadu_di          gmx_simd_load_di
#define gmx_simd_storeu_di         gmx_simd_store_di
#define gmx_simd_setzero_di()      vec_splats((int)0)
#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
/* gcc up to at least version 4.9 is missing intrinsics for double precision
 * to integer conversion, use inline asm instead.
 */
#    define gmx_simd_cvtt_d2i(a)   gmx_simd_cvtt_d2i_ibm_vsx(a)
#    define gmx_simd_cvt_i2d(a)    gmx_simd_cvt_i2d_ibm_vsx(a)
#else
/* IBM xlC */
#    define gmx_simd_cvtt_d2i(a)       vec_cts(a, 0)
#    define gmx_simd_cvt_i2d(a)        vec_ctd(a, 0)
#endif
#define gmx_simd_cvt_d2i(a)         gmx_simd_cvtt_d2i(gmx_simd_round_d(a))
#define gmx_simd_extract_di(a, i)   gmx_simd_extract_fi_ibm_vsx(a, (i)*2)
/* Integer logical ops on gmx_simd_dint32_t */
#define gmx_simd_slli_di(a, i)      vec_sl(a, vec_splats((unsigned int)(i)))
#define gmx_simd_srli_di(a, i)      vec_sr(a, vec_splats((unsigned int)(i)))
#define gmx_simd_and_di(a, b)       vec_and(a, b)
#define gmx_simd_andnot_di(a, b)    vec_andc(b, a)
#define gmx_simd_or_di(a, b)        vec_or(a, b)
#define gmx_simd_xor_di(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on gmx_simd_dint32_t */
#define gmx_simd_add_di(a, b)       vec_add(a, b)
#define gmx_simd_sub_di(a, b)       vec_sub(a, b)
#define gmx_simd_mul_di(a, b)       ((a)*(b))
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t           __vector vsx_bool long long
#define gmx_simd_cmpeq_d(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_d(a, b)     vec_cmplt(a, b)
#define gmx_simd_cmple_d(a, b)     vec_cmple(a, b)
#define gmx_simd_and_db(a, b)      (__vector vsx_bool long long)vec_and((__vector signed int)a, (__vector signed int)b)
#define gmx_simd_or_db(a, b)       (__vector vsx_bool long long)vec_or((__vector signed int)a, (__vector signed int)b)
#define gmx_simd_anytrue_db(a)     vec_any_ne((__vector vsx_bool int)a, (__vector vsx_bool int)vec_splats(0))
#define gmx_simd_blendzero_d(a, sel)    vec_and(a, (__vector double)sel)
#define gmx_simd_blendnotzero_d(a, sel) vec_andc(a, (__vector double)sel)
#define gmx_simd_blendv_d(a, b, sel)    vec_sel(a, b, sel)
#define gmx_simd_reduce_d(a)       gmx_simd_reduce_d_ibm_vsx(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_dibool_t          __vector vsx_bool int
#define gmx_simd_cmpeq_di(a, b)    vec_cmpeq(a, b)
#define gmx_simd_cmplt_di(a, b)    vec_cmplt(a, b)
#define gmx_simd_and_dib(a, b)     vec_and(a, b)
#define gmx_simd_or_dib(a, b)      vec_or(a, b)
/* Since we have applied all operations to pairs of elements we can work on all elements here */
#define gmx_simd_anytrue_dib(a)          vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define gmx_simd_blendzero_di(a, sel)    vec_and(a, (__vector signed int)sel)
#define gmx_simd_blendnotzero_di(a, sel) vec_andc(a, (__vector signed int)sel)
#define gmx_simd_blendv_di(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define gmx_simd_cvt_db2dib(x)     (__vector vsx_bool int)(x)
#define gmx_simd_cvt_dib2db(x)     (__vector vsx_bool long long)(x)
/* Float/double conversion */
#define gmx_simd_cvt_f2dd(f, d0, d1)  gmx_simd_cvt_f2dd_ibm_vsx(f, d0, d1)
#define gmx_simd_cvt_dd2f(d0, d1)     gmx_simd_cvt_dd2f_ibm_vsx(d0, d1)

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
static gmx_inline gmx_simd_dint32_t
gmx_simd_load_di_ibm_vsx(const int *m)
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

static gmx_inline void
gmx_simd_store_di_ibm_vsx(int *m, gmx_simd_dint32_t x)
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
static gmx_inline gmx_simd_dint32_t
gmx_simd_cvtt_d2i_ibm_vsx(gmx_simd_double_t x)
{
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
    gmx_simd_dint32_t            ix;

    __asm__ ("xvcvdpsxws %0,%1" : "=ww" (ix) : "ww" ((__vector double) (x)));

    return vec_perm(ix, ix, perm);
}

static gmx_inline gmx_simd_double_t
gmx_simd_cvt_i2d_ibm_vsx(gmx_simd_dint32_t ix)
{
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
    gmx_simd_double_t            x;
    ix = vec_perm(ix, ix, perm);
    __asm__ ("xvcvsxwdp %0,%1" : "=ww" (x) : "ww" ((__vector signed int) (ix)));
    return x;
}
#endif


static gmx_inline gmx_simd_double_t
gmx_simd_get_exponent_d_ibm_vsx(gmx_simd_double_t x)
{
#ifndef __BIG_ENDIAN__
    const __vector unsigned char perm = {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11};
#endif
    gmx_simd_double_t            expmask = (__vector double)vec_splats(0x7FF0000000000000ULL);
    gmx_simd_dint32_t            i1023   = vec_splats(1023);
    gmx_simd_dint32_t            iexp;

    iexp = (__vector signed int)gmx_simd_and_d(x, expmask);
    /* The data is in the upper half of each double (corresponding to elements 1/3).
     * First shift 52-32=20bits, and then permute to swap 0 with 1 and 2 with 3
     * For big endian they are in opposite order, so we avoid the swap.
     */
    iexp = gmx_simd_srli_fi(iexp, 20);
#ifndef __BIG_ENDIAN__
    iexp = vec_perm(iexp, iexp, perm);
#endif
    iexp = vec_sub(iexp, i1023);
    /* Now we have the correct integer in elements 0 & 2. Never mind about elements 1,3 */
    return gmx_simd_cvt_i2d(iexp);
}

static gmx_inline gmx_simd_double_t
gmx_simd_get_mantissa_d_ibm_vsx(gmx_simd_double_t x)
{
    gmx_simd_double_t  expmask = (__vector double)vec_splats(0x7FF0000000000000ULL);

    x = gmx_simd_andnot_d(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return gmx_simd_or_d(x, vec_splats(1.0));
}

static gmx_inline gmx_simd_double_t
gmx_simd_set_exponent_d_ibm_vsx(gmx_simd_double_t x)
{
    gmx_simd_dint32_t            iexp  = gmx_simd_cvt_d2i(x);
    gmx_simd_dint32_t            i1023 = vec_splats(1023);
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

static gmx_inline double
gmx_simd_reduce_d_ibm_vsx(gmx_simd_double_t x)
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
static gmx_inline void
gmx_simd_cvt_f2dd_ibm_vsx(gmx_simd_float_t f0,
                          gmx_simd_double_t * d0, gmx_simd_double_t * d1)
{
    __vector float fA, fB;
    fA  = vec_mergel(f0, f0); /* 0011 */
    fB  = vec_mergeh(f0, f0); /* 2233 */
    *d0 = gmx_vsx_f2d(fA);    /* 01 */
    *d1 = gmx_vsx_f2d(fB);    /* 23 */
}


static gmx_inline gmx_simd_float_t
gmx_simd_cvt_dd2f_ibm_vsx(gmx_simd_double_t d0, gmx_simd_double_t d1)
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
