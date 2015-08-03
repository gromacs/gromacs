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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_H

#include <math.h>

#include <altivec.h>

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

/* Since we've had to use quite a few compiler and endian defines in this code
 * it might not 'just work' when e.g. clang provides VSX support. If you are
 * reading this because you saw the error below for a new compiler, try removing
 * the checks, but then make sure you run 'make test'.
 */
#if !(defined __GNUC__) && !(defined __ibmxl__) && !(defined __xlC__)
#    error VSX acceleration is very compiler-dependent, and only tested for gcc & xlc.
#endif

/* Capability definitions for IBM VSX */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#define GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#define GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS

/* With GCC, only version 4.9 or later supports all parts of double precision VSX.
 * We check explicitly for xlc, since that compiler appears to like pretending it is gcc,
 * but there double precision seems to work fine.
 */
#if defined(__ibmxl__) || defined(__xlC__) || !(defined(__GNUC__) && ((__GNUC__ < 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ < 9))))
#    define GMX_SIMD_HAVE_DOUBLE
#    define GMX_SIMD_HAVE_DINT32
#    define GMX_SIMD_HAVE_DINT32_EXTRACT
#    define GMX_SIMD_HAVE_DINT32_LOGICAL
#    define GMX_SIMD_HAVE_DINT32_ARITHMETICS
#endif

#define GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         4
#define GMX_SIMD_DOUBLE_WIDTH        2
#define GMX_SIMD_FINT32_WIDTH        4
#define GMX_SIMD_DINT32_WIDTH        2
#define GMX_SIMD_RSQRT_BITS         14
#define GMX_SIMD_RCP_BITS           14

/* The VSX load & store operations are a slight mess. The interface is different
 * for xlc version 12, xlc version 13, and gcc. Long-term IBM recommends
 * simply using pointer dereferencing both for aligned and unaligned loads.
 * That's nice, but unfortunately xlc still bugs out when the pointer is
 * not aligned. Sticking to vec_xl/vec_xst isn't a solution either, since
 * that appears to be buggy for some _aligned_ loads :-)
 *
 * For now, we use pointer dereferencing for all aligned load/stores, and
 * for unaligned ones with gcc. On xlc we use vec_xlw4/vec_xstw4 for
 * unaligned memory operations. The latest docs recommend using the overloaded
 * vec_xl/vec_xst, but that is not supported on xlc version 12. We'll
 * revisit things once xlc is a bit more stable - for now you probably want
 * to stick to gcc...
 */

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __vector float
#define gmx_simd_load_f(m)         (*(const gmx_simd_float_t *)(m))
#define gmx_simd_store_f(m, x)      { *(gmx_simd_float_t *)(m) = (x); }
#define gmx_simd_load1_f(m)        vec_splats((float)(*m))
#define gmx_simd_set1_f(x)         vec_splats((float)(x))
#if defined(__ibmxl__) || defined(__xlC__)
#    define gmx_simd_loadu_f(m)    vec_xlw4(0, (float *)(m))
#    define gmx_simd_storeu_f(m, x) vec_xstw4(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define gmx_simd_loadu_f       gmx_simd_load_f
#    define gmx_simd_storeu_f      gmx_simd_store_f
#endif
#define gmx_simd_setzero_f()       vec_splats(0.0f)
#define gmx_simd_add_f(a, b)       vec_add(a, b)
#define gmx_simd_sub_f(a, b)       vec_sub(a, b)
#define gmx_simd_mul_f(a, b)       vec_mul(a, b)
#define gmx_simd_fmadd_f(a, b, c)  vec_madd(a, b, c)
#define gmx_simd_fmsub_f(a, b, c)  vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_f(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_f(a, b, c) vec_nmadd(a, b, c)
#define gmx_simd_and_f(a, b)       vec_and(a, b)
#define gmx_simd_andnot_f(a, b)    vec_andc(b, a)
#define gmx_simd_or_f(a, b)        vec_or(a, b)
#define gmx_simd_xor_f(a, b)       vec_xor(a, b)
#define gmx_simd_rsqrt_f(a)        vec_rsqrte(a)
#define gmx_simd_rcp_f(a)          vec_re(a)
#define gmx_simd_fabs_f(a)         vec_abs(a)
#define gmx_simd_fneg_f(a)         (-(a))
#define gmx_simd_max_f(a, b)       vec_max(a, b)
#define gmx_simd_min_f(a, b)       vec_min(a, b)
#define gmx_simd_round_f(a)        vec_round(a)
#define gmx_simd_trunc_f(a)        vec_trunc(a)
#define gmx_simd_fraction_f(x)     vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_f(a) gmx_simd_get_exponent_f_ibm_vsx(a)
#define gmx_simd_get_mantissa_f(a) gmx_simd_get_mantissa_f_ibm_vsx(a)
#define gmx_simd_set_exponent_f(a) gmx_simd_set_exponent_f_ibm_vsx(a)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __vector signed int
#define gmx_simd_load_fi(m)        (*(const gmx_simd_fint32_t *)(m))
#define gmx_simd_store_fi(m, x)     { *(gmx_simd_fint32_t *)(m) = (x); }
#define gmx_simd_set1_fi(i)        vec_splats((int)(i))
#if defined(__ibmxl__) || defined(__xlC__)
#    define gmx_simd_loadu_fi(m)    vec_xlw4(0, (int *)(m))
#    define gmx_simd_storeu_fi(m, x) vec_xstw4(x, 0, (m))
#else
/* GCC can handle unaligned load/store as pointer dereference */
#    define gmx_simd_loadu_fi      gmx_simd_load_fi
#    define gmx_simd_storeu_fi     gmx_simd_store_fi
#endif
#define gmx_simd_setzero_fi()      vec_splats((int)0)
#define gmx_simd_cvt_f2i(a)        vec_cts(vec_round(a), 0)
#define gmx_simd_cvtt_f2i(a)       vec_cts(a, 0)
#define gmx_simd_cvt_i2f(a)        vec_ctf(a, 0)
#define gmx_simd_extract_fi(a, i)   gmx_simd_extract_fi_ibm_vsx(a, i)
/* Integer logical ops on gmx_simd_fint32_t */
#define gmx_simd_slli_fi(a, i)      vec_sl(a, vec_splats((unsigned int)(i)))
#define gmx_simd_srli_fi(a, i)      vec_sr(a, vec_splats((unsigned int)(i)))
#define gmx_simd_and_fi(a, b)       vec_and(a, b)
#define gmx_simd_andnot_fi(a, b)    vec_andc(b, a)
#define gmx_simd_or_fi(a, b)        vec_or(a, b)
#define gmx_simd_xor_fi(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi(a, b)       vec_add(a, b)
#define gmx_simd_sub_fi(a, b)       vec_sub(a, b)
#define gmx_simd_mul_fi(a, b)      ((a)*(b))
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __vector vsx_bool int
#define gmx_simd_cmpeq_f(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_f(a, b)     vec_cmplt(a, b)
#define gmx_simd_cmple_f(a, b)     vec_cmple(a, b)
#define gmx_simd_and_fb(a, b)      vec_and(a, b)
#define gmx_simd_or_fb(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fb(a)     vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define gmx_simd_blendzero_f(a, sel)    vec_and(a, (__vector float)sel)
#define gmx_simd_blendnotzero_f(a, sel) vec_andc(a, (__vector float)sel)
#define gmx_simd_blendv_f(a, b, sel)    vec_sel(a, b, sel)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_ibm_vsx(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __vector vsx_bool int
#define gmx_simd_cmpeq_fi(a, b)    vec_cmpeq(a, b)
#define gmx_simd_cmplt_fi(a, b)    vec_cmplt(a, b)
#define gmx_simd_and_fib(a, b)     vec_and(a, b)
#define gmx_simd_or_fib(a, b)      vec_or(a, b)
#define gmx_simd_anytrue_fib(a)          vec_any_ne(a, (__vector vsx_bool int)vec_splats(0))
#define gmx_simd_blendzero_fi(a, sel)    vec_and(a, (__vector signed int)sel)
#define gmx_simd_blendnotzero_fi(a, sel) vec_andc(a, (__vector signed int)sel)
#define gmx_simd_blendv_fi(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

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
 * SINGLE PREC. IMPLEMENTATION HELPER FUNCTIONS     *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f_ibm_vsx(gmx_simd_float_t x)
{
    /* Generate 0x7F800000 without memory operations */
    gmx_simd_float_t     expmask = (__vector float)vec_splats(0x7F800000);
    gmx_simd_fint32_t    i127    = vec_splats(127);
    gmx_simd_fint32_t    iexp;

    iexp = (__vector signed int)gmx_simd_and_f(x, expmask);
    iexp = vec_sub(gmx_simd_srli_fi(iexp, 23), i127);
    return vec_ctf(iexp, 0);
}

static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f_ibm_vsx(gmx_simd_float_t x)
{
    gmx_simd_float_t     expmask = (__vector float)vec_splats(0x7F800000);

    x = gmx_simd_andnot_f(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return gmx_simd_or_f(x, vec_splats(1.0f));
}

static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f_ibm_vsx(gmx_simd_float_t x)
{
    gmx_simd_fint32_t  iexp = gmx_simd_cvt_f2i(x);
    gmx_simd_fint32_t  i127 = vec_splats(127);

    iexp = gmx_simd_slli_fi(vec_add(iexp, i127), 23);
    return (__vector float)iexp;
}

static gmx_inline float
gmx_simd_reduce_f_ibm_vsx(gmx_simd_float_t x)
{
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };

    x = vec_add(x, vec_perm(x, x, perm1));
    x = vec_add(x, vec_perm(x, x, perm2));
    return vec_extract(x, 0);
}

/* gcc-4.9 does not detect that vec_extract() uses its argument */
static gmx_inline int
gmx_simd_extract_fi_ibm_vsx(gmx_simd_fint32_t gmx_unused x, unsigned int i)
{
    return vec_extract(x, i);
}

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

/* Single precision VSX is 4 elements wide, use for SIMD4 */
#define gmx_simd4_float_t                gmx_simd_float_t
#define gmx_simd4_load_f                 gmx_simd_load_f
#define gmx_simd4_load1_f                gmx_simd_load1_f
#define gmx_simd4_set1_f                 gmx_simd_set1_f
#define gmx_simd4_store_f                gmx_simd_store_f
#define gmx_simd4_loadu_f                gmx_simd_loadu_f
#define gmx_simd4_storeu_f               gmx_simd_storeu_f
#define gmx_simd4_setzero_f              gmx_simd_setzero_f
#define gmx_simd4_add_f                  gmx_simd_add_f
#define gmx_simd4_sub_f                  gmx_simd_sub_f
#define gmx_simd4_mul_f                  gmx_simd_mul_f
#define gmx_simd4_fmadd_f                gmx_simd_fmadd_f
#define gmx_simd4_fmsub_f                gmx_simd_fmsub_f
#define gmx_simd4_fnmadd_f               gmx_simd_fnmadd_f
#define gmx_simd4_fnmsub_f               gmx_simd_fnmsub_f
#define gmx_simd4_and_f                  gmx_simd_and_f
#define gmx_simd4_andnot_f               gmx_simd_andnot_f
#define gmx_simd4_or_f                   gmx_simd_or_f
#define gmx_simd4_xor_f                  gmx_simd_xor_f
#define gmx_simd4_rsqrt_f                gmx_simd_rsqrt_f
#define gmx_simd4_rcp_f                  gmx_simd_rcp_f
#define gmx_simd4_fabs_f                 gmx_simd_fabs_f
#define gmx_simd4_fneg_f                 gmx_simd_fneg_f
#define gmx_simd4_max_f                  gmx_simd_max_f
#define gmx_simd4_min_f                  gmx_simd_min_f
#define gmx_simd4_round_f                gmx_simd_round_f
#define gmx_simd4_trunc_f                gmx_simd_trunc_f
#define gmx_simd4_fraction_f             gmx_simd_fraction_f
#define gmx_simd4_get_exponent_f         gmx_simd_get_exponent_f
#define gmx_simd4_get_mantissa_f         gmx_simd_get_mantissa_f
#define gmx_simd4_set_exponent_f         gmx_simd_set_exponent_f
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_ibm_vsx
#define gmx_simd4_fint32_t               gmx_simd_fint32_t
#define gmx_simd4_load_fi                gmx_simd_load_fi
#define gmx_simd4_load1_fi               gmx_simd_load1_fi
#define gmx_simd4_set1_fi                gmx_simd_set1_fi
#define gmx_simd4_store_fi               gmx_simd_store_fi
#define gmx_simd4_loadu_fi               gmx_simd_loadu_fi
#define gmx_simd4_storeu_fi              gmx_simd_storeu_fi
#define gmx_simd4_setzero_fi             gmx_simd_setzero_fi
#define gmx_simd4_cvt_f2i                gmx_simd_cvt_f2i
#define gmx_simd4_cvtt_f2i               gmx_simd_cvtt_f2i
#define gmx_simd4_cvt_i2f                gmx_simd_cvt_i2f
#define gmx_simd4_fbool_t                gmx_simd_fbool_t
#define gmx_simd4_cmpeq_f                gmx_simd_cmpeq_f
#define gmx_simd4_cmplt_f                gmx_simd_cmplt_f
#define gmx_simd4_cmple_f                gmx_simd_cmple_f
#define gmx_simd4_and_fb                 gmx_simd_and_fb
#define gmx_simd4_or_fb                  gmx_simd_or_fb
#define gmx_simd4_anytrue_fb             gmx_simd_anytrue_fb
#define gmx_simd4_blendzero_f            gmx_simd_blendzero_f
#define gmx_simd4_blendnotzero_f         gmx_simd_blendnotzero_f
#define gmx_simd4_blendv_f               gmx_simd_blendv_f
#define gmx_simd4_reduce_f               gmx_simd_reduce_f

static gmx_inline float
gmx_simd4_dotproduct3_f_ibm_vsx(gmx_simd4_float_t a, gmx_simd4_float_t b)
{
    gmx_simd4_float_t            c     = gmx_simd_mul_f(a, b);
    const __vector unsigned char perm1 = { 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7 };
    const __vector unsigned char perm2 = { 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3 };
    __vector float               sum;
    sum = vec_add(c, vec_perm(c, c, (__vector unsigned char)perm1));
    sum = vec_add(sum, vec_perm(c, c, (__vector unsigned char)perm2));
    return vec_extract(sum, 0);
}

/* Undefine our temporary work-arounds so they are not used by mistake */
#undef gmx_vsx_f2d
#undef gmx_vsx_d2f

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_H */
