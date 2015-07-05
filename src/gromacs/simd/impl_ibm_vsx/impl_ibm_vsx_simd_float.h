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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H

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

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VSX_SIMD_FLOAT_H */
