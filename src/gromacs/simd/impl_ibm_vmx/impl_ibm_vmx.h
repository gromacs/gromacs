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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_VMX_H
#define GMX_SIMD_IMPLEMENTATION_IBM_VMX_H

#include <math.h>

#include <altivec.h>

/* Make sure we do not screw up c++ - undefine vector/bool, and rely on __vector and __bool */
#undef vector
#undef bool

/* IBM VMX SIMD instruction wrappers. Power6 and later.
 *
 * Please see documentation in gromacs/simd/simd.h for the available
 * defines.
 */
/* Capability definitions for IBM VMX */
#define GMX_SIMD_HAVE_FLOAT
#undef  GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_HARDWARE
#undef  GMX_SIMD_HAVE_LOADU
#undef  GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_LOGICAL
/* VMX only provides fmadd/fnmadd (our definitions), but not fmsub/fnmsub.
 * However, fnmadd is what we need for 1/sqrt(x).
 */
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#undef  GMX_SIMD_HAVE_FINT32_EXTRACT
#define GMX_SIMD_HAVE_FINT32_LOGICAL
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS
#undef  GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
#undef  GMX_SIMD_HAVE_DINT32_LOGICAL
#undef  GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD4_HAVE_FLOAT
#undef  GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         4
#undef  GMX_SIMD_DOUBLE_WIDTH
#define GMX_SIMD_FINT32_WIDTH        4
#undef  GMX_SIMD_DINT32_WIDTH
#define GMX_SIMD_RSQRT_BITS         14
#define GMX_SIMD_RCP_BITS           14

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t           __vector float
#define gmx_simd_load_f(m)         vec_ld(0, (const __vector float *)(m))
#define gmx_simd_load1_f(m)        gmx_simd_load1_f_ibm_vmx(m)
#define gmx_simd_set1_f(x)         gmx_simd_set1_f_ibm_vmx(x)
#define gmx_simd_store_f(m, x)     vec_st(x, 0, (__vector float *)(m))
#undef  gmx_simd_loadu_f
#undef  gmx_simd_storeu_f
#define gmx_simd_setzero_f()       ((__vector float)vec_splat_u32(0))
#define gmx_simd_add_f(a, b)       vec_add(a, b)
#define gmx_simd_sub_f(a, b)       vec_sub(a, b)
#define gmx_simd_mul_f(a, b)       vec_mul(a, b)
#define gmx_simd_fmadd_f(a, b, c)  vec_madd(a, b, c)
#define gmx_simd_fmsub_f(a, b, c)  vec_sub(vec_mul(a, b), c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_f(a, b, c) vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_f(a, b, c) vec_sub(gmx_simd_setzero_f(), vec_madd(a, b, c))
#define gmx_simd_and_f(a, b)       vec_and(a, b)
#define gmx_simd_andnot_f(a, b)    vec_andc(b, a)
#define gmx_simd_or_f(a, b)        vec_or(a, b)
#define gmx_simd_xor_f(a, b)       vec_xor(a, b)
#define gmx_simd_rsqrt_f(a)        vec_rsqrte(a)
#define gmx_simd_rcp_f(a)          vec_re(a)
#define gmx_simd_fabs_f(a)         vec_abs(a)
#define gmx_simd_fneg_f(a)         vec_xor(a, (__vector float)vec_sl(vec_splat_u32(-1), vec_splat_u32(-1)))
#define gmx_simd_max_f(a, b)       vec_max(a, b)
#define gmx_simd_min_f(a, b)       vec_min(a, b)
#define gmx_simd_round_f(a)        vec_round(a)
#define gmx_simd_trunc_f(a)        vec_trunc(a)
#define gmx_simd_fraction_f(x)     vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_f(a) gmx_simd_get_exponent_f_ibm_vmx(a)
#define gmx_simd_get_mantissa_f(a) gmx_simd_get_mantissa_f_ibm_vmx(a)
#define gmx_simd_set_exponent_f(a) gmx_simd_set_exponent_f_ibm_vmx(a)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t          __vector int
#define gmx_simd_load_fi(m)        vec_ld(0, (const __vector int *)(m))
#define gmx_simd_set1_fi(i)        gmx_simd_set1_fi_ibm_vmx((int)(i))
#define gmx_simd_store_fi(m, x)    vec_st(x, 0, (__vector int *)(m))
#undef  gmx_simd_loadu_fi
#undef  gmx_simd_storeu_fi
#define gmx_simd_setzero_fi()      vec_splat_s32(0)
#define gmx_simd_cvt_f2i(a)        vec_cts(vec_round(a), 0)
#define gmx_simd_cvtt_f2i(a)       vec_cts(a, 0)
#define gmx_simd_cvt_i2f(a)        vec_ctf(a, 0)
#undef  gmx_simd_extract_fi
/* Integer logical ops on gmx_simd_fint32_t */
/* The shift constant magic requires an explanation:
 * VMX only allows literals up to 15 to be created directly with vec_splat_u32,
 * and we need to be able to shift up to 31 bits. The code on the right hand
 * side splits the constant in three parts with values in the range 0..15.
 * Since the argument has to be a constant (but our and VMX requirement), these
 * constants will be evaluated at compile-time, and if one or two parts evaluate
 * to zero they will be removed with -O2 or higher optimization (checked).
 */
#define gmx_simd_slli_fi(a, i)      vec_sl(a, vec_add(vec_add(vec_splat_u32( (((i&0xF)+(i/16))&0xF)+i/31 ), vec_splat_u32( (i/16)*15 )), vec_splat_u32( (i/31)*15 )))
#define gmx_simd_srli_fi(a, i)      vec_sr(a, vec_add(vec_add(vec_splat_u32( (((i&0xF)+(i/16))&0xF)+i/31 ), vec_splat_u32( (i/16)*15 )), vec_splat_u32( (i/31)*15 )))
#define gmx_simd_and_fi(a, b)       vec_and(a, b)
#define gmx_simd_andnot_fi(a, b)   vec_andc(b, a)
#define gmx_simd_or_fi(a, b)        vec_or(a, b)
#define gmx_simd_xor_fi(a, b)       vec_xor(a, b)
/* Integer arithmetic ops on gmx_simd_fint32_t */
#define gmx_simd_add_fi(a, b)       vec_add(a, b)
#define gmx_simd_sub_fi(a, b)       vec_sub(a, b)
#define gmx_simd_mul_fi(a, b)       vec_mule((__vector short)a, (__vector short)b)
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __vector __bool int
#define gmx_simd_cmpeq_f(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_f(a, b)     vec_cmplt(a, b)
#define gmx_simd_cmple_f(a, b)     vec_cmple(a, b)
#define gmx_simd_and_fb(a, b)      vec_and(a, b)
#define gmx_simd_or_fb(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fb(a)     vec_any_ne(a, (__vector __bool int)vec_splat_u32(0))
#define gmx_simd_blendzero_f(a, sel)    vec_and(a, (__vector float)sel)
#define gmx_simd_blendnotzero_f(a, sel) vec_andc(a, (__vector float)sel)
#define gmx_simd_blendv_f(a, b, sel)    vec_sel(a, b, sel)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_ibm_vmx(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __vector __bool int
#define gmx_simd_cmpeq_fi(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_fi(a, b)     vec_cmplt(a, b)
#define gmx_simd_and_fib(a, b)      vec_and(a, b)
#define gmx_simd_or_fib(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fib(a)          vec_any_ne(a, (__vector __bool int)vec_splat_u32(0))
#define gmx_simd_blendzero_fi(a, sel)    vec_and(a, (__vector int)sel)
#define gmx_simd_blendnotzero_fi(a, sel) vec_andc(a, (__vector int)sel)
#define gmx_simd_blendv_fi(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)

/* Double is not available with VMX SIMD */

/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_set1_f_ibm_vmx(const float x)
{
    /* In the old days when PPC was all big endian we could
     * use the vec_lvsl() instruction to permute bytes based on
     * a load adress. However, at least with gcc-4.8.2 the bytes
     * end up reversed on Power8 running little endian (Linux).
     * Since this is not a critical instruction we work around
     * it by first putting the data in an aligned position before
     * loading, so we can avoid vec_lvsl() entirely. We can
     * do this slightly faster on GCC with alignment attributes.
     */
    __vector float vx;
#ifdef __GNUC__
    float alignedx __attribute ((aligned (16)));
    alignedx = x;
    vx       = vec_lde(0, &alignedx);
#else
    struct {
        vector float vx; float x[4];
    } conv;
    conv.x[0] = x;
    vx        = vec_lde(0, conv.x);
#endif
    return vec_splat(vx, 0);
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_f_ibm_vmx(const float * m)
{
    return gmx_simd_set1_f_ibm_vmx(*m);
}

static gmx_inline gmx_simd_fint32_t
gmx_simd_set1_fi_ibm_vmx(const int x)
{
    /* See comment in gmx_simd_set1_f_ibm_vmx why we
     * cannot use vec_lvsl().
     */
    __vector int vx;
#ifdef __GNUC__
    int alignedx __attribute ((aligned (16)));
    alignedx = x;
    vx       = vec_lde(0, &alignedx);
#else
    struct {
        vector int vx; int x[4];
    } conv;
    conv.x[0] = x;
    vx        = vec_lde(0, conv.x);
#endif
    return vec_splat(vx, 0);
}


static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f_ibm_vmx(gmx_simd_float_t x)
{
    /* Generate 0x7F800000 without memory operations */
    gmx_simd_float_t     expmask = (__vector float)gmx_simd_slli_fi(vec_add(vec_splat_s32(15), vec_sl(vec_splat_s32(15), vec_splat_u32(4))), 23);
    gmx_simd_fint32_t    i127    = vec_sub(vec_sl(vec_splat_s32(1), vec_splat_u32(7)), vec_splat_s32(1));
    gmx_simd_fint32_t    iexp;

    iexp = (__vector int)gmx_simd_and_f(x, expmask);
    iexp = vec_sub(gmx_simd_srli_fi(iexp, 23), i127);
    return vec_ctf(iexp, 0);
}

static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f_ibm_vmx(gmx_simd_float_t x)
{
    gmx_simd_float_t     expmask = (__vector float)gmx_simd_slli_fi(vec_add(vec_splat_s32(15), vec_sl(vec_splat_s32(15), vec_splat_u32(4))), 23);

    /* Get mantissa. By taking the absolute value (to get rid of the sign bit) we can
     * use the same mask as for gmx_simd_get_exponent_f() (but complement it). Since
     * these two routines are typically called together, this will save a few operations.
     */
    x = gmx_simd_andnot_f(expmask, vec_abs(x));
    /* Reset zero (but correctly biased) exponent */
    return gmx_simd_or_f(x, vec_ctf(vec_splat_s32(1), 0));
}

static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f_ibm_vmx(gmx_simd_float_t x)
{
    gmx_simd_fint32_t  iexp = gmx_simd_cvt_f2i(x);
    gmx_simd_fint32_t  i127 = vec_sub(vec_sl(vec_splat_s32(1), vec_splat_u32(7)), vec_splat_s32(1));

    iexp = gmx_simd_slli_fi(vec_add(iexp, i127), 23);
    return (__vector float)iexp;
}

static gmx_inline float
gmx_simd_reduce_f_ibm_vmx(gmx_simd_float_t x)
{
    float res;
    x = vec_add(x, vec_sld(x, x, 8));
    x = vec_add(x, vec_sld(x, x, 4));
    vec_ste(x, 0, &res);
    return res;
}



/* SINGLE */
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
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_ibm_vmx
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
gmx_simd4_dotproduct3_f_ibm_vmx(gmx_simd4_float_t a, gmx_simd4_float_t b)
{
    gmx_simd4_float_t c = vec_mul(a, b);
    /* Keep only elements 0,1,2 by shifting in zero from right */
    c = vec_sld(c, gmx_simd_setzero_f(), 4);
    /* calculate sum */
    return gmx_simd_reduce_f_ibm_vmx(c);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VMX_H */
