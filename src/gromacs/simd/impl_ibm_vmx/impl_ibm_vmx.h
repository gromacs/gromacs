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

#include "config.h"

#include <math.h>

#include <altivec.h>

#if defined(__GNUC__) && !defined(__ibmxl__) && !defined(__xlC__)
#    define gmx_vsx_bool __bool
#    undef  bool
#else
#    define gmx_vsx_bool bool
#endif

/* IBM VMX SIMD instruction wrappers. Power6 and later.
 *
 * Please see documentation in gromacs/simd/simd.h for the available
 * defines.
 */

#define GMX_SIMD_V2

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
#define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT
#undef  GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT  /* No need for half-simd, width is 4 */

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
#define gmx_simd_load_f(m)         vec_ld(0, (float *)(m))
#define gmx_simd_load1_f(m)        gmx_simd_load1_f_ibm_vmx(m)
#define gmx_simd_set1_f(x)         gmx_simd_set1_f_ibm_vmx(x)
#define gmx_simd_store_f(m, x)     vec_st(x, 0, (__vector float *)(m))
#undef  gmx_simd_loadu_f
#undef  gmx_simd_storeu_f
#define gmx_simd_setzero_f()       ((__vector float)vec_splat_u32(0))
#define gmx_simd_add_f(a, b)       vec_add(a, b)
#define gmx_simd_sub_f(a, b)       vec_sub(a, b)
#define gmx_simd_mul_f(a, b)       vec_madd(a, b, (__vector float)vec_splat_u32(0))
#define gmx_simd_fmadd_f(a, b, c)  vec_madd(a, b, c)
#define gmx_simd_fmsub_f(a, b, c)  vec_madd(a, b, -c)
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
#define gmx_simd_mul_mask_f(a, b, m)       vec_and(gmx_simd_mul_f(a, b), m)
#define gmx_simd_fmadd_mask_f(a, b, c, m)  vec_and(vec_madd(a, b, c), m)
#ifdef NDEBUG
#    define gmx_simd_rcp_mask_f(a, m)      vec_and(vec_re(a), m)
#    define gmx_simd_rsqrt_mask_f(a, m)    vec_and(vec_rsqrte(a), m)
#else
/* For masked rcp/rsqrt we need to make sure we do not use the masked-out arguments if FP exceptions are enabled */
#    define gmx_simd_rcp_mask_f(a, m)      vec_and(vec_re(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), a, m)), m)
#    define gmx_simd_rsqrt_mask_f(a, m)    vec_and(vec_rsqrte(gmx_simd_blendv_f(gmx_simd_set1_f(1.0f), a, m)), m)
#endif
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
#define gmx_simd_fint32_t          __vector signed int
#define gmx_simd_load_fi(m)        vec_ld(0, (const __vector signed int *)(m))
#define gmx_simd_set1_fi(i)        gmx_simd_set1_fi_ibm_vmx((int)(i))
#define gmx_simd_store_fi(m, x)    vec_st(x, 0, (__vector signed int *)(m))
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
#define gmx_simd_mul_fi(a, b)       ((a)*(b))
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t           __vector gmx_vsx_bool int
#define gmx_simd_cmpeq_f(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmpnz_f(a)        (gmx_simd_fbool_t)vec_cmpgt( (__vector unsigned int)a, vec_splat_u32(0))
#define gmx_simd_cmplt_f(a, b)     vec_cmplt(a, b)
#define gmx_simd_cmple_f(a, b)     vec_cmple(a, b)
#define gmx_simd_and_fb(a, b)      vec_and(a, b)
#define gmx_simd_or_fb(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fb(a)     vec_any_ne(a, (__vector gmx_vsx_bool int)vec_splat_u32(0))
#define gmx_simd_blendzero_f(a, sel)    vec_and(a, (__vector float)sel)
#define gmx_simd_blendnotzero_f(a, sel) vec_andc(a, (__vector float)sel)
#define gmx_simd_blendv_f(a, b, sel)    vec_sel(a, b, sel)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_f_ibm_vmx(a)
/* Boolean & comparison operations on gmx_simd_fint32_t */
#define gmx_simd_fibool_t          __vector gmx_vsx_bool int
#define gmx_simd_cmpeq_fi(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_fi(a, b)     vec_cmplt(a, b)
#define gmx_simd_and_fib(a, b)      vec_and(a, b)
#define gmx_simd_or_fib(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fib(a)          vec_any_ne(a, (__vector gmx_vsx_bool int)vec_splat_u32(0))
#define gmx_simd_blendzero_fi(a, sel)    vec_and(a, (__vector signed int)sel)
#define gmx_simd_blendnotzero_fi(a, sel) vec_andc(a, (__vector signed int)sel)
#define gmx_simd_blendv_fi(a, b, sel)    vec_sel(a, b, sel)
/* Conversions between different booleans */
#define gmx_simd_cvt_fb2fib(x)     (x)
#define gmx_simd_cvt_fib2fb(x)     (x)
/* Higher-level utility functions */
#define gmx_simd_gather_load_transpose_f             gmx_simd_gather_load_transpose_f_ibm_vmx
#define gmx_simd_best_pair_alignment_f               gmx_simd_best_pair_alignment_f_ibm_vmx
#define gmx_simd_gather_loadu_transpose_f            gmx_simd_gather_loadu_transpose_f_ibm_vmx
#define gmx_simd_transpose_scatter_storeu_f          gmx_simd_transpose_scatter_storeu_f_ibm_vmx
#define gmx_simd_transpose_scatter_incru_f           gmx_simd_transpose_scatter_incru_f_ibm_vmx
#define gmx_simd_transpose_scatter_decru_f           gmx_simd_transpose_scatter_decru_f_ibm_vmx
#define gmx_simd_expand_scalars_to_triplets_f        gmx_simd_expand_scalars_to_triplets_f_ibm_vmx
#define gmx_simd_gather_load_bysimdint_transpose_f   gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx
#define gmx_simd_gather_loadu_bysimdint_transpose_f  gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vmx
#define gmx_simd_reduce_incr_4_return_sum_f          gmx_simd_reduce_incr_4_return_sum_f_ibm_vmx

/* Double is not available with VMX SIMD */

/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static gmx_inline gmx_simd_float_t
gmx_simd_set1_f_ibm_vmx(const float x)
{
    __vector float         vx;
    __vector unsigned char perm;

    vx   = vec_lde(0, (float *)&x);
    perm = vec_lvsl(0, (float *)&x);
    vx   = vec_perm(vx, vx, perm);
    return vec_splat(vx, 0);
}

static gmx_inline gmx_simd_float_t
gmx_simd_load1_f_ibm_vmx(const float * m)
{
    __vector float         vx;
    __vector unsigned char perm;

    vx   = vec_lde(0, (float *)m);
    perm = vec_lvsl(0, (float *)m);
    vx   = vec_perm(vx, vx, perm);
    return vec_splat(vx, 0);
}

static gmx_inline gmx_simd_fint32_t
gmx_simd_set1_fi_ibm_vmx(const int x)
{
    __vector signed int    vx;
    __vector unsigned char perm;

    vx   = vec_lde(0, (int *)&x);
    perm = vec_lvsl(0, (int *)&x);
    vx   = vec_perm(vx, vx, perm);
    return vec_splat(vx, 0);
}


static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f_ibm_vmx(gmx_simd_float_t x)
{
    /* Generate 0x7F800000 without memory operations */
    gmx_simd_float_t     expmask = (__vector float)gmx_simd_slli_fi(vec_add(vec_splat_s32(15), vec_sl(vec_splat_s32(15), vec_splat_u32(4))), 23);
    gmx_simd_fint32_t    i127    = vec_sub(vec_sl(vec_splat_s32(1), vec_splat_u32(7)), vec_splat_s32(1));
    gmx_simd_fint32_t    iexp;

    iexp = (__vector signed int)gmx_simd_and_f(x, expmask);
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

/****************************************************
 * Single precision higher-level utility functions  *
 ****************************************************/
#ifdef __cplusplus

/* Internal macro for transposes */
#define GMX_VMX_TRANSPOSE4(v0, v1, v2, v3)               \
    {                                                    \
        __vector float gmx_vmx_t0 = vec_mergeh(v0, v2);  \
        __vector float gmx_vmx_t1 = vec_mergel(v0, v2);  \
        __vector float gmx_vmx_t2 = vec_mergeh(v1, v3);  \
        __vector float gmx_vmx_t3 = vec_mergel(v1, v3);  \
        v0 = vec_mergeh(gmx_vmx_t0, gmx_vmx_t2);         \
        v1 = vec_mergel(gmx_vmx_t0, gmx_vmx_t2);         \
        v2 = vec_mergeh(gmx_vmx_t1, gmx_vmx_t3);         \
        v3 = vec_mergel(gmx_vmx_t1, gmx_vmx_t3);         \
    }


template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vmx(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1,
                                         gmx_simd_float_t    &v2,
                                         gmx_simd_float_t    &v3)
{
    v0  = gmx_simd_load_f( base + align * offset[0] );
    v1  = gmx_simd_load_f( base + align * offset[1] );
    v2  = gmx_simd_load_f( base + align * offset[2] );
    v3  = gmx_simd_load_f( base + align * offset[3] );
    GMX_VMX_TRANSPOSE4(v0, v1, v2, v3);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_transpose_f_ibm_vmx(const float *        base,
                                         const gmx_int32_t    offset[],
                                         gmx_simd_float_t    &v0,
                                         gmx_simd_float_t    &v1)
{
    gmx_simd_float_t       t0, t1, t2, t3, t4, t5, t6, t7;
    __vector unsigned char p0, p1, p2, p3;

    if ((align & 0x3) == 0)
    {
        gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, offset, v0, v1, t2, t3);
    }
    else
    {
        /* This is REALLY slow, since we have no choice but to load individual
         * elements when we cannot guarantee that we can access beyond the end of
         * the memory. Fortunately, 99% of the usage should be the aligned-to-4
         * case above instead.
         */
        t0 = vec_lde(0, base + align * offset[0]);
        t1 = vec_lde(0, base + align * offset[1]);
        t2 = vec_lde(0, base + align * offset[2]);
        t3 = vec_lde(0, base + align * offset[3]);
        p0 = vec_lvsl(0, base + align * offset[0]);
        p1 = vec_lvsl(0, base + align * offset[1]);
        p2 = vec_lvsl(0, base + align * offset[2]);
        p3 = vec_lvsl(0, base + align * offset[3]);
        t0 = vec_perm(t0, t0, p0);
        t1 = vec_perm(t1, t1, p1);
        t2 = vec_perm(t2, t2, p2);
        t3 = vec_perm(t3, t3, p3);
        t0 = vec_mergeh(t0, t2);
        t1 = vec_mergeh(t1, t3);
        v0 = vec_mergeh(t0, t1);

        t4 = vec_lde(0, base + align * offset[0] + 1);
        t5 = vec_lde(0, base + align * offset[1] + 1);
        t6 = vec_lde(0, base + align * offset[2] + 1);
        t7 = vec_lde(0, base + align * offset[3] + 1);
        p0 = vec_lvsl(0, base + align * offset[0] + 1);
        p1 = vec_lvsl(0, base + align * offset[1] + 1);
        p2 = vec_lvsl(0, base + align * offset[2] + 1);
        p3 = vec_lvsl(0, base + align * offset[3] + 1);
        t4 = vec_perm(t4, t4, p0);
        t5 = vec_perm(t5, t5, p1);
        t6 = vec_perm(t6, t6, p2);
        t7 = vec_perm(t7, t7, p3);
        t4 = vec_mergeh(t4, t6);
        t5 = vec_mergeh(t5, t7);
        v1 = vec_mergeh(t4, t5);
    }
}

static const int gmx_simd_best_pair_alignment_f_ibm_vmx = 4;

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_transpose_f_ibm_vmx(const float *        base,
                                          const gmx_int32_t    offset[],
                                          gmx_simd_float_t    &v0,
                                          gmx_simd_float_t    &v1,
                                          gmx_simd_float_t    &v2)
{
    gmx_simd_float_t       t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11;
    __vector unsigned char p0, p1, p2, p3;

    if ((align & 0x3) == 0)
    {
        gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, offset, v0, v1, v2, t3);
    }
    else
    {
        /* This is REALLY slow, since we have no choice but to load individual
         * elements when we cannot guarantee that we can access beyond the end of
         * the memory. Unfortunately this is likely the most common case.
         */
        t0  = vec_lde(0, base + align * offset[0]);
        t1  = vec_lde(0, base + align * offset[1]);
        t2  = vec_lde(0, base + align * offset[2]);
        t3  = vec_lde(0, base + align * offset[3]);
        p0  = vec_lvsl(0, base + align * offset[0]);
        p1  = vec_lvsl(0, base + align * offset[1]);
        p2  = vec_lvsl(0, base + align * offset[2]);
        p3  = vec_lvsl(0, base + align * offset[3]);
        t0  = vec_perm(t0, t0, p0);
        t1  = vec_perm(t1, t1, p1);
        t2  = vec_perm(t2, t2, p2);
        t3  = vec_perm(t3, t3, p3);
        t0  = vec_mergeh(t0, t2);
        t1  = vec_mergeh(t1, t3);
        v0  = vec_mergeh(t0, t1);

        t4  = vec_lde(0, base + align * offset[0] + 1);
        t5  = vec_lde(0, base + align * offset[1] + 1);
        t6  = vec_lde(0, base + align * offset[2] + 1);
        t7  = vec_lde(0, base + align * offset[3] + 1);
        p0  = vec_lvsl(0, base + align * offset[0] + 1);
        p1  = vec_lvsl(0, base + align * offset[1] + 1);
        p2  = vec_lvsl(0, base + align * offset[2] + 1);
        p3  = vec_lvsl(0, base + align * offset[3] + 1);
        t4  = vec_perm(t4, t4, p0);
        t5  = vec_perm(t5, t5, p1);
        t6  = vec_perm(t6, t6, p2);
        t7  = vec_perm(t7, t7, p3);
        t4  = vec_mergeh(t4, t6);
        t5  = vec_mergeh(t5, t7);
        v1  = vec_mergeh(t4, t5);

        t8  = vec_lde(0, base + align * offset[0] + 2);
        t9  = vec_lde(0, base + align * offset[1] + 2);
        t10 = vec_lde(0, base + align * offset[2] + 2);
        t11 = vec_lde(0, base + align * offset[3] + 2);
        p0  = vec_lvsl(0, base + align * offset[0] + 2);
        p1  = vec_lvsl(0, base + align * offset[1] + 2);
        p2  = vec_lvsl(0, base + align * offset[2] + 2);
        p3  = vec_lvsl(0, base + align * offset[3] + 2);
        t8  = vec_perm(t8, t8, p0);
        t9  = vec_perm(t9, t9, p1);
        t10 = vec_perm(t10, t10, p2);
        t11 = vec_perm(t11, t11, p3);
        t8  = vec_mergeh(t8, t10);
        t9  = vec_mergeh(t9, t11);
        v2  = vec_mergeh(t8, t9);
    }
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_storeu_f_ibm_vmx(float *                        base,
                                            const gmx_int32_t              offset[],
                                            gmx_simd_float_t               v0,
                                            gmx_simd_float_t               v1,
                                            gmx_simd_float_t gmx_unused    v2)
{
    __vector float         t3;
    __vector unsigned char p0, p1, p2, p3;

    t3 = (__vector float)vec_splat_u32(0);
    GMX_VMX_TRANSPOSE4(v0, v1, v2, t3);

    p0 = vec_lvsr(0, base + align * offset[0]);
    p1 = vec_lvsr(0, base + align * offset[1]);
    p2 = vec_lvsr(0, base + align * offset[2]);
    p3 = vec_lvsr(0, base + align * offset[3]);

    v0 = vec_perm(v0, v0, p0);
    v1 = vec_perm(v1, v1, p1);
    v2 = vec_perm(v2, v2, p2);
    t3 = vec_perm(t3, t3, p3);

    vec_ste(v0, 0, base + align * offset[0]);
    vec_ste(v0, 4, base + align * offset[0]);
    vec_ste(v0, 8, base + align * offset[0]);
    vec_ste(v1, 0, base + align * offset[1]);
    vec_ste(v1, 4, base + align * offset[1]);
    vec_ste(v1, 8, base + align * offset[1]);
    vec_ste(v2, 0, base + align * offset[2]);
    vec_ste(v2, 4, base + align * offset[2]);
    vec_ste(v2, 8, base + align * offset[2]);
    vec_ste(t3, 0, base + align * offset[3]);
    vec_ste(t3, 4, base + align * offset[3]);
    vec_ste(t3, 8, base + align * offset[3]);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_incru_f_ibm_vmx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    gmx_simd_float_t             t0, t1, t2;

    gmx_simd_gather_loadu_transpose_f_ibm_vmx<align>(base, offset, t0, t1, t2);

    t0 = vec_add(t0, v0);
    t1 = vec_add(t1, v1);
    t2 = vec_add(t2, v2);

    gmx_simd_transpose_scatter_storeu_f_ibm_vmx<align>(base, offset, t0, t1, t2);
}

template <int align>
static gmx_inline void
gmx_simd_transpose_scatter_decru_f_ibm_vmx(float *                       base,
                                           const gmx_int32_t             offset[],
                                           gmx_simd_float_t              v0,
                                           gmx_simd_float_t              v1,
                                           gmx_simd_float_t gmx_unused   v2)
{
    gmx_simd_float_t             t0, t1, t2;

    gmx_simd_gather_loadu_transpose_f_ibm_vmx<align>(base, offset, t0, t1, t2);

    t0 = vec_sub(t0, v0);
    t1 = vec_sub(t1, v1);
    t2 = vec_sub(t2, v2);

    gmx_simd_transpose_scatter_storeu_f_ibm_vmx<align>(base, offset, t0, t1, t2);
}

static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f_ibm_vmx(gmx_simd_float_t    scalar,
                                              gmx_simd_float_t   &triplets0,
                                              gmx_simd_float_t   &triplets1,
                                              gmx_simd_float_t   &triplets2)
{
    __vector unsigned char perm0 = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 7 };
    __vector unsigned char perm1 = { 4, 5, 6, 7, 4, 5, 6, 7, 8, 9, 10, 11, 8, 9, 10, 11 };
    __vector unsigned char perm2 = { 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15, 12, 13, 14, 15 };

    triplets0 = vec_perm(scalar, scalar, perm0);
    triplets1 = vec_perm(scalar, scalar, perm1);
    triplets2 = vec_perm(scalar, scalar, perm2);
}

template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1,
                                                   gmx_simd_float_t   &v2,
                                                   gmx_simd_float_t   &v3)
{
#ifdef __GNUC__
    gmx_int32_t __attribute__ ((aligned (16))) ioffset[4];
    gmx_simd_store_fi(ioffset, offset);
#else
    struct {
        __vector signed int simd; gmx_int32_t i[4];
    } conv;
    gmx_int32_t *ioffset = conv.i;
    gmx_simd_store_fi((int *)&conv.simd, offset);
#endif
    gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, ioffset, v0, v1, v2, v3);
}


template <int align>
static gmx_inline void
gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                   gmx_simd_fint32_t   offset,
                                                   gmx_simd_float_t   &v0,
                                                   gmx_simd_float_t   &v1)
{
#ifdef __GNUC__
    gmx_int32_t __attribute__ ((aligned (16))) ioffset[4];
    gmx_simd_store_fi(ioffset, offset);
#else
    struct {
        __vector signed int simd; gmx_int32_t i[4];
    } conv;
    gmx_int32_t *ioffset = conv.i;
    gmx_simd_store_fi((int *)&conv.simd, offset);
#endif
    gmx_simd_gather_load_transpose_f_ibm_vmx<align>(base, ioffset, v0, v1);
}

template <int align>
static gmx_inline void
gmx_simd_gather_loadu_bysimdint_transpose_f_ibm_vmx(const float *       base,
                                                    gmx_simd_fint32_t   offset,
                                                    gmx_simd_float_t   &v0,
                                                    gmx_simd_float_t   &v1)
{
    gmx_simd_gather_load_bysimdint_transpose_f_ibm_vmx<align>(base, offset, v0, v1);
}


static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f_ibm_vmx(float *           m,
                                            gmx_simd_float_t  v0,
                                            gmx_simd_float_t  v1,
                                            gmx_simd_float_t  v2,
                                            gmx_simd_float_t  v3)
{
    GMX_VMX_TRANSPOSE4(v0, v1, v2, v3);
    v0 = vec_add(v0, v1);
    v2 = vec_add(v2, v3);
    v0 = vec_add(v0, v2);
    v2 = vec_add(v0, gmx_simd_load_f(m));
    gmx_simd_store_f(m, v2);

    return gmx_simd_reduce_f_ibm_vmx(v0);
}
#endif /* __cplusplus */


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
#define gmx_simd4_transpose_f            gmx_simd4_transpose_f_ibm_vmx
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
    gmx_simd4_float_t c = gmx_simd_mul_f(a, b);
    /* Keep only elements 0,1,2 by shifting in zero from right */
    c = vec_sld(c, gmx_simd_setzero_f(), 4);
    /* calculate sum */
    return gmx_simd_reduce_f_ibm_vmx(c);
}

#ifdef __cplusplus
static gmx_inline void gmx_simdcall
gmx_simd4_transpose_f_ibm_vmx(gmx_simd4_float_t &v0, gmx_simd4_float_t &v1,
                              gmx_simd4_float_t &v2, gmx_simd4_float_t &v3)
{
    GMX_VMX_TRANSPOSE4(v0, v1, v2, v3);
}
#endif

/* Function to check whether SIMD operations have resulted in overflow.
 * For now, this is unfortunately a dummy for this architecture.
 */
static int
gmx_simd_check_and_reset_overflow(void)
{
    return 0;
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_VMX_H */
