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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H

#include "config.h"

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"

/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t         vector4double
#ifdef NDEBUG
#    define gmx_simd_load_d(m)    vec_ld(0, (double *)(m))
#    define gmx_simd_store_d(m, a) vec_st(a, 0, (double *)(m))
#else
#    define gmx_simd_load_d(m)    vec_lda(0, (double *)(m))
#    define gmx_simd_store_d(m, a) vec_sta(a, 0, (double *)(m))
#endif
#    define gmx_simd_load1_d(m)   vec_lds(0, (double *)(m))
#define gmx_simd_set1_d(x)        vec_splats(x)
/* No support for unaligned load/store */
#define gmx_simd_setzero_d        gmx_simd_setzero_ibm_qpx
#define gmx_simd_add_d(a, b)       vec_add(a, b)
#define gmx_simd_sub_d(a, b)       vec_sub(a, b)
#define gmx_simd_mul_d(a, b)       vec_mul(a, b)
#define gmx_simd_fmadd_d(a, b, c)   vec_madd(a, b, c)
#define gmx_simd_fmsub_d(a, b, c)   vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_d(a, b, c)  vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_d(a, b, c)  vec_nmadd(a, b, c)
/* gmx_simd_and_d not supported - no bitwise logical ops */
/* gmx_simd_andnot_d not supported - no bitwise logical ops */
/* gmx_simd_or_d not supported - no bitwise logical ops */
/* gmx_simd_xor_d not supported - no bitwise logical ops */
#define gmx_simd_rsqrt_d(a)       vec_rsqrte(a)
#define gmx_simd_rcp_d(a)         vec_re(a)
#define gmx_simd_fabs_d(a)        vec_abs(a)
#define gmx_simd_fneg_d           gmx_simd_fneg_ibm_qpx
#define gmx_simd_max_d(a, b)       vec_sel(b, a, vec_sub(a, b))
#define gmx_simd_min_d(a, b)       vec_sel(b, a, vec_sub(b, a))
/* Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
 * of gmx_simd_round_f(), since vec_round() does not adhere to the FP control
 * word rounding scheme. We rely on float-to-float and float-to-integer
 * rounding being the same for half-way values in a few algorithms.
 */
#define gmx_simd_round_d(a)       vec_cfid(vec_ctid(a))
#define gmx_simd_trunc_d(a)       vec_trunc(a)
#define gmx_simd_fraction_d(x)    vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_d(a) gmx_simd_get_exponent_ibm_qpx(a)
#define gmx_simd_get_mantissa_d(a) gmx_simd_get_mantissa_ibm_qpx(a)
#define gmx_simd_set_exponent_d(a) gmx_simd_set_exponent_ibm_qpx(a)
/* integer datatype corresponding to double: gmx_simd_dint32_t */
#define gmx_simd_dint32_t         vector4double
#ifdef NDEBUG
#    define gmx_simd_load_di(m)   vec_ldia(0, (int *)(m))
#else
#    define gmx_simd_load_di(m)   vec_ldiaa(0, (int *)(m))
#endif
#define gmx_simd_set1_di(i)       gmx_simd_set1_int_ibm_qpx(i)
#define gmx_simd_store_di(m, x)   vec_st(x, 0, (int *)(m))
#define gmx_simd_setzero_di       gmx_simd_setzero_ibm_qpx
#define gmx_simd_cvt_d2i(a)       vec_ctiw(a)
#define gmx_simd_cvtt_d2i(a)      vec_ctiwz(a)
#define gmx_simd_cvt_i2d(a)       vec_cfid(a)
/* Integer simd extract not available */
/* Integer logical ops on gmx_simd_dint32_t not supported */
/* Integer arithmetic ops on gmx_simd_dint32_t not supported */
/* Boolean & comparison operations on gmx_simd_double_t */
#define gmx_simd_dbool_t          vector4double
#define gmx_simd_cmpeq_d(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_d(a, b)     vec_cmplt((a), (b))
#define gmx_simd_cmple_d(a, b)     gmx_simd_or_fb(vec_cmpeq(a, b), vec_cmplt(a, b))
#define gmx_simd_and_db(a, b)      vec_and(a, b)
#define gmx_simd_or_db(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_db(a)    gmx_simd_anytrue_bool_ibm_qpx(a)
#define gmx_simd_blendzero_d(a, sel) vec_sel(vec_splats(0.0), a, sel)
#define gmx_simd_blendnotzero_d(a, sel) vec_sel(a, vec_splats(0.0), sel)
#define gmx_simd_blendv_d(a, b, sel)  vec_sel(a, b, sel)
#define gmx_simd_reduce_d(a)      gmx_simd_reduce_ibm_qpx(a)

/* Boolean & comparison operations on gmx_simd_dint32_t not supported */
/* Conversions between different booleans not supported */


/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static __attribute__((always_inline)) vector4double
gmx_simd_fneg_ibm_qpx(vector4double a)
{
    return vec_neg(a);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_setzero_ibm_qpx(void)
{
    return vec_splats(0.0);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_get_exponent_ibm_qpx(vector4double x)
{
    const gmx_int64_t    expmask   = 0x7ff0000000000000LL;
    const gmx_int64_t    expbase   = 1023;
    gmx_int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory */
    vec_st(x, 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = ((idata[0] & expmask) >> 52) - expbase;
    idata[1] = ((idata[1] & expmask) >> 52) - expbase;
    idata[2] = ((idata[2] & expmask) >> 52) - expbase;
    idata[3] = ((idata[3] & expmask) >> 52) - expbase;
    /* Reload from memory */
    return vec_cfid(vec_ld(0, idata));
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_get_mantissa_ibm_qpx(vector4double x)
{
    const gmx_int64_t    exp_and_sign_mask = 0xfff0000000000000LL;
    const gmx_int64_t    ione              = 0x3ff0000000000000LL;
    gmx_int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory */
    vec_st(x, 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = (idata[0] & (~exp_and_sign_mask)) | ione;
    idata[1] = (idata[1] & (~exp_and_sign_mask)) | ione;
    idata[2] = (idata[2] & (~exp_and_sign_mask)) | ione;
    idata[3] = (idata[3] & (~exp_and_sign_mask)) | ione;
    /* Reload from memory */
    return vec_ld(0, idata);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_set_exponent_ibm_qpx(vector4double x)
{
    const gmx_int64_t    expbase = 1023;
    gmx_int64_t          idata[4] __attribute__((aligned(32)));

    /* Store to memory for shifts. It is REALLY critical that we use the same
     * rounding mode as for gmx_simd_round_r() here. In particular, for QPX
     * this means we implement gmx_simd_round_r(a) as vec_cfid(vec_ctid(a)),
     * since vec_round() uses a different rounding scheme.
     */
    vec_st(vec_ctid(x), 0, idata);
    /* Perform integer arithmetics in general registers. */
    idata[0] = (idata[0] + expbase) << 52;
    idata[1] = (idata[1] + expbase) << 52;
    idata[2] = (idata[2] + expbase) << 52;
    idata[3] = (idata[3] + expbase) << 52;
    /* Reload from memory */
    return vec_ld(0, idata);
}

static __attribute__((always_inline)) double gmx_simdcall
gmx_simd_reduce_ibm_qpx(vector4double x)
{
    vector4double y = vec_sldw(x, x, 2);
    vector4double z;

    y = vec_add(y, x);
    z = vec_sldw(y, y, 1);
    y = vec_add(y, z);
    return vec_extract(y, 0);
}

static __attribute__((always_inline)) vector4double gmx_simdcall
gmx_simd_set1_int_ibm_qpx(int i)
{
    int idata[4] __attribute__((aligned(32)));

    idata[0] = i;

    /* Reload from memory */
    return vec_splat(vec_ldia(0, idata), 0);
}

/* This works in both single and double */
static __attribute__((always_inline)) int gmx_simdcall
gmx_simd_anytrue_bool_ibm_qpx(vector4double a)
{
    vector4double b = vec_sldw(a, a, 2);

    a = vec_or(a, b);
    b = vec_sldw(a, a, 1);
    a = vec_or(a, b);
    return (vec_extract(a, 0) > 0);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_DOUBLE_H */
