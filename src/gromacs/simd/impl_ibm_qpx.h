/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_H

#ifndef GMX_SIMD_SIMD_H
# error "Never include architecture-specific simd headers directly; use simd.h."
#endif

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

/* IBM QPX SIMD instruction wrappers
 *
 * Please see documentation in gromacs/simd/simd.h for the available
 * defines.
 */
/* Capability definitions for IBM QPX */
#define GMX_SIMD_HAVE_FLOAT
#define GMX_SIMD_HAVE_DOUBLE
#define GMX_SIMD_HAVE_HARDWARE
#define GMX_SIMD_HAVE_LOADU
#undef  GMX_SIMD_HAVE_LOADU
#undef  GMX_SIMD_HAVE_STOREU
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FRACTION
#define GMX_SIMD_HAVE_FINT32
#undef  GMX_SIMD_HAVE_FINT32_EXTRACT
#undef  GMX_SIMD_HAVE_FINT32_LOGICAL
#undef  GMX_SIMD_HAVE_FINT32_ARITHMETICS
#define GMX_SIMD_HAVE_DINT32
#undef  GMX_SIMD_HAVE_DINT32_EXTRACT
#undef  GMX_SIMD_HAVE_DINT32_LOGICAL
#undef  GMX_SIMD_HAVE_DINT32_ARITHMETICS
#define GMX_SIMD4_HAVE_FLOAT
#define GMX_SIMD4_HAVE_DOUBLE

/* Implementation details */
#define GMX_SIMD_FLOAT_WIDTH         4
#define GMX_SIMD_DOUBLE_WIDTH        4
#define GMX_SIMD_FINT32_WIDTH        4
#define GMX_SIMD_DINT32_WIDTH        4
#define GMX_SIMD_RSQRT_BITS          5
#define GMX_SIMD_RCP_BITS            8

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t          vector4double
#ifdef NDEBUG
#    define gmx_simd_load_f(m)    vec_ld(0, (float *)(m))
#    define gmx_simd_load1_f(m)   vec_lds(0, (float *)(m))
#    define gmx_simd_store_f(m, a) vec_st(a, 0, (float *)(m))
#else
#    define gmx_simd_load_f(m)    vec_lda(0, (float *)(m))
#    define gmx_simd_load1_f(m)   vec_ldsa(0, (float *)(m))
#    define gmx_simd_store_f(m, a) vec_sta(a, 0, (float *)(m))
#endif
#define gmx_simd_set1_f(x)        vec_splats(x)
/* No support for unaligned load/store */
#define gmx_simd_setzero_f        gmx_simd_setzero_ibm_qpx
#define gmx_simd_add_f(a, b)       vec_add(a, b)
#define gmx_simd_sub_f(a, b)       vec_sub(a, b)
#define gmx_simd_mul_f(a, b)       vec_mul(a, b)
#define gmx_simd_fmadd_f(a, b, c)   vec_madd(a, b, c)
#define gmx_simd_fmsub_f(a, b, c)   vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_f(a, b, c)  vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_f(a, b, c)  vec_nmadd(a, b, c)
#define gmx_simd_and_f(a, b)       vec_and(a, b)
#define gmx_simd_andnot_f(a, b)    vec_andc(b, a)
#define gmx_simd_or_f(a, b)        vec_or(a, b)
#define gmx_simd_xor_f(a, b)       vec_xor(a, b)
#define gmx_simd_rsqrt_f(a)       vec_rsqrte(a)
#define gmx_simd_rcp_f(a)         vec_re(a)
#define gmx_simd_fabs_f(a)        vec_abs(a)
#define gmx_simd_fneg_f           gmx_simd_fneg_ibm_qpx
#define gmx_simd_max_f(a, b)       vec_sel(b, a, vec_sub(a, b))
#define gmx_simd_min_f(a, b)       vec_sel(b, a, vec_sub(b, a))
#define gmx_simd_round_f(a)       vec_round(a)
#define gmx_simd_trunc_f(a)       vec_trunc(a)
#define gmx_simd_fraction_f(x)    vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_f(a) gmx_simd_get_exponent_ibm_qpx(a)
#define gmx_simd_get_mantissa_f(a) gmx_simd_get_mantissa_ibm_qpx(a)
#define gmx_simd_set_exponent_f(a) gmx_simd_set_exponent_ibm_qpx(a)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t         vector4double
#ifdef NDEBUG
#    define gmx_simd_load_fi(m)   vec_ldi(0, (int *)(m))
#else
#    define gmx_simd_load_fi(m)   vec_ldia(0, (int *)(m))
#endif
#define gmx_simd_set1_fi(i)       gmx_simd_set1_int_ibm_qpx(i)
#define gmx_simd_store_fi(m, x)    vec_st(x, 0, (int *)(m))
#define gmx_simd_setzero_fi       gmx_simd_setzero_ibm_qpx
#define gmx_simd_cvt_f2i(a, b)     vec_ctiw(a, b)
#define gmx_simd_cvtt_f2i(a, b)    vec_ctiwz(a, b)
#define gmx_simd_cvt_i2f(a, b)     vec_cfid(a, b)
/* Integer simd extract not available */
/* Integer logical ops on gmx_simd_fint32_t not supported */
/* Integer arithmetic ops on gmx_simd_fint32_t not supported */
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t          vector4double
#define gmx_simd_cmpeq_f(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_f(a, b)     vec_cmplt((a), (b))
#define gmx_simd_cmple_f(a, b)     gmx_simd_or_fb(vec_cmpeq(a, b), vec_cmplt(a, b))
#define gmx_simd_and_fb(a, b)      vec_logical(a, b, 0x1)
#define gmx_simd_or_fb(a, b)       vec_logical(a, b, 0x7)
#define gmx_simd_anytrue_fb(a)    gmx_simd_anytrue_bool_ibm_qpx(a)
#define gmx_simd_blendzero_f(a, sel) vec_sel(vec_splats(0.0), a, sel)
#define gmx_simd_blendv_f(a, b, sel)  vec_sel(a, b, sel)
/* Boolean & comparison operations on gmx_simd_fint32_t not supported */
/* Conversions between different booleans not supported */

static __attribute__((always_inline)) vector4double
gmx_simd_fneg_ibm_qpx(vector4double a)
{
    return vec_neg(a);
}
/****************************************************
 *      DOUBLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_double_t         vector4double
#ifdef NDEBUG
#    define gmx_simd_load_d(m)    vec_ld(0, (double *)(m))
#    define gmx_simd_load1_d(m)   vec_lds(0, (double *)(m))
#    define gmx_simd_store_d(m, a) vec_st(a, 0, (double *)(m))
#else
#    define gmx_simd_load_d(m)    vec_lda(0, (double *)(m))
#    define gmx_simd_load1_d(m)   vec_ldsa(0, (double *)(m))
#    define gmx_simd_store_d(m, a) vec_sta(a, 0, (double *)(m))
#endif
#define gmx_simd_set1_d           vec_splats
/* No unaligned load/store */
#define gmx_simd_setzero_d        gmx_simd_setzero_ibm_qpx
#define gmx_simd_add_d            vec_add
#define gmx_simd_sub_d            vec_sub
#define gmx_simd_mul_d            vec_mul
#define gmx_simd_fmadd_d(a, b, c)   vec_madd(a, b, c)
#define gmx_simd_fmsub_d(a, b, c)   vec_msub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b+c=-(a*b-c) is "nmsub" */
#define gmx_simd_fnmadd_d(a, b, c)  vec_nmsub(a, b, c)
/* IBM uses an alternative FMA definition, so -a*b-c=-(a*b+c) is "nmadd" */
#define gmx_simd_fnmsub_d(a, b, c)  vec_nmadd(a, b, c)
#define gmx_simd_and_d            vec_and
#define gmx_simd_andnot_d(a, b)    vec_andc(b, a)
#define gmx_simd_or_d             vec_or
#define gmx_simd_xor_d            vec_xor
#define gmx_simd_rsqrt_d(a)       vec_rsqrte(a)
#define gmx_simd_rcp_d            vec_re
#define gmx_simd_fabs_d           vec_abs
#define gmx_simd_fneg_d           gmx_simd_fneg_ibm_qpx
#define gmx_simd_max_d            vec_sel(b, a, vec_sub(a, b))
#define gmx_simd_min_d            vec_sel(b, a, vec_sub(b, a))
#define gmx_simd_round_d          vec_round
#define gmx_simd_trunc_d          vec_trunc
#define gmx_simd_fraction_d(x)    vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_d   gmx_simd_get_exponent_ibm_qpx
#define gmx_simd_get_mantissa_d   gmx_simd_get_mantissa_ibm_qpx
#define gmx_simd_set_exponent_d   gmx_simd_set_exponent_ibm_qpx
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_dint32_t         vector4double
#define gmx_simd_load_di(m)       vec_ldia(0, (int *)(m))
#define gmx_simd_set1_di          gmx_simd_set1_int_ibm_qpx
#define gmx_simd_store_di(m, x)    vec_st(x, 0, (int *)(m))
#define gmx_simd_setzero_di       gmx_simd_setzero_ibm_qpx
#define gmx_simd_cvt_d2i          vec_ctiw
#define gmx_simd_cvtt_d2i         vec_ctiwz
#define gmx_simd_cvt_i2d          vec_cfid
/* Integer simd extract not available */
/* Integer logical ops on gmx_simd_fint32_t not supported */
/* Integer arithmetic ops on gmx_simd_fint32_t not supported */
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_dbool_t          vector4double
#define gmx_simd_cmpeq_d          vec_cmpeq
#define gmx_simd_cmplt_d          vec_cmplt
#define gmx_simd_cmple_d(a, b)     gmx_simd_or_fb(vec_cmpeq(a, b), vec_cmplt(a, b))
#define gmx_simd_and_db(a, b)      vec_logical(a, b, 0x1)
#define gmx_simd_or_db(a, b)       vec_logical(a, b, 0x7)
#define gmx_simd_anytrue_db       gmx_simd_anytrue_bool_ibm_qpx
#define gmx_simd_blendzero_d(a, sel) vec_sel(vec_splats(0.0), a, sel)
#define gmx_simd_blendv_d         vec_sel
/* Boolean & comparison operations on gmx_simd_fint32_t not supported */
/* Conversions between different booleans not supported */


/****************************************************
 * IMPLEMENTATION HELPER FUNCTIONS                  *
 ****************************************************/
static __attribute__((always_inline)) vector4double
gmx_simd_setzero_ibm_qpx(void)
{
    return vec_splats(0.0);
}

static __attribute__((always_inline)) vector4double
gmx_simd_get_exponent_ibm_qpx(vector4double x)
{
    const long long int    expmask = 0x7ff0000000000000LL;
    const long long int    expbase = 1023;

    union {
        vector4double v; long long i[4];
    } conv;

    /* Store to memory */
    vec_st(x, 0, (double *)&conv.v);
    /* Perform integer arithmetics in general registers. */
    conv.i[0] = ((conv.i[0] & expmask) >> 52) - expbase;
    conv.i[1] = ((conv.i[1] & expmask) >> 52) - expbase;
    conv.i[2] = ((conv.i[2] & expmask) >> 52) - expbase;
    conv.i[3] = ((conv.i[3] & expmask) >> 52) - expbase;
    /* Reload from memory */
    return vec_ld(0, (double *)&conv.v);
}

static __attribute__((always_inline)) vector4double
gmx_simd_get_mantissa_ibm_qpx(vector4double x)
{
    /* We cannot represent the mask 0x7ff0000000000000 directly as a double,
     * but the expression (0x5550000000000000 | 0x2aa0000000000000) works.
     * This will work for single-precision QPX too, since the floating-point
     * numbers will have been converted to double upon loading.
     */
    const vector4double expmask = vec_or(vec_splats(8.9589789687112168422e+102),
                                         vec_splats(2.2323972485981933246e-103));
    const vector4double one     = vec_splats(1.0);

    /* Get mantissa */
    x = vec_andc(x, expmask);
    /* Reset zero (but correctly biased) exponent */
    return vec_or(x, one);
}

static __attribute__((always_inline)) vector4double
gmx_simd_set_exponent_ibm_qpx(vector4double x)
{
    const long long int    expbase = 1023;
    union {
        vector4double v; long long i[4];
    } conv;

    /* Store to memory for shifts */
    vec_st(x, 0, (double *)&conv.v);
    /* Perform integer arithmetics in general registers. */
    conv.i[0] = (conv.i[0] + expbase) << 52;
    conv.i[1] = (conv.i[1] + expbase) << 52;
    conv.i[2] = (conv.i[2] + expbase) << 52;
    conv.i[3] = (conv.i[3] + expbase) << 52;
    /* Reload from memory */
    return vec_ld(0, (double *)&conv.v);
}

static __attribute__((always_inline)) vector4double
gmx_simd_set1_int_ibm_qpx(int i)
{
    union {
        vector4double v; long long i[4];
    } conv;

    conv.i[0] = i;
    /* Reload from memory */
    return vec_splat(vec_ldia(0, (int *)&conv.v), 0);
}

static __attribute__((always_inline)) int
gmx_simd_anytrue_db_ibm_qpx(vector4double a)
{
    vector4double b = vec_sldw(a, a, 2);

    a = vec_logical(a, b, 0x7);
    b = vec_sldw(a, a, 1);
    a = vec_logical(a, b, 0x7);
    return (vec_extract(a, 0) > 0);
}

static __attribute__((always_inline)) int
gmx_simd_anytrue_ibm_qpx(vector4double a)
{
    vector4double b = vec_sldw(a, a, 2);

    a = vec_logical(a, b, 0x7);
    b = vec_sldw(a, a, 1);
    a = vec_logical(a, b, 0x7);
    return (vec_extract(a, 0) > 0);
}

/* QPX is already 4-wide both in single and double, so just reuse for SIMD4 */

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
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_ibm_qpx
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
#define gmx_simd4_blendv_f               gmx_simd_blendv_f
/* DOUBLE */
#define gmx_simd4_double_t               gmx_simd_double_t
#define gmx_simd4_load_d                 gmx_simd_load_d
#define gmx_simd4_load1_d                gmx_simd_load1_d
#define gmx_simd4_set1_d                 gmx_simd_set1_d
#define gmx_simd4_store_d                gmx_simd_store_d
#define gmx_simd4_loadu_d                gmx_simd_loadu_d
#define gmx_simd4_storeu_d               gmx_simd_storeu_d
#define gmx_simd4_setzero_d              gmx_simd_setzero_d
#define gmx_simd4_add_d                  gmx_simd_add_d
#define gmx_simd4_sub_d                  gmx_simd_sub_d
#define gmx_simd4_mul_d                  gmx_simd_mul_d
#define gmx_simd4_fmadd_d                gmx_simd_fmadd_d
#define gmx_simd4_fmsub_d                gmx_simd_fmsub_d
#define gmx_simd4_fnmadd_d               gmx_simd_fnmadd_d
#define gmx_simd4_fnmsub_d               gmx_simd_fnmsub_d
#define gmx_simd4_and_d                  gmx_simd_and_d
#define gmx_simd4_andnot_d               gmx_simd_andnot_d
#define gmx_simd4_or_d                   gmx_simd_or_d
#define gmx_simd4_xor_d                  gmx_simd_xor_d
#define gmx_simd4_rsqrt_d                gmx_simd_rsqrt_d
#define gmx_simd4_rcp_d                  gmx_simd_rcp_d
#define gmx_simd4_fabs_d                 gmx_simd_fabs_d
#define gmx_simd4_fneg_d                 gmx_simd_fneg_d
#define gmx_simd4_max_d                  gmx_simd_max_d
#define gmx_simd4_min_d                  gmx_simd_min_d
#define gmx_simd4_round_d                gmx_simd_round_d
#define gmx_simd4_trunc_d                gmx_simd_trunc_d
#define gmx_simd4_fraction_d             gmx_simd_fraction_d
#define gmx_simd4_get_exponent_d         gmx_simd_get_exponent_d
#define gmx_simd4_get_mantissa_d         gmx_simd_get_mantissa_d
#define gmx_simd4_set_exponent_d         gmx_simd_set_exponent_d
#define gmx_simd4_dotproduct3_d          gmx_simd4_dotproduct3_d_ibm_qpx
#define gmx_simd4_dint32_t               gmx_simd_dint32_t
#define gmx_simd4_load_di                gmx_simd_load_di
#define gmx_simd4_load1_di               gmx_simd_load1_di
#define gmx_simd4_set1_di                gmx_simd_set1_di
#define gmx_simd4_store_di               gmx_simd_store_di
#define gmx_simd4_loadu_di               gmx_simd_loadu_di
#define gmx_simd4_storeu_di              gmx_simd_storeu_di
#define gmx_simd4_setzero_di             gmx_simd_setzero_di
#define gmx_simd4_cvt_d2i                gmx_simd_cvt_d2i
#define gmx_simd4_cvtt_d2i               gmx_simd_cvtt_d2i
#define gmx_simd4_cvt_i2f                gmx_simd_cvt_i2f
#define gmx_simd4_dbool_t                gmx_simd_dbool_t
#define gmx_simd4_cmpeq_d                gmx_simd_cmpeq_d
#define gmx_simd4_cmplt_d                gmx_simd_cmplt_d
#define gmx_simd4_cmple_d                gmx_simd_cmple_d
#define gmx_simd4_and_db                 gmx_simd_and_db
#define gmx_simd4_or_db                  gmx_simd_or_db
#define gmx_simd4_anytrue_db             gmx_simd_anytrue_db
#define gmx_simd4_blendzero_d            gmx_simd_blendzero_d
#define gmx_simd4_blendv_d               gmx_simd_blendv_d

static __attribute__((always_inline)) float
gmx_simd4_dotproduct3_f_ibm_qpx(vector4double a, vector4double b)
{
    float   data[GMX_SIMD_FLOAT_WIDTH*2];
    float * pdata = gmx_simd_align_f(data);
    a = vec_mul(a, b);
    vec_st(a, 0, pdata);
    return pdata[0]+pdata[1]+pdata[2];
}

static __attribute__((always_inline)) double
gmx_simd4_dotproduct3_d_ibm_qpx(vector4double a, vector4double b)
{
    double   data[GMX_SIMD_DOUBLE_WIDTH*2];
    double * pdata = gmx_simd_align_d(data);
    a = vec_mul(a, b);
    vec_st(a, 0, pdata);
    return pdata[0]+pdata[1]+pdata[2];
}

/* Function to check whether SIMD operations have resulted in overflow.
 * For now, this is unfortunately a dummy for this architecture.
 */
static int
gmx_simd_check_and_reset_overflow(void)
{
    return 0;
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_H */
