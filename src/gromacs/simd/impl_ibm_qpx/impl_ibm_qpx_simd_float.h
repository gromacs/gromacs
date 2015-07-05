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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"

/****************************************************
 *      SINGLE PRECISION SIMD IMPLEMENTATION        *
 ****************************************************/
#define gmx_simd_float_t          vector4double
#ifdef NDEBUG
#    define gmx_simd_load_f(m)    vec_ld(0, (float *)(m))
#    define gmx_simd_store_f(m, a) vec_st(a, 0, (float *)(m))
#else
#    define gmx_simd_load_f(m)    vec_lda(0, (float *)(m))
#    define gmx_simd_store_f(m, a) vec_sta(a, 0, (float *)(m))
#endif
#    define gmx_simd_load1_f(m)   vec_lds(0, (float *)(m))
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
/* gmx_simd_and_f not supported - no bitwise logical ops */
/* gmx_simd_andnot_f not supported - no bitwise logical ops */
/* gmx_simd_or_f not supported - no bitwise logical ops */
/* gmx_simd_xor_f not supported - no bitwise logical ops */
#define gmx_simd_rsqrt_f(a)       vec_rsqrte(a)
#define gmx_simd_rcp_f(a)         vec_re(a)
#define gmx_simd_fabs_f(a)        vec_abs(a)
#define gmx_simd_fneg_f           gmx_simd_fneg_ibm_qpx
#define gmx_simd_max_f(a, b)       vec_sel(b, a, vec_sub(a, b))
#define gmx_simd_min_f(a, b)       vec_sel(b, a, vec_sub(b, a))
/* Note: It is critical to use vec_cfid(vec_ctid(a)) for the implementation
 * of gmx_simd_round_f(), since vec_round() does not adhere to the FP control
 * word rounding scheme. We rely on float-to-float and float-to-integer
 * rounding being the same for half-way values in a few algorithms.
 */
#define gmx_simd_round_f(a)       vec_cfid(vec_ctid(a))
#define gmx_simd_trunc_f(a)       vec_trunc(a)
#define gmx_simd_fraction_f(x)    vec_sub(x, vec_trunc(x))
#define gmx_simd_get_exponent_f(a) gmx_simd_get_exponent_ibm_qpx(a)
#define gmx_simd_get_mantissa_f(a) gmx_simd_get_mantissa_ibm_qpx(a)
#define gmx_simd_set_exponent_f(a) gmx_simd_set_exponent_ibm_qpx(a)
/* integer datatype corresponding to float: gmx_simd_fint32_t */
#define gmx_simd_fint32_t         vector4double
#ifdef NDEBUG
#    define gmx_simd_load_fi(m)   vec_ldia(0, (int *)(m))
#else
#    define gmx_simd_load_fi(m)   vec_ldiaa(0, (int *)(m))
#endif
#define gmx_simd_set1_fi(i)       gmx_simd_set1_int_ibm_qpx(i)
#define gmx_simd_store_fi(m, x)    vec_st(x, 0, (int *)(m))
#define gmx_simd_setzero_fi       gmx_simd_setzero_ibm_qpx
#define gmx_simd_cvt_f2i(a)       vec_ctiw(a)
#define gmx_simd_cvtt_f2i(a)      vec_ctiwz(a)
#define gmx_simd_cvt_i2f(a)       vec_cfid(a)
/* Integer simd extract not available */
/* Integer logical ops on gmx_simd_fint32_t not supported */
/* Integer arithmetic ops on gmx_simd_fint32_t not supported */
/* Boolean & comparison operations on gmx_simd_float_t */
#define gmx_simd_fbool_t          vector4double
#define gmx_simd_cmpeq_f(a, b)     vec_cmpeq(a, b)
#define gmx_simd_cmplt_f(a, b)     vec_cmplt((a), (b))
#define gmx_simd_cmple_f(a, b)     gmx_simd_or_fb(vec_cmpeq(a, b), vec_cmplt(a, b))
#define gmx_simd_and_fb(a, b)      vec_and(a, b)
#define gmx_simd_or_fb(a, b)       vec_or(a, b)
#define gmx_simd_anytrue_fb(a)    gmx_simd_anytrue_bool_ibm_qpx(a)
#define gmx_simd_blendzero_f(a, sel) vec_sel(vec_splats(0.0), a, sel)
#define gmx_simd_blendnotzero_f(a, sel) vec_sel(a, vec_splats(0.0), sel)
#define gmx_simd_blendv_f(a, b, sel)  vec_sel(a, b, sel)
#define gmx_simd_reduce_f(a)       gmx_simd_reduce_ibm_qpx(a)


/* Boolean & comparison operations on gmx_simd_fint32_t not supported */
/* Conversions between different booleans not supported */

/* Note: Since QPX registers are always double internally, the single
 * precision defines use several double precision helper functions.
 */

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD_FLOAT_H */
