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

#ifndef GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H

#include "config.h"

#include <math.h>
#ifdef __clang__
#include <qpxmath.h>
#endif

#include "impl_ibm_qpx_common.h"
#include "impl_ibm_qpx_simd_double.h"

/* QPX is already 4-wide both in single and double, so just reuse for SIMD4 */

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
#define gmx_simd4_blendnotzero_d         gmx_simd_blendnotzero_d
#define gmx_simd4_blendv_d               gmx_simd_blendv_d
#define gmx_simd4_reduce_d               gmx_simd_reduce_d

static __attribute__((always_inline)) double gmx_simdcall
gmx_simd4_dotproduct3_d_ibm_qpx(vector4double a, vector4double b)
{
    vector4double dp_sh0 = vec_mul(a, b);
    vector4double dp_sh1 = vec_sldw(dp_sh0, dp_sh0, 1);
    vector4double dp_sh2 = vec_sldw(dp_sh0, dp_sh0, 2);
    vector4double dp     = vec_add(dp_sh2, vec_add(dp_sh0, dp_sh1));

    return vec_extract(dp, 0);
}

#endif /* GMX_SIMD_IMPLEMENTATION_IBM_QPX_SIMD4_DOUBLE_H */
