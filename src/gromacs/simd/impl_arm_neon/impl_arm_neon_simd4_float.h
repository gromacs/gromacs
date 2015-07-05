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

#ifndef GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H

#include <math.h>

#include <arm_neon.h>

#include "impl_arm_neon_common.h"
#include "impl_arm_neon_simd_float.h"

/* ARM 32-bit Neon is already 4-wide in single, so just reuse float type for SIMD4 */
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
#define gmx_simd4_fabs_f                 gmx_simd_fabs_f
#define gmx_simd4_fneg_f                 gmx_simd_fneg_f
#define gmx_simd4_max_f                  gmx_simd_max_f
#define gmx_simd4_min_f                  gmx_simd_min_f
#define gmx_simd4_round_f                gmx_simd_round_f
#define gmx_simd4_trunc_f                gmx_simd_trunc_f
#define gmx_simd4_dotproduct3_f          gmx_simd4_dotproduct3_f_arm_neon
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

/* SIMD4 Dotproduct helper function */
static gmx_inline float
gmx_simd4_dotproduct3_f_arm_neon(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    c = gmx_simd_mul_f(a, b);
    /* set 4th element to 0, then add all of them */
    c = vsetq_lane_f32(0.0f, c, 3);
    return gmx_simd_reduce_f_arm_neon(c);
}

#endif /* GMX_SIMD_IMPL_ARM_NEON_SIMD4_FLOAT_H */
