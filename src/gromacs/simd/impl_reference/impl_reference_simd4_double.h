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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD4 double precision.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <math.h>

#include "impl_reference_common.h"
#include "impl_reference_simd_double.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name SIMD4. Constant width-4 double precision SIMD types and instructions
 * \{
 */

#if (GMX_SIMD_DOUBLE_WIDTH == 4) || defined DOXYGEN

/*! \brief SIMD4 double type. Available with \ref GMX_SIMD4_HAVE_DOUBLE.
 *
 * Unless you specifically want a double-precision type you should check
 * \ref gmx_simd4_real_t instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define gmx_simd4_double_t   gmx_simd_double_t

/*! \brief Double precision SIMD4 load aligned.
 * \copydetails gmx_simd_load_d
 */
#    define gmx_simd4_load_d     gmx_simd_load_d

/*! \brief Double precision SIMD4 load single value to all elements.
 * \copydetails gmx_simd_load1_d
 */
#    define gmx_simd4_load1_d    gmx_simd_load1_d

/*! \brief Double precision SIMD4 set all elements from value.
 * \copydetails gmx_simd_set1_d
 */
#    define gmx_simd4_set1_d     gmx_simd_set1_d

/*! \brief Double precision SIMD4 store to aligned memory.
 * \copydetails gmx_simd_store_d
 */
#    define gmx_simd4_store_d   gmx_simd_store_d

/*! \brief Load unaligned SIMD4 double.
 * \copydetails gmx_simd_loadu_d
 */
#    define gmx_simd4_loadu_d   gmx_simd_loadu_d

/*! \brief Store unaligned SIMD4 double.
 * \copydetails gmx_simd_storeu_d
 */
#    define gmx_simd4_storeu_d  gmx_simd_storeu_d

/*! \brief Set all elements in SIMD4 double to 0.0.
 * \copydetails gmx_simd_setzero_d
 */
#    define gmx_simd4_setzero_d gmx_simd_setzero_d

/*! \brief Bitwise and for two SIMD4 double variables.
 * \copydetails gmx_simd_and_d
 */
#    define gmx_simd4_and_d     gmx_simd_and_d

/*! \brief Bitwise andnot for SIMD4 double. c=(~a) & b.
 * \copydetails gmx_simd_andnot_d
 */
#    define gmx_simd4_andnot_d  gmx_simd_andnot_d

/*! \brief Bitwise or for SIMD4 double.
 * \copydetails gmx_simd_or_d
 */
#    define gmx_simd4_or_d      gmx_simd_or_d

/*! \brief Bitwise xor for SIMD4 double.
 * \copydetails gmx_simd_xor_d
 */
#    define gmx_simd4_xor_d     gmx_simd_xor_d

/*! \brief Add two SIMD4 double values.
 * \copydetails gmx_simd_add_d
 */
#    define gmx_simd4_add_d     gmx_simd_add_d

/*! \brief Subtract two SIMD4 double values.
 * \copydetails gmx_simd_sub_d
 */
#    define gmx_simd4_sub_d     gmx_simd_sub_d

/*! \brief Multiply two SIMD4 double values.
 * \copydetails gmx_simd_mul_d
 */
#    define gmx_simd4_mul_d     gmx_simd_mul_d

/*! \brief Fused-multiply-add for SIMD4 double. Result is a*b+c.
 * \copydetails gmx_simd_fmadd_d
 */
#    define gmx_simd4_fmadd_d   gmx_simd_fmadd_d

/*! \brief Fused-multiply-subtract for SIMD4 double. Result is a*b-c.
 * \copydetails gmx_simd_fmsub_d
 */
#    define gmx_simd4_fmsub_d   gmx_simd_fmsub_d

/*! \brief Fused-negated-multiply-add for SIMD4 double. Result is -a*b+c.
 * \copydetails gmx_simd_fnmadd_d
 */
#    define gmx_simd4_fnmadd_d  gmx_simd_fnmadd_d

/*! \brief Fused-negated-multiply-sub for SIMD4 double. Result is -a*b-c.
 * \copydetails gmx_simd_fnmsub_d
 */
#    define gmx_simd4_fnmsub_d  gmx_simd_fnmsub_d

/*! \brief SIMD4 double 1.0/sqrt(x) lookup.
 * \copydetails gmx_simd_rsqrt_d
 */
#    define gmx_simd4_rsqrt_d   gmx_simd_rsqrt_d

/*! \brief SIMD4 double Floating-point fabs().
 * \copydetails gmx_simd_fabs_d
 */
#    define gmx_simd4_fabs_d    gmx_simd_fabs_d

/*! \brief SIMD4 double floating-point negate.
 * \copydetails gmx_simd_fneg_d
 */
#    define gmx_simd4_fneg_d    gmx_simd_fneg_d

/*! \brief Set each SIMD4 element to the largest from two variables.
 * \copydetails gmx_simd_max_d
 */
#    define gmx_simd4_max_d     gmx_simd_max_d

/*! \brief Set each SIMD4 element to the smallest from two variables.
 * \copydetails gmx_simd_min_d
 */
#    define gmx_simd4_min_d     gmx_simd_min_d

/*!  \brief Round SIMD4 double to nearest integer value (in floating-point format).
 * \copydetails gmx_simd_round_d
 */
#    define gmx_simd4_round_d   gmx_simd_round_d

/*! \brief Truncate SIMD4 double, i.e. round towards zero.
 * \copydetails gmx_simd_trunc_d
 */
#    define gmx_simd4_trunc_d   gmx_simd_trunc_d

/*! \brief Return dot product of two double precision SIMD4 variables.
 * \copydetails gmx_simd_setzero_f
 */
static gmx_inline double
gmx_simd4_dotproduct3_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/*! \brief SIMD4 variable type to use for logical comparisons on doubles.
 * \copydetails gmx_simd_dbool_t
 */
#    define gmx_simd4_dbool_t   gmx_simd_dbool_t

/*! \brief Equality comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmpeq_d
 */
#    define gmx_simd4_cmpeq_d   gmx_simd_cmpeq_d

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmplt_d
 */
#    define gmx_simd4_cmplt_d   gmx_simd_cmplt_d

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmple_d
 */
#    define gmx_simd4_cmple_d   gmx_simd_cmple_d

/*! \brief Logical AND on double SIMD4 booleans.
 * \copydetails gmx_simd_and_db
 */
#    define gmx_simd4_and_db gmx_simd_and_db

/*! \brief Logical OR on double SIMD4 booleans.
 * \copydetails gmx_simd_or_db
 */
#    define gmx_simd4_or_db gmx_simd_or_db

/*! \brief Returns non-zero if any of the SIMD4 booleans in x is True.
 * \copydetails gmx_simd_anytrue_db
 */
#    define gmx_simd4_anytrue_db gmx_simd_anytrue_db

/*! \brief Select from double precision SIMD4 variable where boolean is true.
 * \copydetails gmx_simd_blendzero_d
 */
#    define gmx_simd4_blendzero_d gmx_simd_blendzero_d

/*! \brief Select from double precision SIMD4 variable where boolean is false.
 * \copydetails gmx_simd_blendnotzero_d
 */
#    define gmx_simd4_blendnotzero_d gmx_simd_blendnotzero_d

/*! \brief Vector-blend instruction for SIMD4 double.
 * \copydetails gmx_simd_blendv_d
 */
#    define gmx_simd4_blendv_d  gmx_simd_blendv_d

/*! \brief Return sum of all elements in SIMD4 double.
 * \copydetails gmx_simd_reduce_d
 */
#    define gmx_simd4_reduce_d  gmx_simd_reduce_d

#else /* GMX_SIMD4_DOUBLE_WIDTH!=4 */
#    define GMX_SIMD4_HAVE_DOUBLE      0
#endif

/*! \} */


/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H */
