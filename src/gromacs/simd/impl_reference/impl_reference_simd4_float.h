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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD4 single precision.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <math.h>

#include "impl_reference_common.h"
#include "impl_reference_simd_float.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Constant width-4 single precision SIMD types and instructions
 * \{
 */

#if (GMX_SIMD_FLOAT_WIDTH == 4) || defined DOXYGEN


/*! \brief SIMD4 float type. Available if \ref GMX_SIMD4_HAVE_FLOAT.
 *
 * Unless you specifically want a single-precision type you should check
 * \ref gmx_simd4_real_t instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define gmx_simd4_float_t    gmx_simd_float_t

/*! \brief Load SIMD4 float from aligned memory.
 *  \copydetails gmx_simd_load_f
 */
#    define gmx_simd4_load_f     gmx_simd_load_f

/*! \brief Set all elements of SIMD4 float from single pointer.
 *  \copydetails gmx_simd_load1_f
 */
#    define gmx_simd4_load1_f    gmx_simd_load1_f

/*! \brief Set all SIMD4 float elements to the value r.
 *  \copydetails gmx_simd_set1_f
 */
#    define gmx_simd4_set1_f     gmx_simd_set1_f

/*! \brief Store the contents of SIMD4 float pr to aligned memory m.
 *  \copydetails gmx_simd_store_f
 */
#    define gmx_simd4_store_f    gmx_simd_store_f

/*! \brief Load SIMD4 float from unaligned memory.
 * \copydetails gmx_simd_loadu_f
 */
#    define gmx_simd4_loadu_f    gmx_simd_loadu_f

/*! \brief Store SIMD4 float to unaligned memory.
 * \copydetails gmx_simd_storeu_f
 */
#    define gmx_simd4_storeu_f   gmx_simd_storeu_f

/*! \brief Set all SIMD4 float elements to 0.
 * \copydetails gmx_simd_setzero_f
 */
#    define gmx_simd4_setzero_f  gmx_simd_setzero_f

/*! \brief Bitwise and for two SIMD4 float variables.
 * \copydetails gmx_simd_and_f
 */
#    define gmx_simd4_and_f      gmx_simd_and_f

/*! \brief Bitwise andnot for two SIMD4 float variables. c=(~a) & b.
 * \copydetails gmx_simd_andnot_f
 */
#    define gmx_simd4_andnot_f   gmx_simd_andnot_f

/*! \brief Bitwise or for two SIMD4 float variables.
 * \copydetails gmx_simd_or_f
 */
#    define gmx_simd4_or_f       gmx_simd_or_f

/*! \brief Bitwise xor for two SIMD4 float variables.
 * \copydetails gmx_simd_xor_f
 */
#    define gmx_simd4_xor_f      gmx_simd_xor_f

/*! \brief Add two SIMD4 float variables.
 * \copydetails gmx_simd_add_f
 */
#    define gmx_simd4_add_f      gmx_simd_add_f

/*! \brief Subtract two SIMD4 float variables.
 * \copydetails gmx_simd_sub_f
 */
#    define gmx_simd4_sub_f      gmx_simd_sub_f

/*! \brief Multiply two SIMD4 float variables.
 * \copydetails gmx_simd_mul_f
 */
#    define gmx_simd4_mul_f      gmx_simd_mul_f

/*! \brief Fused-multiply-add for SIMD4 float. Result is a*b+c.
 * \copydetails gmx_simd_fmadd_f
 */
#    define gmx_simd4_fmadd_f    gmx_simd_fmadd_f

/*! \brief Fused-multiply-subtract for SIMD4 float. Result is a*b-c.
 * \copydetails gmx_simd_fmsub_f
 */
#    define gmx_simd4_fmsub_f    gmx_simd_fmsub_f

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b+c.
 * \copydetails gmx_simd_fnmadd_f
 */
#    define gmx_simd4_fnmadd_f   gmx_simd_fnmadd_f

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b-c.
 * \copydetails gmx_simd_fnmsub_f
 */
#    define gmx_simd4_fnmsub_f   gmx_simd_fnmsub_f

/*! \brief Lookup of approximate 1/sqrt(x) for SIMD4 float.
 * \copydetails gmx_simd_rsqrt_f
 */
#    define gmx_simd4_rsqrt_f    gmx_simd_rsqrt_f

/*! \brief Floating-point absolute value for SIMD4 float.
 * \copydetails gmx_simd_fabs_f
 */
#    define gmx_simd4_fabs_f     gmx_simd_fabs_f

/*! \brief Floating-point negate for SIMD4 float.
 * \copydetails gmx_simd_fneg_f
 */
#    define gmx_simd4_fneg_f     gmx_simd_fneg_f

/*! \brief Set each SIMD4 float element to the largest from two variables.
 * \copydetails gmx_simd_max_f
 */
#    define gmx_simd4_max_f      gmx_simd_max_f

/*! \brief Set each SIMD4 float element to the smallest from two variables.
 * \copydetails gmx_simd_min_f
 */
#    define gmx_simd4_min_f      gmx_simd_min_f

/*! \brief Round to nearest integer value for SIMD4 float.
 * \copydetails gmx_simd_round_f
 */
#    define gmx_simd4_round_f    gmx_simd_round_f

/*! \brief Round to largest integral value for SIMD4 float.
 * \copydetails gmx_simd_trunc_f
 */
#    define gmx_simd4_trunc_f    gmx_simd_trunc_f

/*! \brief Return dot product of two single precision SIMD4 variables.
 *
 * The dot product is calculated between the first three elements in the two
 * vectors, while the fourth is ignored. The result is returned as a scalar.
 *
 * \param a vector1
 * \param b vector2
 * \result a[0]*b[0]+a[1]*b[1]+a[2]*b[2], returned as scalar. Last element is ignored.
 */
static gmx_inline float
gmx_simd4_dotproduct3_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/*! \brief SIMD4 variable type to use for logical comparisons on floats.
 * \copydetails gmx_simd_fbool_t
 */
#    define gmx_simd4_fbool_t   gmx_simd_fbool_t

/*! \brief Equality comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmpeq_f
 */
#    define gmx_simd4_cmpeq_f   gmx_simd_cmpeq_f

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmplt_f
 */
#    define gmx_simd4_cmplt_f   gmx_simd_cmplt_f

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmple_f
 */
#    define gmx_simd4_cmple_f   gmx_simd_cmple_f

/*! \brief Logical AND on float SIMD4 booleans.
 * \copydetails gmx_simd_and_fb
 */
#    define gmx_simd4_and_fb gmx_simd_and_fb

/*! \brief Logical OR on float SIMD4 booleans.
 * \copydetails gmx_simd_or_fb
 */
#    define gmx_simd4_or_fb gmx_simd_or_fb

/*! \brief Returns non-zero if any of the SIMD4 boolean in x is True.
 * \copydetails gmx_simd_anytrue_fb
 */
#    define gmx_simd4_anytrue_fb gmx_simd_anytrue_fb

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 * \copydetails gmx_simd_blendzero_f
 */
#    define gmx_simd4_blendzero_f gmx_simd_blendzero_f

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 * \copydetails gmx_simd_blendnotzero_f
 */
#    define gmx_simd4_blendnotzero_f gmx_simd_blendnotzero_f

/*! \brief Vector-blend instruction form SIMD4 float.
 * \copydetails gmx_simd_blendv_f
 */
#    define gmx_simd4_blendv_f  gmx_simd_blendv_f

/*! \brief Return sum of all elements in SIMD4 float.
 * \copydetails gmx_simd_reduce_f
 */
#    define gmx_simd4_reduce_f  gmx_simd_reduce_f

#else /* GMX_SIMD_FLOAT_WIDTH!=4 */
#    define GMX_SIMD4_HAVE_FLOAT    0
#endif

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H */
