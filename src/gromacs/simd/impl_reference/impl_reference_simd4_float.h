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
 * \ref gmx::Simd4Real instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define Simd4Float    SimdFloat

/*! \brief Load SIMD4 float from aligned memory.
 *  \copydetails simdLoadF
 */
#    define simd4LoadF     simdLoadF

/*! \brief Set all elements of SIMD4 float from single pointer.
 *  \copydetails simdLoad1F
 */
#    define simd4Load1F    simdLoad1F

/*! \brief Set all SIMD4 float elements to the value r.
 *  \copydetails simdSet1F
 */
#    define simd4Set1F     simdSet1F

/*! \brief Store the contents of SIMD4 float pr to aligned memory m.
 *  \copydetails simdStoreF
 */
#    define simd4StoreF    simdStoreF

/*! \brief Load SIMD4 float from unaligned memory.
 * \copydetails simdLoadUF
 */
#    define simd4LoadUF    simdLoadUF

/*! \brief Store SIMD4 float to unaligned memory.
 * \copydetails simdStoreUF
 */
#    define simd4StoreUF   simdStoreUF

/*! \brief Set all SIMD4 float elements to 0.
 * \copydetails simdSetZeroF
 */
#    define simd4SetZeroF  simdSetZeroF

/*! \brief Bitwise and for two SIMD4 float variables.
 * \copydetails simdAndF
 */
#    define simd4AndF      simdAndF

/*! \brief Bitwise andnot for two SIMD4 float variables. c=(~a) & b.
 * \copydetails simdAndNotF
 */
#    define simd4AndNotF   simdAndNotF

/*! \brief Bitwise or for two SIMD4 float variables.
 * \copydetails simdOrF
 */
#    define simd4OrF       simdOrF

/*! \brief Bitwise xor for two SIMD4 float variables.
 * \copydetails simdXorF
 */
#    define simd4XorF      simdXorF

/*! \brief Add two SIMD4 float variables.
 * \copydetails simdAddF
 */
#    define simd4AddF      simdAddF

/*! \brief Subtract two SIMD4 float variables.
 * \copydetails simdSubF
 */
#    define simd4SubF      simdSubF

/*! \brief Multiply two SIMD4 float variables.
 * \copydetails simdMulF
 */
#    define simd4MulF      simdMulF

/*! \brief Fused-multiply-add for SIMD4 float. Result is a*b+c.
 * \copydetails simdFmaddF
 */
#    define simd4FmaddF    simdFmaddF

/*! \brief Fused-multiply-subtract for SIMD4 float. Result is a*b-c.
 * \copydetails simdFmsubF
 */
#    define simd4FmsubF    simdFmsubF

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b+c.
 * \copydetails simdFnmaddF
 */
#    define simd4FnmaddF   simdFnmaddF

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b-c.
 * \copydetails simdFnmsubF
 */
#    define simd4FnmsubF   simdFnmsubF

/*! \brief Lookup of approximate 1/sqrt(x) for SIMD4 float.
 * \copydetails simdRsqrtF
 */
#    define simd4RsqrtF    simdRsqrtF

/*! \brief Floating-point absolute value for SIMD4 float.
 * \copydetails simdAbsF
 */
#    define simd4AbsF     simdAbsF

/*! \brief Floating-point negate for SIMD4 float.
 * \copydetails simdNegF
 */
#    define simd4NegF     simdNegF

/*! \brief Set each SIMD4 float element to the largest from two variables.
 * \copydetails simdMaxF
 */
#    define simd4MaxF      simdMaxF

/*! \brief Set each SIMD4 float element to the smallest from two variables.
 * \copydetails simdMinF
 */
#    define simd4MinF      simdMinF

/*! \brief Round to nearest integer value for SIMD4 float.
 * \copydetails simdRoundF
 */
#    define simd4RoundF    simdRoundF

/*! \brief Round to largest integral value for SIMD4 float.
 * \copydetails simdTruncF
 */
#    define simd4TruncF    simdTruncF

/*! \brief Return dot product of two single precision SIMD4 variables.
 *
 * The dot product is calculated between the first three elements in the two
 * vectors, while the fourth is ignored. The result is returned as a scalar.
 *
 * \param a vector1
 * \param b vector2
 * \result a[0]*b[0]+a[1]*b[1]+a[2]*b[2], returned as scalar. Last element is ignored.
 */
static inline float
simd4DotProductF(SimdFloat a, SimdFloat b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/*! \brief SIMD4 variable type to use for logical comparisons on floats.
 * \copydetails SimdFBool
 */
#    define Simd4FBool   SimdFBool

/*! \brief Equality comparison of two single precision SIMD4.
 * \copydetails simdCmpEqF
 */
#    define simd4CmpEqF   simdCmpEqF

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLtF
 */
#    define simd4CmpLtF   simdCmpLtF

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLeF
 */
#    define simd4CmpLeF   simdCmpLeF

/*! \brief Logical AND on float SIMD4 booleans.
 * \copydetails simdAndFB
 */
#    define simd4AndFB simdAndFB

/*! \brief Logical OR on float SIMD4 booleans.
 * \copydetails simdOrFB
 */
#    define simd4OrFB simdOrFB

/*! \brief Returns non-zero if any of the SIMD4 boolean in x is True.
 * \copydetails simdAnyTrueFB
 */
#    define simd4AnyTrueFB simdAnyTrueFB

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 * \copydetails simdMaskF
 */
#    define simd4MaskF simdMaskF

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 * \copydetails simdMaskNotF
 */
#    define simd4MaskNotF simdMaskNotF

/*! \brief Vector-blend instruction form SIMD4 float.
 * \copydetails simdBlendF
 */
#    define simd4BlendF  simdBlendF

/*! \brief Return sum of all elements in SIMD4 float.
 * \copydetails simdReduceF
 */
#    define simd4ReduceF  simdReduceF

#else /* GMX_SIMD_FLOAT_WIDTH!=4 */
#    define GMX_SIMD4_HAVE_FLOAT    0
#endif

/*! \} */

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H */
