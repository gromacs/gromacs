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
 * \ref gmx::Simd4Real instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define Simd4Double    SimdDouble

/*! \brief Double precision SIMD4 load aligned.
 * \copydetails simdLoadD
 */
#    define simd4LoadD     simdLoadD

/*! \brief Double precision SIMD4 load single value to all elements.
 * \copydetails simdLoad1D
 */
#    define simd4Load1D    simdLoad1D

/*! \brief Double precision SIMD4 set all elements from value.
 * \copydetails simdSet1D
 */
#    define simd4Set1D     simdSet1D

/*! \brief Double precision SIMD4 store to aligned memory.
 * \copydetails simdStoreD
 */
#    define simd4StoreD   simdStoreD

/*! \brief Load unaligned SIMD4 double.
 * \copydetails simdLoadUD
 */
#    define simd4LoadUD   simdLoadUD

/*! \brief Store unaligned SIMD4 double.
 * \copydetails simdStoreUD
 */
#    define simd4StoreUD  simdStoreUD

/*! \brief Set all elements in SIMD4 double to 0.0.
 * \copydetails simdSetZeroD
 */
#    define simd4SetZeroD simdSetZeroD

/*! \brief Bitwise and for two SIMD4 double variables.
 * \copydetails simdAndD
 */
#    define simd4AndD     simdAndD

/*! \brief Bitwise andnot for SIMD4 double. c=(~a) & b.
 * \copydetails simdAndNotD
 */
#    define simd4AndNotD  simdAndNotD

/*! \brief Bitwise or for SIMD4 double.
 * \copydetails gmx::simdOrD
 */
#    define simd4OrD      simdOrD

/*! \brief Bitwise xor for SIMD4 double.
 * \copydetails gmx::simdXorD
 */
#    define simd4XorD     simdXorD

/*! \brief Add two SIMD4 double values.
 * \copydetails gmx::simdAddD
 */
#    define simd4AddD     simdAddD

/*! \brief Subtract two SIMD4 double values.
 * \copydetails gmx::simdSubD
 */
#    define simd4SubD     simdSubD

/*! \brief Multiply two SIMD4 double values.
 * \copydetails gmx::simdMulD
 */
#    define simd4MulD     simdMulD

/*! \brief Fused-multiply-add for SIMD4 double. Result is a*b+c.
 * \copydetails gmx::simdFmaddD
 */
#    define simd4FmaddD   simdFmaddD

/*! \brief Fused-multiply-subtract for SIMD4 double. Result is a*b-c.
 * \copydetails gmx::simdFmsubD
 */
#    define simd4FmsubD   simdFmsubD

/*! \brief Fused-negated-multiply-add for SIMD4 double. Result is -a*b+c.
 * \copydetails gmx::simdFnmaddD
 */
#    define simd4FnmaddD  simdFnmaddD

/*! \brief Fused-negated-multiply-sub for SIMD4 double. Result is -a*b-c.
 * \copydetails gmx::simdFnmsubD
 */
#    define simd4FnmsubD  simdFnmsubD

/*! \brief SIMD4 double 1.0/sqrt(x) lookup.
 * \copydetails gmx::simdRsqrtD
 */
#    define simd4RsqrtD   simdRsqrtD

/*! \brief SIMD4 double Floating-point fabs().
 * \copydetails gmx::simdAbsD
 */
#    define simd4AbsD    simdAbsD

/*! \brief SIMD4 double floating-point negate.
 * \copydetails gmx::simdNegD
 */
#    define simd4NegD    simdNegD

/*! \brief Set each SIMD4 element to the largest from two variables.
 * \copydetails gmx::simdMaxD
 */
#    define simd4MaxD     simdMaxD

/*! \brief Set each SIMD4 element to the smallest from two variables.
 * \copydetails gmx::simdMinD
 */
#    define simd4MinD     simdMinD

/*!  \brief Round SIMD4 double to nearest integer value (in floating-point format).
 * \copydetails gmx::simdRoundD
 */
#    define simd4RoundD   simdRoundD

/*! \brief Truncate SIMD4 double, i.e. round towards zero.
 * \copydetails gmx::simdTruncD
 */
#    define simd4TruncD   simdTruncD

/*! \brief Return dot product of two double precision SIMD4 variables.
 *
 * \copydetails simd4DotProductF
 */
static inline double
simd4DotProductD(SimdDouble a, SimdDouble b)
{
    return a.r[0]*b.r[0]+a.r[1]*b.r[1]+a.r[2]*b.r[2];
}

/*! \brief SIMD4 variable type to use for logical comparisons on doubles.
 * \copydetails SimdDBool
 */
#    define Simd4DBool   SimdDBool

/*! \brief Equality comparison of two double precision SIMD4 values.
 * \copydetails simdCmpEqD
 */
#    define simd4CmpEqD   simdCmpEqD

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails simdCmpLtD
 */
#    define simd4CmpLtD   simdCmpLtD

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails simdCmpLeD
 */
#    define simd4CmpLeD   simdCmpLeD

/*! \brief Logical AND on double SIMD4 booleans.
 * \copydetails simdAndDB
 */
#    define simd4AndDB simdAndDB

/*! \brief Logical OR on double SIMD4 booleans.
 * \copydetails simdOrDB
 */
#    define simd4OrDB simdOrDB

/*! \brief Returns non-zero if any of the SIMD4 booleans in x is True.
 * \copydetails simdAnyTrueDB
 */
#    define simd4AnyTrueDB simdAnyTrueDB

/*! \brief Select from double precision SIMD4 variable where boolean is true.
 * \copydetails simdMaskD
 */
#    define simd4MaskD simdMaskD

/*! \brief Select from double precision SIMD4 variable where boolean is false.
 * \copydetails simdMaskNotD
 */
#    define simd4MaskNotD simdMaskNotD

/*! \brief Vector-blend instruction for SIMD4 double.
 * \copydetails simdBlendD
 */
#    define simd4BlendD  simdBlendD

/*! \brief Return sum of all elements in SIMD4 double.
 * \copydetails simdReduceD
 */
#    define simd4ReduceD  simdReduceD

#else /* GMX_SIMD4_DOUBLE_WIDTH!=4 */
#    undef  GMX_SIMD4_HAVE_DOUBLE
#    define GMX_SIMD4_HAVE_DOUBLE      0
#endif

/*! \} */


/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H */
