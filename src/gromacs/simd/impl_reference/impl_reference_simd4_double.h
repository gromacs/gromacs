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

#include <cmath>

#include "impl_reference_common.h"
#include "impl_reference_simd_double.h"

namespace gmx
{

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
typedef SimdDouble         Simd4Double;

/*! \brief Double precision SIMD4 load aligned.
 * \copydetails simdLoadD
 */
static inline Simd4Double
simd4LoadD(const double *m) { return simdLoadD(m); }

/*! \brief Double precision SIMD4 load single value to all elements.
 * \copydetails simdLoad1D
 */
static inline Simd4Double
simd4Load1D(const double *m) { return simdLoad1D(m); }

/*! \brief Double precision SIMD4 set all elements from value.
 * \copydetails simdSet1D
 */
static inline Simd4Double
simd4Set1D(double r) { return simdSet1D(r); }

/*! \brief Double precision SIMD4 store to aligned memory.
 * \copydetails simdStoreD
 */
static inline void
simd4StoreD(double *m, Simd4Double a) { simdStoreD(m, a); }

/*! \brief Load unaligned SIMD4 double.
 * \copydetails simdLoadUD
 */
static inline Simd4Double
simd4LoadUD(const double *m) { return simdLoadD(m); }

/*! \brief Store unaligned SIMD4 double.
 * \copydetails simdStoreUD
 */
static inline void
simd4StoreUD(double *m, Simd4Double a) { simdStoreD(m, a); }

/*! \brief Set all elements in SIMD4 double to 0.0.
 * \copydetails simdSetZeroD
 */
static inline Simd4Double
simd4SetZeroD() { return simdSetZeroD(); }

/*! \brief Bitwise and for two SIMD4 double variables.
 * \copydetails simdAndD
 */
static inline Simd4Double
simd4AndD(Simd4Double a, Simd4Double b) { return simdAndD(a, b); }

/*! \brief Bitwise andnot for SIMD4 double. c=(~a) & b.
 * \copydetails simdAndNotD
 */
static inline Simd4Double
simd4AndNotD(Simd4Double a, Simd4Double b) { return simdAndNotD(a, b); }

/*! \brief Bitwise or for SIMD4 double.
 * \copydetails gmx::simdOrD
 */
static inline Simd4Double
simd4OrD(Simd4Double a, Simd4Double b) { return simdOrD(a, b); }

/*! \brief Bitwise xor for SIMD4 double.
 * \copydetails gmx::simdXorD
 */
static inline Simd4Double
simd4XorD(Simd4Double a, Simd4Double b) { return simdXorD(a, b); }

/*! \brief Add two SIMD4 double values.
 * \copydetails gmx::simdAddD
 */
static inline Simd4Double
simd4AddD(Simd4Double a, Simd4Double b) { return simdAddD(a, b); }

/*! \brief Subtract two SIMD4 double values.
 * \copydetails gmx::simdSubD
 */
static inline Simd4Double
simd4SubD(Simd4Double a, Simd4Double b) { return simdSubD(a, b); }

/*! \brief Multiply two SIMD4 double values.
 * \copydetails gmx::simdMulD
 */
static inline Simd4Double
simd4MulD(Simd4Double a, Simd4Double b) { return simdMulD(a, b); }

/*! \brief Fused-multiply-add for SIMD4 double. Result is a*b+c.
 * \copydetails gmx::simdFmaddD
 */
static inline Simd4Double
simd4FmaddD(Simd4Double a, Simd4Double b, Simd4Double c) { return simdFmaddD(a, b, c); }

/*! \brief Fused-multiply-subtract for SIMD4 double. Result is a*b-c.
 * \copydetails gmx::simdFmsubD
 */
static inline Simd4Double
simd4FmsubD(Simd4Double a, Simd4Double b, Simd4Double c) { return simdFmsubD(a, b, c); }

/*! \brief Fused-negated-multiply-add for SIMD4 double. Result is -a*b+c.
 * \copydetails gmx::simdFnmaddD
 */
static inline Simd4Double
simd4FnmaddD(Simd4Double a, Simd4Double b, Simd4Double c) { return simdFnmaddD(a, b, c); }

/*! \brief Fused-negated-multiply-sub for SIMD4 double. Result is -a*b-c.
 * \copydetails gmx::simdFnmsubD
 */
static inline Simd4Double
simd4FnmsubD(Simd4Double a, Simd4Double b, Simd4Double c) { return simdFnmsubD(a, b, c); }

/*! \brief SIMD4 double 1.0/sqrt(x) lookup.
 * \copydetails gmx::simdRsqrtD
 */
static inline Simd4Double
simd4RsqrtD(Simd4Double x) { return simdRsqrtD(x); }

/*! \brief SIMD4 double Floating-point fabs().
 * \copydetails gmx::simdAbsD
 */
static inline Simd4Double
simd4AbsD(Simd4Double a) { return simdAbsD(a); }

/*! \brief SIMD4 double floating-point negate.
 * \copydetails gmx::simdNegD
 */
static inline Simd4Double
simd4NegD(Simd4Double a) { return simdNegD(a); }

/*! \brief Set each SIMD4 element to the largest from two variables.
 * \copydetails gmx::simdMaxD
 */
static inline Simd4Double
simd4MaxD(Simd4Double a, Simd4Double b) { return simdMaxD(a, b); }

/*! \brief Set each SIMD4 element to the smallest from two variables.
 * \copydetails gmx::simdMinD
 */
static inline Simd4Double
simd4MinD(Simd4Double a, Simd4Double b) { return simdMinD(a, b); }

/*!  \brief Round SIMD4 double to nearest integer value (in floating-point format).
 * \copydetails gmx::simdRoundD
 */
static inline Simd4Double
simd4RoundD(Simd4Double a) { return simdRoundD(a); }

/*! \brief Truncate SIMD4 double, i.e. round towards zero.
 * \copydetails gmx::simdTruncD
 */
static inline Simd4Double
simd4TruncD(Simd4Double a) { return simdTruncD(a); }

/*! \brief Return dot product of two double precision SIMD4 variables.
 *
 * \copydetails simd4DotProductF
 */
static inline double
simd4DotProductD(SimdDouble a, SimdDouble b)
{
    return a.r[0] * b.r[0] + a.r[1] * b.r[1] + a.r[2] * b.r[2];
}

/*! \brief SIMD4 double transpose
 *
 * \param[in,out] v0  Row 0 on input, column 0 on output
 * \param[in,out] v1  Row 1 on input, column 1 on output
 * \param[in,out] v2  Row 2 on input, column 2 on output
 * \param[in,out] v3  Row 3 on input, column 3 on output
 *
 * This is only available in C++.
 */
static inline void
simd4Transpose(Simd4Double &v0, Simd4Double &v1,
               Simd4Double &v2, Simd4Double &v3)
{
    Simd4Double t0, t1, t2, t3;
    t0      = v0;
    t1      = v1;
    t2      = v2;
    t3      = v3;
    v0.r[0] = t0.r[0];
    v0.r[1] = t1.r[0];
    v0.r[2] = t2.r[0];
    v0.r[3] = t3.r[0];
    v1.r[0] = t0.r[1];
    v1.r[1] = t1.r[1];
    v1.r[2] = t2.r[1];
    v1.r[3] = t3.r[1];
    v2.r[0] = t0.r[2];
    v2.r[1] = t1.r[2];
    v2.r[2] = t2.r[2];
    v2.r[3] = t3.r[2];
    v3.r[0] = t0.r[3];
    v3.r[1] = t1.r[3];
    v3.r[2] = t2.r[3];
    v3.r[3] = t3.r[3];
}

/*! \brief SIMD4 variable type to use for logical comparisons on doubles.
 * \copydetails SimdDBool
 */
typedef SimdDBool        Simd4DBool;

/*! \brief Equality comparison of two double precision SIMD4 values.
 * \copydetails simdCmpEqD
 */
static inline Simd4Double
simd4CmpEqD(Simd4Double a, Simd4Double b) { return simdCmpEqD(a, b); }

/*! \brief Return true for nonzero SIMD4 elements (i.e., with bits set)
 * \copydetails simdCmpNzD
 */
static inline Simd4Double
simd4CmpNzD(Simd4Double a) { return simdCmpNzD(a); }

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails simdCmpLtD
 */
static inline Simd4Double
simd4CmpLtD(Simd4Double a, Simd4Double b) { return simdCmpLtD(a, b); }

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails simdCmpLeD
 */
static inline Simd4Double
simd4CmpLeD(Simd4Double a, Simd4Double b) { return simdCmpLeD(a, b); }

/*! \brief Logical AND on double SIMD4 booleans.
 * \copydetails simdAndDB
 */
static inline Simd4Double
simd4AndFB(Simd4FDool a, Simd4DBool b) { return simdAndFB(a, b); }

/*! \brief Logical OR on double SIMD4 booleans.
 * \copydetails simdOrDB
 */
static inline Simd4Double
simd4OrFB(Simd4DBool a, Simd4DBool b) { return simdOrFB(a, b); }

/*! \brief Returns non-zero if any of the SIMD4 booleans in x is True.
 * \copydetails simdAnyTrueDB
 */
static inline int
simd4AnyTrueFB(Simd4DBool a) { return simdAnyTrueFB(a); }

/*! \brief Select from double precision SIMD4 variable where boolean is true.
 * \copydetails simdMaskD
 */
static inline Simd4Double
simd4MaskD(Simd4Double a, Simd4DBool mask) { return simdMaskD(a, mask); }

/*! \brief Select from double precision SIMD4 variable where boolean is false.
 * \copydetails simdMaskNotD
 */
static inline Simd4Double
simd4MaskNotD(Simd4Double a, Simd4DBool mask) { return simdMaskNotD(a, mask); }

/*! \brief Vector-blend instruction for SIMD4 double.
 * \copydetails simdBlendD
 */
static inline Simd4Double
simd4BlendD(Simd4Double a, Simd4Double b, Simd4DBool sel) { return simdBlendD(a, b, sel); }

/*! \brief Return sum of all elements in SIMD4 double.
 * \copydetails simdReduceD
 */
static inline double
simd4ReduceD(Simd4Double a) { return simdReduceD(a); }

#endif  /* GMX_SIMD_DOUBLE_WIDTH == 4 */

/*! \} */


/*! \} */
/*! \endcond */

}      // namespace gmx

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H */
