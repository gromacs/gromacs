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

#include <cmath>

#include "impl_reference_common.h"
#include "impl_reference_simd_float.h"

namespace gmx
{

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
typedef SimdFloat         Simd4Float;

/*! \brief Load SIMD4 float from aligned memory.
 *  \copydetails simdLoadF
 */
static inline Simd4Float
simd4LoadF(const float *m) { return simdLoadF(m); }

/*! \brief Set all elements of SIMD4 float from single pointer.
 *  \copydetails simdLoad1F
 */
static inline Simd4Float
simd4Load1F(const float *m) { return simdLoad1F(m); }

/*! \brief Set all SIMD4 float elements to the value r.
 *  \copydetails simdSet1F
 */
static inline Simd4Float
simd4Set1F(float r) { return simdSet1F(r); }

/*! \brief Store the contents of SIMD4 float pr to aligned memory m.
 *  \copydetails simdStoreF
 */
static inline void
simd4StoreF(float *m, Simd4Float a) { simdStoreF(m, a); }

/*! \brief Load SIMD4 float from unaligned memory.
 * \copydetails simdLoadUF
 */
static inline Simd4Float
simd4LoadUF(const float *m) { return simdLoadF(m); }

/*! \brief Store SIMD4 float to unaligned memory.
 * \copydetails simdStoreUF
 */
static inline void
simd4StoreUF(float *m, Simd4Float a) { simdStoreF(m, a); }

/*! \brief Set all SIMD4 float elements to 0.
 * \copydetails simdSetZeroF
 */
static inline Simd4Float
simd4SetZeroF() { return simdSetZeroF(); }

/*! \brief Bitwise and for two SIMD4 float variables.
 * \copydetails simdAndF
 */
static inline Simd4Float
simd4AndF(Simd4Float a, Simd4Float b) { return simdAndF(a, b); }

/*! \brief Bitwise andnot for two SIMD4 float variables. c=(~a) & b.
 * \copydetails simdAndNotF
 */
static inline Simd4Float
simd4AndNotF(Simd4Float a, Simd4Float b) { return simdAndNotF(a, b); }

/*! \brief Bitwise or for two SIMD4 float variables.
 * \copydetails simdOrF
 */
static inline Simd4Float
simd4OrF(Simd4Float a, Simd4Float b) { return simdOrF(a, b); }

/*! \brief Bitwise xor for two SIMD4 float variables.
 * \copydetails simdXorF
 */
static inline Simd4Float
simd4XorF(Simd4Float a, Simd4Float b) { return simdXorF(a, b); }

/*! \brief Add two SIMD4 float variables.
 * \copydetails simdAddF
 */
static inline Simd4Float
simd4AddF(Simd4Float a, Simd4Float b) { return simdAddF(a, b); }

/*! \brief Subtract two SIMD4 float variables.
 * \copydetails simdSubF
 */
static inline Simd4Float
simd4SubF(Simd4Float a, Simd4Float b) { return simdSubF(a, b); }

/*! \brief Multiply two SIMD4 float variables.
 * \copydetails simdMulF
 */
static inline Simd4Float
simd4MulF(Simd4Float a, Simd4Float b) { return simdMulF(a, b); }

/*! \brief Fused-multiply-add for SIMD4 float. Result is a*b+c.
 * \copydetails simdFmaddF
 */
static inline Simd4Float
simd4FmaddF(Simd4Float a, Simd4Float b, Simd4Float c) { return simdFmaddF(a, b, c); }

/*! \brief Fused-multiply-subtract for SIMD4 float. Result is a*b-c.
 * \copydetails simdFmsubF
 */
static inline Simd4Float
simd4FmsubF(Simd4Float a, Simd4Float b, Simd4Float c) { return simdFmsubF(a, b, c); }

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b+c.
 * \copydetails simdFnmaddF
 */
static inline Simd4Float
simd4FnmaddF(Simd4Float a, Simd4Float b, Simd4Float c) { return simdFnmaddF(a, b, c); }

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b-c.
 * \copydetails simdFnmsubF
 */
static inline Simd4Float
simd4FnmsubF(Simd4Float a, Simd4Float b, Simd4Float c) { return simdFnmsubF(a, b, c); }

/*! \brief Lookup of approximate 1/sqrt(x) for SIMD4 float.
 * \copydetails simdRsqrtF
 */
static inline Simd4Float
simd4RsqrtF(Simd4Float x) { return simdRsqrtF(x); }

/*! \brief Floating-point absolute value for SIMD4 float.
 * \copydetails simdAbsF
 */
static inline Simd4Float
simd4AbsF(Simd4Float a) { return simdAbsF(a); }

/*! \brief Floating-point negate for SIMD4 float.
 * \copydetails simdNegF
 */
static inline Simd4Float
simd4NegF(Simd4Float a) { return simdNegF(a); }

/*! \brief Set each SIMD4 float element to the largest from two variables.
 * \copydetails simdMaxF
 */
static inline Simd4Float
simd4MaxF(Simd4Float a, Simd4Float b) { return simdMaxF(a, b); }

/*! \brief Set each SIMD4 float element to the smallest from two variables.
 * \copydetails simdMinF
 */
static inline Simd4Float
simd4MinF(Simd4Float a, Simd4Float b) { return simdMinF(a, b); }

/*! \brief Round to nearest integer value for SIMD4 float.
 * \copydetails simdRoundF
 */
static inline Simd4Float
simd4RoundF(Simd4Float a) { return simdRoundF(a); }

/*! \brief Round to largest integral value for SIMD4 float.
 * \copydetails simdTruncF
 */
static inline Simd4Float
simd4TruncF(Simd4Float a) { return simdTruncF(a); }

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
    return a.r[0] * b.r[0] + a.r[1] * b.r[1] + a.r[2] * b.r[2];
}

/*! \brief SIMD4 float transpose
 *
 * \param[in,out] v0  Row 0 on input, column 0 on output
 * \param[in,out] v1  Row 1 on input, column 1 on output
 * \param[in,out] v2  Row 2 on input, column 2 on output
 * \param[in,out] v3  Row 3 on input, column 3 on output
 *
 * This is only available in C++.
 */
static inline void
simd4Transpose(Simd4Float &v0, Simd4Float &v1,
               Simd4Float &v2, Simd4Float &v3)
{
    Simd4Float t0, t1, t2, t3;
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

/*! \brief SIMD4 variable type to use for logical comparisons on floats.
 * \copydetails SimdFBool
 */
typedef SimdFBool        Simd4FBool;

/*! \brief Equality comparison of two single precision SIMD4.
 * \copydetails simdCmpEqF
 */
static inline Simd4FBool
simd4CmpEqF(Simd4Float a, Simd4Float b) { return simdCmpEqF(a, b); }

/*! \brief Return true for nonzero SIMD4 elements (i.e., with bits set)
 * \copydetails simdCmpNzF
 */
static inline Simd4FBool
simd4CmpNzF(Simd4Float a) { return simdCmpNzF(a); }

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLtF
 */
static inline Simd4FBool
simd4CmpLtF(Simd4Float a, Simd4Float b) { return simdCmpLtF(a, b); }

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLeF
 */
static inline Simd4FBool
simd4CmpLeF(Simd4Float a, Simd4Float b) { return simdCmpLeF(a, b); }

/*! \brief Logical AND on float SIMD4 booleans.
 * \copydetails simdAndFB
 */
static inline Simd4FBool
simd4AndFB(Simd4FBool a, Simd4FBool b) { return simdAndFB(a, b); }

/*! \brief Logical OR on float SIMD4 booleans.
 * \copydetails simdOrFB
 */
static inline Simd4FBool
simd4OrFB(Simd4FBool a, Simd4FBool b) { return simdOrFB(a, b); }

/*! \brief Returns non-zero if any of the SIMD4 boolean in x is True.
 * \copydetails simdAnyTrueFB
 */
static inline int
simd4AnyTrueFB(Simd4FBool a) { return simdAnyTrueFB(a); }

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 * \copydetails simdMaskF
 */
static inline Simd4Float
simd4MaskF(Simd4Float a, Simd4FBool mask) { return simdMaskF(a, mask); }

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 * \copydetails simdMaskNotF
 */
static inline Simd4Float
simd4MaskNotF(Simd4Float a, Simd4FBool mask) { return simdMaskNotF(a, mask); }

/*! \brief Vector-blend instruction form SIMD4 float.
 * \copydetails simdBlendF
 */
static inline Simd4Float
simd4BlendF(Simd4Float a, Simd4Float b, Simd4FBool sel) { return simdBlendF(a, b, sel); }

/*! \brief Return sum of all elements in SIMD4 float.
 * \copydetails simdReduceF
 */
static inline float
simd4ReduceF(Simd4Float a) { return simdReduceF(a); }

#endif  /* GMX_SIMD_FLOAT_WIDTH==4 */

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif /* GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H */
