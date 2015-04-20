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

#include <algorithm>

#include "impl_reference_definitions.h"
#include "impl_reference_simd_float.h"

namespace gmx
{

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Constant width-4 single precision SIMD types and instructions
 * \{
 */

/*! \libinternal \brief SIMD4 float type. Available if \ref GMX_SIMD4_HAVE_FLOAT.
 *
 * Unless you specifically want a single-precision type you should check
 * \ref gmx::Simd4Real instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
struct Simd4Float
{
    float r[GMX_SIMD4_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \brief Load SIMD4 float from aligned memory.
 *  \copydetails simdLoadF
 */
static inline Simd4Float
simd4LoadF(const float *m)
{
    Simd4Float        a;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Set all elements of SIMD4 float from single pointer.
 *  \copydetails simdLoad1F
 */
static inline Simd4Float
simd4Load1F(const float *m)
{
    Simd4Float        a;
    float             f = *m;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = f;
    }
    return a;
}


/*! \brief Set all SIMD4 float elements to the value r.
 *  \copydetails simdSet1F
 */
static inline Simd4Float
simd4Set1F(float r)
{
    Simd4Float        a;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}


/*! \brief Store the contents of SIMD4 float pr to aligned memory m.
 *  \copydetails simdStoreF
 */
static inline void
simd4StoreF(float *m, Simd4Float a)
{
    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD4 float from unaligned memory.
 * \copydetails simdLoadUF
 */
static inline Simd4Float
simd4LoadUF(const float *m) { return simd4LoadF(m); }

/*! \brief Store SIMD4 float to unaligned memory.
 * \copydetails simdStoreUF
 */
static inline void
simd4StoreUF(float *m, Simd4Float a) { simd4StoreF(m, a); }

/*! \brief Set all SIMD4 float elements to 0.
 * \copydetails simdSetZeroF
 */
static inline Simd4Float
simd4SetZeroF()
{
    Simd4Float a;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        a.r[i] = 0.0f;
    }
    return a;
}


/*! \brief Bitwise and for two SIMD4 float variables.
 * \copydetails simdAndF
 */
static inline Simd4Float
simd4AndF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}


/*! \brief Bitwise andnot for two SIMD4 float variables. c=(~a) & b.
 * \copydetails simdAndNotF
 */
static inline Simd4Float
simd4AndNotF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = (~conv1.i) & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}


/*! \brief Bitwise or for two SIMD4 float variables.
 * \copydetails simdOrF
 */
static inline Simd4Float
simd4OrF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i | conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise xor for two SIMD4 float variables.
 * \copydetails simdXorF
 */
static inline Simd4Float
simd4XorF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i ^ conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Add two SIMD4 float variables.
 * \copydetails simdAddF
 */
static inline Simd4Float
simd4AddF(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }

    return c;
}

/*! \brief Subtract two SIMD4 float variables.
 * \copydetails simdSubF
 */
static inline Simd4Float
simd4SubF(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }

    return c;
}

/*! \brief Multiply two SIMD4 float variables.
 * \copydetails simdMulF
 */
static inline Simd4Float
simd4MulF(Simd4Float a, Simd4Float b)
{
    Simd4Float c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = a.r[i] * b.r[i];
    }

    return c;
}

/*! \brief Fused-multiply-add for SIMD4 float. Result is a*b+c.
 * \copydetails simdFmaddF
 */
static inline Simd4Float
simd4FmaddF(Simd4Float a, Simd4Float b, Simd4Float c) { return simd4AddF(simd4MulF(a, b), c); }

/*! \brief Fused-multiply-subtract for SIMD4 float. Result is a*b-c.
 * \copydetails simdFmsubF
 */
static inline Simd4Float
simd4FmsubF(Simd4Float a, Simd4Float b, Simd4Float c) { return simd4SubF(simd4MulF(a, b), c); }

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b+c.
 * \copydetails simdFnmaddF
 */
static inline Simd4Float
simd4FnmaddF(Simd4Float a, Simd4Float b, Simd4Float c) { return simd4SubF(c, simd4MulF(a, b)); }

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b-c.
 * \copydetails simdFnmsubF
 */
static inline Simd4Float
simd4FnmsubF(Simd4Float a, Simd4Float b, Simd4Float c)
{
    return simd4SubF(simd4SetZeroF(), simd4FmaddF(a, b, c));
}

/*! \brief Lookup of approximate 1/sqrt(x) for SIMD4 float.
 * \copydetails simdRsqrtF
 */
static inline Simd4Float
simd4RsqrtF(Simd4Float x)
{
    Simd4Float        b;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        b.r[i] = 1.0f / std::sqrt(x.r[i]);
    }
    return b;
};


/*! \brief Floating-point absolute value for SIMD4 float.
 * \copydetails simdAbsF
 */
static inline Simd4Float
simd4AbsF(Simd4Float a)
{
    Simd4Float        c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = std::abs(a.r[i]);
    }
    return c;
}

/*! \brief Floating-point negate for SIMD4 float.
 * \copydetails simdNegF
 */
static inline Simd4Float
simd4NegF(Simd4Float a)
{
    Simd4Float        c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/*! \brief Set each SIMD4 float element to the largest from two variables.
 * \copydetails simdMaxF
 */
static inline Simd4Float
simd4MaxF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = std::max(a.r[i], b.r[i]);
    }
    return c;
}


/*! \brief Set each SIMD4 float element to the smallest from two variables.
 * \copydetails simdMinF
 */
static inline Simd4Float
simd4MinF(Simd4Float a, Simd4Float b)
{
    Simd4Float        c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = std::min(a.r[i], b.r[i]);
    }
    return c;
}


/*! \brief Round to nearest integer value for SIMD4 float.
 * \copydetails simdRoundF
 */
static inline Simd4Float
simd4RoundF(Simd4Float a)
{
    Simd4Float        b;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        b.r[i] = std::round(a.r[i]);
    }
    return b;
}


/*! \brief Round to largest integral value for SIMD4 float.
 * \copydetails simdTruncF
 */
static inline Simd4Float
simd4TruncF(Simd4Float a)
{
    Simd4Float        b;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        b.r[i] = std::trunc(a.r[i]);
    }
    return b;
}

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
simd4DotProductF(Simd4Float a, Simd4Float b)
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
simd4Transpose(Simd4Float * v0, Simd4Float * v1,
               Simd4Float * v2, Simd4Float * v3)
{
    Simd4Float t0 = *v0;
    Simd4Float t1 = *v1;
    Simd4Float t2 = *v2;
    Simd4Float t3 = *v3;
    v0->r[0] = t0.r[0];
    v0->r[1] = t1.r[0];
    v0->r[2] = t2.r[0];
    v0->r[3] = t3.r[0];
    v1->r[0] = t0.r[1];
    v1->r[1] = t1.r[1];
    v1->r[2] = t2.r[1];
    v1->r[3] = t3.r[1];
    v2->r[0] = t0.r[2];
    v2->r[1] = t1.r[2];
    v2->r[2] = t2.r[2];
    v2->r[3] = t3.r[2];
    v3->r[0] = t0.r[3];
    v3->r[1] = t1.r[3];
    v3->r[2] = t2.r[3];
    v3->r[3] = t3.r[3];
}

/*! \libinternal  \brief SIMD4 variable type to use for logical comparisons on floats.
 * \copydetails SimdFBool
 */
struct Simd4FBool
{
    std::int32_t b[GMX_SIMD4_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \brief Equality comparison of two single precision SIMD4.
 * \copydetails simdCmpEqF
 */
static inline Simd4FBool
simd4CmpEqF(Simd4Float a, Simd4Float b)
{
    Simd4FBool        c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLtF
 */
static inline Simd4FBool
simd4CmpLtF(Simd4Float a, Simd4Float b)
{
    Simd4FBool         c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}


/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails simdCmpLeF
 */
static inline Simd4FBool
simd4CmpLeF(Simd4Float a, Simd4Float b)
{
    Simd4FBool         c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}

/*! \brief Logical AND on float SIMD4 booleans.
 * \copydetails simdAndFB
 */
static inline Simd4FBool
simd4AndFB(Simd4FBool a, Simd4FBool b)
{
    Simd4FBool c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.b[i] = a.b[i] && b.b[i];
    }

    return c;
}

/*! \brief Logical OR on float SIMD4 booleans.
 * \copydetails simdOrFB
 */
static inline Simd4FBool
simd4OrFB(Simd4FBool a, Simd4FBool b)
{
    Simd4FBool c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.b[i] = a.b[i] || b.b[i];
    }

    return c;
}

/*! \brief Returns true if any of the SIMD4 boolean in x is True.
 * \copydetails simdAnyTrueFB
 */
static inline bool
simd4AnyTrueFB(Simd4FBool a)
{
    return a.b[0] || a.b[1] || a.b[2] || a.b[3];
}

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 * \copydetails simdMaskF
 */
static inline Simd4Float
simd4MaskF(Simd4Float a, Simd4FBool mask)
{
    Simd4Float         c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? a.r[i] : 0.0f;
    }
    return c;
}

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 * \copydetails simdMaskNotF
 */
static inline Simd4Float
simd4MaskNotF(Simd4Float a, Simd4FBool mask)
{
    Simd4Float         c;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? 0.0f : a.r[i];
    }
    return c;
}


/*! \brief Vector-blend instruction form SIMD4 float.
 * \copydetails simdBlendF
 */
static inline Simd4Float
simd4BlendF(Simd4Float a, Simd4Float b, Simd4FBool sel)
{
    Simd4Float        d;

    for (int i = 0; i < GMX_SIMD4_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}


/*! \brief Return sum of all elements in SIMD4 float.
 * \copydetails simdReduceF
 */
static inline float
simd4ReduceF(Simd4Float a)
{
    return a.r[0]+a.r[1]+a.r[2]+a.r[3];
}

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_SIMD4_FLOAT_H
