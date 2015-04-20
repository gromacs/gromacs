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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD single precision.

 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <algorithm>

#include "impl_reference_definitions.h"

namespace gmx
{

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name SIMD implementation data types
 * \{
 */
/*! \libinternal \brief Float SIMD variable. Supported if GMX_SIMD_HAVE_FLOAT.
 */
struct SimdFloat
{
    float r[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from float.
 *
 * This is also the widest integer SIMD type. Available if GMX_SIMD_HAVE_FLOAT.
 *
 * \note The integer SIMD type will always be available, but on architectures
 * that do not have any real integer SIMD support it might be defined as the
 * floating-point type. This will work fine, since there are separate defines
 * for whether the implementation can actually do any operations on integer
 * SIMD types.
 */
struct SimdFInt32
{
    std::int32_t i[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \libinternal \brief Boolean type for float SIMD data.
 *
 * You should likely use SimdBool
 * (for SimdReal) instead, unless you really know what you are doing.
 */
struct SimdFBool
{
    std::int32_t b[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \libinternal \brief Boolean type for integer datatypes corresponding to float SIMD. */
struct SimdFIBool
{
    std::int32_t b[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
};

/*! \}
 *
 * \name SIMD implementation load/store operations for single precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_FLOAT_WIDTH numbers from aligned memory.
 *
 * \param m Pointer to memory aligned to the SIMD width.
 * \return SIMD variable with data loaded.
 */
static inline SimdFloat
simdLoadF(const float *m)
{
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH*sizeof(float)) == 0);

    SimdFloat         a;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Set all SIMD variable elements to float pointed to by m (unaligned).
 *
 * \param m Pointer to single value in memory.
 * \return SIMD variable with all elements set to *m.
 */
static inline SimdFloat
simdLoad1F(const float *m)
{
    SimdFloat         a;
    float             f = *m;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = f;
    }
    return a;
}

/*! \brief Set all SIMD float variable elements to the value r.
 *
 *  \param r floating-point constant
 *  \return SIMD variable with all elements set to r.
 */
static inline SimdFloat
simdSet1F(float r)
{
    SimdFloat         a;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}

/*! \brief Set all SIMD float variable elements to 0.0f.
 *
 *  \return The value 0.0 in all elements of a SIMD variable.
 */
static inline SimdFloat
simdSetZeroF()
{
    SimdFloat         a;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = 0.0f;
    }
    return a;
}

/*! \brief Store the contents of the SIMD float variable pr to aligned memory m.
 *
 * \param[out] m Pointer to memory, aligned to SIMD width.
 * \param a SIMD variable to store
 */
static inline void
simdStoreF(float *m, SimdFloat a)
{
    assert(std::size_t(m) % (GMX_SIMD_FLOAT_WIDTH*sizeof(float)) == 0);

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD float from unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirement.
 * \return SIMD variable with data loaded.
 */
static inline SimdFloat
simdLoadUF(const float *m)
{
    SimdFloat         a;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Store SIMD float to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a SIMD variable to store.
 */
static inline void
simdStoreUF(float *m, SimdFloat a)
{
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \} */


/*!
 * \name SIMD implementation load/store operations for integers (corresponding to float)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdLoadI.
 *
 * \param m Pointer to memory, aligned to integer SIMD width.
 * \return SIMD integer variable.
 */
static inline SimdFInt32
simdLoadFI(const std::int32_t * m)
{
    SimdFInt32         a;

    assert(std::size_t(m) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdSet1I.
 *
 *  \param b integer value to set variable to.
 *  \return SIMD variable with all elements set to b.
 */
static inline SimdFInt32
simdSet1FI(std::int32_t b)
{
    SimdFInt32         a;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdSetZeroI.
 *
 * \return SIMD integer variable with all bits set to zero.
 */
static inline SimdFInt32
simdSetZeroFI()
{
    SimdFInt32         a;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdStoreI.
 *
 * \param m Memory aligned to integer SIMD width.
 * \param a SIMD variable to store.
 */
static inline void
simdStoreFI(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % (GMX_SIMD_FINT32_WIDTH*sizeof(std::int32_t)) == 0);

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdLoadUI.
 *
 * Supported with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirements.
 * \return SIMD integer variable.
 */
static inline SimdFInt32
simdLoadUFI(const std::int32_t *m)
{
    SimdFInt32         a;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
}

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically call the real-precision \ref gmx::simdStoreUI.
 *
 * Supported with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param m Memory pointer, no alignment requirements.
 * \param a SIMD variable to store.
 */
static inline void
simdStoreUFI(std::int32_t * m, SimdFInt32 a)
{
    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
}

/*! \brief Extract element with index i from \ref gmx::SimdFInt32.
 *
 * You should typically call the real-precision \ref gmx::simdExtractI.
 *
 * Available with \ref GMX_SIMD_HAVE_FINT32_EXTRACT.
 *
 * \tparam index Compile-time constant, position to extract (first position is 0)
 * \param  a     SIMD variable from which to extract value.
 * \return Single integer from position index in SIMD variable.
 */
template<int index>
static inline std::int32_t
simdExtractFI(SimdFInt32 a)
{
    return a.i[index];
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point bitwise logical operations
 * \{
 */

/*! \brief Bitwise and for two SIMD float variables. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx::simdAnd.
 *
 * \param a data1
 * \param b data2
 * \return data1 & data2
 */
static inline SimdFloat
simdAndF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise andnot for SIMD float. c=(~a) & b. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx::simdAndNot.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 */
static inline SimdFloat
simdAndNotF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = (~conv1.i) & conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise or for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx::simdOr.
 *
 * \param a data1
 * \param b data2
 * \return data1 | data2
 */
static inline SimdFloat
simdOrF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i | conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \brief Bitwise xor for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx::simdXor.
 *
 * \param a data1
 * \param b data2
 * \return data1 ^ data2
 */
static inline SimdFloat
simdXorF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i ^ conv2.i;
        c.r[i]  = conv1.r;
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point arithmetics
 * \{
 */
/*! \brief Add two float SIMD variables.
 *
 * You should typically call the real-precision \ref gmx::simdAdd.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static inline SimdFloat
simdAddF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/*! \brief Subtract two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx::simdSub.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static inline SimdFloat
simdSubF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/*! \brief Multiply two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx::simdMul.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 */
static inline SimdFloat
simdMulF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] * b.r[i];
    }
    return c;
}

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * You should typically call the real-precision \ref gmx::simdFmadd.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static inline SimdFloat
simdFmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return simdAddF(simdMulF(a, b), c);
}



/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * You should typically call the real-precision \ref gmx::simdFmsub.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static inline SimdFloat
simdFmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return simdSubF(simdMulF(a, b), c);
}


/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * You should typically call the real-precision \ref gmx::simdFnmadd.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static inline SimdFloat
simdFnmaddF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return simdSubF(c, simdMulF(a, b));
}


/*! \brief Fused-negated-multiply-sub. Result is -a*b-c.
 *
 * You should typically call the real-precision \ref gmx::simdFnmsub.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is 1 this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static inline SimdFloat
simdFnmsubF(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return simdSubF(simdSetZeroF(), simdFmaddF(a, b, c));
}


/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * You should typically call the real-precision \ref gmx::simdRsqrt.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 */
static inline SimdFloat
simdRsqrtF(SimdFloat x)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = 1.0f / std::sqrt(x.r[i]);
    }
    return b;
};

/*! \brief SIMD 1.0/x lookup.
 *
 * You should typically call the real-precision \ref gmx::simdRcp.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x!=0
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 */
static inline SimdFloat
simdRcpF(SimdFloat x)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = 1.0f / x.r[i];
    }
    return b;
};

/*! \brief Multiply two SIMD variables, masked version.
 *
 * You should typically call the real-precision simdMulMask().
 *
 * \param a factor1
 * \param b factor2
 * \param m mask
 * \return a*b where mask is true, 0.0 otherwise.
 */
static inline SimdFloat
simdMulMaskF(SimdFloat a, SimdFloat b, SimdFBool m)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = m.b[i] ? (a.r[i] * b.r[i]) : 0.0f;
    }
    return c;
}

/*! \brief Fused-multiply-add. Result is a*b+c, masked version.
 *
 * You should typically call the real-precision simdFmaddMask().
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \param m mask
 * \return a*b+c where mask is true, 0.0 otherwise.
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static inline SimdFloat
simdFmaddMaskF(SimdFloat a, SimdFloat b, SimdFloat c,
               SimdFBool m)
{
    SimdFloat         d;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d.r[i] = m.b[i] ? (a.r[i] * b.r[i] + c.r[i]) : 0.0f;
    }
    return d;
}

/*! \brief SIMD 1.0/sqrt(x) lookup, masked version.
 *
 * You should typically call the real-precision \ref gmx::simdRsqrtMask().
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0 for entries where mask is true.
 * \param m Mask
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 *         The result for masked-out entries will be 0.0.
 */
static inline SimdFloat
simdRsqrtMaskF(SimdFloat x, SimdFBool m)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (m.b[i] != 0) ? 1.0f / std::sqrt(x.r[i]) : 0.0f;
    }
    return b;
}

/*! \brief SIMD 1.0/x lookup, masked version.
 *
 * You should typically call the real-precision \ref gmx::simdRcpMask().
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x>0 for entries where mask is true.
 * \param m Mask
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 *         The result for masked-out entries will be 0.0.
 */
static inline SimdFloat
simdRcpMaskF(SimdFloat x, SimdFBool m)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (m.b[i] != 0) ? 1.0f / x.r[i] : 0.0f;
    }
    return b;
}

/*! \brief SIMD Floating-point fabs().
 *
 * You should typically call the real-precision \ref gmx::simdAbs.
 *
 * \param a any floating point values
 * \return fabs(a) for each element.
 */
static inline SimdFloat
simdAbsF(SimdFloat a)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = std::abs(a.r[i]);
    }
    return c;
}

/*! \brief SIMD floating-point negate.
 *
 * You should typically call the real-precision \ref gmx::simdNeg.
 *
 * \param a Any floating-point value
 * \return -a
 */
static inline SimdFloat
simdNegF(SimdFloat a)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * You should typically call the real-precision \ref gmx::simdMax.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
static inline SimdFloat
simdMaxF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = std::max(a.r[i], b.r[i]);
    }
    return c;
}

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * You should typically call the real-precision \ref gmx::simdMin.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return min(a,b) for each element.
 */
static inline SimdFloat
simdMinF(SimdFloat a, SimdFloat b)
{
    SimdFloat         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = std::min(a.r[i], b.r[i]);
    }
    return c;
}

/*! \brief Round to nearest integer value (in floating-point format).
 *
 * You should typically call the real-precision \ref gmx::simdRound.
 *
 * \param a Any floating-point value
 * \return The nearest integer, represented in floating-point format.
 *
 * \note The reference implementation rounds exact half-way cases
 * away from zero, whereas most SIMD intrinsics will round to nearest even.
 * This could be fixed by using rint/rintf, but the bigger problem is that
 * MSVC does not support full C99, and none of the round or rint
 * functions are defined. It's much easier to approximately implement
 * round() than rint(), so we do that and hope we never get bitten in
 * testing. (Thanks, Microsoft.)
 */
static inline SimdFloat
simdRoundF(SimdFloat a)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = std::round(a.r[i]);
    }
    return b;
}

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * You should typically call the real-precision \ref gmx::simdTrunc.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This is truncation towards zero, not floor(). The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but floor() frequently isn't.
 */
static inline SimdFloat
simdTruncF(SimdFloat a)
{
    SimdFloat         b;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = std::trunc(a.r[i]);
    }
    return b;
}


/*! \brief Fraction of the SIMD floating point number.
 *
 * You should typically call the real-precision \ref gmx::simdFraction.
 *
 * \param a Any floating-point value
 * \return a-trunc(r)
 *
 * To maximize compatibility, we use the same definition of fractions as used
 * e.g. for the AMD64 hardware instructions. This relies on truncation towards
 * zero for the integer part, and the remaining fraction can thus be either
 * positive or negative. As an example, -1.42 would return the fraction -0.42.
 *
 * Hardware support with \ref GMX_SIMD_HAVE_FRACTION, otherwise emulated.
 */
static inline SimdFloat
simdFractionF(SimdFloat a)
{
    return simdSubF(a, simdTruncF(a));
}

/*! \brief Extract (integer) exponent from single precision SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdGetExponent.
 *
 * \param a Any floating-point value
 * \return Exponent value, represented in floating-point format.
 *
 * The IEEE754 exponent field is selected, the bias removed, and it is converted to
 * a normal floating-point SIMD.
 */
static inline SimdFloat
simdGetExponentF(SimdFloat a)
{
    /* Mask with ones for the exponent field of single precision fp */
    const std::int32_t  exponentMask = 0x7f800000;
    SimdFloat           b;

    union
    {
        float         f;
        std::int32_t  i;
    }
    conv;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* Keep exponent, shift 23 right (float mantissa), remove bias (127) */
        b.r[i] = ((conv.i & exponentMask) >> 23) - 127;
    }
    return b;
}

/*! \brief Get SIMD mantissa.
 *
 * You should typically call the real-precision \ref gmx::simdGetMantissa.
 *
 * \param a Any floating-point value
 * \return Mantissa, represented in floating-point format.
 *
 * The mantissa field is selected, and a new neutral exponent created.
 */
static inline SimdFloat
simdGetMantissaF(SimdFloat a)
{
    const std::int32_t  mantissaMask = 0x007fffff;
    const std::int32_t  one          = 0x3f800000;
    SimdFloat           b;

    union
    {
        float         f;
        std::int32_t  i;
    }
    conv;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* remove current exponent, add a biased exponent for 1.0 (i.e., 2^0=1) */
        conv.i = (conv.i & (mantissaMask)) | one;
        b.r[i] = conv.f;
    }
    return b;
}

/*! \brief Set (integer) exponent from single precision floating-point SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdSetExponent.
 *
 * \param a A floating point value that will not overflow as 2^a.
 * \return 2^(round(a)).
 *
 * The input is \a rounded to the nearest integer, the exponent bias is added
 * to this integer, and the bits are shifted to the IEEE754 exponent part of the number.
 *
 * \note The argument will be \a rounded to nearest integer since that is what
 * we need for the exponential functions, and this integer x will be set as the
 * exponent so the new floating-point number will be 2^x.
 */
static inline SimdFloat
simdSetExponentF(SimdFloat a)
{
    SimdFloat           b;
    std::int32_t        iExponent;

    union
    {
        float         f;
        std::int32_t  i;
    }
    conv;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        /* Critical to use same algorithm as for simdRoundF() */
        iExponent = static_cast<std::int32_t>(std::round(a.r[i]));
        /* Add bias (127), and shift 23 bits left (mantissa size) */
        conv.i = (iExponent + 127) << 23;
        b.r[i] = conv.f;
    }
    return b;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point comparisons, boolean, selection.
 * \{
 */
/*! \brief SIMD a==b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdCmpEq.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a==b.
 *
 * Beware that exact floating-point comparisons are difficult.
 */
static inline SimdFBool
simdCmpEqF(SimdFloat a, SimdFloat b)
{
    SimdFBool         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/*! \brief SIMD a!=0 for single SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdCmpNz().
 *
 * \param a value
 * \return Each element of the boolean will be true if any bit in a is nonzero.
 *         The behaviour for negative zero is undefined, and should not be
 *         relied on - it will depend on the architecture.
 */
static inline SimdFBool
simdCmpNzF(SimdFloat a)
{
    SimdFBool         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] != 0.0f);
    }
    return c;
}

/*! \brief SIMD a<b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdCmpLt.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<b.
 */
static inline SimdFBool
simdCmpLtF(SimdFloat a, SimdFloat b)
{
    SimdFBool          c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/*! \brief SIMD a<=b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx::simdCmpLe.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<=b.
 */
static inline SimdFBool
simdCmpLeF(SimdFloat a, SimdFloat b)
{
    SimdFBool          c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}

/*! \brief Logical \a and on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx::simdAnd.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a \& b are true.
 *
 * \note This is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa simdAndIB
 */
static inline SimdFBool
simdAndFB(SimdFBool a, SimdFBool b)
{
    SimdFBool         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical \a or on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx::simdOr.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a or b is true.
 *
 * Note that this is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa simdOrIB
 */
static inline SimdFBool
simdOrFB(SimdFBool a, SimdFBool b)
{
    SimdFBool         c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx::simdAnyTrueB.
 *
 * \param a Logical variable.
 * \return true if any element in a is true, otherwise false.
 *
 * The actual return value for truth will depend on the architecture,
 * so any non-zero value is considered truth.
 */
static inline bool
simdAnyTrueFB(SimdFBool a)
{
    bool anyTrue = false;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        anyTrue = anyTrue || a.b[i];
    }
    return anyTrue;
}

/*! \brief Select from single precision SIMD variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx::simdMask.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for true, 0 for false.
 */
static inline SimdFloat
simdMaskF(SimdFloat a, SimdFBool mask)
{
    SimdFloat          c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? a.r[i] : 0.0f;
    }
    return c;
}

/*! \brief Select from single precision SIMD variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx::simdMaskNot.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for false, 0 for true (sic).
 */
static inline SimdFloat
simdMaskNotF(SimdFloat a, SimdFBool mask)
{
    SimdFloat          c;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = mask.b[i] ? 0.0f : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx::simdBlend.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static inline SimdFloat
simdBlendF(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    SimdFloat         d;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}

/*! \brief Return sum of all elements in SIMD float variable.
 *
 * You should typically call the real-precision \ref gmx::simdReduce.
 *
 * \param a SIMD variable to reduce/sum.
 * \return The sum of all elements in the argument variable.
 *
 */
static inline float
simdReduceF(SimdFloat a)
{
    float     sum = 0.0f;

    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        sum += a.r[i];
    }
    return sum;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) bitwise logical operations
 * \{
 */

/*! \brief SIMD integer shift left logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx::simdSlliI.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  left, while zeros are shifted in from the right. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is 1.
 *
 * \param a integer data to shift
 * \param n number of positions to shift left. n<=32.
 * \return shifted values
 */
static inline SimdFInt32
simdSlliFI(SimdFInt32 a, int n)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/*! \brief SIMD integer shift right logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx::simdSrliI.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  right, while zeros are shifted in from the left. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is 1.
 *
 * \param a integer data to shift
 * \param n number of positions to shift right. n<=32.
 * \return shifted values
 */
static inline SimdFInt32
simdSrliFI(SimdFInt32 a, int n)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/*! \brief Integer SIMD bitwise and.
 *
 * You should typically call the real-precision \ref gmx::simdAndI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \note You can \a not use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx::simdMaskI instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \& b (bitwise and)
 */
static inline SimdFInt32
simdAndFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise not-and.
 *
 * You should typically call the real-precision \ref gmx::simdAndNotI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * Note that you can NOT use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx::simdMaskNotI instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return (~a) \& b (bitwise andnot)
 */
static inline SimdFInt32
simdAndNotFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise or.
 *
 * You should typically call the real-precision \ref gmx::simdOrI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \| b (bitwise or)
 */
static inline SimdFInt32
simdOrFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise xor.
 *
 * You should typically call the real-precision \ref gmx::simdXorI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a ^ b (bitwise xor)
 */
static inline SimdFInt32
simdXorFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] ^ b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) arithmetics
 * \{
 */
/*! \brief Add SIMD integers.
 *
 * You should typically call the real-precision \ref gmx::simdXorI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static inline SimdFInt32
simdAddFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/*! \brief Subtract SIMD integers.
 *
 * You should typically call the real-precision \ref gmx::simdXorI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static inline SimdFInt32
simdSubFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/*! \brief Multiply SIMD integers.
 *
 * You should typically call the real-precision \ref gmx::simdXorI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 *
 * \note Only the low 32 bits are retained, so this can overflow.
 */
static inline SimdFInt32
simdMulFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] * b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) comparisons, boolean, selection
 * \{
 */

/*! \brief Equality comparison of two integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx::simdCmpEqI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a==b
 */
static inline SimdFIBool
simdCmpEqFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFIBool         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/*! \brief Less-than comparison of two SIMD integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx::simdCmpLtI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a<b
 */
static inline SimdFIBool
simdCmpLtFI(SimdFInt32 a, SimdFInt32 b)
{
    SimdFIBool         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/*! \brief Logical AND on SimdFIBool.
 *
 * You should typically call the real-precision \ref gmx::simdAndIB.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static inline SimdFIBool
simdAndFIB(SimdFIBool a, SimdFIBool b)
{
    SimdFIBool        c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical OR on SimdFIBool.
 *
 * You should typically call the real-precision \ref gmx::simdOrIB.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static inline SimdFIBool
simdOrFIB(SimdFIBool a, SimdFIBool b)
{
    SimdFIBool         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns true if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx::simdAnyTrueIB.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * The actual return value for "any true" will depend on the architecture.
 * Any non-zero value should be considered truth.
 *
 * \param a SIMD boolean
 * \return True if any of the elements in a is true, otherwise 0.
 */
static inline bool
simdAnyTrueFIB(SimdFIBool a)
{
    bool anyTrue = false;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        anyTrue = anyTrue || a.b[i];
    }
    return anyTrue;
}

/*! \brief Select from \ref gmx::SimdFInt32 variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx::simdMaskI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param mask Boolean selector
 * \return Elements from a where sel is true, 0 otherwise.
 */
static inline SimdFInt32
simdMaskFI(SimdFInt32 a, SimdFIBool mask)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = mask.b[i] ? a.i[i] : 0.0f;
    }
    return c;
}

/*! \brief Select from \ref gmx::SimdFInt32 variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx::simdMaskNotI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param mask Boolean selector
 * \return Elements from a where sel is false, 0 otherwise (sic).
 */
static inline SimdFInt32
simdMaskNotFI(SimdFInt32 a, SimdFIBool mask)
{
    SimdFInt32         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = mask.b[i] ? 0.0f : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx::simdBlendI.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static inline SimdFInt32
simdBlendFI(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    SimdFInt32        d;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        d.i[i] = sel.b[i] ? b.i[i] : a.i[i];
    }
    return d;
}

/*! \}
 *
 * \name SIMD implementation conversion operations
 * \{
 */

/*! \brief Round single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx::simdCvtR2I.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, rounded to nearest integer.
 */
static inline SimdFInt32
simdCvtF2I(SimdFloat a)
{
    SimdFInt32         b;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.i[i] = std::round(a.r[i]);
    }
    return b;
};

/*! \brief Truncate single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx::simdCvttR2I.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, truncated towards zero.
 */
static inline SimdFInt32
simdCvttF2I(SimdFloat a)
{
    SimdFInt32         b;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/*! \brief Convert integer to single precision floating-point.
 *
 * You should typically call the real-precision \ref gmx::simdCvtI2R.
 *
 * \param a SIMD integer
 * \return SIMD floating-pint
 */
static inline SimdFloat
simdCvtI2F(SimdFInt32 a)
{
    SimdFloat          b;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/*! \brief Convert from float boolean to corresponding integer boolean.
 *
 * You should typically call the real-precision \ref gmx::simdCvtB2IB.
 *
 * \param a Boolean corresponding to SIMD floating-point
 * \return Boolean that can be applied to SIMD integer operations.
 */
static inline SimdFIBool
simdCvtFB2FIB(SimdFBool a)
{
    SimdFIBool         b;

    /* Integer width >= float width */
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \brief Convert from integer boolean (corresponding to float) to float boolean.
 *
 * You should typically call the real-precision \ref gmx::simdCvtIB2B.
 *
 * \param a Boolean corresponding to SIMD integer
 * \return Boolean that can be applied to SIMD floating-point.
 */
static inline SimdFBool
simdCvtFIB2FB(SimdFIBool a)
{
    SimdFBool         b;

    /* Integer width >= float width */
    for (int i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H
