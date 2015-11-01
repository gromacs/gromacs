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
#include <array>

#include "impl_reference_definitions.h"

namespace gmx
{

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/* \name SIMD implementation data types and built-in conversions between types
 * \{
 */

class SimdFloat;
class SimdFInt32;
class SimdFBool;
class SimdFIBool;


/*! \libinternal \brief Represent of a SIMD float result, either data or intermediate
 *
 *  \tparam E  Type of the expression, such as SimdFloat or an intermediate result.
 *
 *  This class is used to represent intermediate results in some compound
 *  expressions. The class will not really do anything, but by using it we get
 *  a way to specify that an intermediate value is actually the result of e.g.
 *  a multiply operation, and this makes it possible to automatically turn adjacent
 *  multiply and add operations into a single fused multiply-add instruction.
 */
template <typename E>
class SimdFloatExpression
{
    public:
        /*! \brief Return a constant reference to the expression itself */
        operator E const &() const { return static_cast<const E &>(*this); }
};

/*! \libinternal \brief Float SIMD variable. Supported if GMX_SIMD_HAVE_FLOAT.
 */
class SimdFloat : public SimdFloatExpression<SimdFloat>
{
    public:
        SimdFloat() {}

        /*! \brief Construct from scalar */
        SimdFloat(float f) { simdInternal_.fill(f); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<float, GMX_SIMD_FLOAT_WIDTH>  simdInternal_;
};

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from float.
 *
 * This is also the widest integer SIMD type. Available if GMX_SIMD_HAVE_FLOAT is true.
 *
 * \note The integer SIMD type will always be available, but on architectures
 * that do not have any real integer SIMD support it might be defined as the
 * floating-point type. This will work fine, since there are separate defines
 * for whether the implementation can actually do any operations on integer
 * SIMD types.
 */
class SimdFInt32
{
    public:
        SimdFInt32() {}

        /*! \brief Construct from scalar signed int */
        SimdFInt32(std::int32_t i) { simdInternal_.fill(i); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<std::int32_t, GMX_SIMD_FINT32_WIDTH>  simdInternal_;
};

/*! \libinternal \brief Boolean type for float SIMD data. */
class SimdFBool
{
    public:
        SimdFBool() {}

        /*! \brief Construct from scalar bool */
        SimdFBool(bool b) { simdInternal_.fill(b); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<bool, GMX_SIMD_FLOAT_WIDTH>  simdInternal_;
};

/*! \libinternal \brief Boolean type for integer datatypes corresponding to float SIMD. */
class SimdFIBool
{
    public:
        SimdFIBool() {}

        /*! \brief Construct from scalar bool */
        SimdFIBool(bool b) { simdInternal_.fill(b); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<bool, GMX_SIMD_FINT32_WIDTH>  simdInternal_;
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
load(const float *m)
{
    SimdFloat a;

    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(float)) == 0);

    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

/*! \brief Store the contents of the SIMD float variable pr to aligned memory m.
 *
 * \param[out] m Pointer to memory, aligned to SIMD width.
 * \param a SIMD variable to store
 */
static inline void
store(float *m, SimdFloat a)
{
    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(float)) == 0);

    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
}

/*! \brief Load SIMD float from unaligned memory.
 *
 * Available if \ref GMX_SIMD_HAVE_LOADU is true.
 *
 * \param m Pointer to memory, no alignment requirement.
 * \return SIMD variable with data loaded.
 */
static inline SimdFloat
loadU(const float *m)
{
    SimdFloat a;
    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

/*! \brief Store SIMD float to unaligned memory.
 *
 * Available if \ref GMX_SIMD_HAVE_STOREU is true.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a SIMD variable to store.
 */
static inline void
storeU(float *m, SimdFloat a)
{
    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
}

/*! \} */


/*!
 * \name SIMD implementation load/store operations for integers (corresponding to float)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically just call \ref gmx::load(), which uses proxy objects
 * internally to handle all types rather than adding the suffix used here.
 *
 * \param m Pointer to memory, aligned to integer SIMD width.
 * \return SIMD integer variable.
 */
static inline SimdFInt32
loadFI(const std::int32_t * m)
{
    SimdFInt32 a;

    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(std::int32_t)) == 0);

    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
};

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * \param m Memory aligned to integer SIMD width.
 * \param a SIMD variable to store.
 */
static inline void
store(std::int32_t * m, SimdFInt32 a)
{
    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(std::int32_t)) == 0);

    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx::SimdFloat.
 *
 * You should typically just call \ref gmx::loadU(), which uses proxy objects
 * internally to handle all types rather than adding the suffix used here.
 *
 * Available if \ref GMX_SIMD_HAVE_LOADU is true.
 *
 * \param m Pointer to memory, no alignment requirements.
 * \return SIMD integer variable.
 */
static inline SimdFInt32
loadUFI(const std::int32_t *m)
{
    SimdFInt32 a;
    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx::SimdFloat.
 *
 * Available if \ref GMX_SIMD_HAVE_STOREU is true.
 *
 * \param m Memory pointer, no alignment requirements.
 * \param a SIMD variable to store.
 */
static inline void
storeU(std::int32_t * m, SimdFInt32 a)
{
    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
}

/*! \brief Extract element with index i from \ref gmx::SimdFInt32.
 *
 * Available if \ref GMX_SIMD_HAVE_FINT32_EXTRACT is true.
 *
 * \tparam index Compile-time constant, position to extract (first position is 0)
 * \param  a     SIMD variable from which to extract value.
 * \return Single integer from position index in SIMD variable.
 */
template<int index>
static inline std::int32_t
extract(SimdFInt32 a)
{
    return a.simdInternal_[index];
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point bitwise logical operations
 * \{
 */

/*! \brief Bitwise and for two SIMD float variables.
 *
 * Supported if \ref GMX_SIMD_HAVE_LOGICAL is true.
 *
 * \param a data1
 * \param b data2
 * \return data1 & data2
 */
static inline SimdFloat
operator&(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i & conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \brief Bitwise andnot for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 */
static inline SimdFloat
andNot(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = ~conv1.i & conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \brief Bitwise or for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \param a data1
 * \param b data2
 * \return data1 | data2
 */
static inline SimdFloat
operator|(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i | conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \brief Bitwise xor for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \param a data1
 * \param b data2
 * \return data1 ^ data2
 */
static inline SimdFloat
operator^(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    union
    {
        float         r;
        std::int32_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i ^ conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point arithmetics
 * \{
 */

/*! \brief Add two float SIMD variables.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static inline SimdFloat
operator+(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] + b.simdInternal_[i];
    }
    return res;
}

/*! \brief Subtract two SIMD variables.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static inline SimdFloat
operator-(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] - b.simdInternal_[i];
    }
    return res;
}

/*! \brief SIMD floating-point negate.
 *
 * \param a Any floating-point value
 * \return -a
 */
static inline SimdFloat
operator-(SimdFloat a)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = -a.simdInternal_[i];
    }
    return res;
}

/*! \brief Multiply two SIMD variables.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 */
static inline SimdFloat
operator*(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] * b.simdInternal_[i];
    }
    return res;
}

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b+c
 */
static inline SimdFloat
fma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return a*b+c;
}

/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b-c
 */
static inline SimdFloat
fms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return a*b-c;
}

/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b+c
 */
static inline SimdFloat
fnma(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return c-a*b;
}

/*! \brief Fused-negated-multiply-add. Result is -a*b-c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b-c
 */
static inline SimdFloat
fnms(SimdFloat a, SimdFloat b, SimdFloat c)
{
    return -a*b-c;
}

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 */
static inline SimdFloat
rsqrt(SimdFloat x)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = 1.0f / std::sqrt(x.simdInternal_[i]);
    }
    return res;
};

/*! \brief SIMD 1.0/x lookup.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x!=0
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 */
static inline SimdFloat
rcp(SimdFloat x)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = 1.0f / x.simdInternal_[i];
    }
    return res;
};

/*! \brief Add two SIMD variables, masked version.
 *
 * \param a term1
 * \param b term2
 * \param m mask
 * \return a+b where mask is true, 0.0 otherwise.
 */
static inline SimdFloat
addMask(SimdFloat a, SimdFloat b, SimdFBool m)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = m.simdInternal_[i] ? (a.simdInternal_[i] + b.simdInternal_[i]) : 0.0f;
    }
    return res;
}

/*! \brief Multiply two SIMD variables, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param m mask
 * \return a*b where mask is true, 0.0 otherwise.
 */
static inline SimdFloat
mulMask(SimdFloat a, SimdFloat b, SimdFBool m)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = m.simdInternal_[i] ? (a.simdInternal_[i] * b.simdInternal_[i]) : 0.0f;
    }
    return res;
}

/*! \brief SIMD fused multiply-add, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param b term
 * \param m mask
 * \return a*b+c where mask is true, 0.0 otherwise.
 */
static inline SimdFloat
fmaMask(SimdFloat a, SimdFloat b, SimdFloat c, SimdFBool m)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = m.simdInternal_[i] ? (a.simdInternal_[i] * b.simdInternal_[i] + c.simdInternal_[i]) : 0.0f;
    }
    return res;
}

/*! \brief SIMD 1.0/sqrt(x) lookup, masked version.
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
rsqrtMask(SimdFloat x, SimdFBool m)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (m.simdInternal_[i] != 0) ? 1.0f / std::sqrt(x.simdInternal_[i]) : 0.0f;
    }
    return res;
}

/*! \brief SIMD 1.0/x lookup, masked version.
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
rcpMask(SimdFloat x, SimdFBool m)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (m.simdInternal_[i] != 0) ? 1.0f / x.simdInternal_[i] : 0.0f;
    }
    return res;
}

/*! \brief SIMD Floating-point fabs().
 *
 * \param a any floating point values
 * \return fabs(a) for each element.
 */
static inline SimdFloat
abs(SimdFloat a)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::abs(a.simdInternal_[i]);
    }
    return res;
}

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
static inline SimdFloat
max(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::max(a.simdInternal_[i], b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return min(a,b) for each element.
 */
static inline SimdFloat
min(SimdFloat a, SimdFloat b)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::min(a.simdInternal_[i], b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Round to nearest integer value (in floating-point format).
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
round(SimdFloat a)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::round(a.simdInternal_[i]);
    }
    return res;
}

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This is truncation towards zero, not floor(). The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but floor() frequently isn't.
 */
static inline SimdFloat
trunc(SimdFloat a)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::trunc(a.simdInternal_[i]);
    }
    return res;
}

/*! \brief Extract (integer) exponent from single precision SIMD.
 *
 * \param a Any floating-point value
 * \return Exponent value, represented in floating-point format.
 *
 * The IEEE754 exponent field is selected, the bias removed, and it is converted to
 * a normal floating-point SIMD.
 */
static inline SimdFloat
frexp(SimdFloat value, SimdFInt32 * exponent)
{
    SimdFloat fraction;

    for (std::size_t i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        fraction.simdInternal_[i] = std::frexp(value.simdInternal_[i], &exponent->simdInternal_[i]);
    }
    return fraction;
}

/*! \brief Multiplies a floating point value by the number 2 raised to an exp power.
 *
 * \param value Floating-point number to multiply with new exponent
 * \param exponent Integer that will not overflow as 2^exponent.
 * \return value*2^exponent
 */
static inline SimdFloat
ldexp(SimdFloat value, SimdFInt32 exponent)
{
    SimdFloat           res;

    for (std::size_t i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        res.simdInternal_[i] = std::ldexp(value.simdInternal_[i], exponent.simdInternal_[i]);
    }
    return res;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point comparisons, boolean, selection.
 * \{
 */

/*! \brief SIMD a==b for single SIMD.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a==b.
 *
 * Beware that exact floating-point comparisons are difficult.
 */
static inline SimdFBool
operator==(SimdFloat a, SimdFloat b)
{
    SimdFBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] == b.simdInternal_[i]);
    }
    return res;
}

/*! \brief SIMD a!=0 for single SIMD.
 *
 * \param a value
 * \return Each element of the boolean will be true if any bit in a is nonzero.
 *         The behaviour for negative zero is undefined, and should not be
 *         relied on - it will depend on the architecture.
 */
static inline SimdFBool
isNonZero(SimdFloat a)
{
    SimdFBool         res;

    for (std::size_t i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] != 0.0f);
    }
    return res;
}

/*! \brief SIMD a<b for single SIMD.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<b.
 */
static inline SimdFBool
operator<(SimdFloat a, SimdFloat b)
{
    SimdFBool          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] < b.simdInternal_[i]);
    }
    return res;
}

/*! \brief SIMD a<=b for single SIMD.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<=b.
 */
static inline SimdFBool
operator<=(SimdFloat a, SimdFloat b)
{
    SimdFBool          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] <= b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Logical \a and on single precision SIMD booleans.
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
operator&&(SimdFBool a, SimdFBool b)
{
    SimdFBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] && b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Logical \a or on single precision SIMD booleans.
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
operator||(SimdFBool a, SimdFBool b)
{
    SimdFBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] || b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * \param a Logical variable.
 * \return true if any element in a is true, otherwise false.
 *
 * The actual return value for truth will depend on the architecture,
 * so any non-zero value is considered truth.
 */
static inline bool
anyTrue(SimdFBool a)
{
    bool res = false;

    for (std::size_t i = 0; i < a.simdInternal_.size(); i++)
    {
        res = res || a.simdInternal_[i];
    }
    return res;
}

/*! \brief Select from single precision SIMD variable where boolean is true.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for true, 0 for false.
 */
static inline SimdFloat
selectMask(SimdFloat a, SimdFBool mask)
{
    SimdFloat          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? a.simdInternal_[i] : 0.0f;
    }
    return res;
}

/*! \brief Select from single precision SIMD variable where boolean is false.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for false, 0 for true (sic).
 */
static inline SimdFloat
selectNotMask(SimdFloat a, SimdFBool mask)
{
    SimdFloat          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? 0.0f : a.simdInternal_[i];
    }
    return res;
}

/*! \brief Vector-blend SIMD selection.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static inline SimdFloat
blend(SimdFloat a, SimdFloat b, SimdFBool sel)
{
    SimdFloat         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = sel.simdInternal_[i] ? b.simdInternal_[i] : a.simdInternal_[i];
    }
    return res;
}

/*! \brief Return sum of all elements in SIMD float variable.
 *
 * \param a SIMD variable to reduce/sum.
 * \return The sum of all elements in the argument variable.
 *
 */
static inline float
reduce(SimdFloat a)
{
    float sum = 0.0f;

    for (std::size_t i = 0; i < a.simdInternal_.size(); i++)
    {
        sum += a.simdInternal_[i];
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
operator<<(SimdFInt32 a, int n)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] << n;
    }
    return res;
}

/*! \brief SIMD integer shift right logical, based on immediate value.
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
operator>>(SimdFInt32 a, int n)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] >> n;
    }
    return res;
}

/*! \brief Integer SIMD bitwise and.
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
operator&(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] & b.simdInternal_[i];
    }
    return res;
}

/*! \brief Integer SIMD bitwise not/complement.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \note You can \a not use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx::simdMaskI instead.
 *
 * \param a integer SIMD
 * \return ~a
 */
static inline SimdFInt32
operator~(SimdFInt32 a)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = ~a.simdInternal_[i];
    }
    return res;
}

/*! \brief Integer SIMD bitwise not/complement.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \note You can \a not use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx::simdMaskI instead.
 *
 * \param a integer SIMD
 * \return ~a
 */
static inline SimdFInt32
andNot(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = ~a.simdInternal_[i] & b.simdInternal_[i];
    }
    return res;
}

/*! \brief Integer SIMD bitwise or.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \| b (bitwise or)
 */
static inline SimdFInt32
operator|(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] | b.simdInternal_[i];
    }
    return res;
}

/*! \brief Integer SIMD bitwise xor.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is 1.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a ^ b (bitwise xor)
 */
static inline SimdFInt32
operator^(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] ^ b.simdInternal_[i];
    }
    return res;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) arithmetics
 * \{
 */

/*! \brief Add SIMD integers.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static inline SimdFInt32
operator+(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] + b.simdInternal_[i];
    }
    return res;
}

/*! \brief Subtract SIMD integers.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static inline SimdFInt32
operator-(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] - b.simdInternal_[i];
    }
    return res;
}

/*! \brief Multiply SIMD integers.
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
operator*(SimdFInt32 a, SimdFInt32 b)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] * b.simdInternal_[i];
    }
    return res;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) comparisons, boolean, selection
 * \{
 */

/*! \brief Equality comparison of two integers corresponding to float values.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a==b
 */
static inline SimdFIBool
operator==(SimdFInt32 a, SimdFInt32 b)
{
    SimdFIBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] == b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Less-than comparison of two SIMD integers corresponding to float values.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a<b
 */
static inline SimdFIBool
operator<(SimdFInt32 a, SimdFInt32 b)
{
    SimdFIBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] < b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Check if any bit is set in each element
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer
 * \return SIMD integer boolean with true for elements where any bit is set
 */
static inline SimdFIBool
isNonZero(SimdFInt32 a)
{
    SimdFIBool         c;

    for (int i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.simdInternal_[i] = (a.simdInternal_[i] != 0);
    }
    return c;
}

/*! \brief Logical AND on SimdFIBool.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static inline SimdFIBool
operator&&(SimdFIBool a, SimdFIBool b)
{
    SimdFIBool        res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] && b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Logical OR on SimdFIBool.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static inline SimdFIBool
operator||(SimdFIBool a, SimdFIBool b)
{
    SimdFIBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] || b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Returns true if any of the boolean in x is True, otherwise 0.
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
anyTrue(SimdFIBool a)
{
    bool res = false;

    for (std::size_t i = 0; i < a.simdInternal_.size(); i++)
    {
        res = res || a.simdInternal_[i];
    }
    return res;
}

/*! \brief Select from \ref gmx::SimdFInt32 variable where boolean is true.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param mask Boolean selector
 * \return Elements from a where sel is true, 0 otherwise.
 */
static inline SimdFInt32
selectMask(SimdFInt32 a, SimdFIBool mask)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? a.simdInternal_[i] : 0.0f;
    }
    return res;
}

/*! \brief Select from \ref gmx::SimdFInt32 variable where boolean is false.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is 1.
 *
 * \param a SIMD integer to select from
 * \param mask Boolean selector
 * \return Elements from a where sel is false, 0 otherwise (sic).
 */
static inline SimdFInt32
selectNotMask(SimdFInt32 a, SimdFIBool mask)
{
    SimdFInt32         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? 0.0f : a.simdInternal_[i];
    }
    return res;
}

/*! \brief Vector-blend SIMD selection.
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
blend(SimdFInt32 a, SimdFInt32 b, SimdFIBool sel)
{
    SimdFInt32        res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = sel.simdInternal_[i] ? b.simdInternal_[i] : a.simdInternal_[i];
    }
    return res;
}

/*! \}
 *
 * \name SIMD implementation conversion operations
 * \{
 */

/*! \brief Round single precision floating point to integer.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, rounded to nearest integer.
 */
static inline SimdFInt32
cvtR2I(SimdFloat a)
{
    SimdFInt32         b;

    for (std::size_t i = 0; i < b.simdInternal_.size(); i++)
    {
        b.simdInternal_[i] = std::round(a.simdInternal_[i]);
    }
    return b;
};

/*! \brief Truncate single precision floating point to integer.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, truncated to nearest integer.
 */
static inline SimdFInt32
cvttR2I(SimdFloat a)
{
    SimdFInt32         b;

    for (std::size_t i = 0; i < b.simdInternal_.size(); i++)
    {
        b.simdInternal_[i] = std::trunc(a.simdInternal_[i]);
    }
    return b;
};

/*! \brief Convert integer to single precision floating point.
 *
 * \param a SIMD integer
 * \return SIMD floating-point
 */
static inline SimdFloat
cvtI2R(SimdFInt32 a)
{
    SimdFloat         b;

    for (std::size_t i = 0; i < b.simdInternal_.size(); i++)
    {
        b.simdInternal_[i] = a.simdInternal_[i];
    }
    return b;
};

/*! \brief Convert from single precision boolean to corresponding integer boolean
 *
 * \param a SIMD floating-point boolean
 * \return SIMD integer boolean
 */
static inline SimdFIBool
cvtB2IB(SimdFBool a)
{
    SimdFIBool         b;

    for (std::size_t i = 0; i < b.simdInternal_.size(); i++)
    {
        b.simdInternal_[i] = a.simdInternal_[i];
    }
    return b;
};

/*! \brief Convert from integer boolean to corresponding single precision boolean
 *
 * \param a SIMD integer boolean
 * \return SIMD floating-point boolean
 */
static inline SimdFBool
cvtIB2B(SimdFIBool a)
{
    SimdFBool         b;

    for (std::size_t i = 0; i < b.simdInternal_.size(); i++)
    {
        b.simdInternal_[i] = a.simdInternal_[i];
    }
    return b;
};

/*! \} */

/*! \} */
/*! \endcond */

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_SIMD_FLOAT_H
