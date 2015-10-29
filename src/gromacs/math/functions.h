/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares simple math functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_FUNCTIONS_H
#define GMX_MATH_FUNCTIONS_H

#include <cmath>
#include <cstdint>

#include <limits>

#include <array>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Evaluate log2(n) for integer n statically at compile time.
 *
 * Use as staticLog2<n>::value, where n must be a positive integer.
 * Negative n will be reinterpreted as the corresponding unsigned integer,
 * and you will get a compile-time error if n==0.
 * The calculation is done by recursively dividing n by 2 (until it is 1),
 * and incrementing the result by 1 in each step.
 *
 * \tparam n Value to recursively calculate log2(n) for
 */
template<std::size_t n>
struct StaticLog2
{
    static const int value = StaticLog2<n/2>::value+1; //!< Variable value used for recursive static calculation of Log2(int)
};

/*! \brief Specialization of StaticLog2<n> for n==1.
 *
 *  This specialization provides the final value in the recursion; never
 *  call it directly, but use StaticLog2<n>::value.
 */
template<>
struct StaticLog2<1>
{
    static const int value = 0; //!< Base value for recursive static calculation of Log2(int)
};

/*! \brief Specialization of StaticLog2<n> for n==0.
 *
 *  This specialization should never actually be used since log2(0) is
 *  negative infinity. However, since Log2() is often used to calculate the number
 *  of bits needed for a number, we might be using the value 0 with a conditional
 *  statement around the logarithm. Depending on the compiler the expansion of
 *  the template can occur before the conditional statement, so to avoid infinite
 *  recursion we need a specialization for the case n==0.
 */
template<>
struct StaticLog2<0>
{
    static const int value = -1; //!< Base value for recursive static calculation of Log2(int)
};


/*! \brief Compute floor of logarithm to base 2
 *
 *  \param x Argument, must be positive
 *
 *  \return log2(x)
 */
unsigned int
log2I(std::size_t x);

/*! \brief Find greatest common divisor of two numbers
 *
 *  \param p First number, positive
 *  \param q Second number, positive
 *
 * \return Greatest common divisor of p and q
 */
std::size_t
greatestCommonDivisor(std::size_t p, std::size_t q);
    
    
/*! \brief Declaration of external 8-bit table for 1.0/sqrt(x) exponent lookup */
extern const
std::array<std::uint32_t,256>  invsqrtExponentTable;
    
/*! \brief Declaration of external 12-bit table for 1.0/sqrt(x) fraction lookup */
extern const
std::array<std::uint32_t,4096> invsqrtFractionTable;

/*! \brief Calculate 1.0/sqrt(x) in single precision
 *
 * \param  x  Positive value to calculate inverse square root for
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 * \return 1.0/sqrt(x)
 */
static inline float
invsqrt(float x)
{
    const int     fractionBits  = std::numeric_limits<float>::digits - 1;
    const int     lookupBits    = 12;
    
    union
    {
        std::uint32_t  i;
        float          f;
    }
    convert;

    GMX_ASSERT(x>=0,"gmx::invsqrt() called with negative argument");

    convert.f                   = x;
    std::uint32_t ix            = convert.i;
    std::uint32_t exponentIndex = (ix & 0x7f800000) >> fractionBits;
    std::uint32_t fractionIndex = (ix & 0x00ffffff) >> (fractionBits+1-lookupBits);
    convert.i                   = invsqrtExponentTable[exponentIndex] | invsqrtFractionTable[fractionIndex];
    float         lookup        = convert.f;
    
    return 0.5f*lookup*(3.0f-x*lookup*lookup);
}

/*! \brief Calculate 1.0/sqrt(x) in double precision, but single range
 *
 * \param  x  Positive value to calculate inverse square root for, must be
 *            in the input domain valid for single precision.
 *
 * This routine performs the first table lookup step in single precision, after
 * which the Newton-Rhapson iterations are performed in double. For this reason,
 * the result will only be valid for arguments that can be represented in single
 * precision. If this does not hold for your data, use the system 1.0/sqrt(x).
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 * \return 1.0/sqrt(x)
 */
static inline double
invsqrt(double x)
{
    double lookup = invsqrt(static_cast<float>(x));
    return 0.5*lookup*(3.0-x*lookup*lookup);
}

/*! \brief Calculate 1.0/sqrt(x) for integer x in double precision.
 *
 * \param  x  Positive value to calculate inverse square root for, must be
 *            in the input domain valid for single precision.
 *
 * This routine performs the first table lookup step in single precision, after
 * which the Newton-Rhapson iterations are performed in double. For this reason,
 * the result will only be valid for arguments that can be represented in single
 * precision. If this does not hold for your data, use the system 1.0/sqrt(x).
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 * \return 1.0/sqrt(x)
 */
static inline double
invsqrt(int x)
{
    return invsqrt(static_cast<double>(x));
}
    
/*! \brief GROMACS version of sqrt(x) that uses invsqrt(x), single precision
 *
 * \param  x  Positive value to calculate square root of, cannot be negative.
 *
 * This function achieves higher performance by forgoing the usual checks on
 * the input values, it does not try to get the rounding effects on the least
 * significant bit right, and it does not care about e.g. returning NaN when
 * the argument is negative.
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 * To avoid unintentional mixing with sqrt() from the standard library we use
 * a slightly different name.
 *
 * \return sqrt(x)
 */
static inline float
fastsqrt(float x)
{
    GMX_ASSERT(x>=0,"gmx::sqrt() called with negative or zero argument");
    return (x>0.0f) ? x*invsqrt(x) : 0.0f;
}

/*! \brief GROMACS version of sqrt(x) that uses invsqrt(x), double precision
 *
 * \param  x  Positive value to calculate square root of, cannot be negative.
 *
 * * This function achieves higher performance by forgoing the usual checks on
 * the input values, it only gets 44 bits correct (twice that of single
 * precision), and it does not care about e.g. returning NaN when
 * the argument is negative.
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 * To avoid unintentional mixing with sqrt() from the standard library we use
 * a slightly different name.
 *
 * \return sqrt(x)
 */
static inline double
fastsqrt(double x)
{
    GMX_ASSERT(x>=0,"gmx::sqrt() called with negative or zero argument");
    return (x>0.0) ? x*invsqrt(x) : 0.0;
}

/*! \brief GROMACS version of sqrt(x) that uses invsqrt(x), integer x
 *
 * \param  x  Positive value to calculate square root of, cannot be negative.
 *
 * * This function achieves higher performance by forgoing the usual checks on
 * the input values, it only gets 44 bits correct (twice that of single
 * precision), and it does not care about e.g. returning NaN when
 * the argument is negative.
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 * To avoid unintentional mixing with sqrt() from the standard library we use
 * a slightly different name.
 *
 * \return sqrt(x)
 */
static inline double
fastsqrt(int x)
{
    return fastsqrt(static_cast<double>(x));
}

/*! \brief Calculate inverse cube root of x in single precision
 *
 *  \param  x  Argument
 *
 *  \return x^(-1/3)
 *
 *  This routine is typically faster than using std::pow().
 */
static inline float
invcbrt(float x)
{
    return 1.0/std::cbrt(x);
}
    
/*! \brief Calculate inverse sixth root of x in double precision
 *
 *  \param  x  Argument
 *
 *  \return x^(-1/3)
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
invcbrt(double x)
{
    return 1.0/std::cbrt(x);
}
    
/*! \brief Calculate inverse sixth root of integer x in double precision
 *
 *  \param  x  Argument
 *
 *  \return x^(-1/3)
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
invcbrt(int x)
{
    return 1.0/std::cbrt(x);
}

/*! \brief Calculate sixth root of x in single precision.
 *
 *  \param  x  Argument, must be greater than or equal to zero.
 *
 *  \return x^(1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline float
sixthroot(float x)
{
    GMX_ASSERT(x>=0,"gmx::sixthroot() called with negative argument");
    return std::sqrt(std::cbrt(x));
}

/*! \brief Calculate sixth root of x in double precision.
 *
 *  \param  x  Argument, must be greater than or equal to zero.
 *
 *  \return x^(1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
sixthroot(double x)
{
    GMX_ASSERT(x>=0,"gmx::sixthroot() called with negative argument");
    return std::sqrt(std::cbrt(x));
}

/*! \brief Calculate sixth root of integer x, return double.
 *
 *  \param  x  Argument, must be greater than or equal to zero.
 *
 *  \return x^(1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is non-negative in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
sixthroot(int x)
{
    GMX_ASSERT(x>=0,"gmx::sixthroot() called with negative argument");
    return std::sqrt(std::cbrt(x));
}

/*! \brief Calculate inverse sixth root of x in single precision
 *
 *  \param  x  Argument, must be greater than zero.
 *
 *  \return x^(-1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline float
invsixthroot(float x)
{
    GMX_ASSERT(x>0,"gmx::invsixthroot() called with negative or zero argument");
    return invsqrt(std::cbrt(x));
}

/*! \brief Calculate inverse sixth root of x in double precision
 *
 *  \param  x  Argument, must be greater than zero.
 *
 *  \return x^(-1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
invsixthroot(double x)
{
    GMX_ASSERT(x>0,"gmx::invsixthroot() called with negative or zero argument");
    return invsqrt(std::cbrt(x));
}

/*! \brief Calculate inverse sixth root of integer x in double precision
 *
 *  \param  x  Argument, must be greater than zero.
 *
 *  \return x^(-1/6)
 *
 * Since this is a performance-critical function, we only use an assert to
 * check that the argument is positive in debug builds; release code will
 * not check the value.
 *
 *  This routine is typically faster than using std::pow().
 */
static inline double
invsixthroot(int x)
{
    return invsqrt(std::cbrt(x));
}

/*! \brief calculate x^2
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^2
 */
template <typename T>
T
square(T x)
{
    return x*x;
}

/*! \brief calculate x^3
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^3
 */
template <typename T>
T
power3(T x)
{
    return x*square(x);
}

/*! \brief calculate x^4
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^4
 */
template <typename T>
T
power4(T x)
{
    return square(square(x));
}

/*! \brief calculate x^5
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^5
 */
template <typename T>
T
power5(T x)
{
    return x*power4(x);
}

/*! \brief calculate x^6
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^6
 */
template <typename T>
T
power6(T x)
{
    return square(power3(x));
}

/*! \brief calculate x^12
 *
 *  \tparam T  Type of argument and return value
 *  \param  x  argument
 *
 *  \return x^12
 */
template <typename T>
T
power12(T x)
{
    return square(power6(x));
}

} // namespace gmx


#endif // GMX_MATH_FUNCTIONS_H
