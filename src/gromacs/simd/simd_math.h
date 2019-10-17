/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_SIMD_MATH_H
#define GMX_SIMD_SIMD_MATH_H

/*! \libinternal \file
 *
 * \brief Math functions for SIMD datatypes.
 *
 * \attention This file is generic for all SIMD architectures, so you cannot
 * assume that any of the optional SIMD features (as defined in simd.h) are
 * present. In particular, this means you cannot assume support for integers,
 * logical operations (neither on floating-point nor integer values), shifts,
 * and the architecture might only have SIMD for either float or double.
 * Second, to keep this file clean and general, any additions to this file
 * must work for all possible SIMD architectures in both single and double
 * precision (if they support it), and you cannot make any assumptions about
 * SIMD width.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#include "config.h"

#include <cmath>

#include <limits>

#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

#if GMX_SIMD

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Implementation accuracy settings
 *  \{
 */

/*! \} */

#    if GMX_SIMD_HAVE_FLOAT

/*! \name Single precision SIMD math functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/****************************************
 * SINGLE PRECISION SIMD MATH FUNCTIONS *
 ****************************************/

#        if !GMX_SIMD_HAVE_NATIVE_COPYSIGN_FLOAT
/*! \brief Composes floating point value with the magnitude of x and the sign of y.
 *
 * \param x Values to set sign for
 * \param y Values used to set sign
 * \return  Magnitude of x, sign of y
 */
static inline SimdFloat gmx_simdcall copysign(SimdFloat x, SimdFloat y)
{
#            if GMX_SIMD_HAVE_LOGICAL
    return abs(x) | (SimdFloat(GMX_FLOAT_NEGZERO) & y);
#            else
    return blend(abs(x), -abs(x), y < setZero());
#            endif
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_FLOAT
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdFloat gmx_simdcall rsqrtIter(SimdFloat lu, SimdFloat x)
{
    SimdFloat tmp1 = x * lu;
    SimdFloat tmp2 = SimdFloat(-0.5F) * lu;
    tmp1           = fma(tmp1, lu, SimdFloat(-3.0F));
    return tmp1 * tmp2;
}
#        endif

/*! \brief Calculate 1/sqrt(x) for SIMD float.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall invsqrt(SimdFloat x)
{
    SimdFloat lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD floats.
 *
 * \param x0  First set of arguments, x0 must be in single range (see below).
 * \param x1  Second set of arguments, x1 must be in single range (see below).
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 *
 * \note Both arguments must be larger than GMX_FLOAT_MIN and smaller than
 *       GMX_FLOAT_MAX, i.e. within the range of single precision.
 *       For the single precision implementation this is obviously always
 *       true for positive values, but for double precision it adds an
 *       extra restriction since the first lookup step might have to be
 *       performed in single precision on some architectures. Note that the
 *       responsibility for checking falls on you - this routine does not
 *       check arguments.
 */
static inline void gmx_simdcall invsqrtPair(SimdFloat x0, SimdFloat x1, SimdFloat* out0, SimdFloat* out1)
{
    *out0 = invsqrt(x0);
    *out1 = invsqrt(x1);
}

#        if !GMX_SIMD_HAVE_NATIVE_RCP_ITER_FLOAT
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the reciprocal.
 *
 *  \param lu Approximation of 1/x, typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/x.
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdFloat gmx_simdcall rcpIter(SimdFloat lu, SimdFloat x)
{
    return lu * fnma(lu, x, SimdFloat(2.0F));
}
#        endif

/*! \brief Calculate 1/x for SIMD float.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall inv(SimdFloat x)
{
    SimdFloat lu = rcp(x);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}

/*! \brief Division for SIMD floats
 *
 * \param nom    Nominator
 * \param denom  Denominator, with magnitude in range (GMX_FLOAT_MIN,GMX_FLOAT_MAX).
 *               For single precision this is equivalent to a nonzero argument,
 *               but in double precision it adds an extra restriction since
 *               the first lookup step might have to be performed in single
 *               precision on some architectures. Note that the responsibility
 *               for checking falls on you - this routine does not check arguments.
 *
 * \return nom/denom
 *
 * \note This function does not use any masking to avoid problems with
 *       zero values in the denominator.
 */
static inline SimdFloat gmx_simdcall operator/(SimdFloat nom, SimdFloat denom)
{
    return nom * inv(denom);
}

/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD float.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX for masked-in entries.
 *           See \ref invsqrt for the discussion about argument restrictions.
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdFloat maskzInvsqrt(SimdFloat x, SimdFBool m)
{
    SimdFloat lu = maskzRsqrt(x, m);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate 1/x for SIMD float, masked version.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX for masked-in entries.
 *           See \ref invsqrt for the discussion about argument restrictions.
 *  \param m Mask
 *  \return 1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdFloat gmx_simdcall maskzInv(SimdFloat x, SimdFBool m)
{
    SimdFloat lu = maskzRcp(x, m);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate sqrt(x) for SIMD floats
 *
 *  \tparam opt By default, this function checks if the input value is 0.0
 *              and masks this to return the correct result. If you are certain
 *              your argument will never be zero, and you know you need to save
 *              every single cycle you can, you can alternatively call the
 *              function as sqrt<MathOptimization::Unsafe>(x).
 *
 *  \param  x   Argument that must be in range 0 <=x <= GMX_FLOAT_MAX, since the
 *              lookup step often has to be implemented in single precision.
 *              Arguments smaller than GMX_FLOAT_MIN will always lead to a zero
 *              result, even in double precision. If you are using the unsafe
 *              math optimization parameter, the argument must be in the range
 *              GMX_FLOAT_MIN <= x <= GMX_FLOAT_MAX.
 *
 *  \return sqrt(x). The result is undefined if the input value does not fall
 *          in the allowed range specified for the argument.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall sqrt(SimdFloat x)
{
    if (opt == MathOptimization::Safe)
    {
        SimdFloat res = maskzInvsqrt(x, setZero() < x);
        return res * x;
    }
    else
    {
        return x * invsqrt(x);
    }
}

/*! \brief Cube root for SIMD floats
 *
 * \param x      Argument to calculate cube root of. Can be negative or zero,
 *               but NaN or Inf values are not supported. Denormal values will
 *               be treated as 0.0.
 * \return       Cube root of x.
 */
static inline SimdFloat gmx_simdcall cbrt(SimdFloat x)
{
    const SimdFloat signBit(GMX_FLOAT_NEGZERO);
    const SimdFloat minFloat(std::numeric_limits<float>::min());
    // Bias is 128-1 = 127, which is not divisible by 3. Since the largest-magnitude
    // negative exponent from frexp() is -126, we can subtract one more unit to get 126
    // as offset, which is divisible by 3 (result 42). To avoid clang warnings about fragile integer
    // division mixed with FP, we let the divided value (42) be the original constant.
    const std::int32_t offsetDiv3(42);
    const SimdFloat    c2(-0.191502161678719066F);
    const SimdFloat    c1(0.697570460207922770F);
    const SimdFloat    c0(0.492659620528969547F);
    const SimdFloat    one(1.0F);
    const SimdFloat    two(2.0F);
    const SimdFloat    three(3.0F);
    const SimdFloat    oneThird(1.0F / 3.0F);
    const SimdFloat    cbrt2(1.2599210498948731648F);
    const SimdFloat    sqrCbrt2(1.5874010519681994748F);

    // To calculate cbrt(x) we first take the absolute value of x but save the sign,
    // since cbrt(-x) = -cbrt(x). Then we only need to consider positive values for
    // the main step.
    // A number x is represented in IEEE754 as fraction*2^e. We rewrite this as
    // x=fraction*2^(3*n)*2^m, where e=3*n+m, and m is a remainder.
    // The cube root can the be evaluated by calculating the cube root of the fraction
    // limited to the mantissa range, multiplied by 2^mod (which is either 1, +/-2^(1/3) or
    // +/-2^(2/3), and then we load this into a new IEEE754 fp number with the exponent 2^n, where
    // n is the integer part of the original exponent divided by 3.

    SimdFloat xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdFloat xAbs     = andNot(signBit, x);   // select everthing but the sign bit => abs(x)
    SimdFBool xIsNonZero = (minFloat <= xAbs); // treat denormals as 0

    SimdFInt32 exponent;
    SimdFloat  y = frexp(xAbs, &exponent);
    // For the mantissa (y) we will use a limited-range approximation of cbrt(y),
    // by first using a polynomial and then evaluating
    // Transform y to z = c2*y^2 + c1*y + c0, then w = z^3, and finally
    // evaluate the quotient q = z * (w + 2 * y) / (2 * w + y).
    SimdFloat z        = fma(fma(y, c2, c1), y, c0);
    SimdFloat w        = z * z * z;
    SimdFloat nom      = z * fma(two, y, w);
    SimdFloat invDenom = inv(fma(two, w, y));

    // Handle the exponent. In principle there are beautiful ways to do this with custom 16-bit
    // division converted to multiplication... but we can't do that since our SIMD layer cannot
    // assume the presence of integer shift operations!
    // However, when I first worked with the integer algorithm I still came up with a neat
    // optimization, so I'll describe the full algorithm here in case we ever want to use it
    // in the future:
    //
    // Our dividend is signed, which is a complication, but let's consider the unsigned case
    // first: Division by 3 corresponds to multiplication by 1010101... Since we also know
    // our dividend is less than 16 bits (exponent range) we can accomplish this by
    // multiplying with 21845 (which is almost 2^16/3 - 21845.333 would be exact) and then
    // right-shifting by 16 bits to divide out the 2^16 part.
    // If we add 1 to the dividend to handle the extra 0.333, the integer result will be correct.
    // To handle the signed exponent one alternative would be to take absolute values, saving
    // signs, etc - but that gets a bit complicated with 2-complement integers.
    // Instead, we remember that we don't really want the exact division per se - what we're
    // really after is only rewriting e = 3*n+m. That will actually be *easier* to handle if
    // we require that m must be positive (fewer cases to handle) instead of having n as the
    // strict e/3.
    // To handle this we start by adding 127 to the exponent. This value corresponds to the
    // exponent bias, minus 1 because frexp() has a different standard for the value it returns,
    // but then we add 1 back to handle the extra 0.333 in 21845. So, we have offsetExp = e+127
    // and then multiply by 21845 to get a division result offsetExpDiv3.
    // A (signed) value for n is then recovered by subtracting 42 (bias-1)/3 from k.
    // To calculate a strict remainder we should evaluate offsetExp - 3*offsetExpDiv3 - 1, where
    // the extra 1 corrects for the value we added to the exponent to get correct division.
    // This remainder would have the value 0,1, or 2, but since we only use it to select
    // other numbers we can skip the last step and just handle the cases as 1,2 or 3 instead.
    //
    // OK; end of long detour. Here's how we actually do it in our implementation by using
    // floating-point for the exponent instead to avoid needing integer shifts:
    //
    // 1) Convert the exponent (obtained from frexp) to a float
    // 2) Calculate offsetExp = exp + offset. Note that we should not add the extra 1 here since we
    //    do floating-point division instead of our integer hack, so it's the exponent bias-1, or
    //    the largest exponent minus 2.
    // 3) Divide the float by 3 by multiplying with 1/3
    // 4) Truncate it to an integer to get the division result. This is potentially dangerous in
    //    combination with floating-point, because many integers cannot be represented exactly in
    //    floating point, and if we are just epsilon below the result might be truncated to a lower
    //    integer. I have not observed this on x86, but to have a safety margin we can add a small
    //    fraction - since we already know the fraction part should be either 0, 0.333..., or 0.666...
    //    We can even save this extra floating-point addition by adding a small fraction (0.1) when
    //    we introduce the exponent offset - that will correspond to a safety margin of 0.1/3, which is plenty.
    // 5) Get the remainder part by subtracting the truncated floating-point part.
    //    Here too we will have a plain division, so the remainder is a strict modulus
    //    and will have the values 0, 1 or 2.
    //
    // Before worrying about the few wasted cycles due to longer fp latency, this has the
    // additional advantage that we don't use a single integer operation, so the algorithm
    // will work just A-OK on all SIMD implementations, which avoids diverging code paths.

    // The  0.1 here is the safety margin due to  truncation described in item 4 in the comments above.
    SimdFloat offsetExp = cvtI2R(exponent) + SimdFloat(static_cast<float>(3 * offsetDiv3) + 0.1);

    SimdFloat offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division

    SimdFInt32 expDiv3 = cvtR2I(offsetExpDiv3 - SimdFloat(static_cast<float>(offsetDiv3)));

    SimdFloat remainder = offsetExp - offsetExpDiv3 * three;

    // If remainder is 0 we should just have the factor 1.0,
    // so first pick 1.0 if it is below 0.5, and 2^(1/3) if it's above 0.5 (i.e., 1 or 2)
    SimdFloat factor = blend(one, cbrt2, SimdFloat(0.5) < remainder);
    // Second, we overwrite with 2^(2/3) if rem>1.5 (i.e., 2)
    factor = blend(factor, sqrCbrt2, SimdFloat(1.5) < remainder);

    // Assemble the non-signed fraction, and add the sign back by xor
    SimdFloat fraction = (nom * invDenom * factor) ^ xSignBit;
    // Load to IEEE754 number, and set result to 0.0 if x was 0.0 or denormal
    SimdFloat result = selectByMask(ldexp(fraction, expDiv3), xIsNonZero);

    return result;
}

/*! \brief Inverse cube root for SIMD floats
 *
 * \param x      Argument to calculate cube root of. Can be positive or
 *               negative, but the magnitude cannot be lower than
 *               the smallest normal number.
 * \return       Cube root of x. Undefined for values that don't
 *               fulfill the restriction of abs(x) > minFloat.
 */
static inline SimdFloat gmx_simdcall invcbrt(SimdFloat x)
{
    const SimdFloat signBit(GMX_FLOAT_NEGZERO);
    const SimdFloat minFloat(std::numeric_limits<float>::min());
    // Bias is 128-1 = 127, which is not divisible by 3. Since the largest-magnitude
    // negative exponent from frexp() is -126, we can subtract one more unit to get 126
    // as offset, which is divisible by 3 (result 42). To avoid clang warnings about fragile integer
    // division mixed with FP, we let the divided value (42) be the original constant.
    const std::int32_t offsetDiv3(42);
    const SimdFloat    c2(-0.191502161678719066F);
    const SimdFloat    c1(0.697570460207922770F);
    const SimdFloat    c0(0.492659620528969547F);
    const SimdFloat    one(1.0F);
    const SimdFloat    two(2.0F);
    const SimdFloat    three(3.0F);
    const SimdFloat    oneThird(1.0F / 3.0F);
    const SimdFloat    invCbrt2(1.0F / 1.2599210498948731648F);
    const SimdFloat    invSqrCbrt2(1.0F / 1.5874010519681994748F);

    // We use pretty much exactly the same implementation as for cbrt(x),
    // but to compute the inverse we swap the nominator/denominator
    // in the quotient, and also swap the sign of the exponent parts.

    SimdFloat xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdFloat xAbs     = andNot(signBit, x); // select everthing but the sign bit => abs(x)

    SimdFInt32 exponent;
    SimdFloat  y = frexp(xAbs, &exponent);
    // For the mantissa (y) we will use a limited-range approximation of cbrt(y),
    // by first using a polynomial and then evaluating
    // Transform y to z = c2*y^2 + c1*y + c0, then w = z^3, and finally
    // evaluate the quotient q = z * (w + 2 * y) / (2 * w + y).
    SimdFloat z        = fma(fma(y, c2, c1), y, c0);
    SimdFloat w        = z * z * z;
    SimdFloat nom      = fma(two, w, y);
    SimdFloat invDenom = inv(z * fma(two, y, w));

    // The  0.1 here is the safety margin due to  truncation described in item 4 in the comments above.
    SimdFloat offsetExp = cvtI2R(exponent) + SimdFloat(static_cast<float>(3 * offsetDiv3) + 0.1);
    SimdFloat offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division

    // We should swap the sign here, so we change order of the terms in the subtraction
    SimdFInt32 expDiv3 = cvtR2I(SimdFloat(static_cast<float>(offsetDiv3)) - offsetExpDiv3);

    // Swap sign here too, so remainder is either 0, -1 or -2
    SimdFloat remainder = offsetExpDiv3 * three - offsetExp;

    // If remainder is 0 we should just have the factor 1.0,
    // so first pick 1.0 if it is above -0.5, and 2^(-1/3) if it's below -0.5 (i.e., -1 or -2)
    SimdFloat factor = blend(one, invCbrt2, remainder < SimdFloat(-0.5));
    // Second, we overwrite with 2^(-2/3) if rem<-1.5 (i.e., -2)
    factor = blend(factor, invSqrCbrt2, remainder < SimdFloat(-1.5));

    // Assemble the non-signed fraction, and add the sign back by xor
    SimdFloat fraction = (nom * invDenom * factor) ^ xSignBit;
    // Load to IEEE754 number, and set result to 0.0 if x was 0.0 or denormal
    SimdFloat result = ldexp(fraction, expDiv3);

    return result;
}

/*! \brief SIMD float log2(x). This is the base-2 logarithm.
 *
 * \param x Argument, should be >0.
 * \result The base-2 logarithm of x. Undefined if argument is invalid.
 */
static inline SimdFloat gmx_simdcall log2(SimdFloat x)
{
    // This implementation computes log2 by
    // 1) Extracting the exponent and adding it to...
    // 2) A 9th-order minimax approximation using only odd
    //    terms of (x-1)/(x+1), where x is the mantissa.

#        if GMX_SIMD_HAVE_NATIVE_LOG_FLOAT
    // Just rescale if native log2() is not present, but log() is.
    return log(x) * SimdFloat(std::log2(std::exp(1.0)));
#        else
    const SimdFloat one(1.0F);
    const SimdFloat two(2.0F);
    const SimdFloat invsqrt2(1.0F / std::sqrt(2.0F));
    const SimdFloat CL9(0.342149508897807708152F);
    const SimdFloat CL7(0.411570606888219447939F);
    const SimdFloat CL5(0.577085979152320294183F);
    const SimdFloat CL3(0.961796550607099898222F);
    const SimdFloat CL1(2.885390081777926774009F);
    SimdFloat       fExp, x2, p;
    SimdFBool       m;
    SimdFInt32      iExp;

    x    = frexp(x, &iExp);
    fExp = cvtI2R(iExp);

    m = x < invsqrt2;
    // Adjust to non-IEEE format for x<1/sqrt(2): exponent -= 1, mantissa *= 2.0
    fExp = fExp - selectByMask(one, m);
    x    = x * blend(one, two, m);

    x  = (x - one) * inv(x + one);
    x2 = x * x;

    p = fma(CL9, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, fExp);

    return p;
#        endif
}

#        if !GMX_SIMD_HAVE_NATIVE_LOG_FLOAT
/*! \brief SIMD float log(x). This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
static inline SimdFloat gmx_simdcall log(SimdFloat x)
{
    const SimdFloat one(1.0F);
    const SimdFloat two(2.0F);
    const SimdFloat invsqrt2(1.0F / std::sqrt(2.0F));
    const SimdFloat corr(0.693147180559945286226764F);
    const SimdFloat CL9(0.2371599674224853515625F);
    const SimdFloat CL7(0.285279005765914916992188F);
    const SimdFloat CL5(0.400005519390106201171875F);
    const SimdFloat CL3(0.666666567325592041015625F);
    const SimdFloat CL1(2.0F);
    SimdFloat       fExp, x2, p;
    SimdFBool       m;
    SimdFInt32      iExp;

    x    = frexp(x, &iExp);
    fExp = cvtI2R(iExp);

    m = x < invsqrt2;
    // Adjust to non-IEEE format for x<1/sqrt(2): exponent -= 1, mantissa *= 2.0
    fExp = fExp - selectByMask(one, m);
    x    = x * blend(one, two, m);

    x  = (x - one) * inv(x + one);
    x2 = x * x;

    p = fma(CL9, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, corr * fExp);

    return p;
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_EXP2_FLOAT
/*! \brief SIMD float 2^x
 *
 * \tparam opt If this is changed from the default (safe) into the unsafe
 *             option, input values that would otherwise lead to zero-clamped
 *             results are not allowed and will lead to undefined results.
 *
 * \param x Argument. For the default (safe) function version this can be
 *          arbitrarily small value, but the routine might clamp the result to
 *          zero for arguments that would produce subnormal IEEE754-2008 results.
 *          This corresponds to inputs below -126 in single or -1022 in double,
 *          and it might overflow for arguments reaching 127 (single) or
 *          1023 (double). If you enable the unsafe math optimization,
 *          very small arguments will not necessarily be zero-clamped, but
 *          can produce undefined results.
 *
 * \result 2^x. The result is undefined for very large arguments that cause
 *          internal floating-point overflow. If unsafe optimizations are enabled,
 *          this is also true for very small values.
 *
 * \note    The definition range of this function is just-so-slightly smaller
 *          than the allowed IEEE exponents for many architectures. This is due
 *          to the implementation, which will hopefully improve in the future.
 *
 * \warning You cannot rely on this implementation returning inf for arguments
 *          that cause overflow. If you have some very large
 *          values and need to rely on getting a valid numerical output,
 *          take the minimum of your variable and the largest valid argument
 *          before calling this routine.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall exp2(SimdFloat x)
{
    const SimdFloat CC6(0.0001534581200287996416911311F);
    const SimdFloat CC5(0.001339993121934088894618990F);
    const SimdFloat CC4(0.009618488957115180159497841F);
    const SimdFloat CC3(0.05550328776964726865751735F);
    const SimdFloat CC2(0.2402264689063408646490722F);
    const SimdFloat CC1(0.6931472057372680777553816F);
    const SimdFloat one(1.0F);

    SimdFloat intpart;
    SimdFloat fexppart;
    SimdFloat p;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -127, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer.
    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdFloat(std::numeric_limits<std::int32_t>::lowest()));
    }

    fexppart = ldexp<opt>(one, cvtR2I(x));
    intpart  = round(x);
    x        = x - intpart;

    p = fma(CC6, x, CC5);
    p = fma(p, x, CC4);
    p = fma(p, x, CC3);
    p = fma(p, x, CC2);
    p = fma(p, x, CC1);
    p = fma(p, x, one);
    x = p * fexppart;
    return x;
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_EXP_FLOAT
/*! \brief SIMD float exp(x).
 *
 * In addition to scaling the argument for 2^x this routine correctly does
 * extended precision arithmetics to improve accuracy.
 *
 * \tparam opt If this is changed from the default (safe) into the unsafe
 *             option, input values that would otherwise lead to zero-clamped
 *             results are not allowed and will lead to undefined results.
 *
 * \param x Argument. For the default (safe) function version this can be
 *          arbitrarily small value, but the routine might clamp the result to
 *          zero for arguments that would produce subnormal IEEE754-2008 results.
 *          This corresponds to input arguments reaching
 *          -126*ln(2)=-87.3 in single, or -1022*ln(2)=-708.4 (double).
 *          Similarly, it might overflow for arguments reaching
 *          127*ln(2)=88.0 (single) or 1023*ln(2)=709.1 (double). If the
 *          unsafe math optimizations are enabled, small input values that would
 *          result in zero-clamped output are not allowed.
 *
 * \result exp(x). Overflowing arguments are likely to either return 0 or inf,
 *         depending on the underlying implementation. If unsafe optimizations
 *         are enabled, this is also true for very small values.
 *
 * \note    The definition range of this function is just-so-slightly smaller
 *          than the allowed IEEE exponents for many architectures. This is due
 *          to the implementation, which will hopefully improve in the future.
 *
 * \warning You cannot rely on this implementation returning inf for arguments
 *          that cause overflow. If you have some very large
 *          values and need to rely on getting a valid numerical output,
 *          take the minimum of your variable and the largest valid argument
 *          before calling this routine.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall exp(SimdFloat x)
{
    const SimdFloat argscale(1.44269504088896341F);
    const SimdFloat invargscale0(-0.693145751953125F);
    const SimdFloat invargscale1(-1.428606765330187045e-06F);
    const SimdFloat CC4(0.00136324646882712841033936F);
    const SimdFloat CC3(0.00836596917361021041870117F);
    const SimdFloat CC2(0.0416710823774337768554688F);
    const SimdFloat CC1(0.166665524244308471679688F);
    const SimdFloat CC0(0.499999850988388061523438F);
    const SimdFloat one(1.0F);
    SimdFloat       fexppart;
    SimdFloat       intpart;
    SimdFloat       y, p;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -127, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer
    // after being multiplied by argscale.

    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdFloat(std::numeric_limits<std::int32_t>::lowest()) / argscale);
    }

    y = x * argscale;


    fexppart = ldexp<opt>(one, cvtR2I(y));
    intpart  = round(y);

    // Extended precision arithmetics
    x = fma(invargscale0, intpart, x);
    x = fma(invargscale1, intpart, x);

    p = fma(CC4, x, CC3);
    p = fma(p, x, CC2);
    p = fma(p, x, CC1);
    p = fma(p, x, CC0);
    p = fma(x * x, p, x);
#            if GMX_SIMD_HAVE_FMA
    x = fma(p, fexppart, fexppart);
#            else
    x = (p + one) * fexppart;
#            endif
    return x;
}
#        endif

/*! \brief SIMD float pow(x,y)
 *
 * This returns x^y for SIMD values.
 *
 * \tparam opt If this is changed from the default (safe) into the unsafe
 *             option, there are no guarantees about correct results for x==0.
 *
 * \param x Base.
 *
 * \param y exponent.

 * \result x^y. Overflowing arguments are likely to either return 0 or inf,
 *         depending on the underlying implementation. If unsafe optimizations
 *         are enabled, this is also true for x==0.
 *
 * \warning You cannot rely on this implementation returning inf for arguments
 *          that cause overflow. If you have some very large
 *          values and need to rely on getting a valid numerical output,
 *          take the minimum of your variable and the largest valid argument
 *          before calling this routine.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall pow(SimdFloat x, SimdFloat y)
{
    SimdFloat xcorr;

    if (opt == MathOptimization::Safe)
    {
        xcorr = max(x, SimdFloat(std::numeric_limits<float>::min()));
    }
    else
    {
        xcorr = x;
    }

    SimdFloat result = exp2<opt>(y * log2(xcorr));

    if (opt == MathOptimization::Safe)
    {
        // if x==0 and y>0 we explicitly set the result to 0.0
        // For any x with y==0, the result will already be 1.0 since we multiply by y (0.0) and call exp().
        result = blend(result, setZero(), x == setZero() && setZero() < y);
    }

    return result;
}


/*! \brief SIMD float erf(x).
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to full precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdFloat gmx_simdcall erf(SimdFloat x)
{
    // Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1]
    const SimdFloat CA6(7.853861353153693e-5F);
    const SimdFloat CA5(-8.010193625184903e-4F);
    const SimdFloat CA4(5.188327685732524e-3F);
    const SimdFloat CA3(-2.685381193529856e-2F);
    const SimdFloat CA2(1.128358514861418e-1F);
    const SimdFloat CA1(-3.761262582423300e-1F);
    const SimdFloat CA0(1.128379165726710F);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2]
    const SimdFloat CB9(-0.0018629930017603923F);
    const SimdFloat CB8(0.003909821287598495F);
    const SimdFloat CB7(-0.0052094582210355615F);
    const SimdFloat CB6(0.005685614362160572F);
    const SimdFloat CB5(-0.0025367682853477272F);
    const SimdFloat CB4(-0.010199799682318782F);
    const SimdFloat CB3(0.04369575504816542F);
    const SimdFloat CB2(-0.11884063474674492F);
    const SimdFloat CB1(0.2732120154030589F);
    const SimdFloat CB0(0.42758357702025784F);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19]
    const SimdFloat CC10(-0.0445555913112064F);
    const SimdFloat CC9(0.21376355144663348F);
    const SimdFloat CC8(-0.3473187200259257F);
    const SimdFloat CC7(0.016690861551248114F);
    const SimdFloat CC6(0.7560973182491192F);
    const SimdFloat CC5(-1.2137903600145787F);
    const SimdFloat CC4(0.8411872321232948F);
    const SimdFloat CC3(-0.08670413896296343F);
    const SimdFloat CC2(-0.27124782687240334F);
    const SimdFloat CC1(-0.0007502488047806069F);
    const SimdFloat CC0(0.5642114853803148F);
    const SimdFloat one(1.0F);
    const SimdFloat two(2.0F);

    SimdFloat x2, x4, y;
    SimdFloat t, t2, w, w2;
    SimdFloat pA0, pA1, pB0, pB1, pC0, pC1;
    SimdFloat expmx2;
    SimdFloat res_erf, res_erfc, res;
    SimdFBool m, maskErf;

    // Calculate erf()
    x2 = x * x;
    x4 = x2 * x2;

    pA0 = fma(CA6, x4, CA4);
    pA1 = fma(CA5, x4, CA3);
    pA0 = fma(pA0, x4, CA2);
    pA1 = fma(pA1, x4, CA1);
    pA0 = pA0 * x4;
    pA0 = fma(pA1, x2, pA0);
    // Constant term must come last for precision reasons
    pA0 = pA0 + CA0;

    res_erf = x * pA0;

    // Calculate erfc
    y       = abs(x);
    maskErf = SimdFloat(0.75F) <= y;
    t       = maskzInv(y, maskErf);
    w       = t - one;
    t2      = t * t;
    w2      = w * w;

    // No need for a floating-point sieve here (as in erfc), since erf()
    // will never return values that are extremely small for large args.
    expmx2 = exp(-y * y);

    pB1 = fma(CB9, w2, CB7);
    pB0 = fma(CB8, w2, CB6);
    pB1 = fma(pB1, w2, CB5);
    pB0 = fma(pB0, w2, CB4);
    pB1 = fma(pB1, w2, CB3);
    pB0 = fma(pB0, w2, CB2);
    pB1 = fma(pB1, w2, CB1);
    pB0 = fma(pB0, w2, CB0);
    pB0 = fma(pB1, w, pB0);

    pC0 = fma(CC10, t2, CC8);
    pC1 = fma(CC9, t2, CC7);
    pC0 = fma(pC0, t2, CC6);
    pC1 = fma(pC1, t2, CC5);
    pC0 = fma(pC0, t2, CC4);
    pC1 = fma(pC1, t2, CC3);
    pC0 = fma(pC0, t2, CC2);
    pC1 = fma(pC1, t2, CC1);

    pC0 = fma(pC0, t2, CC0);
    pC0 = fma(pC1, t, pC0);
    pC0 = pC0 * t;

    // Select pB0 or pC0 for erfc()
    m        = two < y;
    res_erfc = blend(pB0, pC0, m);
    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    m        = x < setZero();
    res_erfc = blend(res_erfc, two - res_erfc, m);

    // Select erf() or erfc()
    res = blend(res_erf, one - res_erfc, maskErf);

    return res;
}

/*! \brief SIMD float erfc(x).
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves full precision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-30 for single precision)
 * since that is not relevant for MD.
 */
static inline SimdFloat gmx_simdcall erfc(SimdFloat x)
{
    // Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1]
    const SimdFloat CA6(7.853861353153693e-5F);
    const SimdFloat CA5(-8.010193625184903e-4F);
    const SimdFloat CA4(5.188327685732524e-3F);
    const SimdFloat CA3(-2.685381193529856e-2F);
    const SimdFloat CA2(1.128358514861418e-1F);
    const SimdFloat CA1(-3.761262582423300e-1F);
    const SimdFloat CA0(1.128379165726710F);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2]
    const SimdFloat CB9(-0.0018629930017603923F);
    const SimdFloat CB8(0.003909821287598495F);
    const SimdFloat CB7(-0.0052094582210355615F);
    const SimdFloat CB6(0.005685614362160572F);
    const SimdFloat CB5(-0.0025367682853477272F);
    const SimdFloat CB4(-0.010199799682318782F);
    const SimdFloat CB3(0.04369575504816542F);
    const SimdFloat CB2(-0.11884063474674492F);
    const SimdFloat CB1(0.2732120154030589F);
    const SimdFloat CB0(0.42758357702025784F);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19]
    const SimdFloat CC10(-0.0445555913112064F);
    const SimdFloat CC9(0.21376355144663348F);
    const SimdFloat CC8(-0.3473187200259257F);
    const SimdFloat CC7(0.016690861551248114F);
    const SimdFloat CC6(0.7560973182491192F);
    const SimdFloat CC5(-1.2137903600145787F);
    const SimdFloat CC4(0.8411872321232948F);
    const SimdFloat CC3(-0.08670413896296343F);
    const SimdFloat CC2(-0.27124782687240334F);
    const SimdFloat CC1(-0.0007502488047806069F);
    const SimdFloat CC0(0.5642114853803148F);
    // Coefficients for expansion of exp(x) in [0,0.1]
    // CD0 and CD1 are both 1.0, so no need to declare them separately
    const SimdFloat CD2(0.5000066608081202F);
    const SimdFloat CD3(0.1664795422874624F);
    const SimdFloat CD4(0.04379839977652482F);
    const SimdFloat one(1.0F);
    const SimdFloat two(2.0F);

    /* We need to use a small trick here, since we cannot assume all SIMD
     * architectures support integers, and the flag we want (0xfffff000) would
     * evaluate to NaN (i.e., it cannot be expressed as a floating-point num).
     * Instead, we represent the flags 0xf0f0f000 and 0x0f0f0000 as valid
     * fp numbers, and perform a logical or. Since the expression is constant,
     * we can at least hope it is evaluated at compile-time.
     */
#        if GMX_SIMD_HAVE_LOGICAL
    const SimdFloat sieve(SimdFloat(-5.965323564e+29F) | SimdFloat(7.05044434e-30F));
#        else
    const int                         isieve = 0xFFFFF000;
    alignas(GMX_SIMD_ALIGNMENT) float mem[GMX_SIMD_FLOAT_WIDTH];

    union {
        float f;
        int   i;
    } conv;
    int i;
#        endif

    SimdFloat x2, x4, y;
    SimdFloat q, z, t, t2, w, w2;
    SimdFloat pA0, pA1, pB0, pB1, pC0, pC1;
    SimdFloat expmx2, corr;
    SimdFloat res_erf, res_erfc, res;
    SimdFBool m, msk_erf;

    // Calculate erf()
    x2 = x * x;
    x4 = x2 * x2;

    pA0 = fma(CA6, x4, CA4);
    pA1 = fma(CA5, x4, CA3);
    pA0 = fma(pA0, x4, CA2);
    pA1 = fma(pA1, x4, CA1);
    pA1 = pA1 * x2;
    pA0 = fma(pA0, x4, pA1);
    // Constant term must come last for precision reasons
    pA0 = pA0 + CA0;

    res_erf = x * pA0;

    // Calculate erfc
    y       = abs(x);
    msk_erf = SimdFloat(0.75F) <= y;
    t       = maskzInv(y, msk_erf);
    w       = t - one;
    t2      = t * t;
    w2      = w * w;
    /*
     * We cannot simply calculate exp(-y2) directly in single precision, since
     * that will lose a couple of bits of precision due to the multiplication.
     * Instead, we introduce y=z+w, where the last 12 bits of precision are in w.
     * Then we get exp(-y2) = exp(-z2)*exp((z-y)*(z+y)).
     *
     * The only drawback with this is that it requires TWO separate exponential
     * evaluations, which would be horrible performance-wise. However, the argument
     * for the second exp() call is always small, so there we simply use a
     * low-order minimax expansion on [0,0.1].
     *
     * However, this neat idea requires support for logical ops (and) on
     * FP numbers, which some vendors decided isn't necessary in their SIMD
     * instruction sets (Hi, IBM VSX!). In principle we could use some tricks
     * in double, but we still need memory as a backup when that is not available,
     * and this case is rare enough that we go directly there...
     */
#        if GMX_SIMD_HAVE_LOGICAL
    z = y & sieve;
#        else
    store(mem, y);
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = mem[i];
        conv.i = conv.i & isieve;
        mem[i] = conv.f;
    }
    z = load<SimdFloat>(mem);
#        endif
    q    = (z - y) * (z + y);
    corr = fma(CD4, q, CD3);
    corr = fma(corr, q, CD2);
    corr = fma(corr, q, one);
    corr = fma(corr, q, one);

    expmx2 = exp(-z * z);
    expmx2 = expmx2 * corr;

    pB1 = fma(CB9, w2, CB7);
    pB0 = fma(CB8, w2, CB6);
    pB1 = fma(pB1, w2, CB5);
    pB0 = fma(pB0, w2, CB4);
    pB1 = fma(pB1, w2, CB3);
    pB0 = fma(pB0, w2, CB2);
    pB1 = fma(pB1, w2, CB1);
    pB0 = fma(pB0, w2, CB0);
    pB0 = fma(pB1, w, pB0);

    pC0 = fma(CC10, t2, CC8);
    pC1 = fma(CC9, t2, CC7);
    pC0 = fma(pC0, t2, CC6);
    pC1 = fma(pC1, t2, CC5);
    pC0 = fma(pC0, t2, CC4);
    pC1 = fma(pC1, t2, CC3);
    pC0 = fma(pC0, t2, CC2);
    pC1 = fma(pC1, t2, CC1);

    pC0 = fma(pC0, t2, CC0);
    pC0 = fma(pC1, t, pC0);
    pC0 = pC0 * t;

    // Select pB0 or pC0 for erfc()
    m        = two < y;
    res_erfc = blend(pB0, pC0, m);
    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    m        = x < setZero();
    res_erfc = blend(res_erfc, two - res_erfc, m);

    // Select erf() or erfc()
    res = blend(one - res_erf, res_erfc, msk_erf);

    return res;
}

/*! \brief SIMD float sin \& cos.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * This version achieves close to machine precision, but for very large
 * magnitudes of the argument we inherently begin to lose accuracy due to the
 * argument reduction, despite using extended precision arithmetics internally.
 */
static inline void gmx_simdcall sincos(SimdFloat x, SimdFloat* sinval, SimdFloat* cosval)
{
    // Constants to subtract Pi/4*x from y while minimizing precision loss
    const SimdFloat argred0(-1.5703125);
    const SimdFloat argred1(-4.83751296997070312500e-04F);
    const SimdFloat argred2(-7.54953362047672271729e-08F);
    const SimdFloat argred3(-2.56334406825708960298e-12F);
    const SimdFloat two_over_pi(static_cast<float>(2.0F / M_PI));
    const SimdFloat const_sin2(-1.9515295891e-4F);
    const SimdFloat const_sin1(8.3321608736e-3F);
    const SimdFloat const_sin0(-1.6666654611e-1F);
    const SimdFloat const_cos2(2.443315711809948e-5F);
    const SimdFloat const_cos1(-1.388731625493765e-3F);
    const SimdFloat const_cos0(4.166664568298827e-2F);
    const SimdFloat half(0.5F);
    const SimdFloat one(1.0F);
    SimdFloat       ssign, csign;
    SimdFloat       x2, y, z, psin, pcos, sss, ccc;
    SimdFBool       m;

#        if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdFInt32 ione(1);
    const SimdFInt32 itwo(2);
    SimdFInt32       iy;

    z  = x * two_over_pi;
    iy = cvtR2I(z);
    y  = round(z);

    m     = cvtIB2B((iy & ione) == SimdFInt32(0));
    ssign = selectByMask(SimdFloat(GMX_FLOAT_NEGZERO), cvtIB2B((iy & itwo) == itwo));
    csign = selectByMask(SimdFloat(GMX_FLOAT_NEGZERO), cvtIB2B(((iy + ione) & itwo) == itwo));
#        else
    const SimdFloat quarter(0.25f);
    const SimdFloat minusquarter(-0.25f);
    SimdFloat       q;
    SimdFBool       m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign = x;
    x     = abs(x);
    // It is critical that half-way cases are rounded down
    z = fma(x, two_over_pi, half);
    y = trunc(z);
    q = z * quarter;
    q = q - trunc(q);
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q = q - half;
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#            if GMX_SIMD_HAVE_LOGICAL
    ssign = ssign & SimdFloat(GMX_FLOAT_NEGZERO);
    csign = andNot(q, SimdFloat(GMX_FLOAT_NEGZERO));
    ssign = ssign ^ csign;
#            else
    ssign = copysign(SimdFloat(1.0f), ssign);
    csign = copysign(SimdFloat(1.0f), q);
    csign = -csign;
    ssign = ssign * csign; // swap ssign if csign was set.
#            endif
    // Check if y had value 1 or 3 (remember we subtracted 0.5 from q)
    m1 = (q < minusquarter);
    m2 = (setZero() <= q);
    m3 = (q < quarter);
    m2 = m2 && m3;
    m  = m1 || m2;
    // where mask is FALSE, swap sign.
    csign   = csign * blend(SimdFloat(-1.0f), one, m);
#        endif
    x  = fma(y, argred0, x);
    x  = fma(y, argred1, x);
    x  = fma(y, argred2, x);
    x  = fma(y, argred3, x);
    x2 = x * x;

    psin = fma(const_sin2, x2, const_sin1);
    psin = fma(psin, x2, const_sin0);
    psin = fma(psin, x * x2, x);
    pcos = fma(const_cos2, x2, const_cos1);
    pcos = fma(pcos, x2, const_cos0);
    pcos = fms(pcos, x2, half);
    pcos = fma(pcos, x2, one);

    sss = blend(pcos, psin, m);
    ccc = blend(psin, pcos, m);
    // See comment for GMX_SIMD_HAVE_LOGICAL section above.
#        if GMX_SIMD_HAVE_LOGICAL
    *sinval = sss ^ ssign;
    *cosval = ccc ^ csign;
#        else
    *sinval = sss * ssign;
    *cosval = ccc * csign;
#        endif
}

/*! \brief SIMD float sin(x).
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall sin(SimdFloat x)
{
    SimdFloat s, c;
    sincos(x, &s, &c);
    return s;
}

/*! \brief SIMD float cos(x).
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall cos(SimdFloat x)
{
    SimdFloat s, c;
    sincos(x, &s, &c);
    return c;
}

/*! \brief SIMD float tan(x).
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdFloat gmx_simdcall tan(SimdFloat x)
{
    const SimdFloat argred0(-1.5703125);
    const SimdFloat argred1(-4.83751296997070312500e-04F);
    const SimdFloat argred2(-7.54953362047672271729e-08F);
    const SimdFloat argred3(-2.56334406825708960298e-12F);
    const SimdFloat two_over_pi(static_cast<float>(2.0F / M_PI));
    const SimdFloat CT6(0.009498288995810566122993911);
    const SimdFloat CT5(0.002895755790837379295226923);
    const SimdFloat CT4(0.02460087336161924491836265);
    const SimdFloat CT3(0.05334912882656359828045988);
    const SimdFloat CT2(0.1333989091464957704418495);
    const SimdFloat CT1(0.3333307599244198227797507);

    SimdFloat x2, p, y, z;
    SimdFBool m;

#        if GMX_SIMD_HAVE_FINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdFInt32 iy;
    SimdFInt32 ione(1);

    z  = x * two_over_pi;
    iy = cvtR2I(z);
    y  = round(z);
    m  = cvtIB2B((iy & ione) == ione);

    x = fma(y, argred0, x);
    x = fma(y, argred1, x);
    x = fma(y, argred2, x);
    x = fma(y, argred3, x);
    x = selectByMask(SimdFloat(GMX_FLOAT_NEGZERO), m) ^ x;
#        else
    const SimdFloat quarter(0.25f);
    const SimdFloat half(0.5f);
    const SimdFloat threequarter(0.75f);
    SimdFloat       w, q;
    SimdFBool       m1, m2, m3;

    w     = abs(x);
    z     = fma(w, two_over_pi, half);
    y     = trunc(z);
    q     = z * quarter;
    q     = q - trunc(q);
    m1    = quarter <= q;
    m2    = q < half;
    m3    = threequarter <= q;
    m1    = m1 && m2;
    m     = m1 || m3;
    w     = fma(y, argred0, w);
    w     = fma(y, argred1, w);
    w     = fma(y, argred2, w);
    w     = fma(y, argred3, w);
    w     = blend(w, -w, m);
    x     = w * copysign(SimdFloat(1.0), x);
#        endif
    x2 = x * x;
    p  = fma(CT6, x2, CT5);
    p  = fma(p, x2, CT4);
    p  = fma(p, x2, CT3);
    p  = fma(p, x2, CT2);
    p  = fma(p, x2, CT1);
    p  = fma(x2 * p, x, x);

    p = blend(p, maskzInv(p, m), m);
    return p;
}

/*! \brief SIMD float asin(x).
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdFloat gmx_simdcall asin(SimdFloat x)
{
    const SimdFloat limitlow(1e-4F);
    const SimdFloat half(0.5F);
    const SimdFloat one(1.0F);
    const SimdFloat halfpi(static_cast<float>(M_PI / 2.0F));
    const SimdFloat CC5(4.2163199048E-2F);
    const SimdFloat CC4(2.4181311049E-2F);
    const SimdFloat CC3(4.5470025998E-2F);
    const SimdFloat CC2(7.4953002686E-2F);
    const SimdFloat CC1(1.6666752422E-1F);
    SimdFloat       xabs;
    SimdFloat       z, z1, z2, q, q1, q2;
    SimdFloat       pA, pB;
    SimdFBool       m, m2;

    xabs = abs(x);
    m    = half < xabs;
    z1   = half * (one - xabs);
    m2   = xabs < one;
    q1   = z1 * maskzInvsqrt(z1, m2);
    q2   = xabs;
    z2   = q2 * q2;
    z    = blend(z2, z1, m);
    q    = blend(q2, q1, m);

    z2 = z * z;
    pA = fma(CC5, z2, CC3);
    pB = fma(CC4, z2, CC2);
    pA = fma(pA, z2, CC1);
    pA = pA * z;
    z  = fma(pB, z2, pA);
    z  = fma(z, q, q);
    q2 = halfpi - z;
    q2 = q2 - z;
    z  = blend(z, q2, m);

    m = limitlow < xabs;
    z = blend(xabs, z, m);
    z = copysign(z, x);

    return z;
}

/*! \brief SIMD float acos(x).
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdFloat gmx_simdcall acos(SimdFloat x)
{
    const SimdFloat one(1.0F);
    const SimdFloat half(0.5F);
    const SimdFloat pi(static_cast<float>(M_PI));
    const SimdFloat halfpi(static_cast<float>(M_PI / 2.0F));
    SimdFloat       xabs;
    SimdFloat       z, z1, z2, z3;
    SimdFBool       m1, m2, m3;

    xabs = abs(x);
    m1   = half < xabs;
    m2   = setZero() < x;

    z  = fnma(half, xabs, half);
    m3 = xabs < one;
    z  = z * maskzInvsqrt(z, m3);
    z  = blend(x, z, m1);
    z  = asin(z);

    z2 = z + z;
    z1 = pi - z2;
    z3 = halfpi - z;
    z  = blend(z1, z2, m2);
    z  = blend(z3, z, m1);

    return z;
}

/*! \brief SIMD float asin(x).
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdFloat gmx_simdcall atan(SimdFloat x)
{
    const SimdFloat halfpi(static_cast<float>(M_PI / 2.0F));
    const SimdFloat CA17(0.002823638962581753730774F);
    const SimdFloat CA15(-0.01595690287649631500244F);
    const SimdFloat CA13(0.04250498861074447631836F);
    const SimdFloat CA11(-0.07489009201526641845703F);
    const SimdFloat CA9(0.1063479334115982055664F);
    const SimdFloat CA7(-0.1420273631811141967773F);
    const SimdFloat CA5(0.1999269574880599975585F);
    const SimdFloat CA3(-0.3333310186862945556640F);
    const SimdFloat one(1.0F);
    SimdFloat       x2, x3, x4, pA, pB;
    SimdFBool       m, m2;

    m  = x < setZero();
    x  = abs(x);
    m2 = one < x;
    x  = blend(x, maskzInv(x, m2), m2);

    x2 = x * x;
    x3 = x2 * x;
    x4 = x2 * x2;
    pA = fma(CA17, x4, CA13);
    pB = fma(CA15, x4, CA11);
    pA = fma(pA, x4, CA9);
    pB = fma(pB, x4, CA7);
    pA = fma(pA, x4, CA5);
    pB = fma(pB, x4, CA3);
    pA = fma(pA, x2, pB);
    pA = fma(pA, x3, x);

    pA = blend(pA, halfpi - pA, m2);
    pA = blend(pA, -pA, m);

    return pA;
}

/*! \brief SIMD float atan2(y,x).
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static inline SimdFloat gmx_simdcall atan2(SimdFloat y, SimdFloat x)
{
    const SimdFloat pi(static_cast<float>(M_PI));
    const SimdFloat halfpi(static_cast<float>(M_PI / 2.0));
    SimdFloat       xinv, p, aoffset;
    SimdFBool       mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = x != setZero();
    mask_ynz  = y != setZero();
    mask_xlt0 = x < setZero();
    mask_ylt0 = y < setZero();

    aoffset = selectByNotMask(halfpi, mask_xnz);
    aoffset = selectByMask(aoffset, mask_ynz);

    aoffset = blend(aoffset, pi, mask_xlt0);
    aoffset = blend(aoffset, -aoffset, mask_ylt0);

    xinv = maskzInv(x, mask_xnz);
    p    = y * xinv;
    p    = atan(p);
    p    = p + aoffset;

    return p;
}

/*! \brief Calculate the force correction due to PME analytically in SIMD float.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb force - see below for details.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be \f$ \mbox{erfc}(\beta r)/r\f$, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of \f$(\beta r)\f$ that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since \f$\mbox{erfc}(z)=1-\mbox{erf}(z)\f$, i.e.
 * \f[
 * V = \frac{1}{r} - \frac{\mbox{erf}(\beta r)}{r}
 * \f]
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument \f$beta r\f$ will be in the range 0.15 to ~4, which is
 * the range used for the minimax fit. Use your favorite plotting program
 * to realize how well-behaved \f$\frac{\mbox{erf}(z)}{z}\f$ is in this range!
 *
 * We approximate \f$f(z)=\mbox{erf}(z)/z\f$ with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate \f$f(z)/z\f$ and
 * then only use even powers. This is another minor optimization, since
 * we actually \a want \f$f(z)/z\f$, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=(\beta r)^2\f$.
 * 3. Evaluate this routine with \f$z^2\f$ as the argument.
 * 4. The return value is the expression:
 *
 * \f[
 *    \frac{2 \exp{-z^2}}{\sqrt{\pi} z^2}-\frac{\mbox{erf}(z)}{z^3}
 * \f]
 *
 * 5. Multiply the entire expression by \f$\beta^3\f$. This will get you
 *
 *  \f[
 *    \frac{2 \beta^3 \exp(-z^2)}{\sqrt{\pi} z^2} - \frac{\beta^3 \mbox{erf}(z)}{z^3}
 *  \f]
 *
 *    or, switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{2 \beta \exp(-r^2 \beta^2)}{\sqrt{\pi} r^2} - \frac{\mbox{erf}(r \beta)}{r^3}
 *  \f]
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly
 *
 *  \f[
 *   \frac{\frac{d}{dr}\left( \frac{\mbox{erf}(\beta r)}{r} \right)}{r}
 *  \f]
 *
 * 6. Add the result to \f$r^{-3}\f$, multiply by the product of the charges,
 *    and you have your force (divided by \f$r\f$). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 * This approximation achieves an error slightly lower than 1e-6
 * in single precision and 1e-11 in double precision
 * for arguments smaller than 16 (\f$\beta r \leq 4 \f$);
 * when added to \f$1/r\f$ the error will be insignificant.
 * For \f$\beta r \geq 7206\f$ the return value can be inf or NaN.
 *
 */
static inline SimdFloat gmx_simdcall pmeForceCorrection(SimdFloat z2)
{
    const SimdFloat FN6(-1.7357322914161492954e-8F);
    const SimdFloat FN5(1.4703624142580877519e-6F);
    const SimdFloat FN4(-0.000053401640219807709149F);
    const SimdFloat FN3(0.0010054721316683106153F);
    const SimdFloat FN2(-0.019278317264888380590F);
    const SimdFloat FN1(0.069670166153766424023F);
    const SimdFloat FN0(-0.75225204789749321333F);

    const SimdFloat FD4(0.0011193462567257629232F);
    const SimdFloat FD3(0.014866955030185295499F);
    const SimdFloat FD2(0.11583842382862377919F);
    const SimdFloat FD1(0.50736591960530292870F);
    const SimdFloat FD0(1.0F);

    SimdFloat z4;
    SimdFloat polyFN0, polyFN1, polyFD0, polyFD1;

    z4 = z2 * z2;

    polyFD0 = fma(FD4, z4, FD2);
    polyFD1 = fma(FD3, z4, FD1);
    polyFD0 = fma(polyFD0, z4, FD0);
    polyFD0 = fma(polyFD1, z2, polyFD0);

    polyFD0 = inv(polyFD0);

    polyFN0 = fma(FN6, z4, FN4);
    polyFN1 = fma(FN5, z4, FN3);
    polyFN0 = fma(polyFN0, z4, FN2);
    polyFN1 = fma(polyFN1, z4, FN1);
    polyFN0 = fma(polyFN0, z4, FN0);
    polyFN0 = fma(polyFN1, z2, polyFN0);

    return polyFN0 * polyFD0;
}


/*! \brief Calculate the potential correction due to PME analytically in SIMD float.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * See \ref pmeForceCorrection for details about the approximation.
 *
 * This routine calculates \f$\mbox{erf}(z)/z\f$, although you should provide \f$z^2\f$
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=\beta^2*r^2\f$.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *  \f[
 *   \frac{\mbox{erf}(z)}{z}
 *  \f]
 *
 * 5. Multiply the entire expression by beta and switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{\mbox{erf}(r \beta)}{r}
 *  \f]
 *
 * 6. Subtract the result from \f$1/r\f$, multiply by the product of the charges,
 *    and you have your potential.
 *
 * This approximation achieves an error slightly lower than 1e-6
 * in single precision and 4e-11 in double precision
 * for arguments smaller than 16 (\f$ 0.15 \leq \beta r \leq 4 \f$);
 * for \f$ \beta r \leq 0.15\f$ the error can be twice as high;
 * when added to \f$1/r\f$ the error will be insignificant.
 * For \f$\beta r \geq 7142\f$ the return value can be inf or NaN.
 */
static inline SimdFloat gmx_simdcall pmePotentialCorrection(SimdFloat z2)
{
    const SimdFloat VN6(1.9296833005951166339e-8F);
    const SimdFloat VN5(-1.4213390571557850962e-6F);
    const SimdFloat VN4(0.000041603292906656984871F);
    const SimdFloat VN3(-0.00013134036773265025626F);
    const SimdFloat VN2(0.038657983986041781264F);
    const SimdFloat VN1(0.11285044772717598220F);
    const SimdFloat VN0(1.1283802385263030286F);

    const SimdFloat VD3(0.0066752224023576045451F);
    const SimdFloat VD2(0.078647795836373922256F);
    const SimdFloat VD1(0.43336185284710920150F);
    const SimdFloat VD0(1.0F);

    SimdFloat z4;
    SimdFloat polyVN0, polyVN1, polyVD0, polyVD1;

    z4 = z2 * z2;

    polyVD1 = fma(VD3, z4, VD1);
    polyVD0 = fma(VD2, z4, VD0);
    polyVD0 = fma(polyVD1, z2, polyVD0);

    polyVD0 = inv(polyVD0);

    polyVN0 = fma(VN6, z4, VN4);
    polyVN1 = fma(VN5, z4, VN3);
    polyVN0 = fma(polyVN0, z4, VN2);
    polyVN1 = fma(polyVN1, z4, VN1);
    polyVN0 = fma(polyVN0, z4, VN0);
    polyVN0 = fma(polyVN1, z2, polyVN0);

    return polyVN0 * polyVD0;
}
#    endif

/*! \} */

#    if GMX_SIMD_HAVE_DOUBLE


/*! \name Double precision SIMD math functions
 *
 *  \note In most cases you should use the real-precision functions instead.
 *  \{
 */

/****************************************
 * DOUBLE PRECISION SIMD MATH FUNCTIONS *
 ****************************************/

#        if !GMX_SIMD_HAVE_NATIVE_COPYSIGN_DOUBLE
/*! \brief Composes floating point value with the magnitude of x and the sign of y.
 *
 * \param x Values to set sign for
 * \param y Values used to set sign
 * \return  Magnitude of x, sign of y
 */
static inline SimdDouble gmx_simdcall copysign(SimdDouble x, SimdDouble y)
{
#            if GMX_SIMD_HAVE_LOGICAL
    return abs(x) | (SimdDouble(GMX_DOUBLE_NEGZERO) & y);
#            else
    return blend(abs(x), -abs(x), (y < setZero()));
#            endif
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_RSQRT_ITER_DOUBLE
/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD double.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdDouble gmx_simdcall rsqrtIter(SimdDouble lu, SimdDouble x)
{
    SimdDouble tmp1 = x * lu;
    SimdDouble tmp2 = SimdDouble(-0.5) * lu;
    tmp1            = fma(tmp1, lu, SimdDouble(-3.0));
    return tmp1 * tmp2;
}
#        endif

/*! \brief Calculate 1/sqrt(x) for SIMD double.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall invsqrt(SimdDouble x)
{
    SimdDouble lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles.
 *
 * \param x0  First set of arguments, x0 must be in single range (see below).
 * \param x1  Second set of arguments, x1 must be in single range (see below).
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 *
 * \note Both arguments must be larger than GMX_FLOAT_MIN and smaller than
 *       GMX_FLOAT_MAX, i.e. within the range of single precision.
 *       For the single precision implementation this is obviously always
 *       true for positive values, but for double precision it adds an
 *       extra restriction since the first lookup step might have to be
 *       performed in single precision on some architectures. Note that the
 *       responsibility for checking falls on you - this routine does not
 *       check arguments.
 */
static inline void gmx_simdcall invsqrtPair(SimdDouble x0, SimdDouble x1, SimdDouble* out0, SimdDouble* out1)
{
#        if GMX_SIMD_HAVE_FLOAT && (GMX_SIMD_FLOAT_WIDTH == 2 * GMX_SIMD_DOUBLE_WIDTH) \
                && (GMX_SIMD_RSQRT_BITS < 22)
    SimdFloat  xf  = cvtDD2F(x0, x1);
    SimdFloat  luf = rsqrt(xf);
    SimdDouble lu0, lu1;
    // Intermediate target is single - mantissa+1 bits
#            if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
#            if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
#            if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
    cvtF2DD(luf, &lu0, &lu1);
    // Last iteration(s) performed in double - if we had 22 bits, this gets us to 44 (~1e-15)
#            if (GMX_SIMD_ACCURACY_BITS_SINGLE < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = rsqrtIter(lu0, x0);
    lu1 = rsqrtIter(lu1, x1);
#            endif
#            if (GMX_SIMD_ACCURACY_BITS_SINGLE * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu0 = rsqrtIter(lu0, x0);
    lu1 = rsqrtIter(lu1, x1);
#            endif
    *out0 = lu0;
    *out1 = lu1;
#        else
    *out0 = invsqrt(x0);
    *out1 = invsqrt(x1);
#        endif
}

#        if !GMX_SIMD_HAVE_NATIVE_RCP_ITER_DOUBLE
/*! \brief Perform one Newton-Raphson iteration to improve 1/x for SIMD double.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the reciprocal.
 *
 *  \param lu Approximation of 1/x, typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/x.
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline SimdDouble gmx_simdcall rcpIter(SimdDouble lu, SimdDouble x)
{
    return lu * fnma(lu, x, SimdDouble(2.0));
}
#        endif

/*! \brief Calculate 1/x for SIMD double.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall inv(SimdDouble x)
{
    SimdDouble lu = rcp(x);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}

/*! \brief Division for SIMD doubles
 *
 * \param nom    Nominator
 * \param denom  Denominator, with magnitude in range (GMX_FLOAT_MIN,GMX_FLOAT_MAX).
 *               For single precision this is equivalent to a nonzero argument,
 *               but in double precision it adds an extra restriction since
 *               the first lookup step might have to be performed in single
 *               precision on some architectures. Note that the responsibility
 *               for checking falls on you - this routine does not check arguments.
 *
 * \return nom/denom
 *
 * \note This function does not use any masking to avoid problems with
 *       zero values in the denominator.
 */
static inline SimdDouble gmx_simdcall operator/(SimdDouble nom, SimdDouble denom)
{
    return nom * inv(denom);
}


/*! \brief Calculate 1/sqrt(x) for masked entries of SIMD double.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX for masked-in entries.
 *           See \ref invsqrt for the discussion about argument restrictions.
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdDouble maskzInvsqrt(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = maskzRsqrt(x, m);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate 1/x for SIMD double, masked version.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX for masked-in entries.
 *           See \ref invsqrt for the discussion about argument restrictions.
 *  \param m Mask
 *  \return 1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdDouble gmx_simdcall maskzInv(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = maskzRcp(x, m);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}


/*! \brief Calculate sqrt(x) for SIMD doubles.
 *
 *  \copydetails sqrt(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall sqrt(SimdDouble x)
{
    if (opt == MathOptimization::Safe)
    {
        // As we might use a float version of rsqrt, we mask out small values
        SimdDouble res = maskzInvsqrt(x, SimdDouble(GMX_FLOAT_MIN) < x);
        return res * x;
    }
    else
    {
        return x * invsqrt(x);
    }
}

/*! \brief Cube root for SIMD doubles
 *
 * \param x      Argument to calculate cube root of. Can be negative or zero,
 *               but NaN or Inf values are not supported. Denormal values will
 *               be treated as 0.0.
 * \return       Cube root of x.
 */
static inline SimdDouble gmx_simdcall cbrt(SimdDouble x)
{
    const SimdDouble signBit(GMX_DOUBLE_NEGZERO);
    const SimdDouble minDouble(std::numeric_limits<double>::min());
    // Bias is 1024-1 = 1023, which is divisible by 3, so no need to change it more.
    // To avoid clang warnings about fragile integer division mixed with FP, we let
    // the divided value (1023/3=341) be the original constant.
    const std::int32_t offsetDiv3(341);
    const SimdDouble   c6(-0.145263899385486377);
    const SimdDouble   c5(0.784932344976639262);
    const SimdDouble   c4(-1.83469277483613086);
    const SimdDouble   c3(2.44693122563534430);
    const SimdDouble   c2(-2.11499494167371287);
    const SimdDouble   c1(1.50819193781584896);
    const SimdDouble   c0(0.354895765043919860);
    const SimdDouble   one(1.0);
    const SimdDouble   two(2.0);
    const SimdDouble   three(3.0);
    const SimdDouble   oneThird(1.0 / 3.0);
    const SimdDouble   cbrt2(1.2599210498948731648);
    const SimdDouble   sqrCbrt2(1.5874010519681994748);

    // See the single precision routines for documentation of the algorithm

    SimdDouble xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdDouble xAbs     = andNot(signBit, x);    // select everthing but the sign bit => abs(x)
    SimdDBool  xIsNonZero = (minDouble <= xAbs); // treat denormals as 0

    SimdDInt32 exponent;
    SimdDouble y        = frexp(xAbs, &exponent);
    SimdDouble z        = fma(y, c6, c5);
    z                   = fma(z, y, c4);
    z                   = fma(z, y, c3);
    z                   = fma(z, y, c2);
    z                   = fma(z, y, c1);
    z                   = fma(z, y, c0);
    SimdDouble w        = z * z * z;
    SimdDouble nom      = z * fma(two, y, w);
    SimdDouble invDenom = inv(fma(two, w, y));

    SimdDouble offsetExp = cvtI2R(exponent) + SimdDouble(static_cast<double>(3 * offsetDiv3) + 0.1);
    SimdDouble offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division
    SimdDInt32 expDiv3   = cvtR2I(offsetExpDiv3 - SimdDouble(static_cast<double>(offsetDiv3)));
    SimdDouble remainder = offsetExp - offsetExpDiv3 * three;
    SimdDouble factor    = blend(one, cbrt2, SimdDouble(0.5) < remainder);
    factor               = blend(factor, sqrCbrt2, SimdDouble(1.5) < remainder);
    SimdDouble fraction  = (nom * invDenom * factor) ^ xSignBit;
    SimdDouble result    = selectByMask(ldexp(fraction, expDiv3), xIsNonZero);
    return result;
}

/*! \brief Inverse cube root for SIMD doubles.
 *
 * \param x      Argument to calculate cube root of. Can be positive or
 *               negative, but the magnitude cannot be lower than
 *               the smallest normal number.
 * \return       Cube root of x. Undefined for values that don't
 *               fulfill the restriction of abs(x) > minDouble.
 */
static inline SimdDouble gmx_simdcall invcbrt(SimdDouble x)
{
    const SimdDouble signBit(GMX_DOUBLE_NEGZERO);
    // Bias is 1024-1 = 1023, which is divisible by 3, so no need to change it more.
    // To avoid clang warnings about fragile integer division mixed with FP, we let
    // the divided value (1023/3=341) be the original constant.
    const std::int32_t offsetDiv3(341);
    const SimdDouble   c6(-0.145263899385486377);
    const SimdDouble   c5(0.784932344976639262);
    const SimdDouble   c4(-1.83469277483613086);
    const SimdDouble   c3(2.44693122563534430);
    const SimdDouble   c2(-2.11499494167371287);
    const SimdDouble   c1(1.50819193781584896);
    const SimdDouble   c0(0.354895765043919860);
    const SimdDouble   one(1.0);
    const SimdDouble   two(2.0);
    const SimdDouble   three(3.0);
    const SimdDouble   oneThird(1.0 / 3.0);
    const SimdDouble   invCbrt2(1.0 / 1.2599210498948731648);
    const SimdDouble   invSqrCbrt2(1.0F / 1.5874010519681994748);

    // See the single precision routines for documentation of the algorithm

    SimdDouble xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdDouble xAbs     = andNot(signBit, x); // select everthing but the sign bit => abs(x)

    SimdDInt32 exponent;
    SimdDouble y         = frexp(xAbs, &exponent);
    SimdDouble z         = fma(y, c6, c5);
    z                    = fma(z, y, c4);
    z                    = fma(z, y, c3);
    z                    = fma(z, y, c2);
    z                    = fma(z, y, c1);
    z                    = fma(z, y, c0);
    SimdDouble w         = z * z * z;
    SimdDouble nom       = fma(two, w, y);
    SimdDouble invDenom  = inv(z * fma(two, y, w));
    SimdDouble offsetExp = cvtI2R(exponent) + SimdDouble(static_cast<double>(3 * offsetDiv3) + 0.1);
    SimdDouble offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division
    SimdDInt32 expDiv3   = cvtR2I(SimdDouble(static_cast<double>(offsetDiv3)) - offsetExpDiv3);
    SimdDouble remainder = offsetExpDiv3 * three - offsetExp;
    SimdDouble factor    = blend(one, invCbrt2, remainder < SimdDouble(-0.5));
    factor               = blend(factor, invSqrCbrt2, remainder < SimdDouble(-1.5));
    SimdDouble fraction  = (nom * invDenom * factor) ^ xSignBit;
    SimdDouble result    = ldexp(fraction, expDiv3);
    return result;
}

/*! \brief SIMD double log2(x). This is the base-2 logarithm.
 *
 * \param x Argument, should be >0.
 * \result The base-2 logarithm of x. Undefined if argument is invalid.
 */
static inline SimdDouble gmx_simdcall log2(SimdDouble x)
{
#        if GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE
    // Just rescale if native log2() is not present, but log is.
    return log(x) * SimdDouble(std::log2(std::exp(1.0)));
#        else
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble invsqrt2(1.0 / std::sqrt(2.0));
    const SimdDouble CL15(0.2138031565795550370534528);
    const SimdDouble CL13(0.2208884091496370882801159);
    const SimdDouble CL11(0.2623358279761824340958754);
    const SimdDouble CL9(0.3205984930182496084327681);
    const SimdDouble CL7(0.4121985864521960363227038);
    const SimdDouble CL5(0.5770780163410746954610886);
    const SimdDouble CL3(0.9617966939260027547931031);
    const SimdDouble CL1(2.885390081777926774009302);
    SimdDouble       fExp, x2, p;
    SimdDBool        m;
    SimdDInt32       iExp;

    x    = frexp(x, &iExp);
    fExp = cvtI2R(iExp);

    m = x < invsqrt2;
    // Adjust to non-IEEE format for x<1/sqrt(2): exponent -= 1, mantissa *= 2.0
    fExp = fExp - selectByMask(one, m);
    x    = x * blend(one, two, m);

    x  = (x - one) * inv(x + one);
    x2 = x * x;

    p = fma(CL15, x2, CL13);
    p = fma(p, x2, CL11);
    p = fma(p, x2, CL9);
    p = fma(p, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, fExp);

    return p;
#        endif
}

#        if !GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE
/*! \brief SIMD double log(x). This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
static inline SimdDouble gmx_simdcall log(SimdDouble x)
{
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble invsqrt2(1.0 / std::sqrt(2.0));
    const SimdDouble corr(0.693147180559945286226764);
    const SimdDouble CL15(0.148197055177935105296783);
    const SimdDouble CL13(0.153108178020442575739679);
    const SimdDouble CL11(0.181837339521549679055568);
    const SimdDouble CL9(0.22222194152736701733275);
    const SimdDouble CL7(0.285714288030134544449368);
    const SimdDouble CL5(0.399999999989941956712869);
    const SimdDouble CL3(0.666666666666685503450651);
    const SimdDouble CL1(2.0);
    SimdDouble       fExp, x2, p;
    SimdDBool        m;
    SimdDInt32       iExp;

    x    = frexp(x, &iExp);
    fExp = cvtI2R(iExp);

    m = x < invsqrt2;
    // Adjust to non-IEEE format for x<1/sqrt(2): exponent -= 1, mantissa *= 2.0
    fExp = fExp - selectByMask(one, m);
    x    = x * blend(one, two, m);

    x  = (x - one) * inv(x + one);
    x2 = x * x;

    p = fma(CL15, x2, CL13);
    p = fma(p, x2, CL11);
    p = fma(p, x2, CL9);
    p = fma(p, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, corr * fExp);

    return p;
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_EXP2_DOUBLE
/*! \brief SIMD double 2^x.
 *
 * \copydetails exp2(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall exp2(SimdDouble x)
{
    const SimdDouble CE11(4.435280790452730022081181e-10);
    const SimdDouble CE10(7.074105630863314448024247e-09);
    const SimdDouble CE9(1.017819803432096698472621e-07);
    const SimdDouble CE8(1.321543308956718799557863e-06);
    const SimdDouble CE7(0.00001525273348995851746990884);
    const SimdDouble CE6(0.0001540353046251466849082632);
    const SimdDouble CE5(0.001333355814678995257307880);
    const SimdDouble CE4(0.009618129107588335039176502);
    const SimdDouble CE3(0.05550410866481992147457793);
    const SimdDouble CE2(0.2402265069591015620470894);
    const SimdDouble CE1(0.6931471805599453304615075);
    const SimdDouble one(1.0);

    SimdDouble intpart;
    SimdDouble fexppart;
    SimdDouble p;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -1023, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer.
    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdDouble(std::numeric_limits<std::int32_t>::lowest()));
    }

    fexppart = ldexp<opt>(one, cvtR2I(x));
    intpart  = round(x);
    x        = x - intpart;

    p = fma(CE11, x, CE10);
    p = fma(p, x, CE9);
    p = fma(p, x, CE8);
    p = fma(p, x, CE7);
    p = fma(p, x, CE6);
    p = fma(p, x, CE5);
    p = fma(p, x, CE4);
    p = fma(p, x, CE3);
    p = fma(p, x, CE2);
    p = fma(p, x, CE1);
    p = fma(p, x, one);
    x = p * fexppart;
    return x;
}
#        endif

#        if !GMX_SIMD_HAVE_NATIVE_EXP_DOUBLE
/*! \brief SIMD double exp(x).
 *
 * \copydetails exp(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall exp(SimdDouble x)
{
    const SimdDouble argscale(1.44269504088896340735992468100);
    const SimdDouble invargscale0(-0.69314718055966295651160180568695068359375);
    const SimdDouble invargscale1(-2.8235290563031577122588448175013436025525412068e-13);
    const SimdDouble CE12(2.078375306791423699350304e-09);
    const SimdDouble CE11(2.518173854179933105218635e-08);
    const SimdDouble CE10(2.755842049600488770111608e-07);
    const SimdDouble CE9(2.755691815216689746619849e-06);
    const SimdDouble CE8(2.480158383706245033920920e-05);
    const SimdDouble CE7(0.0001984127043518048611841321);
    const SimdDouble CE6(0.001388888889360258341755930);
    const SimdDouble CE5(0.008333333332907368102819109);
    const SimdDouble CE4(0.04166666666663836745814631);
    const SimdDouble CE3(0.1666666666666796929434570);
    const SimdDouble CE2(0.5);
    const SimdDouble one(1.0);
    SimdDouble       fexppart;
    SimdDouble       intpart;
    SimdDouble       y, p;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -1023, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer
    // after being multiplied by argscale.

    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdDouble(std::numeric_limits<std::int32_t>::lowest()) / argscale);
    }

    y = x * argscale;

    fexppart = ldexp<opt>(one, cvtR2I(y));
    intpart  = round(y);

    // Extended precision arithmetics
    x = fma(invargscale0, intpart, x);
    x = fma(invargscale1, intpart, x);

    p = fma(CE12, x, CE11);
    p = fma(p, x, CE10);
    p = fma(p, x, CE9);
    p = fma(p, x, CE8);
    p = fma(p, x, CE7);
    p = fma(p, x, CE6);
    p = fma(p, x, CE5);
    p = fma(p, x, CE4);
    p = fma(p, x, CE3);
    p = fma(p, x, CE2);
    p = fma(p, x * x, x);
#            if GMX_SIMD_HAVE_FMA
    x = fma(p, fexppart, fexppart);
#            else
    x = (p + one) * fexppart;
#            endif

    return x;
}
#        endif

/*! \brief SIMD double pow(x,y)
 *
 * This returns x^y for SIMD values.
 *
 * \tparam opt If this is changed from the default (safe) into the unsafe
 *             option, there are no guarantees about correct results for x==0.
 *
 * \param x Base.
 *
 * \param y exponent.
 *
 * \result x^y. Overflowing arguments are likely to either return 0 or inf,
 *         depending on the underlying implementation. If unsafe optimizations
 *         are enabled, this is also true for x==0.
 *
 * \warning You cannot rely on this implementation returning inf for arguments
 *          that cause overflow. If you have some very large
 *          values and need to rely on getting a valid numerical output,
 *          take the minimum of your variable and the largest valid argument
 *          before calling this routine.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall pow(SimdDouble x, SimdDouble y)
{
    SimdDouble xcorr;

    if (opt == MathOptimization::Safe)
    {
        xcorr = max(x, SimdDouble(std::numeric_limits<double>::min()));
    }
    else
    {
        xcorr = x;
    }

    SimdDouble result = exp2<opt>(y * log2(xcorr));

    if (opt == MathOptimization::Safe)
    {
        // if x==0 and y>0 we explicitly set the result to 0.0
        // For any x with y==0, the result will already be 1.0 since we multiply by y (0.0) and call exp().
        result = blend(result, setZero(), x == setZero() && setZero() < y);
    }

    return result;
}


/*! \brief SIMD double erf(x).
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to full precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdDouble gmx_simdcall erf(SimdDouble x)
{
    // Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75]
    const SimdDouble CAP4(-0.431780540597889301512e-4);
    const SimdDouble CAP3(-0.00578562306260059236059);
    const SimdDouble CAP2(-0.028593586920219752446);
    const SimdDouble CAP1(-0.315924962948621698209);
    const SimdDouble CAP0(0.14952975608477029151);

    const SimdDouble CAQ5(-0.374089300177174709737e-5);
    const SimdDouble CAQ4(0.00015126584532155383535);
    const SimdDouble CAQ3(0.00536692680669480725423);
    const SimdDouble CAQ2(0.0668686825594046122636);
    const SimdDouble CAQ1(0.402604990869284362773);
    // CAQ0 == 1.0
    const SimdDouble CAoffset(0.9788494110107421875);

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5]
    const SimdDouble CBP6(2.49650423685462752497647637088e-10);
    const SimdDouble CBP5(0.00119770193298159629350136085658);
    const SimdDouble CBP4(0.0164944422378370965881008942733);
    const SimdDouble CBP3(0.0984581468691775932063932439252);
    const SimdDouble CBP2(0.317364595806937763843589437418);
    const SimdDouble CBP1(0.554167062641455850932670067075);
    const SimdDouble CBP0(0.427583576155807163756925301060);
    const SimdDouble CBQ7(0.00212288829699830145976198384930);
    const SimdDouble CBQ6(0.0334810979522685300554606393425);
    const SimdDouble CBQ5(0.2361713785181450957579508850717);
    const SimdDouble CBQ4(0.955364736493055670530981883072);
    const SimdDouble CBQ3(2.36815675631420037315349279199);
    const SimdDouble CBQ2(3.55261649184083035537184223542);
    const SimdDouble CBQ1(2.93501136050160872574376997993);
    // CBQ0 == 1.0

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf]
    const SimdDouble CCP6(-2.8175401114513378771);
    const SimdDouble CCP5(-3.22729451764143718517);
    const SimdDouble CCP4(-2.5518551727311523996);
    const SimdDouble CCP3(-0.687717681153649930619);
    const SimdDouble CCP2(-0.212652252872804219852);
    const SimdDouble CCP1(0.0175389834052493308818);
    const SimdDouble CCP0(0.00628057170626964891937);

    const SimdDouble CCQ6(5.48409182238641741584);
    const SimdDouble CCQ5(13.5064170191802889145);
    const SimdDouble CCQ4(22.9367376522880577224);
    const SimdDouble CCQ3(15.930646027911794143);
    const SimdDouble CCQ2(11.0567237927800161565);
    const SimdDouble CCQ1(2.79257750980575282228);
    // CCQ0 == 1.0
    const SimdDouble CCoffset(0.5579090118408203125);

    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble minFloat(std::numeric_limits<float>::min());

    SimdDouble xabs, x2, x4, t, t2, w, w2;
    SimdDouble PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    SimdDouble PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    SimdDouble PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    SimdDouble res_erf, res_erfcB, res_erfcC, res_erfc, res;
    SimdDouble expmx2;
    SimdDBool  mask, mask_erf, notmask_erf;

    // Calculate erf()
    xabs        = abs(x);
    mask_erf    = (xabs < one);
    notmask_erf = (one <= xabs);
    x2          = x * x;
    x4          = x2 * x2;

    PolyAP0 = fma(CAP4, x4, CAP2);
    PolyAP1 = fma(CAP3, x4, CAP1);
    PolyAP0 = fma(PolyAP0, x4, CAP0);
    PolyAP0 = fma(PolyAP1, x2, PolyAP0);

    PolyAQ1 = fma(CAQ5, x4, CAQ3);
    PolyAQ0 = fma(CAQ4, x4, CAQ2);
    PolyAQ1 = fma(PolyAQ1, x4, CAQ1);
    PolyAQ0 = fma(PolyAQ0, x4, one);
    PolyAQ0 = fma(PolyAQ1, x2, PolyAQ0);

    res_erf = PolyAP0 * maskzInv(PolyAQ0, mask_erf && (minFloat <= abs(PolyAQ0)));
    res_erf = CAoffset + res_erf;
    res_erf = x * res_erf;

    // Calculate erfc() in range [1,4.5]
    t  = xabs - one;
    t2 = t * t;

    PolyBP0 = fma(CBP6, t2, CBP4);
    PolyBP1 = fma(CBP5, t2, CBP3);
    PolyBP0 = fma(PolyBP0, t2, CBP2);
    PolyBP1 = fma(PolyBP1, t2, CBP1);
    PolyBP0 = fma(PolyBP0, t2, CBP0);
    PolyBP0 = fma(PolyBP1, t, PolyBP0);

    PolyBQ1 = fma(CBQ7, t2, CBQ5);
    PolyBQ0 = fma(CBQ6, t2, CBQ4);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ3);
    PolyBQ0 = fma(PolyBQ0, t2, CBQ2);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ1);
    PolyBQ0 = fma(PolyBQ0, t2, one);
    PolyBQ0 = fma(PolyBQ1, t, PolyBQ0);

    // The denominator polynomial can be zero outside the range
    res_erfcB = PolyBP0 * maskzInv(PolyBQ0, notmask_erf && (minFloat <= abs(PolyBQ0)));

    res_erfcB = res_erfcB * xabs;

    // Calculate erfc() in range [4.5,inf]
    w  = maskzInv(xabs, notmask_erf && (minFloat <= xabs));
    w2 = w * w;

    PolyCP0 = fma(CCP6, w2, CCP4);
    PolyCP1 = fma(CCP5, w2, CCP3);
    PolyCP0 = fma(PolyCP0, w2, CCP2);
    PolyCP1 = fma(PolyCP1, w2, CCP1);
    PolyCP0 = fma(PolyCP0, w2, CCP0);
    PolyCP0 = fma(PolyCP1, w, PolyCP0);

    PolyCQ0 = fma(CCQ6, w2, CCQ4);
    PolyCQ1 = fma(CCQ5, w2, CCQ3);
    PolyCQ0 = fma(PolyCQ0, w2, CCQ2);
    PolyCQ1 = fma(PolyCQ1, w2, CCQ1);
    PolyCQ0 = fma(PolyCQ0, w2, one);
    PolyCQ0 = fma(PolyCQ1, w, PolyCQ0);

    expmx2 = exp(-x2);

    // The denominator polynomial can be zero outside the range
    res_erfcC = PolyCP0 * maskzInv(PolyCQ0, notmask_erf && (minFloat <= abs(PolyCQ0)));
    res_erfcC = res_erfcC + CCoffset;
    res_erfcC = res_erfcC * w;

    mask     = (SimdDouble(4.5) < xabs);
    res_erfc = blend(res_erfcB, res_erfcC, mask);

    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    mask     = (x < setZero());
    res_erfc = blend(res_erfc, two - res_erfc, mask);

    // Select erf() or erfc()
    res = blend(one - res_erfc, res_erf, mask_erf);

    return res;
}

/*! \brief SIMD double erfc(x).
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves full precision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-200 for double)
 * since that is not relevant for MD.
 */
static inline SimdDouble gmx_simdcall erfc(SimdDouble x)
{
    // Coefficients for minimax approximation of erf(x)=x*(CAoffset + P(x^2)/Q(x^2)) in range [-0.75,0.75]
    const SimdDouble CAP4(-0.431780540597889301512e-4);
    const SimdDouble CAP3(-0.00578562306260059236059);
    const SimdDouble CAP2(-0.028593586920219752446);
    const SimdDouble CAP1(-0.315924962948621698209);
    const SimdDouble CAP0(0.14952975608477029151);

    const SimdDouble CAQ5(-0.374089300177174709737e-5);
    const SimdDouble CAQ4(0.00015126584532155383535);
    const SimdDouble CAQ3(0.00536692680669480725423);
    const SimdDouble CAQ2(0.0668686825594046122636);
    const SimdDouble CAQ1(0.402604990869284362773);
    // CAQ0 == 1.0
    const SimdDouble CAoffset(0.9788494110107421875);

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)*x*(P(x-1)/Q(x-1)) in range [1.0,4.5]
    const SimdDouble CBP6(2.49650423685462752497647637088e-10);
    const SimdDouble CBP5(0.00119770193298159629350136085658);
    const SimdDouble CBP4(0.0164944422378370965881008942733);
    const SimdDouble CBP3(0.0984581468691775932063932439252);
    const SimdDouble CBP2(0.317364595806937763843589437418);
    const SimdDouble CBP1(0.554167062641455850932670067075);
    const SimdDouble CBP0(0.427583576155807163756925301060);
    const SimdDouble CBQ7(0.00212288829699830145976198384930);
    const SimdDouble CBQ6(0.0334810979522685300554606393425);
    const SimdDouble CBQ5(0.2361713785181450957579508850717);
    const SimdDouble CBQ4(0.955364736493055670530981883072);
    const SimdDouble CBQ3(2.36815675631420037315349279199);
    const SimdDouble CBQ2(3.55261649184083035537184223542);
    const SimdDouble CBQ1(2.93501136050160872574376997993);
    // CBQ0 == 1.0

    // Coefficients for minimax approximation of erfc(x)=exp(-x^2)/x*(P(1/x)/Q(1/x)) in range [4.5,inf]
    const SimdDouble CCP6(-2.8175401114513378771);
    const SimdDouble CCP5(-3.22729451764143718517);
    const SimdDouble CCP4(-2.5518551727311523996);
    const SimdDouble CCP3(-0.687717681153649930619);
    const SimdDouble CCP2(-0.212652252872804219852);
    const SimdDouble CCP1(0.0175389834052493308818);
    const SimdDouble CCP0(0.00628057170626964891937);

    const SimdDouble CCQ6(5.48409182238641741584);
    const SimdDouble CCQ5(13.5064170191802889145);
    const SimdDouble CCQ4(22.9367376522880577224);
    const SimdDouble CCQ3(15.930646027911794143);
    const SimdDouble CCQ2(11.0567237927800161565);
    const SimdDouble CCQ1(2.79257750980575282228);
    // CCQ0 == 1.0
    const SimdDouble CCoffset(0.5579090118408203125);

    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble minFloat(std::numeric_limits<float>::min());

    SimdDouble xabs, x2, x4, t, t2, w, w2;
    SimdDouble PolyAP0, PolyAP1, PolyAQ0, PolyAQ1;
    SimdDouble PolyBP0, PolyBP1, PolyBQ0, PolyBQ1;
    SimdDouble PolyCP0, PolyCP1, PolyCQ0, PolyCQ1;
    SimdDouble res_erf, res_erfcB, res_erfcC, res_erfc, res;
    SimdDouble expmx2;
    SimdDBool  mask, mask_erf, notmask_erf;

    // Calculate erf()
    xabs        = abs(x);
    mask_erf    = (xabs < one);
    notmask_erf = (one <= xabs);
    x2          = x * x;
    x4          = x2 * x2;

    PolyAP0 = fma(CAP4, x4, CAP2);
    PolyAP1 = fma(CAP3, x4, CAP1);
    PolyAP0 = fma(PolyAP0, x4, CAP0);
    PolyAP0 = fma(PolyAP1, x2, PolyAP0);
    PolyAQ1 = fma(CAQ5, x4, CAQ3);
    PolyAQ0 = fma(CAQ4, x4, CAQ2);
    PolyAQ1 = fma(PolyAQ1, x4, CAQ1);
    PolyAQ0 = fma(PolyAQ0, x4, one);
    PolyAQ0 = fma(PolyAQ1, x2, PolyAQ0);

    res_erf = PolyAP0 * maskzInv(PolyAQ0, mask_erf && (minFloat <= abs(PolyAQ0)));
    res_erf = CAoffset + res_erf;
    res_erf = x * res_erf;

    // Calculate erfc() in range [1,4.5]
    t  = xabs - one;
    t2 = t * t;

    PolyBP0 = fma(CBP6, t2, CBP4);
    PolyBP1 = fma(CBP5, t2, CBP3);
    PolyBP0 = fma(PolyBP0, t2, CBP2);
    PolyBP1 = fma(PolyBP1, t2, CBP1);
    PolyBP0 = fma(PolyBP0, t2, CBP0);
    PolyBP0 = fma(PolyBP1, t, PolyBP0);

    PolyBQ1 = fma(CBQ7, t2, CBQ5);
    PolyBQ0 = fma(CBQ6, t2, CBQ4);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ3);
    PolyBQ0 = fma(PolyBQ0, t2, CBQ2);
    PolyBQ1 = fma(PolyBQ1, t2, CBQ1);
    PolyBQ0 = fma(PolyBQ0, t2, one);
    PolyBQ0 = fma(PolyBQ1, t, PolyBQ0);

    // The denominator polynomial can be zero outside the range
    res_erfcB = PolyBP0 * maskzInv(PolyBQ0, notmask_erf && (minFloat <= abs(PolyBQ0)));

    res_erfcB = res_erfcB * xabs;

    // Calculate erfc() in range [4.5,inf]
    // Note that 1/x can only handle single precision!
    w  = maskzInv(xabs, minFloat <= xabs);
    w2 = w * w;

    PolyCP0 = fma(CCP6, w2, CCP4);
    PolyCP1 = fma(CCP5, w2, CCP3);
    PolyCP0 = fma(PolyCP0, w2, CCP2);
    PolyCP1 = fma(PolyCP1, w2, CCP1);
    PolyCP0 = fma(PolyCP0, w2, CCP0);
    PolyCP0 = fma(PolyCP1, w, PolyCP0);

    PolyCQ0 = fma(CCQ6, w2, CCQ4);
    PolyCQ1 = fma(CCQ5, w2, CCQ3);
    PolyCQ0 = fma(PolyCQ0, w2, CCQ2);
    PolyCQ1 = fma(PolyCQ1, w2, CCQ1);
    PolyCQ0 = fma(PolyCQ0, w2, one);
    PolyCQ0 = fma(PolyCQ1, w, PolyCQ0);

    expmx2 = exp(-x2);

    // The denominator polynomial can be zero outside the range
    res_erfcC = PolyCP0 * maskzInv(PolyCQ0, notmask_erf && (minFloat <= abs(PolyCQ0)));
    res_erfcC = res_erfcC + CCoffset;
    res_erfcC = res_erfcC * w;

    mask     = (SimdDouble(4.5) < xabs);
    res_erfc = blend(res_erfcB, res_erfcC, mask);

    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    mask     = (x < setZero());
    res_erfc = blend(res_erfc, two - res_erfc, mask);

    // Select erf() or erfc()
    res = blend(res_erfc, one - res_erf, mask_erf);

    return res;
}

/*! \brief SIMD double sin \& cos.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * This version achieves close to machine precision, but for very large
 * magnitudes of the argument we inherently begin to lose accuracy due to the
 * argument reduction, despite using extended precision arithmetics internally.
 */
static inline void gmx_simdcall sincos(SimdDouble x, SimdDouble* sinval, SimdDouble* cosval)
{
    // Constants to subtract Pi/4*x from y while minimizing precision loss
    const SimdDouble argred0(-2 * 0.78539816290140151978);
    const SimdDouble argred1(-2 * 4.9604678871439933374e-10);
    const SimdDouble argred2(-2 * 1.1258708853173288931e-18);
    const SimdDouble argred3(-2 * 1.7607799325916000908e-27);
    const SimdDouble two_over_pi(2.0 / M_PI);
    const SimdDouble const_sin5(1.58938307283228937328511e-10);
    const SimdDouble const_sin4(-2.50506943502539773349318e-08);
    const SimdDouble const_sin3(2.75573131776846360512547e-06);
    const SimdDouble const_sin2(-0.000198412698278911770864914);
    const SimdDouble const_sin1(0.0083333333333191845961746);
    const SimdDouble const_sin0(-0.166666666666666130709393);

    const SimdDouble const_cos7(-1.13615350239097429531523e-11);
    const SimdDouble const_cos6(2.08757471207040055479366e-09);
    const SimdDouble const_cos5(-2.75573144028847567498567e-07);
    const SimdDouble const_cos4(2.48015872890001867311915e-05);
    const SimdDouble const_cos3(-0.00138888888888714019282329);
    const SimdDouble const_cos2(0.0416666666666665519592062);
    const SimdDouble half(0.5);
    const SimdDouble one(1.0);
    SimdDouble       ssign, csign;
    SimdDouble       x2, y, z, psin, pcos, sss, ccc;
    SimdDBool        mask;
#        if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdDInt32 ione(1);
    const SimdDInt32 itwo(2);
    SimdDInt32       iy;

    z  = x * two_over_pi;
    iy = cvtR2I(z);
    y  = round(z);

    mask  = cvtIB2B((iy & ione) == setZero());
    ssign = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), cvtIB2B((iy & itwo) == itwo));
    csign = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), cvtIB2B(((iy + ione) & itwo) == itwo));
#        else
    const SimdDouble quarter(0.25);
    const SimdDouble minusquarter(-0.25);
    SimdDouble       q;
    SimdDBool        m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign = x;
    x     = abs(x);
    // It is critical that half-way cases are rounded down
    z = fma(x, two_over_pi, half);
    y = trunc(z);
    q = z * quarter;
    q = q - trunc(q);
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q = q - half;
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#            if GMX_SIMD_HAVE_LOGICAL
    ssign = ssign & SimdDouble(GMX_DOUBLE_NEGZERO);
    csign = andNot(q, SimdDouble(GMX_DOUBLE_NEGZERO));
    ssign = ssign ^ csign;
#            else
    ssign = copysign(SimdDouble(1.0), ssign);
    csign = copysign(SimdDouble(1.0), q);
    csign = -csign;
    ssign = ssign * csign; // swap ssign if csign was set.
#            endif
    // Check if y had value 1 or 3 (remember we subtracted 0.5 from q)
    m1   = (q < minusquarter);
    m2   = (setZero() <= q);
    m3   = (q < quarter);
    m2   = m2 && m3;
    mask = m1 || m2;
    // where mask is FALSE, swap sign.
    csign   = csign * blend(SimdDouble(-1.0), one, mask);
#        endif
    x  = fma(y, argred0, x);
    x  = fma(y, argred1, x);
    x  = fma(y, argred2, x);
    x  = fma(y, argred3, x);
    x2 = x * x;

    psin = fma(const_sin5, x2, const_sin4);
    psin = fma(psin, x2, const_sin3);
    psin = fma(psin, x2, const_sin2);
    psin = fma(psin, x2, const_sin1);
    psin = fma(psin, x2, const_sin0);
    psin = fma(psin, x2 * x, x);

    pcos = fma(const_cos7, x2, const_cos6);
    pcos = fma(pcos, x2, const_cos5);
    pcos = fma(pcos, x2, const_cos4);
    pcos = fma(pcos, x2, const_cos3);
    pcos = fma(pcos, x2, const_cos2);
    pcos = fms(pcos, x2, half);
    pcos = fma(pcos, x2, one);

    sss = blend(pcos, psin, mask);
    ccc = blend(psin, pcos, mask);
    // See comment for GMX_SIMD_HAVE_LOGICAL section above.
#        if GMX_SIMD_HAVE_LOGICAL
    *sinval = sss ^ ssign;
    *cosval = ccc ^ csign;
#        else
    *sinval = sss * ssign;
    *cosval = ccc * csign;
#        endif
}

/*! \brief SIMD double sin(x).
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall sin(SimdDouble x)
{
    SimdDouble s, c;
    sincos(x, &s, &c);
    return s;
}

/*! \brief SIMD double cos(x).
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall cos(SimdDouble x)
{
    SimdDouble s, c;
    sincos(x, &s, &c);
    return c;
}

/*! \brief SIMD double tan(x).
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdDouble gmx_simdcall tan(SimdDouble x)
{
    const SimdDouble argred0(-2 * 0.78539816290140151978);
    const SimdDouble argred1(-2 * 4.9604678871439933374e-10);
    const SimdDouble argred2(-2 * 1.1258708853173288931e-18);
    const SimdDouble argred3(-2 * 1.7607799325916000908e-27);
    const SimdDouble two_over_pi(2.0 / M_PI);
    const SimdDouble CT15(1.01419718511083373224408e-05);
    const SimdDouble CT14(-2.59519791585924697698614e-05);
    const SimdDouble CT13(5.23388081915899855325186e-05);
    const SimdDouble CT12(-3.05033014433946488225616e-05);
    const SimdDouble CT11(7.14707504084242744267497e-05);
    const SimdDouble CT10(8.09674518280159187045078e-05);
    const SimdDouble CT9(0.000244884931879331847054404);
    const SimdDouble CT8(0.000588505168743587154904506);
    const SimdDouble CT7(0.00145612788922812427978848);
    const SimdDouble CT6(0.00359208743836906619142924);
    const SimdDouble CT5(0.00886323944362401618113356);
    const SimdDouble CT4(0.0218694882853846389592078);
    const SimdDouble CT3(0.0539682539781298417636002);
    const SimdDouble CT2(0.133333333333125941821962);
    const SimdDouble CT1(0.333333333333334980164153);
    const SimdDouble minFloat(std::numeric_limits<float>::min());

    SimdDouble x2, p, y, z;
    SimdDBool  m;

#        if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdDInt32 iy;
    SimdDInt32 ione(1);

    z  = x * two_over_pi;
    iy = cvtR2I(z);
    y  = round(z);
    m  = cvtIB2B((iy & ione) == ione);

    x = fma(y, argred0, x);
    x = fma(y, argred1, x);
    x = fma(y, argred2, x);
    x = fma(y, argred3, x);
    x = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), m) ^ x;
#        else
    const SimdDouble quarter(0.25);
    const SimdDouble half(0.5);
    const SimdDouble threequarter(0.75);
    const SimdDouble minFloat(std::numeric_limits<float>::min());
    SimdDouble       w, q;
    SimdDBool        m1, m2, m3;

    w  = abs(x);
    z  = fma(w, two_over_pi, half);
    y  = trunc(z);
    q  = z * quarter;
    q  = q - trunc(q);
    m1 = (quarter <= q);
    m2 = (q < half);
    m3 = (threequarter <= q);
    m1 = m1 && m2;
    m  = m1 || m3;
    w  = fma(y, argred0, w);
    w  = fma(y, argred1, w);
    w  = fma(y, argred2, w);
    w  = fma(y, argred3, w);

    w     = blend(w, -w, m);
    x     = w * copysign(SimdDouble(1.0), x);
#        endif
    x2 = x * x;
    p  = fma(CT15, x2, CT14);
    p  = fma(p, x2, CT13);
    p  = fma(p, x2, CT12);
    p  = fma(p, x2, CT11);
    p  = fma(p, x2, CT10);
    p  = fma(p, x2, CT9);
    p  = fma(p, x2, CT8);
    p  = fma(p, x2, CT7);
    p  = fma(p, x2, CT6);
    p  = fma(p, x2, CT5);
    p  = fma(p, x2, CT4);
    p  = fma(p, x2, CT3);
    p  = fma(p, x2, CT2);
    p  = fma(p, x2, CT1);
    p  = fma(x2, p * x, x);

    p = blend(p, maskzInv(p, m && (minFloat < abs(p))), m);
    return p;
}

/*! \brief SIMD double asin(x).
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdDouble gmx_simdcall asin(SimdDouble x)
{
    // Same algorithm as cephes library
    const SimdDouble limit1(0.625);
    const SimdDouble limit2(1e-8);
    const SimdDouble one(1.0);
    const SimdDouble quarterpi(M_PI / 4.0);
    const SimdDouble morebits(6.123233995736765886130e-17);

    const SimdDouble P5(4.253011369004428248960e-3);
    const SimdDouble P4(-6.019598008014123785661e-1);
    const SimdDouble P3(5.444622390564711410273e0);
    const SimdDouble P2(-1.626247967210700244449e1);
    const SimdDouble P1(1.956261983317594739197e1);
    const SimdDouble P0(-8.198089802484824371615e0);

    const SimdDouble Q4(-1.474091372988853791896e1);
    const SimdDouble Q3(7.049610280856842141659e1);
    const SimdDouble Q2(-1.471791292232726029859e2);
    const SimdDouble Q1(1.395105614657485689735e2);
    const SimdDouble Q0(-4.918853881490881290097e1);

    const SimdDouble R4(2.967721961301243206100e-3);
    const SimdDouble R3(-5.634242780008963776856e-1);
    const SimdDouble R2(6.968710824104713396794e0);
    const SimdDouble R1(-2.556901049652824852289e1);
    const SimdDouble R0(2.853665548261061424989e1);

    const SimdDouble S3(-2.194779531642920639778e1);
    const SimdDouble S2(1.470656354026814941758e2);
    const SimdDouble S1(-3.838770957603691357202e2);
    const SimdDouble S0(3.424398657913078477438e2);

    SimdDouble xabs;
    SimdDouble zz, ww, z, q, w, zz2, ww2;
    SimdDouble PA, PB;
    SimdDouble QA, QB;
    SimdDouble RA, RB;
    SimdDouble SA, SB;
    SimdDouble nom, denom;
    SimdDBool  mask, mask2;

    xabs = abs(x);

    mask = (limit1 < xabs);

    zz  = one - xabs;
    ww  = xabs * xabs;
    zz2 = zz * zz;
    ww2 = ww * ww;

    // R
    RA = fma(R4, zz2, R2);
    RB = fma(R3, zz2, R1);
    RA = fma(RA, zz2, R0);
    RA = fma(RB, zz, RA);

    // S, SA = zz2
    SB = fma(S3, zz2, S1);
    SA = zz2 + S2;
    SA = fma(SA, zz2, S0);
    SA = fma(SB, zz, SA);

    // P
    PA = fma(P5, ww2, P3);
    PB = fma(P4, ww2, P2);
    PA = fma(PA, ww2, P1);
    PB = fma(PB, ww2, P0);
    PA = fma(PA, ww, PB);

    // Q, QA = ww2
    QB = fma(Q4, ww2, Q2);
    QA = ww2 + Q3;
    QA = fma(QA, ww2, Q1);
    QB = fma(QB, ww2, Q0);
    QA = fma(QA, ww, QB);

    RA = RA * zz;
    PA = PA * ww;

    nom   = blend(PA, RA, mask);
    denom = blend(QA, SA, mask);

    mask2 = (limit2 < xabs);
    q     = nom * maskzInv(denom, mask2);

    zz = zz + zz;
    zz = sqrt(zz);
    z  = quarterpi - zz;
    zz = fms(zz, q, morebits);
    z  = z - zz;
    z  = z + quarterpi;

    w = xabs * q;
    w = w + xabs;

    z = blend(w, z, mask);

    z = blend(xabs, z, mask2);

    z = copysign(z, x);

    return z;
}

/*! \brief SIMD double acos(x).
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdDouble gmx_simdcall acos(SimdDouble x)
{
    const SimdDouble one(1.0);
    const SimdDouble half(0.5);
    const SimdDouble quarterpi0(7.85398163397448309616e-1);
    const SimdDouble quarterpi1(6.123233995736765886130e-17);

    SimdDBool  mask1;
    SimdDouble z, z1, z2;

    mask1 = (half < x);
    z1    = half * (one - x);
    z1    = sqrt(z1);
    z     = blend(x, z1, mask1);

    z = asin(z);

    z1 = z + z;

    z2 = quarterpi0 - z;
    z2 = z2 + quarterpi1;
    z2 = z2 + quarterpi0;

    z = blend(z2, z1, mask1);

    return z;
}

/*! \brief SIMD double asin(x).
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdDouble gmx_simdcall atan(SimdDouble x)
{
    // Same algorithm as cephes library
    const SimdDouble limit1(0.66);
    const SimdDouble limit2(2.41421356237309504880);
    const SimdDouble quarterpi(M_PI / 4.0);
    const SimdDouble halfpi(M_PI / 2.0);
    const SimdDouble mone(-1.0);
    const SimdDouble morebits1(0.5 * 6.123233995736765886130E-17);
    const SimdDouble morebits2(6.123233995736765886130E-17);

    const SimdDouble P4(-8.750608600031904122785E-1);
    const SimdDouble P3(-1.615753718733365076637E1);
    const SimdDouble P2(-7.500855792314704667340E1);
    const SimdDouble P1(-1.228866684490136173410E2);
    const SimdDouble P0(-6.485021904942025371773E1);

    const SimdDouble Q4(2.485846490142306297962E1);
    const SimdDouble Q3(1.650270098316988542046E2);
    const SimdDouble Q2(4.328810604912902668951E2);
    const SimdDouble Q1(4.853903996359136964868E2);
    const SimdDouble Q0(1.945506571482613964425E2);

    SimdDouble y, xabs, t1, t2;
    SimdDouble z, z2;
    SimdDouble P_A, P_B, Q_A, Q_B;
    SimdDBool  mask1, mask2;

    xabs = abs(x);

    mask1 = (limit1 < xabs);
    mask2 = (limit2 < xabs);

    t1 = (xabs + mone) * maskzInv(xabs - mone, mask1);
    t2 = mone * maskzInv(xabs, mask2);

    y    = selectByMask(quarterpi, mask1);
    y    = blend(y, halfpi, mask2);
    xabs = blend(xabs, t1, mask1);
    xabs = blend(xabs, t2, mask2);

    z  = xabs * xabs;
    z2 = z * z;

    P_A = fma(P4, z2, P2);
    P_B = fma(P3, z2, P1);
    P_A = fma(P_A, z2, P0);
    P_A = fma(P_B, z, P_A);

    // Q_A = z2
    Q_B = fma(Q4, z2, Q2);
    Q_A = z2 + Q3;
    Q_A = fma(Q_A, z2, Q1);
    Q_B = fma(Q_B, z2, Q0);
    Q_A = fma(Q_A, z, Q_B);

    z = z * P_A;
    z = z * inv(Q_A);
    z = fma(z, xabs, xabs);

    t1 = selectByMask(morebits1, mask1);
    t1 = blend(t1, morebits2, mask2);

    z = z + t1;
    y = y + z;

    y = copysign(y, x);

    return y;
}

/*! \brief SIMD double atan2(y,x).
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static inline SimdDouble gmx_simdcall atan2(SimdDouble y, SimdDouble x)
{
    const SimdDouble pi(M_PI);
    const SimdDouble halfpi(M_PI / 2.0);
    const SimdDouble minFloat(std::numeric_limits<float>::min());
    SimdDouble       xinv, p, aoffset;
    SimdDBool        mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = x != setZero();
    mask_ynz  = y != setZero();
    mask_xlt0 = (x < setZero());
    mask_ylt0 = (y < setZero());

    aoffset = selectByNotMask(halfpi, mask_xnz);
    aoffset = selectByMask(aoffset, mask_ynz);

    aoffset = blend(aoffset, pi, mask_xlt0);
    aoffset = blend(aoffset, -aoffset, mask_ylt0);

    xinv = maskzInv(x, mask_xnz && (minFloat <= abs(x)));
    p    = y * xinv;
    p    = atan(p);
    p    = p + aoffset;

    return p;
}


/*! \brief Calculate the force correction due to PME analytically in SIMD double.
 *
 * \param z2 This should be the value \f$(r \beta)^2\f$, where r is your
 *           interaction distance and beta the ewald splitting parameters.
 * \result Correction factor to coulomb force.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables. For details, see the
 * single precision function.
 */
static inline SimdDouble gmx_simdcall pmeForceCorrection(SimdDouble z2)
{
    const SimdDouble FN10(-8.0072854618360083154e-14);
    const SimdDouble FN9(1.1859116242260148027e-11);
    const SimdDouble FN8(-8.1490406329798423616e-10);
    const SimdDouble FN7(3.4404793543907847655e-8);
    const SimdDouble FN6(-9.9471420832602741006e-7);
    const SimdDouble FN5(0.000020740315999115847456);
    const SimdDouble FN4(-0.00031991745139313364005);
    const SimdDouble FN3(0.0035074449373659008203);
    const SimdDouble FN2(-0.031750380176100813405);
    const SimdDouble FN1(0.13884101728898463426);
    const SimdDouble FN0(-0.75225277815249618847);

    const SimdDouble FD5(0.000016009278224355026701);
    const SimdDouble FD4(0.00051055686934806966046);
    const SimdDouble FD3(0.0081803507497974289008);
    const SimdDouble FD2(0.077181146026670287235);
    const SimdDouble FD1(0.41543303143712535988);
    const SimdDouble FD0(1.0);

    SimdDouble z4;
    SimdDouble polyFN0, polyFN1, polyFD0, polyFD1;

    z4 = z2 * z2;

    polyFD1 = fma(FD5, z4, FD3);
    polyFD1 = fma(polyFD1, z4, FD1);
    polyFD1 = polyFD1 * z2;
    polyFD0 = fma(FD4, z4, FD2);
    polyFD0 = fma(polyFD0, z4, FD0);
    polyFD0 = polyFD0 + polyFD1;

    polyFD0 = inv(polyFD0);

    polyFN0 = fma(FN10, z4, FN8);
    polyFN0 = fma(polyFN0, z4, FN6);
    polyFN0 = fma(polyFN0, z4, FN4);
    polyFN0 = fma(polyFN0, z4, FN2);
    polyFN0 = fma(polyFN0, z4, FN0);
    polyFN1 = fma(FN9, z4, FN7);
    polyFN1 = fma(polyFN1, z4, FN5);
    polyFN1 = fma(polyFN1, z4, FN3);
    polyFN1 = fma(polyFN1, z4, FN1);
    polyFN0 = fma(polyFN1, z2, polyFN0);


    return polyFN0 * polyFD0;
}


/*! \brief Calculate the potential correction due to PME analytically in SIMD double.
 *
 * \param z2 This should be the value \f$(r \beta)^2\f$, where r is your
 *           interaction distance and beta the ewald splitting parameters.
 * \result Correction factor to coulomb force.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic potential to avoid tables. For details, see the
 * single precision function.
 */
static inline SimdDouble gmx_simdcall pmePotentialCorrection(SimdDouble z2)
{
    const SimdDouble VN9(-9.3723776169321855475e-13);
    const SimdDouble VN8(1.2280156762674215741e-10);
    const SimdDouble VN7(-7.3562157912251309487e-9);
    const SimdDouble VN6(2.6215886208032517509e-7);
    const SimdDouble VN5(-4.9532491651265819499e-6);
    const SimdDouble VN4(0.00025907400778966060389);
    const SimdDouble VN3(0.0010585044856156469792);
    const SimdDouble VN2(0.045247661136833092885);
    const SimdDouble VN1(0.11643931522926034421);
    const SimdDouble VN0(1.1283791671726767970);

    const SimdDouble VD5(0.000021784709867336150342);
    const SimdDouble VD4(0.00064293662010911388448);
    const SimdDouble VD3(0.0096311444822588683504);
    const SimdDouble VD2(0.085608012351550627051);
    const SimdDouble VD1(0.43652499166614811084);
    const SimdDouble VD0(1.0);

    SimdDouble z4;
    SimdDouble polyVN0, polyVN1, polyVD0, polyVD1;

    z4 = z2 * z2;

    polyVD1 = fma(VD5, z4, VD3);
    polyVD0 = fma(VD4, z4, VD2);
    polyVD1 = fma(polyVD1, z4, VD1);
    polyVD0 = fma(polyVD0, z4, VD0);
    polyVD0 = fma(polyVD1, z2, polyVD0);

    polyVD0 = inv(polyVD0);

    polyVN1 = fma(VN9, z4, VN7);
    polyVN0 = fma(VN8, z4, VN6);
    polyVN1 = fma(polyVN1, z4, VN5);
    polyVN0 = fma(polyVN0, z4, VN4);
    polyVN1 = fma(polyVN1, z4, VN3);
    polyVN0 = fma(polyVN0, z4, VN2);
    polyVN1 = fma(polyVN1, z4, VN1);
    polyVN0 = fma(polyVN0, z4, VN0);
    polyVN0 = fma(polyVN1, z2, polyVN0);

    return polyVN0 * polyVD0;
}

/*! \} */


/*! \name SIMD math functions for double prec. data, single prec. accuracy
 *
 *  \note In some cases we do not need full double accuracy of individual
 *        SIMD math functions, although the data is stored in double precision
 *        SIMD registers. This might be the case for special algorithms, or
 *        if the architecture does not support single precision.
 *        Since the full double precision evaluation of math functions
 *        typically require much more expensive polynomial approximations
 *        these functions implement the algorithms used in the single precision
 *        SIMD math functions, but they operate on double precision
 *        SIMD variables.
 *
 *  \{
 */

/*********************************************************************
 * SIMD MATH FUNCTIONS WITH DOUBLE PREC. DATA, SINGLE PREC. ACCURACY *
 *********************************************************************/

/*! \brief Calculate 1/sqrt(x) for SIMD double, but in single accuracy.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall invsqrtSingleAccuracy(SimdDouble x)
{
    SimdDouble lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief 1/sqrt(x) for masked-in entries of SIMD double, but in single accuracy.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdDouble maskzInvsqrtSingleAccuracy(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = maskzRsqrt(x, m);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}

/*! \brief Calculate 1/sqrt(x) for two SIMD doubles, but single accuracy.
 *
 * \param x0  First set of arguments, x0 must be in single range (see below).
 * \param x1  Second set of arguments, x1 must be in single range (see below).
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 *
 * \note Both arguments must be larger than GMX_FLOAT_MIN and smaller than
 *       GMX_FLOAT_MAX, i.e. within the range of single precision.
 *       For the single precision implementation this is obviously always
 *       true for positive values, but for double precision it adds an
 *       extra restriction since the first lookup step might have to be
 *       performed in single precision on some architectures. Note that the
 *       responsibility for checking falls on you - this routine does not
 *       check arguments.
 */
static inline void gmx_simdcall invsqrtPairSingleAccuracy(SimdDouble  x0,
                                                          SimdDouble  x1,
                                                          SimdDouble* out0,
                                                          SimdDouble* out1)
{
#        if GMX_SIMD_HAVE_FLOAT && (GMX_SIMD_FLOAT_WIDTH == 2 * GMX_SIMD_DOUBLE_WIDTH) \
                && (GMX_SIMD_RSQRT_BITS < 22)
    SimdFloat  xf  = cvtDD2F(x0, x1);
    SimdFloat  luf = rsqrt(xf);
    SimdDouble lu0, lu1;
    // Intermediate target is single - mantissa+1 bits
#            if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
#            if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
#            if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    luf = rsqrtIter(luf, xf);
#            endif
    cvtF2DD(luf, &lu0, &lu1);
    // We now have single-precision accuracy values in lu0/lu1
    *out0 = lu0;
    *out1 = lu1;
#        else
    *out0 = invsqrtSingleAccuracy(x0);
    *out1 = invsqrtSingleAccuracy(x1);
#        endif
}

/*! \brief Calculate 1/x for SIMD double, but in single accuracy.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \return 1/x. Result is undefined if your argument was invalid.
 */
static inline SimdDouble gmx_simdcall invSingleAccuracy(SimdDouble x)
{
    SimdDouble lu = rcp(x);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}

/*! \brief 1/x for masked entries of SIMD double, single accuracy.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *
 *  \param m Mask
 *  \return 1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdDouble gmx_simdcall maskzInvSingleAccuracy(SimdDouble x, SimdDBool m)
{
    SimdDouble lu = maskzRcp(x, m);
#        if (GMX_SIMD_RCP_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
#        if (GMX_SIMD_RCP_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rcpIter(lu, x);
#        endif
    return lu;
}


/*! \brief Calculate sqrt(x) (correct for 0.0) for SIMD double, with single accuracy.
 *
 *  \copydetails sqrt(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall sqrtSingleAccuracy(SimdDouble x)
{
    if (opt == MathOptimization::Safe)
    {
        SimdDouble res = maskzInvsqrt(x, SimdDouble(GMX_FLOAT_MIN) < x);
        return res * x;
    }
    else
    {
        return x * invsqrtSingleAccuracy(x);
    }
}

/*! \brief Cube root for SIMD doubles, single accuracy.
 *
 * \param x      Argument to calculate cube root of. Can be negative or zero,
 *               but NaN or Inf values are not supported. Denormal values will
 *               be treated as 0.0.
 * \return       Cube root of x.
 */
static inline SimdDouble gmx_simdcall cbrtSingleAccuracy(SimdDouble x)
{
    const SimdDouble signBit(GMX_DOUBLE_NEGZERO);
    const SimdDouble minDouble(std::numeric_limits<double>::min());
    // Bias is 1024-1 = 1023, which is divisible by 3, so no need to change it more.
    // Use the divided value as original constant to avoid division warnings.
    const std::int32_t offsetDiv3(341);
    const SimdDouble   c2(-0.191502161678719066);
    const SimdDouble   c1(0.697570460207922770);
    const SimdDouble   c0(0.492659620528969547);
    const SimdDouble   one(1.0);
    const SimdDouble   two(2.0);
    const SimdDouble   three(3.0);
    const SimdDouble   oneThird(1.0 / 3.0);
    const SimdDouble   cbrt2(1.2599210498948731648);
    const SimdDouble   sqrCbrt2(1.5874010519681994748);

    // See the single precision routines for documentation of the algorithm

    SimdDouble xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdDouble xAbs     = andNot(signBit, x);    // select everthing but the sign bit => abs(x)
    SimdDBool  xIsNonZero = (minDouble <= xAbs); // treat denormals as 0

    SimdDInt32 exponent;
    SimdDouble y        = frexp(xAbs, &exponent);
    SimdDouble z        = fma(fma(y, c2, c1), y, c0);
    SimdDouble w        = z * z * z;
    SimdDouble nom      = z * fma(two, y, w);
    SimdDouble invDenom = inv(fma(two, w, y));

    SimdDouble offsetExp = cvtI2R(exponent) + SimdDouble(static_cast<double>(3 * offsetDiv3) + 0.1);
    SimdDouble offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division
    SimdDInt32 expDiv3   = cvtR2I(offsetExpDiv3 - SimdDouble(static_cast<double>(offsetDiv3)));
    SimdDouble remainder = offsetExp - offsetExpDiv3 * three;
    SimdDouble factor    = blend(one, cbrt2, SimdDouble(0.5) < remainder);
    factor               = blend(factor, sqrCbrt2, SimdDouble(1.5) < remainder);
    SimdDouble fraction  = (nom * invDenom * factor) ^ xSignBit;
    SimdDouble result    = selectByMask(ldexp(fraction, expDiv3), xIsNonZero);
    return result;
}

/*! \brief Inverse cube root for SIMD doubles, single accuracy.
 *
 * \param x      Argument to calculate cube root of. Can be positive or
 *               negative, but the magnitude cannot be lower than
 *               the smallest normal number.
 * \return       Cube root of x. Undefined for values that don't
 *               fulfill the restriction of abs(x) > minDouble.
 */
static inline SimdDouble gmx_simdcall invcbrtSingleAccuracy(SimdDouble x)
{
    const SimdDouble signBit(GMX_DOUBLE_NEGZERO);
    // Bias is 1024-1 = 1023, which is divisible by 3, so no need to change it more.
    // Use the divided value as original constant to avoid division warnings.
    const std::int32_t offsetDiv3(341);
    const SimdDouble   c2(-0.191502161678719066);
    const SimdDouble   c1(0.697570460207922770);
    const SimdDouble   c0(0.492659620528969547);
    const SimdDouble   one(1.0);
    const SimdDouble   two(2.0);
    const SimdDouble   three(3.0);
    const SimdDouble   oneThird(1.0 / 3.0);
    const SimdDouble   invCbrt2(1.0 / 1.2599210498948731648);
    const SimdDouble   invSqrCbrt2(1.0F / 1.5874010519681994748);

    // See the single precision routines for documentation of the algorithm

    SimdDouble xSignBit = x & signBit; // create bit mask where the sign bit is 1 for x elements < 0
    SimdDouble xAbs     = andNot(signBit, x); // select everthing but the sign bit => abs(x)

    SimdDInt32 exponent;
    SimdDouble y         = frexp(xAbs, &exponent);
    SimdDouble z         = fma(fma(y, c2, c1), y, c0);
    SimdDouble w         = z * z * z;
    SimdDouble nom       = fma(two, w, y);
    SimdDouble invDenom  = inv(z * fma(two, y, w));
    SimdDouble offsetExp = cvtI2R(exponent) + SimdDouble(static_cast<double>(3 * offsetDiv3) + 0.1);
    SimdDouble offsetExpDiv3 =
            trunc(offsetExp * oneThird); // important to truncate here to mimic integer division
    SimdDInt32 expDiv3   = cvtR2I(SimdDouble(static_cast<double>(offsetDiv3)) - offsetExpDiv3);
    SimdDouble remainder = offsetExpDiv3 * three - offsetExp;
    SimdDouble factor    = blend(one, invCbrt2, remainder < SimdDouble(-0.5));
    factor               = blend(factor, invSqrCbrt2, remainder < SimdDouble(-1.5));
    SimdDouble fraction  = (nom * invDenom * factor) ^ xSignBit;
    SimdDouble result    = ldexp(fraction, expDiv3);
    return result;
}

/*! \brief SIMD log2(x). Double precision SIMD data, single accuracy.
 *
 * \param x Argument, should be >0.
 * \result The base 2 logarithm of x. Undefined if argument is invalid.
 */
static inline SimdDouble gmx_simdcall log2SingleAccuracy(SimdDouble x)
{
#        if GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE
    return log(x) * SimdDouble(std::log2(std::exp(1.0)));
#        else
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble sqrt2(std::sqrt(2.0));
    const SimdDouble CL9(0.342149508897807708152F);
    const SimdDouble CL7(0.411570606888219447939F);
    const SimdDouble CL5(0.577085979152320294183F);
    const SimdDouble CL3(0.961796550607099898222F);
    const SimdDouble CL1(2.885390081777926774009F);
    SimdDouble       fexp, x2, p;
    SimdDInt32       iexp;
    SimdDBool        mask;

    x    = frexp(x, &iexp);
    fexp = cvtI2R(iexp);

    mask = (x < sqrt2);
    // Adjust to non-IEEE format for x<sqrt(2): exponent -= 1, mantissa *= 2.0
    fexp = fexp - selectByMask(one, mask);
    x    = x * blend(one, two, mask);

    x  = (x - one) * invSingleAccuracy(x + one);
    x2 = x * x;

    p = fma(CL9, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, fexp);

    return p;
#        endif
}

/*! \brief SIMD log(x). Double precision SIMD data, single accuracy.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
static inline SimdDouble gmx_simdcall logSingleAccuracy(SimdDouble x)
{
#        if GMX_SIMD_HAVE_NATIVE_LOG_DOUBLE
    return log(x);
#        else
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);
    const SimdDouble invsqrt2(1.0 / std::sqrt(2.0));
    const SimdDouble corr(0.693147180559945286226764);
    const SimdDouble CL9(0.2371599674224853515625);
    const SimdDouble CL7(0.285279005765914916992188);
    const SimdDouble CL5(0.400005519390106201171875);
    const SimdDouble CL3(0.666666567325592041015625);
    const SimdDouble CL1(2.0);
    SimdDouble       fexp, x2, p;
    SimdDInt32       iexp;
    SimdDBool        mask;

    x    = frexp(x, &iexp);
    fexp = cvtI2R(iexp);

    mask = x < invsqrt2;
    // Adjust to non-IEEE format for x<1/sqrt(2): exponent -= 1, mantissa *= 2.0
    fexp = fexp - selectByMask(one, mask);
    x    = x * blend(one, two, mask);

    x  = (x - one) * invSingleAccuracy(x + one);
    x2 = x * x;

    p = fma(CL9, x2, CL7);
    p = fma(p, x2, CL5);
    p = fma(p, x2, CL3);
    p = fma(p, x2, CL1);
    p = fma(p, x, corr * fexp);

    return p;
#        endif
}

/*! \brief SIMD 2^x. Double precision SIMD, single accuracy.
 *
 * \copydetails exp2(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall exp2SingleAccuracy(SimdDouble x)
{
#        if GMX_SIMD_HAVE_NATIVE_EXP2_DOUBLE
    return exp2(x);
#        else
    const SimdDouble CC6(0.0001534581200287996416911311);
    const SimdDouble CC5(0.001339993121934088894618990);
    const SimdDouble CC4(0.009618488957115180159497841);
    const SimdDouble CC3(0.05550328776964726865751735);
    const SimdDouble CC2(0.2402264689063408646490722);
    const SimdDouble CC1(0.6931472057372680777553816);
    const SimdDouble one(1.0);

    SimdDouble intpart;
    SimdDouble p;
    SimdDInt32 ix;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -1023, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer.
    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdDouble(std::numeric_limits<std::int32_t>::lowest()));
    }

    ix      = cvtR2I(x);
    intpart = round(x);
    x       = x - intpart;

    p = fma(CC6, x, CC5);
    p = fma(p, x, CC4);
    p = fma(p, x, CC3);
    p = fma(p, x, CC2);
    p = fma(p, x, CC1);
    p = fma(p, x, one);
    x = ldexp<opt>(p, ix);

    return x;
#        endif
}


/*! \brief SIMD exp(x). Double precision SIMD, single accuracy.
 *
 * \copydetails exp(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall expSingleAccuracy(SimdDouble x)
{
#        if GMX_SIMD_HAVE_NATIVE_EXP_DOUBLE
    return exp(x);
#        else
    const SimdDouble argscale(1.44269504088896341);
    // Lower bound: Clamp args that would lead to an IEEE fp exponent below -1023.
    const SimdDouble smallArgLimit(-709.0895657128);
    const SimdDouble invargscale(-0.69314718055994528623);
    const SimdDouble CC4(0.00136324646882712841033936);
    const SimdDouble CC3(0.00836596917361021041870117);
    const SimdDouble CC2(0.0416710823774337768554688);
    const SimdDouble CC1(0.166665524244308471679688);
    const SimdDouble CC0(0.499999850988388061523438);
    const SimdDouble one(1.0);
    SimdDouble       intpart;
    SimdDouble       y, p;
    SimdDInt32       iy;

    // Large negative values are valid arguments to exp2(), so there are two
    // things we need to account for:
    // 1. When the exponents reaches -1023, the (biased) exponent field will be
    //    zero and we can no longer multiply with it. There are special IEEE
    //    formats to handle this range, but for now we have to accept that
    //    we cannot handle those arguments. If input value becomes even more
    //    negative, it will start to loop and we would end up with invalid
    //    exponents. Thus, we need to limit or mask this.
    // 2. For VERY large negative values, we will have problems that the
    //    subtraction to get the fractional part loses accuracy, and then we
    //    can end up with overflows in the polynomial.
    //
    // For now, we handle this by forwarding the math optimization setting to
    // ldexp, where the routine will return zero for very small arguments.
    //
    // However, before doing that we need to make sure we do not call cvtR2I
    // with an argument that is so negative it cannot be converted to an integer
    // after being multiplied by argscale.

    if (opt == MathOptimization::Safe)
    {
        x = max(x, SimdDouble(std::numeric_limits<std::int32_t>::lowest()) / argscale);
    }

    y = x * argscale;

    iy      = cvtR2I(y);
    intpart = round(y); // use same rounding algorithm here

    // Extended precision arithmetics not needed since
    // we have double precision and only need single accuracy.
    x = fma(invargscale, intpart, x);

    p = fma(CC4, x, CC3);
    p = fma(p, x, CC2);
    p = fma(p, x, CC1);
    p = fma(p, x, CC0);
    p = fma(x * x, p, x);
    p = p + one;
    x = ldexp<opt>(p, iy);
    return x;
#        endif
}

/*! \brief SIMD pow(x,y). Double precision SIMD data, single accuracy.
 *
 * This returns x^y for SIMD values.
 *
 * \tparam opt If this is changed from the default (safe) into the unsafe
 *             option, there are no guarantees about correct results for x==0.
 *
 * \param x Base.
 *
 * \param y exponent.

 * \result x^y. Overflowing arguments are likely to either return 0 or inf,
 *         depending on the underlying implementation. If unsafe optimizations
 *         are enabled, this is also true for x==0.
 *
 * \warning You cannot rely on this implementation returning inf for arguments
 *          that cause overflow. If you have some very large
 *          values and need to rely on getting a valid numerical output,
 *          take the minimum of your variable and the largest valid argument
 *          before calling this routine.
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdDouble gmx_simdcall powSingleAccuracy(SimdDouble x, SimdDouble y)
{
    SimdDouble xcorr;

    if (opt == MathOptimization::Safe)
    {
        xcorr = max(x, SimdDouble(std::numeric_limits<double>::min()));
    }
    else
    {
        xcorr = x;
    }

    SimdDouble result = exp2SingleAccuracy<opt>(y * log2SingleAccuracy(xcorr));

    if (opt == MathOptimization::Safe)
    {
        // if x==0 and y>0 we explicitly set the result to 0.0
        // For any x with y==0, the result will already be 1.0 since we multiply by y (0.0) and call exp().
        result = blend(result, setZero(), x == setZero() && setZero() < y);
    }

    return result;
}

/*! \brief SIMD erf(x). Double precision SIMD data, single accuracy.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to single precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdDouble gmx_simdcall erfSingleAccuracy(SimdDouble x)
{
    // Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1]
    const SimdDouble CA6(7.853861353153693e-5);
    const SimdDouble CA5(-8.010193625184903e-4);
    const SimdDouble CA4(5.188327685732524e-3);
    const SimdDouble CA3(-2.685381193529856e-2);
    const SimdDouble CA2(1.128358514861418e-1);
    const SimdDouble CA1(-3.761262582423300e-1);
    const SimdDouble CA0(1.128379165726710);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2]
    const SimdDouble CB9(-0.0018629930017603923);
    const SimdDouble CB8(0.003909821287598495);
    const SimdDouble CB7(-0.0052094582210355615);
    const SimdDouble CB6(0.005685614362160572);
    const SimdDouble CB5(-0.0025367682853477272);
    const SimdDouble CB4(-0.010199799682318782);
    const SimdDouble CB3(0.04369575504816542);
    const SimdDouble CB2(-0.11884063474674492);
    const SimdDouble CB1(0.2732120154030589);
    const SimdDouble CB0(0.42758357702025784);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19]
    const SimdDouble CC10(-0.0445555913112064);
    const SimdDouble CC9(0.21376355144663348);
    const SimdDouble CC8(-0.3473187200259257);
    const SimdDouble CC7(0.016690861551248114);
    const SimdDouble CC6(0.7560973182491192);
    const SimdDouble CC5(-1.2137903600145787);
    const SimdDouble CC4(0.8411872321232948);
    const SimdDouble CC3(-0.08670413896296343);
    const SimdDouble CC2(-0.27124782687240334);
    const SimdDouble CC1(-0.0007502488047806069);
    const SimdDouble CC0(0.5642114853803148);
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);

    SimdDouble x2, x4, y;
    SimdDouble t, t2, w, w2;
    SimdDouble pA0, pA1, pB0, pB1, pC0, pC1;
    SimdDouble expmx2;
    SimdDouble res_erf, res_erfc, res;
    SimdDBool  mask, msk_erf;

    // Calculate erf()
    x2 = x * x;
    x4 = x2 * x2;

    pA0 = fma(CA6, x4, CA4);
    pA1 = fma(CA5, x4, CA3);
    pA0 = fma(pA0, x4, CA2);
    pA1 = fma(pA1, x4, CA1);
    pA0 = pA0 * x4;
    pA0 = fma(pA1, x2, pA0);
    // Constant term must come last for precision reasons
    pA0 = pA0 + CA0;

    res_erf = x * pA0;

    // Calculate erfc
    y       = abs(x);
    msk_erf = (SimdDouble(0.75) <= y);
    t       = maskzInvSingleAccuracy(y, msk_erf);
    w       = t - one;
    t2      = t * t;
    w2      = w * w;

    expmx2 = expSingleAccuracy(-y * y);

    pB1 = fma(CB9, w2, CB7);
    pB0 = fma(CB8, w2, CB6);
    pB1 = fma(pB1, w2, CB5);
    pB0 = fma(pB0, w2, CB4);
    pB1 = fma(pB1, w2, CB3);
    pB0 = fma(pB0, w2, CB2);
    pB1 = fma(pB1, w2, CB1);
    pB0 = fma(pB0, w2, CB0);
    pB0 = fma(pB1, w, pB0);

    pC0 = fma(CC10, t2, CC8);
    pC1 = fma(CC9, t2, CC7);
    pC0 = fma(pC0, t2, CC6);
    pC1 = fma(pC1, t2, CC5);
    pC0 = fma(pC0, t2, CC4);
    pC1 = fma(pC1, t2, CC3);
    pC0 = fma(pC0, t2, CC2);
    pC1 = fma(pC1, t2, CC1);

    pC0 = fma(pC0, t2, CC0);
    pC0 = fma(pC1, t, pC0);
    pC0 = pC0 * t;

    // Select pB0 or pC0 for erfc()
    mask     = (two < y);
    res_erfc = blend(pB0, pC0, mask);
    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    mask     = (x < setZero());
    res_erfc = blend(res_erfc, two - res_erfc, mask);

    // Select erf() or erfc()
    mask = (y < SimdDouble(0.75));
    res  = blend(one - res_erfc, res_erf, mask);

    return res;
}

/*! \brief SIMD erfc(x). Double precision SIMD data, single accuracy.
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves singleprecision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-30) since that is not
 * relevant for MD.
 */
static inline SimdDouble gmx_simdcall erfcSingleAccuracy(SimdDouble x)
{
    // Coefficients for minimax approximation of erf(x)=x*P(x^2) in range [-1,1]
    const SimdDouble CA6(7.853861353153693e-5);
    const SimdDouble CA5(-8.010193625184903e-4);
    const SimdDouble CA4(5.188327685732524e-3);
    const SimdDouble CA3(-2.685381193529856e-2);
    const SimdDouble CA2(1.128358514861418e-1);
    const SimdDouble CA1(-3.761262582423300e-1);
    const SimdDouble CA0(1.128379165726710);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*P((1/(x-1))^2) in range [0.67,2]
    const SimdDouble CB9(-0.0018629930017603923);
    const SimdDouble CB8(0.003909821287598495);
    const SimdDouble CB7(-0.0052094582210355615);
    const SimdDouble CB6(0.005685614362160572);
    const SimdDouble CB5(-0.0025367682853477272);
    const SimdDouble CB4(-0.010199799682318782);
    const SimdDouble CB3(0.04369575504816542);
    const SimdDouble CB2(-0.11884063474674492);
    const SimdDouble CB1(0.2732120154030589);
    const SimdDouble CB0(0.42758357702025784);
    // Coefficients for minimax approximation of erfc(x)=Exp(-x^2)*(1/x)*P((1/x)^2) in range [2,9.19]
    const SimdDouble CC10(-0.0445555913112064);
    const SimdDouble CC9(0.21376355144663348);
    const SimdDouble CC8(-0.3473187200259257);
    const SimdDouble CC7(0.016690861551248114);
    const SimdDouble CC6(0.7560973182491192);
    const SimdDouble CC5(-1.2137903600145787);
    const SimdDouble CC4(0.8411872321232948);
    const SimdDouble CC3(-0.08670413896296343);
    const SimdDouble CC2(-0.27124782687240334);
    const SimdDouble CC1(-0.0007502488047806069);
    const SimdDouble CC0(0.5642114853803148);
    const SimdDouble one(1.0);
    const SimdDouble two(2.0);

    SimdDouble x2, x4, y;
    SimdDouble t, t2, w, w2;
    SimdDouble pA0, pA1, pB0, pB1, pC0, pC1;
    SimdDouble expmx2;
    SimdDouble res_erf, res_erfc, res;
    SimdDBool  mask, msk_erf;

    // Calculate erf()
    x2 = x * x;
    x4 = x2 * x2;

    pA0 = fma(CA6, x4, CA4);
    pA1 = fma(CA5, x4, CA3);
    pA0 = fma(pA0, x4, CA2);
    pA1 = fma(pA1, x4, CA1);
    pA1 = pA1 * x2;
    pA0 = fma(pA0, x4, pA1);
    // Constant term must come last for precision reasons
    pA0 = pA0 + CA0;

    res_erf = x * pA0;

    // Calculate erfc
    y       = abs(x);
    msk_erf = (SimdDouble(0.75) <= y);
    t       = maskzInvSingleAccuracy(y, msk_erf);
    w       = t - one;
    t2      = t * t;
    w2      = w * w;

    expmx2 = expSingleAccuracy(-y * y);

    pB1 = fma(CB9, w2, CB7);
    pB0 = fma(CB8, w2, CB6);
    pB1 = fma(pB1, w2, CB5);
    pB0 = fma(pB0, w2, CB4);
    pB1 = fma(pB1, w2, CB3);
    pB0 = fma(pB0, w2, CB2);
    pB1 = fma(pB1, w2, CB1);
    pB0 = fma(pB0, w2, CB0);
    pB0 = fma(pB1, w, pB0);

    pC0 = fma(CC10, t2, CC8);
    pC1 = fma(CC9, t2, CC7);
    pC0 = fma(pC0, t2, CC6);
    pC1 = fma(pC1, t2, CC5);
    pC0 = fma(pC0, t2, CC4);
    pC1 = fma(pC1, t2, CC3);
    pC0 = fma(pC0, t2, CC2);
    pC1 = fma(pC1, t2, CC1);

    pC0 = fma(pC0, t2, CC0);
    pC0 = fma(pC1, t, pC0);
    pC0 = pC0 * t;

    // Select pB0 or pC0 for erfc()
    mask     = (two < y);
    res_erfc = blend(pB0, pC0, mask);
    res_erfc = res_erfc * expmx2;

    // erfc(x<0) = 2-erfc(|x|)
    mask     = (x < setZero());
    res_erfc = blend(res_erfc, two - res_erfc, mask);

    // Select erf() or erfc()
    mask = (y < SimdDouble(0.75));
    res  = blend(res_erfc, one - res_erf, mask);

    return res;
}

/*! \brief SIMD sin \& cos. Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 */
static inline void gmx_simdcall sinCosSingleAccuracy(SimdDouble x, SimdDouble* sinval, SimdDouble* cosval)
{
    // Constants to subtract Pi/4*x from y while minimizing precision loss
    const SimdDouble argred0(2 * 0.78539816290140151978);
    const SimdDouble argred1(2 * 4.9604678871439933374e-10);
    const SimdDouble argred2(2 * 1.1258708853173288931e-18);
    const SimdDouble two_over_pi(2.0 / M_PI);
    const SimdDouble const_sin2(-1.9515295891e-4);
    const SimdDouble const_sin1(8.3321608736e-3);
    const SimdDouble const_sin0(-1.6666654611e-1);
    const SimdDouble const_cos2(2.443315711809948e-5);
    const SimdDouble const_cos1(-1.388731625493765e-3);
    const SimdDouble const_cos0(4.166664568298827e-2);

    const SimdDouble half(0.5);
    const SimdDouble one(1.0);
    SimdDouble       ssign, csign;
    SimdDouble       x2, y, z, psin, pcos, sss, ccc;
    SimdDBool        mask;

#        if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    const SimdDInt32 ione(1);
    const SimdDInt32 itwo(2);
    SimdDInt32       iy;

    z  = x * two_over_pi;
    iy = cvtR2I(z);
    y  = round(z);

    mask  = cvtIB2B((iy & ione) == setZero());
    ssign = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), cvtIB2B((iy & itwo) == itwo));
    csign = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), cvtIB2B(((iy + ione) & itwo) == itwo));
#        else
    const SimdDouble quarter(0.25);
    const SimdDouble minusquarter(-0.25);
    SimdDouble       q;
    SimdDBool        m1, m2, m3;

    /* The most obvious way to find the arguments quadrant in the unit circle
     * to calculate the sign is to use integer arithmetic, but that is not
     * present in all SIMD implementations. As an alternative, we have devised a
     * pure floating-point algorithm that uses truncation for argument reduction
     * so that we get a new value 0<=q<1 over the unit circle, and then
     * do floating-point comparisons with fractions. This is likely to be
     * slightly slower (~10%) due to the longer latencies of floating-point, so
     * we only use it when integer SIMD arithmetic is not present.
     */
    ssign = x;
    x     = abs(x);
    // It is critical that half-way cases are rounded down
    z = fma(x, two_over_pi, half);
    y = trunc(z);
    q = z * quarter;
    q = q - trunc(q);
    /* z now starts at 0.0 for x=-pi/4 (although neg. values cannot occur), and
     * then increased by 1.0 as x increases by 2*Pi, when it resets to 0.0.
     * This removes the 2*Pi periodicity without using any integer arithmetic.
     * First check if y had the value 2 or 3, set csign if true.
     */
    q = q - half;
    /* If we have logical operations we can work directly on the signbit, which
     * saves instructions. Otherwise we need to represent signs as +1.0/-1.0.
     * Thus, if you are altering defines to debug alternative code paths, the
     * two GMX_SIMD_HAVE_LOGICAL sections in this routine must either both be
     * active or inactive - you will get errors if only one is used.
     */
#            if GMX_SIMD_HAVE_LOGICAL
    ssign = ssign & SimdDouble(GMX_DOUBLE_NEGZERO);
    csign = andNot(q, SimdDouble(GMX_DOUBLE_NEGZERO));
    ssign = ssign ^ csign;
#            else
    ssign = copysign(SimdDouble(1.0), ssign);
    csign = copysign(SimdDouble(1.0), q);
    csign = -csign;
    ssign = ssign * csign; // swap ssign if csign was set.
#            endif
    // Check if y had value 1 or 3 (remember we subtracted 0.5 from q)
    m1   = (q < minusquarter);
    m2   = (setZero() <= q);
    m3   = (q < quarter);
    m2   = m2 && m3;
    mask = m1 || m2;
    // where mask is FALSE, swap sign.
    csign   = csign * blend(SimdDouble(-1.0), one, mask);
#        endif
    x  = fnma(y, argred0, x);
    x  = fnma(y, argred1, x);
    x  = fnma(y, argred2, x);
    x2 = x * x;

    psin = fma(const_sin2, x2, const_sin1);
    psin = fma(psin, x2, const_sin0);
    psin = fma(psin, x * x2, x);
    pcos = fma(const_cos2, x2, const_cos1);
    pcos = fma(pcos, x2, const_cos0);
    pcos = fms(pcos, x2, half);
    pcos = fma(pcos, x2, one);

    sss = blend(pcos, psin, mask);
    ccc = blend(psin, pcos, mask);
    // See comment for GMX_SIMD_HAVE_LOGICAL section above.
#        if GMX_SIMD_HAVE_LOGICAL
    *sinval = sss ^ ssign;
    *cosval = ccc ^ csign;
#        else
    *sinval = sss * ssign;
    *cosval = ccc * csign;
#        endif
}

/*! \brief SIMD sin(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall sinSingleAccuracy(SimdDouble x)
{
    SimdDouble s, c;
    sinCosSingleAccuracy(x, &s, &c);
    return s;
}

/*! \brief SIMD cos(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdDouble gmx_simdcall cosSingleAccuracy(SimdDouble x)
{
    SimdDouble s, c;
    sinCosSingleAccuracy(x, &s, &c);
    return c;
}

/*! \brief SIMD tan(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdDouble gmx_simdcall tanSingleAccuracy(SimdDouble x)
{
    const SimdDouble argred0(2 * 0.78539816290140151978);
    const SimdDouble argred1(2 * 4.9604678871439933374e-10);
    const SimdDouble argred2(2 * 1.1258708853173288931e-18);
    const SimdDouble two_over_pi(2.0 / M_PI);
    const SimdDouble CT6(0.009498288995810566122993911);
    const SimdDouble CT5(0.002895755790837379295226923);
    const SimdDouble CT4(0.02460087336161924491836265);
    const SimdDouble CT3(0.05334912882656359828045988);
    const SimdDouble CT2(0.1333989091464957704418495);
    const SimdDouble CT1(0.3333307599244198227797507);

    SimdDouble x2, p, y, z;
    SimdDBool  mask;

#        if GMX_SIMD_HAVE_DINT32_ARITHMETICS && GMX_SIMD_HAVE_LOGICAL
    SimdDInt32 iy;
    SimdDInt32 ione(1);

    z    = x * two_over_pi;
    iy   = cvtR2I(z);
    y    = round(z);
    mask = cvtIB2B((iy & ione) == ione);

    x = fnma(y, argred0, x);
    x = fnma(y, argred1, x);
    x = fnma(y, argred2, x);
    x = selectByMask(SimdDouble(GMX_DOUBLE_NEGZERO), mask) ^ x;
#        else
    const SimdDouble quarter(0.25);
    const SimdDouble half(0.5);
    const SimdDouble threequarter(0.75);
    SimdDouble       w, q;
    SimdDBool        m1, m2, m3;

    w    = abs(x);
    z    = fma(w, two_over_pi, half);
    y    = trunc(z);
    q    = z * quarter;
    q    = q - trunc(q);
    m1   = (quarter <= q);
    m2   = (q < half);
    m3   = (threequarter <= q);
    m1   = m1 && m2;
    mask = m1 || m3;
    w    = fnma(y, argred0, w);
    w    = fnma(y, argred1, w);
    w    = fnma(y, argred2, w);

    w = blend(w, -w, mask);
    x = w * copysign(SimdDouble(1.0), x);
#        endif
    x2 = x * x;
    p  = fma(CT6, x2, CT5);
    p  = fma(p, x2, CT4);
    p  = fma(p, x2, CT3);
    p  = fma(p, x2, CT2);
    p  = fma(p, x2, CT1);
    p  = fma(x2, p * x, x);

    p = blend(p, maskzInvSingleAccuracy(p, mask), mask);
    return p;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdDouble gmx_simdcall asinSingleAccuracy(SimdDouble x)
{
    const SimdDouble limitlow(1e-4);
    const SimdDouble half(0.5);
    const SimdDouble one(1.0);
    const SimdDouble halfpi(M_PI / 2.0);
    const SimdDouble CC5(4.2163199048E-2);
    const SimdDouble CC4(2.4181311049E-2);
    const SimdDouble CC3(4.5470025998E-2);
    const SimdDouble CC2(7.4953002686E-2);
    const SimdDouble CC1(1.6666752422E-1);
    SimdDouble       xabs;
    SimdDouble       z, z1, z2, q, q1, q2;
    SimdDouble       pA, pB;
    SimdDBool        mask, mask2;

    xabs  = abs(x);
    mask  = (half < xabs);
    z1    = half * (one - xabs);
    mask2 = (xabs < one);
    q1    = z1 * maskzInvsqrtSingleAccuracy(z1, mask2);
    q2    = xabs;
    z2    = q2 * q2;
    z     = blend(z2, z1, mask);
    q     = blend(q2, q1, mask);

    z2 = z * z;
    pA = fma(CC5, z2, CC3);
    pB = fma(CC4, z2, CC2);
    pA = fma(pA, z2, CC1);
    pA = pA * z;
    z  = fma(pB, z2, pA);
    z  = fma(z, q, q);
    q2 = halfpi - z;
    q2 = q2 - z;
    z  = blend(z, q2, mask);

    mask = (limitlow < xabs);
    z    = blend(xabs, z, mask);
    z    = copysign(z, x);

    return z;
}

/*! \brief SIMD acos(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdDouble gmx_simdcall acosSingleAccuracy(SimdDouble x)
{
    const SimdDouble one(1.0);
    const SimdDouble half(0.5);
    const SimdDouble pi(M_PI);
    const SimdDouble halfpi(M_PI / 2.0);
    SimdDouble       xabs;
    SimdDouble       z, z1, z2, z3;
    SimdDBool        mask1, mask2, mask3;

    xabs  = abs(x);
    mask1 = (half < xabs);
    mask2 = (setZero() < x);

    z     = half * (one - xabs);
    mask3 = (xabs < one);
    z     = z * maskzInvsqrtSingleAccuracy(z, mask3);
    z     = blend(x, z, mask1);
    z     = asinSingleAccuracy(z);

    z2 = z + z;
    z1 = pi - z2;
    z3 = halfpi - z;
    z  = blend(z1, z2, mask2);
    z  = blend(z3, z, mask1);

    return z;
}

/*! \brief SIMD asin(x). Double precision SIMD data, single accuracy.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdDouble gmx_simdcall atanSingleAccuracy(SimdDouble x)
{
    const SimdDouble halfpi(M_PI / 2);
    const SimdDouble CA17(0.002823638962581753730774);
    const SimdDouble CA15(-0.01595690287649631500244);
    const SimdDouble CA13(0.04250498861074447631836);
    const SimdDouble CA11(-0.07489009201526641845703);
    const SimdDouble CA9(0.1063479334115982055664);
    const SimdDouble CA7(-0.1420273631811141967773);
    const SimdDouble CA5(0.1999269574880599975585);
    const SimdDouble CA3(-0.3333310186862945556640);
    SimdDouble       x2, x3, x4, pA, pB;
    SimdDBool        mask, mask2;

    mask  = (x < setZero());
    x     = abs(x);
    mask2 = (SimdDouble(1.0) < x);
    x     = blend(x, maskzInvSingleAccuracy(x, mask2), mask2);

    x2 = x * x;
    x3 = x2 * x;
    x4 = x2 * x2;
    pA = fma(CA17, x4, CA13);
    pB = fma(CA15, x4, CA11);
    pA = fma(pA, x4, CA9);
    pB = fma(pB, x4, CA7);
    pA = fma(pA, x4, CA5);
    pB = fma(pB, x4, CA3);
    pA = fma(pA, x2, pB);
    pA = fma(pA, x3, x);

    pA = blend(pA, halfpi - pA, mask2);
    pA = blend(pA, -pA, mask);

    return pA;
}

/*! \brief SIMD atan2(y,x). Double precision SIMD data, single accuracy.
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static inline SimdDouble gmx_simdcall atan2SingleAccuracy(SimdDouble y, SimdDouble x)
{
    const SimdDouble pi(M_PI);
    const SimdDouble halfpi(M_PI / 2.0);
    SimdDouble       xinv, p, aoffset;
    SimdDBool        mask_xnz, mask_ynz, mask_xlt0, mask_ylt0;

    mask_xnz  = x != setZero();
    mask_ynz  = y != setZero();
    mask_xlt0 = (x < setZero());
    mask_ylt0 = (y < setZero());

    aoffset = selectByNotMask(halfpi, mask_xnz);
    aoffset = selectByMask(aoffset, mask_ynz);

    aoffset = blend(aoffset, pi, mask_xlt0);
    aoffset = blend(aoffset, -aoffset, mask_ylt0);

    xinv = maskzInvSingleAccuracy(x, mask_xnz);
    p    = y * xinv;
    p    = atanSingleAccuracy(p);
    p    = p + aoffset;

    return p;
}

/*! \brief Analytical PME force correction, double SIMD data, single accuracy.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb force - see below for details.
 *
 * This routine is meant to enable analytical evaluation of the
 * direct-space PME electrostatic force to avoid tables.
 *
 * The direct-space potential should be \f$ \mbox{erfc}(\beta r)/r\f$, but there
 * are some problems evaluating that:
 *
 * First, the error function is difficult (read: expensive) to
 * approxmiate accurately for intermediate to large arguments, and
 * this happens already in ranges of \f$(\beta r)\f$ that occur in simulations.
 * Second, we now try to avoid calculating potentials in Gromacs but
 * use forces directly.
 *
 * We can simply things slight by noting that the PME part is really
 * a correction to the normal Coulomb force since \f$\mbox{erfc}(z)=1-\mbox{erf}(z)\f$, i.e.
 * \f[
 * V = \frac{1}{r} - \frac{\mbox{erf}(\beta r)}{r}
 * \f]
 * The first term we already have from the inverse square root, so
 * that we can leave out of this routine.
 *
 * For pme tolerances of 1e-3 to 1e-8 and cutoffs of 0.5nm to 1.8nm,
 * the argument \f$beta r\f$ will be in the range 0.15 to ~4. Use your
 * favorite plotting program to realize how well-behaved \f$\frac{\mbox{erf}(z)}{z}\f$ is
 * in this range!
 *
 * We approximate \f$f(z)=\mbox{erf}(z)/z\f$ with a rational minimax polynomial.
 * However, it turns out it is more efficient to approximate \f$f(z)/z\f$ and
 * then only use even powers. This is another minor optimization, since
 * we actually \a want \f$f(z)/z\f$, because it is going to be multiplied by
 * the vector between the two atoms to get the vectorial force. The
 * fastest flops are the ones we can avoid calculating!
 *
 * So, here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=(\beta r)^2\f$.
 * 3. Evaluate this routine with \f$z^2\f$ as the argument.
 * 4. The return value is the expression:
 *
 * \f[
 *    \frac{2 \exp{-z^2}}{\sqrt{\pi} z^2}-\frac{\mbox{erf}(z)}{z^3}
 * \f]
 *
 * 5. Multiply the entire expression by \f$\beta^3\f$. This will get you
 *
 *  \f[
 *    \frac{2 \beta^3 \exp(-z^2)}{\sqrt{\pi} z^2} - \frac{\beta^3 \mbox{erf}(z)}{z^3}
 *  \f]
 *
 *    or, switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{2 \beta \exp(-r^2 \beta^2)}{\sqrt{\pi} r^2} - \frac{\mbox{erf}(r \beta)}{r^3}
 *  \f]
 *
 *    With a bit of math exercise you should be able to confirm that
 *    this is exactly
 *
 *  \f[
 *   \frac{\frac{d}{dr}\left( \frac{\mbox{erf}(\beta r)}{r} \right)}{r}
 *  \f]
 *
 * 6. Add the result to \f$r^{-3}\f$, multiply by the product of the charges,
 *    and you have your force (divided by \f$r\f$). A final multiplication
 *    with the vector connecting the two particles and you have your
 *    vectorial force to add to the particles.
 *
 * This approximation achieves an accuracy slightly lower than 1e-6; when
 * added to \f$1/r\f$ the error will be insignificant.
 *
 */
static inline SimdDouble gmx_simdcall pmeForceCorrectionSingleAccuracy(SimdDouble z2)
{
    const SimdDouble FN6(-1.7357322914161492954e-8);
    const SimdDouble FN5(1.4703624142580877519e-6);
    const SimdDouble FN4(-0.000053401640219807709149);
    const SimdDouble FN3(0.0010054721316683106153);
    const SimdDouble FN2(-0.019278317264888380590);
    const SimdDouble FN1(0.069670166153766424023);
    const SimdDouble FN0(-0.75225204789749321333);

    const SimdDouble FD4(0.0011193462567257629232);
    const SimdDouble FD3(0.014866955030185295499);
    const SimdDouble FD2(0.11583842382862377919);
    const SimdDouble FD1(0.50736591960530292870);
    const SimdDouble FD0(1.0);

    SimdDouble z4;
    SimdDouble polyFN0, polyFN1, polyFD0, polyFD1;

    z4 = z2 * z2;

    polyFD0 = fma(FD4, z4, FD2);
    polyFD1 = fma(FD3, z4, FD1);
    polyFD0 = fma(polyFD0, z4, FD0);
    polyFD0 = fma(polyFD1, z2, polyFD0);

    polyFD0 = invSingleAccuracy(polyFD0);

    polyFN0 = fma(FN6, z4, FN4);
    polyFN1 = fma(FN5, z4, FN3);
    polyFN0 = fma(polyFN0, z4, FN2);
    polyFN1 = fma(polyFN1, z4, FN1);
    polyFN0 = fma(polyFN0, z4, FN0);
    polyFN0 = fma(polyFN1, z2, polyFN0);

    return polyFN0 * polyFD0;
}


/*! \brief Analytical PME potential correction, double SIMD data, single accuracy.
 *
 * \param z2 \f$(r \beta)^2\f$ - see below for details.
 * \result Correction factor to coulomb potential - see below for details.
 *
 * This routine calculates \f$\mbox{erf}(z)/z\f$, although you should provide \f$z^2\f$
 * as the input argument.
 *
 * Here's how it should be used:
 *
 * 1. Calculate \f$r^2\f$.
 * 2. Multiply by \f$\beta^2\f$, so you get \f$z^2=\beta^2*r^2\f$.
 * 3. Evaluate this routine with z^2 as the argument.
 * 4. The return value is the expression:
 *
 *  \f[
 *   \frac{\mbox{erf}(z)}{z}
 *  \f]
 *
 * 5. Multiply the entire expression by beta and switching back to \f$r\f$ (since \f$z=r \beta\f$):
 *
 *  \f[
 *    \frac{\mbox{erf}(r \beta)}{r}
 *  \f]
 *
 * 6. Subtract the result from \f$1/r\f$, multiply by the product of the charges,
 *    and you have your potential.
 *
 * This approximation achieves an accuracy slightly lower than 1e-6; when
 * added to \f$1/r\f$ the error will be insignificant.
 */
static inline SimdDouble gmx_simdcall pmePotentialCorrectionSingleAccuracy(SimdDouble z2)
{
    const SimdDouble VN6(1.9296833005951166339e-8);
    const SimdDouble VN5(-1.4213390571557850962e-6);
    const SimdDouble VN4(0.000041603292906656984871);
    const SimdDouble VN3(-0.00013134036773265025626);
    const SimdDouble VN2(0.038657983986041781264);
    const SimdDouble VN1(0.11285044772717598220);
    const SimdDouble VN0(1.1283802385263030286);

    const SimdDouble VD3(0.0066752224023576045451);
    const SimdDouble VD2(0.078647795836373922256);
    const SimdDouble VD1(0.43336185284710920150);
    const SimdDouble VD0(1.0);

    SimdDouble z4;
    SimdDouble polyVN0, polyVN1, polyVD0, polyVD1;

    z4 = z2 * z2;

    polyVD1 = fma(VD3, z4, VD1);
    polyVD0 = fma(VD2, z4, VD0);
    polyVD0 = fma(polyVD1, z2, polyVD0);

    polyVD0 = invSingleAccuracy(polyVD0);

    polyVN0 = fma(VN6, z4, VN4);
    polyVN1 = fma(VN5, z4, VN3);
    polyVN0 = fma(polyVN0, z4, VN2);
    polyVN1 = fma(polyVN1, z4, VN1);
    polyVN0 = fma(polyVN0, z4, VN0);
    polyVN0 = fma(polyVN1, z2, polyVN0);

    return polyVN0 * polyVD0;
}

#    endif


/*! \name SIMD4 math functions
 *
 * \note Only a subset of the math functions are implemented for SIMD4.
 *  \{
 */


#    if GMX_SIMD4_HAVE_FLOAT

/*************************************************************************
 * SINGLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 float.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline Simd4Float gmx_simdcall rsqrtIter(Simd4Float lu, Simd4Float x)
{
    Simd4Float tmp1 = x * lu;
    Simd4Float tmp2 = Simd4Float(-0.5F) * lu;
    tmp1            = fma(tmp1, lu, Simd4Float(-3.0F));
    return tmp1 * tmp2;
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 float.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return  1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline Simd4Float gmx_simdcall invsqrt(Simd4Float x)
{
    Simd4Float lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}


#    endif // GMX_SIMD4_HAVE_FLOAT


#    if GMX_SIMD4_HAVE_DOUBLE
/*************************************************************************
 * DOUBLE PRECISION SIMD4 MATH FUNCTIONS - JUST A SMALL SUBSET SUPPORTED *
 *************************************************************************/

/*! \brief Perform one Newton-Raphson iteration to improve 1/sqrt(x) for SIMD4 double.
 *
 * This is a low-level routine that should only be used by SIMD math routine
 * that evaluates the inverse square root.
 *
 *  \param lu Approximation of 1/sqrt(x), typically obtained from lookup.
 *  \param x  The reference (starting) value x for which we want 1/sqrt(x).
 *  \return   An improved approximation with roughly twice as many bits of accuracy.
 */
static inline Simd4Double gmx_simdcall rsqrtIter(Simd4Double lu, Simd4Double x)
{
    Simd4Double tmp1 = x * lu;
    Simd4Double tmp2 = Simd4Double(-0.5F) * lu;
    tmp1             = fma(tmp1, lu, Simd4Double(-3.0F));
    return tmp1 * tmp2;
}

/*! \brief Calculate 1/sqrt(x) for SIMD4 double.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return  1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline Simd4Double gmx_simdcall invsqrt(Simd4Double x)
{
    Simd4Double lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 8 < GMX_SIMD_ACCURACY_BITS_DOUBLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}


/**********************************************************************
 * SIMD4 MATH FUNCTIONS WITH DOUBLE PREC. DATA, SINGLE PREC. ACCURACY *
 **********************************************************************/

/*! \brief Calculate 1/sqrt(x) for SIMD4 double, but in single accuracy.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return  1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline Simd4Double gmx_simdcall invsqrtSingleAccuracy(Simd4Double x)
{
    Simd4Double lu = rsqrt(x);
#        if (GMX_SIMD_RSQRT_BITS < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 2 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
#        if (GMX_SIMD_RSQRT_BITS * 4 < GMX_SIMD_ACCURACY_BITS_SINGLE)
    lu = rsqrtIter(lu, x);
#        endif
    return lu;
}


#    endif // GMX_SIMD4_HAVE_DOUBLE

/*! \} */

#    if GMX_SIMD_HAVE_FLOAT
/*! \brief Calculate 1/sqrt(x) for SIMD float, only targeting single accuracy.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return  1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall invsqrtSingleAccuracy(SimdFloat x)
{
    return invsqrt(x);
}

/*! \brief Calculate 1/sqrt(x) for masked SIMD floats, only targeting single accuracy.
 *
 *  This routine only evaluates 1/sqrt(x) for elements for which mask is true.
 *  Illegal values in the masked-out elements will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \param m Mask
 *  \return  1/sqrt(x). Result is undefined if your argument was invalid or
 *           entry was not masked, and 0.0 for masked-out entries.
 */
static inline SimdFloat maskzInvsqrtSingleAccuracy(SimdFloat x, SimdFBool m)
{
    return maskzInvsqrt(x, m);
}

/*! \brief Calculate 1/sqrt(x) for two SIMD floats, only targeting single accuracy.
 *
 * \param x0  First set of arguments, x0 must be in single range (see below).
 * \param x1  Second set of arguments, x1 must be in single range (see below).
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 *  In particular for double precision we can sometimes calculate square root
 *  pairs slightly faster by using single precision until the very last step.
 *
 * \note Both arguments must be larger than GMX_FLOAT_MIN and smaller than
 *       GMX_FLOAT_MAX, i.e. within the range of single precision.
 *       For the single precision implementation this is obviously always
 *       true for positive values, but for double precision it adds an
 *       extra restriction since the first lookup step might have to be
 *       performed in single precision on some architectures. Note that the
 *       responsibility for checking falls on you - this routine does not
 *       check arguments.
 */
static inline void gmx_simdcall invsqrtPairSingleAccuracy(SimdFloat x0, SimdFloat x1, SimdFloat* out0, SimdFloat* out1)
{
    return invsqrtPair(x0, x1, out0, out1);
}

/*! \brief Calculate 1/x for SIMD float, only targeting single accuracy.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return  1/x. Result is undefined if your argument was invalid.
 */
static inline SimdFloat gmx_simdcall invSingleAccuracy(SimdFloat x)
{
    return inv(x);
}


/*! \brief Calculate 1/x for masked SIMD floats, only targeting single accuracy.
 *
 *  \param x Argument with magnitude larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \param m Mask
 *  \return  1/x for elements where m is true, or 0.0 for masked-out entries.
 */
static inline SimdFloat maskzInvSingleAccuracy(SimdFloat x, SimdFBool m)
{
    return maskzInv(x, m);
}

/*! \brief Calculate sqrt(x) for SIMD float, always targeting single accuracy.
 *
 * \copydetails sqrt(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall sqrtSingleAccuracy(SimdFloat x)
{
    return sqrt<opt>(x);
}

/*! \brief Calculate cbrt(x) for SIMD float, always targeting single accuracy.
 *
 * \copydetails cbrt(SimdFloat)
 */
static inline SimdFloat gmx_simdcall cbrtSingleAccuracy(SimdFloat x)
{
    return cbrt(x);
}

/*! \brief Calculate 1/cbrt(x) for SIMD float, always targeting single accuracy.
 *
 * \copydetails cbrt(SimdFloat)
 */
static inline SimdFloat gmx_simdcall invcbrtSingleAccuracy(SimdFloat x)
{
    return invcbrt(x);
}

/*! \brief SIMD float log2(x), only targeting single accuracy. This is the base-2 logarithm.
 *
 * \param x Argument, should be >0.
 * \result The base-2 logarithm of x. Undefined if argument is invalid.
 */
static inline SimdFloat gmx_simdcall log2SingleAccuracy(SimdFloat x)
{
    return log2(x);
}

/*! \brief SIMD float log(x), only targeting single accuracy. This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 */
static inline SimdFloat gmx_simdcall logSingleAccuracy(SimdFloat x)
{
    return log(x);
}

/*! \brief SIMD float 2^x, only targeting single accuracy.
 *
 * \copydetails exp2(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall exp2SingleAccuracy(SimdFloat x)
{
    return exp2<opt>(x);
}

/*! \brief SIMD float e^x, only targeting single accuracy.
 *
 * \copydetails exp(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall expSingleAccuracy(SimdFloat x)
{
    return exp<opt>(x);
}

/*! \brief SIMD pow(x,y), only targeting single accuracy.
 *
 * \copydetails pow(SimdFloat)
 */
template<MathOptimization opt = MathOptimization::Safe>
static inline SimdFloat gmx_simdcall powSingleAccuracy(SimdFloat x, SimdFloat y)
{
    return pow<opt>(x, y);
}

/*! \brief SIMD float erf(x), only targeting single accuracy.
 *
 * \param x The value to calculate erf(x) for.
 * \result erf(x)
 *
 * This routine achieves very close to single precision, but we do not care about
 * the last bit or the subnormal result range.
 */
static inline SimdFloat gmx_simdcall erfSingleAccuracy(SimdFloat x)
{
    return erf(x);
}

/*! \brief SIMD float erfc(x), only targeting single accuracy.
 *
 * \param x The value to calculate erfc(x) for.
 * \result erfc(x)
 *
 * This routine achieves singleprecision (bar the last bit) over most of the
 * input range, but for large arguments where the result is getting close
 * to the minimum representable numbers we accept slightly larger errors
 * (think results that are in the ballpark of 10^-30) since that is not
 * relevant for MD.
 */
static inline SimdFloat gmx_simdcall erfcSingleAccuracy(SimdFloat x)
{
    return erfc(x);
}

/*! \brief SIMD float sin \& cos, only targeting single accuracy.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 */
static inline void gmx_simdcall sinCosSingleAccuracy(SimdFloat x, SimdFloat* sinval, SimdFloat* cosval)
{
    sincos(x, sinval, cosval);
}

/*! \brief SIMD float sin(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall sinSingleAccuracy(SimdFloat x)
{
    return sin(x);
}

/*! \brief SIMD float cos(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \attention Do NOT call both sin & cos if you need both results, since each of them
 * will then call \ref sincos and waste a factor 2 in performance.
 */
static inline SimdFloat gmx_simdcall cosSingleAccuracy(SimdFloat x)
{
    return cos(x);
}

/*! \brief SIMD float tan(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 */
static inline SimdFloat gmx_simdcall tanSingleAccuracy(SimdFloat x)
{
    return tan(x);
}

/*! \brief SIMD float asin(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 */
static inline SimdFloat gmx_simdcall asinSingleAccuracy(SimdFloat x)
{
    return asin(x);
}

/*! \brief SIMD float acos(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 */
static inline SimdFloat gmx_simdcall acosSingleAccuracy(SimdFloat x)
{
    return acos(x);
}

/*! \brief SIMD float atan(x), only targeting single accuracy.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x), same argument/value range as standard math library.
 */
static inline SimdFloat gmx_simdcall atanSingleAccuracy(SimdFloat x)
{
    return atan(x);
}

/*! \brief SIMD float atan2(y,x), only targeting single accuracy.
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x), same argument/value range as standard math library.
 *
 * \note This routine should provide correct results for all finite
 * non-zero or positive-zero arguments. However, negative zero arguments will
 * be treated as positive zero, which means the return value will deviate from
 * the standard math library atan2(y,x) for those cases. That should not be
 * of any concern in Gromacs, and in particular it will not affect calculations
 * of angles from vectors.
 */
static inline SimdFloat gmx_simdcall atan2SingleAccuracy(SimdFloat y, SimdFloat x)
{
    return atan2(y, x);
}

/*! \brief SIMD Analytic PME force correction, only targeting single accuracy.
 *
 * \param z2 \f$(r \beta)^2\f$ - see default single precision version for details.
 * \result Correction factor to coulomb force.
 */
static inline SimdFloat gmx_simdcall pmeForceCorrectionSingleAccuracy(SimdFloat z2)
{
    return pmeForceCorrection(z2);
}

/*! \brief SIMD Analytic PME potential correction, only targeting single accuracy.
 *
 * \param z2 \f$(r \beta)^2\f$ - see default single precision version for details.
 * \result Correction factor to coulomb force.
 */
static inline SimdFloat gmx_simdcall pmePotentialCorrectionSingleAccuracy(SimdFloat z2)
{
    return pmePotentialCorrection(z2);
}
#    endif // GMX_SIMD_HAVE_FLOAT

#    if GMX_SIMD4_HAVE_FLOAT
/*! \brief Calculate 1/sqrt(x) for SIMD4 float, only targeting single accuracy.
 *
 *  \param x Argument that must be larger than GMX_FLOAT_MIN and smaller than
 *           GMX_FLOAT_MAX, i.e. within the range of single precision.
 *           For the single precision implementation this is obviously always
 *           true for positive values, but for double precision it adds an
 *           extra restriction since the first lookup step might have to be
 *           performed in single precision on some architectures. Note that the
 *           responsibility for checking falls on you - this routine does not
 *           check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 */
static inline Simd4Float gmx_simdcall invsqrtSingleAccuracy(Simd4Float x)
{
    return invsqrt(x);
}
#    endif // GMX_SIMD4_HAVE_FLOAT

/*! \}   end of addtogroup module_simd */
/*! \endcond  end of condition libabl */

#endif // GMX_SIMD

} // namespace gmx

#endif // GMX_SIMD_SIMD_MATH_H
