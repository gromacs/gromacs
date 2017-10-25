/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_SIMD_SCALAR_H
#define GMX_SIMD_SCALAR_H

#include <cmath>
#include <cstdint>
#include <cstdlib>

#include <algorithm>

/*! \libinternal \file
 *
 * \brief Scalar float functions corresponding to GROMACS SIMD functions
 *
 * These versions make it possible to write functions that are templated with
 * either a SIMD or scalar type. While some of these functions might not appear
 * SIMD-specific, we have placed them here because the only reason to use these
 * instead of generic function is in templated combined SIMD/non-SIMD code.
 *
 * There are a handful of limitations, in particular that it is impossible
 * to overload the bitwise logical operators on built-in types.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

namespace gmx
{

/************************************************************************
 *   Single-precision floating point functions mimicking SIMD versions  *
 ************************************************************************/

/*! \brief Store contents of float variable to aligned memory m.
 *
 * \param[out] m Pointer to memory.
 * \param a float variable to store
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
store(float *m, float a)
{
    *m = a;
}

/*! \brief Store contents of float variable to unaligned memory m.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a float variable to store.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
storeU(float *m, float a)
{
    *m = a;
}

// We cannot overload the logical operators and, or, andNot, xor for
// built-in types.

/*! \brief Float Fused-multiply-add. Result is a*b + c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b + c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
fma(float a, float b, float c)
{
    // Note that we purposely do not use the single-rounding std::fma
    // as that can be very slow without hardware support
    return a*b + c;
}

/*! \brief Float Fused-multiply-subtract. Result is a*b - c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b - c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
fms(float a, float b, float c)
{
    return a*b - c;
}

/*! \brief Float Fused-negated-multiply-add. Result is -a*b + c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b + c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
fnma(float a, float b, float c)
{
    return c - a*b;
}

/*! \brief Float Fused-negated-multiply-subtract. Result is -a*b - c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b - c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
fnms(float a, float b, float c)
{
    return -a*b - c;
}

/*! \brief Add two float variables, masked version.
 *
 * \param a term1
 * \param b term2
 * \param m mask
 * \return a+b where mask is true, a otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
maskAdd(float a, float b, float m)
{
    return a + (m ? b : 0.0f);
}

/*! \brief Multiply two float variables, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param m mask
 * \return a*b where mask is true, 0.0 otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
maskzMul(float a, float b, float m)
{
    return m ? (a * b) : 0.0f;
}

/*! \brief Float fused multiply-add, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \param m mask
 * \return a*b+c where mask is true, 0.0 otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
maskzFma(float a, float b, float c, float m)
{
    return m ? (a * b + c) : 0.0f;
}

/*! \brief Float Floating-point abs().
 *
 * \param a any floating point values
 * \return abs(a) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
abs(float a)
{
    return std::abs(a);
}

/*! \brief Set each float element to the largest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
max(float a, float b)
{
    return std::max(a, b);
}

/*! \brief Set each float element to the smallest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return min(a,b) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
min(float a, float b)
{
    return std::min(a, b);
}

/*! \brief Float round to nearest integer value (in floating-point format).
 *
 * \param a Any floating-point value
 * \return The nearest integer, represented in floating-point format.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
round(float a)
{
    return std::round(a);
}

/*! \brief Truncate float, i.e. round towards zero - common hardware instruction.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
trunc(float a)
{
    return std::trunc(a);
}

/*! \brief Return sum of all elements in float variable (i.e., the variable itself).
 *
 * \param a variable to reduce/sum.
 * \return The argument variable itself.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
reduce(float a)
{
    return a;
}

/*! \brief Bitwise andnot for two scalar float variables.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
andNot(float a, float b)
{
    union
    {
        float         r;
        std::uint32_t i;
    } conv1, conv2;

    conv1.r = a;
    conv2.r = b;

    conv1.i = (~conv1.i) & conv2.i;

    return conv1.r;
}

/*! \brief Return true if any bits are set in the float variable.
 *
 * This function is used to handle bitmasks, mainly for exclusions in the
 * inner kernels. Note that it will return true even for -0.0f (sign bit set),
 * so it is not identical to not-equal.
 *
 * \param a value
 * \return True if any bit in a is nonzero.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
testBits(float a)
{
    union
    {
        std::uint32_t i;
        float         f;
    } conv;

    conv.f = a;
    return (conv.i != 0);
}

/*! \brief Returns if the boolean is true.
 *
 * \param a Logical variable.
 * \return true if a is true, otherwise false.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
anyTrue(bool a)
{
    return a;
}

/*! \brief Select from single precision variable where boolean is true.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  a is selected for true, 0 for false.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
selectByMask(float a, bool mask)
{
    return mask ? a : 0.0f;
}

/*! \brief Select from single precision variable where boolean is false.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  a is selected for false, 0 for true.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
selectByNotMask(float a, bool mask)
{
    return mask ? 0.0f : a;
}

/*! \brief Blend float selection.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return Select b if sel is true, a otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
blend(float a, float b, bool sel)
{
    return sel ? b : a;
}

/*! \brief Round single precision floating point to integer.
 *
 * \param a float
 * \return Integer format, a rounded to nearest integer.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
cvtR2I(float a)
{
    return std::round(a);
};

/*! \brief Truncate single precision floating point to integer.
 *
 * \param a float
 * \return Integer format, a truncated to integer.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
cvttR2I(float a)
{
    return std::trunc(a);
};

/*! \brief Return integer.
 *
 * This function mimicks the SIMD integer-to-real conversion routines. By
 * simply returning an integer, we let the compiler sort out whether the
 * conversion should be to float or double rather than using proxy objects.
 *
 * \param a integer
 * \return same value (a)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
cvtI2R(std::int32_t a)
{
    return a;
}

/************************************************************************
 *   Double-precision floating point functions mimicking SIMD versions  *
 ************************************************************************/

/*! \brief Store contents of double variable to aligned memory m.
 *
 * \param[out] m Pointer to memory.
 * \param a double variable to store
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
store(double *m, double a)
{
    *m = a;
}

/*! \brief Store contents of double variable to unaligned memory m.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a double variable to store.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
storeU(double *m, double a)
{
    *m = a;
}

// We cannot overload the logical operators and, or, andNot, xor for
// built-in types.

/*! \brief double Fused-multiply-add. Result is a*b + c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b + c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
fma(double a, double b, double c)
{
    // Note that we purposely do not use the single-rounding std::fma
    // as that can be very slow without hardware support
    return a*b + c;
}

/*! \brief double Fused-multiply-subtract. Result is a*b - c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b - c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
fms(double a, double b, double c)
{
    return a*b - c;
}

/*! \brief double Fused-negated-multiply-add. Result is - a*b + c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b + c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
fnma(double a, double b, double c)
{
    return c - a*b;
}

/*! \brief double Fused-negated-multiply-subtract. Result is -a*b - c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b - c
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
fnms(double a, double b, double c)
{
    return -a*b - c;
}

/*! \brief Add two double variables, masked version.
 *
 * \param a term1
 * \param b term2
 * \param m mask
 * \return a+b where mask is true, a otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskAdd(double a, double b, double m)
{
    return a + (m ? b : 0.0);
}

/*! \brief Multiply two double variables, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param m mask
 * \return a*b where mask is true, 0.0 otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzMul(double a, double b, double m)
{
    return m ? (a * b) : 0.0;
}

/*! \brief double fused multiply-add, masked version.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \param m mask
 * \return a*b+c where mask is true, 0.0 otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzFma(double a, double b, double c, double m)
{
    return m ? (a * b + c) : 0.0;
}

/*! \brief double doubleing-point abs().
 *
 * \param a any doubleing point values
 * \return abs(a) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
abs(double a)
{
    return std::abs(a);
}

/*! \brief Set each double element to the largest from two variables.
 *
 * \param a Any doubleing-point value
 * \param b Any doubleing-point value
 * \return max(a,b) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
max(double a, double b)
{
    return std::max(a, b);
}

/*! \brief Set each double element to the smallest from two variables.
 *
 * \param a Any doubleing-point value
 * \param b Any doubleing-point value
 * \return min(a,b) for each element.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
min(double a, double b)
{
    return std::min(a, b);
}

/*! \brief double round to nearest integer value (in doubleing-point format).
 *
 * \param a Any doubleing-point value
 * \return The nearest integer, represented in doubleing-point format.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
round(double a)
{
    return std::round(a);
}

/*! \brief Truncate double, i.e. round towards zero - common hardware instruction.
 *
 * \param a Any doubleing-point value
 * \return Integer rounded towards zero, represented in doubleing-point format.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
trunc(double a)
{
    return std::trunc(a);
}

/*! \brief Return sum of all elements in double variable (i.e., the variable itself).
 *
 * \param a variable to reduce/sum.
 * \return The argument variable itself.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
reduce(double a)
{
    return a;
}

/*! \brief Bitwise andnot for two scalar double variables.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
andNot(double a, double b)
{
    union
    {
        double        r;
        std::uint64_t i;
    } conv1, conv2;

    conv1.r = a;
    conv2.r = b;

    conv1.i = (~conv1.i) & conv2.i;

    return conv1.r;
}

/*! \brief Return true if any bits are set in the double variable.
 *
 * This function is used to handle bitmasks, mainly for exclusions in the
 * inner kernels. Note that it will return true even for -0.0 (sign bit set),
 * so it is not identical to not-equal.
 *
 * \param a value
 * \return True if any bit in a is nonzero.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
testBits(double a)
{
    union
    {
        std::uint64_t  i;
        double         f;
    } conv;

    conv.f = a;
    return (conv.i != 0);
}

/*! \brief Select from double precision variable where boolean is true.
 *
 * \param a double variable to select from
 * \param mask Boolean selector
 * \return  a is selected for true, 0 for false.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
selectByMask(double a, bool mask)
{
    return mask ? a : 0.0;
}

/*! \brief Select from double precision variable where boolean is false.
 *
 * \param a double variable to select from
 * \param mask Boolean selector
 * \return  a is selected for false, 0 for true.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
selectByNotMask(double a, bool mask)
{
    return mask ? 0.0 : a;
}

/*! \brief Blend double selection.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return Select b if sel is true, a otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
blend(double a, double b, bool sel)
{
    return sel ? b : a;
}

/*! \brief Round single precision doubleing point to integer.
 *
 * \param a double
 * \return Integer format, a rounded to nearest integer.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
cvtR2I(double a)
{
    return std::round(a);
};

/*! \brief Truncate single precision doubleing point to integer.
 *
 * \param a double
 * \return Integer format, a truncated to integer.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
cvttR2I(double a)
{
    return std::trunc(a);
};

// We do not have a separate cvtI2R for double, since that would require
// proxy objects. Instead, the float version returns an integer and lets the
// compiler sort out the conversion type.


/*! \brief Convert float to double (mimicks SIMD conversion)
 *
 * \param a float
 * \return a, as double double
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
cvtF2D(float a)
{
    return a;
}

/*! \brief Convert double to float (mimicks SIMD conversion)
 *
 * \param a double
 * \return a, as float
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
cvtD2F(double a)
{
    return a;
}

/************************************************
 *   Integer functions mimicking SIMD versions  *
 ************************************************/

/*! \brief Store contents of integer variable to aligned memory m.
 *
 * \param[out] m Pointer to memory.
 * \param a integer variable to store
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
store(std::int32_t *m, std::int32_t a)
{
    *m = a;
}

/*! \brief Store contents of integer variable to unaligned memory m.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a integer variable to store.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
storeU(std::int32_t *m, std::int32_t a)
{
    *m = a;
}

/*! \brief Bitwise andnot for two scalar integer variables.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
andNot(std::int32_t a, std::int32_t b)
{
    return ~a & b;
}

/*! \brief Return true if any bits are set in the integer variable.
 *
 * This function is used to handle bitmasks, mainly for exclusions in the
 * inner kernels.
 *
 * \param a value
 * \return True if any bit in a is nonzero.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
testBits(std::int32_t a)
{
    return (a != 0);
}

/*! \brief Select from integer variable where boolean is true.
 *
 * \param a Integer variable to select from
 * \param mask Boolean selector
 * \return  a is selected for true, 0 for false.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
selectByMask(std::int32_t a, bool mask)
{
    return mask ? a : 0;
}

/*! \brief Select from integer variable where boolean is false.
 *
 * \param a Integer variable to select from
 * \param mask Boolean selector
 * \return  a is selected for false, 0 for true.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
selectByNotMask(std::int32_t a, bool mask)
{
    return mask ? 0 : a;
}

/*! \brief Blend integer selection.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return Select b if sel is true, a otherwise.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline std::int32_t
blend(std::int32_t a, std::int32_t b, bool sel)
{
    return sel ? b : a;
}

/*! \brief Just return a boolean (mimicks SIMD real-to-int bool conversions)
 *
 * \param a  boolean
 * \return same boolean
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
cvtB2IB(bool a)
{
    return a;
}

/*! \brief Just return a boolean (mimicks SIMD int-to-real bool conversions)
 *
 * \param a  boolean
 * \return same boolean
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline bool
cvtIB2B(bool a)
{
    return a;
}

} // namespace gmx


#endif // GMX_SIMD_SCALAR_FLOAT_H
