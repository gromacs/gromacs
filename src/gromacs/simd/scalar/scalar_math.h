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
#ifndef GMX_SIMD_SCALAR_MATH_H
#define GMX_SIMD_SCALAR_MATH_H

#include <cmath>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/simd/scalar/scalar.h"

/*! \libinternal \file
 *
 * \brief Scalar math functions mimicking GROMACS SIMD math functions
 *
 * These versions make it possible to write functions that are templated with
 * either a SIMD or scalar type. While some of these functions might not appear
 * SIMD-specific, we have placed them here because the only reason to use these
 * instead of generic function is in templated combined SIMD/non-SIMD code.
 * It is important that these functions match the SIMD versions exactly in their
 * arguments and template arguments so that overload resolution works correctly.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

namespace gmx
{

/*****************************************************************************
 *   Single-precision floating point math functions mimicking SIMD versions  *
 *****************************************************************************/


/*! \brief Composes single value with the magnitude of x and the sign of y.
 *
 * \param x Value to set sign for
 * \param y Value used to set sign
 * \return  Magnitude of x, sign of y
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
copysign(float x, float y)
{
    return std::copysign(x, y);
}

// invsqrt(x) is already defined in math/functions.h

/*! \brief Calculate 1/sqrt(x) for two floats.
 *
 * \param x0  First argument, x0 must be positive - no argument checking.
 * \param x1  Second argument, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
invsqrtPair(float x0,    float x1,
            float *out0, float *out1)
{
    *out0 = invsqrt(x0);
    *out1 = invsqrt(x1);
}

/*! \brief Calculate 1/x for float.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
inv(float x)
{
    return 1.0f/x;
}

/*! \brief Calculate 1/sqrt(x) for masked entry of float.
 *
 *  This routine only evaluates 1/sqrt(x) if mask is true.
 *  Illegal values for a masked-out float will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be >0 if masked-in.
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
maskzInvsqrt(float x, bool m)
{
    return m ? invsqrt(x) : 0.0f;
}

/*! \brief Calculate 1/x for masked entry of float.
 *
 *  This routine only evaluates 1/x if mask is true.
 *  Illegal values for a masked-out float will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be nonzero if masked-in.
 *  \param m Mask
 *  \return 1/x. Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
maskzInv(float x, bool m)
{
    return m ? inv(x) : 0.0f;
}

/*! \brief Float sqrt(x). This is the square root.
 *
 * \param x Argument, should be >= 0.
 * \result The square root of x. Undefined if argument is invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline float
sqrt(float x)
{
    return std::sqrt(x);
}

/*! \brief Float log(x). This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
log(float x)
{
    return std::log(x);
}

/*! \brief Float 2^x.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline float
exp2(float x)
{
    return std::exp2(x);
}

/*! \brief Float exp(x).
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline float
exp(float x)
{
    return std::exp(x);
}

/*! \brief Float erf(x).
 *
 * \param x Argument.
 * \result erf(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
erf(float x)
{
    return std::erf(x);
}

/*! \brief Float erfc(x).
 *
 * \param x Argument.
 * \result erfc(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
erfc(float x)
{
    return std::erfc(x);
}


/*! \brief Float sin \& cos.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
sincos(float x, float *sinval, float *cosval)
{
    *sinval = std::sin(x);
    *cosval = std::cos(x);
}

/*! \brief Float sin.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
sin(float x)
{
    return std::sin(x);
}

/*! \brief Float cos.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
cos(float x)
{
    return std::cos(x);
}

/*! \brief Float tan.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
tan(float x)
{
    return std::tan(x);
}

/*! \brief float asin.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
asin(float x)
{
    return std::asin(x);
}

/*! \brief Float acos.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
acos(float x)
{
    return std::acos(x);
}

/*! \brief Float atan.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
atan(float x)
{
    return std::atan(x);
}

/*! \brief Float atan2(y,x).
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
atan2(float y, float x)
{
    return std::atan2(y, x);
}

/*! \brief Calculate the force correction due to PME analytically in float.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on force
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
pmeForceCorrection(float z2)
{
    const float  FN6(-1.7357322914161492954e-8f);
    const float  FN5(1.4703624142580877519e-6f);
    const float  FN4(-0.000053401640219807709149f);
    const float  FN3(0.0010054721316683106153f);
    const float  FN2(-0.019278317264888380590f);
    const float  FN1(0.069670166153766424023f);
    const float  FN0(-0.75225204789749321333f);

    const float  FD4(0.0011193462567257629232f);
    const float  FD3(0.014866955030185295499f);
    const float  FD2(0.11583842382862377919f);
    const float  FD1(0.50736591960530292870f);
    const float  FD0(1.0f);

    float        z4;
    float        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = z2 * z2;

    polyFD0        = fma(FD4, z4, FD2);
    polyFD1        = fma(FD3, z4, FD1);
    polyFD0        = fma(polyFD0, z4, FD0);
    polyFD0        = fma(polyFD1, z2, polyFD0);

    polyFN0        = fma(FN6, z4, FN4);
    polyFN1        = fma(FN5, z4, FN3);
    polyFN0        = fma(polyFN0, z4, FN2);
    polyFN1        = fma(polyFN1, z4, FN1);
    polyFN0        = fma(polyFN0, z4, FN0);
    polyFN0        = fma(polyFN1, z2, polyFN0);

    return polyFN0 / polyFD0;
}

/*! \brief Calculate the potential correction due to PME analytically in float.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on potential.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline float
pmePotentialCorrection(float z2)
{
    const float  VN6(1.9296833005951166339e-8f);
    const float  VN5(-1.4213390571557850962e-6f);
    const float  VN4(0.000041603292906656984871f);
    const float  VN3(-0.00013134036773265025626f);
    const float  VN2(0.038657983986041781264f);
    const float  VN1(0.11285044772717598220f);
    const float  VN0(1.1283802385263030286f);

    const float  VD3(0.0066752224023576045451f);
    const float  VD2(0.078647795836373922256f);
    const float  VD1(0.43336185284710920150f);
    const float  VD0(1.0f);

    float        z4;
    float        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = z2 * z2;

    polyVD1        = fma(VD3, z4, VD1);
    polyVD0        = fma(VD2, z4, VD0);
    polyVD0        = fma(polyVD1, z2, polyVD0);

    polyVN0        = fma(VN6, z4, VN4);
    polyVN1        = fma(VN5, z4, VN3);
    polyVN0        = fma(polyVN0, z4, VN2);
    polyVN1        = fma(polyVN1, z4, VN1);
    polyVN0        = fma(polyVN0, z4, VN0);
    polyVN0        = fma(polyVN1, z2, polyVN0);

    return polyVN0 / polyVD0;
}

/*****************************************************************************
 *   Double-precision floating point math functions mimicking SIMD versions  *
 *****************************************************************************/


/*! \brief Composes double value with the magnitude of x and the sign of y.
 *
 * \param x Value to set sign for
 * \param y Value used to set sign
 * \return  Magnitude of x, sign of y
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
copysign(double x, double y)
{
    return std::copysign(x, y);
}

// invsqrt(x) is already defined in math/functions.h

/*! \brief Calculate 1/sqrt(x) for two doubles.
 *
 * \param x0  First argument, x0 must be positive - no argument checking.
 * \param x1  Second argument, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
invsqrtPair(double x0,    double x1,
            double *out0, double *out1)
{
    *out0 = invsqrt(x0);
    *out1 = invsqrt(x1);
}

/*! \brief Calculate 1/x for double.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
inv(double x)
{
    return 1.0/x;
}

/*! \brief Calculate 1/sqrt(x) for masked entry of double.
 *
 *  This routine only evaluates 1/sqrt(x) if mask is true.
 *  Illegal values for a masked-out double will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be >0 if masked-in.
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzInvsqrt(double x, bool m)
{
    return m ? invsqrt(x) : 0.0;
}

/*! \brief Calculate 1/x for masked entry of double.
 *
 *  This routine only evaluates 1/x if mask is true.
 *  Illegal values for a masked-out double will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be nonzero if masked-in.
 *  \param m Mask
 *  \return 1/x. Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzInv(double x, bool m)
{
    return m ? inv(x) : 0.0;
}

/*! \brief Double sqrt(x). This is the square root.
 *
 * \param x Argument, should be >= 0.
 * \result The square root of x. Undefined if argument is invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline double
sqrt(double x)
{
    return std::sqrt(x);
}

/*! \brief Double log(x). This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
log(double x)
{
    return std::log(x);
}

/*! \brief Double 2^x.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline double
exp2(double x)
{
    return std::exp2(x);
}

/*! \brief Double exp(x).
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
template <MathOptimization opt = MathOptimization::Safe>
static inline double
exp(double x)
{
    return std::exp(x);
}

/*! \brief Double erf(x).
 *
 * \param x Argument.
 * \result erf(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
erf(double x)
{
    return std::erf(x);
}

/*! \brief Double erfc(x).
 *
 * \param x Argument.
 * \result erfc(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
erfc(double x)
{
    return std::erfc(x);
}


/*! \brief Double sin \& cos.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
sincos(double x, double *sinval, double *cosval)
{
    *sinval = std::sin(x);
    *cosval = std::cos(x);
}

/*! \brief Double sin.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
sin(double x)
{
    return std::sin(x);
}

/*! \brief Double cos.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
cos(double x)
{
    return std::cos(x);
}

/*! \brief Double tan.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
tan(double x)
{
    return std::tan(x);
}

/*! \brief Double asin.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
asin(double x)
{
    return std::asin(x);
}

/*! \brief Double acos.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
acos(double x)
{
    return std::acos(x);
}

/*! \brief Double atan.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
atan(double x)
{
    return std::atan(x);
}

/*! \brief Double atan2(y,x).
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
atan2(double y, double x)
{
    return std::atan2(y, x);
}

/*! \brief Calculate the force correction due to PME analytically in double.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on force
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
pmeForceCorrection(double z2)
{
    const double  FN10(-8.0072854618360083154e-14);
    const double  FN9(1.1859116242260148027e-11);
    const double  FN8(-8.1490406329798423616e-10);
    const double  FN7(3.4404793543907847655e-8);
    const double  FN6(-9.9471420832602741006e-7);
    const double  FN5(0.000020740315999115847456);
    const double  FN4(-0.00031991745139313364005);
    const double  FN3(0.0035074449373659008203);
    const double  FN2(-0.031750380176100813405);
    const double  FN1(0.13884101728898463426);
    const double  FN0(-0.75225277815249618847);

    const double  FD5(0.000016009278224355026701);
    const double  FD4(0.00051055686934806966046);
    const double  FD3(0.0081803507497974289008);
    const double  FD2(0.077181146026670287235);
    const double  FD1(0.41543303143712535988);
    const double  FD0(1.0);

    double        z4;
    double        polyFN0, polyFN1, polyFD0, polyFD1;

    z4             = z2 * z2;

    polyFD1        = fma(FD5, z4, FD3);
    polyFD1        = fma(polyFD1, z4, FD1);
    polyFD1        = polyFD1 * z2;
    polyFD0        = fma(FD4, z4, FD2);
    polyFD0        = fma(polyFD0, z4, FD0);
    polyFD0        = polyFD0 + polyFD1;

    polyFD0        = inv(polyFD0);

    polyFN0        = fma(FN10, z4, FN8);
    polyFN0        = fma(polyFN0, z4, FN6);
    polyFN0        = fma(polyFN0, z4, FN4);
    polyFN0        = fma(polyFN0, z4, FN2);
    polyFN0        = fma(polyFN0, z4, FN0);
    polyFN1        = fma(FN9, z4, FN7);
    polyFN1        = fma(polyFN1, z4, FN5);
    polyFN1        = fma(polyFN1, z4, FN3);
    polyFN1        = fma(polyFN1, z4, FN1);
    polyFN0        = fma(polyFN1, z2, polyFN0);

    return polyFN0 * polyFD0;
}

/*! \brief Calculate the potential correction due to PME analytically in double.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on potential.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
pmePotentialCorrection(double z2)
{
    const double  VN9(-9.3723776169321855475e-13);
    const double  VN8(1.2280156762674215741e-10);
    const double  VN7(-7.3562157912251309487e-9);
    const double  VN6(2.6215886208032517509e-7);
    const double  VN5(-4.9532491651265819499e-6);
    const double  VN4(0.00025907400778966060389);
    const double  VN3(0.0010585044856156469792);
    const double  VN2(0.045247661136833092885);
    const double  VN1(0.11643931522926034421);
    const double  VN0(1.1283791671726767970);

    const double  VD5(0.000021784709867336150342);
    const double  VD4(0.00064293662010911388448);
    const double  VD3(0.0096311444822588683504);
    const double  VD2(0.085608012351550627051);
    const double  VD1(0.43652499166614811084);
    const double  VD0(1.0);

    double        z4;
    double        polyVN0, polyVN1, polyVD0, polyVD1;

    z4             = z2 * z2;

    polyVD1        = fma(VD5, z4, VD3);
    polyVD0        = fma(VD4, z4, VD2);
    polyVD1        = fma(polyVD1, z4, VD1);
    polyVD0        = fma(polyVD0, z4, VD0);
    polyVD0        = fma(polyVD1, z2, polyVD0);

    polyVD0        = inv(polyVD0);

    polyVN1        = fma(VN9, z4, VN7);
    polyVN0        = fma(VN8, z4, VN6);
    polyVN1        = fma(polyVN1, z4, VN5);
    polyVN0        = fma(polyVN0, z4, VN4);
    polyVN1        = fma(polyVN1, z4, VN3);
    polyVN0        = fma(polyVN0, z4, VN2);
    polyVN1        = fma(polyVN1, z4, VN1);
    polyVN0        = fma(polyVN0, z4, VN0);
    polyVN0        = fma(polyVN1, z2, polyVN0);

    return polyVN0 * polyVD0;
}


/*****************************************************************************
 *           Floating point math functions mimicking SIMD versions.          *
 *              Double precision data, single precision accuracy.            *
 *****************************************************************************/


/*! \brief Calculate 1/sqrt(x) for double, but with single accuracy.
 *
 *  \param x Argument that must be >0. This routine does not check arguments.
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
invsqrtSingleAccuracy(double x)
{
    return invsqrt(static_cast<float>(x));
}

/*! \brief Calculate 1/sqrt(x) for two doubles, but with single accuracy.
 *
 * \param x0  First argument, x0 must be positive - no argument checking.
 * \param x1  Second argument, x1 must be positive - no argument checking.
 * \param[out] out0  Result 1/sqrt(x0)
 * \param[out] out1  Result 1/sqrt(x1)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
invsqrtPairSingleAccuracy(double x0,    double x1,
                          double *out0, double *out1)
{
    *out0 = invsqrt(static_cast<float>(x0));
    *out1 = invsqrt(static_cast<float>(x1));
}

/*! \brief Calculate 1/x for double, but with single accuracy.
 *
 *  \param x Argument that must be nonzero. This routine does not check arguments.
 *  \return 1/x. Result is undefined if your argument was invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
invSingleAccuracy(double x)
{
    return 1.0f/x;
}

/*! \brief Calculate 1/sqrt(x) for masked entry of double, but with single accuracy.
 *
 *  This routine only evaluates 1/sqrt(x) if mask is true.
 *  Illegal values for a masked-out double will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be >0 if masked-in.
 *  \param m Mask
 *  \return 1/sqrt(x). Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzInvsqrtSingleAccuracy(double x, bool m)
{
    return m ? invsqrtSingleAccuracy(x) : 0.0;
}

/*! \brief Calculate 1/x for masked entry of double, but with single accuracy.
 *
 *  This routine only evaluates 1/x if mask is true.
 *  Illegal values for a masked-out double will not lead to
 *  floating-point exceptions.
 *
 *  \param x Argument that must be nonzero if masked-in.
 *  \param m Mask
 *  \return 1/x. Result is undefined if your argument was invalid or
 *          entry was not masked, and 0.0 for masked-out entries.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
maskzInvSingleAccuracy(double x, bool m)
{
    return m ? invSingleAccuracy(x) : 0.0;
}

/*! \brief Calculate sqrt(x) for double, but with single accuracy.
 *
 *  \param x Argument that must be >=0.
 *  \return sqrt(x).
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
sqrtSingleAccuracy(double x)
{
    return std::sqrt(static_cast<float>(x));
}

/*! \brief Double log(x), but with single accuracy. This is the natural logarithm.
 *
 * \param x Argument, should be >0.
 * \result The natural logarithm of x. Undefined if argument is invalid.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
logSingleAccuracy(double x)
{
    return std::log(static_cast<float>(x));
}

/*! \brief Double 2^x, but with single accuracy.
 *
 * \param x Argument.
 * \result 2^x. Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
exp2SingleAccuracy(double x)
{
    return std::exp2(static_cast<float>(x));
}

/*! \brief Double exp(x), but with single accuracy.
 *
 * \param x Argument.
 * \result exp(x). Undefined if input argument caused overflow.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
expSingleAccuracy(double x)
{
    return std::exp(static_cast<float>(x));
}

/*! \brief Double erf(x), but with single accuracy.
 *
 * \param x Argument.
 * \result erf(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
erfSingleAccuracy(double x)
{
    return std::erf(static_cast<float>(x));
}

/*! \brief Double erfc(x), but with single accuracy.
 *
 * \param x Argument.
 * \result erfc(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
erfcSingleAccuracy(double x)
{
    return std::erfc(static_cast<float>(x));
}


/*! \brief Double sin \& cos, but with single accuracy.
 *
 * \param x The argument to evaluate sin/cos for
 * \param[out] sinval Sin(x)
 * \param[out] cosval Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline void
sincosSingleAccuracy(double x, double *sinval, double *cosval)
{
    // There is no single-precision sincos guaranteed in C++11, so use
    // separate functions and hope the compiler optimizes it for us.
    *sinval = std::sin(static_cast<float>(x));
    *cosval = std::cos(static_cast<float>(x));
}

/*! \brief Double sin, but with single accuracy.
 *
 * \param x The argument to evaluate sin for
 * \result Sin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
sinSingleAccuracy(double x)
{
    return std::sin(static_cast<float>(x));
}

/*! \brief Double cos, but with single accuracy.
 *
 * \param x The argument to evaluate cos for
 * \result Cos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
cosSingleAccuracy(double x)
{
    return std::cos(static_cast<float>(x));
}

/*! \brief Double tan, but with single accuracy.
 *
 * \param x The argument to evaluate tan for
 * \result Tan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
tanSingleAccuracy(double x)
{
    return std::tan(static_cast<float>(x));
}

/*! \brief Double asin, but with single accuracy.
 *
 * \param x The argument to evaluate asin for
 * \result Asin(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
asinSingleAccuracy(double x)
{
    return std::asin(static_cast<float>(x));
}

/*! \brief Double acos, but with single accuracy.
 *
 * \param x The argument to evaluate acos for
 * \result Acos(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
acosSingleAccuracy(double x)
{
    return std::acos(static_cast<float>(x));
}

/*! \brief Double atan, but with single accuracy.
 *
 * \param x The argument to evaluate atan for
 * \result Atan(x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
atanSingleAccuracy(double x)
{
    return std::atan(static_cast<float>(x));
}

/*! \brief Double atan2(y,x), but with single accuracy.
 *
 * \param y Y component of vector, any quartile
 * \param x X component of vector, any quartile
 * \result Atan(y,x)
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
atan2SingleAccuracy(double y, double x)
{
    return std::atan2(static_cast<float>(y), static_cast<float>(x));
}

/*! \brief Force correction due to PME in double, but with single accuracy.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on force
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
pmeForceCorrectionSingleAccuracy(double z2)
{
    const float  FN6(-1.7357322914161492954e-8f);
    const float  FN5(1.4703624142580877519e-6f);
    const float  FN4(-0.000053401640219807709149f);
    const float  FN3(0.0010054721316683106153f);
    const float  FN2(-0.019278317264888380590f);
    const float  FN1(0.069670166153766424023f);
    const float  FN0(-0.75225204789749321333f);

    const float  FD4(0.0011193462567257629232f);
    const float  FD3(0.014866955030185295499f);
    const float  FD2(0.11583842382862377919f);
    const float  FD1(0.50736591960530292870f);
    const float  FD0(1.0f);

    float        z4;
    float        polyFN0, polyFN1, polyFD0, polyFD1;

    float        z2f = z2;

    z4             = z2f * z2f;

    polyFD0        = fma(FD4, z4, FD2);
    polyFD1        = fma(FD3, z4, FD1);
    polyFD0        = fma(polyFD0, z4, FD0);
    polyFD0        = fma(polyFD1, z2f, polyFD0);

    polyFN0        = fma(FN6, z4, FN4);
    polyFN1        = fma(FN5, z4, FN3);
    polyFN0        = fma(polyFN0, z4, FN2);
    polyFN1        = fma(polyFN1, z4, FN1);
    polyFN0        = fma(polyFN0, z4, FN0);
    polyFN0        = fma(polyFN1, z2f, polyFN0);

    return polyFN0 / polyFD0;
}

/*! \brief Potential correction due to PME in double, but with single accuracy.
 *
 * See the SIMD version of this function for details.
 *
 * \param z2 input parameter
 * \returns Correction to use on potential.
 *
 * \note This function might be superficially meaningless, but it helps us to
 *       write templated SIMD/non-SIMD code. For clarity it should not be used
 *       outside such code.
 */
static inline double
pmePotentialCorrectionSingleAccuracy(double z2)
{
    const float  VN6(1.9296833005951166339e-8f);
    const float  VN5(-1.4213390571557850962e-6f);
    const float  VN4(0.000041603292906656984871f);
    const float  VN3(-0.00013134036773265025626f);
    const float  VN2(0.038657983986041781264f);
    const float  VN1(0.11285044772717598220f);
    const float  VN0(1.1283802385263030286f);

    const float  VD3(0.0066752224023576045451f);
    const float  VD2(0.078647795836373922256f);
    const float  VD1(0.43336185284710920150f);
    const float  VD0(1.0f);

    float        z4;
    float        polyVN0, polyVN1, polyVD0, polyVD1;

    float        z2f = z2;

    z4             = z2f * z2f;

    polyVD1        = fma(VD3, z4, VD1);
    polyVD0        = fma(VD2, z4, VD0);
    polyVD0        = fma(polyVD1, z2f, polyVD0);

    polyVN0        = fma(VN6, z4, VN4);
    polyVN1        = fma(VN5, z4, VN3);
    polyVN0        = fma(polyVN0, z4, VN2);
    polyVN1        = fma(polyVN1, z4, VN1);
    polyVN0        = fma(polyVN0, z4, VN0);
    polyVN0        = fma(polyVN1, z2f, polyVN0);

    return polyVN0 / polyVD0;
}



} // namespace gmx


#endif // GMX_SIMD_SCALAR_MATH_H
