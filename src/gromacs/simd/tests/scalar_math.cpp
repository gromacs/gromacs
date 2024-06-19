/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/simd/simd.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

#include "data.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \cond internal */
/*! \addtogroup module_simd */
/*! \{ */

// Since these functions are mostly wrappers for similar standard library
// functions, we typically just test 1-2 values to make sure we call the right
// function.

/******************************************
 * Default floating-point precision tests *
 ******************************************/

TEST(SimdScalarMathTest, copysign)
{
    EXPECT_EQ(real(-c1), copysign(real(c1), real(-c2)));
    EXPECT_EQ(real(c2), copysign(real(c2), real(c3)));
}

TEST(SimdScalarMathTest, invsqrtPair)
{
    real x0 = c1;
    real x1 = c2;

    real out0, out1;

    invsqrtPair(x0, x1, &out0, &out1);

    EXPECT_EQ(invsqrt(x0), out0);
    EXPECT_EQ(invsqrt(x1), out1);
}

TEST(SimdScalarMathTest, inv)
{
    real x0 = c0;

    EXPECT_EQ(real(1.0) / x0, inv(x0));
}

TEST(SimdScalarMathTest, maskzInvsqrt)
{
    real x0 = c0;

    EXPECT_EQ(invsqrt(x0), maskzInvsqrt(x0, true));
    EXPECT_EQ(real(0), maskzInvsqrt(x0, false));
}

TEST(SimdScalarMathTest, log)
{
    real x0 = c0;

    EXPECT_EQ(std::log(x0), log(x0));
}

TEST(SimdScalarMathTest, exp2)
{
    real x0 = c0;

    EXPECT_EQ(std::exp2(x0), exp2(x0));
}

TEST(SimdScalarMathTest, exp)
{
    real x0 = c0;

    EXPECT_EQ(std::exp(x0), exp(x0));
}

TEST(SimdScalarMathTest, erf)
{
    real x0 = c0;

    EXPECT_EQ(std::erf(x0), erf(x0));
}

TEST(SimdScalarMathTest, erfc)
{
    real x0 = c0;

    EXPECT_EQ(std::erfc(x0), erfc(x0));
}

TEST(SimdScalarMathTest, sincos)
{
    real x0 = c0;
    real s, c;

    sincos(x0, &s, &c);

    EXPECT_EQ(std::sin(x0), s);
    EXPECT_EQ(std::cos(x0), c);
}

TEST(SimdScalarMathTest, sin)
{
    real x0 = c0;

    EXPECT_EQ(std::sin(x0), sin(x0));
}

TEST(SimdScalarMathTest, cos)
{
    real x0 = c0;

    EXPECT_EQ(std::cos(x0), cos(x0));
}

TEST(SimdScalarMathTest, tan)
{
    real x0 = c0;

    EXPECT_EQ(std::tan(x0), tan(x0));
}


TEST(SimdScalarMathTest, asin)
{
    real x0 = c0;

    EXPECT_EQ(std::asin(x0), asin(x0));
}

TEST(SimdScalarMathTest, acos)
{
    real x0 = c0;

    EXPECT_EQ(std::acos(x0), acos(x0));
}

TEST(SimdScalarMathTest, atan)
{
    real x0 = c0;

    EXPECT_EQ(std::atan(x0), atan(x0));
}

TEST(SimdScalarMathTest, atan2)
{
    real x = c0;
    real y = std::sqrt(c0);


    EXPECT_EQ(std::atan2(y, x), atan2(y, x));
}

TEST(SimdScalarMathTest, pmeForceCorrection)
{
    real z2 = c0;

    // Calculate reference value for z2!=0
    real z   = std::sqrt(z2);
    real ref = 2.0 * std::exp(-z2) / (std::sqrt(M_PI) * z2) - std::erf(z) / (z2 * z);

    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#if GMX_DOUBLE
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-10));
#else
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-6));
#endif

    EXPECT_REAL_EQ_TOL(ref, pmeForceCorrection(z2), tolerance);
}

TEST(SimdScalarMathTest, pmePotentialCorrection)
{
    real z2 = c0;

    // Calculate reference value for z2!=0
    real z   = std::sqrt(z2);
    real ref = std::erf(z) / z;

    // Pme correction only needs to be ~1e-6 accuracy single, 1e-10 double
#if GMX_DOUBLE
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-10));
#else
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-6));
#endif

    EXPECT_REAL_EQ_TOL(ref, pmePotentialCorrection(z2), tolerance);
}

/*******************************************
 * Double precision, single accuracy tests *
 *******************************************/

TEST(SimdScalarMathTest, invsqrtPairSingleAccuracy)
{
    double x0 = c1;
    double x1 = c2;

    double out0, out1;

    invsqrtPairSingleAccuracy(x0, x1, &out0, &out1);

    EXPECT_EQ(invsqrt(static_cast<float>(x0)), static_cast<float>(out0));
    EXPECT_EQ(invsqrt(static_cast<float>(x1)), static_cast<float>(out1));
}

TEST(SimdScalarMathTest, invSingleAccuracy)
{
    double x0 = c1;

    EXPECT_EQ(1.0F / static_cast<float>(x0), static_cast<float>(invSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, maskzInvsqrtSingleAccuracy)
{
    double x0 = c1;

    EXPECT_EQ(invsqrt(static_cast<float>(x0)), static_cast<float>(maskzInvsqrtSingleAccuracy(x0, true)));
    EXPECT_EQ(0.0F, static_cast<float>(maskzInvsqrtSingleAccuracy(x0, false)));
}

TEST(SimdScalarMathTest, logSingleAccuracy)
{
    double x0 = c1;

    EXPECT_EQ(std::log(static_cast<float>(x0)), static_cast<float>(logSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, exp2SingleAccuracy)
{
    double x0 = c1;

    EXPECT_EQ(std::exp2(static_cast<float>(x0)), static_cast<float>(exp2SingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, expSingleAccuracy)
{
    double x0 = c1;

    EXPECT_EQ(std::exp(static_cast<float>(x0)), static_cast<float>(expSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, erfSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::erf(static_cast<float>(x0)), static_cast<float>(erfSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, erfcSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::erfc(static_cast<float>(x0)), static_cast<float>(erfcSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, sincosSingleAccuracy)
{
    double x0 = c0;
    double s, c;

    sincosSingleAccuracy(x0, &s, &c);

    EXPECT_EQ(std::sin(static_cast<float>(x0)), static_cast<float>(s));
    EXPECT_EQ(std::cos(static_cast<float>(x0)), static_cast<float>(c));
}

TEST(SimdScalarMathTest, sinSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::sin(static_cast<float>(x0)), static_cast<float>(sinSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, cosSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::cos(static_cast<float>(x0)), static_cast<float>(cosSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, tanSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::tan(static_cast<float>(x0)), static_cast<float>(tanSingleAccuracy(x0)));
}


TEST(SimdScalarMathTest, asinSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::asin(static_cast<float>(x0)), static_cast<float>(asinSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, acosSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::acos(static_cast<float>(x0)), static_cast<float>(acosSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, atanSingleAccuracy)
{
    double x0 = c0;

    EXPECT_EQ(std::atan(static_cast<float>(x0)), static_cast<float>(atanSingleAccuracy(x0)));
}

TEST(SimdScalarMathTest, atan2SingleAccuracy)
{
    double x = c0;
    double y = std::sqrt(c0);


    EXPECT_EQ(std::atan2(static_cast<float>(y), static_cast<float>(x)),
              static_cast<float>(atan2SingleAccuracy(y, x)));
}

TEST(SimdScalarMathTest, pmeForceCorrectionSingleAccuracy)
{
    double z2 = c0;

    // Calculate reference value for z2!=0 in single precision
    float z   = std::sqrt(static_cast<float>(z2));
    float ref = 2.0 * std::exp(static_cast<float>(-z2)) / (std::sqrt(static_cast<float>(M_PI)) * z2)
                - std::erf(z) / (z2 * z);

    // Pme correction only needs to be ~1e-6 accuracy single
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-6));

    EXPECT_REAL_EQ_TOL(ref, static_cast<float>(pmeForceCorrectionSingleAccuracy(z2)), tolerance);
}

TEST(SimdScalarMathTest, pmePotentialCorrectionSingleAccuracy)
{
    double z2 = c0;

    // Calculate reference value for z2!=0 in single precision
    float z   = std::sqrt(static_cast<float>(z2));
    float ref = std::erf(z) / z;

    // Pme correction only needs to be ~1e-6 accuracy single
    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 5e-6));

    EXPECT_REAL_EQ_TOL(ref, static_cast<float>(pmePotentialCorrectionSingleAccuracy(z2)), tolerance);
}

/*! \} */
/*! \endcond internal */

} // namespace
} // namespace test
} // namespace gmx
