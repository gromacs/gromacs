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
/*! \internal \file
 * \brief
 * Implements floating-point comparison routines from testasserts.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testasserts.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <limits>
#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testoptions.h"

namespace gmx
{

namespace test
{

namespace
{

//! Whether to print the message from expected exceptions.
bool g_showExpectedExceptions = false;

//! \cond
GMX_TEST_OPTIONS(ExceptionOptions, options)
{
    options->addOption(BooleanOption("show-error-messages")
                               .store(&g_showExpectedExceptions)
                               .description("Show error messages from expected "
                                            "exceptions"));
}
//! \endcond
} // namespace

namespace internal
{

//! \cond internal
void processExpectedException(const std::exception& ex)
{
    if (g_showExpectedExceptions)
    {
        std::printf("Exception message (from expected exception):\n");
        formatExceptionMessageToFile(stdout, ex);
    }
}
//! \endcond

} // namespace internal

namespace
{

using ::testing::internal::FloatingPoint;

//! \internal \addtogroup module_testutils
//! \{

/*! \name Helper functions for computing floating-point differences
 *
 * These routines are used to initialize FloatingPointDifference.
 * They peek into some internal types from Google Test (gtest-internal.h),
 * and duplicate some other functionality from there, but that is likely
 * a better alternative than just copying all that code here.
 */
//! \{

/*! \brief
 * Computes biased integer representation for a floating-point value.
 *
 * This moves the integer representation from a sign-and-magnitude
 * representation to a biased representation where the 0x8000... represents
 * zero, and the order of the integer values matches the order of the
 * floating-point values.
 */
template<typename FloatType>
typename FloatingPoint<FloatType>::Bits floatingPointToBiasedInteger(const FloatingPoint<FloatType>& value)
{
    if (value.sign_bit())
    {
        return ~value.bits() + 1;
    }
    else
    {
        return value.bits() | FloatingPoint<FloatType>::kSignBitMask;
    }
}

/*! \brief
 * Computes the magnitude of the difference in ULPs between two numbers,
 * treating also values of different sign.
 */
template<typename FloatType>
uint64_t calculateUlpDifference(const FloatingPoint<FloatType>& value1, const FloatingPoint<FloatType>& value2)
{
    typename FloatingPoint<FloatType>::Bits biased1 = floatingPointToBiasedInteger(value1);
    typename FloatingPoint<FloatType>::Bits biased2 = floatingPointToBiasedInteger(value2);
    return biased1 > biased2 ? biased1 - biased2 : biased2 - biased1;
}

/*! \brief
 * Helper to implement the constructors for FloatingPointDifference.
 */
template<typename FloatType>
void initDifference(FloatType raw1, FloatType raw2, double* absoluteDifference, uint64_t* ulpDifference, bool* bSignDifference)
{
    FloatingPoint<FloatType> value1(raw1);
    FloatingPoint<FloatType> value2(raw2);

    if (value1.is_nan() || value2.is_nan())
    {
        *absoluteDifference = std::numeric_limits<double>::quiet_NaN();
        *bSignDifference    = false;
        *ulpDifference      = 0;
        return;
    }
    *absoluteDifference = std::fabs(raw1 - raw2);
    *bSignDifference    = (value1.sign_bit() != value2.sign_bit());
    *ulpDifference      = calculateUlpDifference(value1, value2);
}

/*! \brief
 * Converts a relative tolerance into an ULP difference.
 */
template<typename FloatType>
uint64_t relativeToleranceToUlp(FloatType tolerance)
{
    FloatingPoint<FloatType> m(1.0);
    FloatingPoint<FloatType> t(1.0 + tolerance);
    return calculateUlpDifference<FloatType>(m, t);
}

//! \}
//! \}

} // namespace

/********************************************************************
 * FloatingPointDifference
 */

FloatingPointDifference::FloatingPointDifference(float ref, float value) :
    termMagnitude_(std::abs(ref))
{
    initDifference(ref, value, &absoluteDifference_, &ulpDifference_, &bSignDifference_);
    bDouble_ = false;
}

FloatingPointDifference::FloatingPointDifference(double ref, double value) :
    termMagnitude_(std::abs(ref))
{
    initDifference(ref, value, &absoluteDifference_, &ulpDifference_, &bSignDifference_);
    bDouble_ = true;
}

bool FloatingPointDifference::isNaN() const
{
    return FloatingPoint<double>(absoluteDifference_).is_nan();
}

std::string FloatingPointDifference::toString() const
{
    std::string relDiffStr;

    if (termMagnitude_ > 0)
    {
        // If the reference value is finite we calculate the proper quotient
        relDiffStr = formatString("%.3g", std::abs(absoluteDifference_ / termMagnitude_));
    }
    else if (absoluteDifference_ == 0.0)
    {
        // If the numbers are identical the quotient is strictly NaN here, but
        // there no reason to worry when we have a perfect match.
        relDiffStr = formatString("%.3g", 0.0);
    }
    else
    {
        // If the reference value is zero and numbers are non-identical, relative difference is infinite.
        relDiffStr = formatString("Inf");
    }

    return formatString("%g (%" PRIu64 " %s-prec. ULPs, rel. %s)%s",
                        absoluteDifference_,
                        ulpDifference_,
                        isDouble() ? "double" : "single",
                        relDiffStr.c_str(),
                        bSignDifference_ ? ", signs differ" : "");
}

/********************************************************************
 * FloatingPointTolerance
 */

bool FloatingPointTolerance::isWithin(const FloatingPointDifference& difference) const
{
    if (difference.isNaN())
    {
        return false;
    }

    if (bSignMustMatch_ && difference.signsDiffer())
    {
        return false;
    }

    const double absoluteTolerance =
            difference.isDouble() ? doubleAbsoluteTolerance_ : singleAbsoluteTolerance_;
    if (difference.asAbsolute() < absoluteTolerance)
    {
        return true;
    }

    // By using smaller-than-or-equal below, we allow the test to pass if
    // the numbers are identical, even if the term magnitude is 0, which seems
    // a reasonable thing to do...
    const double relativeTolerance =
            difference.isDouble() ? doubleRelativeTolerance_ : singleRelativeTolerance_;

    if (difference.asAbsolute() <= relativeTolerance * difference.termMagnitude())
    {
        return true;
    }

    const uint64_t ulpTolerance = difference.isDouble() ? doubleUlpTolerance_ : singleUlpTolerance_;
    return ulpTolerance < UINT64_MAX && difference.asUlps() <= ulpTolerance;
}

std::string FloatingPointTolerance::toString(const FloatingPointDifference& difference) const
{
    std::string  result;
    const double absoluteTolerance =
            difference.isDouble() ? doubleAbsoluteTolerance_ : singleAbsoluteTolerance_;
    const double relativeTolerance =
            difference.isDouble() ? doubleRelativeTolerance_ : singleRelativeTolerance_;
    const uint64_t ulpTolerance = difference.isDouble() ? doubleUlpTolerance_ : singleUlpTolerance_;

    if (absoluteTolerance > 0.0)
    {
        result.append(formatString("abs. %g", absoluteTolerance));
    }
    if (relativeTolerance > 0.0)
    {
        if (!result.empty())
        {
            result.append(", ");
        }
        result.append(formatString("rel. %.3g", relativeTolerance));
    }
    if (ulpTolerance < UINT64_MAX)
    {
        if (!result.empty())
        {
            result.append(", ");
        }
        result.append(formatString("%" PRIu64 " ULPs", ulpTolerance));
    }
    if (bSignMustMatch_)
    {
        if (!result.empty())
        {
            result.append(", ");
        }
        result.append("sign must match");
    }
    return result;
}

// Doxygen does not recognize this as the same function as in the header...
//! \cond
FloatingPointTolerance relativeToleranceAsFloatingPoint(double magnitude, double tolerance)
{
    return relativeToleranceAsPrecisionDependentFloatingPoint(magnitude, float(tolerance), tolerance);
}

FloatingPointTolerance relativeToleranceAsPrecisionDependentFloatingPoint(double magnitude,
                                                                          float  singleTolerance,
                                                                          double doubleTolerance)
{
    const float  absoluteSingleTolerance = std::abs(float(magnitude)) * singleTolerance;
    const double absoluteDoubleTolerance = std::abs(magnitude) * doubleTolerance;
    return {
        absoluteSingleTolerance, absoluteDoubleTolerance, singleTolerance, doubleTolerance, UINT64_MAX, UINT64_MAX, false
    };
}
//! \endcond

void checkTestNameLength(std::optional<std::string> testName)
{
    if (!testName.has_value())
    {
        testName = ::testing::UnitTest::GetInstance()->current_test_info()->test_suite_name();
        testName.value() += "_";
        testName.value() += ::testing::UnitTest::GetInstance()->current_test_info()->name();
    }
    const int maxTestLength = 120;
    EXPECT_LE(testName.value().length(), maxTestLength) << formatString(
            "Tests may not use names longer than %d characters\nThis test %s was %zu characters "
            "long",
            maxTestLength,
            testName.value().c_str(),
            testName.value().size());
}

} // namespace test
} // namespace gmx
