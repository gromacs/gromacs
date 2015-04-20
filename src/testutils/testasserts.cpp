/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements floating-point comparison routines from testasserts.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testasserts.h"

#include <cmath>

#include <limits>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace test
{

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
template <typename FloatType>
typename FloatingPoint<FloatType>::Bits
floatingPointToBiasedInteger(const FloatingPoint<FloatType> &value)
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
template <typename FloatType>
gmx_uint64_t calculateUlpDifference(const FloatingPoint<FloatType> &value1,
                                    const FloatingPoint<FloatType> &value2)
{
    typename FloatingPoint<FloatType>::Bits biased1
        = floatingPointToBiasedInteger(value1);
    typename FloatingPoint<FloatType>::Bits biased2
        = floatingPointToBiasedInteger(value2);
    return biased1 > biased2 ? biased1 - biased2 : biased2 - biased1;
}

/*! \brief
 * Helper to implement the constructors for FloatingPointDifference.
 */
template <typename FloatType>
void initDifference(FloatType raw1, FloatType raw2, double *absoluteDifference,
                    gmx_uint64_t *ulpDifference, bool *bSignDifference)
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
template <typename FloatType>
gmx_uint64_t relativeToleranceToUlp(FloatType tolerance)
{
    FloatingPoint<FloatType> m(1.0);
    FloatingPoint<FloatType> t(1.0 + tolerance);
    return calculateUlpDifference<FloatType>(m, t);
}

//! \}
//! \}

}       // namespace

/********************************************************************
 * FloatingPointDifference
 */

FloatingPointDifference::FloatingPointDifference(float value1, float value2)
{
    initDifference(value1, value2,
                   &absoluteDifference_, &ulpDifference_, &bSignDifference_);
    bDouble_ = false;
}

FloatingPointDifference::FloatingPointDifference(double value1, double value2)
{
    initDifference(value1, value2,
                   &absoluteDifference_, &ulpDifference_, &bSignDifference_);
    bDouble_ = true;
}

bool FloatingPointDifference::isNaN() const
{
    return FloatingPoint<double>(absoluteDifference_).is_nan();
}

std::string FloatingPointDifference::toString() const
{
    const double eps = isDouble() ? GMX_DOUBLE_EPS : GMX_FLOAT_EPS;
    return formatString("%g (%" GMX_PRIu64 " %s-prec. ULPs, rel. %.3g)%s",
                        absoluteDifference_, ulpDifference_,
                        isDouble() ? "double" : "single",
                        ulpDifference_ * eps,
                        bSignDifference_ ? ", signs differ" : "");
}

/********************************************************************
 * FloatingPointTolerance
 */

bool FloatingPointTolerance::isWithin(
        const FloatingPointDifference &difference) const
{
    if (difference.isNaN())
    {
        return false;
    }

    if (bSignMustMatch_ && difference.signsDiffer())
    {
        return false;
    }

    const double absoluteTolerance
        = difference.isDouble() ? doubleAbsoluteTolerance_ : singleAbsoluteTolerance_;
    if (difference.asAbsolute() < absoluteTolerance)
    {
        return true;
    }

    const gmx_uint64_t ulpTolerance
        = difference.isDouble() ? doubleUlpTolerance_ : singleUlpTolerance_;
    if (ulpTolerance < GMX_UINT64_MAX && difference.asUlps() <= ulpTolerance)
    {
        return true;
    }
    return false;
}

std::string FloatingPointTolerance::toString(const FloatingPointDifference &difference) const
{
    std::string        result;
    const double       absoluteTolerance
        = difference.isDouble() ? doubleAbsoluteTolerance_ : singleAbsoluteTolerance_;
    const gmx_uint64_t ulpTolerance
        = difference.isDouble() ? doubleUlpTolerance_ : singleUlpTolerance_;
    const double       eps
        = difference.isDouble() ? GMX_DOUBLE_EPS : GMX_FLOAT_EPS;

    if (absoluteTolerance > 0.0)
    {
        result.append(formatString("abs. %g", absoluteTolerance));
    }
    if (ulpTolerance < GMX_UINT64_MAX)
    {
        if (!result.empty())
        {
            result.append(", ");
        }
        result.append(formatString("%" GMX_PRIu64 " ULPs (rel. %.3g)", ulpTolerance, ulpTolerance * eps));
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
FloatingPointTolerance
relativeToleranceAsFloatingPoint(double magnitude, double tolerance)
{
    const double absoluteTolerance = magnitude * tolerance;
    return FloatingPointTolerance(absoluteTolerance, absoluteTolerance,
                                  relativeToleranceToUlp<float>(tolerance),
                                  relativeToleranceToUlp<double>(tolerance),
                                  false);
}
//! \endcond

} // namespace test
} // namespace gmx
