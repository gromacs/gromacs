/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Implements internal utility functions for spline tables
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "splineutil.h"

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace internal
{

void throwUnlessDerivativeIsConsistentWithFunction(const std::function<double(double)>& function,
                                                   const std::function<double(double)>& derivative,
                                                   const std::pair<real, real>&         range)
{
    // Since the numerical derivative will evaluate extra points
    // we shrink the interval slightly to avoid calling the function with values
    // outside the range specified.
    double                    h = std::cbrt(GMX_DOUBLE_EPS); // ideal spacing
    std::pair<double, double> newRange(range.first + h, range.second - h);
    const int                 points       = 1000; // arbitrary
    double                    dx           = (newRange.second - newRange.first) / points;
    bool                      isConsistent = true;
    double                    minFail      = newRange.second;
    double                    maxFail      = newRange.first;

    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
    for (double x = newRange.first; x <= newRange.second; x += dx)
    {
        double analyticalDerivative = derivative(x);
        double numericalDerivative  = (function(x + h) - function(x - h)) / (2 * h);
        double thirdDerivative = (derivative(x + h) - 2.0 * derivative(x) + derivative(x - h)) / (h * h);

        // We make two types of errors in numerical calculation of the derivative:
        // - The truncation error: eps * |f| / h
        // - The rounding error: h * h * |f'''| / 6.0
        double expectedErr =
                GMX_DOUBLE_EPS * std::abs(function(x)) / h + h * h * std::abs(thirdDerivative) / 6.0;

        // To avoid triggering false errors because of compiler optimization or numerical issues
        // in the function evalulation we allow an extra factor of 10 in the expected error
        if (std::abs(analyticalDerivative - numericalDerivative) > 10.0 * expectedErr)
        {
            isConsistent = false;
            minFail      = std::min(minFail, x);
            maxFail      = std::max(maxFail, x);
        }
    }

    if (!isConsistent)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Derivative inconsistent with analytical function in range [%f,%f]", minFail, maxFail)));
    }
}


void throwUnlessDerivativeIsConsistentWithFunction(ArrayRef<const double>       function,
                                                   ArrayRef<const double>       derivative,
                                                   double                       inputSpacing,
                                                   const std::pair<real, real>& range)
{
    std::size_t firstIndex   = static_cast<std::size_t>(range.first / inputSpacing);
    std::size_t lastIndex    = static_cast<std::size_t>(range.second / inputSpacing);
    bool        isConsistent = true;
    std::size_t minFail      = lastIndex;
    std::size_t maxFail      = firstIndex;

    // The derivative will access one extra point before/after each point, so reduce interval
    for (std::size_t i = firstIndex + 1; (i + 1) < lastIndex; i++)
    {
        double inputDerivative     = derivative[i];
        double numericalDerivative = (function[i + 1] - function[i - 1]) / (2.0 * inputSpacing);
        double thirdDerivative     = (derivative[i + 1] - 2.0 * derivative[i] + derivative[i - 1])
                                 / (inputSpacing * inputSpacing);

        // We make two types of errors in numerical calculation of the derivative:
        // - The truncation error: eps * |f| / h
        // - The rounding error: h * h * |f'''| / 6.0
        double expectedErr = GMX_DOUBLE_EPS * std::abs(function[i]) / inputSpacing
                             + inputSpacing * inputSpacing * std::abs(thirdDerivative) / 6.0;

        // To avoid triggering false errors because of compiler optimization or numerical issues
        // in the function evalulation we allow an extra factor of 10 in the expected error
        if (std::abs(inputDerivative - numericalDerivative) > 10.0 * expectedErr)
        {
            isConsistent = false;
            minFail      = std::min(minFail, i);
            maxFail      = std::max(maxFail, i);
        }
    }
    if (!isConsistent)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Derivative inconsistent with numerical vector for elements %zu-%zu", minFail + 1, maxFail + 1)));
    }
}


/*! \brief Calculate absolute quotient of function and its second derivative
 *
 * This is a utility function used in the functions to find the smallest quotient
 * in a range.
 *
 * \param[in]    previousPoint Value of function at x-h.
 * \param[in]    thisPoint     Value of function at x.
 * \param[in]    nextPoint     Value of function at x+h.
 * \param[in]    spacing       Value of h.
 *
 * \return The absolute value of the quotient. If either the function or second
 *         derivative is smaller than sqrt(GMX_REAL_MIN), they will be set to
 *         that value.
 */
static double quotientOfFunctionAndSecondDerivative(double previousPoint,
                                                    double thisPoint,
                                                    double nextPoint,
                                                    double spacing)
{
    double lowerLimit = static_cast<double>(std::sqrt(GMX_REAL_MIN));
    double value      = std::max(std::abs(thisPoint), lowerLimit);
    double secondDerivative =
            std::abs((previousPoint - 2.0 * thisPoint + nextPoint) / (spacing * spacing));

    // Make sure we do not divide by zero. This limit is arbitrary,
    // but it doesnt matter since this point will have a very large value,
    // and the whole routine is searching for the smallest value.
    secondDerivative = std::max(secondDerivative, lowerLimit);

    return (value / secondDerivative);
}


real findSmallestQuotientOfFunctionAndSecondDerivative(const std::function<double(double)>& f,
                                                       const std::pair<real, real>&         range)
{
    // Since the numerical second derivative will evaluate extra points
    // we shrink the interval slightly to avoid calling the function with values
    // outside the range specified.
    double                    h = std::pow(GMX_DOUBLE_EPS, 0.25);
    std::pair<double, double> newRange(range.first + h, range.second - h);
    const int                 points      = 500; // arbitrary
    double                    dx          = (newRange.second - newRange.first) / points;
    double                    minQuotient = GMX_REAL_MAX;

    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
    for (double x = newRange.first; x <= newRange.second; x += dx)
    {
        minQuotient = std::min(minQuotient,
                               quotientOfFunctionAndSecondDerivative(f(x - h), f(x), f(x + h), h));
    }

    return static_cast<real>(minQuotient);
}


real findSmallestQuotientOfFunctionAndSecondDerivative(ArrayRef<const double>       function,
                                                       double                       inputSpacing,
                                                       const std::pair<real, real>& range)
{

    std::size_t firstIndex  = static_cast<std::size_t>(range.first / inputSpacing);
    std::size_t lastIndex   = static_cast<std::size_t>(range.second / inputSpacing);
    double      minQuotient = GMX_REAL_MAX;

    for (std::size_t i = firstIndex + 1; (i + 1) < lastIndex; i++)
    {
        minQuotient = std::min(minQuotient,
                               quotientOfFunctionAndSecondDerivative(
                                       function[i - 1], function[i], function[i + 1], inputSpacing));
    }
    return static_cast<real>(minQuotient);
}


/*! \brief Calculate absolute quotient of function and its third derivative
 *
 * This is a utility function used in the functions to find the smallest quotient
 * in a range.
 *
 * \param[in]    previousPreviousPoint Value of function at x-2h.
 * \param[in]    previousPoint         Value of function at x-h.
 * \param[in]    thisPoint             Value of function at x.
 * \param[in]    nextPoint             Value of function at x+h.
 * \param[in]    nextNextPoint         Value of function at x+2h.
 * \param[in]    spacing               Value of h.
 *
 * \return The absolute value of the quotient. If either the function or third
 *         derivative is smaller than sqrt(GMX_REAL_MIN), they will be set to
 *         that value.
 */
static double quotientOfFunctionAndThirdDerivative(double previousPreviousPoint,
                                                   double previousPoint,
                                                   double thisPoint,
                                                   double nextPoint,
                                                   double nextNextPoint,
                                                   double spacing)
{
    double lowerLimit = static_cast<double>(std::sqrt(GMX_REAL_MIN));
    double value      = std::max(std::abs(thisPoint), lowerLimit);
    double thirdDerivative =
            std::abs((nextNextPoint - 2 * nextPoint + 2 * previousPoint - previousPreviousPoint)
                     / (2 * spacing * spacing * spacing));

    // Make sure we do not divide by zero. This limit is arbitrary,
    // but it doesnt matter since this point will have a very large value,
    // and the whole routine is searching for the smallest value.
    thirdDerivative = std::max(thirdDerivative, lowerLimit);

    return (value / thirdDerivative);
}


real findSmallestQuotientOfFunctionAndThirdDerivative(const std::function<double(double)>& f,
                                                      const std::pair<real, real>&         range)
{
    // Since the numerical second derivative will evaluate extra points
    // we shrink the interval slightly to avoid calling the function with values
    // outside the range specified.
    double h = std::pow(GMX_DOUBLE_EPS, 0.2); // optimal spacing for 3rd derivative
    std::pair<double, double> newRange(range.first + 2 * h, range.second - 2 * h);
    const int                 points      = 500; // arbitrary
    double                    dx          = (newRange.second - newRange.first) / points;
    double                    minQuotient = GMX_REAL_MAX;

    // NOLINTNEXTLINE(clang-analyzer-security.FloatLoopCounter)
    for (double x = newRange.first; x <= newRange.second; x += dx)
    {
        minQuotient = std::min(minQuotient,
                               quotientOfFunctionAndThirdDerivative(
                                       f(x - 2 * h), f(x - h), f(x), f(x + h), f(x + 2 * h), h));
    }
    return static_cast<real>(minQuotient);
}


real findSmallestQuotientOfFunctionAndThirdDerivative(ArrayRef<const double>       function,
                                                      double                       inputSpacing,
                                                      const std::pair<real, real>& range)
{

    std::size_t firstIndex  = static_cast<std::size_t>(range.first / inputSpacing);
    std::size_t lastIndex   = static_cast<std::size_t>(range.second / inputSpacing);
    double      minQuotient = GMX_REAL_MAX;

    for (std::size_t i = firstIndex + 2; (i + 2) < lastIndex; i++)
    {
        minQuotient = std::min(
                minQuotient,
                quotientOfFunctionAndThirdDerivative(
                        function[i - 2], function[i - 1], function[i], function[i + 1], function[i + 2], inputSpacing));
    }
    return static_cast<real>(minQuotient);
}


std::vector<double> vectorSecondDerivative(ArrayRef<const double> f, double spacing)
{
    if (f.size() < 5)
    {
        GMX_THROW(APIError("Too few points in vector for 5-point differentiation"));
    }

    std::vector<double> d(f.size());
    std::size_t         i;

    // 5-point formula evaluated for points 0,1
    i    = 0;
    d[i] = (11 * f[i + 4] - 56 * f[i + 3] + 114 * f[i + 2] - 104 * f[i + 1] + 35 * f[i])
           / (12 * spacing * spacing);
    i = 1;
    d[i] = (-f[i + 3] + 4 * f[i + 2] + 6 * f[i + 1] - 20 * f[i] + 11 * f[i - 1]) / (12 * spacing * spacing);

    for (std::size_t i = 2; i < d.size() - 2; i++)
    {
        // 5-point formula evaluated for central point (2)
        d[i] = (-f[i + 2] + 16 * f[i + 1] - 30 * f[i] + 16 * f[i - 1] - f[i - 2]) / (12 * spacing * spacing);
    }

    // 5-point formula evaluated for points 3,4
    i = d.size() - 2;
    d[i] = (11 * f[i + 1] - 20 * f[i] + 6 * f[i - 1] + 4 * f[i - 2] - f[i - 3]) / (12 * spacing * spacing);
    i    = d.size() - 1;
    d[i] = (35 * f[i] - 104 * f[i - 1] + 114 * f[i - 2] - 56 * f[i - 3] + 11 * f[i - 4])
           / (12 * spacing * spacing);

    return d;
}


} // namespace internal

} // namespace gmx
