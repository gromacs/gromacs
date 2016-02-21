/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements classes for quadratic spline table functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "quadraticsplinetable.h"

#include <cstdint>

#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

#include "gromacs/math/utilities.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace
{


/******************************************************************************
 *                                                                            *
 *                Routines working on analytical functions                    *
 *                                                                            *
 ******************************************************************************/

/*! \brief Check that analytical derivative is the derivative of analytical function.
 *
 *  This routine evaluates the numerical derivative of the function for
 *  a few (1000) points in the interval and checks that the relative difference
 *  between numerical and analytical derivative is within 1e-5.
 *
 *  The main point of this routine is to make sure the user has not made a
 *  mistake or sign error when defining the functions.
 *
 *  \param function   Analytical function to differentiate
 *  \param derivative Analytical derivative to compare with
 *  \param range      Range to test
 *
 *  \return true if the derivative is consistent with the function.
 *
 *  \note The function/derivative are always double-valued to avoid accuracy loss.
 */
bool
derivativeIsConsistentWithFunction(const std::function<double(double)>  &function,
                                   const std::function<double(double)>  &derivative,
                                   const std::pair<real, real>           &range)
{
    // Since the numerical derivative will evaluate extra points
    // we shrink the interval slightly to avoid calling the function with values
    // outside the range specified.
    double                     h            = std::cbrt(GMX_DOUBLE_EPS); // ideal spacing
    std::pair<double, double>  newRange(range.first + h, range.second - h);
    const int                  points       = 1000;                      // arbitrary
    double                     dx           = (newRange.second - newRange.first) / points;
    bool                       isConsistent = true;

    for (double x = newRange.first; x <= newRange.second && isConsistent; x += dx)
    {
        double analyticalDerivative = derivative(x);
        double numericalDerivative  = (function(x+h)-function(x-h))/(2*h);
        double thirdDerivative      = (derivative(x+h)-2.0*derivative(x)+derivative(x-h))/(h*h);

        // We make two types of errors in numerical calculation of the derivative:
        // - The truncation error: eps * |f| / h
        // - The rounding error: h * h * |f'''| / 6.0
        double expectedErr = GMX_DOUBLE_EPS*std::abs(function(x))/h + h*h*std::abs(thirdDerivative)/6.0;

        // To avoid triggering false errors because of compiler optimization or numerical issues
        // in the function evalulation we allow an extra factor of 10 in the expected error
        isConsistent = isConsistent && (std::abs(analyticalDerivative-numericalDerivative) < 10.0*expectedErr);
    }
    return isConsistent;
}


/*! \brief Find smallest quotient between function and 2nd derivative (analytical)
 *
 *  Used to calculate table spacing. This function divides the function value
 *  by the second derivative (or a very small number when that is zero), and
 *  returns the smallest such quotient found in the range.
 *
 *  Our quadratic tables corresponds to linear interpolation of the derivative,
 *  which means the derivative will typically have larger error than the value
 *  when interpolating. The spacing required to reach a particular relative
 *  tolerance in the derivative depends on the quotient between the first
 *  derivative and the third derivative of the function itself.
 *
 *  To avoid excessive fragile high-order numerical differentiation, you should
 *  call this routine with the analytical derivative, and the quotient between
 *  "function and second derivative" will then correspond to the quotient
 *  bewteen "derivative and third derivative" of the actual function.
 *
 *  When the second derivative is very close to zero, that will correspond to
 *  a region where we are not sensitive to spacing. That point will not result in
 *  the smallest quotient anyway, so we just make sure to avoid division by zero.
 *
 *  Since all functions that can be tabulated efficiently are reasonably smooth,
 *  we simply check 1,000 points in the interval rather than bother about
 *  implementing any complicated global optimization scheme.
 *
 *  Numerical differentiation can be very fragile and at best this routine
 *  will achieve a relative accuracy of 1e-8, but in practice it can be 1e-6.
 *  However, this is perfectly sufficient to set the required table spacing.
 *
 *  \param f          Analytical function
 *  \param range      Interval
 *
 *  \return Smallest quotient found in range.
 *
 *  \note The function is always double-valued to avoid accuracy loss.
 */
real
findSmallestQuotientOfFunctionAndSecondDerivative(const std::function<double(double)>  &f,
                                                  const std::pair<real, real>           &range)
{
    // Since the numerical second derivative will evaluate extra points
    // we shrink the interval slightly to avoid calling the function with values
    // outside the range specified.
    double                     h           = std::pow( GMX_DOUBLE_EPS, 0.25 );
    std::pair<double, double>  newRange(range.first + h, range.second - h);
    const int                  points      = 1000; // arbitrary
    double                     dx          = (newRange.second - newRange.first) / points;
    double                     minQuotient = GMX_REAL_MAX;

    for (double x = newRange.first; x <= newRange.second; x += dx)
    {
        double value             = std::abs( f(x) );
        double secondDerivative  = std::abs( f(x+h) - 2.0 * f(x) + f(x-h) )/( h * h );

        // Make sure we do not divide by zero. This limit is arbitrary,
        // but it doesnt matter since this point will have a very large value,
        // and the whole routine is searching for the smallest value.
        secondDerivative = std::max(secondDerivative, static_cast<double>(std::sqrt(GMX_REAL_MIN)));

        minQuotient  = std::min(minQuotient, value / secondDerivative);
    }
    return static_cast<real>(minQuotient);
}




/******************************************************************************
 *                                                                            *
 *                 Routines working on numerical vectors                      *
 *                                                                            *
 ******************************************************************************/




/*! \brief Check that numerical derivative is the derivative of numerical function.
 *
 *  This routine differentiates a vector of numerical values and checks
 *  that the relative difference to a provided vector of numerical derivatives
 *  is small. Since we do not know how the numerical derivatives were obtained,
 *  we are quite lenient and only require the relative difference is within 0.001.
 *  The main point of this routine is to make sure the user has not made a
 *  mistake or sign error when defining the functions.
 *
 *  To avoid problems if the vectors change from zero to finite values at the
 *  start/end of the interval, we only check inside the range requested.
 *
 *  \param function     Numerical function value vector to differentiate
 *  \param derivative   Numerical derivative vector to compare with
 *  \param inputSpacing Distance between input points
 *  \param range        Range to test
 *
 *  \return true if the derivative is consistent with the function.
 *
 *  \note The function/derivative vectors and spacing are always double-valued
 *        to avoid accuracy loss.
 */
bool
derivativeIsConsistentWithFunction(const std::vector<double>    &function,
                                   const std::vector<double>    &derivative,
                                   double                        inputSpacing,
                                   const std::pair<real, real>   &range)
{
    bool            isConsistent = true;

    std::size_t     firstIndex   = range.first / inputSpacing;
    std::size_t     lastIndex    = range.second / inputSpacing;

    // The derivative will access one extra point before/after each point, so reduce interval
    for (std::size_t i = firstIndex + 1; (i + 1) < lastIndex; i++)
    {
        double inputDerivative    = derivative[i];
        double computedDerivative = (function[i+1] - function[i-1]) / (2.0 * inputSpacing);
        double thirdDerivative    = (derivative[i+1] - 2.0*derivative[i] + derivative[i-1])/(inputSpacing * inputSpacing);

        // We make two types of errors in numerical calculation of the derivative:
        // - The truncation error: eps * |f| / h
        // - The rounding error: h * h * |f'''| / 6.0
        double expectedErr = GMX_DOUBLE_EPS*std::abs(function[i])/inputSpacing + inputSpacing*inputSpacing*std::abs(thirdDerivative)/6.0;

        // To avoid triggering false errors because of compiler optimization or numerical issues
        // in the function evalulation we allow an extra factor of 10 in the expected error
        isConsistent = isConsistent && (std::abs(inputDerivative-computedDerivative) < 10.0*expectedErr);
    }
    return isConsistent;
}




/*! \brief Find smallest quotient between function and 2nd derivative (vectors)
 *
 *  Used to calculate table spacing. This function divides the function value
 *  by the second derivative (or a very small number when that is zero), and
 *  returns the smallest such quotient found in the range.
 *
 *  Our quadratic tables corresponds to linear interpolation of the derivative,
 *  which means the derivative will typically have larger error than the value
 *  when interpolating. The spacing required to reach a particular relative
 *  tolerance in the derivative depends on the quotient between the first
 *  derivative and the third derivative of the function itself.
 *
 *  To avoid excessive fragile high-order numerical differentiation, you should
 *  call this routine with the analytical derivative, and the quotient between
 *  "function and second derivative" will then correspond to the quotient
 *  bewteen "derivative and third derivative" of the actual function.
 *
 *  When the second derivative is very close to zero, that will correspond to
 *  a region where we are not sensitive to spacing. That point will not result in
 *  the smallest quotient anyway, so we just make sure to avoid division by zero.
 *
 *  \param function     Vector with function values
 *  \param inputSpacing Spacing between function values
 *  \param range        Interval to check
 *
 *  \return Smallest quotient found in range.
 *
 *  \note The function vector and input spacing are always double-valued to
 *        avoid accuracy loss.
 */
real
findSmallestQuotientOfFunctionAndSecondDerivative(const std::vector<double>    &function,
                                                  double                        inputSpacing,
                                                  const std::pair<real, real>   &range)
{

    std::size_t  firstIndex  = range.first  / inputSpacing;
    std::size_t  lastIndex   = range.second / inputSpacing;
    double       minQuotient = GMX_REAL_MAX;

    for (std::size_t i = firstIndex + 1; (i + 1) < lastIndex; i++)
    {
        double value             = std::abs(function[i]);
        double secondDerivative  = std::abs((function[i+1] - 2.0 * function[i] + function[i-1]) / (inputSpacing * inputSpacing));

        // Make sure we do not divide by zero. This limit is arbitrary,
        // but it doesnt matter since this point will have a very large value,
        // and the whole routine is searching for the smallest value.
        secondDerivative = std::max(secondDerivative, static_cast<double>(std::sqrt(GMX_REAL_MIN)));

        minQuotient  = std::min(minQuotient, value / secondDerivative);
    }
    return static_cast<real>(minQuotient);
}


/*! \brief Calculate derivative of vector and return derivative vector of same length
 *
 *  5-point approximations are used for accuracy, with the endpoints using
 *  non-center interpolation.
 *
 *  \param f       Vector (function) for which to calculate derivative
 *  \param spacing Spacing of input data.
 *
 * \note This function always uses double precision arguments since it is meant
 *       to be used on raw user input data for tables, where we want to avoid
 *       accuracy loss (since differentiation can be numerically fragile).
 */
std::vector<double>
vectorDerivative(const std::vector<double> &f, double spacing)
{
    GMX_ASSERT(f.size() >= 5, "Too few points in vector for 5-point differentiation");

    std::vector<double> d(f.size());
    std::size_t         i;

    // 5-point formula evaluated for points 0,1
    i    = 0;
    d[i] = (-3 * f[i+4] + 16 * f[i+3] - 36 * f[i+2] + 48 * f[i+1] - 25 * f[i]) / ( 12 * spacing);
    i    = 1;
    d[i] = (f[i+3] - 6 * f[i+2] + 18 * f[i+1] - 10 * f[i] - 3 * f[i-1]) / ( 12 * spacing);

    for (std::size_t i = 2; i < d.size() - 2; i++)
    {
        // 5-point formula evaluated for central point (2)
        d[i] = (-f[i+2] + 8 * f[i+1] - 8 * f[i-1] + f[i-2]) / (12 * spacing);
    }

    // 5-point formula evaluated for points 3,4
    i    = d.size() - 2;
    d[i] = (3 * f[i+1] + 10 * f[i] - 18 * f[i-1] + 6 * f[i-2] - f[i-3]) / ( 12 * spacing);
    i    = d.size() - 1;
    d[i] = (25 * f[i] - 48 * f[i-1] + 36 * f[i-2] - 16 * f[i-3] + 3 * f[i-4]) / ( 12 * spacing);

    return d;
}


/*! \brief Calculate second derivative of vector and return vector of same length
 *
 *  5-point approximations are used for accuracy, with the endpoints using
 *  non-center interpolation.
 *
 *  \param f       Vector (function) for which to calculate second derivative
 *  \param spacing Spacing of input data.
 *
 * \note This function always uses double precision arguments since it is meant
 *       to be used on raw user input data for tables, where we want to avoid
 *       accuracy loss (since differentiation can be numerically fragile).
 */
std::vector<double>
vectorSecondDerivative(const std::vector<double> &f, double spacing)
{
    GMX_ASSERT(f.size() >= 5, "Too few points in vector for 5-point differentiation");

    std::vector<double> d(f.size());
    std::size_t         i;

    // 5-point formula evaluated for points 0,1
    i    = 0;
    d[i] = (11 * f[i+4] - 56 * f[i+3] + 114 * f[i+2] - 104 * f[i+1] + 35 * f[i]) / ( 12 * spacing * spacing);
    i    = 1;
    d[i] = (-f[i+3] + 4 * f[i+2] + 6 * f[i+1] - 20 * f[i] + 11 * f[i-1]) / ( 12 * spacing * spacing);

    for (std::size_t i = 2; i < d.size() - 2; i++)
    {
        // 5-point formula evaluated for central point (2)
        d[i] = (-f[i+2] + 16 * f[i+1] - 30 * f[i] + 16 * f[i-1] - f[i-2]) / (12 * spacing * spacing);
    }

    // 5-point formula evaluated for points 3,4
    i    = d.size() - 2;
    d[i] = (11 * f[i+1] - 20 * f[i] + 6 * f[i-1] + 4 * f[i-2] - f[i-3]) / ( 12 * spacing * spacing);
    i    = d.size() - 1;
    d[i] = (35 * f[i] - 104 * f[i-1] + 114 * f[i-2] - 56 * f[i-3] + 11 * f[i-4]) / ( 12 * spacing * spacing);

    return d;
}


/*! \brief Return true if all elements in vectors are finite, otherwise false
 *
 *  \param v1   First vector to test
 *  \param v2   Second vector to test
 *
 *  \return true if all elements in both vectors are smaller than
 *          0.0001 * GMX_REAL_MAX, otherwise false.
 */
bool
checkVectorsAreFinite(std::vector<real>  &v1,
                      std::vector<real>  &v2)
{
    real maxValue = 0.0001 * GMX_REAL_MAX;

    for (real e : v1)
    {
        if (!(std::abs(e) < maxValue) )
        {
            return false;
        }
    }

    for (real e : v2)
    {
        if (!(std::abs(e) < maxValue) )
        {
            return false;
        }
    }
    return true; // all values were finite if we got here
}


/*! \brief Create merged DDFZ vector from function & derivative data
 *
 *  \param functionData     Function values
 *  \param derivativeData   Derivative values. We have already subtracted the
 *                          small third derivative component when calling this
 *                          function, but in practice it is just an arbitrary
 *                          vector here.
 *  \return ddfzData        Vector four times longer, filled with
 *                          the derivative, the difference to the next derivative
 *                          point, the function value, and zero.
 */
std::vector<real, AlignedAllocator<real> >
makeDdfzVector(const std::vector<real>    &functionData,
               const std::vector<real>    &derivativeData)
{
    std::vector<real, AlignedAllocator<real> > ddfzData;

    GMX_ASSERT(functionData.size() == derivativeData.size(), "Mismatching vector lengths");

    for (std::size_t i = 0; i < functionData.size(); i++)
    {
        ddfzData.push_back( derivativeData[i] );
        ddfzData.push_back( derivativeData[i+1] - derivativeData[i] );
        ddfzData.push_back( functionData[i] );
        ddfzData.push_back( 0.0 );
    }
    return ddfzData;
}

}


const real
QuadraticSplineTable::defaultTolerance = 10.0 * GMX_FLOAT_EPS;

QuadraticSplineTable::QuadraticSplineTable(const std::function<double(double)>  &function,
                                           const std::function<double(double)>  &derivative,
                                           const std::pair<real, real>           &range,
                                           real                                  tolerance)
    : range_(range)
{
    // Sanity check on input values
    if (range.first < 0.0 || (range.second-range.first) < 0.001)
    {
        GMX_THROW(APIError("Range to tabulate cannot include negative values and must span at least 0.001"));
    }

    if (tolerance < GMX_REAL_EPS)
    {
        GMX_THROW(ToleranceError("Table tolerance cannot be smaller than GMX_REAL_EPS"));
    }

    if (!derivativeIsConsistentWithFunction(function, derivative, range_))
    {
        GMX_THROW(APIError("Derivative inconsistent with function, or function/derivative values too large"));
    }

    // Calculate the required table spacing h. The error we make with linear interpolation
    // of the derivative will be described by the third-derivative correction term.
    // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
    // where f'/f''' is the first and third derivative of the function, respectively.

    double minQuotient = findSmallestQuotientOfFunctionAndSecondDerivative(derivative, range_);
    double spacing     = std::sqrt(12.0 * tolerance * minQuotient);

    halfSpacing_ = 0.5 * spacing;
    tableScale_  = 1.0 / spacing;

    if (range_.second * tableScale_ > 1e6)
    {
        GMX_THROW(ToleranceError("Over a million points would be required for table; decrease range or increase tolerance"));
    }

    // Fill the actual data. Largest index we will access is range.second * tableScale_ + 1,
    // since we load the next derivative point during interpolation.
    std::size_t  minIndex = range_.first * tableScale_;
    std::size_t  maxIndex = range_.second * tableScale_ + 1;

    for (std::size_t i = 0; i < minIndex; i++)
    {
        functionData_.push_back(0.0);
        derivativeData_.push_back(0.0);
    }

    for (std::size_t i = minIndex; i <= maxIndex; i++)
    {
        double x               = i * spacing;
        functionData_.push_back( function(x) );

        // Calculate third derivative term (2nd derivative of the derivative)
        // Make sure we stay in range. In practice this means we use one-sided
        // interpolation at the interval endpoints (indentical to an offset for 3-point formula)
        const double h                    = std::pow( GMX_DOUBLE_EPS, 0.25 );
        double       y                    = std::min( std::max(x, range_.first + h), range_.second - h);
        double       thirdDerivativeValue = ( derivative(y+h) - 2.0 * derivative(y) + derivative(y-h) ) / ( h * h );

        derivativeData_.push_back( derivative(x) - spacing * spacing * thirdDerivativeValue / 12.0 );
    }

    if (checkVectorsAreFinite(functionData_, derivativeData_) == false)
    {
        GMX_THROW(APIError("Tabluated function has too large value or derivative in requested range"));
    }

    derivativeDeltaFunctionZeroData_ = makeDdfzVector(functionData_, derivativeData_);
}


QuadraticSplineTable::QuadraticSplineTable(const std::vector<double>     &function,
                                           const std::vector<double>     &derivative,
                                           double                         inputSpacing,
                                           const std::pair<real, real>    &range,
                                           real                           tolerance)
    : range_(range)
{
    // Sanity check on input values
    if (range.first < 0.0 || (range.second-range.first) < 0.001)
    {
        GMX_THROW(APIError("Range to tabulate cannot include negative values and must span at least 0.001"));
    }

    // We do not yet know what the margin is, but we need to test that we at least cover
    // the requested range before starting to calculate derivatives
    if (function.size() < range_.second / inputSpacing + 1)
    {
        GMX_THROW(APIError("Table input vectors must cover requested range, and a margin beyond the upper endpoint"));
    }

    if (function.size() != derivative.size())
    {
        GMX_THROW(APIError("Function and derivative vectors have different lengths"));
    }
    //XXX

    if (!derivativeIsConsistentWithFunction(function, derivative, inputSpacing, range_))
    {
        GMX_THROW(APIError("Derivative inconsistent with function, or function/derivative values too large"));
    }

    if (tolerance < GMX_REAL_EPS)
    {
        GMX_THROW(ToleranceError("Table tolerance cannot be smaller than GMX_REAL_EPS"));
    }

    // Calculate the required table spacing h. The error we make with linear interpolation
    // of the derivative will be described by the third-derivative correction term.
    // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
    // where f'/f''' is the first and third derivative of the function, respectively.

    double minQuotient = findSmallestQuotientOfFunctionAndSecondDerivative(derivative, inputSpacing, range_);
    double spacing     = std::sqrt(12.0 * tolerance * minQuotient);

    // During interpolation we load the next table point, so the largest value we need to
    // evaluate during construction will be (range.second+spacing)/inputSpacing.

    // Now we know exactly how large the margin we need is.
    // The largest argument for which to evaluate the table will be epsilon below range.second,
    // but since we access the next point that one is equivalent to adding 'spacing'
    if (function.size() < (range_.second + spacing) / inputSpacing)
    {
        GMX_THROW(APIError("Table input vectors must cover a margin beyond the upper endpoint"));
    }

    if (spacing < inputSpacing)
    {
        GMX_THROW(ToleranceError("Input vector spacing cannot achieve tolerance requested"));
    }

    halfSpacing_ = 0.5 * spacing;
    tableScale_  = 1.0 / spacing;

    if (range_.second * tableScale_ > 1e6)
    {
        GMX_THROW(ToleranceError("Requested tolerance would require over a million points in table"));
    }

    std::vector<double>  thirdDerivative(vectorSecondDerivative(derivative, inputSpacing));

    // Fill the actual data. Largest index we will access is range.second * tableScale_,
    // since we load the next derivative point during interpolation.
    std::size_t  minIndex = range_.first * tableScale_;
    std::size_t  maxIndex = range_.second * tableScale_ + 1;

    for (std::size_t i = 0; i < minIndex; i++)
    {
        functionData_.push_back(0.0);
        derivativeData_.push_back(0.0);
    }

    for (std::size_t i = minIndex; i < maxIndex; i++)
    {
        double x          = i * spacing; // refers to output table

        // Step 1: Interpolate the function value at x from input table.
        double inputXTab  = x / inputSpacing;
        int    inputIndex = inputXTab;
        double inputEps   = inputXTab - inputIndex;

        // Linear interpolation of input derivative and third derivative
        double thirdDerivativeValue = (1.0 - inputEps) * thirdDerivative[inputIndex] + inputEps * thirdDerivative[inputIndex+1];
        double derivativeValue      = (1.0 - inputEps) *      derivative[inputIndex] + inputEps *      derivative[inputIndex+1];

        // Quadratic interpolation for function value
        double functionValue        = function[inputIndex] + 0.5 * (derivative[inputIndex] + derivativeValue) * inputEps * inputSpacing;

        functionData_.push_back( functionValue );
        derivativeData_.push_back( derivativeValue - spacing * spacing * thirdDerivativeValue / 12.0 );
    }

    if (checkVectorsAreFinite(functionData_, derivativeData_) == false)
    {
        GMX_THROW(APIError("Tabluated function has too large value or derivative in requested range"));
    }

    derivativeDeltaFunctionZeroData_ = makeDdfzVector(functionData_, derivativeData_);
}



QuadraticSplineTable::QuadraticSplineTable(const std::vector<double>     &function,
                                           double                         inputSpacing,
                                           const std::pair<real, real>    &range,
                                           real                           tolerance)
    : range_(range)
{
    // Create derivative vector numerically
    if (function.size() < 6)
    {
        GMX_THROW(APIError("Table input vector is too short for differentiation"));
    }

    // Create the missing derivative data
    std::vector<double> derivative = vectorDerivative(function, inputSpacing);

    // Construct a different temporary table using the constructor taking a derivative vector
    QuadraticSplineTable  tmpTable(function, derivative, inputSpacing, range, tolerance);

    // Move data from tmpTable to this table. The range_ object has already been set.
    halfSpacing_                     = tmpTable.halfSpacing_;
    tableScale_                      = tmpTable.tableScale_;
    functionData_                    = std::move(tmpTable.functionData_);
    derivativeData_                  = std::move(tmpTable.derivativeData_);
    derivativeDeltaFunctionZeroData_ = std::move(tmpTable.derivativeDeltaFunctionZeroData_);
}

}
