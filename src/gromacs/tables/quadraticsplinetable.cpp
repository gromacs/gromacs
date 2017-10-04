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

/*! \internal \file
 * \brief
 * Implements classes for quadratic spline table functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "quadraticsplinetable.h"

#include <cmath>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <utility>
#include <vector>

#include "gromacs/tables/tableinput.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

#include "splineutil.h"

namespace gmx
{

namespace
{

/*! \brief Construct the data for a single quadratic table from analytical functions
 *
 * \param[in]  function             Analytical functiojn
 * \param[in]  derivative           Analytical derivative
 * \param[in]  range                Upper/lower limit of region to tabulate
 * \param[in]  spacing              Distance between table points
 * \param[out] functionTableData    Output table with function data
 * \param[out] derivativeTableData  OUtput table with (adjusted) derivative data
 */
void
fillSingleQuadraticSplineTableData(const std::function<double(double)>   &function,
                                   const std::function<double(double)>   &derivative,
                                   const std::pair<real, real>           &range,
                                   double                                 spacing,
                                   std::vector<real>                     *functionTableData,
                                   std::vector<real>                     *derivativeTableData)
{
    std::size_t  endIndex   = range.second / spacing + 2;

    functionTableData->resize(endIndex);
    derivativeTableData->resize(endIndex);

    double       maxMagnitude      = 0.0001*GMX_REAL_MAX;
    bool         functionIsInRange = true;
    std::size_t  lastIndexInRange  = endIndex - 1;

    for (int i = endIndex - 1; i >= 0; i--)
    {
        double x                = i * spacing;
        double tmpFunctionValue;
        double tmpDerivativeValue;

        if (range.first > 0 && i == 0)
        {
            // Avoid x==0 if it is not in the range, since it can lead to
            // singularities even if the value for i==1 was within or required magnitude
            functionIsInRange = false;
        }

        if (functionIsInRange)
        {
            tmpFunctionValue = function(x);

            // Calculate third derivative term (2nd derivative of the derivative)
            // Make sure we stay in range. In practice this means we use one-sided
            // interpolation at the interval endpoints (indentical to an offset for 3-point formula)
            const double h                    = std::pow( GMX_DOUBLE_EPS, 0.25 );
            double       y                    = std::min( std::max(x, range.first + h), range.second - h);
            double       thirdDerivativeValue = ( derivative(y+h) - 2.0 * derivative(y) + derivative(y-h) ) / ( h * h );

            tmpDerivativeValue   = derivative(x) - spacing * spacing * thirdDerivativeValue / 12.0;

            if (std::abs(tmpFunctionValue) > maxMagnitude || std::abs(tmpDerivativeValue) > maxMagnitude)
            {
                functionIsInRange = false; // Once this happens, it never resets to true again
            }
        }

        if (functionIsInRange)
        {
            (*functionTableData)[i]   = tmpFunctionValue;
            (*derivativeTableData)[i] = tmpDerivativeValue;
            lastIndexInRange--;
        }
        else
        {
            // Once the function or derivative (more likely) has reached very large values,
            // we simply make a linear function from the last in-range value of the derivative.
            double lastIndexFunction   = (*functionTableData)[lastIndexInRange];
            double lastIndexDerivative = (*derivativeTableData)[lastIndexInRange];
            (*functionTableData)[i]    = lastIndexFunction + lastIndexDerivative * (i - lastIndexInRange) * spacing;
            (*derivativeTableData)[i]  = lastIndexDerivative;
        }
    }
}


/*! \brief Construct the data for a single quadratic table from vector data
 *
 * \param[in]  function             Input vector with function data
 * \param[in]  derivative           Input vector with derivative data
 * \param[in]  inputSpacing         Distance between points in input vectors
 * \param[in]  range                Upper/lower limit of region to tabulate
 * \param[in]  spacing              Distance between table points
 * \param[out] functionTableData    Output table with function data
 * \param[out] derivativeTableData  OUtput table with (adjusted) derivative data
 */
void
fillSingleQuadraticSplineTableData(ArrayRef<const double>                 function,
                                   ArrayRef<const double>                 derivative,
                                   double                                 inputSpacing,
                                   const std::pair<real, real>           &range,
                                   double                                 spacing,
                                   std::vector<real>                     *functionTableData,
                                   std::vector<real>                     *derivativeTableData)
{
    std::size_t  endIndex   = range.second / spacing + 2;

    functionTableData->resize(endIndex);
    derivativeTableData->resize(endIndex);

    std::vector<double>  thirdDerivative(internal::vectorSecondDerivative(derivative, inputSpacing));

    double               maxMagnitude      = 0.0001*GMX_REAL_MAX;
    bool                 functionIsInRange = true;
    std::size_t          lastIndexInRange  = endIndex - 1;

    for (int i = endIndex - 1; i >= 0; i--)
    {
        double x                = i * spacing;
        double tmpFunctionValue;
        double tmpDerivativeValue;

        if (range.first > 0 && i == 0)
        {
            // Avoid x==0 if it is not in the range, since it can lead to
            // singularities even if the value for i==1 was within or required magnitude
            functionIsInRange = false;
        }

        if (functionIsInRange)
        {
            // Step 1: Interpolate the function value at x from input table.
            double inputXTab  = x / inputSpacing;
            int    inputIndex = inputXTab;
            double inputEps   = inputXTab - inputIndex;

            // Linear interpolation of input derivative and third derivative
            double thirdDerivativeValue = (1.0 - inputEps) * thirdDerivative[inputIndex] + inputEps * thirdDerivative[inputIndex+1];
            double derivativeValue      = (1.0 - inputEps) *      derivative[inputIndex] + inputEps *      derivative[inputIndex+1];

            // Quadratic interpolation for function value
            tmpFunctionValue     = function[inputIndex] + 0.5 * (derivative[inputIndex] + derivativeValue) * inputEps * inputSpacing;
            tmpDerivativeValue   = derivativeValue - spacing * spacing * thirdDerivativeValue / 12.0;

            if (std::abs(tmpFunctionValue) > maxMagnitude || std::abs(tmpDerivativeValue) > maxMagnitude)
            {
                functionIsInRange = false; // Once this happens, it never resets to true again
            }
        }

        if (functionIsInRange)
        {
            (*functionTableData)[i]   = tmpFunctionValue;
            (*derivativeTableData)[i] = tmpDerivativeValue;
            lastIndexInRange--;
        }
        else
        {
            // Once the function or derivative (more likely) has reached very large values,
            // we simply make a linear function from the last in-range value of the derivative.
            double lastIndexFunction   = (*functionTableData)[lastIndexInRange];
            double lastIndexDerivative = (*derivativeTableData)[lastIndexInRange];
            (*functionTableData)[i]    = lastIndexFunction + lastIndexDerivative * (i - lastIndexInRange) * spacing;
            (*derivativeTableData)[i]  = lastIndexDerivative;
        }
    }
}

/*! \brief Create merged DDFZ vector from function & derivative data
 *
 *  \param functionTableData     Function values
 *  \param derivativeTableData   Derivative values. We have already subtracted the
 *                               small third derivative component when calling this
 *                               function, but in practice it is just an arbitrary
 *                               vector here.
 *  \param ddfzTableData         Vector four times longer, filled with
 *                               the derivative, the difference to the next derivative
 *                               point, the function value, and zero.
 *
 *  \throws If the vector lengths do not match.
 */
void
fillDdfzTableData(const std::vector<real>    &functionTableData,
                  const std::vector<real>    &derivativeTableData,
                  std::vector<real>          *ddfzTableData)
{
    GMX_ASSERT(functionTableData.size() == derivativeTableData.size(), "Mismatching vector lengths");

    std::size_t points = functionTableData.size();

    ddfzTableData->resize(4 * points);

    for (std::size_t i = 0; i < points; i++)
    {
        (*ddfzTableData)[4*i]     = derivativeTableData[i];

        double nextDerivative     = ( i < functionTableData.size() - 1 ) ? derivativeTableData[i+1] : 0.0;

        (*ddfzTableData)[4*i + 1] = nextDerivative - derivativeTableData[i];
        (*ddfzTableData)[4*i + 2] = functionTableData[i];
        (*ddfzTableData)[4*i + 3] = 0.0;
    }
}

}   // namespace anonymous



const real
QuadraticSplineTable::defaultTolerance = 10.0 * GMX_FLOAT_EPS;


QuadraticSplineTable::QuadraticSplineTable(std::initializer_list<AnalyticalSplineTableInput>   analyticalInputList,
                                           const std::pair<real, real>                        &range,
                                           real                                                tolerance)
    : numFuncInTable_(analyticalInputList.size()), range_(range)
{
    // Sanity check on input values
    if (range_.first < 0.0 || (range_.second-range_.first) < 0.001)
    {
        GMX_THROW(InvalidInputError("Range to tabulate cannot include negative values and must span at least 0.001"));
    }

    if (tolerance < GMX_REAL_EPS)
    {
        GMX_THROW(ToleranceError("Table tolerance cannot be smaller than GMX_REAL_EPS"));
    }

    double minQuotient = GMX_REAL_MAX;

    // loop over all functions to find smallest spacing
    for (auto thisFuncInput : analyticalInputList)
    {
        try
        {
            internal::throwUnlessDerivativeIsConsistentWithFunction(thisFuncInput.function, thisFuncInput.derivative, range_);
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating quadratic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
        // Calculate the required table spacing h. The error we make with linear interpolation
        // of the derivative will be described by the third-derivative correction term.
        // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
        // where f'/f''' is the first and third derivative of the function, respectively.

        double thisMinQuotient = internal::findSmallestQuotientOfFunctionAndSecondDerivative(thisFuncInput.derivative, range_);

        minQuotient = std::min(minQuotient, thisMinQuotient);
    }

    double spacing = std::sqrt(12.0 * tolerance * minQuotient);

    halfSpacing_ = 0.5 * spacing;
    tableScale_  = 1.0 / spacing;

    if (range_.second * tableScale_ > 1e6)
    {
        GMX_THROW(ToleranceError("Over a million points would be required for table; decrease range or increase tolerance"));
    }

    // Loop over all tables again.
    // Here we create the actual table for each function, and then
    // combine them into a multiplexed table function.
    std::size_t funcIndex = 0;

    for (auto thisFuncInput : analyticalInputList)
    {
        try
        {
            std::vector<real> tmpFuncTableData;
            std::vector<real> tmpDerTableData;
            std::vector<real> tmpDdfzTableData;

            fillSingleQuadraticSplineTableData(thisFuncInput.function,
                                               thisFuncInput.derivative,
                                               range_,
                                               spacing,
                                               &tmpFuncTableData,
                                               &tmpDerTableData);

            fillDdfzTableData(tmpFuncTableData, tmpDerTableData, &tmpDdfzTableData);

            internal::fillMultiplexedTableData(tmpDerTableData, &derivativeMultiTableData_,
                                               1, numFuncInTable_, funcIndex);

            internal::fillMultiplexedTableData(tmpDdfzTableData, &ddfzMultiTableData_,
                                               4, numFuncInTable_, funcIndex);

            funcIndex++;
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating quadratic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
    }
}


QuadraticSplineTable::QuadraticSplineTable(std::initializer_list<NumericalSplineTableInput>   numericalInputList,
                                           const std::pair<real, real>                       &range,
                                           real                                               tolerance)
    : numFuncInTable_(numericalInputList.size()), range_(range)
{
    // Sanity check on input values
    if (range.first < 0.0 || (range.second-range.first) < 0.001)
    {
        GMX_THROW(InvalidInputError("Range to tabulate cannot include negative values and must span at least 0.001"));
    }

    if (tolerance < GMX_REAL_EPS)
    {
        GMX_THROW(ToleranceError("Table tolerance cannot be smaller than GMX_REAL_EPS"));
    }

    double minQuotient = GMX_REAL_MAX;

    // loop over all functions to find smallest spacing
    for (auto thisFuncInput : numericalInputList)
    {
        try
        {
            // We do not yet know what the margin is, but we need to test that we at least cover
            // the requested range before starting to calculate derivatives
            if (thisFuncInput.function.size() < range_.second / thisFuncInput.spacing + 1)
            {
                GMX_THROW(InconsistentInputError("Table input vectors must cover requested range, and a margin beyond the upper endpoint"));
            }

            if (thisFuncInput.function.size() != thisFuncInput.derivative.size())
            {
                GMX_THROW(InconsistentInputError("Function and derivative vectors have different lengths"));
            }

            internal::throwUnlessDerivativeIsConsistentWithFunction(thisFuncInput.function, thisFuncInput.derivative, thisFuncInput.spacing, range_);
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating quadratic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
        // Calculate the required table spacing h. The error we make with linear interpolation
        // of the derivative will be described by the third-derivative correction term.
        // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
        // where f'/f''' is the first and third derivative of the function, respectively.

        double thisMinQuotient = internal::findSmallestQuotientOfFunctionAndSecondDerivative(thisFuncInput.derivative, thisFuncInput.spacing, range_);

        minQuotient = std::min(minQuotient, thisMinQuotient);
    }

    double spacing     = std::sqrt(12.0 * tolerance * minQuotient);

    halfSpacing_ = 0.5 * spacing;
    tableScale_  = 1.0 / spacing;

    if (range_.second * tableScale_ > 1e6)
    {
        GMX_THROW(ToleranceError("Requested tolerance would require over a million points in table"));
    }

    // Loop over all tables again.
    // Here we create the actual table for each function, and then
    // combine them into a multiplexed table function.
    std::size_t funcIndex = 0;

    for (auto thisFuncInput : numericalInputList)
    {
        try
        {
            if (spacing < thisFuncInput.spacing)
            {
                GMX_THROW(ToleranceError("Input vector spacing cannot achieve tolerance requested"));
            }

            std::vector<real> tmpFuncTableData;
            std::vector<real> tmpDerTableData;
            std::vector<real> tmpDdfzTableData;

            fillSingleQuadraticSplineTableData(thisFuncInput.function,
                                               thisFuncInput.derivative,
                                               thisFuncInput.spacing,
                                               range,
                                               spacing,
                                               &tmpFuncTableData,
                                               &tmpDerTableData);

            fillDdfzTableData(tmpFuncTableData, tmpDerTableData, &tmpDdfzTableData);

            internal::fillMultiplexedTableData(tmpDerTableData, &derivativeMultiTableData_,
                                               1, numFuncInTable_, funcIndex);

            internal::fillMultiplexedTableData(tmpDdfzTableData, &ddfzMultiTableData_,
                                               4, numFuncInTable_, funcIndex);

            funcIndex++;
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating quadratic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
    }
}

} // namespace gmx
