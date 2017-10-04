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
 * Implements classes for cubic spline table functions
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "cubicsplinetable.h"

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

/*! \brief Calculate table elements from function/derivative data
 *
 * \param functionValue0   Function value for the present table index
 * \param functionValue1   Function value for the next table index
 * \param derivativeValue0 Derivative value for the present table index
 * \param derivativeValue1 Derivative value for the next table index
 * \param spacing          Distance between table points
 * \param Y                Function value for table index
 * \param F                Component to multiply with offset eps
 * \param G                Component to multiply with eps^2
 * \param H                Component to multiply with eps^3
 */
void
calculateCubicSplineCoefficients(double  functionValue0,
                                 double  functionValue1,
                                 double  derivativeValue0,
                                 double  derivativeValue1,
                                 double  spacing,
                                 double *Y,
                                 double *F,
                                 double *G,
                                 double *H)
{
    *Y  =  functionValue0;
    *F  =  spacing * derivativeValue0;
    *G  =  3.0*( functionValue1 - functionValue0) - spacing * (derivativeValue1 + 2.0 * derivativeValue0);
    *H  = -2.0*( functionValue1 - functionValue0) + spacing * (derivativeValue1 + derivativeValue0);
}

/*! \brief Perform cubic spline interpolation in interval from function/derivative
 *
 * \param      functionValue0               Function value for the present table index
 * \param      functionValue1               Function value for the next table index
 * \param      derivativeValue0             Derivative value for the present table index
 * \param      derivativeValue1             Derivative value for the next table index
 * \param      spacing                      Distance between table points
 * \param      eps                          Offset from lower table point for evaluation
 * \param[out] interpolatedFunctionValue    Output function value
 * \param[out] interpolatedDerivativeValue  Output derivative value
 */
void
cubicSplineInterpolationFromFunctionAndDerivative(double  functionValue0,
                                                  double  functionValue1,
                                                  double  derivativeValue0,
                                                  double  derivativeValue1,
                                                  double  spacing,
                                                  double  eps,
                                                  double *interpolatedFunctionValue,
                                                  double *interpolatedDerivativeValue)
{
    double Y, F, G, H;

    calculateCubicSplineCoefficients(functionValue0, functionValue1,
                                     derivativeValue0, derivativeValue1,
                                     spacing,
                                     &Y, &F, &G, &H);

    double Fp = fma(fma(H, eps, G), eps, F);

    *interpolatedFunctionValue   = fma(Fp, eps, Y);
    *interpolatedDerivativeValue = fma(eps, fma(2.0*eps, H, G), Fp)/spacing;
}



/*! \brief Construct the data for a single cubic table from analytical functions
 *
 * \param[in]  function             Analytical functiojn
 * \param[in]  derivative           Analytical derivative
 * \param[in]  range                Upper/lower limit of region to tabulate
 * \param[in]  spacing              Distance between table points
 * \param[out] yfghTableData        Output cubic spline table with Y,F,G,H entries
 */
void
fillSingleCubicSplineTableData(const std::function<double(double)>   &function,
                               const std::function<double(double)>   &derivative,
                               const std::pair<real, real>           &range,
                               double                                 spacing,
                               std::vector<real>                     *yfghTableData)
{
    int  endIndex   = range.second / spacing + 2;

    yfghTableData->resize(4*endIndex);

    double       maxMagnitude      = 0.0001*GMX_REAL_MAX;
    bool         functionIsInRange = true;
    std::size_t  lastIndexInRange  = endIndex - 1;

    for (int i = endIndex - 1; i >= 0; i--)
    {
        double x                    = i * spacing;
        double tmpFunctionValue;
        double tmpDerivativeValue;
        double nextHigherFunction;
        double nextHigherDerivative;
        double Y, F, G, H;

        if (range.first > 0 && i == 0)
        {
            // Avoid x==0 if it is not in the range, since it can lead to
            // singularities even if the value for i==1 was within or required magnitude
            functionIsInRange = false;
        }

        if (functionIsInRange)
        {
            tmpFunctionValue     = function(x);
            tmpDerivativeValue   = derivative(x);
            nextHigherFunction   = ((i+1) < endIndex) ? function(x+spacing) : 0.0;
            nextHigherDerivative = ((i+1) < endIndex) ? derivative(x+spacing) : 0.0;

            if (std::abs(tmpFunctionValue) > maxMagnitude || std::abs(tmpDerivativeValue) > maxMagnitude)
            {
                functionIsInRange = false; // Once this happens, it never resets to true again
            }
        }

        if (functionIsInRange)
        {
            calculateCubicSplineCoefficients(tmpFunctionValue, nextHigherFunction,
                                             tmpDerivativeValue, nextHigherDerivative,
                                             spacing,
                                             &Y, &F, &G, &H);
            lastIndexInRange--;
        }
        else
        {
            double lastIndexY = (*yfghTableData)[4*lastIndexInRange];
            double lastIndexF = (*yfghTableData)[4*lastIndexInRange + 1];

            Y = lastIndexY + lastIndexF * (i - lastIndexInRange);
            F = lastIndexF;
            G = 0.0;
            H = 0.0;
        }

        (*yfghTableData)[4*i  ] = Y;
        (*yfghTableData)[4*i+1] = F;
        (*yfghTableData)[4*i+2] = G;
        (*yfghTableData)[4*i+3] = H;
    }
}


/*! \brief Construct the data for a single cubic table from vector data
 *
 * \param[in]  function             Input vector with function data
 * \param[in]  derivative           Input vector with derivative data
 * \param[in]  inputSpacing         Distance between points in input vectors
 * \param[in]  range                Upper/lower limit of region to tabulate
 * \param[in]  spacing              Distance between table points
 * \param[out] yfghTableData        Output cubic spline table with Y,F,G,H entries
 */
void
fillSingleCubicSplineTableData(ArrayRef<const double>                 function,
                               ArrayRef<const double>                 derivative,
                               double                                 inputSpacing,
                               const std::pair<real, real>           &range,
                               double                                 spacing,
                               std::vector<real>                     *yfghTableData)
{
    int                 endIndex   = range.second / spacing + 2;

    std::vector<double> tmpFunction(endIndex);
    std::vector<double> tmpDerivative(endIndex);

    double              maxMagnitude      = 0.0001*GMX_REAL_MAX;
    bool                functionIsInRange = true;
    std::size_t         lastIndexInRange  = endIndex - 1;

    // Interpolate function and derivative values in positions needed for output
    for (int i = endIndex - 1; i >= 0; i--)
    {
        double x     = i * spacing;
        double xtab  = x / inputSpacing;
        int    index = xtab;
        double eps   = xtab - index;

        if (range.first > 0 && i == 0)
        {
            // Avoid x==0 if it is not in the range, since it can lead to
            // singularities even if the value for i==1 was within or required magnitude
            functionIsInRange = false;
        }

        if (functionIsInRange && (std::abs(function[index]) > maxMagnitude || std::abs(derivative[index]) > maxMagnitude))
        {
            functionIsInRange = false; // Once this happens, it never resets to true again
        }

        if (functionIsInRange)
        {
            cubicSplineInterpolationFromFunctionAndDerivative(function[index],
                                                              function[index+1],
                                                              derivative[index],
                                                              derivative[index+1],
                                                              inputSpacing,
                                                              eps,
                                                              &(tmpFunction[i]),
                                                              &(tmpDerivative[i]));
            lastIndexInRange--;
        }
        else
        {
            double lastIndexFunction   = tmpFunction[lastIndexInRange];
            double lastIndexDerivative = tmpDerivative[lastIndexInRange];
            tmpFunction[i]             = lastIndexFunction + lastIndexDerivative * (i - lastIndexInRange) * spacing;
            tmpDerivative[i]           = lastIndexDerivative;
        }
    }

    yfghTableData->resize(4*endIndex);

    for (int i = 0; i < endIndex; i++)
    {
        double Y, F, G, H;

        double nextFunction   = ((i+1) < endIndex) ? tmpFunction[i+1] : 0.0;
        double nextDerivative = ((i+1) < endIndex) ? tmpDerivative[i+1] : 0.0;

        calculateCubicSplineCoefficients(tmpFunction[i], nextFunction,
                                         tmpDerivative[i], nextDerivative,
                                         spacing,
                                         &Y, &F, &G, &H);
        (*yfghTableData)[4*i  ] = Y;
        (*yfghTableData)[4*i+1] = F;
        (*yfghTableData)[4*i+2] = G;
        (*yfghTableData)[4*i+3] = H;
    }

}

}   // namespace anonymous



#if GMX_DOUBLE
const real
CubicSplineTable::defaultTolerance = 1e-10;
#else
const real
CubicSplineTable::defaultTolerance = 10.0 * GMX_FLOAT_EPS;
#endif

CubicSplineTable::CubicSplineTable(std::initializer_list<AnalyticalSplineTableInput>   analyticalInputList,
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
            ex.prependContext("Error generating cubic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
        // Calculate the required table spacing h. The error we make with linear interpolation
        // of the derivative will be described by the third-derivative correction term.
        // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
        // where f'/f''' is the first and third derivative of the function, respectively.

        double thisMinQuotient = internal::findSmallestQuotientOfFunctionAndThirdDerivative(thisFuncInput.derivative, range_);

        minQuotient = std::min(minQuotient, thisMinQuotient);
    }

    double spacing = 0.5 * std::cbrt(72.0 * std::sqrt(3.0) * tolerance * minQuotient);

    tableScale_  = 1.0 / spacing;

    if (range_.second * tableScale_ > 2e6)
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
            std::vector<real> tmpYfghTableData;

            fillSingleCubicSplineTableData(thisFuncInput.function,
                                           thisFuncInput.derivative,
                                           range_,
                                           spacing,
                                           &tmpYfghTableData);

            internal::fillMultiplexedTableData(tmpYfghTableData, &yfghMultiTableData_,
                                               4, numFuncInTable_, funcIndex);

            funcIndex++;
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating cubic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
    }
}


CubicSplineTable::CubicSplineTable(std::initializer_list<NumericalSplineTableInput>   numericalInputList,
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
            ex.prependContext("Error generating cubic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
        // Calculate the required table spacing h. The error we make with linear interpolation
        // of the derivative will be described by the third-derivative correction term.
        // This means we can compute the required spacing as h = sqrt(12*tolerance*min(f'/f''')),
        // where f'/f''' is the first and third derivative of the function, respectively.

        double thisMinQuotient = internal::findSmallestQuotientOfFunctionAndThirdDerivative(thisFuncInput.derivative, thisFuncInput.spacing, range_);

        minQuotient = std::min(minQuotient, thisMinQuotient);
    }

    double spacing     = std::cbrt(72.0 * std::sqrt(3.0) * tolerance * minQuotient);

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

            std::vector<real> tmpYfghTableData;

            fillSingleCubicSplineTableData(thisFuncInput.function,
                                           thisFuncInput.derivative,
                                           thisFuncInput.spacing,
                                           range,
                                           spacing,
                                           &tmpYfghTableData);

            internal::fillMultiplexedTableData(tmpYfghTableData, &yfghMultiTableData_,
                                               4, numFuncInTable_, funcIndex);

            funcIndex++;
        }
        catch (gmx::GromacsException &ex)
        {
            ex.prependContext("Error generating cubic spline table for function '" + thisFuncInput.desc + "'");
            throw;
        }
    }
}

} // namespace gmx
