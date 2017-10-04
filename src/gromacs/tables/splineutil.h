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
 * Declares internal utility functions for spline tables
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#ifndef GMX_TABLES_SPLINEUTIL_H
#define GMX_TABLES_SPLINEUTIL_H

#include <functional>
#include <utility>
#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

namespace internal
{

/*! \brief Ensure analytical derivative is the derivative of analytical function.
 *
 *  This routine evaluates the numerical derivative of the function for
 *  a few (1000) points in the interval and checks that the relative difference
 *  between numerical and analytical derivative is within the expected error
 *  for the numerical derivative approximation we use.
 *
 *  The main point of this routine is to make sure the user has not made a
 *  mistake or sign error when defining the functions.
 *
 *  \param function   Analytical function to differentiate
 *  \param derivative Analytical derivative to compare with
 *  \param range      Range to test
 *
 *  \throws If the provided derivative does not seem to match the function.
 *
 *  \note The function/derivative are always double-valued to avoid accuracy loss.
 */
void
throwUnlessDerivativeIsConsistentWithFunction(const std::function<double(double)>  &function,
                                              const std::function<double(double)>  &derivative,
                                              const std::pair<real, real>          &range);

/*! \brief Ensure vector of derivative values is the derivative of function vector.
 *
 *  This routine differentiates a vector of numerical values and checks
 *  that the relative difference to a provided vector of numerical derivatives
 *  is smaller than the expected error from the numerical differentiation.
 *
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
 *  \throws If the provided derivative does not seem to match the function.
 *
 *  \note The function/derivative vectors and spacing are always double-valued
 *        to avoid accuracy loss.
 */
void
throwUnlessDerivativeIsConsistentWithFunction(ArrayRef<const double>        function,
                                              ArrayRef<const double>        derivative,
                                              double                        inputSpacing,
                                              const std::pair<real, real>  &range);




/*! \brief Find smallest quotient between analytical function and its 2nd derivative
 *
 *  Used to calculate spacing for quadratic spline tables. This function divides the
 *  function value by the second derivative (or a very small number when that is zero),
 *  and returns the smallest such quotient found in the range.
 *
 *  Our quadratic tables corresponds to linear interpolation of the derivative,
 *  which means the derivative will typically have larger error than the value
 *  when interpolating. The spacing required to reach a particular relative
 *  tolerance in the derivative depends on the quotient between the first
 *  derivative and the third derivative of the function itself.
 *
 *  You should call this routine with the analytical derivative as the "function"
 *  parameter, and the quotient between "function and second derivative" will
 *  then correspond to the quotient bewteen the derivative and the third derivative
 *  of the actual function we want to tabulate.
 *
 *  Since all functions that can be tabulated efficiently are reasonably smooth,
 *  we simply check 1,000 points in the interval rather than bother about
 *  implementing any complicated global optimization scheme.
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
                                                  const std::pair<real, real>          &range);





/*! \brief Find smallest quotient between vector of values and its 2nd derivative
 *
 *  Used to calculate spacing for quadratic spline tables. This function divides the
 *  function value by the second derivative (or a very small number when that is zero),
 *  and returns the smallest such quotient found in the range.
 *
 *  Our quadratic tables corresponds to linear interpolation of the derivative,
 *  which means the derivative will typically have larger error than the value
 *  when interpolating. The spacing required to reach a particular relative
 *  tolerance in the derivative depends on the quotient between the first
 *  derivative and the third derivative of the function itself.
 *
 *  You should call this routine with the analytical derivative as the "function"
 *  parameter, and the quotient between "function and second derivative" will
 *  then correspond to the quotient bewteen the derivative and the third derivative
 *  of the actual function we want to tabulate.
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
findSmallestQuotientOfFunctionAndSecondDerivative(ArrayRef<const double>        function,
                                                  double                        inputSpacing,
                                                  const std::pair<real, real>  &range);





/*! \brief Find smallest quotient between analytical function and its 3rd derivative
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
 *  You should call this routine with the analytical derivative as the "function"
 *  parameter, and the quotient between "function and second derivative" will
 *  then correspond to the quotient bewteen the derivative and the third derivative
 *  of the actual function we want to tabulate.
 *
 *  Since all functions that can be tabulated efficiently are reasonably smooth,
 *  we simply check 1,000 points in the interval rather than bother about
 *  implementing any complicated global optimization scheme.
 *
 *  \param f          Analytical function
 *  \param range      Interval
 *
 *  \return Smallest quotient found in range.
 *
 *  \note The function is always double-valued to avoid accuracy loss.
 */
real
findSmallestQuotientOfFunctionAndThirdDerivative(const std::function<double(double)>  &f,
                                                 const std::pair<real, real>          &range);




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
 *  You should call this routine with the analytical derivative as the "function"
 *  parameter, and the quotient between "function and second derivative" will
 *  then correspond to the quotient bewteen the derivative and the third derivative
 *  of the actual function we want to tabulate.
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
findSmallestQuotientOfFunctionAndThirdDerivative(ArrayRef<const double>        function,
                                                 double                        inputSpacing,
                                                 const std::pair<real, real>  &range);


/*! \brief Calculate second derivative of vector and return vector of same length
 *
 *  5-point approximations are used, with endpoints using non-center interpolation.
 *
 *  \param f       Vector (function) for which to calculate second derivative
 *  \param spacing Spacing of input data.
 *
 *  \throws If the input vector has fewer than five data points.
 *
 * \note This function always uses double precision arguments since it is meant
 *       to be used on raw user input data for tables, where we want to avoid
 *       accuracy loss (since differentiation can be numerically fragile).
 */
std::vector<double>
vectorSecondDerivative(ArrayRef<const double>       f,
                       double                       spacing);


/*! \brief Copy (temporary) table data into aligned multiplexed vector
 *
 *  This routine takes the temporary data generated for a single table
 *  and writes multiplexed output into a multiple-table-data vector.
 *  If the output vector is empty we will resize it to fit the data, and
 *  otherwise we assert the size is correct to add out input data.
 *
 *  \tparam T     Type of container for input data
 *  \tparam U     Type of container for output data
 *
 *  \param[in]    inputData               Input data for single table
 *  \param[inout] multiplexedOutputData   Multiplexed output vector, many tables.
 *  \param[in]    valuesPerTablePoint     Number of real values for each table point,
 *                                        for instance 4 in DDFZ tables.
 *  \param[in]    numTables               Number of tables mixed into multiplexed output
 *  \param[in]    thisTableIndex          Index of this table in multiplexed output
 *
 *  \note The output container type can be different from the input since the latter
 *        sometimes uses an aligned allocator so the data can be loaded efficiently
 *        in the GROMACS nonbonded kernels.
 */
template<class T, class U>
void
fillMultiplexedTableData(const T           inputData,
                         U *               multiplexedOutputData,
                         std::size_t       valuesPerTablePoint,
                         std::size_t       numTables,
                         std::size_t       thisTableIndex)
{
    if (multiplexedOutputData->size() == 0)
    {
        multiplexedOutputData->resize( inputData.size() * numTables );
    }
    else
    {
        GMX_ASSERT(inputData.size() * numTables == multiplexedOutputData->size(),
                   "Size mismatch when multiplexing table data");
    }

    GMX_ASSERT(inputData.size() % valuesPerTablePoint == 0,
               "Input table size must be a multiple of valuesPerTablePoint");

    std::size_t points = inputData.size() / valuesPerTablePoint;

    for (std::size_t i = 0; i < points; i++)
    {
        std::size_t inputOffset  = valuesPerTablePoint * i;
        std::size_t outputOffset = valuesPerTablePoint * ( numTables * i + thisTableIndex );

        for (std::size_t j = 0; j < valuesPerTablePoint; j++)
        {
            (*multiplexedOutputData)[outputOffset + j] = inputData[inputOffset + j];
        }
    }
}


}      // namespace internal

}      // namespace gmx

#endif // GMX_TABLES_SPLINEUTIL_H
