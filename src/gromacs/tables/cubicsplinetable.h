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

/*! \libinternal \file
 * \brief
 * Declares classes for cubic spline table
 *
 * \inlibraryapi
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#ifndef GMX_TABLES_CUBICSPLINETABLE_H
#define GMX_TABLES_CUBICSPLINETABLE_H

#include <cstddef>

#include <initializer_list>
#include <memory>
#include <utility>
#include <vector>

#include "gromacs/libgromacs_export.h"
#include "gromacs/simd/simd.h"
#include "gromacs/tables/tableinput.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \libinternal \brief Cubic spline interpolation table.
 *
 * This class interpolates a function specified either as an analytical
 * expression or from user-provided table data.
 *
 * At initialization, you provide the reference function of vectors
 * as a list of tuples that contain a brief name, the function, and
 * derivative for each function to tabulate. To create a table with
 * two functions this initializer list can for instance look like
 *
 *     { {"LJ6", lj6Func, lj6Der}, {"LJ12", lj12Func, lj12Der} }
 *
 * The names are only used so exceptions during initialization can
 * be traced to a specific table.
 *
 * When interpolating, there are methods to interpolate either 1, 2, or 3
 * functions in one go. By default these interpolation routines will
 * operate on tables with the same number of functions as specified in
 * the interpolation method (debug builds check that this is consistent with
 * the table). However, it is also possible to use optional template
 * parameters that specify the total number of functions in a table, and
 * what function index to interpolate. For instance, to interpolate the
 * derivative of the second function (i.e., index 1) in a
 * multi-function-table with three functions in total, you can write
 *
 *     table.evaluateDerivative<3,1>(x,&der);
 *
 * Here too, debug builds will check that the template parameters are
 * consistent with the table.
 *
 * This class interpolates a function specified either as an analytical
 * expression or from user-provided table data. The coefficients for each
 * table point are precalculated such that we simply evaluate
 *
 * \f{eqnarray*}{
 * V(x)  = Y + F \epsilon + G \epsilon^2 + H \epsilon^3
 * V'(x) = (F + 2 G \epsilon + 3 H \epsilon^2)/h
 * \f}
 *
 * Where h is the spacing and epsilon the fractional offset from table point.
 *
 * While it is possible to create tables only from function values
 * (i.e., no derivatives), it is recommended to provide derivatives for higher
 * accuracy and to avoid issues with numerical differentiation. Note that the
 * table input should be smooth, i.e. it should not contain noise e.g. from an
 * (iterative) Boltzmann inversion procedure - you have been warned.
 *
 * \note This class is responsible for fundamental interpolation of any function,
 *       which might or might not correspond to a potential. For this reason
 *       both input and output derivatives are proper function derivatives, and
 *       we do not swap any signs to get forces directly from the table.
 *
 * \note There will be a small additional accuracy loss from the internal
 *       operation where we calculate the epsilon offset from the nearest table
 *       point, since the integer part we subtract can get large in those cases.
 *       While this is technically possible to solve with extended precision
 *       arithmetics, that would introduce extra instructions in some highly
 *       performance-sensitive code parts. For typical GROMACS interaction
 *       functions the derivatives will decay faster than the potential, which
 *       means it will never play any role. For other functions it will only
 *       cause a small increase in the relative error for arguments where the
 *       magnitude of the function or derivative is very small.
 *       Since we typically sum several results in GROMACS, this should never
 *       show up in any real cases, and for this reason we choose not to do
 *       the extended precision arithmetics.
 *
 * \note These routines are not suitable for table ranges starting far away
 *       from zero, since we allocate memory and calculate indices starting from
 *       range zero for efficiency reasons.
 */
class CubicSplineTable
{
private:
    /*! \brief Change that function value falls inside range when debugging
     *
     *  \tparam T   Lookup argument floating-point type, typically SimdReal or real.
     *  \param  r   Lookup argument to test
     *
     *  \throws Debug builds will throw gmx::RangeError for values that are
     *          larger than the upper limit of the range, or smaller than 0.
     *          We allow the table to be called with arguments between 0 and
     *          the lower limit of the range, since this might in theory occur
     *          once-in-a-blue-moon with some algorithms.
     */
    template<typename T>
    void rangeCheck(T gmx_unused r) const
    {
#ifndef NDEBUG
        // Check that all values fall in range when debugging
        if (anyTrue(r < T(0.0) || T(range_.second) <= r))
        {
            GMX_THROW(RangeError("Interpolation input value falls outside table definition range"));
        }
#endif
    }

public:
    /*! \brief Default tolerance for cubic spline tables
     *
     * This is 10*GMX_FLOAT_EPS in single precision, and
     * 1e-10 for double precision. It might not be worth setting
     * this tolerance lower than 1e-10 in double precision, both because
     * you will end up with very large tables, and because
     * functions like r^-12 become so large for small values of r the
     * table generation code will lead to some precision loss even
     * in double precision.
     */
    static LIBGROMACS_EXPORT const real defaultTolerance;

    /*! \brief Initialize table data from function
     *
     * \param analyticalInputList Initializer list with one or more functions to tabulate,
     *                            specified as elements with a string description and
     *                            the function as well as derivative. The function will also
     *                            be called for values smaller than the lower limit of the
     *                            range, but we avoid calling it for 0.0 if that value
     *                            is not included in the range.
     *                            Constructor will throw gmx::APIError for negative values.
     *                            Due to the way the numerical derivative evaluation depends
     *                            on machine precision internally, this range must be
     *                            at least 0.001, or the constructor throws gmx::APIError.
     * \param range                Range over which the function will be tabulated.
     *                             Constructor will throw gmx::APIError for negative values,
     *                             or if the value/derivative vector does not cover the
     *                             range.
     * \param tolerance           Requested accuracy of the table. This will be used to
     *                            calculate the required internal spacing. If this cannot
     *                            be achieved (for instance because the table would require
     *                            too much memory) the constructor will throw gmx::ToleranceError.
     *
     * \note The functions are always defined in double precision to avoid
     *       losing accuracy when constructing tables.
     *
     * \note Since we fill the table for values below range.first, you can achieve
     *       a smaller table by using a smaller range where the tolerance has to be
     *       met, and accept that a few function calls below range.first do not
     *       quite reach the tolerance.
     *
     * \warning For efficiency reasons (since this code is used in some inner
     *       (kernels), we always allocate memory and calculate table indices
     *       for the complete interval [0,range.second], although the data will
     *       not be valid outside the definition range to avoid calling the
     *       function there. This means you should \a not use this class
     *       to tabulate functions for small ranges very far away from zero,
     *       since you would both waste a huge amount of memory and incur
     *       truncation errors when calculating the index.
     *
     * \throws gmx::ToleranceError if the requested tolerance cannot be achieved,
     *         and gmx::APIError for other incorrect input.
     */
    CubicSplineTable(std::initializer_list<AnalyticalSplineTableInput> analyticalInputList,
                     const std::pair<real, real>&                      range,
                     real tolerance = defaultTolerance);

    /*! \brief Initialize table data from tabulated values and derivatives
     *
     * \param numericalInputList  Initializer list with one or more functions to tabulate,
     *                            specified as a string description, vectors with function and
     *                            derivative values, and the input spacing. Data points are
     *                            separated by the spacing parameter, starting from 0.
     *                            Values below the lower limit of the range will be used to
     *                            attempt defining the table, but we avoid using index 0
     *                            unless 0.0 is included in the range. Some extra points beyond
     *                            range.second are required to re-interpolate values, so add
     *                            some margin. The constructor will throw gmx::APIError if the
     *                            input vectors are too short to cover the requested range
     *                            (and they must always be at least five points).
     * \param range               Range over which the function will be tabulated.
     *                            Constructor will throw gmx::APIError for negative values,
     *                            or if the value/derivative vector does not cover the
     *                            range.
     * \param tolerance           Requested accuracy of the table. This will be used to
     *                            calculate the required internal spacing and possibly
     *                            re-interpolate. The constructor will throw
     *                            gmx::ToleranceError if the input spacing is too coarse
     *                            to achieve this accuracy.
     *
     * \note The input data vectors are always double precision to avoid
     *       losing accuracy when constructing tables.
     *
     * \note Since we fill the table for values below range.first, you can achieve
     *       a smaller table by using a smaller range where the tolerance has to be
     *       met, and accept that a few function calls below range.first do not
     *       quite reach the tolerance.
     *
     * \warning For efficiency reasons (since this code is used in some inner
     *       (kernels), we always allocate memory and calculate table indices
     *       for the complete interval [0,range.second], although the data will
     *       not be valid outside the definition range to avoid calling the
     *       function there. This means you should \a not use this class
     *       to tabulate functions for small ranges very far away from zero,
     *       since you would both waste a huge amount of memory and incur
     *       truncation errors when calculating the index.
     */
    CubicSplineTable(std::initializer_list<NumericalSplineTableInput> numericalInputList,
                     const std::pair<real, real>&                     range,
                     real                                             tolerance = defaultTolerance);


    /************************************************************
     *           Evaluation methods for single functions        *
     ************************************************************/

    /*! \brief Evaluate both function and derivative, single table function
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable  Number of separate functions in table, default is 1
     *  \tparam     funcIndex       Index of function to evaluate in table, default is 0
     *  \tparam     T               Type (SimdReal or real) of lookup and result
     *  \param      r               Points for which to evaluate function and derivative
     *  \param[out] functionValue   Function value
     *  \param[out] derivativeValue Function derivative
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 1, int funcIndex = 0, typename T>
    void evaluateFunctionAndDerivative(T r, T* functionValue, T* derivativeValue) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex, tabIndex, &Y, &F, &G, &H);
        *functionValue   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /*! \brief Evaluate function value only, single table function
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable  Number of separate functions in table, default is 1
     *  \tparam     funcIndex       Index of function to evaluate in table, default is 0
     *  \tparam     T               Type (SimdReal or real) of lookup and result
     *  \param      r               Points for which to evaluate function value
     *  \param[out] functionValue   Function value
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 1, int funcIndex = 0, typename T>
    void evaluateFunction(T r, T* functionValue) const
    {
        T der gmx_unused;

        evaluateFunctionAndDerivative<numFuncInTable, funcIndex>(r, functionValue, &der);
    }

    /*! \brief Evaluate function derivative only, single table function
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable  Number of separate functions in table, default is 1
     *  \tparam     funcIndex       Index of function to evaluate in table, default is 0
     *  \tparam     T               Type (SimdReal or real) of lookup and result
     *  \param      r               Points for which to evaluate function derivative
     *  \param[out] derivativeValue Function derivative
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 1, int funcIndex = 0, typename T>
    void evaluateDerivative(T r, T* derivativeValue) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex, tabIndex, &Y, &F, &G, &H);
        *derivativeValue = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /************************************************************
     *             Evaluation methods for two functions         *
     ************************************************************/

    /*! \brief Evaluate both function and derivative, two table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable  Number of separate functions in table, default is 2
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function and derivative
     *  \param[out] functionValue0   Interpolated value for first function
     *  \param[out] derivativeValue0 Interpolated derivative for first function
     *  \param[out] functionValue1   Interpolated value for second function
     *  \param[out] derivativeValue1 Interpolated derivative for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateFunctionAndDerivative(T r, T* functionValue0, T* derivativeValue0, T* functionValue1, T* derivativeValue1) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex0, tabIndex, &Y, &F, &G, &H);
        *functionValue0   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue0 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex1, tabIndex, &Y, &F, &G, &H);
        *functionValue1   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue1 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /*! \brief Evaluate function value only, two table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable   Number of separate functions in table, default is 2
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function value
     *  \param[out] functionValue0   Interpolated value for first function
     *  \param[out] functionValue1   Interpolated value for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateFunction(T r, T* functionValue0, T* functionValue1) const
    {
        T der0 gmx_unused;
        T der1 gmx_unused;

        evaluateFunctionAndDerivative<numFuncInTable, funcIndex0, funcIndex1>(
                r, functionValue0, &der0, functionValue1, &der1);
    }

    /*! \brief Evaluate function derivative only, two table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable   Number of separate functions in table, default is 2
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function derivative
     *  \param[out] derivativeValue0 Interpolated derivative for first function
     *  \param[out] derivativeValue1 Interpolated derivative for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateDerivative(T r, T* derivativeValue0, T* derivativeValue1) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex0, tabIndex, &Y, &F, &G, &H);
        *derivativeValue0 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex1, tabIndex, &Y, &F, &G, &H);
        *derivativeValue1 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /************************************************************
     *            Evaluation methods for three functions        *
     ************************************************************/


    /*! \brief Evaluate both function and derivative, three table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable   Number of separate functions in table, default is 3
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     funcIndex2       Index of 3rd function to evaluate in table, default is 2
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function and derivative
     *  \param[out] functionValue0   Interpolated value for first function
     *  \param[out] derivativeValue0 Interpolated derivative for first function
     *  \param[out] functionValue1   Interpolated value for second function
     *  \param[out] derivativeValue1 Interpolated derivative for second function
     *  \param[out] functionValue2   Interpolated value for third function
     *  \param[out] derivativeValue2 Interpolated derivative for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateFunctionAndDerivative(T  r,
                                       T* functionValue0,
                                       T* derivativeValue0,
                                       T* functionValue1,
                                       T* derivativeValue1,
                                       T* functionValue2,
                                       T* derivativeValue2) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable && funcIndex2 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex0, tabIndex, &Y, &F, &G, &H);
        *functionValue0   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue0 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex1, tabIndex, &Y, &F, &G, &H);
        *functionValue1   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue1 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex2, tabIndex, &Y, &F, &G, &H);
        *functionValue2   = fma(fma(fma(H, eps, G), eps, F), eps, Y);
        *derivativeValue2 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /*! \brief Evaluate function value only, three table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable   Number of separate functions in table, default is 3
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     funcIndex2       Index of 3rd function to evaluate in table, default is 2
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function value
     *  \param[out] functionValue0   Interpolated value for first function
     *  \param[out] functionValue1   Interpolated value for second function
     *  \param[out] functionValue2   Interpolated value for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateFunction(T r, T* functionValue0, T* functionValue1, T* functionValue2) const
    {
        T der0 gmx_unused;
        T der1 gmx_unused;
        T der2 gmx_unused;

        evaluateFunctionAndDerivative<numFuncInTable, funcIndex0, funcIndex1, funcIndex2>(
                r, functionValue0, &der0, functionValue1, &der1, functionValue2, &der2);
    }

    /*! \brief Evaluate function derivative only, three table functions
     *
     *  This is a templated method where the template can be either real or SimdReal.
     *
     *  \tparam     numFuncInTable   Number of separate functions in table, default is 3
     *  \tparam     funcIndex0       Index of 1st function to evaluate in table, default is 0
     *  \tparam     funcIndex1       Index of 2nd function to evaluate in table, default is 1
     *  \tparam     funcIndex2       Index of 3rd function to evaluate in table, default is 2
     *  \tparam     T                Type (SimdReal or real) of lookup and result
     *  \param      r                Points for which to evaluate function derivative
     *  \param[out] derivativeValue0 Interpolated derivative for first function
     *  \param[out] derivativeValue1 Interpolated derivative for second function
     *  \param[out] derivativeValue2 Interpolated derivative for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateDerivative(T r, T* derivativeValue0, T* derivativeValue1, T* derivativeValue2) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable && funcIndex2 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    Y, F, G, H;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex0, tabIndex, &Y, &F, &G, &H);
        *derivativeValue0 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex1, tabIndex, &Y, &F, &G, &H);
        *derivativeValue1 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                yfghMultiTableData_.data() + 4 * funcIndex2, tabIndex, &Y, &F, &G, &H);
        *derivativeValue2 = tableScale_ * fma(fma(T(3.0) * H, eps, T(2.0) * G), eps, F);
    }

    /*! \brief Return the table spacing (distance between points)
     *
     *  You should never have to use this for normal code, but due to the
     *  way tables are constructed internally we need this in the unit tests
     *  to check relative tolerances over each interval.
     *
     *  \return table spacing.
     */
    real tableSpacing() const { return 1.0 / tableScale_; }

private:
    std::size_t           numFuncInTable_; //!< Number of separate tabluated functions
    std::pair<real, real> range_;          //!< Range for which table evaluation is allowed
    real                  tableScale_;     //!< Table scale (inverse of spacing between points)

    /*! \brief Vector with combined table data to save calculations after lookup.
     *
     *  For table point i, this vector contains the four coefficients
     *  Y,F,G,H that we use to express the function value as
     *  V(x)  = Y + F e + G e^2 + H e^3, where e is the epsilon offset from
     *  the nearest table point.
     *
     *  To allow aligned SIMD loads we need to use an aligned allocator for
     *  this container.
     */
    std::vector<real, AlignedAllocator<real>> yfghMultiTableData_;

    // There should never be any reason to copy the table since it is read-only
    GMX_DISALLOW_COPY_AND_ASSIGN(CubicSplineTable);
};


} // namespace gmx

#endif // GMX_TABLES_CUBICSPLINETABLE_H
