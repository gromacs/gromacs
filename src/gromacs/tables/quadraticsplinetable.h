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

/*! \libinternal
 * \defgroup module_tables  Classes for table interpolation
 * \ingroup group_utilitymodules
 *
 * \brief Table interpolation from analytical or numerical input
 *
 * This module provides quadratic spline interpolation tables used
 * both for the nonbonded kernels and listed interactions.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 */


/*! \libinternal \file
 * \brief
 * Declares classes for quadratic spline table
 *
 * \inlibraryapi
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#ifndef GMX_TABLES_QUADRATICSPLINETABLE_H
#define GMX_TABLES_QUADRATICSPLINETABLE_H

#include <cstddef>

#include <functional>
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

/*! \libinternal \brief Quadratic spline interpolation table.
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
 * The table data is internally adjusted to guarantee that the interpolated
 * derivative is the true derivative of the interpolated potential, which is
 * important to avoid systematic errors for the common case when the derivative
 * is concave/convex in the entire interval.
 * We do this by expressing the difference in the function value
 * at a small offset h relative to a reference value in position 0 with a forward
 * Taylor series expanded around 0, and then doing the opposite of expressing
 * difference in the function at position 0 relative to a reference value in
 * position h when using a backward Taylor expansion:
 *
 * \f{eqnarray*}{
 *  \Delta V & = & hV'(0) + \frac{1}{2} h^2 V''(0) + \frac{1}{6} h^3 V'''(0) + O(h^4) \\
 *  \Delta V & = & hV'(h) - \frac{1}{2} h^2 V''(h) + \frac{1}{6} h^3 V'''(h) + O(h^4)
 * \f}
 *
 * Summing the equations leads to
 *
 * \f[
 *  2 \Delta V = h(V'(0) + V'(h)) + \frac{1}{2} h^2 (V''(0)-V''(h)) +
 * \frac{1}{6}h^3(V'''(0)+V'''(h)) + O(h^4) \f]
 *
 * To make the second term symmetric too, we can replace it with the average of
 * the Taylor expansion at 0 and h (i.e., using the third derivative). This gives
 *
 * \f[
 *  2 \Delta V = h(V'(0) + V'(h)) - \frac{1}{12} h^3 (V'''(0)+V'''(h)) + O(h^4)
 * \f]
 *
 * Thus, if we replace the derivative in the internal quadratic table data with
 *
 * \f[
 *  V' - \frac{1}{12}h^2 V'''
 * \f]
 *
 * we will cancel the h^3 term in the error. This will make the integral of the
 * forces match the potential much better (The h^4 term actually disappears, so
 * when summing over 1/h points the remaining error will be O(h^4).
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
 *       The absolute such error both in the function and derivative value will
 *       be roughly f''*x*GMX_REAL_EPS, where x is the argument and f'' the
 *       second derivative.
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
class QuadraticSplineTable
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
    /*! \brief Default tolerance for tables is 10*GMX_FLOAT_EPS
     *
     *  \note Even for double precision builds we set the tolerance to
     *        one order of magnitude above the single precision epsilon.
     */
    static LIBGROMACS_EXPORT const real defaultTolerance;

    /*! \brief Initialize table data from function
     *
     * \param analyticalInputList Initializer list with one or more functions to tabulate,
     *                            specified as pairs containing analytical
     *                            functions and their derivatives. The function will also
     *                            be called for values smaller than the lower limit of the
     *                            range, but we avoid calling it for 0.0 if that value
     *                            is not included in the range.
     * \param range               Range over which the function will be tabulated.
     *                            Constructor will throw gmx::APIError for negative values.
     *                            Due to the way the numerical derivative evaluation depends
     *                            on machine precision internally, this range must be
     *                            at least 0.001, or the constructor throws gmx::APIError.
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
    QuadraticSplineTable(std::initializer_list<AnalyticalSplineTableInput> analyticalInputList,
                         const std::pair<real, real>&                      range,
                         real tolerance = defaultTolerance);

    /*! \brief Initialize table data from tabulated values and derivatives
     *
     * \param numericalInputList  Initializer list with one or more functions to tabulate,
     *                            specified as pairs containing containing vectors for the
     *                            function values and their derivatives. Data points are
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
     * \param tolerance           Requested accuracy of the table in the range. This will be
     *                            used to calculate the required internal spacing and possibly
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
    QuadraticSplineTable(std::initializer_list<NumericalSplineTableInput> numericalInputList,
                         const std::pair<real, real>&                     range,
                         real tolerance = defaultTolerance);


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
        T    t0;
        T    t1;
        T    t2;
        T t3 gmx_unused;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex, tabIndex, &t0, &t1, &t2, &t3);

        t1               = t0 + eps * t1;
        *functionValue   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue = t1;
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
        T    t0;
        T    t1;
        T t2 gmx_unused;

        if (numFuncInTable == 1)
        {
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex, tabIndex, &t0, &t1); // works for scalar T too
        }
        else
        {
            // This is not ideal, but we need a version of gatherLoadUBySimdIntTranspose that
            // only loads a single value from memory to implement it better (will be written)
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
        }

        // (1-eps)*t0 + eps*t1
        *derivativeValue = fma(t1 - t0, eps, t0);
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
     *  \param[out] functionValue1   Interpolated value for first function
     *  \param[out] derivativeValue1 Interpolated derivative for first function
     *  \param[out] functionValue2   Interpolated value for second function
     *  \param[out] derivativeValue2 Interpolated derivative for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateFunctionAndDerivative(T r, T* functionValue1, T* derivativeValue1, T* functionValue2, T* derivativeValue2) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    t0;
        T    t1;
        T    t2;
        T t3 gmx_unused;

        // Load Derivative, Delta, Function, and Zero values for each table point.
        // The 4 refers to these four values - not any SIMD width.
        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex0, tabIndex, &t0, &t1, &t2, &t3);
        t1                = t0 + eps * t1;
        *functionValue1   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue1 = t1;

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex1, tabIndex, &t0, &t1, &t2, &t3);
        t1                = t0 + eps * t1;
        *functionValue2   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue2 = t1;
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
     *  \param[out] functionValue1   Interpolated value for first function
     *  \param[out] functionValue2   Interpolated value for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateFunction(T r, T* functionValue1, T* functionValue2) const
    {
        T der1 gmx_unused;
        T der2 gmx_unused;

        evaluateFunctionAndDerivative<numFuncInTable, funcIndex0, funcIndex1>(
                r, functionValue1, &der1, functionValue2, &der2);
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
     *  \param[out] derivativeValue1 Interpolated derivative for first function
     *  \param[out] derivativeValue2 Interpolated derivative for second function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 2, int funcIndex0 = 0, int funcIndex1 = 1, typename T>
    void evaluateDerivative(T r, T* derivativeValue1, T* derivativeValue2) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);

        if (numFuncInTable == 2 && funcIndex0 == 0 && funcIndex1 == 1)
        {
            T t0A, t0B, t1A, t1B;
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data(), tabIndex, &t0A, &t0B); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + 2, tabIndex, &t1A, &t1B); // works for scalar T too
            *derivativeValue1 = fma(t1A - t0A, eps, t0A);
            *derivativeValue2 = fma(t1B - t0B, eps, t0B);
        }
        else
        {
            T t0, t1, t2;
            // This is not ideal, but we need a version of gatherLoadUBySimdIntTranspose that
            // only loads a single value from memory to implement it better (will be written)
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex0, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex0,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
            *derivativeValue1 = fma(t1 - t0, eps, t0);

            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex1, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex1,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
            *derivativeValue2 = fma(t1 - t0, eps, t0);
        }
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
     *  \param[out] functionValue1   Interpolated value for first function
     *  \param[out] derivativeValue1 Interpolated derivative for first function
     *  \param[out] functionValue2   Interpolated value for second function
     *  \param[out] derivativeValue2 Interpolated derivative for second function
     *  \param[out] functionValue3   Interpolated value for third function
     *  \param[out] derivativeValue3 Interpolated derivative for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateFunctionAndDerivative(T  r,
                                       T* functionValue1,
                                       T* derivativeValue1,
                                       T* functionValue2,
                                       T* derivativeValue2,
                                       T* functionValue3,
                                       T* derivativeValue3) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable && funcIndex2 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);
        T    t0;
        T    t1;
        T    t2;
        T t3 gmx_unused;

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex0, tabIndex, &t0, &t1, &t2, &t3);
        t1                = t0 + eps * t1;
        *functionValue1   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue1 = t1;

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex1, tabIndex, &t0, &t1, &t2, &t3);
        t1                = t0 + eps * t1;
        *functionValue2   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue2 = t1;

        gatherLoadBySimdIntTranspose<4 * numFuncInTable>(
                ddfzMultiTableData_.data() + 4 * funcIndex2, tabIndex, &t0, &t1, &t2, &t3);
        t1                = t0 + eps * t1;
        *functionValue3   = fma(eps * T(halfSpacing_), t0 + t1, t2);
        *derivativeValue3 = t1;
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
     *  \param[out] functionValue1   Interpolated value for first function
     *  \param[out] functionValue2   Interpolated value for second function
     *  \param[out] functionValue3   Interpolated value for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateFunction(T r, T* functionValue1, T* functionValue2, T* functionValue3) const
    {
        T der1 gmx_unused;
        T der2 gmx_unused;
        T der3 gmx_unused;

        evaluateFunctionAndDerivative<numFuncInTable, funcIndex0, funcIndex1, funcIndex2>(
                r, functionValue1, &der1, functionValue2, &der2, functionValue3, &der3);
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
     *  \param[out] derivativeValue1 Interpolated derivative for first function
     *  \param[out] derivativeValue2 Interpolated derivative for second function
     *  \param[out] derivativeValue3 Interpolated derivative for third function
     *
     *  For debug builds we assert that the input values fall in the range
     *  specified when constructing the table.
     */
    template<int numFuncInTable = 3, int funcIndex0 = 0, int funcIndex1 = 1, int funcIndex2 = 2, typename T>
    void evaluateDerivative(T r, T* derivativeValue1, T* derivativeValue2, T* derivativeValue3) const
    {
        rangeCheck(r);
        GMX_ASSERT(numFuncInTable == numFuncInTable_,
                   "Evaluation method not matching number of functions in table");
        GMX_ASSERT(funcIndex0 < numFuncInTable && funcIndex1 < numFuncInTable && funcIndex2 < numFuncInTable,
                   "Function index not in range of the number of tables");

        T    rTable   = r * T(tableScale_);
        auto tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
        T    eps      = rTable - trunc(rTable);

        if (numFuncInTable == 3 && funcIndex0 == 0 && funcIndex1 == 1 && funcIndex2 == 2)
        {
            T t0A, t0B, t0C, t1A, t1B, t1C;
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data(), tabIndex, &t0A, &t0B);
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + 2, tabIndex, &t0C, &t1A);
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + 4, tabIndex, &t1B, &t1C);
            *derivativeValue1 = fma(t1A - t0A, eps, t0A);
            *derivativeValue2 = fma(t1B - t0B, eps, t0B);
            *derivativeValue3 = fma(t1C - t0C, eps, t0C);
        }
        else
        {
            T t0, t1, t2;
            // This is not ideal, but we need a version of gatherLoadUBySimdIntTranspose that
            // only loads a single value from memory to implement it better (will be written)
            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex0, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex0,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
            *derivativeValue1 = fma(t1 - t0, eps, t0);

            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex1, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex1,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
            *derivativeValue2 = fma(t1 - t0, eps, t0);

            gatherLoadUBySimdIntTranspose<numFuncInTable>(
                    derivativeMultiTableData_.data() + funcIndex2, tabIndex, &t0, &t2); // works for scalar T too
            gatherLoadUBySimdIntTranspose<numFuncInTable>(derivativeMultiTableData_.data() + funcIndex2,
                                                          tabIndex + T(1),
                                                          &t1,
                                                          &t2); // works for scalar T too
            *derivativeValue3 = fma(t1 - t0, eps, t0);
        }
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
    real                  halfSpacing_;    //!< 0.5*spacing (used for DDFZ table data)

    //!< Derivative values only, with the third-derivative subtraction described in the class documentation.
    std::vector<real> derivativeMultiTableData_;

    /*! \brief Combined derivative, difference to next derivative, value, and zero.
     *
     *  For table point i, this vector contains the four values:
     *  - derivative[i]
     *  - (derivative[i+1]-derivative[i])
     *  - value[i]
     *  - 0.0
     *
     *  For the derivative terms we have subtracted the third-derivative term described
     *  in the main class documentation.
     *
     *  This is typically more efficient than the individual tables, in particular
     *  when using SIMD. The result should be identical outside the class, so this
     *  is merely an internal implementation optimization. However, to allow
     *  aligned SIMD loads we need to use an aligned allocator for this container.
     *  We occasionally abbreviate this data as DDFZ.
     */
    std::vector<real, AlignedAllocator<real>> ddfzMultiTableData_;

    // There should never be any reason to copy the table since it is read-only
    GMX_DISALLOW_COPY_AND_ASSIGN(QuadraticSplineTable);
};


} // namespace gmx

#endif // GMX_TABLES_QUADRATICSPLINETABLE_H
