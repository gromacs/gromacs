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
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#ifndef GMX_TABLES_QUADRATICSPLINETABLE_H
#define GMX_TABLES_QUADRATICSPLINETABLE_H

#include <cstdint>

#include <functional>
#include <vector>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \libinternal \brief Quadratic spline interpolation table.
 *
 * This class interpolates a function specified either as an analytical
 * expression or from user-provided table data. The internal implementation is
 * optimized to guarantee that the interpolated derivative is the true derivative
 * of the interpolated potential, which is important to avoid systematic errors
 * for the common case when the derivative is concave/convex in the entire
 * interval. We do this by expressing the difference in function values
 * between two table locations as forward/backward Taylor series:
 *
 * \f{eqnarray*}{
 *  \Delta V & = & hV'(0) + \frac{1}{2} h^2 V''(0) + \frac{1}{6} h^3 V'''(0) + O(h^4) \\
 *  \Delta V & = & hV'(h) - \frac{1}{2} h^2 V''(0) + \frac{1}{6} h^3 V'''(0) + O(h^4)
 * \f}
 *
 * Summing the equations leads to
 *
 * \f[
 *  2 \Delta V = h(V'(0) + V'(h)) + \frac{1}{2} h^2 (V''(0)-V''(h)) + \frac{1}{6}h^3(V'''(0)+V'''(h)) + O(h^4)
 * \f]
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
 */
class QuadraticSplineTable
{
    private:
        template <typename T>
        void
        rangeCheck(T r) const
        {
#ifndef NDEBUG
            // Check that all values fall in range when debugging
            if (anyTrue( r < T(range_.first) || T(range_.second) <= r ) )
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
        static const real
            defaultTolerance;

        /*! \brief Initialize table data from function
         *
         * \param  function   Function to be tabulated in the interval specified by range.
         *                    To avoid floating-point exceptions the function must be
         *                    defined and produce values in the range
         *                    requested (including the endpoints) with absolute values
         *                    smaller than sqrt(GMX_DOUBLE_EPS).
         * \param  derivative Derivative of function. Same rules as for the
         *                    function object apply. We use an explicit analytical derivative
         *                    to avoid problems with numerical differentiation of
         *                    arbitrary functions, which would be prone to loss of accuracy.
         * \param  range      Range over which the function will be tabulated.
         *                    Constructor will throw gmx::APIError for negative values.
         *                    Due to the way the numerical derivative evaluation depends
         *                    on machine precision internally, this range must be
         *                    at least 0.001, or the constructor throws gmx::APIError.
         * \param  tolerance  Requested accuracy of the table. This will be used to
         *                    calculate the required internal spacing. If this cannot
         *                    be achieved (for instance because the table would require
         *                    too much memory) the constructor will throw gmx::ToleranceError.
         *
         * \note The functions are always defined in double precision to avoid
         *       losing accuracy when constructing tables.
         *
         * \note For efficiency reasons (since this code is used in some inner
         *       (kernels), we always allocate memory and calculate table indices
         *       for the complete interval [0,range.second], although the data will
         *       not be valid outside the definition range to avoid calling the
         *       function there. However, this means you should not use this class
         *       to tabulate functions for small ranges very far away from zero,
         *       since you would waste a huge amount of memory.
         *
         * \throws gmx::ToleranceError if the requested tolerance cannot be achieved,
         *         and gmx::APIError for other incorrect input.
         */
        QuadraticSplineTable(const std::function<double(double)> &function,
                             const std::function<double(double)> &derivative,
                             const std::pair<real, real>          &range,
                             real                                 tolerance = defaultTolerance);

        /*! \brief Initialize table data from tabulated values and derivatives
         *
         * \param  function     Vector with function values. Data points are
         *                      separated by the spacing parameter, starting from 0.
         *                      Values outside the range will not be used, with the
         *                      exception that some extra data points beyond range.second
         *                      are required to re-interpolate values, so add some margin.
         *                      The constructor will throw gmx::APIError if the input
         *                      vectors are too short to cover the requested range
         *                      (and they must always be at least five points).
         * \param  derivative   Vector containing function derivative values,
         *                      for the same points as the function value vector. Same
         *                      rules apply as for the value vector.
         * \param  inputSpacing Distance between input data points. We assume the first
         *                      point corresponds to 0.0.
         * \param  range        Range over which the function will be tabulated.
         *                      Constructor will throw gmx::APIError for negative values,
         *                      or if the value/derivative vector does not cover the
         *                      range.
         * \param  tolerance    Requested accuracy of the table. This will be used to
         *                      calculate the required internal spacing and possibly
         *                      re-interpolate. The constructor will throw
         *                      gmx::ToleranceError if the input spacing is too coarse
         *                      to achieve this accuracy.
         *
         * \note The input data vectors are always double precision to avoid
         *       losing accuracy when constructing tables.
         *
         * \note For efficiency reasons (since this code is used in some inner
         *       (kernels), we always allocate memory and calculate table indices
         *       for the complete interval [0,range.second], although the data will
         *       not be valid outside the definition range to avoid calling the
         *       function there. However, this means you should not use this class
         *       to tabulate functions for small ranges very far away from zero,
         *       since you would waste a huge amount of memory.
         */
        QuadraticSplineTable(const std::vector<double>     &function,
                             const std::vector<double>     &derivative,
                             double                         inputSpacing,
                             const std::pair<real, real>    &range,
                             real                           tolerance = defaultTolerance);

        /*! \brief Initialize table data from tabulated values only (no derivative)
         *
         * \param  function     Vector with function values. Data points are
         *                      separated by the spacing parameter, starting from 0.
         *                      Values outside the range will not be used, with the
         *                      exception that some extra data points beyond range.second
         *                      are required to re-interpolate values, so add some margin.
         *                      The constructor will throw gmx::APIError if the input
         *                      vector is too short to cover the requested range
         *                      (and it must always be at least five points).
         * \param  inputSpacing Distance between input data points. We assume the first
         *                      point corresponds to 0.0.
         * \param  range        Range over which the function will be tabulated.
         *                      Constructor will throw gmx::APIError for negative values,
         *                      or if the value/derivative vector does not cover the
         *                      range. Note that we might need the value/derivative
         *                      at the endpoint of the range to construct the table.
         * \param  tolerance    Requested accuracy of the table. This will be used to
         *                      calculate the required internal spacing and possibly
         *                      re-interpolate. The constructor will throw
         *                      gmx::ToleranceError if the input spacing is too coarse
         *                      to achieve this accuracy.
         *
         * \note The input data vector is always double precision to avoid
         *       losing accuracy when constructing tables.
         *
         * \note For efficiency reasons (since this code is used in some inner
         *       (kernels), we always allocate memory and calculate table indices
         *       for the complete interval [0,range.second], although the data will
         *       not be valid outside the definition range to avoid calling the
         *       function there. However, this means you should not use this class
         *       to tabulate functions for small ranges very far away from zero,
         *       since you would waste a huge amount of memory.
         */
        QuadraticSplineTable(const std::vector<double>    &function,
                             double                        inputSpacing,
                             const std::pair<real, real>   &range,
                             real                          tolerance = defaultTolerance);

        /*! \brief Evaluate both function and derivative values from table
         *
         *  This is a templated method where the template can be either real or SimdReal.
         *
         *  \param      r               Points for which to evaluate function and derivative
         *  \param[out] functionValue   Function values
         *  \param[out] derivativeValue Function derivatives
         *
         *  For debug builds we assert that the input values fall in the range
         *  specified when constructing the table.
         */
        template <typename T>
        void
        evaluateFunctionAndDerivative(T     r,
                                      T *   functionValue,
                                      T *   derivativeValue) const
        {
            rangeCheck(r);

            T     rTable   = r * T(tableScale_);
            auto  tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
            T     eps      = rTable - trunc(rTable);
            T     t0;
            T     t1;
            T     t2;
            T     t3 gmx_unused;

            // Load Derivative, Delta, Function, and Zero values for each table point.
            // The 4 refers to these four values - not any SIMD width.
            gatherLoadBySimdIntTranspose<4>(derivativeDeltaFunctionZeroData_.data(), tabIndex, &t0, &t1, &t2, &t3);

            t1               = t0 + eps * t1;
            *functionValue   = fma( eps * T(halfSpacing_),  t0 + t1, t2);
            *derivativeValue = t1;
        }

        /*! \brief Evaluate function value from table
         *
         *  This is a templated method where the template can be either real or SimdReal.
         *
         *  \param      r          Points for which to evaluate function value
         *  \return Function value
         *
         *  For debug builds we assert that the input values fall in the range
         *  specified when constructing the table.
         */
        template <typename T>
        T
        evaluateFunction(T r) const
        {
            T     func;
            T     der gmx_unused;

            evaluateFunctionAndDerivative(r, &func, &der);

            return func;
        }

        /*! \brief Evaluate function derivative from table
         *
         *  This is a templated method where the template can be either real or SimdReal.
         *
         *  \param      r          Points for which to evaluate function derivative
         *  \return Function derivative
         *
         *  For debug builds we assert that the input values fall in the range
         *  specified when constructing the table.
         */
        template <typename T>
        T
        evaluateDerivative(T r) const
        {
            rangeCheck(r);

            T     rTable   = r * T(tableScale_);
            auto  tabIndex = cvttR2I(rTable); // type is either std::int32_t or SimdInt32
            T     eps      = rTable - trunc(rTable);
            T     t0;
            T     t1;

            gatherLoadUBySimdIntTranspose<1>(derivativeData_.data(), tabIndex, &t0, &t1); // works for scalar T too

            // (1-eps)*t0 + eps*t1
            return fma(t1-t0, eps, t0);
        }

    private:

        std::pair<real, real>   range_;       //!< Range for which table evaluation is allowed
        real                    tableScale_;  //!< Table scale (inverse of spacing between points)
        real                    halfSpacing_; //!< 0.5*spacing (used for DDFZ table data)

        //!< Function values only.
        std::vector<real>      functionData_;
        //!< Derivative values only, with the third-derivative subtraction described in the class documentation.
        std::vector<real>      derivativeData_;

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
        std::vector<real, AlignedAllocator<real> >  derivativeDeltaFunctionZeroData_;
};

}      // namespace gmx

#endif // GMX_TABLES_QUADRATICSPLINETABLE_H
