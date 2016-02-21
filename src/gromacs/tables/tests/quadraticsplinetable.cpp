/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Tests for simple math functions.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include "gromacs/tables/quadraticsplinetable.h"

#include <cmath>
#include <cstdint>

#include <utility>

#include <gtest/gtest.h>

#include "gromacs/simd/simd.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

/*! \brief Test fixture for table comparision with analytical/numerical functions */
class QuadraticSplineTableTest : public ::testing::Test
{
    public:
        QuadraticSplineTableTest() : tolerance_(QuadraticSplineTable::defaultTolerance) {}

        /*! \brief Set a new tolerance to be used in table function comparison
         *
         *  \param tol New tolerance to use
         */
        void
        setTolerance(real tol) { tolerance_ = tol; }

        //! \cond internal
        /*! \internal \brief
         * Assertion predicate formatter for comparing table with function/derivative
         */
        ::testing::AssertionResult
        compareTableWithFunctions(const char *                          refFuncExpr,
                                  const char *                          refDerExpr,
                                  const char *                          tableExpr,
                                  const char gmx_unused *               testRangeExpr,
                                  const std::function<double(double)>  &refFunc,
                                  const std::function<double(double)>  &refDer,
                                  const QuadraticSplineTable           &table,
                                  const std::pair<real, real>           &testRange);
        //! \endcond

    private:
        static int  s_testPoints_; //!< Number of points (in range) to use
        real        tolerance_;    //!< Tolerance to use
};


int
QuadraticSplineTableTest::s_testPoints_ = 100;

/*! \brief Asserts that a table reproduces a function/derivative within tolerance.
 *
 * \hideinitializer
 */
#define GMX_TABLE_TEST_FUNC(refFunc, refDer, table, testRange) \
    EXPECT_PRED_FORMAT4(compareTableWithFunctions, refFunc, refDer, table, testRange)


::testing::AssertionResult
QuadraticSplineTableTest::compareTableWithFunctions(const char *                          refFuncExpr,
                                                    const char *                          refDerExpr,
                                                    const char *                          tableExpr,
                                                    const char gmx_unused *               testRangeExpr,
                                                    const std::function<double(double)>  &refFunc,
                                                    const std::function<double(double)>  &refDer,
                                                    const QuadraticSplineTable           &table,
                                                    const std::pair<real, real>           &testRange)
{
    real dx = (testRange.second - testRange.first) / s_testPoints_;

    for (real x = testRange.first; x < testRange.second; x += dx)
    {
        real refFuncValue     = refFunc(x);
        real refDerValue      = refDer(x);
        real h                = std::sqrt(GMX_REAL_EPS);
        real secondDerivative = (refDer(x+h)-refDer(x))/h;

        real testFuncValue;
        real testDerValue;

        table.evaluateFunctionAndDerivative(x, &testFuncValue, &testDerValue);

        // Interpolation using both Function and Derivative or only one of them should match exactly
        EXPECT_EQ(testFuncValue, table.evaluateFunction(x));
        EXPECT_EQ(testDerValue, table.evaluateDerivative(x));

        // Allowed *absolute* errors due to the loss-of-accuracy when calculating the index through rtable
        real                    allowedAbsFuncErr = 5.0*std::abs(secondDerivative) * x * GMX_REAL_EPS;
        real                    allowedAbsDerErr  = 5.0*std::abs(secondDerivative) * x * GMX_REAL_EPS;

        FloatingPointDifference funcDiff(refFuncValue, testFuncValue);
        FloatingPointDifference derDiff(refDerValue, testDerValue);

        FloatingPointTolerance  funcTol = relativeToleranceAsFloatingPoint(std::abs(refFuncValue), tolerance_);
        FloatingPointTolerance  derTol  = relativeToleranceAsFloatingPoint(std::abs(refDerValue), tolerance_);

        if ((funcDiff.asAbsolute() > allowedAbsFuncErr && !funcTol.isWithin(funcDiff)) ||
            (derDiff.asAbsolute() > allowedAbsDerErr && !derTol.isWithin(derDiff)))
        {
            return ::testing::AssertionFailure()
                   << "Failing comparison with function for table: " << tableExpr << std::endl
                   << "Reference function:   " << refFuncExpr << std::endl
                   << "Reference derivative: " << refDerExpr << std::endl
                   << "Tolerance:            " << tolerance_ << std::endl
                   << "Test range is ( " << testRange.first << " , " << testRange.second << " ) " << std::endl
                   << "First failure at    x = " << x << std::endl
                   << "Reference function    = " << refFuncValue << std::endl
                   << "Test table function   = " << testFuncValue << std::endl
                   << "Reference derivative  = " << refDerValue << std::endl
                   << "Test table derivative = " << testDerValue << std::endl;
        }
    }
    return ::testing::AssertionSuccess();
}


/*! \brief Function similar to power-12 Lennard-Jones repulsion
 *
 *  \param r argument
 *  \return r^-12
 */
double
lj12Function(double r)
{
    return std::pow(r, -12.0);
}

/*! \brief Derivative (not force) of the power-12 Lennard-Jones repulsion
 *
 *  \param r argument
 *  \return -12.0*r^-13
 */
double
lj12Derivative(double r)
{
    return -12.0*std::pow(r, -13.0);
}

/*! \brief The sinc function, sin(r)/r
 *
 *  \param r argument
 *  \return sin(r)/r
 */
double
sincFunction(double r)
{
    return std::sin(r)/r;
}

/*! \brief Derivative of the sinc function
 *
 *  \param r argument
 *  \return derivative of sinc, (r*cos(r)-sin(r))/r^2
 */
double
sincDerivative(double r)
{
    return (r*std::cos(r)-std::sin(r))/(r*r);
}

/*! \brief Function for the direct-space PME correction to 1/r
 *
 *  \param r argument
 *  \return PME correction function, erf(r)/r
 */
double
pmeCorrFunction(double r)
{
    return std::erf(r)/r;
}

/*! \brief Derivative of the direct-space PME correction to 1/r
 *
 *  \param r argument
 *  \return Derivative of the PME correction function.
 */
double
pmeCorrDerivative(double r)
{
    return (2.0*std::exp(-r*r)/std::sqrt(M_PI)*r-erf(r))/(r*r);
}


TEST_F(QuadraticSplineTableTest, IncorrectInput)
{
    // negative range
    EXPECT_THROW_GMX(QuadraticSplineTable(lj12Function, lj12Derivative, std::make_pair(-1.0, 0.0)), gmx::APIError);

    // Too small range
    EXPECT_THROW_GMX(QuadraticSplineTable(lj12Function, lj12Derivative, std::make_pair(1.0, 1.00001)), gmx::APIError);

    // bad tolerance
    EXPECT_THROW_GMX(QuadraticSplineTable(lj12Function, lj12Derivative, std::make_pair(1.0, 2.0), 1e-10), gmx::ToleranceError);

    // Range is so close to 0.0 that table would require >1e6 points
    EXPECT_THROW_GMX(QuadraticSplineTable(lj12Function, lj12Derivative, std::make_pair(1e-4, 2.0)), gmx::ToleranceError);

    // mismatching function/derivative
    EXPECT_THROW_GMX(QuadraticSplineTable(lj12Derivative, lj12Function, std::make_pair(1.0, 2.0)), gmx::APIError);
}

#ifndef NDEBUG
TEST_F(QuadraticSplineTableTest, CatchOutOfRangeValues)
{
    QuadraticSplineTable table(lj12Function, lj12Derivative, std::make_pair(0.1, 1.0));
    real                 x, func, der;

    x = 0.1*(1.0-GMX_REAL_EPS);
    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);

    x = 1.0;
    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);
}
#endif

TEST_F(QuadraticSplineTableTest, Sinc)
{
    std::pair<real, real>  tableRange(0.1, 10);

    QuadraticSplineTable   sincTable(sincFunction, sincDerivative, tableRange);

    GMX_TABLE_TEST_FUNC(sincFunction, sincDerivative, sincTable, tableRange);
}


TEST_F(QuadraticSplineTableTest, LJ12)
{
    std::pair<real, real>  tableRange(0.1, 2.0);

    QuadraticSplineTable   lj12Table(lj12Function, lj12Derivative, tableRange);

    GMX_TABLE_TEST_FUNC(lj12Function, lj12Derivative, lj12Table, tableRange);
}

TEST_F(QuadraticSplineTableTest, PmeCorrection)
{
    std::pair<real, real>  tableRange(0.15, 4.0);
    real                   tolerance = 1e-5;

    QuadraticSplineTable   pmeCorrTable(pmeCorrFunction, pmeCorrDerivative, tableRange, tolerance);

    setTolerance(tolerance);
    GMX_TABLE_TEST_FUNC(pmeCorrFunction, pmeCorrDerivative, pmeCorrTable, tableRange);
}


#ifdef GMX_SIMD_HAVE_REAL
TEST_F(QuadraticSplineTableTest, Simd)
{
    std::pair<real, real>  range(0.1, 1.0);
    QuadraticSplineTable   table(lj12Function, lj12Derivative, range);

    // We already test that the SIMD operations handle the different elements
    // correctly in the SIMD module, so here we only test that interpolation
    // works for a single value in the middle of the interval

    real     x       = 0.5 * (range.first + range.second);
    real     refFunc = lj12Function(x);
    real     refDer  = lj12Derivative(x);
    SimdReal tstFunc, tstDer;
    real     funcErr, derErr;
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) alignedMem[GMX_SIMD_REAL_WIDTH];

    table.evaluateFunctionAndDerivative(SimdReal(x), &tstFunc, &tstDer);

    store(alignedMem, tstFunc);
    funcErr = std::abs(alignedMem[0]-refFunc) / std::abs(refFunc);

    store(alignedMem, tstDer);
    derErr  = std::abs(alignedMem[0]-refDer ) / std::abs(refDer );

    EXPECT_LT(funcErr, QuadraticSplineTable::defaultTolerance);
    EXPECT_LT(derErr, QuadraticSplineTable::defaultTolerance);
}
#endif


#if defined GMX_SIMD_HAVE_REAL && !defined NDEBUG
TEST_F(QuadraticSplineTableTest, CatchOutOfRangeValuesSimd)
{
    std::pair<real, real>  range(0.1, 1.0);
    QuadraticSplineTable   table(lj12Function, lj12Derivative, range);
    SimdReal               x, func, der;

    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) alignedMem[GMX_SIMD_REAL_WIDTH];

    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        alignedMem[i] = range.first;
    }
    // Make position 1 incorrect if width>=2, otherwise position 0
    alignedMem[ (GMX_SIMD_REAL_WIDTH >= 2) ? 1 : 0] = range.first*(1.0-GMX_REAL_EPS);
    x = load(alignedMem);

    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);

    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        alignedMem[i] = range.second*(1.0-GMX_REAL_EPS);
    }
    // Make position 1 incorrect if width>=2, otherwise position 0
    alignedMem[ (GMX_SIMD_REAL_WIDTH >= 2) ? 1 : 0] = range.second;
    x = load(alignedMem);

    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);
}
#endif


TEST_F(QuadraticSplineTableTest, IncorrectNumericalInput)
{
    // Lengths do not match
    std::vector<double>   functionValues(10);
    std::vector<double>   derivativeValues(20);
    EXPECT_THROW_GMX(QuadraticSplineTable(functionValues, derivativeValues, 0.001, std::make_pair(1.0, 2.0)), gmx::APIError);

    // Upper range is 2.0, spacing 0.1. This requires at least 21 points. Make sure we get an error for 20.
    functionValues.resize(20);
    derivativeValues.resize(20);
    EXPECT_THROW_GMX(QuadraticSplineTable(functionValues, derivativeValues, 0.1, std::make_pair(1.0, 2.0)), gmx::APIError);

    // We need at least 5 points for the numerical derivative evaluation (test for constructor w/o derivative vector)
    functionValues.resize(4);
    EXPECT_THROW_GMX(QuadraticSplineTable(functionValues, 0.1, std::make_pair(0.1, 1.0)), gmx::APIError);

    // Create some test data
    functionValues.clear();
    derivativeValues.clear();

    std::vector<double>   badDerivativeValues;
    double                spacing = 1e-3;

    for (std::size_t i = 0; i < 1001; i++)
    {
        double x    = i * spacing;
        double func = (x >= 0.1) ? lj12Function(x)   : 0.0;
        double der  = (x >= 0.1) ? lj12Derivative(x) : 0.0;

        functionValues.push_back(func);
        derivativeValues.push_back(der);
        badDerivativeValues.push_back(-der);
    }

    // Derivatives not consistent with function
    EXPECT_THROW_GMX(QuadraticSplineTable(functionValues, badDerivativeValues, spacing, std::make_pair(0.1, 1.0)), gmx::APIError);

    // Spacing 1e-3 is not sufficient for r^-12 in range [0.1,1.0]
    // Make sure we get a tolerance error
    EXPECT_THROW_GMX(QuadraticSplineTable(functionValues, derivativeValues, spacing, std::make_pair(0.1, 1.0)), gmx::ToleranceError);
}


TEST_F(QuadraticSplineTableTest, NumericalInputBothVectorsPmeCorr)
{
    std::pair<real, real>  tableRange(0.15, 4.0);
    std::vector<double>    functionValues;
    std::vector<double>    derivativeValues;

    double                 inputSpacing = 1e-3;
    real                   tolerance    = 1e-5;

    // We only need data up to the argument 4.0, but add 1% margin
    for (std::size_t i = 0; i < tableRange.second*1.01/inputSpacing; i++)
    {
        double x    = i * inputSpacing;
        double func = (x >= 0.1) ? pmeCorrFunction(x)   : 0.0;
        double der  = (x >= 0.1) ? pmeCorrDerivative(x) : 0.0;

        functionValues.push_back(func);
        derivativeValues.push_back(der);
    }

    QuadraticSplineTable pmeCorrTable(functionValues, derivativeValues, inputSpacing, tableRange, tolerance);

    setTolerance(tolerance);
    GMX_TABLE_TEST_FUNC(pmeCorrFunction, pmeCorrDerivative, pmeCorrTable, tableRange);
}

TEST_F(QuadraticSplineTableTest, NumericalInputFunctionOnlyPmeCorr)
{
    std::pair<real, real>  tableRange(0.15, 4.0);
    std::vector<double>    functionValues;

    double                 inputSpacing = 1e-3;
    real                   tolerance    = 1e-5;

    // We only need data up to the argument 4.0, but add some margin
    for (std::size_t i = 0; i < tableRange.second*1.01/inputSpacing; i++)
    {
        double x    = i * inputSpacing;
        double func = (x >= 0.1) ? pmeCorrFunction(x) : 0.0;

        functionValues.push_back(func);
    }

    QuadraticSplineTable pmeCorrTable(functionValues, inputSpacing, tableRange, tolerance);

    setTolerance(tolerance);
    GMX_TABLE_TEST_FUNC(pmeCorrFunction, pmeCorrDerivative, pmeCorrTable, tableRange);
}

} // namespace

} // namespace test

} // namespace gmx
