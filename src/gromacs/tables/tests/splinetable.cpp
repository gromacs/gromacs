/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 * Tests for simple math functions.eval
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */
#include "gmxpre.h"

#include <cmath>

#include <algorithm>
#include <functional>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/simd/simd.h"
#include "gromacs/tables/cubicsplinetable.h"
#include "gromacs/tables/quadraticsplinetable.h"

#include "testutils/testasserts.h"
#include "testutils/testoptions.h"


namespace gmx
{

namespace test
{

namespace
{

class SplineTableTestBase : public ::testing::Test
{
    public:
        static int  s_testPoints_; //!< Number of points to use. Public so we can set it as option
};

int
SplineTableTestBase::s_testPoints_ = 100;

/*! \cond */
/*! \brief Command-line option to adjust the number of points used to test SIMD math functions. */
GMX_TEST_OPTIONS(SplineTableTestOptions, options)
{
    options->addOption(::gmx::IntegerOption("npoints")
                           .store(&SplineTableTestBase::s_testPoints_)
                           .description("Number of points to test for spline table functions"));
}
/*! \endcond */




/*! \brief Test fixture for table comparision with analytical/numerical functions */
template <typename T>
class SplineTableTest : public SplineTableTestBase
{
    public:
        SplineTableTest() : tolerance_(T::defaultTolerance) {}

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
        template<int numFuncInTable = 1, int funcIndex = 0>
        void
        testSplineTableAgainstFunctions(const std::string                    &desc,
                                        const std::function<double(double)>  &refFunc,
                                        const std::function<double(double)>  &refDer,
                                        const T                              &table,
                                        const std::pair<real, real>          &testRange);
        //! \endcond

    private:
        real        tolerance_;    //!< Tolerance to use
};

template <class T>
template<int numFuncInTable, int funcIndex>
void
SplineTableTest<T>::testSplineTableAgainstFunctions(const std::string                    &desc,
                                                    const std::function<double(double)>  &refFunc,
                                                    const std::function<double(double)>  &refDer,
                                                    const T                              &table,
                                                    const std::pair<real, real>          &testRange)
{
    real                   dx = (testRange.second - testRange.first) / s_testPoints_;

    FloatingPointTolerance funcTolerance(relativeToleranceAsFloatingPoint(0.0, tolerance_));

    for (real x = testRange.first; x < testRange.second; x += dx)
    {
        real h                = std::sqrt(GMX_REAL_EPS);
        real secondDerivative = (refDer(x+h)-refDer(x))/h;

        real testFuncValue;
        real testDerValue;

        table.template evaluateFunctionAndDerivative<numFuncInTable, funcIndex>(x, &testFuncValue, &testDerValue);

        // Check that we get the same values from function/derivative-only methods
        real tmpFunc, tmpDer;

        table.template evaluateFunction<numFuncInTable, funcIndex>(x, &tmpFunc);

        table.template evaluateDerivative<numFuncInTable, funcIndex>(x, &tmpDer);

        if (testFuncValue != tmpFunc)
        {
            ADD_FAILURE()
            << "Interpolation inconsistency for table " << desc << std::endl
            << numFuncInTable << " function(s) in table, testing index " << funcIndex << std::endl
            << "First failure at x = " << x << std::endl
            << "Function value when evaluating function & derivative: " << testFuncValue << std::endl
            << "Function value when evaluating only function:         " << tmpFunc << std::endl;
            return;
        }
        if (testDerValue != tmpDer)
        {
            ADD_FAILURE()
            << "Interpolation inconsistency for table " << desc << std::endl
            << numFuncInTable << " function(s) in table, testing index " << funcIndex << std::endl
            << "First failure at x = " << x << std::endl
            << "Derivative value when evaluating function & derivative: " << testDerValue << std::endl
            << "Derivative value when evaluating only derivative:       " << tmpDer << std::endl;
            return;
        }

        // There are two sources of errors that we need to account for when checking the values,
        // and we only fail the test if both of these tolerances are violated:
        //
        // 1) First, we have the normal relative error of the test vs. reference value. For this
        //    we use the normal testutils relative tolerance checking.
        //    However, there is an additional source of error: When we calculate the forces we
        //    use average higher derivatives over the interval to improve accuracy, but this
        //    also means we won't reproduce values at table points exactly. This is usually not
        //    an issue since the tolerances we have are much larger, but when the reference derivative
        //    value is exactly zero the relative error will be infinite. To account for this, we
        //    extract the spacing from the table and evaluate the reference derivative at a point
        //    this much larger too, and use the largest of the two values as the reference
        //    magnitude for the derivative when setting the relative tolerance.
        //    Note that according to the table function definitions, we should be allowed to evaluate
        //    it one table point beyond the range (this is done already for construction).
        //
        // 2) Second, due to the loss-of-accuracy when calculating the index through rtable
        //    there is an internal absolute tolerance that we can calculate.
        //    The source of this error is the subtraction eps=rtab-[rtab], which leaves an
        //    error proportional to eps_machine*rtab=eps_machine*x*tableScale.
        //    To lowest order, the term in the function and derivative values (respectively) that
        //    are proportional to eps will be the next-higher derivative multiplied by the spacing.
        //    This means the truncation error in the value is derivative*x*eps_machine, and in the
        //    derivative the error is 2nd_derivative*x*eps_machine.

        real refFuncValue     = refFunc(x);
        real refDerValue      = refDer(x);
        real nextRefDerValue  = refDer(x + table.tableSpacing());

        real derMagnitude     = std::max( std::abs(refDerValue), std::abs(nextRefDerValue));

        // Since the reference magnitude will change over each interval we need to re-evaluate
        // the derivative tolerance inside the loop.
        FloatingPointTolerance  derTolerance(relativeToleranceAsFloatingPoint(derMagnitude, tolerance_));

        FloatingPointDifference funcDiff(refFuncValue, testFuncValue);
        FloatingPointDifference derDiff(refDerValue, testDerValue);

        real                    allowedAbsFuncErr = std::abs(refDerValue)      * x * GMX_REAL_EPS;
        real                    allowedAbsDerErr  = std::abs(secondDerivative) * x * GMX_REAL_EPS;

        if ((!funcTolerance.isWithin(funcDiff) && funcDiff.asAbsolute() > allowedAbsFuncErr) ||
            (!derTolerance.isWithin(derDiff)  &&  derDiff.asAbsolute() > allowedAbsDerErr))
        {
            ADD_FAILURE()
            << "Failing comparison with function for table " << desc << std::endl
            << numFuncInTable << " function(s) in table, testing index " << funcIndex << std::endl
            << "Test range is ( " << testRange.first << " , " << testRange.second << " ) " << std::endl
            << "Tolerance             = " << tolerance_ << std::endl
            << "First failure at    x = " << x << std::endl
            << "Reference function    = " << refFuncValue << std::endl
            << "Test table function   = " << testFuncValue << std::endl
            << "Reference derivative  = " << refDerValue << std::endl
            << "Test table derivative = " << testDerValue << std::endl;
            return;
        }
    }
}


/*! \brief Function similar to coulomb electrostatics
 *
 *  \param r argument
 *  \return r^-1
 */
double
coulombFunction(double r)
{
    return 1.0/r;
}

/*! \brief Derivative (not force) of coulomb electrostatics
 *
 *  \param r argument
 *  \return -r^-2
 */
double
coulombDerivative(double r)
{
    return -1.0/(r*r);
}

/*! \brief Function similar to power-6 Lennard-Jones dispersion
 *
 *  \param r argument
 *  \return r^-6
 */
double
lj6Function(double r)
{
    return std::pow(r, -6.0);
}

/*! \brief Derivative (not force) of the power-6 Lennard-Jones dispersion
 *
 *  \param r argument
 *  \return -6.0*r^-7
 */
double
lj6Derivative(double r)
{
    return -6.0*std::pow(r, -7.0);
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
    if (r == 0)
    {
        return 2.0/std::sqrt(M_PI);
    }
    else
    {
        return std::erf(r)/r;
    }
}

/*! \brief Derivative of the direct-space PME correction to 1/r
 *
 *  \param r argument
 *  \return Derivative of the PME correction function.
 */
double
pmeCorrDerivative(double r)
{
    if (r == 0)
    {
        return 0;
    }
    else
    {
        return (2.0*std::exp(-r*r)/std::sqrt(3.14159265358979323846)*r-erf(r))/(r*r);
    }
}

/*! \brief Typed-test list. We test QuadraticSplineTable and CubicSplineTable
 */
typedef ::testing::Types<QuadraticSplineTable, CubicSplineTable> SplineTableTypes;
TYPED_TEST_CASE(SplineTableTest, SplineTableTypes);


TYPED_TEST(SplineTableTest, HandlesIncorrectInput)
{
    // negative range
    EXPECT_THROW_GMX(TypeParam( {{"LJ12", lj12Function, lj12Derivative}}, {-1.0, 0.0}), gmx::InvalidInputError);

    // Too small range
    EXPECT_THROW_GMX(TypeParam( {{"LJ12", lj12Function, lj12Derivative}}, {1.0, 1.00001}), gmx::InvalidInputError);

    // bad tolerance
    EXPECT_THROW_GMX(TypeParam( {{"LJ12", lj12Function, lj12Derivative}}, {1.0, 2.0}, 1e-20), gmx::ToleranceError);

    // Range is so close to 0.0 that table would require >1e6 points
    EXPECT_THROW_GMX(TypeParam( {{"LJ12", lj12Function, lj12Derivative}}, {1e-4, 2.0}), gmx::ToleranceError);

    // mismatching function/derivative
    EXPECT_THROW_GMX(TypeParam( { {"BadLJ12", lj12Derivative, lj12Function}}, {1.0, 2.0}), gmx::InconsistentInputError);
}


#ifndef NDEBUG
TYPED_TEST(SplineTableTest, CatchesOutOfRangeValues)
{
    TypeParam      table( {{"LJ12", lj12Function, lj12Derivative}}, {0.2, 1.0});
    real           x, func, der;

    x = -GMX_REAL_EPS;
    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);

    x = 1.0;
    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);
}
#endif


TYPED_TEST(SplineTableTest, Sinc)
{
    std::pair<real, real>  range(0.1, 10);

    TypeParam              sincTable( {{"Sinc", sincFunction, sincDerivative}}, range);

    TestFixture::testSplineTableAgainstFunctions("Sinc", sincFunction, sincDerivative, sincTable, range);
}


TYPED_TEST(SplineTableTest, LJ12)
{
    std::pair<real, real>  range(0.2, 2.0);

    TypeParam              lj12Table( {{"LJ12", lj12Function, lj12Derivative}}, range);

    TestFixture::testSplineTableAgainstFunctions("LJ12", lj12Function, lj12Derivative, lj12Table, range);
}


TYPED_TEST(SplineTableTest, PmeCorrection)
{
    std::pair<real, real>  range(0.0, 4.0);
    real                   tolerance = 1e-5;

    TypeParam              pmeCorrTable( {{"PMECorr", pmeCorrFunction, pmeCorrDerivative}}, range, tolerance);

    TestFixture::setTolerance(tolerance);
    TestFixture::testSplineTableAgainstFunctions("PMECorr", pmeCorrFunction, pmeCorrDerivative, pmeCorrTable, range);
}



TYPED_TEST(SplineTableTest, HandlesIncorrectNumericalInput)
{
    // Lengths do not match
    std::vector<double>   functionValues(10);
    std::vector<double>   derivativeValues(20);
    EXPECT_THROW_GMX(TypeParam( {{"EmptyVectors", functionValues, derivativeValues, 0.001}},
                                {1.0, 2.0}), gmx::InconsistentInputError);

    // Upper range is 2.0, spacing 0.1. This requires at least 21 points. Make sure we get an error for 20.
    functionValues.resize(20);
    derivativeValues.resize(20);
    EXPECT_THROW_GMX(TypeParam( {{"EmptyVectors", functionValues, derivativeValues, 0.1}},
                                {1.0, 2.0}), gmx::InconsistentInputError);

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
    EXPECT_THROW_GMX(TypeParam( {{"NumericalBadLJ12", functionValues, badDerivativeValues, spacing}},
                                {0.2, 1.0}), gmx::InconsistentInputError);

    // Spacing 1e-3 is not sufficient for r^-12 in range [0.1,1.0]
    // Make sure we get a tolerance error
    EXPECT_THROW_GMX(TypeParam( {{"NumericalLJ12", functionValues, derivativeValues, spacing}},
                                {0.2, 1.0}), gmx::ToleranceError);
}


TYPED_TEST(SplineTableTest, NumericalInputPmeCorr)
{
    std::pair<real, real>  range(0.0, 4.0);
    std::vector<double>    functionValues;
    std::vector<double>    derivativeValues;

    double                 inputSpacing = 1e-3;
    real                   tolerance    = 1e-5;

    // We only need data up to the argument 4.0, but add 1% margin
    for (std::size_t i = 0; i < range.second*1.01/inputSpacing; i++)
    {
        double x    = i * inputSpacing;

        functionValues.push_back(pmeCorrFunction(x));
        derivativeValues.push_back(pmeCorrDerivative(x));
    }

    TypeParam  pmeCorrTable( {{"NumericalPMECorr", functionValues, derivativeValues, inputSpacing}},
                             range, tolerance);

    TestFixture::setTolerance(tolerance);
    TestFixture::testSplineTableAgainstFunctions("NumericalPMECorr", pmeCorrFunction, pmeCorrDerivative, pmeCorrTable, range);
}

TYPED_TEST(SplineTableTest, TwoFunctions)
{
    std::pair<real, real>  range(0.2, 2.0);

    TypeParam              table( {{"LJ6", lj6Function, lj6Derivative}, {"LJ12", lj12Function, lj12Derivative}}, range);

    // Test entire range for each function. This will use the method that interpolates a single function
    TestFixture::template testSplineTableAgainstFunctions<2, 0>("LJ6", lj6Function, lj6Derivative, table, range);
    TestFixture::template testSplineTableAgainstFunctions<2, 1>("LJ12", lj12Function, lj12Derivative, table, range);

    // Test the methods that evaluated both functions for one value
    real     x        = 0.5 * (range.first + range.second);
    real     refFunc0 = lj6Function(x);
    real     refDer0  = lj6Derivative(x);
    real     refFunc1 = lj12Function(x);
    real     refDer1  = lj12Derivative(x);

    real     tstFunc0, tstDer0, tstFunc1, tstDer1;
    real     tmpFunc0, tmpFunc1, tmpDer0, tmpDer1;

    // test that we reproduce the reference functions
    table.evaluateFunctionAndDerivative(x, &tstFunc0, &tstDer0, &tstFunc1, &tstDer1);

    real funcErr0 = std::abs(tstFunc0-refFunc0) / std::abs(refFunc0);
    real funcErr1 = std::abs(tstFunc1-refFunc1) / std::abs(refFunc1);
    real derErr0  = std::abs(tstDer0-refDer0) / std::abs(refDer0);
    real derErr1  = std::abs(tstDer1-refDer1) / std::abs(refDer1);

    // Use asserts, since the following ones compare to these values.
    ASSERT_LT(funcErr0, TypeParam::defaultTolerance);
    ASSERT_LT(derErr0, TypeParam::defaultTolerance);
    ASSERT_LT(funcErr1, TypeParam::defaultTolerance);
    ASSERT_LT(derErr1, TypeParam::defaultTolerance);

    // Test that function/derivative-only interpolation methods work
    table.evaluateFunction(x, &tmpFunc0, &tmpFunc1);
    table.evaluateDerivative(x, &tmpDer0, &tmpDer1);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstDer0, tmpDer0);

    // Test that scrambled order interpolation methods work
    table.template evaluateFunctionAndDerivative<2, 1, 0>(x, &tstFunc1, &tstDer1, &tstFunc0, &tstDer0);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);

    // Test scrambled order for function/derivative-only methods
    table.template evaluateFunction<2, 1, 0>(x, &tmpFunc1, &tmpFunc0);
    table.template evaluateDerivative<2, 1, 0>(x, &tmpDer1, &tmpDer0);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);
}

TYPED_TEST(SplineTableTest, ThreeFunctions)
{
    std::pair<real, real>  range(0.2, 2.0);

    TypeParam              table( {{"Coulomb", coulombFunction, coulombDerivative}, {"LJ6", lj6Function, lj6Derivative}, {"LJ12", lj12Function, lj12Derivative}}, range);

    // Test entire range for each function
    TestFixture::template testSplineTableAgainstFunctions<3, 0>("Coulomb", coulombFunction, coulombDerivative, table, range);
    TestFixture::template testSplineTableAgainstFunctions<3, 1>("LJ6", lj6Function, lj6Derivative, table, range);
    TestFixture::template testSplineTableAgainstFunctions<3, 2>("LJ12", lj12Function, lj12Derivative, table, range);

    // Test the methods that evaluated both functions for one value
    real     x        = 0.5 * (range.first + range.second);
    real     refFunc0 = coulombFunction(x);
    real     refDer0  = coulombDerivative(x);
    real     refFunc1 = lj6Function(x);
    real     refDer1  = lj6Derivative(x);
    real     refFunc2 = lj12Function(x);
    real     refDer2  = lj12Derivative(x);

    real     tstFunc0, tstDer0, tstFunc1, tstDer1, tstFunc2, tstDer2;
    real     tmpFunc0, tmpFunc1, tmpFunc2, tmpDer0, tmpDer1, tmpDer2;

    // test that we reproduce the reference functions
    table.evaluateFunctionAndDerivative(x, &tstFunc0, &tstDer0, &tstFunc1, &tstDer1, &tstFunc2, &tstDer2);

    real funcErr0 = std::abs(tstFunc0-refFunc0) / std::abs(refFunc0);
    real derErr0  = std::abs(tstDer0-refDer0) / std::abs(refDer0);
    real funcErr1 = std::abs(tstFunc1-refFunc1) / std::abs(refFunc1);
    real derErr1  = std::abs(tstDer1-refDer1) / std::abs(refDer1);
    real funcErr2 = std::abs(tstFunc2-refFunc2) / std::abs(refFunc2);
    real derErr2  = std::abs(tstDer2-refDer2) / std::abs(refDer2);

    // Use asserts, since the following ones compare to these values.
    ASSERT_LT(funcErr0, TypeParam::defaultTolerance);
    ASSERT_LT(derErr0, TypeParam::defaultTolerance);
    ASSERT_LT(funcErr1, TypeParam::defaultTolerance);
    ASSERT_LT(derErr1, TypeParam::defaultTolerance);
    ASSERT_LT(funcErr2, TypeParam::defaultTolerance);
    ASSERT_LT(derErr2, TypeParam::defaultTolerance);

    // Test that function/derivative-only interpolation methods work
    table.evaluateFunction(x, &tmpFunc0, &tmpFunc1, &tmpFunc2);
    table.evaluateDerivative(x, &tmpDer0, &tmpDer1, &tmpDer2);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstFunc2, tmpFunc2);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);
    EXPECT_EQ(tstDer2, tmpDer2);

    // Test two functions out of three
    table.template evaluateFunctionAndDerivative<3, 0, 1>(x, &tmpFunc0, &tmpDer0, &tmpFunc1, &tmpDer1);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);

    // two out of three, function/derivative-only
    table.template evaluateFunction<3, 0, 1>(x, &tmpFunc0, &tmpFunc1);
    table.template evaluateDerivative<3, 0, 1>(x, &tmpDer0, &tmpDer1);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);

    // Test that scrambled order interpolation methods work
    table.template evaluateFunctionAndDerivative<3, 2, 1, 0>(x, &tstFunc2, &tstDer2, &tstFunc1, &tstDer1, &tstFunc0, &tstDer0);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstFunc2, tmpFunc2);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);
    EXPECT_EQ(tstDer2, tmpDer2);

    // Test scrambled order for function/derivative-only methods
    table.template evaluateFunction<3, 2, 1, 0>(x, &tmpFunc2, &tmpFunc1, &tmpFunc0);
    table.template evaluateDerivative<3, 2, 1, 0>(x, &tmpDer2, &tmpDer1, &tmpDer0);
    EXPECT_EQ(tstFunc0, tmpFunc0);
    EXPECT_EQ(tstFunc1, tmpFunc1);
    EXPECT_EQ(tstFunc2, tmpFunc2);
    EXPECT_EQ(tstDer0, tmpDer0);
    EXPECT_EQ(tstDer1, tmpDer1);
    EXPECT_EQ(tstDer2, tmpDer2);
}

#if GMX_SIMD_HAVE_REAL
TYPED_TEST(SplineTableTest, Simd)
{
    std::pair<real, real>  range(0.2, 1.0);
    TypeParam              table( {{"LJ12", lj12Function, lj12Derivative}}, range);

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

    EXPECT_LT(funcErr, TypeParam::defaultTolerance);
    EXPECT_LT(derErr, TypeParam::defaultTolerance);
}

TYPED_TEST(SplineTableTest, SimdTwoFunctions)
{
    std::pair<real, real>  range(0.2, 2.0);

    TypeParam              table( {{"LJ6", lj6Function, lj6Derivative}, {"LJ12", lj12Function, lj12Derivative}}, range);

    // We already test that the SIMD operations handle the different elements
    // correctly in the SIMD module, so here we only test that interpolation
    // works for a single value in the middle of the interval

    real     x        = 0.5 * (range.first + range.second);
    real     refFunc0 = lj6Function(x);
    real     refDer0  = lj6Derivative(x);
    real     refFunc1 = lj12Function(x);
    real     refDer1  = lj12Derivative(x);
    SimdReal tstFunc0, tstDer0;
    SimdReal tstFunc1, tstDer1;
    real     funcErr0, derErr0;
    real     funcErr1, derErr1;
    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) alignedMem[GMX_SIMD_REAL_WIDTH];

    table.evaluateFunctionAndDerivative(SimdReal(x), &tstFunc0, &tstDer0, &tstFunc1, &tstDer1);

    store(alignedMem, tstFunc0);
    funcErr0 = std::abs(alignedMem[0]-refFunc0) / std::abs(refFunc0);

    store(alignedMem, tstDer0);
    derErr0  = std::abs(alignedMem[0]-refDer0 ) / std::abs(refDer0 );

    store(alignedMem, tstFunc1);
    funcErr1 = std::abs(alignedMem[0]-refFunc1) / std::abs(refFunc1);

    store(alignedMem, tstDer1);
    derErr1  = std::abs(alignedMem[0]-refDer1 ) / std::abs(refDer1 );

    EXPECT_LT(funcErr0, TypeParam::defaultTolerance);
    EXPECT_LT(derErr0, TypeParam::defaultTolerance);
    EXPECT_LT(funcErr1, TypeParam::defaultTolerance);
    EXPECT_LT(derErr1, TypeParam::defaultTolerance);
}
#endif

#if GMX_SIMD_HAVE_REAL && !defined NDEBUG
TYPED_TEST(SplineTableTest, CatchesOutOfRangeValuesSimd)
{
    std::pair<real, real>                   range(0.2, 1.0);
    TypeParam                               table( {{"LJ12", lj12Function, lj12Derivative}}, range);
    SimdReal                                x, func, der;

    AlignedArray<real, GMX_SIMD_REAL_WIDTH> alignedMem;

    alignedMem.fill(range.first);
    // Make position 1 incorrect if width>=2, otherwise position 0
    // range.first-GMX_REAL_EPS is not invalid. See comment in table.
    alignedMem[ (GMX_SIMD_REAL_WIDTH >= 2) ? 1 : 0] = -GMX_REAL_EPS;
    x = load<SimdReal>(alignedMem);

    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);

    // Make position 1 incorrect if width>=2, otherwise position 0
    alignedMem[ (GMX_SIMD_REAL_WIDTH >= 2) ? 1 : 0] = range.second;
    x = load<SimdReal>(alignedMem);

    EXPECT_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der), gmx::RangeError);
}

TYPED_TEST(SplineTableTest, AcceptsInRangeValuesSimd)
{
    std::pair<real, real>  range(0.2, 1.0);
    TypeParam              table( {{"LJ12", lj12Function, lj12Derivative}}, range);
    SimdReal               x, func, der;

    GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) alignedMem[GMX_SIMD_REAL_WIDTH];

    // Test all values between 0 and range.second
    for (std::size_t i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    {
        alignedMem[i] = range.second*(1.0-GMX_REAL_EPS)*i/(GMX_SIMD_REAL_WIDTH-1);
    }
    x = load<SimdReal>(alignedMem);

    EXPECT_NO_THROW_GMX(table.evaluateFunctionAndDerivative(x, &func, &der));
}

#endif

} // namespace

} // namespace test

} // namespace gmx
