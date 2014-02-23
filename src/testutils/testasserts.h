/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Extra assertions for unit tests.
 *
 * This file provides assertion macros that extend/replace Google Test
 * assertions for:
 *  - exceptions
 *  - floating-point comparison
 *
 * \if internal
 * \todo
 * The implementation is somewhat ugly, and accesses some Google Test
 * internals.  Could be nice to clean it up a bit.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_TESTASSERTS_H
#define GMX_TESTUTILS_TESTASSERTS_H

#include <string>

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/types/simple.h"

#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace test
{

//! \addtogroup module_testutils
//! \{

/*! \name Assertions for exceptions
 *
 * These macros replace `(ASSERT|EXPECT)(_NO)?_THROW` from Google Test.
 * They behave otherwise the same as the Google Test ones, but also print
 * details of any unexpected exceptions using \Gromacs-specific routines.  This
 * makes it much easier to see at one glance what went wrong.
 * See Google Test documentation for details on how to use the macros.
 */
//! \{

//! \cond internal
/*! \brief
 * Internal implementation macro for exception assertations.
 *
 * \param statement          Statements to execute.
 * \param expected_exception Exception type that \p statement should throw.
 * \param fail               Function/macro to call on failure.
 *
 * The implementation is copied and adjusted from
 * include/gtest/internal/gtest-internal.h in Google Test 1.6.0.
 */
#define GMX_TEST_THROW_(statement, expected_exception, fail) \
    GTEST_AMBIGUOUS_ELSE_BLOCKER_ \
    if (::testing::AssertionResult gmx_ar = ::testing::AssertionSuccess()) { \
        bool gmx_caught_expected = false; \
        try { \
            GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement); \
        } \
        catch (expected_exception const &) { \
            gmx_caught_expected = true; \
        } \
        catch (std::exception const &ex) { \
            gmx_ar << "Expected: " #statement " throws an exception of type " \
            << #expected_exception ".\n  Actual: it throws a different type.\n" \
            << "Exception details:\n" << ::gmx::formatExceptionMessageToString(ex); \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__); \
        } \
        catch (...) { \
            gmx_ar << "Expected: " #statement " throws an exception of type " \
            << #expected_exception ".\n  Actual: it throws a different type."; \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__); \
        } \
        if (!gmx_caught_expected) { \
            gmx_ar << "Expected: " #statement " throws an exception of type " \
            << #expected_exception ".\n  Actual: it throws nothing."; \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__); \
        } \
    } else \
        GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__) : \
            fail(gmx_ar.message())

/*! \brief
 * Internal implementation macro for exception assertations.
 *
 * \param statement          Statements to execute.
 * \param fail               Function/macro to call on failure.
 *
 * The implementation is copied and adjusted from
 * include/gtest/internal/gtest-internal.h in Google Test 1.6.0.
 */
#define GMX_TEST_NO_THROW_(statement, fail) \
    GTEST_AMBIGUOUS_ELSE_BLOCKER_ \
    if (::testing::AssertionResult gmx_ar = ::testing::AssertionSuccess()) { \
        try { \
            GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement); \
        } \
        catch (std::exception const &ex) { \
            gmx_ar << "Expected: " #statement " doesn't throw an exception.\n" \
            << "  Actual: it throws.\n" \
            << "Exception details:\n" << ::gmx::formatExceptionMessageToString(ex); \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__); \
        } \
        catch (...) { \
            gmx_ar << "Expected: " #statement " doesn't throw an exception.\n" \
            << "  Actual: it throws."; \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__); \
        } \
    } else \
        GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__) : \
            fail(gmx_ar.message())
//! \endcond

/*! \brief
 * Asserts that a statement throws a given exception.
 *
 * \hideinitializer
 */
#define EXPECT_THROW_GMX(statement, expected_exception) \
    GMX_TEST_THROW_(statement, expected_exception, GTEST_NONFATAL_FAILURE_)
/*! \brief
 * Asserts that a statement does not throw.
 *
 * \hideinitializer
 */
#define EXPECT_NO_THROW_GMX(statement) \
    GMX_TEST_NO_THROW_(statement, GTEST_NONFATAL_FAILURE_)
/*! \brief
 * Asserts that a statement throws a given exception.
 *
 * \hideinitializer
 */
#define ASSERT_THROW_GMX(statement, expected_exception) \
    GMX_TEST_THROW_(statement, expected_exception, GTEST_FATAL_FAILURE_)
/*! \brief
 * Asserts that a statement does not throw.
 *
 * \hideinitializer
 */
#define ASSERT_NO_THROW_GMX(statement) \
    GMX_TEST_NO_THROW_(statement, GTEST_FATAL_FAILURE_)

//! \}

/*! \name Assertions for floating-point comparison
 *
 * These routines extend `(EXPECT|ASSERT)_(FLOAT|DOUBLE)_EQ` and
 * `(EXPECT|ASSERT)_NEAR` from Google Test to provide more flexible assertions
 * for floating-point values.
 *
 * See FloatingPointTolerance for the possible ways to specify the tolerance.
 */
//! \{

/*! \libinternal \brief
 * Computes and represents a floating-point difference value.
 *
 * Methods in this class do not throw.
 */
class FloatingPointDifference
{
    public:
        //! Initializes a single-precision difference.
        FloatingPointDifference(float value1, float value2);
        //! Initializes a double-precision difference.
        FloatingPointDifference(double value1, double value2);

        /*! \brief
         * Whether one or both of the compared values were NaN.
         *
         * If this returns `true`, other accessors return meaningless values.
         */
        bool isNaN() const;
        //! Returns the difference as an absolute number (always non-negative).
        double asAbsolute() const { return absoluteDifference_; }
        /*! \brief
         * Returns the difference as ULPs (always non-negative).
         *
         * The ULPs are calculated for the type that corresponds to the
         * constructor used to initialize the difference.
         * The ULP difference between 0.0 and -0.0 is zero.
         */
        gmx_int64_t asUlps() const { return ulpDifference_; }
        /*! \brief
         * Whether the compared values were of different sign.
         *
         * 0.0 and -0.0 are treated as having the same sign with both positive
         * and negative numbers.
         */
        bool signsDiffer() const { return bSignDifference_; }

        //! Formats the difference as a string for assertion failure messages.
        std::string toString() const;

    private:
        //! Stores the absolute difference, or NaN if one or both values were NaN.
        double      absoluteDifference_;
        gmx_int64_t ulpDifference_;
        bool        bSignDifference_;
        /*! \brief
         * Whether the difference was computed for single or double precision.
         *
         * This sets the units for `ulpDifference_`.
         */
        bool        bDouble_;
};

/*! \libinternal \brief
 * Specifies a floating-point comparison tolerance and checks whether a
 * difference is within the tolerance.
 *
 * Several types of tolerances are possible:
 *  - _absolute tolerance_: difference between the values must be smaller than
 *    the given tolerance for the check to pass.
 *    Setting the absolute tolerance to zero disables the absolute tolerance
 *    check.
 *  - _ULP tolerance_: ULP (units of least precision) difference between the
 *    values must be smaller than the given tolerance for the check to pass.
 *    Setting the ULP tolerance to zero requires exact match.
 *    Setting the ULP tolerance to negative disables the ULP check.
 *    `0.0` and `-0.0` are treated as equal for the ULP check.
 *  - _sign check_: if set, any values that are of different signs fail the
 *    check (except if one of the values is `0.0` or `-0.0`, in which case
 *    either sign is acceptable).
 *
 * Either an absolute or a ULP tolerance must always be specified.
 * If both are specified, then the check passes if either of the tolerances is
 * satisfied.
 *
 * Any combination absolute and ULP tolerance can be combined with the sign
 * check.  In this case, the sign check must succeed for the check to pass,
 * even if other tolerances are satisfied.
 */
class FloatingPointTolerance
{
    public:
        /*! \brief
         * Creates a tolerance with the specified values.
         *
         * \param[in] absolute       Allowed absolute difference.
         * \param[in] ulp            Allowed ULP difference.
         * \param[in] bSignMustMatch Whether sign mismatch fails the comparison.
         */
        FloatingPointTolerance(double absolute, gmx_int64_t ulp,
                               bool bSignMustMatch)
            : absoluteTolerance_(absolute), ulpTolerance_(ulp),
              bSignMustMatch_(bSignMustMatch)
        {
        }

        /*! \brief
         * Checks whether a difference is within the specified tolerance.
         *
         * NaNs are always treated outside the tolerance.
         */
        bool isWithin(const FloatingPointDifference &difference) const;

        //! Formats the tolerance as a string for assertion failure messages.
        std::string toString() const;

    private:
        double      absoluteTolerance_;
        gmx_int64_t ulpTolerance_;
        bool        bSignMustMatch_;
};

/*! \brief
 * Creates a tolerance that only allows a specified ULP difference.
 */
static inline FloatingPointTolerance
ulpTolerance(gmx_int64_t ulpDiff)
{
    return FloatingPointTolerance(0.0, ulpDiff, false);
}

/*! \brief
 * Creates a tolerance that allows a relative difference in a complex computation.
 */
static inline FloatingPointTolerance
relativeRealTolerance(double magnitude, gmx_int64_t ulpDiff)
{
    return FloatingPointTolerance(magnitude*ulpDiff*GMX_REAL_EPS, ulpDiff, false);
}

//! Returns the default tolerance for comparing `real` numbers.
static inline FloatingPointTolerance defaultRealTolerance()
{
    return relativeRealTolerance(1.0, 4);
}

//! \cond internal
/*! \internal \brief
 * Assertion predicate formatter for comparing two floating-point values.
 */
template <typename FloatType>
static inline ::testing::AssertionResult assertEqualWithinTolerance(
        const char *expr1, const char *expr2, const char * /*exprTolerance*/,
        FloatType value1, FloatType value2,
        const FloatingPointTolerance &tolerance)
{
    FloatingPointDifference diff(value1, value2);
    if (tolerance.isWithin(diff))
    {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "  Value of: " << expr2 << std::endl
           << "    Actual: " << value2 << std::endl
           << "  Expected: " << expr1 << std::endl
           << "  Which is: " << value1 << std::endl
           << "Difference: " << diff.toString() << std::endl
           << " Tolerance: " << tolerance.toString();
}
//! \endcond

/*! \brief
 * Asserts that two single-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define EXPECT_FLOAT_EQ_TOL(value1, value2, tolerance) \
    EXPECT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<float>, \
                        value1, value2, tolerance)
/*! \brief
 * Asserts that two double-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define EXPECT_DOUBLE_EQ_TOL(value1, value2, tolerance) \
    EXPECT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<double>, \
                        value1, value2, tolerance)
/*! \def EXPECT_REAL_EQ_TOL
 * \brief
 * Asserts that two `real` values are within the given tolerance.
 *
 * \hideinitializer
 */
/*! \brief
 * Asserts that two single-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define ASSERT_FLOAT_EQ_TOL(value1, value2, tolerance) \
    ASSERT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<float>, \
                        value1, value2, tolerance)
/*! \brief
 * Asserts that two double-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define ASSERT_DOUBLE_EQ_TOL(value1, value2, tolerance) \
    ASSERT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<double>, \
                        value1, value2, tolerance)
/*! \def ASSERT_REAL_EQ_TOL
 * \brief
 * Asserts that two `real` values are within the given tolerance.
 *
 * \hideinitializer
 */

#ifdef GMX_DOUBLE
#define EXPECT_REAL_EQ_TOL(value1, value2, tolerance) \
    EXPECT_DOUBLE_EQ_TOL(value1, value2, tolerance)
#define ASSERT_REAL_EQ_TOL(value1, value2, tolerance) \
    ASSERT_DOUBLE_EQ_TOL(value1, value2, tolerance)
#else
#define EXPECT_REAL_EQ_TOL(value1, value2, tolerance) \
    EXPECT_FLOAT_EQ_TOL(value1, value2, tolerance)
#define ASSERT_REAL_EQ_TOL(value1, value2, tolerance) \
    ASSERT_FLOAT_EQ_TOL(value1, value2, tolerance)
#endif

//! \}
//! \}

} // namespace test
} // namespace gmx

#endif
