/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * This file provides assert macros to replace (ASSERT|EXPECT)(_NO)?_THROW
 * from Google Test.  They behave otherwise the same as the Google Test ones,
 * but also print details of any unexpected exceptions.  This makes it much
 * easier to see at one glance what went wrong.
 *
 * This file also provides extra floating-point assertions.
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

#include <gtest/gtest.h>

#include "gromacs/legacyheaders/maths.h"

#include "gromacs/utility/exceptions.h"

namespace gmx
{
namespace test
{

/*! \cond internal */
/*! \internal
 * \brief
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
            << "Exception details:\n" << ::gmx::formatException(ex); \
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

/*! \internal
 * \brief
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
            << "Exception details:\n" << ::gmx::formatException(ex); \
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
 * See Google Test documentation on EXPECT_THROW.
 * This macro works the same, but additionally prints details of unexpected
 * exceptions.
 */
#define EXPECT_THROW_GMX(statement, expected_exception) \
    GMX_TEST_THROW_(statement, expected_exception, GTEST_NONFATAL_FAILURE_)
/*! \brief
 * Asserts that a statement does not throw.
 *
 * See Google Test documentation on EXPECT_NO_THROW.
 * This macro works the same, but additionally prints details of unexpected
 * exceptions.
 */
#define EXPECT_NO_THROW_GMX(statement) \
    GMX_TEST_NO_THROW_(statement, GTEST_NONFATAL_FAILURE_)
/*! \brief
 * Asserts that a statement throws a given exception.
 *
 * See Google Test documentation on ASSERT_THROW.
 * This macro works the same, but additionally prints details of unexpected
 * exceptions.
 */
#define ASSERT_THROW_GMX(statement, expected_exception) \
    GMX_TEST_THROW_(statement, expected_exception, GTEST_FATAL_FAILURE_)
/*! \brief
 * Asserts that a statement does not throw.
 *
 * See Google Test documentation on ASSERT_NO_THROW.
 * This macro works the same, but additionally prints details of unexpected
 * exceptions.
 */
#define ASSERT_NO_THROW_GMX(statement) \
    GMX_TEST_NO_THROW_(statement, GTEST_FATAL_FAILURE_)

/*! \cond internal */
/*! \internal \brief
 * Assertion predicate formatter for comparing two floating-point values.
 */
static ::testing::AssertionResult assertWithinRelativeTolerance(
        const char *expr1, const char *expr2, const char * /*exprTolerance*/,
        real val1, real val2, real relativeTolerance)
{
    if (gmx_within_tol(val1, val2, relativeTolerance))
    {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "Value of: " << expr2 << " not within tolerance of " << relativeTolerance << "\n"
           << "  Actual: " << val2 << "\n"
           << "Expected: " << expr1 << "\n"
           << "Which is: " << val1;
}
//! \endcond

/*! \brief
 * Asserts that two floating-point values are within the given relative error.
 *
 * This assert works as EXPECT_NEAR from Google Test, except that it uses a
 * relative instead of a absolute tolerance.
 * See gmx_within_tol() for definition of the tolerance.
 */
#define EXPECT_NEAR_REL(val1, val2, rel_error) \
    EXPECT_PRED_FORMAT3(::gmx::test::assertWithinRelativeTolerance, val1, val2, rel_error)
/*! \brief
 * Asserts that two floating-point values are within the given relative error.
 *
 * This assert works as ASSERT_NEAR from Google Test, except that it uses a
 * relative instead of a absolute tolerance.
 * See gmx_within_tol() for definition of the tolerance.
 */
#define ASSERT_NEAR_REL(val1, val2, rel_error) \
    ASSERT_PRED_FORMAT3(::gmx::test::assertWithinRelativeTolerance, val1, val2, rel_error)

} // namespace test
} // namespace gmx

#endif
