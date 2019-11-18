/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  - comparison against NULL
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

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

namespace internal
{
//! \cond internal
/*! \internal
 * \brief
 * Called for an expected exception from EXPECT_THROW_GMX().
 *
 * \param[in] ex  Exception that was thrown.
 */
void processExpectedException(const std::exception& ex);
//! \endcond
} // namespace internal

//! \libinternal \addtogroup module_testutils
//! \{

/*! \name Assertions for exceptions
 *
 * These macros replace `(ASSERT|EXPECT)(_NO)?_THROW` from Google Test.
 * They are used exactly like the Google Test ones, but also print details of
 * any unexpected exceptions using \Gromacs-specific routines.
 * This makes it much easier to see at one glance what went wrong.
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
#define GMX_TEST_THROW_(statement, expected_exception, fail)                           \
    GTEST_AMBIGUOUS_ELSE_BLOCKER_                                                      \
    if (::testing::AssertionResult gmx_ar = ::testing::AssertionSuccess())             \
    {                                                                                  \
        bool gmx_caught_expected = false;                                              \
        try                                                                            \
        {                                                                              \
            GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement);                 \
        }                                                                              \
        catch (expected_exception const& ex)                                           \
        {                                                                              \
            gmx_caught_expected = true;                                                \
            ::gmx::test::internal::processExpectedException(ex);                       \
        }                                                                              \
        catch (std::exception const& ex)                                               \
        {                                                                              \
            gmx_ar << "Expected: " #statement " throws an exception of type "          \
                   << #expected_exception ".\n  Actual: it throws a different type.\n" \
                   << "Exception details:\n"                                           \
                   << ::gmx::formatExceptionMessageToString(ex);                       \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__);                  \
        }                                                                              \
        catch (...)                                                                    \
        {                                                                              \
            gmx_ar << "Expected: " #statement " throws an exception of type "          \
                   << #expected_exception ".\n  Actual: it throws a different type.";  \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__);                  \
        }                                                                              \
        if (!gmx_caught_expected)                                                      \
        {                                                                              \
            gmx_ar << "Expected: " #statement " throws an exception of type "          \
                   << #expected_exception ".\n  Actual: it throws nothing.";           \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__);                  \
        }                                                                              \
    }                                                                                  \
    else                                                                               \
        GTEST_CONCAT_TOKEN_(gmx_label_testthrow_, __LINE__) : fail(gmx_ar.message())

/*! \brief
 * Internal implementation macro for exception assertations.
 *
 * \param statement          Statements to execute.
 * \param fail               Function/macro to call on failure.
 *
 * The implementation is copied and adjusted from
 * include/gtest/internal/gtest-internal.h in Google Test 1.6.0.
 */
#define GMX_TEST_NO_THROW_(statement, fail)                                    \
    GTEST_AMBIGUOUS_ELSE_BLOCKER_                                              \
    if (::testing::AssertionResult gmx_ar = ::testing::AssertionSuccess())     \
    {                                                                          \
        try                                                                    \
        {                                                                      \
            GTEST_SUPPRESS_UNREACHABLE_CODE_WARNING_BELOW_(statement);         \
        }                                                                      \
        catch (std::exception const& ex)                                       \
        {                                                                      \
            gmx_ar << "Expected: " #statement " doesn't throw an exception.\n" \
                   << "  Actual: it throws.\n"                                 \
                   << "Exception details:\n"                                   \
                   << ::gmx::formatExceptionMessageToString(ex);               \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__);        \
        }                                                                      \
        catch (...)                                                            \
        {                                                                      \
            gmx_ar << "Expected: " #statement " doesn't throw an exception.\n" \
                   << "  Actual: it throws.";                                  \
            goto GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__);        \
        }                                                                      \
    }                                                                          \
    else                                                                       \
        GTEST_CONCAT_TOKEN_(gmx_label_testnothrow_, __LINE__) : fail(gmx_ar.message())
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
#define EXPECT_NO_THROW_GMX(statement) GMX_TEST_NO_THROW_(statement, GTEST_NONFATAL_FAILURE_)
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
#define ASSERT_NO_THROW_GMX(statement) GMX_TEST_NO_THROW_(statement, GTEST_FATAL_FAILURE_)

//! \}

/*! \libinternal \brief
 * Computes and represents a floating-point difference value.
 *
 * Methods in this class do not throw, except for toString(), which may throw
 * std::bad_alloc.
 *
 * \see FloatingPointTolerance
 */
class FloatingPointDifference
{
public:
    /*! \brief Initializes a single-precision difference.
     *
     *  \param ref    First term in difference
     *  \param value  Second term in difference
     *
     *  For absolute and ULP differences the two parameters are equivalent,
     *  since the difference is symmetric. For relative differences
     *  the first term is interpreted as the reference value, from which
     *  we extract the magnitude to compare with.
     */
    FloatingPointDifference(float ref, float value);

    /*! \brief Initializes a double-precision difference.
     *
     *  \param ref    First term in difference
     *  \param value  Second term in difference
     *
     *  For absolute and ULP differences the two parameters are equivalent,
     *  since the difference is symmetric. For relative differences
     *  the first term is interpreted as the reference value, from which
     *  we extract the magnitude to compare with.
     */
    FloatingPointDifference(double ref, double value);

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
    uint64_t asUlps() const { return ulpDifference_; }
    /*! \brief
     * Whether the compared values were of different sign.
     *
     * 0.0 and -0.0 are treated as positive and negative, respectively.
     */
    bool signsDiffer() const { return bSignDifference_; }
    /*! \brief
     * Whether the difference is between single- or double-precision
     * numbers.
     */
    bool isDouble() const { return bDouble_; }
    //! Formats the difference as a string for assertion failure messages.
    std::string toString() const;

    //! Returns the magnitude of the original second term of the difference.
    double termMagnitude() const { return termMagnitude_; }

private:
    //! Save the magnitude of the reference value for relative (i.e., not ULP) tolerance
    double termMagnitude_;
    //! Stores the absolute difference, or NaN if one or both values were NaN.
    double   absoluteDifference_;
    uint64_t ulpDifference_;
    bool     bSignDifference_;
    /*! \brief
     * Whether the difference was computed for single or double precision.
     *
     * This sets the units for `ulpDifference_`.
     */
    bool bDouble_;
};

/*! \libinternal \brief
 * Specifies a floating-point comparison tolerance and checks whether a
 * difference is within the tolerance.
 *
 * The related functions section lists methods that can be construct methods
 * using less parameters than the full constructor, and with more obvious
 * semantics.  These should be preferred over using the constructor directly.
 *
 * Several types of tolerances are possible:
 *  - _absolute tolerance_: difference between the values must be smaller than
 *    the given tolerance for the check to pass.
 *    Setting the absolute tolerance to zero disables the absolute tolerance
 *    check.
 *  - _relative tolerance_: the absolute difference between the numbers must
 *    be smaller than the tolerance multiplied by the first number. Setting
 *    the relative tolerance to zero disables this check.
 *  - _ULP tolerance_: ULP (units of least precision) difference between the
 *    values must be smaller than the given tolerance for the check to pass.
 *    Setting the ULP tolerance to zero requires exact match.
 *    Setting the ULP tolerance to UINT64_MAX disables the ULP check.
 *    `0.0` and `-0.0` are treated as equal for the ULP check.
 *  - _sign check_: if set, any values that are of different signs fail the
 *    check (note that this also applies to `0.0` and `-0.0`: a value with a
 *    different sign than the zero will fail the check).
 *
 * Either an absolute, relative, or ULP tolerance must always be specified.
 * If several of them are specified, then the check passes if either of the
 * tolerances is satisfied.
 *
 * Any combination of absolute, relative, and ULP tolerance can be combined with
 * the sign check.  In this case, the sign check must succeed for the check to
 * pass, even if other tolerances are satisfied.
 *
 * The tolerances can be specified separately for single and double precision
 * comparison.  Different initialization functions have different semantics on
 * how the provided tolerance values are interpreted; check their
 * documentation.
 *
 * Methods in this class do not throw, except for toString(), which may throw
 * std::bad_alloc.
 *
 * \todo
 * The factory methods that take ULP difference could be better formulated as
 * methods that take the acceptable number of incorrect bits and/or the number
 * of accurate bits.
 *
 * \see FloatingPointDifference
 */
class FloatingPointTolerance
{
public:
    /*! \brief
     * Creates a tolerance with the specified values.
     *
     * \param[in]  singleAbsoluteTolerance
     *     Allowed absolute difference in a single-precision number.
     * \param[in]  doubleAbsoluteTolerance
     *     Allowed absolute difference in a double-precision number.
     * \param[in]  singleRelativeTolerance
     *     Allowed relative difference in a single-precision number.
     * \param[in]  doubleRelativeTolerance
     *     Allowed relative difference in a double-precision number.
     * \param[in]  singleUlpTolerance
     *     Allowed ULP difference in a single-precision number.
     * \param[in]  doubleUlpTolerance
     *     Allowed ULP difference in a double-precision number.
     * \param[in]  bSignMustMatch
     *     Whether sign mismatch fails the comparison.
     */
    FloatingPointTolerance(float    singleAbsoluteTolerance,
                           double   doubleAbsoluteTolerance,
                           float    singleRelativeTolerance,
                           double   doubleRelativeTolerance,
                           uint64_t singleUlpTolerance,
                           uint64_t doubleUlpTolerance,
                           bool     bSignMustMatch) :
        singleAbsoluteTolerance_(singleAbsoluteTolerance),
        doubleAbsoluteTolerance_(doubleAbsoluteTolerance),
        singleRelativeTolerance_(singleRelativeTolerance),
        doubleRelativeTolerance_(doubleRelativeTolerance),
        singleUlpTolerance_(singleUlpTolerance),
        doubleUlpTolerance_(doubleUlpTolerance),
        bSignMustMatch_(bSignMustMatch)
    {
    }

    /*! \brief
     * Checks whether a difference is within the specified tolerance.
     *
     * NaNs are always treated outside the tolerance.
     */
    bool isWithin(const FloatingPointDifference& difference) const;

    //! Formats the tolerance as a string for assertion failure messages.
    std::string toString(const FloatingPointDifference& difference) const;

private:
    float    singleAbsoluteTolerance_;
    double   doubleAbsoluteTolerance_;
    float    singleRelativeTolerance_;
    double   doubleRelativeTolerance_;
    uint64_t singleUlpTolerance_;
    uint64_t doubleUlpTolerance_;
    bool     bSignMustMatch_;
};

/*! \brief
 * Creates a tolerance that only allows a specified ULP difference.
 *
 * The tolerance uses the given ULP value for both precisions, i.e., double
 * precision will have much stricter tolerance.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance ulpTolerance(uint64_t ulpDiff)
{
    return FloatingPointTolerance(0.0, 0.0, 0.0, 0.0, ulpDiff, ulpDiff, false);
}

/*! \brief
 * Creates a tolerance that allows a difference in two compared values that is
 * relative to the given magnitude.
 *
 * \param[in] magnitude  Magnitude of the numbers the computation operates in.
 * \param[in] tolerance  Relative tolerance permitted (e.g. 1e-4).
 *
 * In addition to setting an relative tolerance for both
 * precisions, this sets the absolute tolerance such that values close to zero
 * (in general, smaller than \p magnitude) do not fail the check if they
 * differ by less than \p tolerance evaluated at \p magnitude.  This accounts
 * for potential loss of precision for small values, and should be used when
 * accuracy of values much less than \p magnitude do not matter for
 * correctness.
 *
 * \related FloatingPointTolerance
 */
FloatingPointTolerance relativeToleranceAsFloatingPoint(double magnitude, double tolerance);

/*! \brief
 * Creates a tolerance that allows a precision-dependent difference in two
 * compared values that is relative to the given magnitude.
 *
 * \param[in] magnitude        Magnitude of the numbers the computation
 *     operates in.
 * \param[in] singleTolerance  Relative tolerance permitted (e.g. 1e-4)
 *     in single precision.
 * \param[in] doubleTolerance  Relative tolerance permitted (e.g. 1e-4)
 *     in double precision.
 *
 * In addition to setting an relative tolerance for both
 * precisions, this sets the absolute tolerance such that values close to zero
 * (in general, smaller than \p magnitude) do not fail the check if they
 * differ by less than \p tolerance evaluated at \p magnitude.  This accounts
 * for potential loss of precision for small values, and should be used when
 * accuracy of values much less than \p magnitude do not matter for
 * correctness.
 *
 * \related FloatingPointTolerance
 */
FloatingPointTolerance relativeToleranceAsPrecisionDependentFloatingPoint(double magnitude,
                                                                          float  singleTolerance,
                                                                          double doubleTolerance);

/*! \brief
 * Creates a tolerance that allows a precision-dependent relative difference in
 * a complex computation.
 *
 * \param[in] magnitude      Magnitude of the numbers the computation operates in.
 * \param[in] singleUlpDiff  Expected accuracy of single-precision
 *     computation (in ULPs).
 * \param[in] doubleUlpDiff  Expected accuracy of double-precision
 *     computation (in ULPs).
 *
 * This works as relativeToleranceAsUlp(), but allows setting the ULP
 * difference separately for the different precisions.  This supports
 * cases where the double-precision calculation can acceptably has a higher ULP
 * difference, but relaxing the single-precision tolerance would lead to an
 * unnecessarily loose test.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance relativeToleranceAsPrecisionDependentUlp(double   magnitude,
                                                                              uint64_t singleUlpDiff,
                                                                              uint64_t doubleUlpDiff)
{
    return FloatingPointTolerance(float(magnitude) * singleUlpDiff * GMX_FLOAT_EPS,
                                  magnitude * doubleUlpDiff * GMX_DOUBLE_EPS, 0.0, 0.0,
                                  singleUlpDiff, doubleUlpDiff, false);
}

/*! \brief
 * Creates a tolerance that allows a specified absolute difference.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance absoluteTolerance(double tolerance)
{
    return FloatingPointTolerance(float(tolerance), tolerance, 0.0, 0.0, UINT64_MAX, UINT64_MAX, false);
}

/*! \brief
 * Creates a tolerance that allows a relative difference in a complex
 * computation.
 *
 * \param[in] magnitude  Magnitude of the numbers the computation operates in.
 * \param[in] ulpDiff    Expected accuracy of the computation (in ULPs).
 *
 * In addition to setting the ULP tolerance as ulpTolerance(), this sets the
 * absolute tolerance such that values close to zero (in general, smaller than
 * \p magnitude) do not fail the check if they differ by less than \p ulpDiff
 * evaluated at \p magnitude.  This accounts for potential loss of precision
 * for small values, and should be used when accuracy of values much less than
 * \p magnitude do not matter for correctness.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance relativeToleranceAsUlp(double magnitude, uint64_t ulpDiff)
{
    return relativeToleranceAsPrecisionDependentUlp(magnitude, ulpDiff, ulpDiff);
}

namespace detail
{
//! Default tolerance in ULPs for two floating-point values to compare equal.
constexpr uint64_t g_defaultUlpTolerance = 4;
} // namespace detail

/*! \brief
 * Returns the default tolerance for comparing `real` numbers.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance defaultRealTolerance()
{
    return relativeToleranceAsUlp(1.0, detail::g_defaultUlpTolerance);
}


/*! \brief
 * Returns the default tolerance for comparing single-precision numbers when
 * compared by \Gromacs built in either precision mode.
 *
 * This permits a checker compiled with any \Gromacs precision to compare
 * equal or not in the same way.
 *
 * \related FloatingPointTolerance
 */
static inline FloatingPointTolerance defaultFloatTolerance()
{
    return relativeToleranceAsPrecisionDependentUlp(
            1.0, detail::g_defaultUlpTolerance,
            static_cast<uint64_t>(detail::g_defaultUlpTolerance * (GMX_FLOAT_EPS / GMX_DOUBLE_EPS)));
}

/*! \name Assertions for floating-point comparison
 *
 * These routines extend `(EXPECT|ASSERT)_(FLOAT|DOUBLE)_EQ` and
 * `(EXPECT|ASSERT)_NEAR` from Google Test to provide more flexible assertions
 * for floating-point values.
 *
 * See gmx::test::FloatingPointTolerance for the possible ways to specify the
 * tolerance, and gmx::test::FloatingPointDifference for some additional
 * details of the difference calculation.
 */
//! \{

//! \cond internal
/*! \internal \brief
 * Assertion predicate formatter for comparing two floating-point values.
 */
template<typename FloatType>
static inline ::testing::AssertionResult assertEqualWithinTolerance(const char* expr1,
                                                                    const char* expr2,
                                                                    const char* /*exprTolerance*/,
                                                                    FloatType value1,
                                                                    FloatType value2,
                                                                    const FloatingPointTolerance& tolerance)
{
    FloatingPointDifference diff(value1, value2);
    if (tolerance.isWithin(diff))
    {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << "  Value of: " << expr2 << std::endl
                                         << "    Actual: " << value2 << std::endl
                                         << "  Expected: " << expr1 << std::endl
                                         << "  Which is: " << value1 << std::endl
                                         << "Difference: " << diff.toString() << std::endl
                                         << " Tolerance: " << tolerance.toString(diff);
}
//! \endcond

/*! \brief
 * Asserts that two single-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define EXPECT_FLOAT_EQ_TOL(value1, value2, tolerance) \
    EXPECT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<float>, value1, value2, tolerance)
/*! \brief
 * Asserts that two double-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define EXPECT_DOUBLE_EQ_TOL(value1, value2, tolerance) \
    EXPECT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<double>, value1, value2, tolerance)
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
    ASSERT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<float>, value1, value2, tolerance)
/*! \brief
 * Asserts that two double-precision values are within the given tolerance.
 *
 * \hideinitializer
 */
#define ASSERT_DOUBLE_EQ_TOL(value1, value2, tolerance) \
    ASSERT_PRED_FORMAT3(::gmx::test::assertEqualWithinTolerance<double>, value1, value2, tolerance)
/*! \def ASSERT_REAL_EQ_TOL
 * \brief
 * Asserts that two `real` values are within the given tolerance.
 *
 * \hideinitializer
 */

#if GMX_DOUBLE
#    define EXPECT_REAL_EQ_TOL(value1, value2, tolerance) \
        EXPECT_DOUBLE_EQ_TOL(value1, value2, tolerance)
#    define ASSERT_REAL_EQ_TOL(value1, value2, tolerance) \
        ASSERT_DOUBLE_EQ_TOL(value1, value2, tolerance)
#else
#    define EXPECT_REAL_EQ_TOL(value1, value2, tolerance) \
        EXPECT_FLOAT_EQ_TOL(value1, value2, tolerance)
#    define ASSERT_REAL_EQ_TOL(value1, value2, tolerance) \
        ASSERT_FLOAT_EQ_TOL(value1, value2, tolerance)
#endif

//! EXPECT_REAL_EQ_TOL with default tolerance
#define EXPECT_REAL_EQ(value1, value2) \
    EXPECT_REAL_EQ_TOL(value1, value2, ::gmx::test::defaultRealTolerance())
//! ASSERT_REAL_EQ_TOL with default tolerance
#define ASSERT_REAL_EQ(value1, value2) \
    ASSERT_REAL_EQ_TOL(value1, value2, ::gmx::test::defaultRealTolerance())
//! \}

//! \cond internal
/*! \internal \brief
 * Helper method for `(EXPECT|ASSERT)_PLAIN`.
 */
static inline ::testing::AssertionResult plainAssertHelper(const char* /*expr*/,
                                                           const ::testing::AssertionResult& expr)
{
    return expr;
}
//! \endcond

/*! \brief
 * Assert for predicates that return AssertionResult and produce a full failure
 * message.
 *
 * `expr` should evaluate to AssertionResult, and on failure the message from
 * the result is used as-is, unlike in EXPECT_TRUE().
 *
 * \hideinitializer
 */
#define EXPECT_PLAIN(expr) EXPECT_PRED_FORMAT1(plainAssertHelper, expr)
/*! \brief
 * Assert for predicates that return AssertionResult and produce a full failure
 * message.
 *
 * `expr` should evaluate to AssertionResult, and on failure the message from
 * the result is used as-is, unlike in ASSERT_TRUE().
 *
 * \hideinitializer
 */
#define ASSERT_PLAIN(expr) ASSERT_PRED_FORMAT1(plainAssertHelper, expr)

//! \}

} // namespace test
} // namespace gmx

#endif
