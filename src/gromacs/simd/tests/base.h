/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

#ifndef GMX_SIMD_TESTS_BASE_H
#define GMX_SIMD_TESTS_BASE_H

/*! \internal \file
 * \brief
 * Declares common base class for testing SIMD and SIMD4.
 *
 * The base class contains the settings for absolute and ulp tolerances,
 * as well as testing ranges used for both SIMD and SIMD4 tests, mainly
 * to keep everything symmetric and clean. The class also defines a couple
 * of generic tests that compare vectors of elements with arbitrary length for
 * either exact or approximate matching (in terms of ulp). These are used in
 * derived classes that convert either SIMD or SIMD4 values to
 * std::vector<real> and then performs the comparison.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 * \ingroup module_simd
 */
#include "config.h"

#include <cstdint>

#include <limits>
#include <ostream>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

//! \internal \brief Test-time utility macro for current precision accuracy
#define GMX_SIMD_ACCURACY_BITS_REAL \
    (GMX_DOUBLE ? GMX_SIMD_ACCURACY_BITS_DOUBLE : GMX_SIMD_ACCURACY_BITS_SINGLE)

/*! \internal
 * \brief
 * Base class for SIMD test fixtures.
 *
 * This class contains settings that are common for SIMD and SIMD4 tests,
 * and it is thus not used directly for any tests, but derived separately
 * in simd.h and simd4.h.
 *
 * \ingroup module_simd
 */
class SimdBaseTest : public ::testing::Test
{
public:
    /*! \brief Return the default ulp tolerance for current precision
     */
    static constexpr std::int64_t defaultRealUlpTol()
    {
        return (1LL << (2 + std::numeric_limits<real>::digits - GMX_SIMD_ACCURACY_BITS_REAL));
    }

    /*! \brief Initialize new SIMD test fixture with default tolerances.
     *
     * The default absolute tolerance is set to 0, which means the we always
     * check the ulp tolerance by default (passing the absolute tolerance
     * test would otherwise mean we approve the test instantly).
     *
     * The default ulp tolerance is set based on the target number of
     * bits requested for single or double precision, depending on what
     * the default Gromacs precision is. We add two bits to avoid
     * tests failing due to corner cases where compiler optimization might
     * cause a slight precision loss e.g. for very small numbers.
     *
     * Most SIMD math functions actually achieve 2-3 ulp accuracy in single,
     * but by being a bit liberal we only catch real errors rather than
     * doing compiler-standard-compliance debugging.
     *
     * The range is used by derived classes to test math functions. The
     * default test range will be [1,10], which is intentionally
     * conservative so it works with (inverse) square root, division,
     * exponentials, logarithms, and error functions.
     */
    SimdBaseTest() : ulpTol_(defaultRealUlpTol()), absTol_(0) {}

    /*! \brief Adjust ulp tolerance from the default 10 (float) or 255 (double). */
    void setUlpTol(std::int64_t newTol) { ulpTol_ = newTol; }

    /*! \brief Adjust ulp tolerance for single accuracy functions. */
    void setUlpTolSingleAccuracy(std::int64_t newTol)
    {
        const int realBits   = std::numeric_limits<real>::digits;
        const int singleBits = std::numeric_limits<float>::digits;
        // In single precision the expression (1LL << 0) evaluates to 1.
        setUlpTol(newTol * (1LL << (realBits - singleBits)));
    }

    /*! \brief Adjust the absolute tolerance from the default 0.
     *
     * If values are closer than the absolute tolerance, the test will pass
     * no matter what their ulp difference is.
     */
    void setAbsTol(real newTol) { absTol_ = newTol; }

    /*! \brief Number of test points to use, settable on command line.
     *
     * \note While this has to be a static non-const variable for the
     *       command-line option to work, you should never change it
     *       manually in any of the tests, because the static storage
     *       class will make the value apply to all subsequent tests
     *       unless you remember to reset it.
     */
    static int s_nPoints;

    /*! \brief Compare two std::vector<real> for approximate equality.
     *
     * This is an internal implementation routine that will be used by
     * routines in derived child classes that first convert SIMD or SIMD4
     * variables to std::vector<real>. Do not call it directly.
     *
     * This routine is designed according to the Google test specs, so the char
     * strings will describe the arguments to the macro.
     *
     * The comparison is applied to each element, and it returns true if each element
     * in the vector test variable is within the class tolerances of the corresponding
     * reference elements.
     */
    ::testing::AssertionResult compareVectorRealUlp(const char*              refExpr,
                                                    const char*              tstExpr,
                                                    const std::vector<real>& ref,
                                                    const std::vector<real>& tst) const;

    /*! \brief Compare std::vectors for exact equality.
     *
     * The template in this class makes it usable for testing both
     * SIMD floating-point and integers variables, after conversion to
     * vectors.
     * This is an internal implementation routine that will be used by
     * routines in derived child classes that first convert SIMD or SIMD4
     * variables to std::vector<real>. Do not call it directly.
     *
     * This routine is designed according to the Google test specs, so the char
     * strings will describe the arguments to the macro.
     *
     * The comparison is applied to each element, and it returns true if each element
     * in the vector test variable is within the class tolerances of the corresponding
     * reference elements.
     */
    template<typename T>
    ::testing::AssertionResult compareVectorEq(const char*           refExpr,
                                               const char*           tstExpr,
                                               const std::vector<T>& ref,
                                               const std::vector<T>& tst)
    {
        if (ref == tst)
        {
            return ::testing::AssertionSuccess();
        }
        else
        {
            return ::testing::AssertionFailure()
                   << "Failing SIMD comparison between " << refExpr << " and " << tstExpr << std::endl
                   << "Ref. values: " << ::testing::PrintToString(ref) << std::endl
                   << "Test values: " << ::testing::PrintToString(tst) << std::endl;
        }
    }

protected:
    std::int64_t ulpTol_; //!< Current tolerance in units-in-last-position.
    real         absTol_; //!< Current absolute tolerance.
};

} // namespace test
} // namespace gmx

#endif // GMX_SIMD_TESTS_BASE_H
