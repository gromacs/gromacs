/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

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
        SimdBaseTest() :
#ifdef GMX_DOUBLE
            ulpTol_((1LL << (2 + std::numeric_limits<double>::digits-GMX_SIMD_ACCURACY_BITS_DOUBLE))),
#else
            ulpTol_((1LL << (2 + std::numeric_limits<float>::digits-GMX_SIMD_ACCURACY_BITS_SINGLE))),
#endif
            absTol_(0), range_(std::pair<real, real>(1, 10))
        {
        }

        /*! \brief Adjust ulp tolerance from the default 10 (float) or 255 (double). */
        void setUlpTol(gmx_int64_t newTol)   { ulpTol_ = newTol; }

        /*! \brief Adjust the absolute tolerance from the default 0.
         *
         * If values are closer than the absolute tolerance, the test will pass
         * no matter what their ulp difference is.
         */
        void setAbsTol(real newTol)          { absTol_ = newTol; }

        /*! \brief Change math function testing range from the default [1,10]. */
        void setRange(real low, real high) { range_.first = low; range_.second = high; }

        static int  s_nPoints;    //!< Number of test points to use, settable on command line.

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
        ::testing::AssertionResult
        compareVectorRealUlp(const char * refExpr,  const char * tstExpr,
                             const std::vector<real> &ref, const std::vector<real> &tst);

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
        template <typename T> ::testing::AssertionResult
        compareVectorEq(const char * refExpr,  const char * tstExpr,
                        const std::vector<T> &ref, const std::vector<T> &tst)
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
        gmx_int64_t            ulpTol_;       //!< Current tolerance in units-in-last-position.
        real                   absTol_;       //!< Current absolute tolerance.
        std::pair<real, real>  range_;        //!< Range for math function tests.
};

}      // namespace
}      // namespace

#endif // GMX_SIMD_TESTS_BASE_H
