/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

/*! \internal \file
 * \brief Implements floating-point matchers from testmatchers.h for
 * use with Google Mock.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testmatchers.h"

#include <memory>
#include <ostream>
#include <string>
#include <tuple>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/*! \brief Implementation class for RealEq matcher.
 *
 * See RealEq().
 *
 * The implementation is templated so that we can support all of real,
 * float and double in the same build without duplication.
 */
template<typename FloatType>
class FloatTypeMatcher : public testing::MatcherInterface<std::tuple<FloatType, FloatType>>
{
public:
    //! Constructor
    FloatTypeMatcher(const FloatingPointTolerance& tolerance) : tolerance_(tolerance) {}
    //! Compare the two elements of \c arg, return whether they are equal, and comment on \c listener when they are not.
    bool MatchAndExplain(std::tuple<FloatType, FloatType> arg, testing::MatchResultListener* listener) const override
    {
        const FloatType&        value1 = std::get<0>(arg);
        const FloatType&        value2 = std::get<1>(arg);
        FloatingPointDifference diff(value1, value2);
        if (tolerance_.isWithin(diff))
        {
            return true;
        }
        if (listener->stream())
        {
            *listener->stream() << "  Actual value: " << value2 << std::endl
                                << "Expected value: " << value1 << std::endl
                                << "    Difference: " << diff.toString() << std::endl
                                << "     Tolerance: " << tolerance_.toString(diff);
        }

        return false;
    }
    //! Describe to a human what matching means.
    void DescribeTo(::std::ostream* os) const override { *os << "matches within tolerance"; }
    //! Describe to a human what failing to match means.
    void DescribeNegationTo(::std::ostream* os) const override
    {
        *os << "does not match within tolerance";
    }

private:
    //! Tolerance used in matching
    FloatingPointTolerance tolerance_;
};

testing::Matcher<std::tuple<float, float>> FloatEq(const FloatingPointTolerance& tolerance)
{
    return testing::MakeMatcher(new FloatTypeMatcher<float>(tolerance));
}

testing::Matcher<std::tuple<double, double>> DoubleEq(const FloatingPointTolerance& tolerance)
{
    return testing::MakeMatcher(new FloatTypeMatcher<double>(tolerance));
}

testing::Matcher<std::tuple<real, real>> RealEq(const FloatingPointTolerance& tolerance)
{
    return testing::MakeMatcher(new FloatTypeMatcher<real>(tolerance));
}

/*! \brief Implementation class for RvecEq matcher
 *
 * See RvecEq().
 */
template<typename FloatType>
class RVecMatcher :
    public testing::MatcherInterface<std::tuple<BasicVector<FloatType>, BasicVector<FloatType>>>
{
public:
    //! Convenience type
    using VectorType = BasicVector<FloatType>;
    //! Constructor
    RVecMatcher(const FloatingPointTolerance& tolerance) : tolerance_(tolerance) {}
    //! Compare the two elements of \c arg, return whether they are equal, and comment on \c listener when they are not.
    bool MatchAndExplain(std::tuple<VectorType, VectorType> arg,
                         testing::MatchResultListener*      listener) const override
    {
        const VectorType&           lhs = std::get<0>(arg);
        const VectorType&           rhs = std::get<1>(arg);
        FloatTypeMatcher<FloatType> floatTypeMatcher(tolerance_);
        bool                        matches = true;
        for (int d = 0; d < DIM; ++d)
        {
            auto floatTuple = std::make_tuple<FloatType, FloatType>(lhs[d], rhs[d]);
            matches         = matches && floatTypeMatcher.MatchAndExplain(floatTuple, listener);
        }
        return matches;
    }
    //! Describe to a human what matching means.
    void DescribeTo(::std::ostream* os) const override
    {
        *os << "matches all elements within tolerance";
    }
    //! Describe to a human what failing to match means.
    void DescribeNegationTo(::std::ostream* os) const override
    {
        *os << "does not match all elements within tolerance";
    }

private:
    //! Tolerance used in matching
    FloatingPointTolerance tolerance_;
};

// Currently there's no need for explicit float or double flavours of
// RVec comparison, but those would be simple to add later.
testing::Matcher<std::tuple<RVec, RVec>> RVecEq(const FloatingPointTolerance& tolerance)
{
    return testing::MakeMatcher(new RVecMatcher<real>(tolerance));
}

} // namespace test
} // namespace gmx
