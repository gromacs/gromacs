/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Tests test-naming functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/naming.h"

#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/utility/strconvert.h"

namespace gmx
{
namespace test
{
namespace
{

//! Example for use in testing
enum class Flavor : int
{
    A,
    B,
    C,
    Count
};
//! Formatter for Flavor
EnumerationArray<Flavor, const char*> sc_flavorName = { "A", "B", "C" };
//! Format function for Flavor
std::string enumValueToString(const Flavor f)
{
    return std::string("Flavor_") + sc_flavorName[f];
}

//! Format functor
class FormatFunctorForInt
{
public:
    //! The call operator
    std::string operator()(const int v) { return intToString(v); }
};

TEST(NameOfTestFromTupleTest, WorksWithEmptyTuple)
{
    using TestParameters = std::tuple<>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple() };
    EXPECT_EQ("", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithFormatFunction)
{
    using TestParameters = std::tuple<int>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(intToString) };
    EXPECT_EQ("3", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(3), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithFormatFunctionOfEnumVariable)
{
    using TestParameters = std::tuple<Flavor>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(enumValueToString) };
    EXPECT_EQ("Flavor_A", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(Flavor::A), 0)));
}

TEST(NameOfTestFromTupleTest, RejectsNullptrFormatFunction)
{
    using TestParameters = std::tuple<int>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(nullptr) };
    EXPECT_THROW_GMX(namer(testing::TestParamInfo<TestParameters>(std::make_tuple(3), 0)), APIError);
}

TEST(NameOfTestFromTupleTest, WorksWithFormatLambda)
{
    using TestParameters = std::tuple<int>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple([](int /* a */) { return "foo"; }) };
    EXPECT_EQ("foo", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(3), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithUseStringFormat)
{
    using TestParameters = std::tuple<std::string>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(useString) };
    EXPECT_EQ("foo_bar", namer(testing::TestParamInfo<TestParameters>(std::make_tuple("foo/bar"), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithPrefixFormatter)
{
    using TestParameters = std::tuple<double>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(
            PrefixFormatter<double, toString>{ "pi_" }) };
    EXPECT_EQ("pi_3_14", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(3.14), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithFormatFunctor)
{
    using TestParameters = std::tuple<int>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(FormatFunctorForInt{}) };
    EXPECT_EQ("4", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(4), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithFormatFromEnumerationArray)
{
    using TestParameters = std::tuple<Flavor>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(sc_flavorName) };
    EXPECT_EQ("B", namer(testing::TestParamInfo<TestParameters>(std::make_tuple(Flavor::B), 0)));
}

TEST(NameOfTestFromTupleTest, WorksWithMixtureOfFormatters)
{
    using TestParameters = std::tuple<Flavor, int, double, int>;
    NameOfTestFromTuple<TestParameters> namer{ std::make_tuple(
            sc_flavorName, intToString, doubleToString, FormatFunctorForInt{}) };
    EXPECT_EQ("B_2_3_2_9",
              namer(testing::TestParamInfo<TestParameters>(std::make_tuple(Flavor::B, 2, 3.2, 9), 0)));
}

TEST(RefDataFilenameMakerTest, WorksWithFormatFunction)
{
    using TestParameters = std::tuple<int>;
    RefDataFilenameMaker<TestParameters> maker{ std::make_tuple(intToString) };
    EXPECT_EQ("RefDataFilenameMakerTest_3.xml", maker(std::make_tuple(3)));
}

TEST(RefDataFilenameMakerTest, WorksWithMixtureOfFormatters)
{
    using TestParameters = std::tuple<Flavor, int, double, int>;
    RefDataFilenameMaker<TestParameters> maker{ std::make_tuple(
            sc_flavorName, intToString, doubleToString, FormatFunctorForInt{}) };
    EXPECT_EQ("RefDataFilenameMakerTest_B_2_3_2_9.xml", maker(std::make_tuple(Flavor::B, 2, 3.2, 9)));
}

TEST(RefDataFilenameMakerTest, WorksWithToEmpty)
{
    using TestParameters = std::tuple<int>;
    RefDataFilenameMaker<TestParameters> maker{ std::make_tuple(toEmptyString<int>) };
    EXPECT_EQ("RefDataFilenameMakerTest.xml", maker(std::make_tuple(3)));
}

} // namespace
} // namespace test
} // namespace gmx
