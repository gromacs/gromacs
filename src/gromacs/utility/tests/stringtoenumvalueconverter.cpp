/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Tests for string-to-enum-value helper functor
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/stringtoenumvalueconverter.h"

#include <optional>
#include <ostream>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/stringcompare.h"

namespace gmx
{
namespace test
{
namespace
{

//! Type to use in testing
enum class Foo : int
{
    Bar,
    Ugh,
    Fooz,
    Weird,
    Empty,
    Count,
    Default = Ugh
};

//! Declare the conventional conversion function
const char* enumValueToString(Foo f)
{
    static constexpr gmx::EnumerationArray<Foo, const char*> toString = {
        "Bar", "Ugh ", "Foo z", "We-i_rd", ""
    };
    return toString[f];
}

//! Declare an atypical conversion function
const char* enumValueToLetterAsString(Foo f)
{
    static constexpr gmx::EnumerationArray<Foo, const char*> toString = { "B", "U", "F", "W", "" };
    return toString[f];
}

//! Pretty-printer for GoogleTest to use
::std::ostream& operator<<(::std::ostream& os, const std::optional<Foo>& f)
{
    if (f)
    {
        return os << enumValueToString(f.value());
    }
    else
    {
        return os << "Out-of-range Foo value";
    }
}

TEST(StringToEnumValueConverterTest, ExactStringComparisonWorksWithoutStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::Exact, StripStrings::No> converter;
    EXPECT_EQ(converter.valueFrom("Bar"), Foo::Bar);
    EXPECT_FALSE(converter.valueFrom("bar"));
    EXPECT_EQ(converter.valueFrom("Ugh "), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("ugh "));
    EXPECT_EQ(converter.valueFrom("Foo z"), Foo::Fooz);
    EXPECT_FALSE(converter.valueFrom("foo z"));
    EXPECT_EQ(converter.valueFrom("We-i_rd"), Foo::Weird);
    EXPECT_FALSE(converter.valueFrom("we-i_rd"));
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("count"));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("BarFoo z"));
    EXPECT_FALSE(converter.valueFrom("Ugh"));
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, CaseInsensitiveStringComparisonWorksWithoutStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::CaseInsensitive, StripStrings::No> converter;
    EXPECT_EQ(converter.valueFrom("Bar"), Foo::Bar);
    EXPECT_EQ(converter.valueFrom("bar"), Foo::Bar);
    EXPECT_EQ(converter.valueFrom("Ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("Foo z"), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("foo z"), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("We-i_rd"), Foo::Weird);
    EXPECT_EQ(converter.valueFrom("we-i_rd"), Foo::Weird);
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("Count"));
    EXPECT_FALSE(converter.valueFrom("count"));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("barfoo z"));
    EXPECT_FALSE(converter.valueFrom("Ugh"));
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, CaseAndDashInsensitiveStringComparisonWorksWithoutStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::CaseAndDashInsensitive, StripStrings::No> converter;
    EXPECT_EQ(converter.valueFrom("Bar"), Foo::Bar);
    EXPECT_EQ(converter.valueFrom("b-ar"), Foo::Bar);
    EXPECT_EQ(converter.valueFrom("Ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("_ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("Foo z"), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("fo_o z"), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("We-i_rd"), Foo::Weird);
    EXPECT_EQ(converter.valueFrom("we-i_rd"), Foo::Weird);
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("Count"));
    EXPECT_FALSE(converter.valueFrom("count"));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("Bar-Foo z"));
    EXPECT_FALSE(converter.valueFrom("Ugh"));
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, ExactStringComparisonWorksWithStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::Exact, StripStrings::Yes> converter;
    EXPECT_EQ(converter.valueFrom("Bar "), Foo::Bar);
    EXPECT_FALSE(converter.valueFrom("Ba r"));
    EXPECT_EQ(converter.valueFrom("Ugh"), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("ugh"));
    EXPECT_EQ(converter.valueFrom("  Foo z "), Foo::Fooz);
    EXPECT_FALSE(converter.valueFrom(" foo z"));
    EXPECT_EQ(converter.valueFrom("We-i_rd  "), Foo::Weird);
    EXPECT_FALSE(converter.valueFrom("  we-i_rd "));
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom(" Count"));
    EXPECT_FALSE(converter.valueFrom("count  "));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("Bar Foo z"));
    EXPECT_EQ(converter.valueFrom(" Ugh"), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, CaseInsensitiveStringComparisonWorksWithStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::CaseInsensitive, StripStrings::Yes> converter;
    EXPECT_EQ(converter.valueFrom("Bar "), Foo::Bar);
    EXPECT_FALSE(converter.valueFrom("ba r"));
    EXPECT_EQ(converter.valueFrom("Ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("  ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom("  Foo z "), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("foo z   "), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("We-i_rd  "), Foo::Weird);
    EXPECT_EQ(converter.valueFrom("  we-i_rd "), Foo::Weird);
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("  Count"));
    EXPECT_FALSE(converter.valueFrom("count "));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("  barfoo z"));
    EXPECT_EQ(converter.valueFrom(" Ugh"), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, CaseAndDashInsensitiveStringComparisonWorksWithStripping)
{
    StringToEnumValueConverter<Foo, enumValueToString, StringCompareType::CaseAndDashInsensitive, StripStrings::Yes> converter;
    EXPECT_EQ(converter.valueFrom("Bar "), Foo::Bar);
    EXPECT_FALSE(converter.valueFrom("b-a r"));
    EXPECT_EQ(converter.valueFrom("Ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom(" _ugh "), Foo::Ugh);
    EXPECT_EQ(converter.valueFrom(" Foo z "), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("fo_o z  "), Foo::Fooz);
    EXPECT_EQ(converter.valueFrom("We-i_rd  "), Foo::Weird);
    EXPECT_EQ(converter.valueFrom("  we-i_rd "), Foo::Weird);
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("  Count"));
    EXPECT_FALSE(converter.valueFrom("coun-t "));
    EXPECT_FALSE(converter.valueFrom("Unknown"));
    EXPECT_FALSE(converter.valueFrom("  Bar-Foo z   "));
    EXPECT_EQ(converter.valueFrom("Ugh  "), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("Default"));
}

TEST(StringToEnumValueConverterTest, CustomConverterWorks)
{
    StringToEnumValueConverter<Foo, enumValueToLetterAsString> converter;
    EXPECT_EQ(converter.valueFrom("B"), Foo::Bar);
    EXPECT_FALSE(converter.valueFrom("b"));
    EXPECT_EQ(converter.valueFrom("U"), Foo::Ugh);
    EXPECT_FALSE(converter.valueFrom("u"));
    EXPECT_EQ(converter.valueFrom("F"), Foo::Fooz);
    EXPECT_FALSE(converter.valueFrom("f"));
    EXPECT_EQ(converter.valueFrom(""), Foo::Empty);
    EXPECT_FALSE(converter.valueFrom("C"));
    EXPECT_FALSE(converter.valueFrom("c"));
    EXPECT_FALSE(converter.valueFrom("X"));
    EXPECT_FALSE(converter.valueFrom("Default"));
}

} // namespace
} // namespace test
} // namespace gmx
