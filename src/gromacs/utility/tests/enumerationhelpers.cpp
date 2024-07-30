/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Tests for enumeration helpers
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/enumerationhelpers.h"

#include <cstddef>

#include <iostream>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Type to use in testing
enum class Foo
{
    Bar,
    Baz,
    Fooz,
    Count
};

TEST(EnumerationHelpersTest, EnumerationWrapperWorks)
{
    EnumerationWrapper<Foo> iter;

    // Range-based for works
    int i = 0;
    for (Foo c : iter)
    {
        EXPECT_EQ(static_cast<int>(c), i++);
    }

    // Normal iterators work
    i = 0;
    for (auto c = iter.begin(); c != iter.end(); ++c)
    {
        EXPECT_EQ(static_cast<int>(*c), i++);
    }

    auto a = std::begin(iter);
    auto b = std::begin(iter);

    ASSERT_EQ(a, b);
    ASSERT_EQ(*(a++), Foo::Bar);
    ASSERT_EQ(*(++b), Foo::Baz);
}

TEST(EnumerationHelpersTest, EnumerationArrayWorks)
{
    using FooArray = EnumerationArray<Foo, std::string>;
    const FooArray fooStrings{ { "Bar", "Baz", "Fooz" } };

    // Keys give you the constants associated with each array index.
    int i = 0;
    for (auto k : FooArray::keys())
    {
        EXPECT_EQ(static_cast<int>(k), i++);
    }

    // Keys give you the constants associated with each array index.
    i = 0;
    for (auto k : keysOf(fooStrings))
    {
        EXPECT_EQ(static_cast<int>(k), i++);
    }

    // Using iterators and operator[] gives the array values.
    i = 0;
    for (const auto& s : fooStrings)
    {
        EXPECT_EQ(s, fooStrings[i++]);
    }

    // Using reverse iterators gives the array values.
    i = 2;
    for (auto s = fooStrings.rbegin(); s != fooStrings.rend(); ++s)
    {
        EXPECT_EQ((*s), fooStrings[i--]);
    }

    // Incrementing iterators works
    const auto* x = std::begin(fooStrings);
    EXPECT_EQ(*x, "Bar");
    ++x;
    EXPECT_EQ(*x, "Baz");
    ++x;
    EXPECT_EQ(*x, "Fooz");

    // Operator[] can be used with enumeration values.
    EXPECT_EQ(fooStrings[Foo::Bar], "Bar");
    EXPECT_EQ(fooStrings[Foo::Baz], "Baz");
    EXPECT_EQ(fooStrings[Foo::Fooz], "Fooz");
}

TEST(EnumerationHelpersTest, EnumerationArrayCountIsSafe)
{
    using FooArray = EnumerationArray<Foo, std::string>;
    const FooArray fooStrings{ { "Bar", "Baz", "Fooz" } };

    // Ensures that the assertions in EnumerationArray::operator[]
    // would fire if an out-range value (including Count) was used.
    EXPECT_LE(fooStrings.size(), size_t(Foo::Count));
#ifndef NDEBUG
    // Tests (where possible) that those assertions do fire in a build
    // with debug behavior.
    GMX_EXPECT_DEATH_IF_SUPPORTED(fooStrings[Foo::Count], "index out of range");
#endif
}

//! Helper function
void func(ArrayRef<const int> a)
{
    EXPECT_EQ(3, a[1]);
}

TEST(EnumerationHelpersTest, ArrayRefOfEnumerationArrayWorks)
{
    using FooArray = EnumerationArray<Foo, int>;

    FooArray counts = { { 2, 3, 1 } };

    // Test that explicit conversion works
    ArrayRef<const int> arrayRef(counts);
    EXPECT_EQ(3, arrayRef[1]);

    // Test that implicit conversion works
    func(counts);

    // Note that ArrayRef<int> arrayRef(counts) does not compile, as
    // expected, but we can't test for that.
}

} // namespace
} // namespace test
} // namespace gmx
