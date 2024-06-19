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
 * \brief
 * Tests for C-style string utility functions.
 *
 * For development, the tests can be run with a '-stdout' command-line option
 * to print out the help to stdout instead of using the XML reference
 * framework.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/cstringutil.h"

#include <cstdlib>
#include <cstring>

#include <string>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{
namespace
{

/********************************************************************
 * Tests for simple c-string utilities
 */

TEST(CStringUtilityTest, CaseInsensitiveComparison)
{
    EXPECT_TRUE(gmx_strcasecmp("foo", "foo") == 0);
    EXPECT_TRUE(gmx_strcasecmp("foo", "foO") == 0);
    EXPECT_FALSE(gmx_strcasecmp("foo", "bar") == 0);

    EXPECT_FALSE(gmx_strcasecmp("foobar", "foo") == 0);
    EXPECT_FALSE(gmx_strcasecmp("foo", "foobar") == 0);
    // This is something that should not work, but the old code allows it.
    EXPECT_TRUE(gmx_strcasecmp("foo", "foo\0bar") == 0);
}

TEST(CStringUtilityTest, CaseInsensitiveComparisonInLength)
{
    // Setting the length to 0 is an automatic pass for the comparison.
    EXPECT_TRUE(gmx_strncasecmp("foo", "foo", 0) == 0);
    EXPECT_TRUE(gmx_strncasecmp("foo", "bar", 0) == 0);
    for (int i = 1; i < 6; i++)
    {
        // The code always checks the length of the first string or up to the user specified length.
        EXPECT_TRUE(gmx_strncasecmp("foo", "foo", i) == 0) << "Fails for number" << i;
        EXPECT_TRUE(gmx_strncasecmp("foo", "foO", i) == 0) << "Fails for number" << i;
    }
    for (int i = 1; i < 4; i++)
    {
        EXPECT_FALSE(gmx_strncasecmp("foo", "bar", i) == 0) << "Fails for number" << i;
    }

    EXPECT_TRUE(gmx_strncasecmp("foo", "foobar", 3) == 0);
    EXPECT_TRUE(gmx_strncasecmp("foo", "foOBar", 3) == 0);
    EXPECT_FALSE(gmx_strncasecmp("foofoo", "foobar", 4) == 0);
    EXPECT_FALSE(gmx_strncasecmp("foo", "foobar", 4) == 0);
    EXPECT_FALSE(gmx_strncasecmp("foobar", "foo", 4) == 0);
}

template<typename F>
void testInplace(F f, const char* input, const char* expectedOutput)
{
    char* tmp = strdup(input);
    f(tmp);
    EXPECT_STREQ(tmp, expectedOutput);
    free(tmp);
}

TEST(CStringUtilityTest, strip_comment)
{
    EXPECT_NO_THROW(strip_comment(nullptr));

    testInplace(strip_comment, "", "");
    testInplace(strip_comment, "foo", "foo");
    testInplace(strip_comment, "foo;", "foo");
    testInplace(strip_comment, "foo;bar", "foo");
    testInplace(strip_comment, ";foo", "");
    testInplace(strip_comment, "foo;bar;baz", "foo");
    testInplace(strip_comment, " foo ; bar", " foo ");
}

TEST(CStringUtilityTest, upstring)
{
    EXPECT_NO_THROW(upstring(nullptr));

    testInplace(upstring, "", "");
    testInplace(upstring, "foo", "FOO");
    testInplace(upstring, "FOO", "FOO");
    testInplace(upstring, "123", "123");
    testInplace(upstring, "aBcDeF\n123", "ABCDEF\n123");
}

TEST(CStringUtilityTest, ltrim)
{
    EXPECT_NO_THROW(ltrim(nullptr));

    testInplace(ltrim, "", "");
    testInplace(ltrim, "    ", "");
    testInplace(ltrim, "   \t ", "");
    testInplace(ltrim, "foo", "foo");
    testInplace(ltrim, "   foo  ", "foo  ");
    testInplace(ltrim, "\tfoo\t", "foo\t");
    testInplace(ltrim, " foo bar ", "foo bar ");
}

TEST(CStringUtilityTest, rtrim)
{
    EXPECT_NO_THROW(rtrim(nullptr));

    testInplace(rtrim, "", "");
    testInplace(rtrim, "    ", "");
    testInplace(rtrim, "   \t ", "");
    testInplace(rtrim, "foo", "foo");
    testInplace(rtrim, "   foo  ", "   foo");
    testInplace(rtrim, "\tfoo\t", "\tfoo");
    testInplace(rtrim, " foo bar ", " foo bar");
}

TEST(CStringUtilityTest, trim)
{
    EXPECT_NO_THROW(trim(nullptr));

    testInplace(trim, "", "");
    testInplace(trim, "    ", "");
    testInplace(trim, "   \t ", "");
    testInplace(trim, "foo", "foo");
    testInplace(trim, "   foo  ", "foo");
    testInplace(trim, "\tfoo\t", "foo");
    testInplace(trim, " foo bar ", "foo bar");
}

} // namespace
} // namespace test
} // namespace gmx
