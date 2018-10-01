/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Tests for gmx::Regex
 *
 * These tests ensure that basic regex operations work. We have
 * two different underlying implementations, so we need to prove
 * to ourselves that these work the same on the range of operations
 * we try to use.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "gromacs/utility/gmxregex.h"

#include <future>
#include <gtest/gtest.h>


namespace gmx
{
namespace test
{
namespace
{

TEST(RegexBasicTest, BasicMatchesWorkWhenSupported)
{
    if (!Regex::isSupported())
    {
        return;
    }

    EXPECT_TRUE(regexMatch("dog", Regex("dog")));
    EXPECT_TRUE(regexMatch("dog", Regex("^dog")));
    EXPECT_TRUE(regexMatch("dog", Regex("dog$")));
    EXPECT_TRUE(regexMatch("dog", Regex("^dog$")));
    EXPECT_TRUE(regexMatch("dog", Regex(".*d.*")));
    EXPECT_TRUE(regexMatch("dog", Regex("d+og")));

    EXPECT_FALSE(regexMatch("dog", Regex("cat")));
    EXPECT_FALSE(regexMatch("dog", Regex("^og")));
    EXPECT_FALSE(regexMatch("dog", Regex("do$")));
    EXPECT_FALSE(regexMatch("dog", Regex("^o$")));
    EXPECT_FALSE(regexMatch("dog", Regex(".*c.*")));
    EXPECT_FALSE(regexMatch("dog", Regex("c+dog")));
}

TEST(RegexBasicTest, MatchesForCharacterClassesWorkWhenSupported)
{
    if (!Regex::isSupported())
    {
        return;
    }

    EXPECT_TRUE(regexMatch("Dog", Regex("[[:alpha:]]+")));
    EXPECT_TRUE(regexMatch("dog", Regex("[[:lower:]]+")));
    EXPECT_TRUE(regexMatch("DOG", Regex("[[:upper:]]+")));
    EXPECT_TRUE(regexMatch("123", Regex("[[:digit:]]+")));
    EXPECT_TRUE(regexMatch("123aE", Regex("[[:xdigit:]]+")));
    EXPECT_TRUE(regexMatch("Dog123", Regex("[[:alnum:]]+")));
    EXPECT_TRUE(regexMatch(" ", Regex("[[:space:]]+")));
    EXPECT_TRUE(regexMatch("\t ", Regex("[[:blank:]]+")));
    EXPECT_TRUE(regexMatch(".,:;-", Regex("[[:punct:]]+")));
    EXPECT_TRUE(regexMatch("Dog,123", Regex("[[:graph:]]+")));
    EXPECT_TRUE(regexMatch("Dog, 123", Regex("[[:print:]]+")));

    EXPECT_FALSE(regexMatch("D0g", Regex("[[:alpha:]]+")));
    EXPECT_FALSE(regexMatch("Dog", Regex("[[:lower:]]+")));
    EXPECT_FALSE(regexMatch("Dog", Regex("[[:upper:]]+")));
    EXPECT_FALSE(regexMatch("D0g", Regex("[[:digit:]]+")));
    EXPECT_FALSE(regexMatch("12xyz3aE", Regex("[[:xdigit:]]+")));
    EXPECT_FALSE(regexMatch("Dog, 123", Regex("[[:alnum:]]+")));
    EXPECT_FALSE(regexMatch("dog", Regex("[[:space:]]+")));
    EXPECT_FALSE(regexMatch("dog\t ", Regex("[[:blank:]]+")));
    EXPECT_FALSE(regexMatch(".,:d;-", Regex("[[:punct:]]+")));
    EXPECT_FALSE(regexMatch("Dog,123 ", Regex("[[:graph:]]+")));

    EXPECT_TRUE(regexMatch(" ", Regex("[[:space:]]")));
    EXPECT_TRUE(regexMatch(" ", Regex("[[:blank:]]")));
    EXPECT_TRUE(regexMatch("\t", Regex("[[:blank:]]")));
    EXPECT_TRUE(regexMatch("\t", Regex("[[:space:]]")));
}

} // namespace
} // namespace test
} // namespace gmx
