/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * \brief
 * Tests for Regex functionality.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/gmxregex.h"

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace
{

/********************************************************************
 * Tests for simple regular expressions.
 */

struct RegexTestHelper
{
    const char *regex_;     // The regular expression to test
    const char *match_;     // A string that will match the regular expression
};

/*! \brief Array of regex and matches to test
 *
 * Only tests of flavours of regular expressions actually used should
 * go here. The purpose is to check that the standard library in use
 * has enough actual support for our needs, since the coverage and
 * quality of the implementations varies a bit. */
RegexTestHelper   regexArray[] = {
    { "test",            "test" },
    { "[BD]",            "B" },
    { ".*",              "ssc139)@#" },
    { "[[:space:]]",     " " },
    { "[[:space:]]{1,}", " \t\n  " },
    { "[[:alpha:]]+",    "sdfosdsc" },
    { "[[:digit:]]*",    "121239" },
};

class RegexTest : public ::testing::TestWithParam<RegexTestHelper>
{
};

TEST_P(RegexTest, CompilesAndMatches)
{
    auto      &regexToTest = GetParam();
    gmx::Regex r;
    EXPECT_NO_THROW_GMX(r = gmx::Regex(regexToTest.regex_));
    EXPECT_TRUE(gmx::regexMatch(regexToTest.match_, r)) << "with regex '" << regexToTest.regex_ << "' and match '" << regexToTest.match_ << "'";
}

INSTANTIATE_TEST_CASE_P(WithInputs, RegexTest, ::testing::ValuesIn(std::begin(regexArray), std::end(regexArray)));

} // namespace
