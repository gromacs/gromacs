/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Tests utilities for routines that parse fields e.g. from grompp input
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "../strdb.h"

#include <stdio.h>

#include <gtest/gtest.h>

#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
namespace testing
{

class get_a_lineTest : public ::testing::Test
{
    public:
        get_a_lineTest() {}

        void runTest(const char *input, const char *expectedOutput)
        {
            std::string lineToParse(input), output(input);

            char *buffer = const_cast<char *>(lineToParse.c_str());
            // Trust fmemopen not to break buffer
            FILE *memstream = fmemopen(buffer, lineToParse.size()+1, "r");

            buffer = const_cast<char *>(output.c_str());
            // Trust get_a_line not to break buffer
            get_a_line(memstream, buffer, output.size()+1);
            fclose(memstream);

            ASSERT_STREQ(expectedOutput, output.c_str());
        }
};

TEST_F(get_a_lineTest, ParsesValidSingleLines)
{
    runTest("mdp-option = 1\n", "mdp-option = 1");
    runTest("mdp-option = 1", "mdp-option = 1");
    runTest("  mdp-option = 1\n", "mdp-option = 1");
    runTest("  mdp-option = 1", "mdp-option = 1");
    runTest("mdp-option = 1;comment\n", "mdp-option = 1");
    runTest("mdp-option = 1; comment\n", "mdp-option = 1");
    runTest("mdp-option = 1  ; comment\n", "mdp-option = 1");
    runTest("  mdp-option = 1  ; comment\n", "mdp-option = 1");
    runTest(";mdp-option = 1  ; comment\n", "");
    runTest("", "");
    runTest("   ", "   "); // Can't skip this line because there is no next line
    runTest("\n", "");
}

TEST_F(get_a_lineTest, ParsesValidMultipleLines)
{
    runTest("mdp-option = 1\nmdp-option = 2\n", "mdp-option = 1");
    runTest("\nmdp-option = 1\n", "mdp-option = 1");
    runTest("   \nmdp-option = 1\n", "mdp-option = 1");
}

} // namespace
} // namespace
