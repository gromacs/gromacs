/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Tests gmx::CommandLineParser.
 *
 * These tests exercise a large fraction of the code, so they may
 * catch errors in other parts than just in command-line parsing.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

namespace
{

using gmx::test::CommandLine;

class CommandLineParserTest : public ::testing::Test
{
    public:
        CommandLineParserTest();

        gmx::Options            options_;
        gmx::CommandLineParser  parser_;
        bool                    flag_;
        std::vector<int>        ivalues_;
        std::vector<double>     dvalues_;
};

CommandLineParserTest::CommandLineParserTest()
    : options_(NULL, NULL), parser_(&options_),
      flag_(false)
{
    using gmx::BooleanOption;
    using gmx::IntegerOption;
    using gmx::DoubleOption;
    options_.addOption(BooleanOption("flag").store(&flag_));
    options_.addOption(IntegerOption("mvi").storeVector(&ivalues_).multiValue());
    options_.addOption(DoubleOption("mvd").storeVector(&dvalues_).allowMultiple());
}

TEST_F(CommandLineParserTest, HandlesSingleValues)
{
    const char *const cmdline[] = {
        "test", "-flag", "yes", "-mvi", "2", "-mvd", "2.7"
    };
    CommandLine       args(CommandLine::create(cmdline));
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());

    EXPECT_TRUE(flag_);
    ASSERT_EQ(1U, ivalues_.size());
    EXPECT_EQ(2, ivalues_[0]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(2.7, dvalues_[0]);
}

TEST_F(CommandLineParserTest, HandlesNegativeNumbers)
{
    const char *const cmdline[] = {
        "test", "-mvi", "1", "-2", "-mvd", "-2.7"
    };
    CommandLine       args(CommandLine::create(cmdline));
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());

    ASSERT_EQ(2U, ivalues_.size());
    EXPECT_EQ(1, ivalues_[0]);
    EXPECT_EQ(-2, ivalues_[1]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(-2.7, dvalues_[0]);
}

} // namespace
