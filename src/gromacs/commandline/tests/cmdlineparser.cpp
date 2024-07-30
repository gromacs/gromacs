/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Tests gmx::CommandLineParser.
 *
 * These tests exercise a large fraction of the code, so they may
 * catch errors in other parts than just in command-line parsing.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/cmdlineparser.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;

class CommandLineParserTest : public ::testing::Test
{
public:
    CommandLineParserTest();

    gmx::Options           options_;
    gmx::CommandLineParser parser_;
    bool                   flag_;
    std::vector<int>       ivalues_;
    std::vector<double>    dvalues_;
    int                    ivalue1p_;
    int                    ivalue12_;
    std::string            stringValue_;
};

CommandLineParserTest::CommandLineParserTest() :
    parser_(&options_), flag_(false), ivalue1p_(0), ivalue12_(0)
{
    using gmx::BooleanOption;
    using gmx::DoubleOption;
    using gmx::IntegerOption;
    using gmx::StringOption;
    options_.addOption(BooleanOption("flag").store(&flag_));
    options_.addOption(IntegerOption("mvi").storeVector(&ivalues_).multiValue());
    options_.addOption(DoubleOption("mvd").storeVector(&dvalues_).allowMultiple());
    options_.addOption(IntegerOption("1p").store(&ivalue1p_));
    options_.addOption(IntegerOption("12").store(&ivalue12_));
    options_.addOption(StringOption("str").store(&stringValue_));
}

TEST_F(CommandLineParserTest, HandlesSingleValues)
{
    const char* const cmdline[] = { "test", "-flag", "yes", "-mvi", "2", "-mvd", "2.7" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(7, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    EXPECT_TRUE(flag_);
    ASSERT_EQ(1U, ivalues_.size());
    EXPECT_EQ(2, ivalues_[0]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(2.7, dvalues_[0]);
}

TEST_F(CommandLineParserTest, HandlesBooleanWithoutArgument)
{
    const char* const cmdline[] = { "test", "-flag" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(2, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    EXPECT_TRUE(flag_);
}

TEST_F(CommandLineParserTest, HandlesBooleanAsNoWithoutArgument)
{
    const char* const cmdline[] = { "test", "-noflag" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(2, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    EXPECT_FALSE(flag_);
}

TEST_F(CommandLineParserTest, ThrowsWithBooleanAsNoWithArgument)
{
    const char* const cmdline[] = { "test", "-noflag", "no" };
    CommandLine       args(cmdline);
    EXPECT_THROW_GMX(parser_.parse(&args.argc(), args.argv()), gmx::InvalidInputError);
}

TEST_F(CommandLineParserTest, HandlesNegativeNumbers)
{
    const char* const cmdline[] = { "test", "-mvi", "1", "-2", "-mvd", "-2.7" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(6, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    ASSERT_EQ(2U, ivalues_.size());
    EXPECT_EQ(1, ivalues_[0]);
    EXPECT_EQ(-2, ivalues_[1]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(-2.7, dvalues_[0]);
}

TEST_F(CommandLineParserTest, HandlesString)
{
    const char* const cmdline[] = { "test", "-str", "text" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(3, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    EXPECT_EQ("text", stringValue_);
}

TEST_F(CommandLineParserTest, RejectsStringWithMultipleValues)
{
    const char* const cmdline[] = { "test", "-str", "text", "excess text" };
    CommandLine       args(cmdline);
    EXPECT_THROW_GMX(parser_.parse(&args.argc(), args.argv()), gmx::InvalidInputError);
}

TEST_F(CommandLineParserTest, HandlesDoubleDashOptionPrefix)
{
    const char* const cmdline[] = { "test", "--mvi", "1", "-2", "--mvd", "-2.7" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(6, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    ASSERT_EQ(2U, ivalues_.size());
    EXPECT_EQ(1, ivalues_[0]);
    EXPECT_EQ(-2, ivalues_[1]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(-2.7, dvalues_[0]);
}

TEST_F(CommandLineParserTest, HandlesOptionsStartingWithNumbers)
{
    const char* const cmdline[] = { "test", "--12", "1", "-1p", "-12" };
    CommandLine       args(cmdline);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    EXPECT_EQ(5, args.argc());
    EXPECT_STREQ("test", args.arg(0));

    EXPECT_EQ(1, ivalue12_);
    EXPECT_EQ(-12, ivalue1p_);
}

TEST_F(CommandLineParserTest, HandlesSkipUnknown)
{
    const char* const cmdline[] = { "test", "-unknown1", "-flag", "-unknown2", "value",
                                    "-mvi", "2",         "-mvd",  "2.7",       "-unknown3" };
    CommandLine       args(cmdline);
    parser_.skipUnknown(true);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());

    ASSERT_EQ(5, args.argc());
    EXPECT_STREQ("test", args.arg(0));
    EXPECT_STREQ("-unknown1", args.arg(1));
    EXPECT_STREQ("-unknown2", args.arg(2));
    EXPECT_STREQ("value", args.arg(3));
    EXPECT_STREQ("-unknown3", args.arg(4));
    EXPECT_TRUE(args.arg(5) == nullptr);

    EXPECT_TRUE(flag_);
    ASSERT_EQ(1U, ivalues_.size());
    EXPECT_EQ(2, ivalues_[0]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(2.7, dvalues_[0]);
}

TEST_F(CommandLineParserTest, RejectsPositionalArgumentsByDefault)
{
    // Ensures that "gmx trjconv f" gets rejected.
    const char* const cmdline[] = { "test", "module", "positional" };
    CommandLine       args(cmdline);
    EXPECT_THROW_GMX(parser_.parse(&args.argc(), args.argv()), gmx::InvalidInputError);
}

TEST_F(CommandLineParserTest, CanAllowPositionalArguments)
{
    // Ensures that "gmx help trjconv" works
    const char* const cmdline[] = { "test", "module", "positional", "-flag" };
    CommandLine       args(cmdline);
    parser_.allowPositionalArguments(true);
    ASSERT_NO_THROW_GMX(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW_GMX(options_.finish());
    ASSERT_EQ(4, args.argc());
    EXPECT_STREQ("test", args.arg(0));
    EXPECT_STREQ("module", args.arg(1));
    EXPECT_STREQ("positional", args.arg(2));
}

TEST_F(CommandLineParserTest, CannotHavePositionalArgumentsAfterOptions)
{
    // Even for the options that can't have arbitrary numbers of
    // values, there's no way to check whether there's been enough
    // values provided, so we can't have positional arguments after
    // any options.
    const char* const cmdline[] = { "test", "module", "-1p", "2", "positional" };
    CommandLine       args(cmdline);
    parser_.allowPositionalArguments(true);
    EXPECT_THROW_GMX(parser_.parse(&args.argc(), args.argv()), gmx::InvalidInputError);
}

} // namespace
} // namespace test
} // namespace gmx
