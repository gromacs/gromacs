/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests gmx::CommandLineParser.
 *
 * These tests exercise a large fraction of the code, so they may
 * catch errors in other parts than just in command-line parsing.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"

#include "testutils/cmdlinetest.h"

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
    ASSERT_NO_THROW(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW(options_.finish());

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
    ASSERT_NO_THROW(parser_.parse(&args.argc(), args.argv()));
    ASSERT_NO_THROW(options_.finish());

    ASSERT_EQ(2U, ivalues_.size());
    EXPECT_EQ(1, ivalues_[0]);
    EXPECT_EQ(-2, ivalues_[1]);
    ASSERT_EQ(1U, dvalues_.size());
    EXPECT_DOUBLE_EQ(-2.7, dvalues_[0]);
}

} // namespace
