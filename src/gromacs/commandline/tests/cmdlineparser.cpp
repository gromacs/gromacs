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
#include <cstdlib>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"

namespace
{

class CommandLineParserTest : public ::testing::Test
{
    public:
        CommandLineParserTest();
        ~CommandLineParserTest();

        void createArguments(const char *cmdline[]);

        gmx::Options            _options;
        gmx::CommandLineParser  _parser;
        bool                    _flag;
        std::vector<int>        _ivalues;
        std::vector<double>     _dvalues;
        int                     _argc;
        char                  **_argv;
};

CommandLineParserTest::CommandLineParserTest()
    : _options(NULL, NULL), _parser(&_options),
      _flag(false),
      _argc(0), _argv(NULL)
{
    using gmx::BooleanOption;
    using gmx::IntegerOption;
    using gmx::DoubleOption;
    _options.addOption(BooleanOption("flag").store(&_flag));
    _options.addOption(IntegerOption("mvi").storeVector(&_ivalues).multiValue());
    _options.addOption(DoubleOption("mvd").storeVector(&_dvalues).allowMultiple());
}

CommandLineParserTest::~CommandLineParserTest()
{
    if (_argv != NULL)
    {
        for (int i = 0; i < _argc; ++i)
        {
            free(_argv[i]);
        }
    }
    delete [] _argv;
}

void CommandLineParserTest::createArguments(const char *cmdline[])
{
    _argc = 0;
    while (cmdline[_argc] != NULL) ++_argc;
    ++_argc;

    _argv = new char *[_argc];
    _argv[0] = strdup("test");
    for (int i = 1; i < _argc; ++i)
    {
        _argv[i] = strdup(cmdline[i - 1]);
    }
}

TEST_F(CommandLineParserTest, HandlesSingleValues)
{
    const char *cmdline[] = {"-flag", "yes", "-mvi", "2", "-mvd", "2.7", NULL};
    createArguments(cmdline);
    ASSERT_NO_THROW(_parser.parse(&_argc, _argv));
    ASSERT_NO_THROW(_options.finish());

    EXPECT_TRUE(_flag);
    ASSERT_EQ(1U, _ivalues.size());
    EXPECT_EQ(2, _ivalues[0]);
    ASSERT_EQ(1U, _dvalues.size());
    EXPECT_DOUBLE_EQ(2.7, _dvalues[0]);
}

TEST_F(CommandLineParserTest, HandlesNegativeNumbers)
{
    const char *cmdline[] = {"-mvi", "1", "-2", "-mvd", "-2.7", NULL};
    createArguments(cmdline);
    ASSERT_NO_THROW(_parser.parse(&_argc, _argv));
    ASSERT_NO_THROW(_options.finish());

    ASSERT_EQ(2U, _ivalues.size());
    EXPECT_EQ(1, _ivalues[0]);
    EXPECT_EQ(-2, _ivalues[1]);
    ASSERT_EQ(1U, _dvalues.size());
    EXPECT_DOUBLE_EQ(-2.7, _dvalues[0]);
}

} // namespace
