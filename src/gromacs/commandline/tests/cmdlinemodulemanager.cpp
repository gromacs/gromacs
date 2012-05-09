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
 * Tests gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
// For GMX_BINARY_SUFFIX
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>

#include <gmock/gmock.h>

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"

#include "cmdlinetest.h"

namespace
{

/********************************************************************
 * MockModule
 */

/*! \internal \brief
 * Mock implementation of gmx::CommandLineModuleInterface.
 *
 * \ingroup module_commandline
 */
class MockModule : public gmx::CommandLineModuleInterface
{
    public:
        //! Creates a mock module with the given name and description.
        MockModule(const char *name, const char *description);

        virtual const char *name() const { return name_; }
        virtual const char *shortDescription() const { return descr_; }

        MOCK_METHOD2(run, int(int argc, char *argv[]));

    private:
        const char             *name_;
        const char             *descr_;
};

MockModule::MockModule(const char *name, const char *description)
    : name_(name), descr_(description)
{
}

/********************************************************************
 * Test fixture for the tests
 */

class CommandLineModuleManagerTest : public ::testing::Test
{
    public:
        CommandLineModuleManagerTest();

        MockModule &addModule(const char *name, const char *description);

        gmx::CommandLineModuleManager manager_;
};

CommandLineModuleManagerTest::CommandLineModuleManagerTest()
    : manager_("g_test")
{
}

MockModule &
CommandLineModuleManagerTest::addModule(const char *name, const char *description)
{
    MockModule *module = new MockModule(name, description);
    manager_.addModule(gmx::CommandLineModulePointer(module));
    return *module;
}

/********************************************************************
 * Actual tests
 */

TEST_F(CommandLineModuleManagerTest, RunsModule)
{
    const char *const cmdline[] = {
        "g_test", "module", "-flag", "yes"
    };
    gmx::test::CommandLine args(cmdline);
    MockModule &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager_.run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleBasedOnBinaryName)
{
    const char *const cmdline[] = {
        "g_module", "-flag", "yes"
    };
    gmx::test::CommandLine args(cmdline);
    MockModule &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv(), args.argc())));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager_.run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleBasedOnBinaryNameWithPathAndSuffix)
{
    const char *const cmdline[] = {
        "/usr/local/gromacs/bin/g_module" GMX_BINARY_SUFFIX ".exe", "-flag", "yes"
    };
    gmx::test::CommandLine args(cmdline);
    MockModule &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv(), args.argc())));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager_.run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, HandlesConflictingBinaryAndModuleNames)
{
    const char *const cmdline[] = {
        "g_test", "test", "-flag", "yes"
    };
    gmx::test::CommandLine args(cmdline);
    MockModule &mod1 = addModule("test", "Test module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager_.run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

} // namespace
