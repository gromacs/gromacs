/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Tests gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/cmdlinemodulemanager.h"

#include <gmock/gmock.h>

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

#include "cmdlinemodulemanagertest.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::MockModule;

//! Test fixture for the tests.
typedef gmx::test::CommandLineModuleManagerTestBase CommandLineModuleManagerTest;

TEST_F(CommandLineModuleManagerTest, RunsModule)
{
    const char *const cmdline[] = {
        "test", "module", "-flag", "yes"
    };
    CommandLine       args(cmdline);
    initManager(args, "test");
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, init(_));
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleHelp)
{
    const char *const cmdline[] = {
        "test", "help", "module"
    };
    CommandLine       args(cmdline);
    initManager(args, "test");
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    EXPECT_CALL(mod1, writeHelp(_));
    mod1.setExpectedDisplayName("test module");
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleHelpWithDashH)
{
    const char *const cmdline[] = {
        "test", "module", "-h"
    };
    CommandLine       args(cmdline);
    initManager(args, "test");
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    EXPECT_CALL(mod1, writeHelp(_));
    mod1.setExpectedDisplayName("test module");
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleHelpWithDashHWithSingleModule)
{
    const char *const cmdline[] = {
        "g_module", "-h"
    };
    CommandLine       args(cmdline);
    initManager(args, "g_module");
    MockModule        mod(NULL, NULL);
    manager().setSingleModule(&mod);
    using ::testing::_;
    EXPECT_CALL(mod, writeHelp(_));
    mod.setExpectedDisplayName("g_module");
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, HandlesConflictingBinaryAndModuleNames)
{
    const char *const cmdline[] = {
        "test", "test", "-flag", "yes"
    };
    CommandLine       args(cmdline);
    initManager(args, "test");
    MockModule       &mod1 = addModule("test", "Test module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, init(_));
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

} // namespace
