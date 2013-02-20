/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
 * Tests gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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
#include "gromacs/utility/programinfo.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mock_helptopic.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::MockHelpTopic;

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
        MOCK_CONST_METHOD1(writeHelp, void(const gmx::HelpWriterContext &context));

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
        void initManager(const CommandLine &args);
        MockModule    &addModule(const char *name, const char *description);
        MockHelpTopic &addHelpTopic(const char *name, const char *title);

        gmx::CommandLineModuleManager &manager() { return *manager_; }

    private:
        boost::scoped_ptr<gmx::ProgramInfo>              programInfo_;
        boost::scoped_ptr<gmx::CommandLineModuleManager> manager_;
};

void CommandLineModuleManagerTest::initManager(const CommandLine &args)
{
    manager_.reset();
    programInfo_.reset(new gmx::ProgramInfo("g_test", args.argc(), args.argv()));
    manager_.reset(new gmx::CommandLineModuleManager(*programInfo_));
}

MockModule &
CommandLineModuleManagerTest::addModule(const char *name, const char *description)
{
    MockModule *module = new MockModule(name, description);
    manager().addModule(gmx::CommandLineModulePointer(module));
    return *module;
}

MockHelpTopic &
CommandLineModuleManagerTest::addHelpTopic(const char *name, const char *title)
{
    MockHelpTopic *topic = new MockHelpTopic(name, title, "Help text");
    manager().addHelpTopic(gmx::HelpTopicPointer(topic));
    return *topic;
}

/********************************************************************
 * Actual tests
 */

TEST_F(CommandLineModuleManagerTest, RunsModule)
{
    const char *const cmdline[] = {
        "g_test", "module", "-flag", "yes"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleHelp)
{
    const char *const cmdline[] = {
        "g_test", "help", "module"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    EXPECT_CALL(mod1, writeHelp(_));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, PrintsHelpOnTopic)
{
    const char *const cmdline[] = {
        "g_test", "help", "topic"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    addModule("module", "First module");
    MockHelpTopic &topic = addHelpTopic("topic", "Test topic");
    using ::testing::_;
    EXPECT_CALL(topic, writeHelp(_));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleBasedOnBinaryName)
{
    const char *const cmdline[] = {
        "g_module", "-flag", "yes"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv(), args.argc())));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, RunsModuleBasedOnBinaryNameWithPathAndSuffix)
{
    const char *const cmdline[] = {
        "/usr/local/gromacs/bin/g_module" GMX_BINARY_SUFFIX ".exe", "-flag", "yes"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    MockModule       &mod1 = addModule("module", "First module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv(), args.argc())));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

TEST_F(CommandLineModuleManagerTest, HandlesConflictingBinaryAndModuleNames)
{
    const char *const cmdline[] = {
        "g_test", "test", "-flag", "yes"
    };
    CommandLine       args(CommandLine::create(cmdline));
    initManager(args);
    MockModule       &mod1 = addModule("test", "Test module");
    addModule("other", "Second module");
    using ::testing::_;
    using ::testing::Args;
    using ::testing::ElementsAreArray;
    EXPECT_CALL(mod1, run(_, _))
        .With(Args<1, 0>(ElementsAreArray(args.argv() + 1, args.argc() - 1)));
    int rc = 0;
    ASSERT_NO_THROW(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
}

} // namespace
