/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests gmx::CommandLineHelpModule through gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include <cstdio>

#include <filesystem>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/onlinehelp/tests/mock_helptopic.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

#include "cmdlinemodulemanagertest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::MockHelpTopic;
using gmx::test::MockOptionsModule;

//! Test fixture for the tests.
typedef gmx::test::CommandLineModuleManagerTestBase CommandLineHelpModuleTest;

TEST_F(CommandLineHelpModuleTest, PrintsGeneralHelp)
{
    const char* const cmdline[] = { "test" };
    CommandLine       args(cmdline);
    initManager(args, "test");
    addModule("module", "First module");
    addModule("other", "Second module");
    addHelpTopic("topic", "Test topic");
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
    checkRedirectedOutput();
}

TEST_F(CommandLineHelpModuleTest, PrintsHelpOnTopic)
{
    const char* const cmdline[] = { "test", "help", "topic" };
    CommandLine       args(cmdline);
    initManager(args, "test");
    addModule("module", "First module");
    MockHelpTopic& topic = addHelpTopic("topic", "Test topic");
    topic.addSubTopic("sub1", "Subtopic 1", "");
    topic.addSubTopic("sub2", "Subtopic 2", "");
    using ::testing::_;
    EXPECT_CALL(topic, writeHelp(_));
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
    checkRedirectedOutput();
}

/*! \brief
 * Initializes Options for help export tests.
 *
 * \ingroup module_commandline
 */
void initOptionsBasic(gmx::IOptionsContainer* options, gmx::ICommandLineOptionsModuleSettings* settings)
{
    const char* const desc[] = { "Sample description", "for testing [THISMODULE]." };
    settings->setHelpText(desc);
    const char* const bug[] = { "Known issue for [THISMODULE].",
                                "With another bug for [THISMODULE]." };
    settings->setBugText(bug);
    options->addOption(gmx::IntegerOption("int").description("Integer option"));
}

TEST_F(CommandLineHelpModuleTest, ExportsHelp)
{
    const char* const cmdline[] = { "test", "help", "-export", "rst" };
    // TODO: Find a more elegant solution, or get rid of the links.dat altogether.
    gmx::TextWriter::writeFileFromString("links.dat", "");
    CommandLine args(cmdline);
    initManager(args, "test");
    MockOptionsModule& mod1 = addOptionsModule("module", "First module");
    MockOptionsModule& mod2 = addOptionsModule("other", "Second module");
    {
        gmx::CommandLineModuleGroup group = manager().addModuleGroup("Group 1");
        group.addModule("module");
    }
    {
        gmx::CommandLineModuleGroup group = manager().addModuleGroup("Group 2");
        group.addModule("other");
    }
    MockHelpTopic& topic1 = addHelpTopic("topic1", "Test topic");
    MockHelpTopic& sub1   = topic1.addSubTopic("sub1", "Subtopic 1", "Sub text");
    MockHelpTopic& sub2   = topic1.addSubTopic("sub2", "Subtopic 2", "Sub text");
    MockHelpTopic& sub3   = topic1.addSubTopic("other", "Out-of-order subtopic", "Sub text");
    MockHelpTopic& topic2 = addHelpTopic("topic2", "Another topic");
    using ::testing::_;
    using ::testing::Invoke;
    EXPECT_CALL(mod1, initOptions(_, _)).WillOnce(Invoke(&initOptionsBasic));
    EXPECT_CALL(mod2, initOptions(_, _));
    EXPECT_CALL(topic1, writeHelp(_));
    EXPECT_CALL(sub1, writeHelp(_));
    EXPECT_CALL(sub2, writeHelp(_));
    EXPECT_CALL(sub3, writeHelp(_));
    EXPECT_CALL(topic2, writeHelp(_));
    int rc = 0;
    ASSERT_NO_THROW_GMX(rc = manager().run(args.argc(), args.argv()));
    ASSERT_EQ(0, rc);
    checkRedirectedOutput();
    std::remove("links.dat");
}

} // namespace
} // namespace test
} // namespace gmx
