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
 * Implements classes from cmdlinemodulemanagertest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinemodulemanagertest.h"

#include <string>

#include <boost/scoped_ptr.hpp>
#include <gmock/gmock.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/utility/stringutil.h"

#include "gromacs/onlinehelp/tests/mock_helptopic.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \brief
 * Helper method to disable nice() calls from CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
void disableNice(gmx::CommandLineModuleSettings *settings)
{
    settings->setDefaultNiceLevel(0);
}

}       // namespace

/********************************************************************
 * MockModule
 */

MockModule::MockModule(const char *name, const char *description)
    : name_(name), descr_(description)
{
    using ::testing::_;
    using ::testing::Invoke;
    ON_CALL(*this, init(_))
        .WillByDefault(Invoke(&disableNice));
    ON_CALL(*this, writeHelp(_))
        .WillByDefault(Invoke(this, &MockModule::checkHelpContext));
}

MockModule::~MockModule()
{
}

void MockModule::checkHelpContext(const gmx::CommandLineHelpContext &context) const
{
    EXPECT_EQ(expectedDisplayName_, context.moduleDisplayName());

    gmx::TextLineWrapperSettings settings;
    std::string                  moduleName =
        context.writerContext().substituteMarkupAndWrapToString(
                settings, "[THISMODULE]");
    EXPECT_EQ(expectedDisplayName_, moduleName);
}

/********************************************************************
 * MockOptionsModule
 */

MockOptionsModule::MockOptionsModule()
{
    using ::testing::_;
    using ::testing::Invoke;
    ON_CALL(*this, init(_))
        .WillByDefault(Invoke(&disableNice));
}

MockOptionsModule::~MockOptionsModule()
{
}

/********************************************************************
 * Test fixture for the tests
 */

class CommandLineModuleManagerTestBase::Impl
{
    public:
        boost::scoped_ptr<CommandLineProgramContext> programContext_;
        boost::scoped_ptr<CommandLineModuleManager>  manager_;
        TestFileManager                              fileManager_;
};

CommandLineModuleManagerTestBase::CommandLineModuleManagerTestBase()
    : impl_(new Impl)
{
}

CommandLineModuleManagerTestBase::~CommandLineModuleManagerTestBase()
{
}

void CommandLineModuleManagerTestBase::initManager(
        const CommandLine &args, const char *realBinaryName)
{
    impl_->manager_.reset();
    impl_->programContext_.reset(
            new gmx::CommandLineProgramContext(args.argc(), args.argv()));
    impl_->manager_.reset(new gmx::CommandLineModuleManager(
                                  realBinaryName, impl_->programContext_.get()));
    impl_->manager_->setQuiet(true);
}

MockModule &
CommandLineModuleManagerTestBase::addModule(const char *name, const char *description)
{
    MockModule *module = new MockModule(name, description);
    manager().addModule(gmx::CommandLineModulePointer(module));
    return *module;
}

MockOptionsModule &
CommandLineModuleManagerTestBase::addOptionsModule(const char *name, const char *description)
{
    MockOptionsModule *module = new MockOptionsModule();
    gmx::CommandLineOptionsModuleInterface::registerModule(
            &manager(), name, description, module);
    return *module;
}

MockHelpTopic &
CommandLineModuleManagerTestBase::addHelpTopic(const char *name, const char *title)
{
    MockHelpTopic *topic = new MockHelpTopic(name, title, "Help text");
    manager().addHelpTopic(gmx::HelpTopicPointer(topic));
    return *topic;
}

CommandLineModuleManager &CommandLineModuleManagerTestBase::manager()
{
    return *impl_->manager_;
}

void CommandLineModuleManagerTestBase::redirectManagerOutput()
{
    impl_->manager_->setOutputRedirector(&initOutputRedirector(&impl_->fileManager_));
}

} // namespace test
} // namespace gmx
