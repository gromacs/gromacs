/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Implements classes from cmdlinemodulemanagertest.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinemodulemanagertest.h"

#include <memory>
#include <string>
#include <utility>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/onlinehelp/ihelptopic.h"
#include "gromacs/onlinehelp/tests/mock_helptopic.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfileredirector.h"

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
void disableNice(gmx::CommandLineModuleSettings* settings)
{
    settings->setDefaultNiceLevel(0);
}

} // namespace

/********************************************************************
 * MockModule
 */

MockModule::MockModule(const char* name, const char* description) : name_(name), descr_(description)
{
    using ::testing::_;
    using ::testing::Invoke;
    ON_CALL(*this, init(_)).WillByDefault(Invoke(&disableNice));
    ON_CALL(*this, writeHelp(_)).WillByDefault(Invoke(this, &MockModule::checkHelpContext));
}

MockModule::~MockModule() {}

void MockModule::checkHelpContext(const gmx::CommandLineHelpContext& context) const
{
    EXPECT_EQ(expectedDisplayName_, context.moduleDisplayName());

    gmx::TextLineWrapperSettings settings;
    std::string                  moduleName =
            context.writerContext().substituteMarkupAndWrapToString(settings, "[THISMODULE]");
    EXPECT_EQ(expectedDisplayName_, moduleName);
}

/********************************************************************
 * MockOptionsModule
 */

MockOptionsModule::MockOptionsModule()
{
    using ::testing::_;
    using ::testing::Invoke;
    ON_CALL(*this, init(_)).WillByDefault(Invoke(&disableNice));
}

MockOptionsModule::~MockOptionsModule() {}

/********************************************************************
 * Test fixture for the tests
 */

class CommandLineModuleManagerTestBase::Impl
{
public:
    Impl(const CommandLine& args, const char* realBinaryName) :
        programContext_(args.argc(), args.argv()), manager_(realBinaryName, &programContext_)
    {
        manager_.setQuiet(true);
        manager_.setOutputRedirector(&redirector_);
    }

    TestFileOutputRedirector  redirector_;
    CommandLineProgramContext programContext_;
    CommandLineModuleManager  manager_;
};

CommandLineModuleManagerTestBase::CommandLineModuleManagerTestBase() : impl_(nullptr) {}

CommandLineModuleManagerTestBase::~CommandLineModuleManagerTestBase() {}

void CommandLineModuleManagerTestBase::initManager(const CommandLine& args, const char* realBinaryName)
{
    GMX_RELEASE_ASSERT(impl_ == nullptr, "initManager() called more than once in test");
    impl_ = std::make_unique<Impl>(args, realBinaryName);
}

MockModule& CommandLineModuleManagerTestBase::addModule(const char* name, const char* description)
{
    MockModule* module = new MockModule(name, description);
    manager().addModule(gmx::CommandLineModulePointer(module));
    return *module;
}

MockOptionsModule& CommandLineModuleManagerTestBase::addOptionsModule(const char* name, const char* description)
{
    std::unique_ptr<MockOptionsModule> modulePtr(new MockOptionsModule());
    MockOptionsModule*                 module = modulePtr.get();
    gmx::ICommandLineOptionsModule::registerModuleDirect(
            &manager(), name, description, std::move(modulePtr));
    return *module;
}

MockHelpTopic& CommandLineModuleManagerTestBase::addHelpTopic(const char* name, const char* title)
{
    MockHelpTopic* topic = new MockHelpTopic(name, title, "Help text");
    manager().addHelpTopic(gmx::HelpTopicPointer(topic));
    return *topic;
}

CommandLineModuleManager& CommandLineModuleManagerTestBase::manager()
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "initManager() not called");
    return impl_->manager_;
}

void CommandLineModuleManagerTestBase::checkRedirectedOutput()
{
    GMX_RELEASE_ASSERT(impl_ != nullptr, "initManager() not called");
    impl_->redirector_.checkRedirectedFiles(&checker());
}

} // namespace test
} // namespace gmx
