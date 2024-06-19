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
 * Test fixture and helper classes for tests using gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGERTEST_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGERTEST_H

#include <memory>
#include <string>

#include <gmock/gmock.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"

#include "testutils/stringtest.h"

namespace gmx
{
class CommandLineModuleManager;
class IOptionsContainer;

namespace test
{

class CommandLine;
class MockHelpTopic;
class TestFileOutputRedirector;

/*! \internal \brief
 * Mock implementation of gmx::ICommandLineModule.
 *
 * \ingroup module_commandline
 */
class MockModule : public gmx::ICommandLineModule
{
public:
    //! Creates a mock module with the given name and description.
    MockModule(const char* name, const char* description);
    ~MockModule() override;

    const char* name() const override { return name_; }
    const char* shortDescription() const override { return descr_; }

    MOCK_METHOD1(init, void(gmx::CommandLineModuleSettings* settings));
    MOCK_METHOD2(run, int(int argc, char* argv[]));
    MOCK_CONST_METHOD1(writeHelp, void(const gmx::CommandLineHelpContext& context));

    //! Sets the expected display name for writeHelp() calls.
    void setExpectedDisplayName(const char* expected) { expectedDisplayName_ = expected; }

private:
    //! Checks the context passed to writeHelp().
    void checkHelpContext(const gmx::CommandLineHelpContext& context) const;

    const char* name_;
    const char* descr_;
    std::string expectedDisplayName_;
};

/*! \internal \brief
 * Mock implementation of gmx::ICommandLineOptionsModule.
 *
 * \ingroup module_commandline
 */
class MockOptionsModule : public gmx::ICommandLineOptionsModule
{
public:
    MockOptionsModule();
    ~MockOptionsModule() override;

    MOCK_METHOD1(init, void(gmx::CommandLineModuleSettings* settings));
    MOCK_METHOD2(initOptions,
                 void(gmx::IOptionsContainer* options, gmx::ICommandLineOptionsModuleSettings* settings));
    MOCK_METHOD0(optionsFinished, void());
    MOCK_METHOD0(run, int());
};

/*! \internal \brief
 * Test fixture for tests using gmx::CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
class CommandLineModuleManagerTestBase : public gmx::test::StringTestBase
{
public:
    CommandLineModuleManagerTestBase();
    ~CommandLineModuleManagerTestBase() override;

    //! Creates the manager to run the given command line.
    void initManager(const CommandLine& args, const char* realBinaryName);
    //! Adds a mock module to the manager.
    MockModule& addModule(const char* name, const char* description);
    //! Adds a mock module using gmx::Options to the manager.
    MockOptionsModule& addOptionsModule(const char* name, const char* description);
    //! Adds a mock help topic to the manager.
    MockHelpTopic& addHelpTopic(const char* name, const char* title);

    /*! \brief
     * Returns the manager for this test.
     *
     * initManager() must have been called.
     */
    CommandLineModuleManager& manager();

    /*! \brief
     * Checks all output from the manager using reference data.
     *
     * Both output to `stdout` and to files is checked.
     *
     * The manager is put into quiet mode by default, so the manager will
     * only print out information if, e.g., help is explicitly requested.
     */
    void checkRedirectedOutput();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
