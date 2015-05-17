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
 * Test fixture and helper classes for tests using gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGERTEST_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGERTEST_H

#include <string>

#include <gmock/gmock.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/utility/classhelpers.h"

#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{

class CommandLine;
class MockHelpTopic;

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
        ~MockModule();

        virtual const char *name() const { return name_; }
        virtual const char *shortDescription() const { return descr_; }

        MOCK_METHOD1(init, void(gmx::CommandLineModuleSettings *settings));
        MOCK_METHOD2(run, int(int argc, char *argv[]));
        MOCK_CONST_METHOD1(writeHelp, void(const gmx::CommandLineHelpContext &context));

        //! Sets the expected display name for writeHelp() calls.
        void setExpectedDisplayName(const char *expected)
        {
            expectedDisplayName_ = expected;
        }

    private:
        //! Checks the context passed to writeHelp().
        void checkHelpContext(const gmx::CommandLineHelpContext &context) const;

        const char             *name_;
        const char             *descr_;
        std::string             expectedDisplayName_;
};

/*! \internal \brief
 * Mock implementation of gmx::CommandLineOptionsModuleInterface.
 *
 * \ingroup module_commandline
 */
class MockOptionsModule : public gmx::CommandLineOptionsModuleInterface
{
    public:
        MockOptionsModule();
        ~MockOptionsModule();

        MOCK_METHOD1(init, void(gmx::CommandLineModuleSettings *settings));
        MOCK_METHOD1(initOptions, void(gmx::Options *options));
        MOCK_METHOD1(optionsFinished, void(gmx::Options *options));
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
        ~CommandLineModuleManagerTestBase();

        //! Creates the manager to run the given command line.
        void initManager(const CommandLine &args, const char *realBinaryName);
        //! Adds a mock module to the manager.
        MockModule               &addModule(const char *name, const char *description);
        //! Adds a mock module using gmx::Options to the manager.
        MockOptionsModule        &addOptionsModule(const char *name, const char *description);
        //! Adds a mock help topic to the manager.
        MockHelpTopic            &addHelpTopic(const char *name, const char *title);

        /*! \brief
         * Returns the manager for this test.
         *
         * initManager() must have been called.
         */
        CommandLineModuleManager &manager();

        /*! \brief
         * Redirects all manager output to files.
         *
         * Can be used to silence tests that would otherwise print out
         * something, and/or checkRedirectedFileContents() from the base class
         * can be used to check the output.
         *
         * The manager is put into quiet mode by default, so the manager will
         * only print out information if, e.g., help is explicitly requested.
         */
        void redirectManagerOutput();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
