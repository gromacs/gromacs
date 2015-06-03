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
 * Implements functions in testinit.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testinit.h"

#include <cstdio>
#include <cstdlib>

#include <boost/scoped_ptr.hpp>
#include <gmock/gmock.h>

#include "buildinfo.h"
#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/math/utilities.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"

#include "testutils/mpi-printer.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \brief
 * Custom program context for test binaries.
 *
 * This context overrides the installationPrefix() implementation to always
 * load data files from the source directory, as the test binaries are never
 * installed.  It also support overriding the directory through a command-line
 * option to the test binary.
 *
 * \ingroup module_testutils
 */
class TestProgramContext : public ProgramContextInterface
{
    public:
        /*! \brief
         * Initializes a test program context.
         *
         * \param[in] context  Current \Gromacs program context.
         */
        explicit TestProgramContext(const ProgramContextInterface &context)
            : context_(context), dataPath_(CMAKE_SOURCE_DIR)
        {
        }

        /*! \brief
         * Sets the source directory root from which to look for data files.
         */
        void overrideSourceRoot(const std::string &sourceRoot)
        {
            dataPath_ = sourceRoot;
        }

        virtual const char *programName() const
        {
            return context_.programName();
        }
        virtual const char *displayName() const
        {
            return context_.displayName();
        }
        virtual const char *fullBinaryPath() const
        {
            return context_.fullBinaryPath();
        }
        virtual InstallationPrefixInfo installationPrefix() const
        {
            return InstallationPrefixInfo(dataPath_.c_str(), true);
        }
        virtual const char *commandLine() const
        {
            return context_.commandLine();
        }

    private:
        const ProgramContextInterface   &context_;
        std::string                      dataPath_;
};

//! Prints the command-line options for the unit test binary.
void printHelp(const Options &options)
{
    const std::string &program = getProgramContext().displayName();
    std::fprintf(stderr,
                 "\nYou can use the following GROMACS-specific command-line flags\n"
                 "to control the behavior of the tests:\n\n");
    CommandLineHelpContext context(&File::standardError(),
                                   eHelpOutputFormat_Console, NULL, program);
    context.setModuleDisplayName(program);
    CommandLineHelpWriter(options).writeHelp(context);
}

//! Global program context instance for test binaries.
boost::scoped_ptr<TestProgramContext> g_testContext;

}       // namespace

//! \cond internal
void initTestUtils(const char *dataPath, const char *tempPath, int *argc, char ***argv)
{
#ifndef NDEBUG
    gmx_feenableexcept();
#endif
    const CommandLineProgramContext &context = initForCommandLine(argc, argv);
    try
    {
        g_testContext.reset(new TestProgramContext(context));
        setProgramContext(g_testContext.get());
        // Use the default finder that does not respect GMXLIB, since the tests
        // generally can only get confused by a different set of data files.
        setLibraryFileFinder(NULL);
        ::testing::InitGoogleMock(argc, *argv);
        if (dataPath != NULL)
        {
            TestFileManager::setInputDataDirectory(
                    Path::join(CMAKE_SOURCE_DIR, dataPath));
        }
        if (tempPath != NULL)
        {
            TestFileManager::setGlobalOutputTempDirectory(tempPath);
        }
        bool        bHelp = false;
        std::string sourceRoot;
        Options     options(NULL, NULL);
        // TODO: A single option that accepts multiple names would be nicer.
        // Also, we recognize -help, but GTest doesn't, which leads to a bit
        // unintuitive behavior.
        options.addOption(BooleanOption("h").store(&bHelp)
                              .description("Print GROMACS-specific unit test options"));
        options.addOption(BooleanOption("help").store(&bHelp).hidden());
        options.addOption(BooleanOption("?").store(&bHelp).hidden());
        // TODO: Make this into a FileNameOption (or a DirectoryNameOption).
        options.addOption(StringOption("src-root").store(&sourceRoot)
                              .description("Override source tree location (for data files)"));
        // TODO: Consider removing this option from test binaries that do not need it.
        initReferenceData(&options);
        initTestOptions(&options);
        try
        {
            CommandLineParser(&options).parse(argc, *argv);
            options.finish();
        }
        catch (const UserInputError &)
        {
            printHelp(options);
            throw;
        }
        if (bHelp)
        {
            printHelp(options);
        }
        if (!sourceRoot.empty())
        {
            g_testContext->overrideSourceRoot(sourceRoot);
            TestFileManager::setInputDataDirectory(
                    Path::join(sourceRoot, dataPath));
        }
        setFatalErrorHandler(NULL);
        initMPIOutput();
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        int retcode = processExceptionAtExitForCommandLine(ex);
        // TODO: It could be nice to destroy things in proper order such that
        // g_testContext would not contain hanging references at this point,
        // but in practice that should not matter.
        g_testContext.reset();
        std::exit(retcode);
    }
}

void finalizeTestUtils()
{
    setProgramContext(NULL);
    g_testContext.reset();
    finalizeForCommandLine();
}
//! \endcond

} // namespace test
} // namespace gmx
