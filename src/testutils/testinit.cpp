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
 * Implements functions in testinit.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testinit.h"

#include <cstdio>
#include <cstdlib>

#include <filesystem>
#include <memory>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/math/utilities.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/test_hardware_environment.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

#include "buildinfo.h"
#include "mpi_printer.h"

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
 * installed.  It also supports overriding the directory through a command-line
 * option to the test binary.
 *
 * \ingroup module_testutils
 */
class TestProgramContext : public IProgramContext
{
public:
    /*! \brief
     * Initializes a test program context.
     *
     * \param[in] context  Current \Gromacs program context.
     */
    explicit TestProgramContext(const IProgramContext& context) :
        context_(context), dataPath_(CMAKE_SOURCE_DIR)
    {
    }

    /*! \brief
     * Sets the source directory root from which to look for data files.
     */
    void overrideSourceRoot(const std::filesystem::path& sourceRoot) { dataPath_ = sourceRoot; }

    const char*            programName() const override { return context_.programName(); }
    const char*            displayName() const override { return context_.displayName(); }
    std::filesystem::path  fullBinaryPath() const override { return context_.fullBinaryPath(); }
    InstallationPrefixInfo installationPrefix() const override
    {
        return InstallationPrefixInfo(dataPath_, true);
    }
    const char* commandLine() const override { return context_.commandLine(); }

private:
    const IProgramContext& context_;
    std::filesystem::path  dataPath_;
};

//! Prints the command-line options for the unit test binary.
void printHelp(const Options& options)
{
    const std::string& program = getProgramContext().displayName();
    std::fprintf(stderr,
                 "\nYou can use the following GROMACS-specific command-line flags\n"
                 "to control the behavior of the tests:\n\n");
    TextWriter             writer(&TextOutputFile::standardError());
    CommandLineHelpContext context(&writer, eHelpOutputFormat_Console, nullptr, program);
    context.setModuleDisplayName(program);
    CommandLineHelpWriter(options).writeHelp(context);
}

//! Global program context instance for test binaries.
// Never releases ownership.
std::unique_ptr<TestProgramContext> g_testContext;

/*! \brief Makes GoogleTest non-failures more verbose
 *
 * By default, GoogleTest does not echo messages appended to explicit
 * assertions of success with SUCCEEDED() e.g.
 *
 *    GTEST_SKIP() << "reason why";
 *
 * produces no output. This test listener changes that behavior, so
 * that the message is echoed.
 *
 * When run with multiple ranks, only the main rank should use this
 * listener, else the output can be very noisy. */
class SuccessListener : public testing::EmptyTestEventListener
{
    void OnTestPartResult(const testing::TestPartResult& result) override
    {
        if (result.type() == testing::TestPartResult::kSuccess)
        {
            printf("%s\n", result.message());
        }
    }
};

} // namespace

//! \cond internal
void initTestUtils(const std::filesystem::path& dataPath,
                   const std::filesystem::path& tempPath,
                   bool                         usesMpi,
                   bool                         usesHardwareDetection,
                   bool                         registersDynamically,
                   int*                         argc,
                   char***                      argv)
{
    if (gmxShouldEnableFPExceptions())
    {
        gmx_feenableexcept();
    }
    const CommandLineProgramContext& context = initForCommandLine(argc, argv);
    try
    {
        if (!usesMpi && gmx_node_num() > 1)
        {
            // We cannot continue, since some tests might be using
            // MPI_COMM_WORLD, which could deadlock if we would only
            // continue with the main rank here.
            if (gmx_node_rank() == 0)
            {
                fprintf(stderr,
                        "NOTE: You are running %s on %d MPI ranks, "
                        "but it does not contain MPI-enabled tests. "
                        "The test will now exit.\n",
                        context.programName(),
                        gmx_node_num());
            }
            finalizeForCommandLine();
            std::exit(1);
        }
        if (registersDynamically || usesHardwareDetection)
        {
            ::gmx::test::TestHardwareEnvironment::gmxSetUp();
        }
        // Note that setting up the hardware environment should
        // precede dynamic registration, in case it is needed for
        // deciding what to register.
        if (registersDynamically)
        {
            ::gmx::test::registerTestsDynamically();
        }
        g_testContext = std::make_unique<TestProgramContext>(context);
        setProgramContext(g_testContext.get());
        // Use the default finder that does not respect GMXLIB, since the tests
        // generally can only get confused by a different set of data files.
        setLibraryFileFinder(nullptr);
        ::testing::InitGoogleMock(argc, *argv);
        if (!dataPath.empty())
        {
            TestFileManager::setInputDataDirectory(
                    std::filesystem::path(CMAKE_SOURCE_DIR).append(dataPath.string()));
        }
        if (!tempPath.empty())
        {
            TestFileManager::setGlobalOutputTempDirectory(tempPath);
        }
        TestFileManager::setTestSimulationDatabaseDirectory(GMX_TESTSIMULATIONDATABASE_DIR);

        bool        bHelp = false;
        std::string sourceRoot;
        bool        echoReasons = false;
        Options     options;
        // TODO: A single option that accepts multiple names would be nicer.
        // Also, we recognize -help, but GTest doesn't, which leads to a bit
        // unintuitive behavior.
        options.addOption(BooleanOption("h").store(&bHelp).description(
                "Print GROMACS-specific unit test options"));
        options.addOption(BooleanOption("help").store(&bHelp).hidden());
        options.addOption(BooleanOption("?").store(&bHelp).hidden());
        // TODO: Make this into a FileNameOption (or a DirectoryNameOption).
        options.addOption(
                StringOption("src-root").store(&sourceRoot).description("Override source tree location (for data files)"));
        options.addOption(
                BooleanOption("echo-reasons").store(&echoReasons).description("When succeeding or skipping a test, echo the reason"));
        // The potential MPI test event listener must be initialized first,
        // because it should appear in the start of the event listener list,
        // before other event listeners that may generate test failures
        // (currently, such an event listener is used by the reference data
        // framework).
        if (usesMpi)
        {
            initMPIOutput();
        }
        // TODO: Consider removing this option from test binaries that do not need it.
        initReferenceData(&options);
        initTestOptions(&options);
        try
        {
            CommandLineParser(&options).parse(argc, *argv);
            options.finish();
        }
        catch (const UserInputError&)
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
                    std::filesystem::path(sourceRoot).append(dataPath.string()));
        }
        // Echo success messages only from the main MPI rank
        if (echoReasons && (gmx_node_rank() == 0))
        {
            testing::UnitTest::GetInstance()->listeners().Append(new SuccessListener);
        }
    }
    catch (const std::exception& ex)
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

void finalizeTestUtils(const bool usesHardwareDetection, const bool registersDynamically)
{
    if (registersDynamically || usesHardwareDetection)
    {
        ::gmx::test::TestHardwareEnvironment::gmxTearDown();
    }
    setProgramContext(nullptr);
    g_testContext.reset();
    finalizeForCommandLine();
}
//! \endcond

} // namespace test
} // namespace gmx
