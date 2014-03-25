/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Implements functions in testoptions.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "testoptions.h"

#include <cstdio>
#include <cstdlib>

#include <list>

#include <gmock/gmock.h>

#include "thread_mpi/mutex.h"

#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/common.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/programcontext.h"

#include "refdata.h"
#include "testfilemanager.h"
#include "mpi-printer.h"

namespace gmx
{
namespace test
{

namespace
{

/*! \brief
 * Singleton registry for test options added with #GMX_TEST_OPTIONS.
 *
 * \ingroup module_testutils
 */
class TestOptionsRegistry
{
    public:
        //! Returns the singleton instance of this class.
        static TestOptionsRegistry &getInstance()
        {
            static TestOptionsRegistry singleton;
            return singleton;
        }

        //! Adds a provider into the registry.
        void add(const char * /*name*/, TestOptionsProvider *provider)
        {
            tMPI::lock_guard<tMPI::mutex> lock(listMutex_);
            providerList_.push_back(provider);
        }

        //! Initializes the options from all the provides.
        void initOptions(Options *options);

    private:
        TestOptionsRegistry() {}

        typedef std::list<TestOptionsProvider *> ProviderList;

        tMPI::mutex             listMutex_;
        ProviderList            providerList_;

        GMX_DISALLOW_COPY_AND_ASSIGN(TestOptionsRegistry);
};

void TestOptionsRegistry::initOptions(Options *options)
{
    // TODO: Have some deterministic order for the options; now it depends on
    // the order in which the global initializers are run.
    tMPI::lock_guard<tMPI::mutex> lock(listMutex_);
    ProviderList::const_iterator  i;
    for (i = providerList_.begin(); i != providerList_.end(); ++i)
    {
        (*i)->initOptions(options);
    }
}

//! Prints the command-line options for the unit test binary.
void printHelp(const Options &options)
{
    std::fprintf(stderr,
                 "\nYou can use the following GROMACS-specific command-line flags\n"
                 "to control the behavior of the tests:\n\n");
    CommandLineHelpContext context(&File::standardError(),
                                   eHelpOutputFormat_Console, NULL);
    context.setModuleDisplayName(getProgramContext().displayName());
    CommandLineHelpWriter(options).writeHelp(context);
}

}       // namespace

void registerTestOptions(const char *name, TestOptionsProvider *provider)
{
    TestOptionsRegistry::getInstance().add(name, provider);
}

void initTestUtils(const char *dataPath, const char *tempPath, int *argc, char ***argv)
{
    gmx::initForCommandLine(argc, argv);
    try
    {
        ::testing::InitGoogleMock(argc, *argv);
        if (dataPath != NULL)
        {
            TestFileManager::setInputDataDirectory(dataPath);
        }
        if (tempPath != NULL)
        {
            TestFileManager::setGlobalOutputTempDirectory(tempPath);
        }
        bool    bHelp = false;
        Options options(NULL, NULL);
        // TODO: A single option that accepts multiple names would be nicer.
        // Also, we recognize -help, but GTest doesn't, which leads to a bit
        // unintuitive behavior.
        options.addOption(BooleanOption("h").store(&bHelp)
                              .description("Print GROMACS-specific unit test options"));
        options.addOption(BooleanOption("help").store(&bHelp).hidden());
        options.addOption(BooleanOption("?").store(&bHelp).hidden());
        // TODO: Consider removing this option from test binaries that do not need it.
        initReferenceData(&options);
        TestOptionsRegistry::getInstance().initOptions(&options);
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
        setFatalErrorHandler(NULL);
        initMPIOutput();
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        std::exit(processExceptionAtExit(ex));
    }
}

void finalizeTestUtils()
{
    gmx::finalizeForCommandLine();
}

} // namespace test
} // namespace gmx
