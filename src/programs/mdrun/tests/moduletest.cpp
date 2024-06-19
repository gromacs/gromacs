/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Implements classes in moduletest.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "moduletest.h"

#include "config.h"

#include <cstdio>

#include <filesystem>
#include <functional>
#include <utility>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

#include "programs/mdrun/mdrun_main.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * MdrunTestFixture
 */

namespace
{

#if GMX_OPENMP || defined(DOXYGEN)
//! Number of OpenMP threads for child mdrun call.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
int g_numOpenMPThreads = 1;
#endif
//! \cond
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
GMX_TEST_OPTIONS(MdrunTestOptions, options)
{
    GMX_UNUSED_VALUE(options);
#if GMX_OPENMP
    options->addOption(
            IntegerOption("ntomp").store(&g_numOpenMPThreads).description("Number of OpenMP threads for child mdrun calls"));
#endif
}
//! \endcond

} // namespace

SimulationRunner::SimulationRunner(TestFileManager* fileManager) :
    fullPrecisionTrajectoryFileName_(fileManager->getTemporaryFilePath(".trr").string()),
    groOutputFileName_(fileManager->getTemporaryFilePath(".gro").string()),
    cptOutputFileName_(fileManager->getTemporaryFilePath(".cpt").string()),
    mdpOutputFileName_(fileManager->getTemporaryFilePath("output.mdp").string()),
    tprFileName_(fileManager->getTemporaryFilePath(".tpr").string()),
    logFileName_(fileManager->getTemporaryFilePath(".log").string()),
    edrFileName_(fileManager->getTemporaryFilePath(".edr").string()),
    mtxFileName_(fileManager->getTemporaryFilePath(".mtx").string()),

    nsteps_(-2),
    maxwarn_(0),
    mdpSource_(SimulationRunnerMdpSource::Undefined),
    fileManager_(*fileManager)
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");

    // It would be better to also detect this in a thread-MPI build,
    // but there is no way to do that currently, and it is also not a
    // problem for such a build. Any code based on such an invalid
    // test fixture will be found in CI testing, however.
    GMX_RELEASE_ASSERT(MdrunTestFixtureBase::s_communicator != MPI_COMM_NULL,
                       "SimulationRunner may only be used from a test fixture that inherits from "
                       "MdrunTestFixtureBase");
#endif
}

// TODO The combination of defaulting to Verlet cut-off scheme, NVE,
// and verlet-buffer-tolerance = -1 gives a grompp error. If we keep
// things that way, this function should be renamed. For now,
// we use the Verlet scheme and hard-code a tolerance.
// TODO There is possible outstanding unexplained behaviour of mdp
// input parsing e.g. Issue #2074, so this particular set of mdp
// contents is also tested with GetIrTest in gmxpreprocess-test.
void SimulationRunner::useEmptyMdpFile()
{
    useStringAsMdpFile("");
}

void SimulationRunner::useStringAsMdpFile(const char* mdpString)
{
    useStringAsMdpFile(std::string(mdpString));
}

void SimulationRunner::useStringAsMdpFile(const std::string& mdpString)
{
    GMX_RELEASE_ASSERT(mdpSource_ != SimulationRunnerMdpSource::File,
                       "Cannot mix .mdp file from database with options set via string.");
    mdpSource_        = SimulationRunnerMdpSource::String;
    mdpInputContents_ = mdpString;
}

void SimulationRunner::useStringAsNdxFile(const char* ndxString) const
{
    gmx::TextWriter::writeFileFromString(ndxFileName_, ndxString);
}

void SimulationRunner::useTopG96AndNdxFromDatabase(const std::string& name)
{
    topFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".top").string();
    groFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".g96").string();
    ndxFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".ndx").string();
}

void SimulationRunner::useTopGroAndNdxFromDatabase(const std::string& name)
{
    topFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".top").string();
    groFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".gro").string();
    ndxFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".ndx").string();
}

void SimulationRunner::useGroFromDatabase(const char* name)
{
    groFileName_ =
            gmx::test::TestFileManager::getInputFilePath((std::string(name) + ".gro").c_str()).string();
}

void SimulationRunner::useNdxFromDatabase(const std::string& name)
{
    ndxFileName_ = gmx::test::TestFileManager::getInputFilePath(name + ".ndx").string();
}

void SimulationRunner::useTopGroAndMdpFromFepTestDatabase(const std::string& name)
{
    GMX_RELEASE_ASSERT(mdpSource_ != SimulationRunnerMdpSource::String,
                       "Cannot mix .mdp file from database with options set via string.");
    mdpSource_ = SimulationRunnerMdpSource::File;
    topFileName_ =
            gmx::test::TestFileManager::getInputFilePath("freeenergy/" + name + "/topol.top").string();
    groFileName_ =
            gmx::test::TestFileManager::getInputFilePath("freeenergy/" + name + "/conf.gro").string();
    mdpFileName_ =
            gmx::test::TestFileManager::getInputFilePath("freeenergy/" + name + "/grompp.mdp").string();
}

void SimulationRunner::setMaxWarn(int maxwarn)
{
    maxwarn_ = maxwarn;
}

int SimulationRunner::callGromppOnThisRank(const CommandLine& callerRef)
{
    std::string mdpInputFileName;
    if (mdpSource_ == SimulationRunnerMdpSource::File)
    {
        mdpInputFileName = mdpFileName_;
    }
    else
    {
        mdpInputFileName = fileManager_.getTemporaryFilePath("input.mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpInputContents_);
    }

    CommandLine caller;
    caller.append("grompp");
    caller.merge(callerRef);
    caller.addOption("-f", mdpInputFileName);
    if (!ndxFileName_.empty())
    {
        caller.addOption("-n", ndxFileName_);
    }
    caller.addOption("-p", topFileName_);
    caller.addOption("-c", groFileName_);
    caller.addOption("-r", groFileName_);

    caller.addOption("-po", mdpOutputFileName_);
    caller.addOption("-o", tprFileName_);
    if (maxwarn_ != 0)
    {
        caller.addOption("-maxwarn", maxwarn_);
    }

    return gmx_grompp(caller.argc(), caller.argv());
}

int SimulationRunner::callGromppOnThisRank()
{
    return callGromppOnThisRank(CommandLine());
}

int SimulationRunner::callGrompp(const CommandLine& callerRef)
{
    int returnValue = 0;
#if GMX_LIB_MPI
    // When compiled with external MPI, we're trying to run mdrun with
    // MPI, but we need to make sure that we only do grompp on one
    // rank
    if (0 == gmx_node_rank())
#endif
    {
        returnValue = callGromppOnThisRank(callerRef);
    }
#if GMX_LIB_MPI
    // Make sure rank zero has written the .tpr file before other
    // ranks try to read it. Thread-MPI and serial do this just fine
    // on their own.
    MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
    return returnValue;
}

int SimulationRunner::callGrompp()
{
    return callGrompp(CommandLine());
}

int SimulationRunner::changeTprNsteps(int nsteps) const
{
    CommandLine caller;
    caller.append("convert-tpr");
    caller.addOption("-nsteps", nsteps);
    // Because the operation is to change the .tpr, we replace the
    // file. TODO Do we need to delete an automatic backup?
    caller.addOption("-s", tprFileName_);
    caller.addOption("-o", tprFileName_);

    return gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &caller);
}

int SimulationRunner::callNmeig() const
{
    /* Conforming to style guide by not passing a non-const reference
       to this function. Passing a non-const reference might make it
       easier to write code that incorrectly re-uses callerRef after
       the call to this function. */

    CommandLine caller;
    caller.append("nmeig");
    caller.addOption("-s", tprFileName_);
    caller.addOption("-f", mtxFileName_);
    // Ignore the overall translation and rotation in the
    // first six eigenvectors.
    caller.addOption("-first", "7");
    // No need to check more than a number of output values.
    caller.addOption("-last", "50");
    caller.addOption("-xvg", "none");

    return gmx_nmeig(caller.argc(), caller.argv());
}

int SimulationRunner::callMdrun(const CommandLine& callerRef)
{
    /* Conforming to style guide by not passing a non-const reference
       to this function. Passing a non-const reference might make it
       easier to write code that incorrectly re-uses callerRef after
       the call to this function. */

    CommandLine caller;
    caller.append("mdrun");
    caller.merge(callerRef);
    caller.addOption("-s", tprFileName_);

    caller.addOption("-g", logFileName_);
    caller.addOption("-e", edrFileName_);
    caller.addOption("-mtx", mtxFileName_);
    caller.addOption("-o", fullPrecisionTrajectoryFileName_);
    caller.addOption("-x", reducedPrecisionTrajectoryFileName_);
    if (!dhdlFileName_.empty())
    {
        caller.addOption("-dhdl", dhdlFileName_);
    }
    caller.addOption("-c", groOutputFileName_);
    caller.addOption("-cpo", cptOutputFileName_);

    caller.addOption("-deffnm", fileManager_.getTemporaryFilePath("state").string());

    if (nsteps_ > -2)
    {
        caller.addOption("-nsteps", nsteps_);
    }

#if GMX_THREAD_MPI
    caller.addOption("-ntmpi", getNumberOfTestMpiRanks());
#endif

#if GMX_OPENMP
    caller.addOption("-ntomp", g_numOpenMPThreads);
#endif

    return gmx_mdrun(MdrunTestFixtureBase::s_communicator,
                     *MdrunTestFixtureBase::s_hwinfo,
                     caller.argc(),
                     caller.argv());
}

int SimulationRunner::callMdrun()
{
    return callMdrun(CommandLine());
}

// ====

// static
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
MPI_Comm MdrunTestFixtureBase::s_communicator = MPI_COMM_NULL;
// static
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::unique_ptr<gmx_hw_info_t> MdrunTestFixtureBase::s_hwinfo;

// static
void MdrunTestFixtureBase::SetUpTestSuite()
{
    s_communicator = MPI_COMM_WORLD;
    auto newHwinfo = gmx_detect_hardware(
            PhysicalNodeCommunicator{ s_communicator, gmx_physicalnode_id_hash() }, s_communicator);
    std::swap(s_hwinfo, newHwinfo);
}

// static
void MdrunTestFixtureBase::TearDownTestSuite()
{
    s_hwinfo.reset(nullptr);
}

MdrunTestFixtureBase::MdrunTestFixtureBase()
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");
#endif
}

MdrunTestFixtureBase::~MdrunTestFixtureBase() {}

// ====

MdrunTestFixture::MdrunTestFixture() : runner_(&fileManager_) {}

MdrunTestFixture::~MdrunTestFixture()
{
#if GMX_LIB_MPI
    // fileManager_ should only clean up after all the ranks are done.
    MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
}

int getNumberOfTestOpenMPThreads()
{
#if GMX_OPENMP
    return g_numOpenMPThreads;
#else
    return 1;
#endif
}

} // namespace test
} // namespace gmx
