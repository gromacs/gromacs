/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements classes in moduletest.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "moduletest.h"

#include "config.h"

#include <cstdio>

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/textwriter.h"
#include "programs/mdrun/mdrun_main.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

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
int g_numOpenMPThreads = 1;
#endif
//! \cond
GMX_TEST_OPTIONS(MdrunTestOptions, options)
{
    GMX_UNUSED_VALUE(options);
#if GMX_OPENMP
    options->addOption(IntegerOption("nt_omp").store(&g_numOpenMPThreads)
                           .description("Number of OpenMP threads for child mdrun calls"));
#endif
}
//! \endcond

}

SimulationRunner::SimulationRunner(TestFileManager *fileManager) :
    topFileName_(),
    groFileName_(),
    fullPrecisionTrajectoryFileName_(),
    ndxFileName_(),
    mdpInputFileName_(fileManager->getTemporaryFilePath("input.mdp")),
    mdpOutputFileName_(fileManager->getTemporaryFilePath("output.mdp")),
    tprFileName_(fileManager->getTemporaryFilePath(".tpr")),
    logFileName_(fileManager->getTemporaryFilePath(".log")),
    edrFileName_(fileManager->getTemporaryFilePath(".edr")),
    nsteps_(-2),
    fileManager_(*fileManager)
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");
#endif
}

// TODO The combination of defaulting to Verlet cut-off scheme, NVE,
// and verlet-buffer-tolerance = -1 gives a grompp error. If we keep
// things that way, this function should be renamed. For now,
// force the use of the group scheme.
// TODO There is possible outstanding unexplained behaviour of mdp
// input parsing e.g. Redmine 2074, so this particular set of mdp
// contents is also tested with GetIrTest in gmxpreprocess-test.
void
SimulationRunner::useEmptyMdpFile()
{
    // TODO When removing the group scheme, update actual and potential users of useEmptyMdpFile
    useStringAsMdpFile("cutoff-scheme = Group\n");
}

void
SimulationRunner::useStringAsMdpFile(const char *mdpString)
{
    useStringAsMdpFile(std::string(mdpString));
}

void
SimulationRunner::useStringAsMdpFile(const std::string &mdpString)
{
    gmx::TextWriter::writeFileFromString(mdpInputFileName_, mdpString);
}

void
SimulationRunner::useStringAsNdxFile(const char *ndxString)
{
    gmx::TextWriter::writeFileFromString(ndxFileName_, ndxString);
}

void
SimulationRunner::useTopGroAndNdxFromDatabase(const char *name)
{
    topFileName_ = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    groFileName_ = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    ndxFileName_ = fileManager_.getInputFilePath((std::string(name) + ".ndx").c_str());
}

void
SimulationRunner::useGroFromDatabase(const char *name)
{
    groFileName_ = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
}

int
SimulationRunner::callGromppOnThisRank(const CommandLine &callerRef)
{
    CommandLine caller;
    caller.append("grompp");
    caller.merge(callerRef);
    caller.addOption("-f", mdpInputFileName_);
    caller.addOption("-n", ndxFileName_);
    caller.addOption("-p", topFileName_);
    caller.addOption("-c", groFileName_);
    caller.addOption("-r", groFileName_);

    caller.addOption("-po", mdpOutputFileName_);
    caller.addOption("-o", tprFileName_);

    return gmx_grompp(caller.argc(), caller.argv());
}

int
SimulationRunner::callGromppOnThisRank()
{
    return callGromppOnThisRank(CommandLine());
}

int
SimulationRunner::callGrompp(const CommandLine &callerRef)
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
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    return returnValue;
}

int
SimulationRunner::callGrompp()
{
    return callGrompp(CommandLine());
}

int
SimulationRunner::callMdrun(const CommandLine &callerRef)
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
    caller.addOption("-o", fullPrecisionTrajectoryFileName_);
    caller.addOption("-x", reducedPrecisionTrajectoryFileName_);

    caller.addOption("-deffnm", fileManager_.getTemporaryFilePath("state"));

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

    return gmx_mdrun(caller.argc(), caller.argv());
}

int
SimulationRunner::callMdrun()
{
    return callMdrun(CommandLine());
}

// ====

MdrunTestFixtureBase::MdrunTestFixtureBase()
{
#if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");
#endif
}

MdrunTestFixtureBase::~MdrunTestFixtureBase()
{
}

// ====

MdrunTestFixture::MdrunTestFixture() : runner_(&fileManager_)
{
}

MdrunTestFixture::~MdrunTestFixture()
{
}

} // namespace test
} // namespace gmx
