/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * \ingroup module_mdrun
 */
#include "moduletest.h"

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "testutils/integrationtests.h"
#include "testutils/testoptions.h"
#include "testutils/cmdlinetest.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/gmxpreprocess/grompp.h"
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

#if defined(GMX_THREAD_MPI) || defined(DOXYGEN)
//! Number of tMPI threads for child mdrun call.
int g_numThreads = 1;
#endif
#if defined(GMX_OPENMP) || defined(DOXYGEN)
//! Number of OpenMP threads for child mdrun call.
int g_numOpenMPThreads = 1;
#endif
//! \cond
GMX_TEST_OPTIONS(MdrunTestOptions, options)
{
    GMX_UNUSED_VALUE(options);
#ifdef GMX_THREAD_MPI
    options->addOption(IntegerOption("nt").store(&g_numThreads)
                           .description("Number of thread-MPI threads/ranks for child mdrun call"));
#endif
#ifdef GMX_OPENMP
    options->addOption(IntegerOption("nt_omp").store(&g_numOpenMPThreads)
                           .description("Number of OpenMP threads for child mdrun call"));
#endif
}
//! \endcond

}

MdrunTestFixture::MdrunTestFixture() :
    topFileName(),
    groFileName(),
    fullPrecisionTrajectoryFileName(),
    ndxFileName(),
    mdpInputFileName(fileManager_.getTemporaryFilePath("input.mdp")),
    mdpOutputFileName(fileManager_.getTemporaryFilePath("output.mdp")),
    tprFileName(fileManager_.getTemporaryFilePath(".tpr")),
    logFileName(fileManager_.getTemporaryFilePath(".log")),
    edrFileName(fileManager_.getTemporaryFilePath(".edr")),
    nsteps(-2)
{
#ifdef GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");
#endif
}

MdrunTestFixture::~MdrunTestFixture()
{
}

// TODO The combination of defaulting to Verlet cut-off scheme, NVE,
// and verlet-buffer-tolerance = -1 gives a grompp error. If we keep
// things that way, this function should be renamed. For now,
// force the use of the group scheme.
void
MdrunTestFixture::useEmptyMdpFile()
{
    useStringAsMdpFile("cutoff-scheme = Group\n");
}

void
MdrunTestFixture::useStringAsMdpFile(const char *mdpString)
{
    useStringAsMdpFile(std::string(mdpString));
}

void
MdrunTestFixture::useStringAsMdpFile(const std::string &mdpString)
{
    gmx::File::writeFileFromString(mdpInputFileName, mdpString);
}

void
MdrunTestFixture::useStringAsNdxFile(const char *ndxString)
{
    gmx::File::writeFileFromString(ndxFileName, ndxString);
}

void
MdrunTestFixture::useTopGroAndNdxFromDatabase(const char *name)
{
    topFileName = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    groFileName = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    ndxFileName = fileManager_.getInputFilePath((std::string(name) + ".ndx").c_str());
}

int
MdrunTestFixture::callGromppOnThisRank()
{
    CommandLine caller;
    caller.append("grompp");
    caller.addOption("-f", mdpInputFileName);
    caller.addOption("-n", ndxFileName);
    caller.addOption("-p", topFileName);
    caller.addOption("-c", groFileName);

    caller.addOption("-po", mdpOutputFileName);
    caller.addOption("-o", tprFileName);

    return gmx_grompp(caller.argc(), caller.argv());
}

int
MdrunTestFixture::callGrompp()
{
#ifdef GMX_LIB_MPI
    // When compiled with external MPI, only call one instance of the
    // grompp function
    if (0 != gmx_node_rank())
    {
        return 0;
    }
#endif
    return callGromppOnThisRank();
}

int
MdrunTestFixture::callMdrun(const CommandLine &callerRef)
{
    /* Conforming to style guide by not passing a non-const reference
       to this function. Passing a non-const reference might make it
       easier to write code that incorrectly re-uses callerRef after
       the call to this function. */
    CommandLine caller(callerRef);
    caller.addOption("-s", tprFileName);

    caller.addOption("-g", logFileName);
    caller.addOption("-e", edrFileName);
    caller.addOption("-o", fullPrecisionTrajectoryFileName);
    caller.addOption("-x", reducedPrecisionTrajectoryFileName);

    caller.addOption("-deffnm", fileManager_.getTemporaryFilePath("state"));

    if (nsteps > -2)
    {
        caller.addOption("-nsteps", nsteps);
    }

#ifdef GMX_THREAD_MPI
    caller.addOption("-nt", g_numThreads);
#endif
#ifdef GMX_OPENMP
    caller.addOption("-ntomp", g_numOpenMPThreads);
#endif

    return gmx_mdrun(caller.argc(), caller.argv());
}

int
MdrunTestFixture::callMdrun()
{
    CommandLine caller;
    caller.append("mdrun");
    return callMdrun(caller);
}

} // namespace test
} // namespace gmx
