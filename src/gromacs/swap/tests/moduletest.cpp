/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * \author Carsten Kutzner <ckutzne@gwdg.de>
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
#include "programs/mdrun/mdrun_main.h"
#include "programs/gmx/legacycmainfunctions.h"

namespace gmx
{
namespace test
{


/********************************************************************
 * SwapTestFixture
 */


/*! /brief Number of tMPI threads for child mdrun call */
static int gmx_unused numThreads = 1;
/*! /brief Number of OpenMP threads for child mdrun call */
static int gmx_unused numOpenMPThreads = 1;


SwapTestFixture::SwapTestFixture() :
    topFileName(),
    groInputFileName(),
    xtcFileName(),
    mdpInputFileName(),
    mdpOutputFileName(fileManager_.getTemporaryFilePath("out.mdp")),
    tprFileName(fileManager_.getTemporaryFilePath(".tpr")),
    logFileName(fileManager_.getTemporaryFilePath(".log")),
    ndxFileName(fileManager_.getTemporaryFilePath(".ndx")),
    edrFileName(fileManager_.getTemporaryFilePath(".edr")),
    cptFileName(fileManager_.getTemporaryFilePath(".cpt")),
    swapFileName(fileManager_.getTemporaryFilePath("swap.xvg")),
    groOutputFileName(fileManager_.getTemporaryFilePath(".gro"))
{
#ifdef GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "MPI system not initialized for mdrun tests");
#endif
}

SwapTestFixture::~SwapTestFixture()
{
    // Clean up checkpoint backup file
    std::string backup = fileManager_.getTemporaryFilePath("prev.cpt");

    std::remove(backup.c_str());
}

int
SwapTestFixture::callGrompp()
{
#ifdef GMX_LIB_MPI
    // When compiled with external MPI, only call one instance of the grompp function
    if (0 != gmx_node_rank())
    {
        return 0;
    }
#endif

    CommandLine caller;
    caller.append("grompp");
    caller.addOption("-f", mdpInputFileName);
    caller.addOption("-p", topFileName);
    caller.addOption("-c", groInputFileName);
    caller.addOption("-n", ndxFileName);
    caller.addOption("-po", mdpOutputFileName);
    caller.addOption("-o", tprFileName);

    return gmx_grompp(caller.argc(), caller.argv());
}

int
SwapTestFixture::callMdrun()
{
    CommandLine caller;
    caller.append("mdrun");
    caller.addOption("-s", tprFileName);
    caller.addOption("-g", logFileName);
    caller.addOption("-cpo", cptFileName);
    caller.addOption("-swap", swapFileName);
    caller.addOption("-e", edrFileName);
    caller.addOption("-c", groOutputFileName);

#ifdef GMX_THREAD_MPI
    caller.addOption("-nt", numThreads);
#endif
#ifdef GMX_OPENMP
    caller.addOption("-ntomp", numOpenMPThreads);
#endif

    return gmx_mdrun(caller.argc(), caller.argv());
}

int
SwapTestFixture::callMdrunAppend()
{
    CommandLine caller;
    caller.append("mdrun");
    caller.addOption("-s", tprFileName);
    caller.addOption("-g", logFileName);
    caller.addOption("-swap", swapFileName);
    caller.addOption("-e", edrFileName);
    caller.addOption("-c", groOutputFileName);
    caller.addOption("-cpo", cptFileName);
    caller.addOption("-cpi", cptFileName);
    caller.addOption("-nsteps", "15");
    caller.append("-append");

#ifdef GMX_THREAD_MPI
    caller.addOption("-nt", numThreads);
#endif
#ifdef GMX_OPENMP
    caller.addOption("-ntomp", numOpenMPThreads);
#endif

    return gmx_mdrun(caller.argc(), caller.argv());
}
} // namespace test
} // namespace gmx
