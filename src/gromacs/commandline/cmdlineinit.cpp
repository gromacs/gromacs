/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Implements functions from cmdlineinit.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlineinit.h"

#include <cstring>

#include <memory>

#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/cmdlineprogramcontext.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/datafilefinder.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

namespace
{

//! \addtogroup module_commandline
//! \{

// These never release ownership.
//! Global context instance initialized in initForCommandLine().
std::unique_ptr<CommandLineProgramContext> g_commandLineContext;
//! Global library data file finder that respects GMXLIB.
std::unique_ptr<DataFileFinder>            g_libFileFinder;

/*! \brief
 * Broadcasts command-line arguments to all ranks.
 *
 * MPI does not ensure that command-line arguments would be passed on any
 * other rank than zero, but our code wants to parse them on each rank
 * separately.
 */
void broadcastArguments(int *argc, char ***argv)
{
    if (gmx_node_num() <= 1)
    {
        return;
    }
    gmx_broadcast_world(sizeof(*argc), argc);

    const bool isMaster = (gmx_node_rank() == 0);
    if (!isMaster)
    {
        snew(*argv, *argc+1);
    }
    for (int i = 0; i < *argc; i++)
    {
        int len;
        if (isMaster)
        {
            len = std::strlen((*argv)[i])+1;
        }
        gmx_broadcast_world(sizeof(len), &len);
        if (!isMaster)
        {
            snew((*argv)[i], len);
        }
        gmx_broadcast_world(len, (*argv)[i]);
    }
}

//! \}

}   // namespace

CommandLineProgramContext &initForCommandLine(int *argc, char ***argv)
{
    gmx::init(argc, argv);
    GMX_RELEASE_ASSERT(!g_commandLineContext,
                       "initForCommandLine() calls cannot be nested");
    // TODO: Consider whether the argument broadcast would better be done
    // in CommandLineModuleManager.
    broadcastArguments(argc, argv);
    try
    {
        g_commandLineContext.reset(new CommandLineProgramContext(*argc, *argv));
        setProgramContext(g_commandLineContext.get());
        g_libFileFinder.reset(new DataFileFinder());
        g_libFileFinder->setSearchPathFromEnv("GMXLIB");
        setLibraryFileFinder(g_libFileFinder.get());
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        std::exit(processExceptionAtExit(ex));
    }
    return *g_commandLineContext;
}

void finalizeForCommandLine()
{
    gmx::finalize();
    setLibraryFileFinder(NULL);
    g_libFileFinder.reset();
    setProgramContext(NULL);
    g_commandLineContext.reset();
}

int processExceptionAtExitForCommandLine(const std::exception &ex)
{
    int rc = processExceptionAtExit(ex); // Currently this aborts for real MPI
    finalizeForCommandLine();            // thus this MPI_Finalize doesn't matter.
    return rc;
}

int runCommandLineModule(int argc, char *argv[],
                         ICommandLineModule *module)
{
    return CommandLineModuleManager::runAsMainSingleModule(argc, argv, module);
}

int runCommandLineModule(
        int argc, char *argv[], const char *name, const char *description,
        std::function<std::unique_ptr<ICommandLineOptionsModule>()> factory)
{
    return ICommandLineOptionsModule::runAsMain(
            argc, argv, name, description, factory);
}

} // namespace gmx

int gmx_run_cmain(int argc, char *argv[], int (*mainFunction)(int, char *[]))
{
    return gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, mainFunction);
}
