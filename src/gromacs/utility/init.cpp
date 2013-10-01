/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements functions from init.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gromacs/utility/init.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstring>

#ifdef GMX_LIB_MPI
#include "gromacs/utility/gmxmpi.h"
#endif

#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/types/commrec.h"

#include "gromacs/utility/programinfo.h"

namespace gmx
{

#ifdef GMX_LIB_MPI
namespace
{

void broadcastArguments(const t_commrec *cr, int *argc, char ***argv)
{
    gmx_bcast(sizeof(*argc), argc, cr);

    if (!MASTER(cr))
    {
        snew(*argv, *argc+1);
    }
    for (int i = 0; i < *argc; i++)
    {
        int len;
        if (MASTER(cr))
        {
            len = std::strlen((*argv)[i])+1;
        }
        gmx_bcast(sizeof(len), &len, cr);
        if (!MASTER(cr))
        {
            snew((*argv)[i], len);
        }
        gmx_bcast(len, (*argv)[i], cr);
    }
}

} // namespace
#endif

ProgramInfo &ProgramInitializer::init(const char *realBinaryName, int *argc, char ***argv)
{
#ifdef GMX_LIB_MPI
    int isInitialized = 0;
    MPI_Initialized(&isInitialized);
    if (isInitialized)
    {
        if (0 == _initializationCounter)
        {
            // Some other code has already initialized MPI, so bump the counter so that
            // we know not to finalize MPI ourselves later.
            _initializationCounter++;
        }
    }
    else
    {
#ifdef GMX_FAHCORE
        (void) fah_MPI_Init(argc, argv);
#else
        (void) MPI_Init(argc, argv);
#endif
    }
    // Bump the counter to record this initialization event
    _initializationCounter++;

    // TODO: Rewrite this to not use t_commrec once there is clarity on
    // the approach for MPI in C++ code.
    // TODO: Consider whether the argument broadcast would better be done
    // in CommandLineModuleManager.
    t_commrec cr;
    std::memset(&cr, 0, sizeof(cr));

    gmx_fill_commrec_from_mpi(&cr);
    if (PAR(&cr))
    {
        broadcastArguments(&cr, argc, argv);
    }
#endif
    return ProgramInfo::init(realBinaryName, *argc, *argv);
}

ProgramInfo &ProgramInitializer::init(int *argc, char ***argv)
{
    return init(NULL, argc, argv);
}

void ProgramInitializer::finalize()
{
#ifdef GMX_LIB_MPI
    assert(0 < _initializationCounter);
    // Bump the counter to record this finalization event
    _initializationCounter--;

    if (0 == _initializationCounter)
    {
        /* We sync the processes here to try to avoid problems
         * with buggy MPI implementations that could cause
         * unfinished processes to terminate.
         */
        MPI_Barrier(MPI_COMM_WORLD);

        /* Apparently certain mpich implementations cause problems
         * with MPI_Finalize. In that case comment out MPI_Finalize.
         */
        MPI_Finalize();
    }
#endif
    return;
}

int ProgramInitializer::_initializationCounter = 0;

} // namespace gmx
