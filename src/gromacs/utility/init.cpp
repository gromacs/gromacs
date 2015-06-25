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
 * Implements functions from init.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "init.h"

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#ifdef GMX_LIB_MPI
#include "gromacs/utility/gmxmpi.h"
#endif

namespace gmx
{

namespace
{
#ifdef GMX_LIB_MPI
//! Maintains global counter of attempts to initialize MPI
int g_initializationCounter = 0;
#endif
}

void init(int *argc, char ***argv)
{
#ifdef GMX_LIB_MPI
    int isInitialized = 0, isFinalized = 0;
    MPI_Finalized(&isFinalized);
    GMX_RELEASE_ASSERT(!isFinalized, "Invalid attempt to initialize MPI after finalization");
    MPI_Initialized(&isInitialized);
    if (isInitialized)
    {
        if (0 == g_initializationCounter)
        {
            // Some other code has already initialized MPI, so bump the counter so that
            // we know not to finalize MPI ourselves later.
            g_initializationCounter++;
        }
    }
    else
    {
#ifdef GMX_FAHCORE
        fah_MPI_Init(argc, argv);
#else
#    ifdef GMX_OPENMP
        /* Formally we need to use MPI_Init_thread and ask for MPI_THREAD_FUNNELED
         * level of thread support when using OpenMP. However, in practice we
         * have never seen any problems with just using MPI_Init(), and some MPI
         * versions that only return MPI_THREAD_SINGLE as their support level
         * still return an OK initialization code when we request MPI_THREAD_FUNNELED.
         *
         * To avoid requiring users to recompile and install new libraries for something
         * that probably isn't a problem, we don't make a fuss about it (no warning)
         * if MPI_Init_thread() returns ok. If it does not return OK we
         * call a second init call without thread support, but since the library
         * apparently tried to prevent this we do issue a warning in this case.
         *
         * support, and if that doesn't work we simply use the simple MPI_Init().
         * Note that MPICH does not allow us to call MPI_Query_thread() before
         * we have initialized the library, and calling MPI_Init twice could
         * be a bit fragile. However, some implementations also advice against
         * doing anything else after calling MPI_Finalize, so at the end of the
         * day we'll have to cross our fingers...
         */
        int provided, rc;
        rc = MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
        if (rc != 0 && provided < MPI_THREAD_FUNNELED)
        {
            gmx_warning("GROMACS was compiled with OpenMP support, but there is no thread support in the MPI library. Keep your fingers crossed.");
            MPI_Init(argc, argv);
        }
#    else
        MPI_Init(argc, argv);
#    endif
#endif
    }
    // Bump the counter to record this initialization event
    g_initializationCounter++;

#else
    GMX_UNUSED_VALUE(argc);
    GMX_UNUSED_VALUE(argv);
#endif
}

void finalize()
{
#ifdef GMX_LIB_MPI
    GMX_RELEASE_ASSERT(0 < g_initializationCounter, "Excess attempt to finalize MPI");
    // Bump the counter to record this finalization event
    g_initializationCounter--;

    if (0 == g_initializationCounter)
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
}

} // namespace gmx
