/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Defines functionality for communicators across physical nodes.
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "physicalnodecommunicator.h"

#include "config.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{

void MPI_Comm_free_wrapper(MPI_Comm* comm)
{
#if GMX_MPI
    // With thread-MPI *comm is shared between ranks which causes issues with
    // freeing. But all thread-mpi data is anyhow freed in tMPI_Finalize()
    // and in practice *comm is always MPI_COMM_WORLD with thread-MPI.
    // Only the thread-affinity test code uses *comm != MPI_COMM_WORLD.
    if (!GMX_THREAD_MPI)
    {
        MPI_Comm_free(comm);
    }
#else
    GMX_UNUSED_VALUE(comm);
#endif
}

PhysicalNodeCommunicator::PhysicalNodeCommunicator(MPI_Comm world, int physicalNodeId)
{
#if GMX_MPI
    int isInitialized;
    MPI_Initialized(&isInitialized);
    if (isInitialized)
    {
        int sizeOfWorld;
        MPI_Comm_size(world, &sizeOfWorld);
        if (sizeOfWorld > 1)
        {
            int rankWithinWorld;
            MPI_Comm_rank(world, &rankWithinWorld);
            MPI_Comm_split(world, physicalNodeId, rankWithinWorld, &comm_);
            auto ptr = MPI_Comm_ptr(&comm_);
            commGuard_.swap(ptr);
            MPI_Comm_size(comm_, &size_);
            MPI_Comm_rank(comm_, &rank_);
        }
        else
        {
            // Handle this trivial case separately, because thread-MPI
            // doesn't have a valid communicator when there is only
            // one rank.
            comm_ = world;
            size_ = 1;
            rank_ = 0;
        }
    }
    else
    {
        comm_ = MPI_COMM_NULL;
        size_ = 1;
        rank_ = 0;
    }
#else
    // Trivial case when there is no MPI support or not initialized
    GMX_UNUSED_VALUE(world);
    GMX_UNUSED_VALUE(physicalNodeId);
    comm_ = nullptr;
    size_ = 1;
    rank_ = 0;
#endif
}

void PhysicalNodeCommunicator::barrier() const
{
#if GMX_MPI
    if (size_ > 1)
    {
        MPI_Barrier(comm_);
    }
#else
    // Nothing to do
#endif
}

} // namespace gmx
