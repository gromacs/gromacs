/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "broadcaststructs.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/state.h"

template<typename AllocatorType>
static void bcastPaddedRVecVector(MPI_Comm                                     communicator,
                                  gmx::PaddedVector<gmx::RVec, AllocatorType>* v,
                                  int                                          numAtoms)
{
    v->resizeWithPadding(numAtoms);
    nblock_bc(communicator, makeArrayRef(*v));
}

void broadcastStateWithoutDynamics(MPI_Comm communicator,
                                   bool     useDomainDecomposition,
                                   bool     isParallelRun,
                                   t_state* state)
{
    GMX_RELEASE_ASSERT(!useDomainDecomposition,
                       "broadcastStateWithoutDynamics should only be used for special cases "
                       "without domain decomposition");

    if (!isParallelRun)
    {
        return;
    }

    /* Broadcasts the state sizes and flags from the master to all ranks
     * in cr->mpi_comm_mygroup.
     */
    block_bc(communicator, state->natoms);
    block_bc(communicator, state->flags);

    for (int i = 0; i < estNR; i++)
    {
        if (state->flags & (1 << i))
        {
            switch (i)
            {
                case estLAMBDA: nblock_bc(communicator, efptNR, state->lambda.data()); break;
                case estFEPSTATE: block_bc(communicator, state->fep_state); break;
                case estBOX: block_bc(communicator, state->box); break;
                case estX: bcastPaddedRVecVector(communicator, &state->x, state->natoms); break;
                default:
                    GMX_RELEASE_ASSERT(false,
                                       "The state has a dynamic entry, while no dynamic entries "
                                       "should be present");
                    break;
            }
        }
    }
}

static void bc_tpxheader(MPI_Comm communicator, TpxFileHeader* tpx)
{
    block_bc(communicator, tpx->bIr);
    block_bc(communicator, tpx->bBox);
    block_bc(communicator, tpx->bTop);
    block_bc(communicator, tpx->bX);
    block_bc(communicator, tpx->bV);
    block_bc(communicator, tpx->bF);
    block_bc(communicator, tpx->natoms);
    block_bc(communicator, tpx->ngtc);
    block_bc(communicator, tpx->lambda);
    block_bc(communicator, tpx->fep_state);
    block_bc(communicator, tpx->sizeOfTprBody);
    block_bc(communicator, tpx->fileVersion);
    block_bc(communicator, tpx->fileGeneration);
    block_bc(communicator, tpx->isDouble);
}

static void bc_tprCharBuffer(MPI_Comm communicator, bool isMasterRank, std::vector<char>* charBuffer)
{
    int elements = charBuffer->size();
    block_bc(communicator, elements);

    nblock_abc(isMasterRank, communicator, elements, charBuffer);
}

void init_parallel(MPI_Comm                    communicator,
                   bool                        isMasterRank,
                   t_inputrec*                 inputrec,
                   gmx_mtop_t*                 mtop,
                   PartialDeserializedTprFile* partialDeserializedTpr)
{
    bc_tpxheader(communicator, &partialDeserializedTpr->header);
    bc_tprCharBuffer(communicator, isMasterRank, &partialDeserializedTpr->body);
    if (!isMasterRank)
    {
        completeTprDeserialization(partialDeserializedTpr, inputrec, mtop);
    }
}
