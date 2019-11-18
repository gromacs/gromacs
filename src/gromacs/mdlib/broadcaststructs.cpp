/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/state.h"

template<typename AllocatorType>
static void bcastPaddedRVecVector(const t_commrec* cr, gmx::PaddedVector<gmx::RVec, AllocatorType>* v, int numAtoms)
{
    v->resizeWithPadding(numAtoms);
    nblock_bc(cr, makeArrayRef(*v));
}

void broadcastStateWithoutDynamics(const t_commrec* cr, t_state* state)
{
    GMX_RELEASE_ASSERT(!DOMAINDECOMP(cr),
                       "broadcastStateWithoutDynamics should only be used for special cases "
                       "without domain decomposition");

    if (!PAR(cr))
    {
        return;
    }

    /* Broadcasts the state sizes and flags from the master to all ranks
     * in cr->mpi_comm_mygroup.
     */
    block_bc(cr, state->natoms);
    block_bc(cr, state->flags);

    for (int i = 0; i < estNR; i++)
    {
        if (state->flags & (1 << i))
        {
            switch (i)
            {
                case estLAMBDA: nblock_bc(cr, efptNR, state->lambda.data()); break;
                case estFEPSTATE: block_bc(cr, state->fep_state); break;
                case estBOX: block_bc(cr, state->box); break;
                case estX: bcastPaddedRVecVector(cr, &state->x, state->natoms); break;
                default:
                    GMX_RELEASE_ASSERT(false,
                                       "The state has a dynamic entry, while no dynamic entries "
                                       "should be present");
                    break;
            }
        }
    }
}

static void bc_tpxheader(const t_commrec* cr, TpxFileHeader* tpx)
{
    block_bc(cr, tpx->bIr);
    block_bc(cr, tpx->bBox);
    block_bc(cr, tpx->bTop);
    block_bc(cr, tpx->bX);
    block_bc(cr, tpx->bV);
    block_bc(cr, tpx->bF);
    block_bc(cr, tpx->natoms);
    block_bc(cr, tpx->ngtc);
    block_bc(cr, tpx->lambda);
    block_bc(cr, tpx->fep_state);
    block_bc(cr, tpx->sizeOfTprBody);
    block_bc(cr, tpx->fileVersion);
    block_bc(cr, tpx->fileGeneration);
    block_bc(cr, tpx->isDouble);
}

static void bc_tprCharBuffer(const t_commrec* cr, std::vector<char>* charBuffer)
{
    int elements = charBuffer->size();
    block_bc(cr, elements);

    nblock_abc(cr, elements, charBuffer);
}

void init_parallel(t_commrec* cr, t_inputrec* inputrec, gmx_mtop_t* mtop, PartialDeserializedTprFile* partialDeserializedTpr)
{
    bc_tpxheader(cr, &partialDeserializedTpr->header);
    bc_tprCharBuffer(cr, &partialDeserializedTpr->body);
    if (!MASTER(cr))
    {
        completeTprDeserialization(partialDeserializedTpr, inputrec, mtop);
    }
}
