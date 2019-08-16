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
 *
 * \brief Implements the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "multisim.h"

#include "config.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

gmx_multisim_t *init_multisystem(MPI_Comm                         comm,
                                 gmx::ArrayRef<const std::string> multidirs)
{
    gmx_multisim_t *ms;
#if GMX_MPI
    MPI_Group       mpi_group_world;
    int            *rank;
#endif

    if (multidirs.empty())
    {
        return nullptr;
    }

    if (!GMX_LIB_MPI && !multidirs.empty())
    {
        gmx_fatal(FARGS, "mdrun -multidir is only supported when GROMACS has been "
                  "configured with a proper external MPI library.");
    }

    if (multidirs.size() == 1)
    {
        /* NOTE: It would be nice if this special case worked, but this requires checks/tests. */
        gmx_fatal(FARGS, "To run mdrun in multiple simulation mode, more then one "
                  "actual simulation is required. The single simulation case is not supported.");
    }

#if GMX_MPI
    int numRanks;
    MPI_Comm_size(comm, &numRanks);
    if (numRanks % multidirs.size() != 0)
    {
        gmx_fatal(FARGS, "The number of ranks (%d) is not a multiple of the number of simulations (%td)", numRanks, multidirs.ssize());
    }

    int numRanksPerSim = numRanks/multidirs.size();
    int rankWithinComm;
    MPI_Comm_rank(comm, &rankWithinComm);

    if (debug)
    {
        fprintf(debug, "We have %td simulations, %d ranks per simulation, local simulation is %d\n", multidirs.ssize(), numRanksPerSim, rankWithinComm/numRanksPerSim);
    }

    ms       = new gmx_multisim_t;
    ms->nsim = multidirs.size();
    ms->sim  = rankWithinComm/numRanksPerSim;
    /* Create a communicator for the master nodes */
    snew(rank, ms->nsim);
    for (int i = 0; i < ms->nsim; i++)
    {
        rank[i] = i*numRanksPerSim;
    }
    MPI_Comm_group(comm, &mpi_group_world);
    MPI_Group_incl(mpi_group_world, ms->nsim, rank, &ms->mpi_group_masters);
    sfree(rank);
    MPI_Comm_create(MPI_COMM_WORLD, ms->mpi_group_masters,
                    &ms->mpi_comm_masters);

#if !MPI_IN_PLACE_EXISTS
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(ms->mpb, 1);
    ms->mpb->ibuf        = NULL;
    ms->mpb->libuf       = NULL;
    ms->mpb->fbuf        = NULL;
    ms->mpb->dbuf        = NULL;
    ms->mpb->ibuf_alloc  = 0;
    ms->mpb->libuf_alloc = 0;
    ms->mpb->fbuf_alloc  = 0;
    ms->mpb->dbuf_alloc  = 0;
#endif

    // TODO This should throw upon error
    gmx_chdir(multidirs[ms->sim].c_str());
#else
    GMX_UNUSED_VALUE(comm);
    ms = nullptr;
#endif

    return ms;
}

void done_multisim(gmx_multisim_t *ms)
{
    if (nullptr != ms)
    {
        done_mpi_in_place_buf(ms->mpb);
        delete ms;
    }
}
