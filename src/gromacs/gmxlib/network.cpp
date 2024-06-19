/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "network.h"

#include "config.h"

#include <cctype>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <limits>

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

std::unique_ptr<t_commrec> init_commrec(MPI_Comm communicator)
{
    int rankInCommunicator, sizeOfCommunicator;
#if GMX_MPI
#    if GMX_LIB_MPI
    GMX_RELEASE_ASSERT(gmx_mpi_initialized(), "Must have initialized MPI before building commrec");
#    endif
    MPI_Comm_rank(communicator, &rankInCommunicator);
    MPI_Comm_size(communicator, &sizeOfCommunicator);
#else
    GMX_UNUSED_VALUE(communicator);
    rankInCommunicator = 0;
    sizeOfCommunicator = 1;
#endif

    std::unique_ptr<t_commrec> cr = std::make_unique<t_commrec>();

    cr->mpiDefaultCommunicator    = communicator;
    cr->sizeOfDefaultCommunicator = sizeOfCommunicator;
    cr->rankInDefaultCommunicator = rankInCommunicator;

    // For now, we want things to go horribly wrong if this is used too early...
    // TODO: Remove when communicators are removed from commrec (#2395)
    cr->nnodes                    = -1;
    cr->sizeOfMyGroupCommunicator = -1;
    cr->nodeid                    = -1;
    cr->sim_nodeid                = -1;
    cr->mpi_comm_mysim            = MPI_COMM_NULL;
    cr->mpi_comm_mygroup          = MPI_COMM_NULL;

    // TODO cr->duty should not be initialized here
    cr->duty = (DUTY_PP | DUTY_PME);

    return cr;
}

void gmx_setup_nodecomm(FILE gmx_unused* fplog, t_commrec* cr)
{
    gmx_nodecomm_t* nc;

    /* Many MPI implementations do not optimize MPI_Allreduce
     * (and probably also other global communication calls)
     * for multi-core nodes connected by a network.
     * We can optimize such communication by using one MPI call
     * within each node and one between the nodes.
     * For MVAPICH2 and Intel MPI this reduces the time for
     * the global_stat communication by 25%
     * for 2x2-core 3 GHz Woodcrest connected by mixed DDR/SDR Infiniband.
     * B. Hess, November 2007
     */

    nc = &cr->nc;

    nc->bUse = FALSE;
#if !GMX_THREAD_MPI
#    if GMX_MPI
    int n, rank;

    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.
    MPI_Comm_size(cr->mpi_comm_mygroup, &n);
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);

    int nodehash = gmx_physicalnode_id_hash();

    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: splitting communicator of size %d\n", n);
    }


    /* The intra-node communicator, split on node number */
    MPI_Comm_split(cr->mpi_comm_mygroup, nodehash, rank, &nc->comm_intra);
    MPI_Comm_rank(nc->comm_intra, &nc->rank_intra);
    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: node ID %d rank within node %d\n", rank, nc->rank_intra);
    }
    /* The inter-node communicator, split on rank_intra.
     * We actually only need the one for rank=0,
     * but it is easier to create them all.
     */
    MPI_Comm_split(cr->mpi_comm_mygroup, nc->rank_intra, rank, &nc->comm_inter);
    /* Check if this really created two step communication */
    int ng, ni;

    MPI_Comm_size(nc->comm_inter, &ng);
    MPI_Comm_size(nc->comm_intra, &ni);
    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: groups %d, my group size %d\n", ng, ni);
    }

    if (getenv("GMX_NO_NODECOMM") == nullptr && ((ng > 1 && ng < n) || (ni > 1 && ni < n)))
    {
        nc->bUse = TRUE;
        if (fplog)
        {
            fprintf(fplog,
                    "Using two step summing over %d groups of on average %.1f ranks\n\n",
                    ng,
                    real(n) / real(ng));
        }
        if (nc->rank_intra > 0)
        {
            MPI_Comm_free(&nc->comm_inter);
        }
    }
    else
    {
        /* One group or all processes in a separate group, use normal summing */
        MPI_Comm_free(&nc->comm_inter);
        MPI_Comm_free(&nc->comm_intra);
        if (debug)
        {
            fprintf(debug,
                    "In gmx_setup_nodecomm: not unsing separate inter- and intra-node "
                    "communicators.\n");
        }
    }
#    endif
#else
    /* tMPI runs only on a single node so just use the nodeid */
    nc->rank_intra = cr->nodeid;
#endif
}

void gmx_barrier(MPI_Comm gmx_unused communicator)
{
    if (communicator == MPI_COMM_NULL)
    {
        return;
    }
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_barrier");
#else
    MPI_Barrier(communicator);
#endif
}

void gmx_bcast(std::size_t gmx_unused nbytes, void gmx_unused* b, MPI_Comm gmx_unused communicator)
{
    // Without MPI we have a single rank, so bcast is a no-op
#if GMX_MPI
    constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();
    char*                 bytePtr      = reinterpret_cast<char*>(b);
    for (std::size_t written = 0, remain = nbytes; remain > 0;)
    {
        std::size_t chunk = std::min(remain, maxSignedInt);
        MPI_Bcast(bytePtr + written, chunk, MPI_BYTE, 0, communicator);
        written += chunk;
        remain -= chunk;
    }
#endif
}

void gmx_sumd(std::size_t gmx_unused nr, double gmx_unused r[], const t_commrec gmx_unused* cr)
{
    // Without MPI we have a single rank, so sum is a no-op
#if GMX_MPI
    if (cr->sizeOfMyGroupCommunicator == 1)
    {
        return;
    }

    constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();
    if (cr->nc.bUse)
    {
        if (cr->nc.rank_intra == 0)
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                /* Use two step summing. */
                MPI_Reduce(MPI_IN_PLACE, r + written, chunk, MPI_DOUBLE, MPI_SUM, 0, cr->nc.comm_intra);
                /* Sum the roots of the internal (intra) buffers. */
                MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_DOUBLE, MPI_SUM, cr->nc.comm_inter);
                written += chunk;
                remain -= chunk;
            }
        }
        else
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                /* This is here because of the silly MPI specification
                   that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
                MPI_Reduce(r + written, nullptr, chunk, MPI_DOUBLE, MPI_SUM, 0, cr->nc.comm_intra);
                written += chunk;
                remain -= chunk;
            }
        }
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Bcast(r + written, chunk, MPI_DOUBLE, 0, cr->nc.comm_intra);
            written += chunk;
            remain -= chunk;
        }
    }
    else
    {
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mygroup);
            written += chunk;
            remain -= chunk;
        }
    }
#endif
}

void gmx_sumf(std::size_t gmx_unused nr, float gmx_unused r[], const t_commrec gmx_unused* cr)
{
    // Without MPI we have a single rank, so sum is a no-op
#if GMX_MPI
    if (cr->sizeOfMyGroupCommunicator == 1)
    {
        return;
    }

    constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();
    if (cr->nc.bUse)
    {
        /* Use two step summing.  */
        if (cr->nc.rank_intra == 0)
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                MPI_Reduce(MPI_IN_PLACE, r + written, chunk, MPI_FLOAT, MPI_SUM, 0, cr->nc.comm_intra);
                /* Sum the roots of the internal (intra) buffers */
                MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_FLOAT, MPI_SUM, cr->nc.comm_inter);
                written += chunk;
                remain -= chunk;
            }
        }
        else
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
                MPI_Reduce(r + written, nullptr, chunk, MPI_FLOAT, MPI_SUM, 0, cr->nc.comm_intra);
                written += chunk;
                remain -= chunk;
            }
        }
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Bcast(r + written, chunk, MPI_FLOAT, 0, cr->nc.comm_intra);
            written += chunk;
            remain -= chunk;
        }
    }
    else
    {
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_FLOAT, MPI_SUM, cr->mpi_comm_mygroup);
            written += chunk;
            remain -= chunk;
        }
    }
#endif
}

void gmx_sumi(std::size_t gmx_unused nr, int gmx_unused r[], const t_commrec gmx_unused* cr)
{
    // Without MPI we have a single rank, so sum is a no-op
#if GMX_MPI
    if (cr->sizeOfMyGroupCommunicator == 1)
    {
        return;
    }

    constexpr std::size_t maxSignedInt = std::numeric_limits<int>::max();
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0)
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                MPI_Reduce(MPI_IN_PLACE, r + written, chunk, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
                /* Sum with the buffers reversed */
                MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_INT, MPI_SUM, cr->nc.comm_inter);
                written += chunk;
                remain -= chunk;
            }
        }
        else
        {
            for (std::size_t written = 0, remain = nr; remain > 0;)
            {
                std::size_t chunk = std::min(remain, maxSignedInt);
                /* This is here because of the silly MPI specification
                   that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
                MPI_Reduce(r + written, nullptr, chunk, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
                written += chunk;
                remain -= chunk;
            }
        }
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Bcast(r + written, chunk, MPI_INT, 0, cr->nc.comm_intra);
            written += chunk;
            remain -= chunk;
        }
    }
    else
    {
        for (std::size_t written = 0, remain = nr; remain > 0;)
        {
            std::size_t chunk = std::min(remain, maxSignedInt);
            MPI_Allreduce(MPI_IN_PLACE, r + written, chunk, MPI_INT, MPI_SUM, cr->mpi_comm_mygroup);
            written += chunk;
            remain -= chunk;
        }
    }
#endif
}

const char* opt2fn_main(const char* opt, int nfile, const t_filenm fnm[], t_commrec* cr)
{
    return SIMMAIN(cr) ? opt2fn(opt, nfile, fnm) : nullptr;
}

void gmx_fatal_collective(int                    f_errno,
                          const char*            file,
                          int                    line,
                          MPI_Comm               comm,
                          gmx_bool               bMain,
                          gmx_fmtstr const char* fmt,
                          ...)
{
    std::va_list ap;
    gmx_bool     bFinalize;
#if GMX_MPI
    int result;
    /* Check if we are calling on all processes in MPI_COMM_WORLD */
    MPI_Comm_compare(comm, MPI_COMM_WORLD, &result);
    /* Any result except MPI_UNEQUAL allows us to call MPI_Finalize */
    bFinalize = (result != MPI_UNEQUAL);
#else
    GMX_UNUSED_VALUE(comm);
    bFinalize = TRUE;
#endif

    va_start(ap, fmt);
    gmx_fatal_mpi_va(f_errno, file, line, bMain, bFinalize, fmt, ap);
    va_end(ap);
}
