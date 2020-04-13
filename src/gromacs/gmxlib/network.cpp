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
#include "gmxpre.h"

#include "network.h"

#include "config.h"

#include <cctype>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/mpiinplacebuffers.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

CommrecHandle init_commrec(MPI_Comm communicator, const gmx_multisim_t* ms)
{
    CommrecHandle handle;
    t_commrec*    cr;

    snew(cr, 1);
    handle.reset(cr);

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

    if (ms != nullptr)
    {
#if GMX_MPI
        cr->nnodes = sizeOfCommunicator / ms->nsim;
        MPI_Comm_split(communicator, ms->sim, rankInCommunicator, &cr->mpi_comm_mysim);
        cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
        MPI_Comm_rank(cr->mpi_comm_mysim, &cr->sim_nodeid);
        MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
#endif
    }
    else
    {
        cr->nnodes           = sizeOfCommunicator;
        cr->nodeid           = rankInCommunicator;
        cr->sim_nodeid       = cr->nodeid;
        cr->mpi_comm_mysim   = communicator;
        cr->mpi_comm_mygroup = communicator;
    }

    // TODO cr->duty should not be initialized here
    cr->duty = (DUTY_PP | DUTY_PME);

#if GMX_MPI && !MPI_IN_PLACE_EXISTS
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(cr->mpb, 1);
    cr->mpb->ibuf        = nullptr;
    cr->mpb->libuf       = nullptr;
    cr->mpb->fbuf        = nullptr;
    cr->mpb->dbuf        = nullptr;
    cr->mpb->ibuf_alloc  = 0;
    cr->mpb->libuf_alloc = 0;
    cr->mpb->fbuf_alloc  = 0;
    cr->mpb->dbuf_alloc  = 0;
#endif

    return handle;
}

void done_commrec(t_commrec* cr)
{
    if (MASTER(cr))
    {
        if (nullptr != cr->dd)
        {
            // TODO: implement
            // done_domdec(cr->dd);
        }
        done_mpi_in_place_buf(cr->mpb);
    }
#if GMX_MPI
    // TODO We need to be able to free communicators, but the
    // structure of the commrec and domdec initialization code makes
    // it hard to avoid both leaks and double frees.
    bool mySimIsMyGroup = (cr->mpi_comm_mysim == cr->mpi_comm_mygroup);
    if (cr->mpi_comm_mysim != MPI_COMM_NULL && cr->mpi_comm_mysim != MPI_COMM_WORLD)
    {
        // TODO see above
        // MPI_Comm_free(&cr->mpi_comm_mysim);
    }
    if (!mySimIsMyGroup && cr->mpi_comm_mygroup != MPI_COMM_NULL && cr->mpi_comm_mygroup != MPI_COMM_WORLD)
    {
        // TODO see above
        // MPI_Comm_free(&cr->mpi_comm_mygroup);
    }
#endif
    sfree(cr);
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
            fprintf(fplog, "Using two step summing over %d groups of on average %.1f ranks\n\n", ng,
                    (real)n / (real)ng);
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
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_barrier");
#else
    MPI_Barrier(communicator);
#endif
}

void gmx_bcast(int gmx_unused nbytes, void gmx_unused* b, MPI_Comm gmx_unused communicator)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_bcast");
#else
    MPI_Bcast(b, nbytes, MPI_BYTE, 0, communicator);
#endif
}

void gmx_sumd(int gmx_unused nr, double gmx_unused r[], const t_commrec gmx_unused* cr)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumd");
#else
#    if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        if (cr->nc.rank_intra == 0)
        {
            /* Use two step summing. */
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, 0, cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers. */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_DOUBLE, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_DOUBLE, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mygroup);
    }
#    else
    int i;

    if (nr > cr->mpb->dbuf_alloc)
    {
        cr->mpb->dbuf_alloc = nr;
        srenew(cr->mpb->dbuf, cr->mpb->dbuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->dbuf, nr, MPI_DOUBLE, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->dbuf, r, nr, MPI_DOUBLE, MPI_SUM, cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_DOUBLE, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->dbuf, nr, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->dbuf[i];
        }
    }
#    endif
#endif
}

void gmx_sumf(int gmx_unused nr, float gmx_unused r[], const t_commrec gmx_unused* cr)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumf");
#else
#    if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing.  */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, 0, cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_FLOAT, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_FLOAT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, cr->mpi_comm_mygroup);
    }
#    else
    int i;

    if (nr > cr->mpb->fbuf_alloc)
    {
        cr->mpb->fbuf_alloc = nr;
        srenew(cr->mpb->fbuf, cr->mpb->fbuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->fbuf, nr, MPI_FLOAT, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->fbuf, r, nr, MPI_FLOAT, MPI_SUM, cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_FLOAT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->fbuf, nr, MPI_FLOAT, MPI_SUM, cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->fbuf[i];
        }
    }
#    endif
#endif
}

void gmx_sumi(int gmx_unused nr, int gmx_unused r[], const t_commrec gmx_unused* cr)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumi");
#else
#    if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_INT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, cr->mpi_comm_mygroup);
    }
#    else
    int i;

    if (nr > cr->mpb->ibuf_alloc)
    {
        cr->mpb->ibuf_alloc = nr;
        srenew(cr->mpb->ibuf, cr->mpb->ibuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->ibuf, nr, MPI_INT, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->ibuf, r, nr, MPI_INT, MPI_SUM, cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_INT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->ibuf, nr, MPI_INT, MPI_SUM, cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->ibuf[i];
        }
    }
#    endif
#endif
}

void gmx_sumli(int gmx_unused nr, int64_t gmx_unused r[], const t_commrec gmx_unused* cr)
{
#if !GMX_MPI
    GMX_RELEASE_ASSERT(false, "Invalid call to gmx_sumli");
#else
#    if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, 0, cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_INT64_T, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_INT64_T, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, cr->mpi_comm_mygroup);
    }
#    else
    int i;

    if (nr > cr->mpb->libuf_alloc)
    {
        cr->mpb->libuf_alloc = nr;
        srenew(cr->mpb->libuf, cr->mpb->libuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->libuf, nr, MPI_INT64_T, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->libuf, r, nr, MPI_INT64_T, MPI_SUM, cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_INT64_T, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->libuf, nr, MPI_INT64_T, MPI_SUM, cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->libuf[i];
        }
    }
#    endif
#endif
}

const char* opt2fn_master(const char* opt, int nfile, const t_filenm fnm[], t_commrec* cr)
{
    return SIMMASTER(cr) ? opt2fn(opt, nfile, fnm) : nullptr;
}

void gmx_fatal_collective(int                    f_errno,
                          const char*            file,
                          int                    line,
                          MPI_Comm               comm,
                          gmx_bool               bMaster,
                          gmx_fmtstr const char* fmt,
                          ...)
{
    va_list  ap;
    gmx_bool bFinalize;
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
    gmx_fatal_mpi_va(f_errno, file, line, bMaster, bFinalize, fmt, ap);
    va_end(ap);
}
