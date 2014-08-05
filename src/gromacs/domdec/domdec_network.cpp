/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008,2009,2010,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * \brief This file defines functions for (mostly) the domdec module
 * to use MPI functionality
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_network.h"

#include "config.h"

#include <string.h>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/utility/gmxmpi.h"


/*! \brief Returns the MPI rank of the domain decomposition master rank */
#define DDMASTERRANK(dd)   (dd->masterrank)


void dd_sendrecv_int(const gmx_domdec_t gmx_unused *dd,
                     int gmx_unused ddimind, int gmx_unused direction,
                     int gmx_unused *buf_s, int gmx_unused n_s,
                     int gmx_unused *buf_r, int gmx_unused n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(int), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(int), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv_real(const gmx_domdec_t gmx_unused *dd,
                      int gmx_unused ddimind, int gmx_unused direction,
                      real gmx_unused *buf_s, int gmx_unused n_s,
                      real gmx_unused *buf_r, int gmx_unused n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s, n_s*sizeof(real), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r, n_r*sizeof(real), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv_rvec(const gmx_domdec_t gmx_unused *dd,
                      int gmx_unused ddimind, int gmx_unused direction,
                      rvec gmx_unused *buf_s, int gmx_unused n_s,
                      rvec gmx_unused *buf_r, int gmx_unused n_r)
{
#ifdef GMX_MPI
    int        rank_s, rank_r;
    MPI_Status stat;

    rank_s = dd->neighbor[ddimind][direction == dddirForward ? 0 : 1];
    rank_r = dd->neighbor[ddimind][direction == dddirForward ? 1 : 0];

    if (n_s && n_r)
    {
        MPI_Sendrecv(buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }
    else if (n_s)
    {
        MPI_Send(    buf_s[0], n_s*sizeof(rvec), MPI_BYTE, rank_s, 0,
                     dd->mpi_comm_all);
    }
    else if (n_r)
    {
        MPI_Recv(    buf_r[0], n_r*sizeof(rvec), MPI_BYTE, rank_r, 0,
                     dd->mpi_comm_all, &stat);
    }

#endif
}

void dd_sendrecv2_rvec(const gmx_domdec_t gmx_unused *dd,
                       int gmx_unused ddimind,
                       rvec gmx_unused *buf_s_fw, int gmx_unused n_s_fw,
                       rvec gmx_unused *buf_r_fw, int gmx_unused n_r_fw,
                       rvec gmx_unused *buf_s_bw, int gmx_unused n_s_bw,
                       rvec gmx_unused *buf_r_bw, int gmx_unused n_r_bw)
{
#ifdef GMX_MPI
    int         rank_fw, rank_bw, nreq;
    MPI_Request req[4];
    MPI_Status  stat[4];

    rank_fw = dd->neighbor[ddimind][0];
    rank_bw = dd->neighbor[ddimind][1];

    if (!dd->bSendRecv2)
    {
        /* Try to send and receive in two directions simultaneously.
         * Should be faster, especially on machines
         * with full 3D communication networks.
         * However, it could be that communication libraries are
         * optimized for MPI_Sendrecv and non-blocking MPI calls
         * are slower.
         * SendRecv2 can be turned on with the env.var. GMX_DD_SENDRECV2
         */
        nreq = 0;
        if (n_r_fw)
        {
            MPI_Irecv(buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_r_bw)
        {
            MPI_Irecv(buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_fw)
        {
            MPI_Isend(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE,
                      rank_fw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_bw)
        {
            MPI_Isend(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE,
                      rank_bw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (nreq)
        {
            MPI_Waitall(nreq, req, stat);
        }
    }
    else
    {
        /* Communicate in two ordered phases.
         * This is slower, even on a dual-core Opteron cluster
         * with a single full-duplex network connection per machine.
         */
        /* Forward */
        MPI_Sendrecv(buf_s_fw[0], n_s_fw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     buf_r_fw[0], n_r_fw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     dd->mpi_comm_all, &stat[0]);
        /* Backward */
        MPI_Sendrecv(buf_s_bw[0], n_s_bw*sizeof(rvec), MPI_BYTE, rank_bw, 0,
                     buf_r_bw[0], n_r_bw*sizeof(rvec), MPI_BYTE, rank_fw, 0,
                     dd->mpi_comm_all, &stat[0]);
    }
#endif
}

/* IBM's BlueGene(/L) MPI_Bcast dereferences the data pointer
 * even when 0 == nbytes, so we protect calls to it on BlueGene.
 * Fortunately dd_bcast() and dd_bcastc() are only
 * called during DD setup and partition.
 */

void dd_bcast(gmx_domdec_t gmx_unused *dd, int gmx_unused nbytes, void gmx_unused *data)
{
#ifdef GMX_MPI
    if (dd->nnodes > 1)
    {
#ifdef GMX_BLUEGENE
        if (nbytes > 0)
        {
#endif
        MPI_Bcast(data, nbytes, MPI_BYTE,
                  DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
    }
#endif
    }
#endif
}

void dd_bcastc(gmx_domdec_t *dd, int nbytes, void *src, void *dest)
{
    if (DDMASTER(dd) || dd->nnodes == 1)
    {
        memcpy(dest, src, nbytes);
    }
#ifdef GMX_MPI
    if (dd->nnodes > 1)
    {
#ifdef GMX_BLUEGENE
        if (nbytes > 0)
        {
#endif
        MPI_Bcast(dest, nbytes, MPI_BYTE,
                  DDMASTERRANK(dd), dd->mpi_comm_all);
#ifdef GMX_BLUEGENE
    }
#endif
    }
#endif
}

void dd_scatter(gmx_domdec_t gmx_unused *dd, int gmx_unused nbytes, void gmx_unused *src, void *dest)
{
#ifdef GMX_MPI
    if (dd->nnodes > 1)
    {
        MPI_Scatter(src, nbytes, MPI_BYTE,
                    dest, nbytes, MPI_BYTE,
                    DDMASTERRANK(dd), dd->mpi_comm_all);
    }
    else
#endif
    {
        /* 1 rank, either we copy everything, or dest=src: nothing to do */
        if (dest != src)
        {
            memcpy(dest, src, nbytes);
        }
    }
}

void dd_gather(gmx_domdec_t gmx_unused *dd, int gmx_unused nbytes, void gmx_unused *src, void gmx_unused *dest)
{
#ifdef GMX_MPI
    MPI_Gather(src, nbytes, MPI_BYTE,
               dest, nbytes, MPI_BYTE,
               DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}

void dd_scatterv(gmx_domdec_t gmx_unused *dd,
                 int gmx_unused *scounts, int gmx_unused *disps, void *sbuf,
                 int rcount, void *rbuf)
{
#ifdef GMX_MPI
    int dum;

    if (dd->nnodes > 1)
    {
        if (rcount == 0)
        {
            /* MPI does not allow NULL pointers */
            rbuf = &dum;
        }
        MPI_Scatterv(sbuf, scounts, disps, MPI_BYTE,
                     rbuf, rcount, MPI_BYTE,
                     DDMASTERRANK(dd), dd->mpi_comm_all);
    }
    else
#endif
    {
        /* 1 rank, either we copy everything, or rbuf=sbuf: nothing to do */
        if (rbuf != sbuf)
        {
            memcpy(rbuf, sbuf, rcount);
        }
    }
}

void dd_gatherv(gmx_domdec_t gmx_unused *dd,
                int gmx_unused scount, void gmx_unused *sbuf,
                int gmx_unused *rcounts, int gmx_unused *disps, void gmx_unused *rbuf)
{
#ifdef GMX_MPI
    int dum;

    if (scount == 0)
    {
        /* MPI does not allow NULL pointers */
        sbuf = &dum;
    }
    MPI_Gatherv(sbuf, scount, MPI_BYTE,
                rbuf, rcounts, disps, MPI_BYTE,
                DDMASTERRANK(dd), dd->mpi_comm_all);
#endif
}
