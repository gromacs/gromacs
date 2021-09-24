/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2008-2019,2020,2021, by the GROMACS development team, led by
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

#include <cstring>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/utility/gmxmpi.h"

#include "domdec_internal.h"


/*! \brief Returns the MPI rank of the domain decomposition master rank */
#define DDMASTERRANK(dd) ((dd)->masterrank)


/*! \brief Move data of type \p T in the communication region one cell along
 * the domain decomposition
 *
 * Moves in the dimension indexed by ddDimensionIndex, either forward
 * (direction=dddirFoward) or backward (direction=dddirBackward).
 */
template<typename T>
void ddSendrecv(const struct gmx_domdec_t* dd,
                int                        ddDimensionIndex,
                int                        direction,
                T*                         sendBuffer,
                int                        numElementsToSend,
                T*                         receiveBuffer,
                int                        numElementsToReceive)
{
#if GMX_MPI
    int sendRank    = dd->neighbor[ddDimensionIndex][direction == dddirForward ? 0 : 1];
    int receiveRank = dd->neighbor[ddDimensionIndex][direction == dddirForward ? 1 : 0];

    constexpr int mpiTag = 0;
    MPI_Status    mpiStatus;
    if (numElementsToSend > 0 && numElementsToReceive > 0)
    {
        MPI_Sendrecv(sendBuffer,
                     numElementsToSend * sizeof(T),
                     MPI_BYTE,
                     sendRank,
                     mpiTag,
                     receiveBuffer,
                     numElementsToReceive * sizeof(T),
                     MPI_BYTE,
                     receiveRank,
                     mpiTag,
                     dd->mpi_comm_all,
                     &mpiStatus);
    }
    else if (numElementsToSend > 0)
    {
        MPI_Send(sendBuffer, numElementsToSend * sizeof(T), MPI_BYTE, sendRank, mpiTag, dd->mpi_comm_all);
    }
    else if (numElementsToReceive > 0)
    {
        MPI_Recv(receiveBuffer, numElementsToReceive * sizeof(T), MPI_BYTE, receiveRank, mpiTag, dd->mpi_comm_all, &mpiStatus);
    }
#else  // GMX_MPI
    GMX_UNUSED_VALUE(dd);
    GMX_UNUSED_VALUE(ddDimensionIndex);
    GMX_UNUSED_VALUE(direction);
    GMX_UNUSED_VALUE(sendBuffer);
    GMX_UNUSED_VALUE(numElementsToSend);
    GMX_UNUSED_VALUE(receiveBuffer);
    GMX_UNUSED_VALUE(numElementsToReceive);
#endif // GMX_MPI
}

//! Specialization of extern template for int
template void ddSendrecv(const gmx_domdec_t*, int, int, int*, int, int*, int);
//! Specialization of extern template for real
template void ddSendrecv(const gmx_domdec_t*, int, int, real*, int, real*, int);
//! Specialization of extern template for gmx::RVec
template void ddSendrecv(const gmx_domdec_t*, int, int, rvec*, int, rvec*, int);

template<typename T>
void ddSendrecv(const gmx_domdec_t* dd,
                int                 ddDimensionIndex,
                int                 direction,
                gmx::ArrayRef<T>    sendBuffer,
                gmx::ArrayRef<T>    receiveBuffer)
{
    ddSendrecv(dd,
               ddDimensionIndex,
               direction,
               sendBuffer.data(),
               sendBuffer.size(),
               receiveBuffer.data(),
               receiveBuffer.size());
}

//! Specialization of extern template for int
template void ddSendrecv(const gmx_domdec_t*, int, int, gmx::ArrayRef<int>, gmx::ArrayRef<int>);
//! Specialization of extern template for real
template void ddSendrecv(const gmx_domdec_t*, int, int, gmx::ArrayRef<real>, gmx::ArrayRef<real>);
//! Specialization of extern template for gmx::RVec
template void ddSendrecv(const gmx_domdec_t*, int, int, gmx::ArrayRef<gmx::RVec>, gmx::ArrayRef<gmx::RVec>);

void dd_sendrecv2_rvec(const struct gmx_domdec_t gmx_unused* dd,
                       int gmx_unused                        ddimind,
                       rvec gmx_unused* buf_s_fw,
                       int gmx_unused   n_s_fw,
                       rvec gmx_unused* buf_r_fw,
                       int gmx_unused   n_r_fw,
                       rvec gmx_unused* buf_s_bw,
                       int gmx_unused   n_s_bw,
                       rvec gmx_unused* buf_r_bw,
                       int gmx_unused   n_r_bw)
{
#if GMX_MPI
    MPI_Request req[4];
    MPI_Status  stat[4];

    int rank_fw = dd->neighbor[ddimind][0];
    int rank_bw = dd->neighbor[ddimind][1];

    if (!dd->comm->ddSettings.useSendRecv2)
    {
        /* Try to send and receive in two directions simultaneously.
         * Should be faster, especially on machines
         * with full 3D communication networks.
         * However, it could be that communication libraries are
         * optimized for MPI_Sendrecv and non-blocking MPI calls
         * are slower.
         * SendRecv2 can be turned on with the env.var. GMX_DD_SENDRECV2
         */
        int nreq = 0;
        if (n_r_fw)
        {
            MPI_Irecv(buf_r_fw[0], n_r_fw * sizeof(rvec), MPI_BYTE, rank_bw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_r_bw)
        {
            MPI_Irecv(buf_r_bw[0], n_r_bw * sizeof(rvec), MPI_BYTE, rank_fw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_fw)
        {
            MPI_Isend(buf_s_fw[0], n_s_fw * sizeof(rvec), MPI_BYTE, rank_fw, 0, dd->mpi_comm_all, &req[nreq++]);
        }
        if (n_s_bw)
        {
            MPI_Isend(buf_s_bw[0], n_s_bw * sizeof(rvec), MPI_BYTE, rank_bw, 1, dd->mpi_comm_all, &req[nreq++]);
        }
        if (nreq)
        {
            MPI_Waitall(nreq, req, stat); //NOLINT(clang-analyzer-optin.mpi.MPI-Checker)
        }
    }
    else
    {
        /* Communicate in two ordered phases.
         * This is slower, even on a dual-core Opteron cluster
         * with a single full-duplex network connection per machine.
         */
        /* Forward */
        MPI_Sendrecv(buf_s_fw[0],
                     n_s_fw * sizeof(rvec),
                     MPI_BYTE,
                     rank_fw,
                     0,
                     buf_r_fw[0],
                     n_r_fw * sizeof(rvec),
                     MPI_BYTE,
                     rank_bw,
                     0,
                     dd->mpi_comm_all,
                     &stat[0]);
        /* Backward */
        MPI_Sendrecv(buf_s_bw[0],
                     n_s_bw * sizeof(rvec),
                     MPI_BYTE,
                     rank_bw,
                     0,
                     buf_r_bw[0],
                     n_r_bw * sizeof(rvec),
                     MPI_BYTE,
                     rank_fw,
                     0,
                     dd->mpi_comm_all,
                     &stat[0]);
    }
#endif
}

void dd_bcast(const gmx_domdec_t gmx_unused* dd, int gmx_unused nbytes, void gmx_unused* data)
{
#if GMX_MPI
    if (dd->nnodes > 1)
    {
        MPI_Bcast(data, nbytes, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
    }
#endif
}

void dd_bcastc(const gmx_domdec_t* dd, int nbytes, void* src, void* dest)
{
    if (DDMASTER(dd) || dd->nnodes == 1)
    {
        memcpy(dest, src, nbytes);
    }
#if GMX_MPI
    if (dd->nnodes > 1)
    {
        MPI_Bcast(dest, nbytes, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
    }
#endif
}

void dd_scatter(const gmx_domdec_t gmx_unused* dd, int gmx_unused nbytes, const void gmx_unused* src, void* dest)
{
#if GMX_MPI
    if (dd->nnodes > 1)
    {
        /* Some MPI implementions don't specify const */
        MPI_Scatter(const_cast<void*>(src), nbytes, MPI_BYTE, dest, nbytes, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
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

void dd_gather(const gmx_domdec_t gmx_unused* dd,
               int gmx_unused                 nbytes,
               const void gmx_unused* src,
               void gmx_unused* dest)
{
#if GMX_MPI
    if (dd->nnodes > 1)
    {
        /* Some MPI implementions don't specify const */
        MPI_Gather(const_cast<void*>(src), nbytes, MPI_BYTE, dest, nbytes, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
    }
    else
#endif
    {
        memcpy(dest, src, nbytes);
    }
}

void dd_scatterv(const gmx_domdec_t gmx_unused* dd,
                 int gmx_unused* scounts,
                 int gmx_unused* disps,
                 const void*     sbuf,
                 int             rcount,
                 void*           rbuf)
{
#if GMX_MPI
    int dum = 0;

    if (dd->nnodes > 1)
    {
        if (rcount == 0)
        {
            /* MPI does not allow NULL pointers */
            rbuf = &dum;
        }
        /* Some MPI implementions don't specify const */
        MPI_Scatterv(
                const_cast<void*>(sbuf), scounts, disps, MPI_BYTE, rbuf, rcount, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
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

void dd_gatherv(const gmx_domdec_t gmx_unused* dd,
                int gmx_unused                 scount,
                const void gmx_unused* sbuf,
                int gmx_unused* rcounts,
                int gmx_unused* disps,
                void gmx_unused* rbuf)
{
#if GMX_MPI
    if (dd->nnodes > 1)
    {
        int dum;

        if (scount == 0)
        {
            /* MPI does not allow NULL pointers */
            sbuf = &dum;
        }
        /* Some MPI implementions don't specify const */
        MPI_Gatherv(
                const_cast<void*>(sbuf), scount, MPI_BYTE, rbuf, rcounts, disps, MPI_BYTE, DDMASTERRANK(dd), dd->mpi_comm_all);
    }
    else
#endif
    {
        memcpy(rbuf, sbuf, rcounts[0]);
    }
}
