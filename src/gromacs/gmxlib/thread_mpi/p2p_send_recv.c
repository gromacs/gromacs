/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1) Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2) Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3) Neither the name of the copyright holders nor the
   names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY US ''AS IS'' AND ANY
   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL WE BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   If you want to redistribute modifications, please consider that
   scientific software is very special. Version control is crucial -
   bugs must be traceable. We will be happy to consider code for
   inclusion in the official distribution, but derived work should not
   be called official thread_mpi. Details are found in the README & COPYING
   files.
 */

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


#include "impl.h"
#include "p2p.h"


/* point-to-point communication exported functions */

int tMPI_Send(void* buf, int count, tMPI_Datatype datatype, int dest,
              int tag, tMPI_Comm comm)
{
    struct envelope    *sev;
    struct tmpi_thread *send_dst;
    struct tmpi_thread *cur = tMPI_Get_current();
    struct tmpi_req_    req;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Send(%p, %d, %p, %d, %d, %p)", buf, count,
                     datatype, dest, tag, comm);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    send_dst = tMPI_Get_thread(comm, dest);
    if (!send_dst)
    {
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST);
    }

    sev = tMPI_Post_send(cur, comm, send_dst, buf, count, datatype, tag, FALSE);
    if (sev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(&req, sev);
    tMPI_Wait_single(cur, &req);

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Send);
#endif
    return req.error;
}




int tMPI_Recv(void* buf, int count, tMPI_Datatype datatype, int source,
              int tag, tMPI_Comm comm, tMPI_Status *status)
{
    struct envelope    *rev;
    struct tmpi_thread *recv_src = 0;
    struct tmpi_thread *cur      = tMPI_Get_current();
    struct tmpi_req_    req;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Recv(%p, %d, %p, %d, %d, %p, %p)", buf, count,
                     datatype, source, tag, comm, status);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (source != TMPI_ANY_SOURCE)
    {
        recv_src = tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC);
        }
    }

    rev = tMPI_Post_match_recv(cur, comm, recv_src, buf, count, datatype, tag,
                               FALSE);
    if (rev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(&req, rev);
    tMPI_Wait_single(cur, &req);

    tMPI_Set_status(&req, status);
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Recv);
#endif
    return req.error;
}





int tMPI_Sendrecv(void *sendbuf, int sendcount, tMPI_Datatype sendtype,
                  int dest, int sendtag, void *recvbuf, int recvcount,
                  tMPI_Datatype recvtype, int source, int recvtag,
                  tMPI_Comm comm, tMPI_Status *status)
{
    struct envelope    *rev, *sev;
    struct tmpi_thread *cur      = tMPI_Get_current();
    struct tmpi_thread *recv_src = 0;
    struct tmpi_thread *send_dst;
    struct tmpi_req_    sreq, rreq;
    int                 ret = TMPI_SUCCESS;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Sendrecv(%p, %d, %p, %d, %d, %p, %d, %p, %d, %d, %p, %p)",
                     sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                     recvcount, recvtype, source, recvtag, comm, status);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    send_dst = tMPI_Get_thread(comm, dest);
    if (!send_dst)
    {
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST);
    }
    if (source != TMPI_ANY_SOURCE)
    {
        recv_src = tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC);
        }
    }

    /* we first prepare to send */
    sev = tMPI_Post_send(cur, comm, send_dst, sendbuf, sendcount,
                         sendtype, sendtag, FALSE);
    if (sev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(&sreq, sev);
    /* the we prepare to receive */
    rev = tMPI_Post_match_recv(cur, comm, recv_src, recvbuf, recvcount,
                               recvtype, recvtag, FALSE);
    if (rev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(&rreq, rev);

    /* fix the pointers */
    sreq.next = &rreq;
    sreq.prev = NULL;
    rreq.prev = &sreq;
    rreq.next = NULL;

    /* and wait for our requests */
    do
    {
        if (tMPI_Test_multi(cur, &sreq, NULL))
        {
            break;
        }
        tMPI_Wait_process_incoming(cur);
    }
    while (TRUE);

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Sendrecv);
#endif

    tMPI_Set_status(&rreq, status);
    ret = sreq.error;
    if (rreq.error != TMPI_SUCCESS)
    {
        ret = rreq.error;
    }

    return ret;
}





/* async */

int tMPI_Isend(void* buf, int count, tMPI_Datatype datatype, int dest,
               int tag, tMPI_Comm comm, tMPI_Request *request)
{
    struct tmpi_thread *cur = tMPI_Get_current();
    struct req_list    *rql = &(cur->rql);
    struct tmpi_req_   *rq  = tMPI_Get_req(rql);
    struct tmpi_thread *send_dst;
    struct envelope    *ev;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Isend(%p, %d, %p, %d, %d, %p, %p)", buf, count,
                     datatype, dest, tag, comm, request);
#endif
    if (!comm)
    {
        tMPI_Return_req(rql, rq);
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    send_dst = tMPI_Get_thread(comm, dest);
    if (!send_dst)
    {
        tMPI_Return_req(rql, rq);
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST);
    }
    ev = tMPI_Post_send(cur, comm, send_dst, buf, count, datatype, tag, TRUE);
    if (ev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(rq, ev);
    *request = rq;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Isend);
#endif
    return ev->error;
}


int tMPI_Irecv(void* buf, int count, tMPI_Datatype datatype, int source,
               int tag, tMPI_Comm comm, tMPI_Request *request)
{
    struct tmpi_thread *cur      = tMPI_Get_current();
    struct req_list    *rql      = &(cur->rql);
    struct tmpi_req_   *rq       = tMPI_Get_req(rql);
    struct tmpi_thread *recv_src = 0;
    struct envelope    *ev;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Irecv(%p, %d, %p, %d, %d, %p, %p)", buf, count,
                     datatype, source, tag, comm, request);
#endif
    if (!comm)
    {
        tMPI_Return_req(rql, rq);
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (source != TMPI_ANY_SOURCE)
    {
        recv_src = tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            tMPI_Return_req(rql, rq);
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC);
        }
    }
    ev = tMPI_Post_match_recv(cur, comm, recv_src, buf, count, datatype, tag,
                              TRUE);
    if (ev == NULL)
    {
        return TMPI_ERR_ENVELOPES;
    }
    tMPI_Req_init(rq, ev);
    *request = rq;
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Irecv);
#endif
    return ev->error;
}
