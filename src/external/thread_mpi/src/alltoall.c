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
#include "collective.h"



int tMPI_Alltoall(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                  void* recvbuf, int recvcount, tMPI_Datatype recvtype,
                  tMPI_Comm comm)
{
    int                 synct;
    struct coll_env    *cev;
    int                 myrank;
    int                 ret = TMPI_SUCCESS;
    int                 i;
    size_t              sendsize = sendtype->size*sendcount;
    size_t              recvsize = recvtype->size*recvcount;
    int                 n_remaining;
    struct tmpi_thread *cur = tMPI_Get_current();

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Alltoall(%p, %d, %p, %p, %d, %p, %p)",
                     sendbuf, sendcount, sendtype,
                     recvbuf, recvcount, recvtype, comm);
#endif

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }

    myrank = tMPI_Comm_seek_rank(comm, cur);

    /* we increase our counter, and determine which coll_env we get */
    cev = tMPI_Get_cev(comm, myrank, &synct);

    /* post our pointers */
    /* we set up multiple posts, so no Post_multi */
    cev->met[myrank].tag      = TMPI_ALLTOALL_TAG;
    cev->met[myrank].datatype = sendtype;
    tMPI_Atomic_set( &(cev->met[myrank].n_remaining), cev->N-1 );
    for (i = 0; i < comm->grp.N; i++)
    {
        cev->met[myrank].bufsize[i]   = sendsize;
        cev->met[myrank].buf[i]       = (char*)sendbuf+sendsize*i;
        cev->met[myrank].read_data[i] = FALSE;
    }
    tMPI_Atomic_memory_barrier_rel();
    tMPI_Atomic_set(&(cev->met[myrank].current_sync), synct);

    /* post availability */
    for (i = 0; i < cev->N; i++)
    {
        if (i != myrank)
        {
            tMPI_Event_signal( &(cev->met[i].recv_ev) );
        }
    }

    /* we don't do the copy buffer thing here because it's pointless:
       the processes have to synchronize anyway, because they all
       send and receive. */

    /* do root transfer */
    tMPI_Coll_root_xfer(comm, sendtype, recvtype,
                        sendsize, recvsize,
                        (char*)sendbuf+sendsize*myrank,
                        (char*)recvbuf+recvsize*myrank, &ret);
    cev->met[myrank].read_data[myrank] = TRUE;
    /* and poll data availability */
    n_remaining = cev->N-1;
    while (n_remaining > 0)
    {
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_start(cur);
#endif
        tMPI_Event_wait( &(cev->met[myrank]).recv_ev );
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_stop(cur, TMPIWAIT_Coll_recv);
#endif
        for (i = 0; i < cev->N; i++)
        {
            if ((!cev->met[myrank].read_data[i]) &&
                (tMPI_Atomic_get(&(cev->met[i].current_sync)) == synct))
            {
                tMPI_Event_process( &(cev->met[myrank]).recv_ev, 1);
                tMPI_Mult_recv(comm, cev, i, myrank, TMPI_ALLTOALL_TAG,
                               recvtype, recvsize, (char*)recvbuf+recvsize*i,
                               &ret);
                if (ret != TMPI_SUCCESS)
                {
                    return ret;
                }
                cev->met[myrank].read_data[i] = TRUE;
                n_remaining--;
            }
        }
    }


    /* and wait until everybody is done copying our data */
    tMPI_Wait_for_others(cev, myrank);

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Alltoall);
#endif
    return ret;
}


int tMPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls,
                   tMPI_Datatype sendtype,
                   void* recvbuf, int *recvcounts, int *rdispls,
                   tMPI_Datatype recvtype,
                   tMPI_Comm comm)

{
    int                 synct;
    struct coll_env    *cev;
    int                 myrank;
    int                 ret = TMPI_SUCCESS;
    int                 i;
    int                 n_remaining;
    struct tmpi_thread *cur = tMPI_Get_current();

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Alltoallv(%p, %p, %p, %p, %p, %p, %p, %p, %p, %p)",
                     sendbuf, sendcounts, sdispls, sendtype,
                     recvbuf, recvcounts, rdispls, recvtype,
                     comm);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }

    myrank = tMPI_Comm_seek_rank(comm, cur);

    /* we increase our counter, and determine which coll_env we get */
    cev = tMPI_Get_cev(comm, myrank, &synct);

    /* post our pointers */
    /* we set up multiple posts, so no Post_multi */
    cev->met[myrank].tag      = TMPI_ALLTOALLV_TAG;
    cev->met[myrank].datatype = sendtype;
    tMPI_Atomic_set( &(cev->met[myrank].n_remaining), cev->N-1 );
    for (i = 0; i < comm->grp.N; i++)
    {
        cev->met[myrank].bufsize[i]   = sendtype->size*sendcounts[i];
        cev->met[myrank].buf[i]       = (char*)sendbuf+sendtype->size*sdispls[i];
        cev->met[myrank].read_data[i] = FALSE;
    }
    tMPI_Atomic_memory_barrier_rel();
    tMPI_Atomic_set(&(cev->met[myrank].current_sync), synct);

    /* post availability */
    for (i = 0; i < cev->N; i++)
    {
        if (i != myrank)
        {
            tMPI_Event_signal( &(cev->met[i].recv_ev) );
        }
    }

    /* we don't do the copy buffer thing here because it's pointless:
       the processes have to synchronize anyway, because they all
       send and receive. */

    /* do root transfer */
    tMPI_Coll_root_xfer(comm, sendtype, recvtype,
                        sendtype->size*sendcounts[myrank],
                        recvtype->size*recvcounts[myrank],
                        (char*)sendbuf+sendtype->size*sdispls[myrank],
                        (char*)recvbuf+recvtype->size*rdispls[myrank], &ret);
    cev->met[myrank].read_data[myrank] = TRUE;

    /* and poll data availability */
    n_remaining = cev->N-1;
    while (n_remaining > 0)
    {
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_start(cur);
#endif
        tMPI_Event_wait( &(cev->met[myrank]).recv_ev );
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_stop(cur, TMPIWAIT_Coll_recv);
#endif
        for (i = 0; i < cev->N; i++)
        {
            if ((!cev->met[myrank].read_data[i]) &&
                (tMPI_Atomic_get(&(cev->met[i].current_sync)) == synct) )
            {
                tMPI_Event_process( &(cev->met[myrank]).recv_ev, 1);
                tMPI_Mult_recv(comm, cev, i, myrank, TMPI_ALLTOALLV_TAG,
                               recvtype, recvtype->size*recvcounts[i],
                               (char*)recvbuf+recvtype->size*rdispls[i], &ret);
                if (ret != TMPI_SUCCESS)
                {
                    return ret;
                }
                cev->met[myrank].read_data[i] = TRUE;
                n_remaining--;
            }
        }
    }

    /* and wait until everybody is done copying our data */
    tMPI_Wait_for_others(cev, myrank);

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Alltoallv);
#endif
    return ret;
}
