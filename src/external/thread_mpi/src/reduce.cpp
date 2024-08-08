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


int tMPI_Reduce_run_op(void *dest, const void *src_a, const void *src_b,
                       tMPI_Datatype datatype, int count, tMPI_Op op,
                       tMPI_Comm comm)
{
    tMPI_Op_fn fn = datatype->op_functions[op];

    if (src_a == src_b)
    {
        return tMPI_Error(comm, TMPI_ERR_XFER_BUF_OVERLAP);
    }
    fn(dest, src_a, src_b, count);
    return TMPI_SUCCESS;
}

int tMPI_Reduce_fast(const void* sendbuf, void* recvbuf, int count,
                     tMPI_Datatype datatype, tMPI_Op op, int root,
                     tMPI_Comm comm)
{
    struct tmpi_thread *cur    = tMPI_Get_current();
    int                 myrank = tMPI_Comm_seek_rank(comm, cur);

    /* this function uses a binary tree-like reduction algorithm: */
    int N          = tMPI_Comm_N(comm);
    int myrank_rtr = (N+myrank-root)%N; /* my rank relative to root */
    int Nred       = N;                 /* number of neighbours that still communicate
                                           (decreases exponentially) */
    int nbr_dist   = 1;                 /* distance between communicating neighbours
                                           (increases exponentially) */
    int stepping   = 2;                 /* distance between non-communicating neighbours
                                           (increases exponentially) */
    int iteration  = 0;

    if (count == 0)
    {
        return TMPI_SUCCESS;
    }
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!recvbuf)
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }
    if ( (!datatype->op_functions) || (!datatype->op_functions[op]) )
    {
        return tMPI_Error(comm, TMPI_ERR_OP_FN);
    }

    if (sendbuf == TMPI_IN_PLACE) /* i.e. sendbuf == tMPI_IN_PLACE */
    {
        sendbuf = recvbuf;
    }
    /* we set our send and recv buffer s*/
    tMPI_Atomic_ptr_set(&(comm->reduce_sendbuf[myrank]), (void*)(sendbuf));
    tMPI_Atomic_ptr_set(&(comm->reduce_recvbuf[myrank]), recvbuf);

    while (Nred > 1)
    {
        /* calculate neighbour rank (here we use the real rank) */
        int nbr = (myrank_rtr%stepping == 0) ?
            (N+myrank+nbr_dist)%N :
            (N+myrank-nbr_dist)%N;

#ifdef TMPI_DEBUG
        printf("%d: iteration %d: myrank_rtr=%d, stepping=%d\n",
               myrank, iteration, myrank_rtr, stepping);
        fflush(stdout);
#endif
        /* check if I'm the reducing thread in this iteration's pair: */
        if (myrank_rtr%stepping == 0)
        {
            const void *a, *b;
            int   ret;

            /* now wait for my neighbor's data to become ready.
               First check if I actually have a neighbor. */
            if (myrank_rtr+nbr_dist < N)
            {
#ifdef TMPI_DEBUG
                printf("%d: waiting to reduce with %d, iteration=%d\n",
                       myrank, nbr, iteration);
                fflush(stdout);
#endif

#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
                tMPI_Profile_wait_start(cur);
#endif
                tMPI_Event_wait( &(comm->csync[myrank].events[nbr]) );
                tMPI_Event_process( &(comm->csync[myrank].events[nbr]), 1);
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
                tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);
#endif

#ifdef TMPI_DEBUG
                printf("%d: reducing with %d, iteration=%d\n",
                       myrank, nbr, iteration);
                fflush(stdout);
#endif
                /* we reduce with our neighbour*/
                if (iteration == 0)
                {
                    /* for the first iteration, the inputs are in the
                       sendbuf*/
                    a = sendbuf;
                    b = (void*)tMPI_Atomic_ptr_get(&(comm->reduce_sendbuf[nbr]));
                }
                else
                {
                    /* after the first operation, they're already in
                       the recvbuf */
                    a = recvbuf;
                    b = (void*)tMPI_Atomic_ptr_get(&(comm->reduce_recvbuf[nbr]));
                }

                if ((ret = tMPI_Reduce_run_op(recvbuf, a, b, datatype,
                                              count, op, comm)) != TMPI_SUCCESS)
                {
                    return ret;
                }

                /* signal to my neighbour that I'm ready. */
                tMPI_Event_signal( &(comm->csync[nbr].events[myrank]) );
            }
            else
            {
#ifdef TMPI_DEBUG
                printf("%d: not waiting copying buffer\n", myrank);
                fflush(stdout);
#endif
                /* we still need to put things in the right buffer for the next
                   iteration. We need to check for overlapping buffers
                   here because MPI_IN_PLACE might cause recvbuf to be the
                   same as sendbuf. */
                if (iteration == 0 && (recvbuf != sendbuf))
                {
                    memcpy(recvbuf, sendbuf, datatype->size*count);
                }
            }

        }
        else
        {
            /* the other thread is doing the reducing; we can just
               wait and break when ready */
            /* Awake our neighbour */
            tMPI_Event_signal( &(comm->csync[nbr].events[myrank]) );


#ifdef TMPI_DEBUG
            printf("%d: signalled %d, now waiting: iteration=%d\n",
                   nbr, myrank,  iteration);
            fflush(stdout);
#endif

            /* And wait for an incoming event from out neighbour */
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
            tMPI_Profile_wait_start(cur);
#endif
            tMPI_Event_wait( &(comm->csync[myrank].events[nbr]) );
            tMPI_Event_process( &(comm->csync[myrank].events[nbr]), 1);
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
            tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);
#endif
            /* now we can break because our data is reduced, and
               our neighbour goes on reducing it further. */
            break;
        }

#ifdef TMPI_DEBUG
        printf("%d: iteration over, iteration=%d\n", myrank,  iteration);
        fflush(stdout);
#endif

        Nred      = Nred/2 + Nred%2;
        nbr_dist *= 2;
        stepping *= 2;
        iteration++;
    }

    return TMPI_SUCCESS;
}

int tMPI_Reduce(const void* sendbuf, void* recvbuf, int count,
                tMPI_Datatype datatype, tMPI_Op op, int root, tMPI_Comm comm)
{
    struct tmpi_thread *cur    = tMPI_Get_current();
    int                 myrank = tMPI_Comm_seek_rank(comm, cur);
    int                 ret;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Reduce(%p, %p, %d, %p, %p, %d, %p)",
                     sendbuf, recvbuf, count, datatype, op, root, comm);
#endif

    if (myrank == root)
    {
        if (sendbuf == TMPI_IN_PLACE) /* i.e. sendbuf == TMPI_IN_PLACE */
        {
            sendbuf = recvbuf;
        }
    }
    else
    {
#ifdef TMPI_WARN_MALLOC
        fprintf(stderr, "Warning: malloc during tMPI_Reduce\n");
#endif
        recvbuf = (void*)tMPI_Malloc(datatype->size*count);
    }
    ret = tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, root, comm);
    if (myrank != root)
    {
        free(recvbuf);
    }
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Reduce);
#endif
    return ret;
}

int tMPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
                   tMPI_Datatype datatype, tMPI_Op op, tMPI_Comm comm)
{
    void               *rootbuf = NULL; /* root process' receive buffer */
    struct tmpi_thread *cur     = tMPI_Get_current();
    int                 myrank  = tMPI_Comm_seek_rank(comm, cur);
    int                 ret;

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Allreduce(%p, %p, %d, %p, %p, %p)",
                     sendbuf, recvbuf, count, datatype, op, comm);
#endif
    if (count == 0)
    {
        return TMPI_SUCCESS;
    }
    if (!recvbuf)
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }
    if (sendbuf == TMPI_IN_PLACE) /* i.e. sendbuf == TMPI_IN_PLACE */
    {
        sendbuf = recvbuf;
    }

    ret = tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, 0, comm);
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_start(cur);
#endif
    tMPI_Barrier_wait( &(comm->barrier));
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);
#endif
    /* distribute rootbuf */
    rootbuf = (void*)tMPI_Atomic_ptr_get(&(comm->reduce_recvbuf[0]));

    /* and now we just copy things back. We know that the root thread
       arrives last, so there's no point in using tMPI_Scatter with
       copy buffers, etc. */
    if (myrank != 0)
    {
        if (rootbuf == recvbuf)
        {
            return tMPI_Error(comm, TMPI_ERR_XFER_BUF_OVERLAP);
        }
        memcpy(recvbuf, rootbuf, datatype->size*count );
    }

#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
    tMPI_Profile_wait_start(cur);
#endif
    tMPI_Barrier_wait( &(comm->barrier));
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);
    tMPI_Profile_count_stop(cur, TMPIFN_Allreduce);
#endif
    return ret;
}
