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



int tMPI_Scan(void* sendbuf, void* recvbuf, int count,
              tMPI_Datatype datatype, tMPI_Op op, tMPI_Comm comm)
{
    struct tmpi_thread *cur    = tMPI_Get_current();
    int                 myrank = tMPI_Comm_seek_rank(comm, cur);
    int                 N      = tMPI_Comm_N(comm);
    int                 prev   = myrank - 1; /* my previous neighbor */
    int                 next   = myrank + 1; /* my next neighbor */

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Scan(%p, %p, %d, %p, %p, %p)",
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
    if (sendbuf == TMPI_IN_PLACE)
    {
        sendbuf = recvbuf;
    }

    /* we set our send and recv buffers */
    tMPI_Atomic_ptr_set(&(comm->reduce_sendbuf[myrank]), sendbuf);
    tMPI_Atomic_ptr_set(&(comm->reduce_recvbuf[myrank]), recvbuf);

    /* now wait for the previous rank to finish */
    if (myrank > 0)
    {
        void *a, *b;
        int   ret;

#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_start(cur);
#endif
        /* wait for the previous neighbor's data to be ready */
        tMPI_Event_wait( &(comm->csync[myrank].events[prev]) );
        tMPI_Event_process( &(comm->csync[myrank].events[prev]), 1);
#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
        tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);
#endif
#ifdef TMPI_DEBUG
        printf("%d: scanning with %d \n", myrank, prev, iteration);
        fflush(stdout);
#endif
        /* now do the reduction */
        if (prev > 0)
        {
            a = (void*)tMPI_Atomic_ptr_get(&(comm->reduce_recvbuf[prev]));
        }
        else
        {
            a = (void*)tMPI_Atomic_ptr_get(&(comm->reduce_sendbuf[prev]));
        }
        b = sendbuf;

        if ((ret = tMPI_Reduce_run_op(recvbuf, a, b, datatype,
                                      count, op, comm)) != TMPI_SUCCESS)
        {
            return ret;
        }

        /* signal to my previous neighbor that I'm done with the data */
        tMPI_Event_signal( &(comm->csync[prev].events[prev]) );
    }
    else
    {
        if (sendbuf != recvbuf)
        {
            /* copy the data if this is rank 0, and not MPI_IN_PLACE */
            memcpy(recvbuf, sendbuf, count*datatype->size);
        }
    }

    if (myrank < N-1)
    {
        /* signal to my next neighbor that I have the data */
        tMPI_Event_signal( &(comm->csync[next].events[myrank]) );
        /* and wait for my next neighbor to finish */
        tMPI_Event_wait( &(comm->csync[myrank].events[myrank]) );
        tMPI_Event_process( &(comm->csync[myrank].events[myrank]), 1);
    }


#if defined(TMPI_PROFILE) && defined(TMPI_CYCLE_COUNT)
    tMPI_Profile_wait_start(cur);
#endif
    /*tMPI_Barrier_wait( &(comm->barrier));*/
#if defined(TMPI_PROFILE)
    /*tMPI_Profile_wait_stop(cur, TMPIWAIT_Reduce);*/
    tMPI_Profile_count_stop(cur, TMPIFN_Scan);
#endif
    return TMPI_SUCCESS;
}
