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



int tMPI_Scatter(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                 void* recvbuf, int recvcount, tMPI_Datatype recvtype,
                 int root, tMPI_Comm comm)
{
    int                 synct;
    struct coll_env    *cev;
    int                 myrank;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();

#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Scatter(%p, %d, %p, %p, %d, %p, %d, %p)",
                     sendbuf, sendcount, sendtype,
                     recvbuf, recvcount, recvtype, root, comm);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank = tMPI_Comm_seek_rank(comm, cur);

    /* we increase our counter, and determine which coll_env we get */
    cev = tMPI_Get_cev(comm, myrank, &synct);

    if (myrank == root)
    {
        int       i;
        size_t    sendsize        = sendtype->size*sendcount;
        size_t    total_send_size = 0;
#ifdef USE_COLLECTIVE_COPY_BUFFER
        tmpi_bool using_cb;
#endif

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }

        /* we set up multiple posts, so no Post_multi */
        cev->met[myrank].tag      = TMPI_SCATTER_TAG;
        cev->met[myrank].datatype = sendtype;
        tMPI_Atomic_memory_barrier_rel();
        tMPI_Atomic_set( &(cev->met[myrank].n_remaining), cev->N-1 );
        for (i = 0; i < comm->grp.N; i++)
        {
            total_send_size            += sendtype->size*sendcount;
            cev->met[myrank].bufsize[i] = sendsize;
            cev->met[myrank].buf[i]     = (char*)sendbuf+sendsize*i;
        }
#ifdef USE_COLLECTIVE_COPY_BUFFER
        /* we must copy our own data too, unfortunately (otherwise there's
           a hole) */
        using_cb                  = (total_send_size < (size_t)((cev->N)*COPY_BUFFER_SIZE));
        cev->met[myrank].using_cb = using_cb;
        if (using_cb)
        {
            /* we set cpbuf stuff to NULL initially */
            for (i = 0; i < cev->N; i++)
            {
                /*cev->met[myrank].cpbuf[i]=NULL;*/
                tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[i]), NULL);
            }
            tMPI_Atomic_set(&(cev->met[myrank].buf_readcount), 0);
        }
#endif

        /* post availability */
        for (i = 0; i < cev->N; i++)
        {
            if (i != myrank)
            {
                tMPI_Event_signal( &(cev->met[i].recv_ev) );
            }
        }




#ifdef USE_COLLECTIVE_COPY_BUFFER
        /* do the copy_buffer thing */
        if (using_cb)
        {
            /* copy the buffer locally. First allocate */
            cev->met[myrank].cb = tMPI_Copy_buffer_list_get(
                        &(tMPI_Get_thread(comm, myrank)->cbl_multi));
            if (cev->met[myrank].cb->size < total_send_size)
            {
                return tMPI_Error(comm, TMPI_ERR_COPY_BUFFER_SIZE);
            }
            /* copy to the new buf */
            memcpy(cev->met[myrank].cb->buf, sendbuf, total_send_size);

            /* post the new buf */
            for (i = 0; i < cev->N; i++)
            {
                tMPI_Atomic_memory_barrier_rel();
                tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[i]),
                                    (char*)cev->met[myrank].cb->buf+sendsize*i);
                /*cev->met[myrank].cpbuf[i] = (char*)cev->met[myrank].cb->buf +
                                            sendsize*i ;*/
            }
        }
#endif

        /* do root transfer */
        if (recvbuf != TMPI_IN_PLACE)
        {
            tMPI_Coll_root_xfer(comm, sendtype, recvtype,
                                sendsize, recvtype->size*recvcount,
                                (char*)sendbuf+sendsize*myrank,
                                recvbuf, &ret);
        }

        /* and wait until everybody is done copying */
        tMPI_Wait_for_others(cev, myrank);
    }
    else
    {
        /* get the root cev */
        size_t bufsize = recvcount*recvtype->size;
        /* wait until root becomes available */
        tMPI_Wait_for_data(cur, cev, myrank);
        tMPI_Mult_recv(comm, cev, root, myrank, TMPI_SCATTER_TAG, recvtype,
                       bufsize, recvbuf, &ret);
    }
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Scatter);
#endif
    return ret;
}



int tMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs,
                  tMPI_Datatype sendtype, void* recvbuf, int recvcount,
                  tMPI_Datatype recvtype, int root, tMPI_Comm comm)
{
    int                 synct;
    struct coll_env    *cev;
    int                 myrank;
    int                 ret = TMPI_SUCCESS;
    struct tmpi_thread *cur = tMPI_Get_current();

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Scatterv(%p, %p, %p, %p, %p, %d, %p, %d, %p)",
                     sendbuf, sendcounts, displs, sendtype, recvbuf,
                     recvcount, recvtype, root, comm);
#endif
#ifdef TMPI_PROFILE
    tMPI_Profile_count_start(cur);
#endif


    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank = tMPI_Comm_seek_rank(comm, cur);

    /* we increase our counter, and determine which coll_env we get */
    cev = tMPI_Get_cev(comm, myrank, &synct);

    if (myrank == root)
    {
        int       i;
        size_t    total_send_size = 0;
#ifdef USE_COLLECTIVE_COPY_BUFFER
        tmpi_bool using_cb;
#endif

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }

        /* we set up multiple posts, so no Post_multi */
        cev->met[myrank].tag      = TMPI_SCATTERV_TAG;
        cev->met[myrank].datatype = sendtype;
        tMPI_Atomic_memory_barrier_rel();
        tMPI_Atomic_set( &(cev->met[myrank].n_remaining), cev->N-1 );
        for (i = 0; i < cev->N; i++)
        {
            total_send_size            += sendtype->size*sendcounts[i];
            cev->met[myrank].bufsize[i] = sendtype->size*sendcounts[i];
            cev->met[myrank].buf[i]     = (char*)sendbuf+sendtype->size*displs[i];
        }
#ifdef USE_COLLECTIVE_COPY_BUFFER
        /* we must copy our own data too, unfortunately (otherwise there's
           a hole) */
        using_cb                  = (total_send_size < (size_t)((cev->N)*COPY_BUFFER_SIZE));
        cev->met[myrank].using_cb = using_cb;
        if (using_cb)
        {
            /* we set cpbuf stuff to NULL initially */
            for (i = 0; i < cev->N; i++)
            {
                /*cev->met[myrank].cpbuf[i]=NULL;*/
                tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[i]), NULL);
            }
            tMPI_Atomic_set(&(cev->met[myrank].buf_readcount), 0);
        }
#endif

        /* post availability */
        for (i = 0; i < cev->N; i++)
        {
            if (i != myrank)
            {
                tMPI_Event_signal( &(cev->met[i].recv_ev) );
            }
        }



#ifdef USE_COLLECTIVE_COPY_BUFFER
        /* do the copy_buffer thing */
        if (using_cb)
        {
            /* copy the buffer locally. First allocate */
            cev->met[myrank].cb = tMPI_Copy_buffer_list_get(
                        &(tMPI_Get_thread(comm, myrank)->cbl_multi));
            if (cev->met[myrank].cb->size < total_send_size)
            {
                return tMPI_Error(comm, TMPI_ERR_COPY_BUFFER_SIZE);
            }
            /* copy to the new buf */
            memcpy(cev->met[myrank].cb->buf, sendbuf, total_send_size);
            /* post the new buf */
            for (i = 0; i < cev->N; i++)
            {
                tMPI_Atomic_memory_barrier_rel();
                tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[i]),
                                    (char*)cev->met[myrank].cb->buf +
                                    sendtype->size*displs[i]);
                /*cev->met[myrank].cpbuf[i]=(char*)cev->met[myrank].cb->buf +
                          sendtype->size*displs[i];*/
            }
        }
#endif

        /* do root transfer */
        if (recvbuf != TMPI_IN_PLACE)
        {
            tMPI_Coll_root_xfer(comm, sendtype, recvtype,
                                sendtype->size*sendcounts[myrank],
                                recvtype->size*recvcount,
                                (char*)sendbuf+sendtype->size*displs[myrank],
                                recvbuf,
                                &ret);
        }

        /* and wait until everybody is done copying */
        tMPI_Wait_for_others(cev, myrank);
    }
    else
    {
        /* get the root cev */
        size_t bufsize = recvcount*recvtype->size;
        /* wait until root becomes available */
        tMPI_Wait_for_data(cur, cev, myrank);
        tMPI_Mult_recv(comm, cev, root, myrank, TMPI_SCATTERV_TAG,
                       recvtype, bufsize, recvbuf, &ret);
    }
#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Scatterv);
#endif
    return ret;
}
