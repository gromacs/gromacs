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
#include "unused.h"

#ifdef USE_COLLECTIVE_COPY_BUFFER
/* initialize a copy buffer */
int tMPI_Copy_buffer_init(struct copy_buffer *cb, size_t size)
{
    cb->buf = tMPI_Malloc(size);
    if (cb->buf == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    cb->size = size;
    return TMPI_SUCCESS;
}

/* destroy a copy buffer */
void tMPI_Copy_buffer_destroy(struct copy_buffer *cb)
{
    free(cb->buf);
}

int tMPI_Copy_buffer_list_init(struct copy_buffer_list *cbl, int Nbufs,
                               size_t size)
{
    int i;
    int ret;

    cbl->size     = size;
    cbl->cb_alloc = (struct copy_buffer*)
        tMPI_Malloc(sizeof(struct copy_buffer)*Nbufs);
    if (cbl->cb_alloc == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    cbl->cb    = cbl->cb_alloc; /* the first one */
    cbl->Nbufs = Nbufs;
    for (i = 0; i < Nbufs; i++)
    {
        ret = tMPI_Copy_buffer_init( &(cbl->cb_alloc[i]), size );
        if (ret != TMPI_SUCCESS)
        {
            return ret;
        }
        if (i < Nbufs-1)
        {
            cbl->cb_alloc[i].next = &(cbl->cb_alloc[i+1]);
        }
        else
        {
            cbl->cb_alloc[i].next = NULL;
        }
    }
    return TMPI_SUCCESS;
}

void tMPI_Copy_buffer_list_destroy(struct copy_buffer_list *cbl)
{
    int i;

    for (i = 0; i < cbl->Nbufs; i++)
    {
        tMPI_Copy_buffer_destroy( &(cbl->cb_alloc[i]) );
    }
    free(cbl->cb_alloc);
}

struct copy_buffer *tMPI_Copy_buffer_list_get(struct copy_buffer_list *cbl)
{
    struct copy_buffer *ret = cbl->cb;
    if (!ret)
    {
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COPY_NBUFFERS);
        return NULL;
    }
    cbl->cb = ret->next;

    return ret;
}

void tMPI_Copy_buffer_list_return(struct copy_buffer_list *cbl,
                                  struct copy_buffer      *cb)
{
    cb->next = cbl->cb;
    cbl->cb  = cb;
}
#endif

int tMPI_Coll_envt_init(struct coll_env_thread *met, int N)
{
    tMPI_Atomic_set(&(met->current_sync), 0);
    tMPI_Atomic_set(&(met->n_remaining), 0);
    met->buf = (void**)tMPI_Malloc(sizeof(void*)*N);
    if (met->buf == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    met->bufsize = (size_t*)tMPI_Malloc(sizeof(size_t)*N);
    if (met->bufsize == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    met->read_data = (tmpi_bool*)tMPI_Malloc(sizeof(tmpi_bool)*N);
    if (met->read_data == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
#ifdef USE_COLLECTIVE_COPY_BUFFER
    met->cpbuf = (tMPI_Atomic_ptr_t*)tMPI_Malloc(sizeof(tMPI_Atomic_ptr_t)*
                                                 N);
    if (met->read_data == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    met->cb       = NULL;
    met->using_cb = FALSE;
#endif
    tMPI_Event_init( &(met->send_ev) );
    tMPI_Event_init( &(met->recv_ev) );
    return TMPI_SUCCESS;
}

void tMPI_Coll_envt_destroy(struct coll_env_thread *met)
{
    free( (void*)met->buf );
    free( (void*)met->bufsize );
    free( (void*)met->read_data );

#ifdef USE_COLLECTIVE_COPY_BUFFER
    free( (void*)met->cpbuf );
#endif
}

int tMPI_Coll_env_init(struct coll_env *cev, int N)
{
    int i;
    int ret;

    cev->met = (struct coll_env_thread*)tMPI_Malloc(
                sizeof(struct coll_env_thread)*N);
    if (cev->met == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    cev->N = N;
    tMPI_Atomic_set(&(cev->coll.current_sync), 0);
    tMPI_Atomic_set(&(cev->coll.n_remaining), 0);
    for (i = 0; i < N; i++)
    {
        ret = tMPI_Coll_envt_init(&(cev->met[i]), N);
        if (ret != TMPI_SUCCESS)
        {
            return ret;
        }
    }
    return TMPI_SUCCESS;
}

void tMPI_Coll_env_destroy(struct coll_env *cev)
{
    int i;
    for (i = 0; i < cev->N; i++)
    {
        tMPI_Coll_envt_destroy(&(cev->met[i]));
    }
    free( (void*)cev->met );
}

int tMPI_Coll_sync_init(struct coll_sync *csync, int N)
{
    int i;

    csync->synct = 0;
    csync->syncs = 0;
    csync->N     = N;

    csync->events = (tMPI_Event*)tMPI_Malloc(sizeof(tMPI_Event)*N);
    if (csync->events == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    for (i = 0; i < N; i++)
    {
        tMPI_Event_init( &(csync->events[i]) );
    }
    return TMPI_SUCCESS;
}

void tMPI_Coll_sync_destroy(struct coll_sync *csync)
{
    int i;

    csync->synct = 0;
    csync->syncs = 0;

    for (i = 0; i < csync->N; i++)
    {
        tMPI_Event_destroy( &(csync->events[i]) );
    }
    free(csync->events);
}

/* get a pointer the next coll_env once it's ready. */
struct coll_env *tMPI_Get_cev(tMPI_Comm comm, int myrank, int *counter)
{
    struct coll_sync *csync = &(comm->csync[myrank]);
    struct coll_env  *cev;
#ifdef USE_COLLECTIVE_COPY_BUFFER
    int               N;
#endif

    /* we increase our counter, and determine which coll_env we get */
    csync->synct++;
    *counter = csync->synct;
    cev      = &(comm->cev[csync->synct % N_COLL_ENV]);


#ifdef USE_COLLECTIVE_COPY_BUFFER
    if (cev->met[myrank].using_cb)
    {
        if (cev->N > 1)
        {
            N = tMPI_Event_wait( &(cev->met[myrank].send_ev));
            tMPI_Event_process( &(cev->met[myrank].send_ev), 1);
        }
    }
#endif
#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* clean up old copy_buffer pointers */
    if (cev->met[myrank].cb)
    {
        tMPI_Copy_buffer_list_return(&(tMPI_Get_current()->cbl_multi),
                                     cev->met[myrank].cb);
        cev->met[myrank].cb       = NULL;
        cev->met[myrank].using_cb = FALSE;
    }
#endif

    return cev;
}

void tMPI_Mult_recv(tMPI_Comm comm, struct coll_env *cev, int rank,
                    int index, int expected_tag, tMPI_Datatype recvtype,
                    size_t recvsize, void *recvbuf, int *ret)
{
    size_t sendsize = cev->met[rank].bufsize[index];

    /* check tags, types */
    if ((cev->met[rank].datatype != recvtype ) ||
        (cev->met[rank].tag != expected_tag))
    {
        *ret = tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
    }

    if (sendsize) /* we allow NULL ptrs if there's nothing to xmit */
    {
        void     *srcbuf;
#ifdef USE_COLLECTIVE_COPY_BUFFER
        tmpi_bool decrease_ctr = FALSE;
#endif

        if (sendsize > recvsize)
        {
            *ret = tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            return;
        }

        if (cev->met[rank].buf == recvbuf)
        {
            *ret = tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
            return;
        }
        /* get source buffer */
#ifdef USE_COLLECTIVE_COPY_BUFFER
        if (!(cev->met[rank].using_cb))
#endif
        {
            srcbuf = cev->met[rank].buf[index];
        }
#ifdef USE_COLLECTIVE_COPY_BUFFER
        else
        {
            srcbuf = tMPI_Atomic_ptr_get(&(cev->met[rank].cpbuf[index]));
            tMPI_Atomic_memory_barrier_acq();

            if (!srcbuf)
            {   /* there was (as of yet) no copied buffer */
                void *try_again_srcbuf;
                /* we need to try checking the pointer again after we increase
                   the read counter, signaling that one more thread
                   is reading. */
                tMPI_Atomic_fetch_add(&(cev->met[rank].buf_readcount), 1);
                /* a full memory barrier */
                tMPI_Atomic_memory_barrier();
                try_again_srcbuf = tMPI_Atomic_ptr_get(
                            &(cev->met[rank].cpbuf[index]));
                if (!try_again_srcbuf)
                {
                    /* apparently the copied buffer is not ready yet. We
                       just use the real source buffer. We have already
                       indicated we're reading from the regular buf. */
                    srcbuf       = cev->met[rank].buf[index];
                    decrease_ctr = TRUE;

                }
                else
                {
                    /* We tried again, and this time there was a copied buffer.
                       We use that, and indicate that we're not reading from the
                       regular buf. This case should be pretty rare.  */
                    tMPI_Atomic_fetch_add(&(cev->met[rank].buf_readcount), -1);
                    tMPI_Atomic_memory_barrier_acq();
                    srcbuf = try_again_srcbuf;
                }
            }

#ifdef TMPI_PROFILE
            if (srcbuf)
            {
                tMPI_Profile_count_buffered_coll_xfer(tMPI_Get_current());
            }
#endif
        }
#endif
        /* copy data */
        memcpy((char*)recvbuf, srcbuf, sendsize);
#ifdef TMPI_PROFILE
        tMPI_Profile_count_coll_xfer(tMPI_Get_current());
#endif

#ifdef USE_COLLECTIVE_COPY_BUFFER
        if (decrease_ctr)
        {
            /* we decrement the read count; potentially releasing the buffer. */
            tMPI_Atomic_memory_barrier_rel();
            tMPI_Atomic_fetch_add( &(cev->met[rank].buf_readcount), -1);
        }
#endif
    }
    /* signal one thread ready */
    {
        int reta;
        tMPI_Atomic_memory_barrier_rel();
        reta = tMPI_Atomic_fetch_add( &(cev->met[rank].n_remaining), -1);
        if (reta <= 1) /* n_remaining == 0 now. */
        {
            tMPI_Event_signal( &(cev->met[rank].send_ev) );
        }
    }
}

void tMPI_Coll_root_xfer(tMPI_Comm comm, tMPI_Datatype sendtype,
                         tMPI_Datatype recvtype,
                         size_t sendsize, size_t recvsize,
                         void* sendbuf, void* recvbuf, int *ret)
{
    /* do root transfer */
    if (recvsize < sendsize)
    {
        *ret = tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
        return;
    }
    if (recvtype != sendtype)
    {
        *ret = tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
        return;
    }
    if (sendbuf == recvbuf)
    {
        *ret = tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
        return;
    }

    memcpy(recvbuf, sendbuf, sendsize);
}

int tMPI_Post_multi(struct coll_env *cev, int myrank, int index,
                    int tag, tMPI_Datatype datatype, size_t bufsize,
                    void *buf, int n_remaining, int synct, int dest)
{
    int i;
#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* decide based on the number of waiting threads */
    tmpi_bool using_cb = (bufsize < (size_t)(n_remaining*COPY_BUFFER_SIZE));

    cev->met[myrank].using_cb = using_cb;
    if (using_cb)
    {
        /* we set it to NULL initially */
        /*cev->met[myrank].cpbuf[index]=NULL;*/
        tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[index]), NULL);

        tMPI_Atomic_set(&(cev->met[myrank].buf_readcount), 0);
    }
#endif
    cev->met[myrank].tag            = tag;
    cev->met[myrank].datatype       = datatype;
    cev->met[myrank].buf[index]     = buf;
    cev->met[myrank].bufsize[index] = bufsize;
    tMPI_Atomic_set(&(cev->met[myrank].n_remaining), n_remaining);
    tMPI_Atomic_memory_barrier_rel();
    tMPI_Atomic_set(&(cev->met[myrank].current_sync), synct);

    /* publish availability. */
    if (dest < 0)
    {
        for (i = 0; i < cev->N; i++)
        {
            if (i != myrank)
            {
                tMPI_Event_signal( &(cev->met[i].recv_ev) );
            }
        }
    }
    else
    {
        tMPI_Event_signal( &(cev->met[dest].recv_ev) );
    }

#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* becase we've published availability, we can start copying --
       possibly in parallel with the receiver */
    if (using_cb)
    {
        struct tmpi_thread *cur = tMPI_Get_current();
        /* copy the buffer locally. First allocate */
        cev->met[myrank].cb = tMPI_Copy_buffer_list_get( &(cur->cbl_multi) );
        if (cev->met[myrank].cb == NULL)
        {
            return TMPI_ERR_COPY_NBUFFERS;
        }
        if (cev->met[myrank].cb->size < bufsize)
        {
            return TMPI_ERR_COPY_BUFFER_SIZE;
        }
        /* copy to the new buf */
        memcpy(cev->met[myrank].cb->buf, buf, bufsize);

        /* post the new buf */
        tMPI_Atomic_memory_barrier_rel();
        /*cev->met[myrank].cpbuf[index]=cev->met[myrank].cb->buf;*/
        tMPI_Atomic_ptr_set(&(cev->met[myrank].cpbuf[index]),
                            cev->met[myrank].cb->buf);
    }
#endif
    return TMPI_SUCCESS;
}

void tMPI_Wait_for_others(struct coll_env *cev, int myrank)
{
#if defined(TMPI_PROFILE)
    struct tmpi_thread *cur = tMPI_Get_current();
    tMPI_Profile_wait_start(cur);
#endif

    if (cev->N > 1)
    {
#ifdef USE_COLLECTIVE_COPY_BUFFER
        if (!(cev->met[myrank].using_cb) )
#endif
        {
            /* wait until everybody else is done copying the buffer */
            tMPI_Event_wait( &(cev->met[myrank].send_ev));
            tMPI_Event_process( &(cev->met[myrank].send_ev), 1);
        }
#ifdef USE_COLLECTIVE_COPY_BUFFER
        else
        {
            /* wait until everybody else is done copying the original buffer.
               This wait is bound to be very short (otherwise it wouldn't
               be double-buffering) so we always spin here. */
#if 0
            /* dummy compare-and-swap to a value that is non-zero. The
               atomic read with barrier below is simpler, but we keep this
               code here commented out for if there is ever a platform
               where the simple read doesn't work because of, say, cache
               coherency issues. */
            while (!tMPI_Atomic_cas( &(cev->met[rank].buf_readcount), 0,
                                     -100000))
#endif
#if 1
            tMPI_Atomic_memory_barrier();         /* a full barrier to make
                                                     sure that the sending
                                                     doesn't interfere with the
                                                     waiting */
            while (tMPI_Atomic_get( &(cev->met[myrank].buf_readcount) ) > 0)
#endif
            {
                tMPI_Atomic_memory_barrier_acq();
            }
            tMPI_Atomic_memory_barrier_acq();
        }
#endif
    }
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_Coll_send);
#endif
}

void tMPI_Wait_for_data(struct tmpi_thread tmpi_unused *cur, struct coll_env *cev,
                        int myrank)
{
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_start(cur);
#endif
    tMPI_Event_wait( &(cev->met[myrank].recv_ev));
    tMPI_Event_process( &(cev->met[myrank].recv_ev), 1);
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_Coll_recv);
#endif
}

int tMPI_Barrier(tMPI_Comm comm)
{
#ifdef TMPI_PROFILE
    struct tmpi_thread *cur = tMPI_Get_current();
    tMPI_Profile_count_start(cur);
#endif

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Barrier(%p, %d, %p, %d, %d, %p, %p)", comm);
#endif

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (comm->grp.N > 1)
    {
#if defined(TMPI_PROFILE)
        tMPI_Profile_wait_start(cur);
#endif

        tMPI_Barrier_wait( &(comm->barrier) );
#if defined(TMPI_PROFILE)
        tMPI_Profile_wait_stop(cur, TMPIWAIT_Barrier);
#endif
    }

#ifdef TMPI_PROFILE
    tMPI_Profile_count_stop(cur, TMPIFN_Barrier);
#endif
    return TMPI_SUCCESS;
}
