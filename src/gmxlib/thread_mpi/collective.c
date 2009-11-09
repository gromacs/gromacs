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

#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "thread_mpi/collective.h"
#include "tmpi_impl.h"


/* get a pointer the next coll_env once it's ready */
static struct coll_env *tMPI_Get_cev(tMPI_Comm comm, int myrank, int *synct);

/* post the availability of data in a cev */
static void tMPI_Post_multi(struct coll_env *cev, int rank, int index, 
                            int tag, tMPI_Datatype datatype, 
                            size_t bufsize, void *buf, int n_remaining, 
                            int synct);

/* transfer data from cev->met[rank] to recvbuf */
static void tMPI_Mult_recv(tMPI_Comm comm, struct coll_env *cev, int rank,
                           int index, int expected_tag, tMPI_Datatype recvtype,
                           size_t recvsize, void *recvbuf, int *ret);

/* do a root transfer (from root send buffer to root recv buffer) */
static void tMPI_Coll_root_xfer(tMPI_Comm comm, 
                                tMPI_Datatype sendtype, tMPI_Datatype recvtype, 
                                size_t sendsize, size_t recvsize, 
                                void* sendbuf, void* recvbuf, int *ret);

/* wait for other processes to copy data from my cev */
static void tMPI_Wait_for_others(struct coll_env *cev, int rank);
/* wait for data to become available from a specific rank */
static void tMPI_Wait_for_data(struct coll_env *cev, int rank, int synct);

/* run a single binary reduce operation on src_a and src_b, producing dest. 
      dest and src_a may be identical */
static int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b,
                              tMPI_Datatype datatype, int count, tMPI_Op op,
                              tMPI_Comm comm);




#ifdef USE_COLLECTIVE_COPY_BUFFER
/* initialize a copy buffer */
void tMPI_Copy_buffer_init(struct copy_buffer *cb, size_t size)
{
    cb->buf=tMPI_Malloc(size);
    cb->size=size;
}

/* destroy a copy buffer */
void tMPI_Copy_buffer_destroy(struct copy_buffer *cb)
{
    free(cb->buf);
}

void tMPI_Copy_buffer_list_init(struct copy_buffer_list *cbl, int Nbufs,
                                size_t size)
{
    int i;

    cbl->size=size;
    cbl->cb_alloc=(struct copy_buffer*)
                  tMPI_Malloc(sizeof(struct copy_buffer)*Nbufs);
    cbl->cb=cbl->cb_alloc; /* the first one */
    cbl->Nbufs = Nbufs;
    for(i=0;i<Nbufs;i++)
    {
        tMPI_Copy_buffer_init( &(cbl->cb_alloc[i]), size );
        if (i<Nbufs-1)
            cbl->cb_alloc[i].next=&(cbl->cb_alloc[i+1]);
        else
            cbl->cb_alloc[i].next=NULL;
    }
}

void tMPI_Copy_buffer_list_destroy(struct copy_buffer_list *cbl)
{
    int i;

    for(i=0;i<cbl->Nbufs;i++)
    {
        tMPI_Copy_buffer_destroy( &(cbl->cb_alloc[i]) );
    }
    free(cbl->cb_alloc);
}

struct copy_buffer *tMPI_Copy_buffer_list_get(struct copy_buffer_list *cbl)
{
    struct copy_buffer *ret=cbl->cb;
    if (!ret)
    {
        fprintf(stderr,"out of copy buffers!!");
        exit(1);
    }
    cbl->cb=ret->next;

    return ret;
}

void tMPI_Copy_buffer_list_return(struct copy_buffer_list *cbl, 
                                  struct copy_buffer *cb)
{
    cb->next=cbl->cb;
    cbl->cb=cb;
}
#endif








static void tMPI_Coll_envt_init(struct coll_env_thread *met, int N)
{
    tMPI_Atomic_set(&(met->current_sync), 0);
    tMPI_Atomic_set(&(met->n_remaining), 0);
    met->buf=(void**)tMPI_Malloc(sizeof(void*)*N);
    met->bufsize=(size_t*)tMPI_Malloc(sizeof(size_t)*N);
    met->read_data=(bool*)tMPI_Malloc(sizeof(bool)*N);
#ifdef USE_COLLECTIVE_COPY_BUFFER
    met->cpbuf=(tMPI_Atomic_ptr_t*)tMPI_Malloc(sizeof(tMPI_Atomic_ptr_t)*N);
    met->cb=NULL;
    met->using_cb=FALSE;
#endif
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_init( &(met->wait_mutex) );
    tMPI_Thread_cond_init( &(met->send_cond) );
    tMPI_Thread_cond_init( &(met->recv_cond) );
#endif
}


static void tMPI_Coll_envt_destroy(struct coll_env_thread *met)
{
    free( (void*)met->buf );
    free( (void*)met->bufsize );
    free( (void*)met->read_data );

#ifdef USE_COLLECTIVE_COPY_BUFFER
    free( (void*)met->cpbuf );
#endif
}

void tMPI_Coll_env_init(struct coll_env *cev, int N)
{
    int i;

    cev->met=(struct coll_env_thread*)tMPI_Malloc(
                        sizeof(struct coll_env_thread)*N);
    cev->N=N;
    tMPI_Atomic_set(&(cev->coll.current_sync), 0);
    tMPI_Atomic_set(&(cev->coll.n_remaining), 0);
    for(i=0;i<N;i++)
    {
        tMPI_Coll_envt_init(&(cev->met[i]), N);
    }
}

void tMPI_Coll_env_destroy(struct coll_env *cev)
{
    int i;
    for(i=0;i<cev->N;i++)
    {
        tMPI_Coll_envt_destroy(&(cev->met[i]));
    }
    free( (void*)cev->met );
}


void tMPI_Coll_sync_init(struct coll_sync *csync)
{
    csync->synct=0;
    csync->syncs=0;
}

void tMPI_Coll_sync_destroy(struct coll_sync *csync)
{
    csync->synct=0;
    csync->syncs=0;
}






/* get a pointer the next coll_env once it's ready. */
static struct coll_env *tMPI_Get_cev(tMPI_Comm comm, int myrank, int *counter)
{
    struct coll_sync *csync=&(comm->csync[myrank]);
    struct coll_env *cev;

    /* we increase our counter, and determine which coll_env we get */
    csync->synct++;
    *counter=csync->synct;
    cev=&(comm->cev[csync->synct % N_COLL_ENV]);

#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_lock( &(cev->met[myrank].wait_mutex ) );
#endif

    /* wait until the thread's cev->met becomes available */
    while (tMPI_Atomic_get( &(cev->met[myrank].n_remaining) ) > 0)
    {
#ifdef TMPI_NO_BUSY_WAIT
        tMPI_Thread_cond_wait(&(cev->met[myrank].send_cond), 
                              &(cev->met[myrank].wait_mutex) );
#else
        tMPI_Atomic_memory_barrier();
#endif
    }
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_unlock( &(cev->met[myrank].wait_mutex ) );
#endif
#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* clean up old copy_buffer pointers */
    if (cev->met[myrank].cb)  
    {
        tMPI_Copy_buffer_list_return(&(tMPI_Get_current()->cbl_multi),
                                     cev->met[myrank].cb);
        cev->met[myrank].cb=NULL;
        cev->met[myrank].using_cb=FALSE;
    }
#endif

    return cev;
}





static void tMPI_Mult_recv(tMPI_Comm comm, struct coll_env *cev, int rank,
                           int index, int expected_tag, tMPI_Datatype recvtype, 
                           size_t recvsize, void *recvbuf, int *ret)
{
    size_t sendsize=cev->met[rank].bufsize[index];

    /* check tags, types */
    if ((cev->met[rank].datatype != recvtype ) || 
        (cev->met[rank].tag != expected_tag))
    {
        *ret=tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
    }
  
    if (sendsize) /* we allow NULL ptrs if there's nothing to xmit */
    {
        void *srcbuf;
#ifdef USE_COLLECTIVE_COPY_BUFFER
        bool decrease_ctr=FALSE;
#endif

        if ( sendsize > recvsize ) 
        {
            *ret=tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            return;
        }

        if ( cev->met[rank].buf == recvbuf )
        {
            *ret=tMPI_Error(TMPI_COMM_WORLD,TMPI_ERR_XFER_BUF_OVERLAP);
            return;
        }
        /* get source buffer */
#ifdef USE_COLLECTIVE_COPY_BUFFER
        if ( !(cev->met[rank].using_cb)) 
#endif
        {
            srcbuf=cev->met[rank].buf[index];
        }
#ifdef USE_COLLECTIVE_COPY_BUFFER
        else
        {
            tMPI_Atomic_memory_barrier();
            /*srcbuf=(char*) (cev->met[rank].cpbuf[index]);*/
            srcbuf=tMPI_Atomic_ptr_get(&(cev->met[rank].cpbuf[index]));

            if(!srcbuf)
            { /* there was (as of yet) no copied buffer */
                char *try_again_srcbuf;
                /* we need to try checking the pointer again after we increase
                   the read counter, signaling that one more thread
                   is reading. */
                tMPI_Atomic_add_return(&(cev->met[rank].buf_readcount), 1);
                tMPI_Atomic_memory_barrier();
                /*try_again_srcbuf=(char*) (cev->met[rank].cpbuf[index]);*/
                try_again_srcbuf=tMPI_Atomic_ptr_get(
                                         &(cev->met[rank].cpbuf[index]));
                if (!try_again_srcbuf)
                {
                    /* apparently the copied buffer is not ready yet. We
                       just use the real source buffer. We have already
                       indicated we're reading from the regular buf. */
                    srcbuf=cev->met[rank].buf[index];
                    decrease_ctr=TRUE;
                }
                else
                {
                    /* We tried again, and this time there was a copied buffer. 
                       We use that, and indicate that we're not reading from the
                       regular buf. This case should be pretty rare.  */
                    tMPI_Atomic_fetch_add(&(cev->met[rank].buf_readcount),-1);
                    srcbuf=try_again_srcbuf;
                }
            }
        }
#endif
        /* copy data */
        memcpy((char*)recvbuf, srcbuf, sendsize);

#ifdef USE_COLLECTIVE_COPY_BUFFER
        if (decrease_ctr)
        {
            /* we decrement the read count; potentially releasing the buffer. */
            tMPI_Atomic_fetch_add( &(cev->met[rank].buf_readcount), -1);
        }
#endif
    }
    /* signal one thread ready */
#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Atomic_fetch_add( &(cev->met[rank].n_remaining), -1);
#else
    tMPI_Thread_mutex_lock( &(cev->met[rank].wait_mutex ) );
    if (tMPI_Atomic_add_return( &(cev->met[rank].n_remaining), -1) <= 0)
    {
        tMPI_Thread_cond_broadcast( &(cev->met[rank].send_cond) );
    }
    tMPI_Thread_mutex_unlock( &(cev->met[rank].wait_mutex ) );
#endif
}

static void tMPI_Coll_root_xfer(tMPI_Comm comm, tMPI_Datatype sendtype, 
                                tMPI_Datatype recvtype, 
                                size_t sendsize, size_t recvsize, 
                                void* sendbuf, void* recvbuf, int *ret)
{
    /* do root transfer */
    if (recvsize < sendsize)
    {
        *ret=tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
        return;
    }
    if (recvtype != sendtype)
    {
        *ret=tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
        return;
    }
    if ( sendbuf == recvbuf )
    {
        *ret=tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
        return;
    }

    memcpy(recvbuf, sendbuf, sendsize);
}

static void tMPI_Post_multi(struct coll_env *cev, int rank, int index, 
                            int tag, tMPI_Datatype datatype, size_t bufsize, 
                            void *buf, int n_remaining, int synct)
{
#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* decide based on the number of waiting threads */
    bool using_cb=(bufsize < (size_t)(n_remaining*COPY_BUFFER_SIZE));

    cev->met[rank].using_cb=using_cb;
    if (using_cb)
    {
        /* we set it to NULL initially */
        /*cev->met[rank].cpbuf[index]=NULL;*/
        tMPI_Atomic_ptr_set(&(cev->met[rank].cpbuf[index]), NULL);

        tMPI_Atomic_set(&(cev->met[rank].buf_readcount), 0);
    }
#endif
    cev->met[rank].tag=tag;
    cev->met[rank].datatype=datatype;
    cev->met[rank].buf[index]=buf;
    cev->met[rank].bufsize[index]=bufsize;
    tMPI_Atomic_set(&(cev->met[rank].n_remaining), n_remaining);

    /* publish availability. */
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_lock( &(cev->met[rank].wait_mutex) );
#else
    tMPI_Atomic_memory_barrier();
#endif
    tMPI_Atomic_set(&(cev->met[rank].current_sync), synct);
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_cond_broadcast( &(cev->met[rank].recv_cond) );
    tMPI_Thread_mutex_unlock( &(cev->met[rank].wait_mutex) );
/*#else
    tMPI_Atomic_memory_barrier();*/
#endif

#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* becase we've published availability, we can start copying -- 
       possibly in parallel with the receiver */
    if (using_cb)
    {
         /* copy the buffer locally. First allocate */
        cev->met[rank].cb=tMPI_Copy_buffer_list_get(
                          &(tMPI_Get_current()->cbl_multi) );
        if (cev->met[rank].cb->size < bufsize)
        {
            fprintf(stderr, "ERROR: cb size too small\n");
            exit(1);
        }
        /* copy to the new buf */
        memcpy(cev->met[rank].cb->buf, buf, bufsize);
        /* post the new buf */
        tMPI_Atomic_memory_barrier();
        /*cev->met[rank].cpbuf[index]=cev->met[rank].cb->buf;*/
        tMPI_Atomic_ptr_set(&(cev->met[rank].cpbuf[index]), 
                            cev->met[rank].cb->buf);
    }
#endif
}


static void tMPI_Wait_for_others(struct coll_env *cev, int rank)
{
#ifdef USE_COLLECTIVE_COPY_BUFFER
    if (! (cev->met[rank].using_cb) )
#endif
    {
#ifdef TMPI_NO_BUSY_WAIT
        tMPI_Thread_mutex_lock( &(cev->met[rank].wait_mutex ) );
#endif
        /* wait until everybody else is done copying the buffer */
        while (tMPI_Atomic_get( &(cev->met[rank].n_remaining) ) > 0)
        {
#ifdef TMPI_NO_BUSY_WAIT
            tMPI_Thread_cond_wait(&(cev->met[rank].send_cond), 
                                  &(cev->met[rank].wait_mutex) );
#else
            tMPI_Atomic_memory_barrier();
#endif
        }
#ifdef TMPI_NO_BUSY_WAIT
        tMPI_Thread_mutex_unlock( &(cev->met[rank].wait_mutex ) );
#endif
#if 0
        /* wait until everybody else is done copying the buffer */
        while(tMPI_Atomic_get( &(cev->met[rank].n_remaining) ) > 0)
        {
        }
#endif
    }
#ifdef USE_COLLECTIVE_COPY_BUFFER
    else
    {
        /* wait until everybody else is done copying the original buffer. 
           We use fetch_add because we want to be sure of coherency.
           This wait is bound to be very short (otherwise it wouldn't 
           be double-buffering) so we always spin here. */
        tMPI_Atomic_memory_barrier();
#if 0
        while (tMPI_Atomic_cas( &(cev->met[rank].buf_readcount), 0,
                                    -100000) != 0)
#endif
#if 1
        while (tMPI_Atomic_fetch_add( &(cev->met[rank].buf_readcount), 0) != 0)
#endif
#if 0
        while (tMPI_Atomic_get( &(cev->met[rank].buf_readcount) )>0)
#endif
        {
            tMPI_Atomic_memory_barrier();
        }
    }
#endif
}

static void tMPI_Wait_for_data(struct coll_env *cev, int rank, int synct)
{
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_lock( &(cev->met[rank].wait_mutex ) );
#endif
    /* wait until the source posts availability by updating its sync value */
    while( tMPI_Atomic_get( &(cev->met[rank].current_sync)) != synct)
    {
#ifdef TMPI_NO_BUSY_WAIT
        tMPI_Thread_cond_wait(&(cev->met[rank].recv_cond), 
                              &(cev->met[rank].wait_mutex) );
#else
        tMPI_Atomic_memory_barrier();
#endif
    }
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_unlock( &(cev->met[rank].wait_mutex ) );
#endif

#if 0
    while( tMPI_Atomic_get( &(cev->met[rank].current_sync)) != synct)
    {
        tMPI_Atomic_memory_barrier();
    }
#endif
}






int tMPI_Barrier(tMPI_Comm comm) 
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Barrier(%p, %d, %p, %d, %d, %p, %p)", comm);
#endif

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (comm->grp.N>1)
    {
#ifndef TMPI_NO_BUSY_WAIT
        tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
#else
        tMPI_Thread_barrier_wait( &(comm->multicast_barrier[0]) );
#endif
    }
    return TMPI_SUCCESS;
}




/* The actual collective functions are #included, so that the static
   functions above are available to them and can get inlined if the
   compiler deems it appropriate. */
#include "bcast.c"
#include "scatter.c"
#include "gather.c"
#include "alltoall.c"
#include "reduce.c"

