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

To help us fund development, we humbly ask that you cite
any papers on the package - you can find them in the top README file.

*/


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#if 0
#include <sys/time.h>
#endif

#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "tmpi_impl.h"


/* do a mulitcast transfer, with checks */
static void tMPI_Mult_xfer(tMPI_Comm comm, int rank, struct multi_env *rmev,
                           void *recvbuf, size_t recvsize, 
                           tMPI_Datatype recvtype, int expected_tag, int *ret);



/* run a single binary reduce operation on src_a and src_b, producing dest. 
      dest and src_a may be identical */
static int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b,
                              tMPI_Datatype datatype, int count, tMPI_Op op,
                              tMPI_Comm comm);





void tMPI_Multi_sync_init(struct multi_sync *msc, int N)
{
    int i;

    msc->counter=0;
    for(i=0;i<N_MULTI_SYNC;i++)
    {
        tMPI_Atomic_set( &(msc->mev[i].current_counter), 0);
        tMPI_Atomic_set( &(msc->mev[i].n_remaining), 0);
        msc->mev[i].buf=(void*)tMPI_Malloc(sizeof(void*)*N);
        msc->mev[i].bufsize=(size_t*)tMPI_Malloc(sizeof(size_t)*N);
        msc->mev[i].read_data=(bool*)tMPI_Malloc(sizeof(bool)*N);
    }
}

void tMPI_Multi_sync_destroy(struct multi_sync *msc)
{
    int i;

    for(i=0;i<N_MULTI_SYNC;i++)
    {
        free((void*)msc->mev[i].buf);
        free((void*)msc->mev[i].bufsize);
        free((void*)msc->mev[i].read_data);
    }
}

static void tMPI_Mult_xfer(tMPI_Comm comm, int rank, struct multi_env *rmev, 
                           void *recvbuf, size_t recvsize, 
                           tMPI_Datatype recvtype, int expected_tag, int *ret)
{
    size_t sendsize=rmev->bufsize[rank];
    /* check tags, types */
    if ( (rmev->datatype != recvtype ) || (rmev->tag != expected_tag) )
    {
        *ret=tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
    }
  
    if (sendsize) /* we allow NULL ptrs if there's nothing to xmit */
    {
        if ( sendsize > recvsize ) 
        {
            *ret=tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
        }

        if ( rmev->buf[rank] == recvbuf )
        {
            *ret=tMPI_Error(TMPI_COMM_WORLD,TMPI_ERR_XFER_BUF_OVERLAP);
        }
        /* copy data */
        memcpy((char*)recvbuf, (void*)(rmev->buf[rank]), sendsize);
    }
    /* signal one thread ready */
    tMPI_Atomic_fetch_add( &(rmev->n_remaining), -1);
}


int tMPI_Barrier(tMPI_Comm comm) 
{
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (comm->grp.N>1)
    {
        tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
    }
    return TMPI_SUCCESS;
}


/* multi */
int tMPI_Bcast(void* buffer, int count, tMPI_Datatype datatype, int root,
              tMPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=TMPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        struct multi_env *mev=&(msc->mev[mevi]);

        /* first set up the data */
        mev->tag=TMPI_BCAST_TAG;
        mev->datatype=datatype;
        mev->buf[0]=buffer;
        mev->bufsize[0]=count*datatype->size;
        tMPI_Atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);
        /* and wait until everybody is done copying */
        while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=count*datatype->size;
        /* wait until root becomes available */
        while( tMPI_Atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }

        tMPI_Mult_xfer(comm, 0, rmev, buffer, bufsize, datatype, 
                       TMPI_BCAST_TAG, &ret);
    }
    return ret;
}





int tMPI_Gather(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                void* recvbuf, int recvcount, tMPI_Datatype recvtype, 
                int root, tMPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=TMPI_SUCCESS;
    struct multi_env *mev;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);
    if (myrank==root)
    {
        int i;
        int n_remaining=comm->grp.N-1;
        /* do root transfer */
        {
            size_t recvsize=recvtype->size*recvcount;
            size_t sendsize=sendtype->size*sendcount;
            if (recvsize < sendsize)
            {
                return tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
            }
            if ( sendbuf == (char*)recvbuf+myrank*recvcount*recvtype->size )
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy((char*)recvbuf+myrank*recvcount*recvtype->size, sendbuf, 
                   sendsize);
        }
        /* poll data availability */
        while(n_remaining>0)
        {
            for(i=0;i<comm->grp.N;i++)
            {
                struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
                if ((i!=myrank) && 
                    (tMPI_Atomic_get(&(rmev->current_counter))==msc->counter)&& 
                    (tMPI_Atomic_get(&(rmev->n_remaining)) > 0) )
                {
                    tMPI_Mult_xfer(comm, 0, rmev, 
                                   (char*)recvbuf+i*recvcount*recvtype->size,
                                   recvcount*recvtype->size,
                                   recvtype, TMPI_GATHERV_TAG,&ret);
                    if (ret!=TMPI_SUCCESS)
                        return ret;
                    n_remaining--;
                }
            }
        }
    }
    else
    {
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }

        /* first set up the data */
        mev->tag=TMPI_GATHERV_TAG;
        mev->datatype=sendtype;
        mev->buf[0]=sendbuf;
        mev->bufsize[0]=bufsize;

        tMPI_Atomic_set(&(mev->n_remaining), 1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);

        /* and wait until root is done copying */
        while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }

    return ret;
}


int tMPI_Gatherv(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                void* recvbuf, int *recvcounts, int *displs,
                tMPI_Datatype recvtype, int root, tMPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=TMPI_SUCCESS;
    struct multi_env *mev;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);
    if (myrank==root)
    {
        int i;
        int n_remaining=comm->grp.N-1;
        /* do root transfer */
        {
            size_t recvsize=recvtype->size*recvcounts[myrank];
            size_t sendsize=sendtype->size*sendcount;
            if (recvsize < sendsize)
            {
                return tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
            }
            if ( sendbuf == (char*)recvbuf+displs[myrank]*recvtype->size )
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy((char*)recvbuf+displs[myrank]*recvtype->size, 
                   sendbuf, sendsize);
        }
        /* poll data availability */
        while(n_remaining>0)
        {
            for(i=0;i<comm->grp.N;i++)
            {
                struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
                if ((i!=myrank) && 
                    (tMPI_Atomic_get(&(rmev->current_counter))==msc->counter)&& 
                    (tMPI_Atomic_get(&(rmev->n_remaining)) > 0) )
                {
                    tMPI_Mult_xfer(comm, 0, rmev,
                                   (char*)recvbuf+displs[i]*recvtype->size,
                                   recvcounts[i]*recvtype->size,
                                   recvtype, TMPI_GATHERV_TAG,&ret);
                    if (ret!=TMPI_SUCCESS)
                        return ret;
                    n_remaining--;
                }
            }
        }
    }
    else
    {
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }

        /* first set up the data */
        mev->tag=TMPI_GATHERV_TAG;
        mev->datatype=sendtype;
        mev->buf[0]=sendbuf;
        mev->bufsize[0]=bufsize;

        tMPI_Atomic_set(&(mev->n_remaining), 1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);

        /* and wait until root is done copying */
        while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }

    return ret;
}


int tMPI_Scatter(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                 void* recvbuf, int recvcount, tMPI_Datatype recvtype, 
                 int root, tMPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=TMPI_SUCCESS;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        int i;
        struct multi_env *mev=&(msc->mev[mevi]);
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }
        /* first set up the data */
        mev->tag=TMPI_SCATTER_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf + bufsize*i;
            mev->bufsize[i]=bufsize;
        }
        tMPI_Atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);

        /* do the root transfer */
        {
            size_t recvsize=recvtype->size*recvcount;
            if (recvsize < bufsize)
            {
                return tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
            }
            if ( recvbuf == mev->buf[myrank] )
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy(recvbuf, (void*)(mev->buf[myrank]), bufsize);
        }

        /* and wait until everybody is done copying */
        while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=recvcount*recvtype->size;
        /* wait until root becomes available */
        while( tMPI_Atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }
        tMPI_Mult_xfer(comm, myrank, rmev, recvbuf, bufsize, recvtype,
                       TMPI_SCATTER_TAG, &ret);
    }
    return ret;
}



int tMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs,
                 tMPI_Datatype sendtype, void* recvbuf, int recvcount,
                 tMPI_Datatype recvtype, int root, tMPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=TMPI_SUCCESS;
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        int i;
        struct multi_env *mev=&(msc->mev[mevi]);

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, TMPI_ERR_BUF);
        }
        /* first set up the data */
        mev->tag=TMPI_SCATTERV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf + displs[i]*sendtype->size;
            mev->bufsize[i]=sendtype->size*sendcounts[i];
        }
        tMPI_Atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);

        /* do the root transfer */
        {
            size_t recvsize=recvtype->size*recvcount;
            if (recvsize < mev->bufsize[myrank])
            {
                return tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, TMPI_ERR_MULTI_MISMATCH);
            }
            if ( recvbuf == mev->buf[myrank] )
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy(recvbuf, (void*)(mev->buf[myrank]), mev->bufsize[myrank]);
        }

        /* and wait until everybody is done copying */
        while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=recvcount*recvtype->size;
        /* wait until root becomes available */
        while( tMPI_Atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }
        tMPI_Mult_xfer(comm, myrank, rmev, recvbuf, bufsize, recvtype,
                       TMPI_SCATTERV_TAG, &ret);
    }
    return ret;

}



int tMPI_Alltoall(void* sendbuf, int sendcount, tMPI_Datatype sendtype,
                 void* recvbuf, int recvcount, tMPI_Datatype recvtype,
                 tMPI_Comm comm)
{
    struct multi_sync *msc;
    struct multi_env *mev;
    int mevi;
    int i;
    int myrank;
    int n_remaining;
    int ret=TMPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }

    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);

    /* post our pointers */
    {
        /* first set up the data */
        mev->tag=TMPI_ALLTOALLV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf+i*sendcount*sendtype->size;
            mev->bufsize[i]=sendcount*sendtype->size;
            mev->read_data[i]=FALSE;
        }

        tMPI_Atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);
    }
    /* do root transfer */
    {
        size_t recvsize=recvtype->size*recvcount;
        tMPI_Mult_xfer(comm, myrank, mev,
                       (char*)recvbuf+myrank*recvcount*recvtype->size, 
                       recvsize, recvtype, TMPI_ALLTOALLV_TAG, &ret);
        if (ret!=TMPI_SUCCESS)
            return ret;
                    
        mev->read_data[myrank]=TRUE;
    }

    n_remaining=comm->grp.N-1;
    /* poll data availability */
    while(n_remaining>0)
    {
        for(i=0;i<comm->grp.N;i++)
        {
            struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
            if ((!mev->read_data[i]) &&
                (tMPI_Atomic_get(&(rmev->current_counter))==msc->counter) )
            {
                tMPI_Mult_xfer(comm, myrank, rmev, 
                               (char*)recvbuf+i*recvcount*recvtype->size,
                               recvcount*recvtype->size, 
                               recvtype, TMPI_ALLTOALLV_TAG,&ret);
                if (ret!=TMPI_SUCCESS)
                    return ret;

                mev->read_data[i]=TRUE;
                n_remaining--;
            }
        }
    }
    /* and wait until everybody is done copying our data */
    while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
    {
    }

    return ret;
}


int tMPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls,
                  tMPI_Datatype sendtype, 
                  void* recvbuf, int *recvcounts, int *rdispls, 
                  tMPI_Datatype recvtype, 
                  tMPI_Comm comm)

{
    struct multi_sync *msc;
    struct multi_env *mev;
    int mevi;
    int i;
    int myrank;
    int n_remaining;
    int ret=TMPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }

    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);

    /* post our pointers */
    {
        /* first set up the data */
        mev->tag=TMPI_ALLTOALLV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf+sdispls[i]*sendtype->size;
            mev->bufsize[i]=sendcounts[i]*sendtype->size;
            mev->read_data[i]=FALSE;
        }

        tMPI_Atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        tMPI_Atomic_set(&(mev->current_counter), msc->counter);
    }
    /* do root transfer */
    {
        size_t recvsize=recvtype->size*recvcounts[myrank];
        tMPI_Mult_xfer(comm, myrank, mev,
                       (char*)recvbuf+rdispls[myrank]*recvtype->size, 
                       recvsize, recvtype, TMPI_ALLTOALLV_TAG, &ret);
        if (ret!=TMPI_SUCCESS)
            return ret;
                    
        mev->read_data[myrank]=TRUE;
    }

    n_remaining=comm->grp.N-1;
    /* poll data availability */
    while(n_remaining>0)
    {
        for(i=0;i<comm->grp.N;i++)
        {
            struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
            if ((!mev->read_data[i]) &&
                (tMPI_Atomic_get(&(rmev->current_counter))==msc->counter) )
            {
                tMPI_Mult_xfer(comm, myrank, rmev, 
                               (char*)recvbuf+rdispls[i]*recvtype->size,
                               recvcounts[i]*recvtype->size, 
                               recvtype, TMPI_ALLTOALLV_TAG,&ret);
                if (ret!=TMPI_SUCCESS)
                    return ret;

                mev->read_data[i]=TRUE;
                n_remaining--;
            }
        }
    }
    /* and wait until everybody is done copying our data */
    while (tMPI_Atomic_get( &(mev->n_remaining) ) > 0)
    {
    }

    return ret;
}






int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
                        tMPI_Datatype datatype, int count, tMPI_Op op, 
                        tMPI_Comm comm)
{
    tMPI_Op_fn fn=datatype->op_functions[op];

    if (src_a==src_b)
    {
        return tMPI_Error(comm, TMPI_ERR_XFER_BUF_OVERLAP);
    }
    fn(dest, src_a, src_b, count);
    return TMPI_SUCCESS;
}

int tMPI_Reduce_fast(void* sendbuf, void* recvbuf, int count,
               tMPI_Datatype datatype, tMPI_Op op, int root, tMPI_Comm comm)
{
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int i;

    /* this function uses a binary tree-like reduction algorithm: */
    int N=tMPI_Comm_N(comm);
    int myrank_rtr=(N+myrank-root)%N; /* my rank relative to root */
    int Nnbrs=N; /* number of neighbours to communicate with 
                    (decreases exponentially) */
    int nbr_dist=1; /* distance between communicating neighbours 
                       (increases exponentially) */
    int stepping=2; /* distance between non-communicating neighbours
                       (increases exponentially) */
    int iteration=0;
  
    if (count==0)
        return TMPI_SUCCESS;
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

    if (!sendbuf)/* i.e. sendbuf == tMPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }
    comm->sendbuf[myrank]=sendbuf;
    comm->recvbuf[myrank]=recvbuf;
    /* there's a barrier to wait for all the processes to put their 
       send/recvbuf in the global list */
    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));

    /* check the buffers */
    for(i=0;i<N;i++)
    {
        if ( (i!=myrank) && ( (comm->recvbuf[i]==recvbuf) || 
                              (comm->sendbuf[i]==sendbuf) ) )
        {
            return tMPI_Error(comm, TMPI_ERR_XFER_BUF_OVERLAP);
        }
    }

    while(Nnbrs>1)
    {
        int nbr;
        /* reduce between myrank and myrank+nbr_dist, if there is such
            a neighbour. */
#ifdef TMPI_DEBUG
        printf("%d: iteration=%d, Nnbrs=%d, barrier size=%d and I'm still in here\n", 
                myrank, iteration, comm->N_multicast_barrier[iteration], Nnbrs);
        fflush(0);
#endif
        if (myrank_rtr%stepping == 0 )
        {
            if (myrank_rtr+nbr_dist<N) 
            {
                void *a,*b;
                int ret;
#ifdef TMPI_DEBUG
                printf("%d: reducing with %d, iteration=%d\n", 
                        myrank, myrank+nbr_dist, iteration);
                fflush(0);
#endif
               /* we reduce with our neighbour*/
                nbr=(N+myrank+nbr_dist)%N; /* here we use the real rank */
                if (iteration==0)
                {
                    a=sendbuf;
                    b=(void*)(comm->sendbuf[nbr]);
                }
                else
                {
                    a=recvbuf;
                    b=(void*)(comm->recvbuf[nbr]);
                }
                if ((ret=tMPI_Reduce_run_op(recvbuf, a, b, datatype, 
                                            count, op, comm)) != TMPI_SUCCESS)
                    return ret;
            }
            else
            {
                /* we still need to put things in the right buffer for the next
                   iteration */
                if (iteration==0)
                    memcpy(recvbuf, sendbuf, datatype->size*count);
            }
            /* split barrier */
            tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
        }
        else 
        {
            /* this barrier is split because the threads not actually 
               calculating still need to be waiting for the thread using its
               data to finish calculating */
#ifdef TMPI_DEBUG
            printf("%d: barrier waiting, iteration=%d\n", myrank,  iteration);
            fflush(0);
#endif
            tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
            break;
        }
#ifdef TMPI_DEBUG
        printf("%d: barrier over, iteration=%d\n", myrank,  iteration);
        fflush(0);
#endif
        Nnbrs = Nnbrs/2 + Nnbrs%2;
        nbr_dist*=2;
        stepping*=2;
        iteration++;
    }
    return TMPI_SUCCESS;
}

int tMPI_Reduce(void* sendbuf, void* recvbuf, int count,
               tMPI_Datatype datatype, tMPI_Op op, int root, tMPI_Comm comm)
{
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int ret;

    if (myrank==root)
    {
        if (!sendbuf) /* i.e. sendbuf == TMPI_IN_PLACE */
        {
            sendbuf=recvbuf;
        }
    }
    else
    {
#ifdef TMPI_WARN_MALLOC
        fprintf(stderr, "Warning: malloc during tMPI_Reduce\n");
#endif
        recvbuf=(void*)tMPI_Malloc(datatype->size*count);
    }
    ret=tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, root, comm);
    if (myrank!=root)
    {
        free(recvbuf);
    }
    return ret;
}

int tMPI_Allreduce(void* sendbuf, void* recvbuf, int count,
                  tMPI_Datatype datatype, tMPI_Op op, tMPI_Comm comm)
{
    void *rootbuf=NULL; /* root process' receive buffer */
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int ret;

    if (count==0)
        return TMPI_SUCCESS;
    if (!recvbuf)
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }
    if (!sendbuf) /* i.e. sendbuf == TMPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }

    ret=tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, 0, comm);
    /* distribute rootbuf */
    rootbuf=(void*)comm->recvbuf[0];

    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
    /* and now we just copy things back inefficiently. We should make
       a better tMPI_Scatter, and use that. */
    if (myrank != 0)
    {
        if (rootbuf==recvbuf)
        {
            return tMPI_Error(comm, TMPI_ERR_XFER_BUF_OVERLAP);
        }
        memcpy(recvbuf, rootbuf, datatype->size*count );
    }
    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));

    return ret;
}

