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

/* this file is #included from collective.c */


/* run a single binary reduce operation on src_a and src_b, producing dest. 
      dest and src_a may be identical */
static int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b,
                              tMPI_Datatype datatype, int count, tMPI_Op op,
                              tMPI_Comm comm);


static int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
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
                     tMPI_Datatype datatype, tMPI_Op op, int root, 
                     tMPI_Comm comm)
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

    if (sendbuf==TMPI_IN_PLACE)/* i.e. sendbuf == tMPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }
    comm->sendbuf[myrank]=sendbuf;
    comm->recvbuf[myrank]=recvbuf;
    /* there's a barrier to wait for all the processes to put their 
       send/recvbuf in the global list */
#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
#else
    tMPI_Thread_barrier_wait( &(comm->multicast_barrier[0]));
#endif

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
                   iteration. We need to check for overlapping buffers
                   here because MPI_IN_PLACE might cause recvbuf to be the
                   same as sendbuf. */
                if (iteration==0 && (recvbuf!=sendbuf))
                    memcpy(recvbuf, sendbuf, datatype->size*count);
            }
            /* split barrier */
#ifndef TMPI_NO_BUSY_WAIT
            tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
#else
            tMPI_Thread_barrier_wait( &(comm->multicast_barrier[iteration]));
#endif
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
#ifndef TMPI_NO_BUSY_WAIT
            tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
#else
            tMPI_Thread_barrier_wait( &(comm->multicast_barrier[iteration]));
#endif
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

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Reduce(%p, %p, %d, %p, %p, %d, %p)",
                       sendbuf, recvbuf, count, datatype, op, root, comm);
#endif

    if (myrank==root)
    {
        if (sendbuf==TMPI_IN_PLACE) /* i.e. sendbuf == TMPI_IN_PLACE */
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

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Allreduce(%p, %p, %d, %p, %p, %p)",
                     sendbuf, recvbuf, count, datatype, op, comm);
#endif
    if (count==0)
        return TMPI_SUCCESS;
    if (!recvbuf)
    {
        return tMPI_Error(comm, TMPI_ERR_BUF);
    }
    if (sendbuf==TMPI_IN_PLACE) /* i.e. sendbuf == TMPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }

    ret=tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, 0, comm);
    /* distribute rootbuf */
    rootbuf=(void*)comm->recvbuf[0];


#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
#else
    tMPI_Thread_barrier_wait( &(comm->multicast_barrier[0]));
#endif

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
#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Spinlock_barrier_wait( &(comm->multicast_barrier[0]));
#else
    tMPI_Thread_barrier_wait( &(comm->multicast_barrier[0]));
#endif

    return ret;
}

