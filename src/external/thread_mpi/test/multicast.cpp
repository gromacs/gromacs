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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef MPICC
#include "tmpi.h"
#else
#include <mpi.h>
#endif

#include "send_recv.h"
#include "multicast.h"
#include "ntests.h"

#define MSG_SIZE 128

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif

static int msgsize(int rank)
{
    return rank + MSG_SIZE/2;
}


void barrier_data_tester(void)
{
    int                   myrank, N;
    int                   i, j, k;
    int                   testnr = 0;
    static tMPI_Atomic_t *shared_buf;

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("Barrier data test:\n");
        fflush(stdout);
        shared_buf = (tMPI_Atomic_t*)malloc(sizeof(tMPI_Atomic_t)*N);
        for (i = 0; i < N; i++)
        {
            tMPI_Atomic_set(&shared_buf[i], testnr);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    for (i = 0; i < BARRIER_NTRIES; i++)
    {
        for (j = 0; j < BARRIER_TRYSIZE; j++)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            /* add one, non-atomically */
            testnr = tMPI_Atomic_get(&shared_buf[myrank])+1;
            tMPI_Atomic_set(&shared_buf[myrank], testnr);

            MPI_Barrier(MPI_COMM_WORLD);

            for (k = 0; k < N; k++)
            {
                if (tMPI_Atomic_get(&shared_buf[k]) != testnr)
                {
                    fprintf(stderr, "ERROR in MPI_Barrier\n");
                    fprintf(stderr, "shared_buf[%d]=%d != %d (%d)\n",
                            k, tMPI_Atomic_get(&shared_buf[k]), testnr, myrank);
                    exit(0);
                }
            }
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
        free(shared_buf);
    }
}




void bcast_data_tester(void)
{
    int myrank, N;
    int i, j;
    int testnr = 0;
    int buf[MSG_SIZE];
    int buf_recv[MSG_SIZE];
    int buf_expt[MSG_SIZE];

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Broadcast data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < BCAST_NTRIES; i++)
    {
        for (j = 0; j < BCAST_TRYSIZE; j++)
        {
            int root = (testnr/3)%N;
            if (N < 2)
            {
                root = 0;    /* because otherwise we'll get stuck */

            }
            if (myrank == root)
            {
                unique_array(myrank, N, testnr, MSG_SIZE, buf);
                if (MPI_Bcast(buf, MSG_SIZE, MPI_INT, root,
                              MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Bcast"); exit(0);
                }
            }
            else
            {
                unique_array(root, N, testnr, MSG_SIZE, buf_expt);
                if (MPI_Bcast(buf_recv, MSG_SIZE, MPI_INT, root,
                              MPI_COMM_WORLD) != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Bcast"); exit(0);
                }

                if (!check_arrays(buf_recv, buf_expt, MSG_SIZE))
                {
                    print_match_err(myrank, buf_recv, buf_expt, MSG_SIZE);
                }
            }
            testnr++;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}



void gather_data_tester(void)
{
    int  myrank, N;
    int  i, j, k;
    int  testnr = 0;
    int *buf;
    int  buf_send[MSG_SIZE];
    int  buf_expt[MSG_SIZE];

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Gather data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    buf = (int*)malloc(sizeof(int)*MSG_SIZE*N);
    for (i = 0; i < GATHER_NTRIES; i++)
    {
        for (j = 0; j < GATHER_TRYSIZE; j++)
        {
            int root = (testnr/3)%N; /* we alter root every third test */
            if (N < 2)
            {
                root = 0;    /* because otherwise we'll get stuck */

            }
            unique_array(myrank, N, testnr, MSG_SIZE, buf_send);

            if (MPI_Gather(buf_send, MSG_SIZE, MPI_INT,
                           buf, MSG_SIZE, MPI_INT, root, MPI_COMM_WORLD)
                != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_Gather"); exit(0);
            }

            if (myrank == root)
            {
                for (k = 0; k < N; k++)
                {
                    unique_array(k, N, testnr, MSG_SIZE, buf_expt);
                    if (!check_arrays(&buf[k*MSG_SIZE], buf_expt, MSG_SIZE))
                    {
                        print_match_err(myrank, &buf[k*MSG_SIZE], buf_expt,
                                        MSG_SIZE);
                    }
                    else
                    {
                        /*printf("OK "); fflush(stdout);*/
                    }
                }
            }
            testnr++;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    free(buf);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}

void gatherv_data_tester(void)
{
    int  myrank, N;
    int  i, j, k;
    int  testnr = 0;
    int *buf_recv;
    int  buf_send[2*MSG_SIZE];
    int  buf_expt[2*MSG_SIZE];
    int *displs;
    int *counts;

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Gatherv data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    buf_recv = (int*)malloc(sizeof(int)*2*MSG_SIZE*N);
    displs   = (int*)malloc(sizeof(int)*N);
    counts   = (int*)malloc(sizeof(int)*N);

    for (i = 0; i < GATHERV_NTRIES; i++)
    {
        for (j = 0; j < GATHERV_TRYSIZE; j++)
        {
            int root = (testnr/3)%N; /* we alter root every third test */
            if (N < 2)
            {
                root = 0;    /* because otherwise we'll get stuck */

            }
            unique_array(myrank, N, testnr+myrank, msgsize(myrank), buf_send);

            if (myrank == root)
            {
                int c = 0;
                for (k = 0; k < N; k++)
                {
                    displs[k] = c;
                    counts[k] = msgsize(k);
                    c        += counts[k];
                }
            }

            if (MPI_Gatherv(buf_send, msgsize(myrank), MPI_INT,
                            buf_recv, counts, displs, MPI_INT, root,
                            MPI_COMM_WORLD)
                != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_Gatherv"); exit(0);
            }

            if (myrank == root)
            {
                for (k = 0; k < N; k++)
                {
                    unique_array(k, N, testnr+k, msgsize(k), buf_expt);
                    if (!check_arrays(buf_recv+displs[k], buf_expt, counts[k]))
                    {
                        print_match_err(myrank, buf_recv + displs[k],
                                        buf_expt, counts[k]);
                    }
                }
            }
            testnr += N;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    free(buf_recv);
    free(displs);
    free(counts);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}




void scatter_data_tester(void)
{
    int  myrank, N;
    int  i, j, k;
    int  testnr = 0;
    int *buf;
    int  buf_recv[MSG_SIZE];
    int  buf_expt[MSG_SIZE];

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Scatter data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    buf = (int*)malloc(sizeof(int)*MSG_SIZE*N);
    for (i = 0; i < SCATTER_NTRIES; i++)
    {
        for (j = 0; j < SCATTER_TRYSIZE; j++)
        {
            int root = (testnr/3)%N; /* we alter root every third test */
            if (N < 2)
            {
                root = 0;    /* because otherwise we'll get stuck */

            }
            if (myrank == root)
            {
                /* prepare data */
                for (k = 0; k < N; k++)
                {
                    unique_array(myrank, N, testnr+k, MSG_SIZE, buf+k*MSG_SIZE);
                }

                if (MPI_Scatter(buf, MSG_SIZE, MPI_INT,
                                buf_recv, MSG_SIZE, MPI_INT,
                                root, MPI_COMM_WORLD)
                    != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Scatter"); exit(0);
                }
            }
            else
            {
                if (MPI_Scatter(buf, MSG_SIZE, MPI_INT,
                                buf_recv, MSG_SIZE, MPI_INT,
                                root, MPI_COMM_WORLD)
                    != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Scatter"); exit(0);
                }
            }

            unique_array(root, N, testnr+myrank, MSG_SIZE, buf_expt);
            if (!check_arrays(buf_recv, buf_expt, MSG_SIZE))
            {
                print_match_err(myrank, buf_recv, buf_expt, MSG_SIZE);
            }

            testnr += N;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    free(buf);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}

void scatterv_data_tester(void)
{
    int  myrank, N;
    int  i, j, k;
    int  testnr = 0;
    int *sendbuf;
    int *displs;
    int *counts;
    int  buf_recv[2*MSG_SIZE];
    int  buf_expt[2*MSG_SIZE];

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Scatterv data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    sendbuf = (int*)malloc(sizeof(int)*2*MSG_SIZE*N);
    displs  = (int*)malloc(sizeof(int)*N);
    counts  = (int*)malloc(sizeof(int)*N);
    if (MSG_SIZE < N+2)
    {
        fprintf(stderr, "ERROR: N > MSG_SIZE; can't do scatterv test\n");
        return;
    }

    for (i = 0; i < SCATTERV_NTRIES; i++)
    {
        for (j = 0; j < SCATTERV_TRYSIZE; j++)
        {
            int root = (testnr/3)%N; /* we alter root every third test */
            if (N < 2)
            {
                root = 0;    /* because otherwise we'll get stuck */

            }
            if (myrank == root)
            {
                int dcount = 0;
                /* prepare data */
                for (k = 0; k < N; k++)
                {
                    displs[k] = dcount;
                    counts[k] = msgsize(k); /*k+MSG_SIZE/2;*/
                    unique_array(myrank, N, testnr+k, counts[k],
                                 sendbuf+displs[k]);
                    dcount += counts[k];
                }

                if (MPI_Scatterv(sendbuf, counts, displs, MPI_INT,
                                 buf_recv, msgsize(myrank), MPI_INT,
                                 root, MPI_COMM_WORLD)
                    != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Scatterv"); exit(0);
                }
            }
            else
            {
                if (MPI_Scatterv(NULL, NULL, NULL, MPI_INT,
                                 buf_recv, msgsize(myrank), MPI_INT,
                                 root, MPI_COMM_WORLD)
                    != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR in MPI_Scatterv"); exit(0);
                }
            }

            unique_array(root, N, testnr+myrank, msgsize(myrank), buf_expt);
            if (!check_arrays(buf_recv, buf_expt, msgsize(myrank)))
            {
                print_match_err(myrank, buf_recv, buf_expt, msgsize(myrank));
            }
            testnr += N;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    free(sendbuf);
    free(displs);
    free(counts);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}



void alltoallv_data_tester(void)
{
    int  myrank, N;
    int  i, j, k;
    int  testnr = 0;
    int *buf_recv;
    int *buf_send;
    int  buf_expt[2*MSG_SIZE];
    int *rdispls, *rcounts;
    int *sdispls, *scounts;

    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Alltoallv data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    buf_send = (int*)malloc(sizeof(int)*2*MSG_SIZE*N);
    buf_recv = (int*)malloc(sizeof(int)*2*MSG_SIZE*N);
    sdispls  = (int*)malloc(sizeof(int)*N);
    scounts  = (int*)malloc(sizeof(int)*N);
    rdispls  = (int*)malloc(sizeof(int)*N);
    rcounts  = (int*)malloc(sizeof(int)*N);

    for (i = 0; i < ALLTOALLV_NTRIES; i++)
    {
        for (j = 0; j < ALLTOALLV_TRYSIZE; j++)
        {
            int rcount = 0;
            int scount = 0;

            /* prepare data */
            for (k = 0; k < N; k++)
            {
                /* send data */
                sdispls[k] = scount;
                scounts[k] = msgsize(myrank);
                unique_array(myrank, N, testnr+myrank+k, scounts[k],
                             buf_send+sdispls[k]);
                scount += scounts[k];
                /* receive size */
                rdispls[k] = rcount;
                rcounts[k] = msgsize(k);
                rcount    += rcounts[k];
            }

            if (MPI_Alltoallv(buf_send, scounts, sdispls, MPI_INT,
                              buf_recv, rcounts, rdispls, MPI_INT,
                              MPI_COMM_WORLD)
                != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_Alltoallv"); exit(0);
            }

            for (k = 0; k < N; k++)
            {
                unique_array(k, N, testnr+k+myrank, msgsize(k), buf_expt);

                if (!check_arrays(buf_recv + rdispls[k], buf_expt, rcounts[k]))
                {
                    print_match_err(myrank, buf_recv + rdispls[k], buf_expt,
                                    rcounts[k]);
                }
            }

            testnr += N;
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    free(buf_recv);
    free(rdispls);
    free(rcounts);
    free(buf_send);
    free(sdispls);
    free(scounts);
    if (myrank == 0)
    {
        printf("\nOK.\n");
        fflush(stdout);
    }
}
