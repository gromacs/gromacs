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
/*#include <unistd.h>*/


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "tmpi.h"

#define MSG_SIZE 80
/*#define N 8*/

#define MAX_THREADS 512

#if defined( _WIN32 ) || defined( _WIN64 )
#define snprintf sprintf_s
#endif

static void tester(const void* arg)
{
    char       message[MSG_SIZE];
    char       recv_message[MSG_SIZE];
    char       send_message[MSG_SIZE];
    int        i, j, k;
    MPI_Status status;
    int        myrank;
    int        left;
    int        right;
    int        N;

    int        argp = *((int*)arg);

    /* now check whether we're in a multithreaded  environment */
    {
        int flag;
        MPI_Initialized(&flag);
        if (!flag)
        {
            printf("MPI not initialized:\n");
            printf("was '-np threadnr' specified on cmd line?\n");
            /*exit(1);*/
        }
    }

    MPI_Comm_size( MPI_COMM_WORLD, &N );
    /*MPI_Init(N);*/
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank);

    printf("thread %d: got arg %d; total number of hw threads: %d\n",
           myrank, argp, tMPI_Thread_get_hw_number());
    if (N > MAX_THREADS)
    {
        printf("ERROR: %d is bigger than max # of threads: %d\n", N,
               MAX_THREADS);
        exit(1);
    }

    left  = ((myrank-1)+N)%N;
    right = (myrank+1)%N;


    printf("my rank = %d\n", myrank);


    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Sendrecv tests\n");
    }
    /* the message we're going to sendrecv */
    snprintf(send_message, MSG_SIZE, "Sendrecv hello, from thread %d", myrank);
    for (i = 0; i < N; i++)
    {
        /* we now shift everything, by sending to rank-1, and receiving
           from rank+1 */
        if (MPI_Sendrecv(send_message, (int)strlen(send_message)+1, MPI_CHAR,
                         left, i+20,
                         recv_message, MSG_SIZE, MPI_CHAR, right,
                         i+20, MPI_COMM_WORLD,
                         &status) != MPI_SUCCESS)
        {
            printf("ERROR: MPI_Sendrecv error");
            fflush(0);
        }
        if (status.MPI_ERROR != MPI_SUCCESS)
        {
            printf("ERROR: got receive status error on rank %d", myrank);
            fflush(0);
        }
        if (myrank == 0)
        {
            printf("Received: '%s'\n", recv_message);
        }
        /* now send what we received in the next step */
        memcpy(send_message, recv_message, MSG_SIZE);
        /*sleep(1);*/
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Bcast test\n");
    }
    {
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < N; j++)
            {
                snprintf(message, MSG_SIZE, "Broadcast message %d from %d", i,
                         myrank);
                MPI_Bcast(message, MSG_SIZE, MPI_CHAR, j, MPI_COMM_WORLD);
                if (myrank == 0)
                {
                    printf("Received: '%s'\n", message);
                }
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Gather test\n");
    }
    for (k = 0; k < 2; k++)
    {
        for (j = 0; j < 2; j++)
        {
            char msg_recv[MAX_THREADS*MSG_SIZE];
            int  msglen;
            int  root = k%N;
            if (N < 2)
            {
                root = 0;
            }

            for (i = 0; i < N*MSG_SIZE; i++)
            {
                msg_recv[i] = 0;
            }
            snprintf(message, MSG_SIZE, "Gather message %02d from %03d\n", j,
                     myrank);
            msglen = (int)strlen(message);
            if (MPI_Gather(message, msglen, MPI_CHAR,
                           msg_recv, msglen, MPI_CHAR,
                           root, MPI_COMM_WORLD) !=
                MPI_SUCCESS)
            {
                printf("MPI_Gather returned error\n"); fflush(0);
            }

            if (myrank == k)
            {
                printf("Received: '%s'\n", msg_recv);
            }
        }
    }


    if (myrank == 0)
    {
        printf("Finishing..\n");
    }

    /*MPI_Finalize();*/

}


int main(int argc, char *argv[])
{
    int arg = 10;
    int n;

    /*if (tMPI_Get_N(&argc, &argv, NULL, &n) != MPI_SUCCESS)*/
    if (tMPI_Get_N(&argc, &argv, "-nt", &n) != MPI_SUCCESS)
    {
        fprintf(stderr, "alloc_check test program.\n");
        fprintf(stderr, "Usage: alloc_check -nt <nthreads>\n");
        exit(0);
    }
    else
    {
        printf("\nMemory (de)allocation tester.\n\n");
        printf("Number of threads: %d\n\n", n);
    }
    tMPI_Init_fn(1, n, TMPI_AFFINITY_ALL_CORES, tester, &arg);
    /* only the main thread goes here: we've said main_thread_returns=false,
       so we can safely run the tester. */
    tester(&arg);
    tMPI_Finalize();
    return 0;
}
