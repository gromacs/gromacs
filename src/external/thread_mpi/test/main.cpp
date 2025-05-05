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
#include <errno.h>


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "tmpi.h"

#include "cond-var.h"
#include "atomic.h"
#include "reduce.h"
#include "once.h"
#include "send_recv.h"
#include "multicast.h"
#include "cycles.h"
#include "split.h"
#include "ntests.h"

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

#ifdef TEST_AFFINITY
    if (tMPI_Thread_setaffinity_support() == TMPI_SETAFFINITY_SUPPORT_YES)
    {
        int ret;

        if (myrank == 0)
        {
            printf("Setting thread affinity\n");
        }
        ret = tMPI_Thread_setaffinity_single(tMPI_Thread_self(), myrank);
        if (ret != 0)
        {
            printf("ERROR setting thread affinity: %d (%s)", ret,
                   strerror(errno));
        }
    }
#endif



#ifdef TEST_ATOMICS
    atomic_fn_tester();
#endif

#if 0
#ifdef CYCLE_COUNT
    cycles_tester();
#endif

#endif
#ifdef TEST_COND_VARS
    cond_var_tester();
#endif
#ifdef TEST_ATOMICS
    atomic_tester();
#endif

#ifdef TEST_BARRIERS
    barrier_data_tester();
#endif

    /*sleep(100);*/

    MPI_Barrier(MPI_COMM_WORLD);
#ifdef TEST_P2P
    if (myrank == 0)
    {
        printf("Starting MPI_Send and MPI_Recv tests\n");
    }
    for (i = 0; i < 2; i++)
    {
        if (myrank != 0)
        {
            snprintf(message, MSG_SIZE, "From thread %d to %d: hello.",
                     myrank, i);
            /*strcpy(message, "Hello, there\n");*/
            if (MPI_Send(message, (int)strlen(message)+1, MPI_CHAR,
                         0, i,
                         MPI_COMM_WORLD) != MPI_SUCCESS)
            {
                printf("ERROR: MPI_Send error");
                fflush(0);
            }
        }
        else
        {
            /*MPI_Recv(message, 20, MPI_CHAR, 0, 99, MPI_COMM_WORLD, &status);*/

            for (j = 0; j < N-1; j++)
            {
#if 1
                if (MPI_Recv(message, MSG_SIZE, MPI_CHAR, MPI_ANY_SOURCE,
                             i, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
#else
                if (MPI_Recv(message, MSG_SIZE, MPI_CHAR, MPI_ANY_SOURCE,
                             i, NULL, &status) != MPI_SUCCESS)
#endif
                {
                    printf("ERROR: MPI_Recv error");
                    fflush(0);
                }
                if (status.MPI_ERROR != MPI_SUCCESS)
                {
                    printf("ERROR: got receive status error on rank %d",
                           myrank);
                    fflush(0);
                }
                printf("Received: '%s'\n", message);
            }
        }
    }


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

    if (myrank == 0)
    {
        printf("\nStarting MPI_Isend/MPI_Irecv tests\n");
    }
#define PRINTMSG 10
    k = 0;
    for (j = 0; j < SEND_RECV_INITIAL_NTRIES*SEND_RECV_INITIAL_TRYSIZE; j++)
    {
        const int   recv_left  = 0;
        const int   recv_right = 1;
        const int   send_left  = 2;
        const int   send_right = 3;

        char        buf[4][MSG_SIZE];
        MPI_Request reqs[4];


        if (j%(SEND_RECV_INITIAL_TRYSIZE) == 0)
        {
            snprintf(buf[send_left], MSG_SIZE,  "Msg %d to  left (%d) from %d",
                     k, left, myrank);
            snprintf(buf[send_right], MSG_SIZE, "Msg %d to right (%d) from %d",
                     k, right, myrank);
            k++;
        }

        for (i = 0; i < N*2; i++)
        {
            MPI_Isend(buf[send_left],
                      (int)strlen(buf[send_left])+1, MPI_CHAR,
                      left, i+100, MPI_COMM_WORLD, &(reqs[send_left]) );
            MPI_Irecv(buf[recv_left], MSG_SIZE, MPI_CHAR,
                      left, i+100, MPI_COMM_WORLD, &(reqs[recv_left]) );

            MPI_Irecv(buf[recv_right], MSG_SIZE, MPI_CHAR,
                      right, i+100, MPI_COMM_WORLD, &(reqs[recv_right]) );
            MPI_Isend(buf[send_right],
                      (int)strlen(buf[send_right])+1, MPI_CHAR,
                      right, i+100, MPI_COMM_WORLD, &(reqs[send_right]) );

            /*MPI_Wait(&(reqs[recv_left]), 0);
               MPI_Wait(&(reqs[recv_right]), 0);
               MPI_Wait(&(reqs[send_left]), 0);
               MPI_Wait(&(reqs[send_right]), 0);*/
            MPI_Waitall(4, reqs, 0);

            if (myrank == 0 && j%(SEND_RECV_INITIAL_TRYSIZE) == 0)
            {
                printf("From  left: '%s', ", buf[recv_left]);
                printf("From right: '%s'\n", buf[recv_right]);
            }
            memcpy(buf[send_right], buf[recv_left], MSG_SIZE);
            memcpy(buf[send_left], buf[recv_right], MSG_SIZE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    send_recv_data_tester();

    MPI_Barrier(MPI_COMM_WORLD);
#endif


#ifdef TEST_SPLIT
    split_tester(N, myrank);
#endif


#ifdef TEST_MULTICAST
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

    bcast_data_tester();


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


    gather_data_tester();
    gatherv_data_tester();

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Scatter test\n");
    }
    {
        char msg_send[MAX_THREADS*MSG_SIZE];
        for (i = 0; i < N*MSG_SIZE; i++)
        {
            msg_send[i] = 0;
        }
        if (myrank == 0)
        {
            /* create scatter message */
            for (i = 0; i < N; i++)
            {
                snprintf(msg_send+i*MSG_SIZE, MSG_SIZE,
                         "Scatter message from %03d to %03d",
                         myrank, i);
            }
        }

        if (MPI_Scatter(msg_send, MSG_SIZE, MPI_CHAR,
                        message, MSG_SIZE, MPI_CHAR,
                        0, MPI_COMM_WORLD) !=
            MPI_SUCCESS)
        {
            printf("MPI_Scatter returned error\n"); fflush(0);
        }

        /* now receive messages back */
        if (myrank != 0)
        {
            MPI_Send(message, MSG_SIZE, MPI_CHAR, 0, 200, MPI_COMM_WORLD);
        }
        else
        {
            printf("Received from %d: '%s'\n", myrank, message);
            for (i = 0; i < N; i++)
            {
                if (i != myrank)
                {
                    MPI_Recv(message, MSG_SIZE, MPI_CHAR, i, 200,
                             MPI_COMM_WORLD, 0);
                    printf("Received back from %d: '%s'\n", i, message);
                }
            }
        }
    }

    scatter_data_tester();
    scatterv_data_tester();

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Alltoall test\n");
    }
    for (k = 0; k < 4; k++)
    {
        char msg_send[MAX_THREADS*MSG_SIZE];
        char msg_recv[MAX_THREADS*MSG_SIZE];
        for (i = 0; i < N*MSG_SIZE; i++)
        {
            msg_send[i] = 0;
        }
        /* create alltoall message */
        for (i = 0; i < N; i++)
        {
            snprintf(msg_send+i*MSG_SIZE, MSG_SIZE, "%d->%d(%d)",
                     myrank, i, k);
        }

        if (MPI_Alltoall(msg_send, MSG_SIZE, MPI_CHAR,
                         msg_recv, MSG_SIZE, MPI_CHAR,
                         MPI_COMM_WORLD) !=
            MPI_SUCCESS)
        {
            printf("MPI_Alltoall returned error\n"); fflush(0);
        }

        /* now receive messages back */
        if (myrank != 0)
        {
            MPI_Send(msg_recv, N*MSG_SIZE, MPI_CHAR, 0, 300, MPI_COMM_WORLD);
        }
        else
        {
            printf("From %d: ", myrank);
            for (j = 0; j < N; j++)
            {
                printf("%s ", msg_recv+j*MSG_SIZE);
            }
            printf("\n");

            for (i = 0; i < N; i++)
            {
                if (i != myrank)
                {
                    MPI_Recv(msg_recv, N*MSG_SIZE, MPI_CHAR, i, 300,
                             MPI_COMM_WORLD, 0);

                    printf("From %d: ", i);
                    for (j = 0; j < N; j++)
                    {
                        printf("%s ", msg_recv+j*MSG_SIZE);
                    }
                    printf("\n");
                }
            }
        }
    }

    alltoallv_data_tester();
#endif

#ifdef TEST_REDUCE

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting MPI_Reduce test\n");
    }
    {
#define REDUCE_SIZE 10
        int source[REDUCE_SIZE];
        int dest[REDUCE_SIZE];

        for (j = 0; j < REDUCE_SIZE; j++)
        {
            source[j] = myrank + N*j;
        }

        MPI_Reduce(source, dest, REDUCE_SIZE, MPI_INT, MPI_MAX, 1,
                   MPI_COMM_WORLD);

        if (myrank == 1)
        {
            int l;
            int answer = 0;
            printf("Results from mpi_reduce:\n");

            for (i = 0; i < REDUCE_SIZE; i++)
            {
                answer = 0;

                for (l = 0; l < N; l++)
                {
                    answer = (answer > l+N*i) ? answer : l+N*i;
                }
                /*answer += l + N*i;*/
                printf(" (%d: %d == %d): %s\n", i, dest[i], answer,
                       (answer == dest[i]) ? "OK" : "NOT OK!!");
                if (answer != dest[i])
                {
                    abort();
                }
            }
        }
    }

    allreduce_tester();
    allreduce_tester_fn();
    scan_tester();
    reduce_tester();
#endif

    once_tester();

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
        fprintf(stderr, "thread_mpi test program.\n");
        fprintf(stderr, "Usage: mpithreads -nt <nthreads>\n");
        exit(0);
    }
    else
    {
        printf("\nthread_mpi tester.\n\n");
        printf("Number of threads: %d\n\n", n);
    }
    tMPI_Init_fn(1, n, TMPI_AFFINITY_ALL_CORES, tester, &arg);
    /* only the main thread goes here: we've said main_thread_returns=false,
       so we can safely run the tester. */
    tester(&arg);
    tMPI_Finalize();
    return 0;
}
