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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "thread_mpi.h"
#include "cond-var.h"

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif

static void signal_tester(void);
static void bcast_tester(void);

void cond_var_tester(void)
{
    signal_tester();
    bcast_tester();
}

tMPI_Thread_cond_t  cond_recv = TMPI_THREAD_COND_INITIALIZER;
tMPI_Thread_cond_t  cond_send = TMPI_THREAD_COND_INITIALIZER;
tMPI_Thread_mutex_t cdt_mtx   = TMPI_THREAD_MUTEX_INITIALIZER;
volatile int        test_receivers;
volatile int        test_send_status;
volatile int        test_count;
volatile int        test_add;
volatile int       *test_processed;

#define NSIGNAL_TESTS 1000
static void signal_tester(void)
{
    int i, j, k;
    int N, myrank;
    int N_receivers;

    tMPI_Comm_size(TMPI_COMM_WORLD, &N);
    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    N_receivers = N-1;

    /*printf("my rank  = %d\n",myrank);*/

    if (myrank == 0)
    {
#if 0
        tMPI_Thread_mutex_init(&cdt_mtx);
        tMPI_Thread_cond_init(&cond_recv);
        tMPI_Thread_cond_init(&cond_send);
#endif
        test_send_status = 0;
        test_receivers   = 0;
        test_count       = 0;
        test_add         = 0;
        test_processed   = (int*)malloc(sizeof(int)*N);
        for (k = 0; k < N; k++)
        {
            test_processed[k] = 0;
        }
        printf("Starting tMPI_Thread_cond_signal test:\n");
    }

    for (i = 0; i < N; i++)
    {
        tMPI_Barrier(TMPI_COMM_WORLD);
        /*printf("i=%d, N=%d, myrank=%d\n", i, N, myrank);*/
        if (i%N == myrank)
        {
            for (j = 0; j < NSIGNAL_TESTS; j++)
            {
                /* i am the signaler */
                tMPI_Thread_mutex_lock(&cdt_mtx);
                /* wait until the first thread is ready */
                while (test_receivers < N_receivers)
                {
                    tMPI_Thread_cond_wait(&cond_send, &cdt_mtx);
                }
                test_send_status = 1;
                test_add = test_add + 1;
                tMPI_Thread_cond_signal(&cond_recv);
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
        }
        else
        {
            int end   = 0;
            int first = 1;
            do
            {
                /* i wait for a signal */
                tMPI_Thread_mutex_lock(&cdt_mtx);
                if (first)
                {
                    test_receivers = test_receivers + 1;
                    first = 0;
                }
                tMPI_Thread_cond_signal(&cond_send);
                while (test_send_status == 0 && test_add == 0)
                {
                    tMPI_Thread_cond_wait(&cond_recv, &cdt_mtx);
                }
                test_count             = test_count + test_add;
                test_processed[myrank] = test_processed[myrank] + test_add;
                test_add                = 0;
                if (test_send_status < 0)
                {
                    end = 1;
                    test_receivers = test_receivers - 1;
                }
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
            while (!end);
        }

        if (i%N == myrank)
        {
            int recvrs;

            /* now send everybody a signal */
            do
            {
                tMPI_Thread_mutex_lock(&cdt_mtx);
                test_send_status = -1;
                recvrs           = test_receivers;
                if (recvrs > 0)
                {
                    tMPI_Thread_cond_broadcast(&cond_recv);
                }
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
            while (recvrs > 0);

            tMPI_Thread_mutex_lock(&cdt_mtx);
            if ( (test_count != NSIGNAL_TESTS) && (N_receivers > 0) )
            {
                printf("ERROR: signals processed: %d, expected: %d\n",
                       test_count, NSIGNAL_TESTS);
                fflush(stdout);
                exit(1);
            }
            else
            {
                printf("tMPI_Thread_cond_signal test: OK\n");
                fflush(stdout);
            }
            test_send_status = 0;
            test_count       = 0;
            test_add         = 0;
            for (k = 0; k < N; k++)
            {
                if (k != myrank)
                {
                    printf("    processed by thread %d: %d\n",
                           k, test_processed[k]);
                }
                test_processed[k] = 0;
            }
            tMPI_Thread_mutex_unlock(&cdt_mtx);
        }
    }
    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\n");
        free((void*)test_processed);
#if 0
        tMPI_Thread_mutex_destroy(&cdt_mtx);
        tMPI_Thread_cond_destroy(&cond_recv);
        tMPI_Thread_cond_destroy(&cond_send);
#endif
    }
}

static void bcast_tester(void)
{
    int i, j, k;
    int N, myrank;
    int N_receivers;

    tMPI_Comm_size(TMPI_COMM_WORLD, &N);
    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    N_receivers = N-1;

    if (myrank == 0)
    {
#if 0
        tMPI_Thread_mutex_init(&cdt_mtx);
        tMPI_Thread_cond_init(&cond_recv);
        tMPI_Thread_cond_init(&cond_send);
#endif
        test_send_status = 0;
        test_receivers   = 0;
        test_count       = 0;
        test_add         = 0;
        test_processed   = (int*)malloc(sizeof(int)*N);
        for (k = 0; k < N; k++)
        {
            test_processed[k] = 0;
        }
        printf("Starting tMPI_Thread_cond_broadcast test:\n");
    }

    for (i = 0; i < N; i++)
    {
        tMPI_Barrier(TMPI_COMM_WORLD);
        if (i%N == myrank)
        {
            for (j = 0; j < NSIGNAL_TESTS; j++)
            {
                /* i am the signaler */
                tMPI_Thread_mutex_lock(&cdt_mtx);
                /* wait until the first thread is ready */
                while (test_receivers < N_receivers)
                {
                    tMPI_Thread_cond_wait(&cond_send, &cdt_mtx);
                }
                test_send_status = 1;
                test_add         = test_add + (N-1);
                tMPI_Thread_cond_broadcast(&cond_recv);
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
        }
        else
        {
            int end   = 0;
            int first = 1;
            do
            {
                /* i wait for a signal */
                tMPI_Thread_mutex_lock(&cdt_mtx);
                if (first)
                {
                    test_receivers = test_receivers + 1;
                    first = 0;
                }
                tMPI_Thread_cond_signal(&cond_send);
                while (test_send_status == 0 && test_add == 0)
                {
                    tMPI_Thread_cond_wait(&cond_recv, &cdt_mtx);
                }
                if (test_add > 0)
                {
                    test_count = test_count + 1;
                    test_processed[myrank] = test_processed[myrank] + 1;
                    test_add = test_add - 1;
                }
                if (test_send_status < 0 && test_add == 0)
                {
                    end = 1;
                    test_receivers = test_receivers - 1;
                }
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
            while (!end);
        }

        if (i%N == myrank)
        {
            int recvrs;

            /* now send everybody a signal */
            do
            {
                tMPI_Thread_mutex_lock(&cdt_mtx);
                test_send_status = -1;
                recvrs           = (test_receivers > 0 || test_add > 0);
                if (recvrs)
                {
                    tMPI_Thread_cond_broadcast(&cond_recv);
                }
                tMPI_Thread_mutex_unlock(&cdt_mtx);
            }
            while (recvrs > 0);

            tMPI_Thread_mutex_lock(&cdt_mtx);
            if (test_count != (N-1)*NSIGNAL_TESTS)
            {
                printf("ERROR: signals processed: %d, expected: %d\n",
                       test_count, (N-1)*NSIGNAL_TESTS);
                fflush(stdout);
                exit(1);
            }
            else
            {
                printf("tMPI_Thread_cond_broadcast test: OK\n");
                fflush(stdout);
            }
            test_send_status = 0;
            test_count       = 0;
            test_add         = 0;
            for (k = 0; k < N; k++)
            {
                if (k != myrank)
                {
                    printf("    processed by thread %d: %d\n",
                           k, test_processed[k]);
                }
                test_processed[k] = 0;
            }
            tMPI_Thread_mutex_unlock(&cdt_mtx);
        }
    }
    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\n");
#if 0
        tMPI_Thread_mutex_destroy(&cdt_mtx);
        tMPI_Thread_cond_destroy(&cond_recv);
        tMPI_Thread_cond_destroy(&cond_send);
#endif
        free((void*)test_processed);
    }
}
