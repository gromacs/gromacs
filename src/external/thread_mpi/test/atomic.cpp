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
#include <math.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "thread_mpi/tmpi.h"
#include "thread_mpi/lock.h"
#include "thread_mpi/atomic.h"
#include "atomic.h"

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif

void spinlock_tester(int locktype, int n_tests);

#define N_TESTS 100000
void atomic_tester(void)
{
    spinlock_tester(0, N_TESTS);
    spinlock_tester(1, N_TESTS);
}

int            *bins      = NULL;
volatile int    curthread = -1;
int             finished;
tMPI_Spinlock_t slock = TMPI_SPINLOCK_INITIALIZER;
tMPI_Lock_t     lock;
unsigned int    sum;

void spinlock_tester(int locktype, int n_tests)
{
    int i, j;
    int N, myrank;

    tMPI_Comm_size(TMPI_COMM_WORLD, &N);
    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        bins     = (int*)malloc(sizeof(int)*N);
        finished = 0;
        if (locktype != 0)
        {
            tMPI_Lock_init(&lock);
        }
    }
    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        /* we cast to int because %zd is not platform-independent, and
           ints should be big enoug anyway */
        if (locktype == 0)
        {
            printf("\nsizeof(atomic_t) = %d, sizeof(atomic_t_ptr)=%d, sizeof(spinlock)=%d\n",
                   (int)sizeof(tMPI_Atomic_t),
                   (int)sizeof(tMPI_Atomic_ptr_t),
                   (int)sizeof(tMPI_Spinlock_t));
            printf("Testing spinlocks..\n");
        }
        else
        {
            printf("Testing locks..\n");
        }
    }
    for (i = 0; i < n_tests; i++)
    {
        int succeeded = 0;
        do
        {
            if (locktype == 0)
            {
                tMPI_Spinlock_lock(&slock);
            }
            else
            {
                tMPI_Lock_lock(&lock);
            }

            /*if (curthread!=myrank && !finished)*/
            {
                curthread = myrank;
                succeeded = 1;
                for (j = 0; j < 50; j++)
                {
                    /*sum+=sin(i*0.1+j*0.01);*/
                    sum += i *12 + j *20;
                }
                bins[myrank]++;
                if (curthread != myrank)
                {
                    printf("ERROR ERROR locks not exclusive\n");
                    exit(1);
                }
            }
            if (locktype == 0)
            {
                tMPI_Spinlock_unlock(&slock);
            }
            else
            {
                tMPI_Lock_unlock(&lock);
            }
        }
        while (!succeeded);
    }
    finished = 1;
    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("OK.\n\n");
        /* for(i=0;i<N;i++)
            printf("  bins[%d]=%d\n", i, bins[i]);*/
    }

    if (myrank == 0)
    {
        free(bins);
    }
}


static tMPI_Atomic_t attest;

void atomic_fn_tester(void)
{
    tMPI_Atomic_t at;
    int           N, myrank;
    int           expect;
    int           i;

    tMPI_Comm_size(TMPI_COMM_WORLD, &N);
    tMPI_Comm_rank(TMPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        printf("Testing atomic functions..\n");
        fflush(stdout);
    }


    tMPI_Atomic_set(&at, myrank);
    expect = tMPI_Atomic_get(&at);

    /* test fetch_add */
    for (i = 0; i < N_TESTS; i++)
    {
        int got;

        got = tMPI_Atomic_fetch_add(&at, myrank);

        if (got != expect)
        {
            printf("ERROR ERROR Atomic_fetch_add doesn't work\n");
            exit(1);
        }

        expect += myrank;
    }
    /* test add_return */
    for (i = 0; i < N_TESTS; i++)
    {
        int got;

        got     = tMPI_Atomic_add_return(&at, myrank);
        expect += myrank;

        if (got != expect)
        {
            printf("ERROR ERROR Atomic_add_return doesn't work\n");
            exit(1);
        }

    }

    tMPI_Atomic_set(&attest, 0);
#if 0
    for (i = 0; i < n_tests; i++)
    {
    }
#endif

    if (myrank == 0)
    {
        printf("OK.\n");
        fflush(stdout);
    }
}
