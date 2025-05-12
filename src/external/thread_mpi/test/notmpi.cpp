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
#include <stdarg.h>
#include <vector>


#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "thread_mpi/atomic.h"
#include "thread_mpi/threads.h"
#include "thread_mpi/tmpi.h"


#if defined(THREAD_PTHREADS)
#include <pthread.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#elif defined(THREAD_WINDOWS)
#define snprintf sprintf_s
#endif



tMPI_Atomic_t at = {0};
int          *id_array;

/* Test for tMPI functions without having started thread_mpi tests, for
   'alien' threads*/


static void *thread_fn(void *arg)
{
    int           id =  (*((int *)arg));
    int           i;
    tMPI_Thread_t tid;
    int           ret;

    tid = tMPI_Thread_self();
    ret = tMPI_Thread_setaffinity_single(tid, id);
    if (ret != 0)
    {
        printf("My ID=%d, ERROR on thread_setaffinity!!\n", id);
    }

    for (i = 0; i < 1000; i++)
    {
        tMPI_Atomic_fetch_add(&at, 1);
    }
    printf("My ID=%d, i=%d, tid=%p\n", id, i, (void *)tid);

    return NULL;
}


#if defined(THREAD_WINDOWS)
static DWORD WINAPI thread_starter( LPVOID lpParam )
{
    thread_fn(lpParam);
    return 0;
}
#endif


int main(int argc, char *argv[])
{
    int n;
    int i;

    if (tMPI_Get_N(&argc, &argv, "-nt", &n) != TMPI_SUCCESS)
    {
        fprintf(stderr, "notmpi program.\n");
        fprintf(stderr, "Usage: notmpi -nt <nthreads>\n");
        exit(0);
    }
    else
    {
        printf("\nnotmpi.\n\n");
        printf("Number of threads: %d\n\n", n);
    }

    id_array = (int*)malloc(n*sizeof(int));
    for (i = 0; i < n; i++)
    {
        id_array[i] = i;
    }

    if (tMPI_Thread_setaffinity_support() == TMPI_SETAFFINITY_SUPPORT_YES)
    {
        printf("This platform supports setting thread affinity\n");
    }
    else
    {
        printf("This platform does NOT support setting thread affinity\n");
    }

#if defined(THREAD_PTHREADS)
    /* we assume POSIX threads here */
    {
        std::vector<pthread_t> th(n);

        /* create n-1 threads */
        for (i = 1; i < n; i++)
        {
            pthread_create(&th[i], NULL, thread_fn, (void*)(id_array+i));
        }
        thread_fn(id_array+0);

        for (i = 1; i < n; i++)
        {
            pthread_join(th[i], NULL);
        }
    }
#elif defined(THREAD_WINDOWS)
    /* Windows threads here */
    {
        DWORD *th;

        th = (DWORD*)malloc(sizeof(DWORD)*n);

        for (i = 1; i < n; i++)
        {
            CreateThread(NULL, 0, thread_starter, (void*)(id_array+i), 0, th+i);
        }
        thread_fn(id_array+0);

        for (i = 1; i < n; i++)
        {
            WaitForSingleObject(th+i, INFINITE);
        }
    }
#endif

    printf("Total count: %d\n", tMPI_Atomic_get(&at));
    free(id_array);
    return 0;
}
