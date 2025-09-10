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

#include "thread_mpi/tmpi.h"
#include "thread_mpi/collective.h"

#include "once.h"

#define MSG_SIZE 80

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif


static void once_fn(void *arg)
{
    int nr = *((int*)arg);
    int myrank;
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);

    printf("Printing this number once: %d (my rank=%d)\n", nr, myrank);
}
#if 0
static void* once_fn_wait(void *arg)
{
    int nr = *((int*)arg);
    int myrank;
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);

    printf("Returning this number once: %d (my rank=%d)\n", nr, myrank);
    return (void*)(size_t)myrank;
}
#endif

void once_tester(void)
{
    int i;
    int myrank;
    int N;

    tMPI_Comm_size( TMPI_COMM_WORLD, &N );
    /*tMPI_Init(N);*/
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);


    printf("again: my rank = %d\n", myrank);

    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting tMPI_Once test\n");
    }
    {
        for (i = 0; i < 10; i++)
        {
            int arg       = i;
            int was_first = 0;
            tMPI_Once(TMPI_COMM_WORLD, once_fn, &arg, &was_first);
            if ((i%10 == 0) && was_first)
            {
                printf("I was first for %d, my rank=%d\n", arg, myrank);
            }
        }
    }
    if (myrank == 0)
    {
        printf("Done.\n\n");
    }
    tMPI_Barrier(TMPI_COMM_WORLD);
#if 0
    if (myrank == 0)
    {
        printf("\nStarting tMPI_Once_wait test\n");
    }
    {
        for (i = 0; i < 10; i++)
        {
            void *ret;
            int   arg       = i;
            int   was_first = 0;

            ret = tMPI_Once_wait(TMPI_COMM_WORLD, once_fn_wait, &arg, &was_first);
            printf("i=%d, my rank=%d, once returned: %d\n", i, myrank,
                   (int)((size_t)(ret)));
        }
    }
    if (myrank == 0)
    {
        printf("Done.\n\n");
    }
#endif
}
