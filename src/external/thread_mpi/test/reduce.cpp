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

#include "reduce.h"


#include "ntests.h"

#define MSG_SIZE 80

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif


#define REDUCE_SIZE 10


void reduce_tester(void)
{
#if 0
    int i, j, k;
    int myrank;
    int N;
    tMPI_Comm_size( TMPI_COMM_WORLD, &N );
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);

    for (i = 0; i < REDUCE_NTRIES; i++)
    {
        for (j = 0; j < REDUCE_TRYSIZE; j++)
        {
            int source[REDUCE_SIZE];
            /*int dest[REDUCE_SIZE];*/

            for (k = 0; k < REDUCE_SIZE; k++)
            {
                source[k] = myrank + N*k;
            }

            tMPI_Reduce(TMPI_IN_PLACE, source, REDUCE_SIZE, TMPI_INT,
                        TMPI_SUM, j%N, TMPI_COMM_WORLD);

            if (j%N == myrank)
            {
                for (k = 0; k < REDUCE_SIZE; k++)
                {
                    int answer = 0;
                    int l;

                    for (l = 0; l < N; l++)
                    {
                        answer += l + N*k;
                    }
                    if (answer != source[k])
                    {
                        printf(" (%d: %d == %d): %s\n", i, source[k], answer,
                               (answer == source[k]) ? "OK" : "NOT OK!!");
                        abort();
                    }
                }

            }
        }
    }
#endif
}

void allreduce_tester(void)
{
    int i, j, k;
    int myrank;
    int N;

    tMPI_Comm_size( TMPI_COMM_WORLD, &N );
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);


    printf("again: my rank = %d\n", myrank);

    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting tMPI_Allreduce test\n");
    }
    {
        int source[REDUCE_SIZE];
        /*int dest[REDUCE_SIZE];*/

        for (j = 0; j < REDUCE_SIZE; j++)
        {
            source[j] = myrank + N*j;
        }

        tMPI_Allreduce(TMPI_IN_PLACE, source, REDUCE_SIZE, TMPI_INT, TMPI_PROD,
                       TMPI_COMM_WORLD);

        for (k = 0; k < N; k++)
        {
            if (k == myrank)
            {
                int l;
                int answer = 0;
                printf("Results from mpi_allreduce process %d:\n", myrank);

                for (i = 5; i < REDUCE_SIZE; i++)
                {
                    answer = 1;

                    for (l = 0; l < N; l++)
                    {
                        answer *= l + N*i;
                    }
                    printf(" (%d: %d == %d): %s\n", i, source[i], answer,
                           (answer == source[i]) ? "OK" : "NOT OK!!");
                    if (answer != source[i])
                    {
                        abort();
                    }
                }
            }
            tMPI_Barrier(TMPI_COMM_WORLD);
        }
    }
}


static int max_int(int a, int b)
{
    return (a > b) ? a : b;
}

void allreduce_tester_fn(void)
{
    int i, j, k;
    int myrank;
    int N;

    tMPI_Comm_size( TMPI_COMM_WORLD, &N );
    /*tMPI_Init(N);*/
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);


    printf("again: my rank = %d\n", myrank);

    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting tMPI_Allreduce test\n");
    }
    {
#define REDUCE_SIZE 10
        int source[REDUCE_SIZE];
        /*int dest[REDUCE_SIZE];*/

        for (j = 0; j < REDUCE_SIZE; j++)
        {
            source[j] = myrank + N*j;
        }

        tMPI_Allreduce(TMPI_IN_PLACE, source, REDUCE_SIZE, TMPI_INT, TMPI_MAX,
                       TMPI_COMM_WORLD);

        for (k = 0; k < N; k++)
        {
            if (k == myrank)
            {
                int l;
                int answer = 0;
                printf("Results from mpi_allreduce process %d:\n", myrank);

                for (i = 5; i < REDUCE_SIZE; i++)
                {
                    answer = 1;

                    for (l = 0; l < N; l++)
                    {
                        answer = max_int(answer, l+N*i);
                    }
                    /*answer *= l + N*i;*/
                    printf(" (%d: %d == %d): %s\n", i, source[i], answer,
                           (answer == source[i]) ? "OK" : "NOT OK!!");
                    if (answer != source[i])
                    {
                        abort();
                    }
                }
            }
            tMPI_Barrier(TMPI_COMM_WORLD);
        }
    }
}


void scan_tester(void)
{
    int i, j, k;
    int myrank;
    int N;

    tMPI_Comm_size( TMPI_COMM_WORLD, &N );
    tMPI_Comm_rank( TMPI_COMM_WORLD, &myrank);

    printf("again: my rank = %d\n", myrank);

    tMPI_Barrier(TMPI_COMM_WORLD);
    if (myrank == 0)
    {
        printf("\nStarting tMPI_Scan test\n");
    }

    {
#define REDUCE_SIZE 10
        int source[REDUCE_SIZE];
        int dest[REDUCE_SIZE];

        for (j = 0; j < REDUCE_SIZE; j++)
        {
            source[j] = myrank + N*j;
        }

        tMPI_Scan(source, dest, REDUCE_SIZE, TMPI_INT, TMPI_SUM,
                  TMPI_COMM_WORLD);

        for (k = 0; k < N; k++)
        {
            if (k == myrank)
            {
                int l;
                int answer = 0;
                if (myrank == 0)
                {
                    printf("Results from mpi_scan process %d:\n", myrank);
                }

                for (i = 5; i < REDUCE_SIZE; i++)
                {
                    answer = 0;

                    for (l = 0; l < myrank+1; l++)
                    {
                        answer += l + N*i;
                    }
                    if (answer != dest[i])
                    {
                        printf("%d (%d: %d == %d): %s\n", myrank, i,
                               dest[i], answer,
                               (answer == dest[i]) ? "OK" : "NOT OK!!");
                        abort();
                    }
                }
                if (myrank == 0)
                {
                    printf("OK \n");
                }
            }
            tMPI_Barrier(TMPI_COMM_WORLD);
        }
    }
}
