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
/* only include tmpi.h if we're not compiling with mpicc */
#include "tmpi.h"
#else
#include <mpi.h>
#endif

#include "send_recv.h"
#include "ntests.h"

#define MSG_SIZE 20

#if defined( _WIN32 ) ||  defined( _WIN64 )
#define snprintf sprintf_s
#endif


void unique_array(int myrank, int N, int testnr, int size, int *array)
{
    int i;
    unsigned int z = myrank + (N<<5) + (testnr<<10);

    array[0] = myrank;
    array[1] = testnr;
    for (i = 2; i < size; i++)
    {
        z        = 36969*(z&65535)+(z>>16);
        array[i] = static_cast<int>(z);
    }
}

int check_arrays(int *a, int *b, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (a[i] != b[i])
        {
            return 0;
        }
    }
    return 1;
}

void print_match_err(int myrank, int *buf_recv, int *buf_expt, int size)
{
    int  i;
#define MATCH_ERR_SIZE 20480
    char buf[MATCH_ERR_SIZE];

    fprintf(stderr, "\n%d: Array did not match!\n", myrank);
    snprintf(buf, MATCH_ERR_SIZE, "\n%d: Received: ", myrank);
    for (i = 0; i < size; i++)
    {
        snprintf(buf+strlen(buf), MATCH_ERR_SIZE, " %d", buf_recv[i]);
    }
    snprintf(buf+strlen(buf), MATCH_ERR_SIZE, "\n%d: Expected: ", myrank);
    for (i = 0; i < size; i++)
    {
        snprintf(buf+strlen(buf), MATCH_ERR_SIZE, " %d", buf_expt[i]);
    }
    fprintf(stderr, "%s\n", buf);
    abort();
}

void send_recv_data_tester(void)
{
    int         myrank, N;
    int         i, j;
    int         testnr = 0;
    int         buf[MSG_SIZE];
    int         buf_left[MSG_SIZE], buf_right[MSG_SIZE];
    int         expect_left[MSG_SIZE], expect_right[MSG_SIZE];
    int         left, right;
    MPI_Request reqs[4];
    const int   recv_left  = 0;
    const int   recv_right = 1;
    const int   send_left  = 2;
    const int   send_right = 3;


    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    left  = ((myrank-1)+N)%N;
    right =  (myrank+1)%N;

    if (myrank == 0)
    {
        printf("Async data test:\n");
        fflush(stdout);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < SEND_RECV_NTRIES; i++)
    {
        for (j = 0; j < SEND_RECV_TRYSIZE; j++)
        {
            unique_array(myrank, N, testnr, MSG_SIZE, buf);
            /* now send it left and right */
            if (MPI_Isend(buf, MSG_SIZE, MPI_INT, left, testnr,
                          MPI_COMM_WORLD, &(reqs[send_left])) != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_send"); exit(0);
            }
            if (MPI_Isend(buf, MSG_SIZE, MPI_INT, right, testnr,
                          MPI_COMM_WORLD, &(reqs[send_right])) != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_send"); exit(0);
            }
            /* and receive from left and right */
            if (MPI_Irecv(buf_left, MSG_SIZE, MPI_INT, left, testnr,
                          MPI_COMM_WORLD, &(reqs[recv_left])) != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_send"); exit(0);
            }
            if (MPI_Irecv(buf_right, MSG_SIZE, MPI_INT, right, testnr,
                          MPI_COMM_WORLD, &(reqs[recv_right])) != MPI_SUCCESS)
            {
                fprintf(stderr, "ERROR in MPI_send"); exit(0);
            }

            unique_array(left, N, testnr, MSG_SIZE, expect_left);
            unique_array(right, N, testnr, MSG_SIZE, expect_right);

            MPI_Waitall(4, reqs, 0);
            /*
               MPI_Wait(&reqs[0], 0);
               MPI_Wait(&reqs[1], 0);
               MPI_Wait(&reqs[2], 0);
               MPI_Wait(&reqs[3], 0);*/

            if (!check_arrays(expect_left, buf_left, MSG_SIZE))
            {
                print_match_err(myrank, buf_left, expect_left, MSG_SIZE);
                /*fprintf(stderr, "%d: Left array did not match!\n",myrank);
                   exit(0);*/
            }
            if (!check_arrays(expect_right, buf_right, MSG_SIZE))
            {
                print_match_err(myrank, buf_right, expect_right, MSG_SIZE);
                /*fprintf(stderr, "%d: Right array did not match!\n",myrank);
                   exit(0);*/
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
