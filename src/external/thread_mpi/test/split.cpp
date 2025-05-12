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


#include "tmpi.h"

#include "ntests.h"


void split_tester(int N, int myrank)
{
    int i, j, k;

    if (myrank == 0)
    {
        printf("\nStarting MPI_Split test\n");
    }
    {
        int      newrank, newsize;
        MPI_Comm comm;

        MPI_Comm_split(MPI_COMM_WORLD, myrank%2, -myrank, &comm);
        MPI_Comm_rank(comm, &newrank);
        MPI_Comm_size(comm, &newsize);
        printf("I'm in world rank %d, and have new rank %d in group %p of size %d\n",
               myrank, newrank, ((void *)comm), newsize);
        MPI_Comm_free(&comm);
    }

    /*MPI_Barrier(MPI_COMM_WORLD);*/
    if (myrank == 0)
    {
        printf("\nStarting MPI_Split test\n");
    }
    {
        int      newrank, newsize;
        MPI_Comm comm;

        MPI_Comm_split(MPI_COMM_WORLD, myrank%3, -myrank, &comm);
        MPI_Comm_rank(comm, &newrank);
        MPI_Comm_size(comm, &newsize);
        printf("I'm in world rank %d, and have new rank %d in group %p of size %d\n", myrank, newrank, ((void *)comm), newsize);

        {
            int compare_res;

            MPI_Comm_compare(comm, MPI_COMM_WORLD, &compare_res);
            printf("my rank is %d, MPI_Comm_compare said: %d\n",
                   myrank, compare_res);
        }
        MPI_Comm_free(&comm);
    }

    if (myrank == 0)
    {
        printf("\nStarting MPI_Split test\n");
    }
    {
        int      newrank, newsize;
        MPI_Comm comm;

        MPI_Comm_split(MPI_COMM_WORLD, myrank/8, -myrank, &comm);
        MPI_Comm_rank(comm, &newrank);
        MPI_Comm_size(comm, &newsize);
        printf("I'm in world rank %d, and have new rank %d in group %p of size %d\n", myrank, newrank, ((void *)comm), newsize);

        {
            int compare_res;

            MPI_Comm_compare(comm, MPI_COMM_WORLD, &compare_res);
            printf("my rank is %d, MPI_Comm_compare said: %d\n",
                   myrank, compare_res);
        }
        MPI_Comm_free(&comm);
    }

    if (myrank == 0)
    {
        printf("Starting repeated MPI_Split() test.\n");
    }
    for (k = 0; k < NTRIES; k++)
    {
        for (i = 0; i < NSPLITTEST; i++)
        {
            int      newrank, newsize;
            MPI_Comm comm;
            int      ncolors = i%16+1;
            int      mycolor = (myrank+i) % ncolors;
            int      expected_rank, expected_size; /* the expected rank and size of
                                                      the new communicator*/

            MPI_Comm_split(MPI_COMM_WORLD, mycolor, myrank, &comm);
            MPI_Comm_rank(comm, &newrank);
            MPI_Comm_size(comm, &newsize);

            expected_size = 0;
            expected_rank = 0;

            for (j = 0; j < N; j++)
            {
                int thiscolor = (j + i) % ncolors;
                if (mycolor == thiscolor)
                {
                    expected_size += 1;
                    if (j < myrank)
                    {
                        expected_rank += 1;
                    }
                }
            }
            if (expected_size != newsize)
            {
                fprintf(stderr,
                        "ERROR: expected communicator size %d, was %d\n",
                        expected_size, newsize);
                exit(1);
            }
            if (expected_rank != newrank)
            {
                fprintf(stderr,
                        "ERROR: expected communicator rank %d, was %d\n",
                        expected_rank, newrank);
                exit(1);
            }
            MPI_Comm_free(&comm);
        }
        if (myrank == 0)
        {
            printf(">");
            fflush(stdout);
        }
    }
    if (myrank == 0)
    {
        printf("\nDone.\n");
    }


    if (N >= 2*3)
    {
        /* maybe we should test more dimensions */
#define Ndims 2
        int      dims[Ndims], periods[Ndims];
        MPI_Comm comm;

        if (myrank == 0)
        {
            printf("\nStarting Cartesian topology test\n");
        }

        for (i = 0; i < Ndims; i++)
        {
            periods[i] = 1;
        }
        dims[0] = 2; dims[1] = 3; /*dims[2]=4;*/
        MPI_Cart_create(MPI_COMM_WORLD, Ndims, dims, periods, 0, &comm);
        if (comm)
        {
            int          newrank;
            int          coord[Ndims];
            float        fcoord[Ndims];
            MPI_Datatype dt;

            MPI_Type_contiguous(Ndims, MPI_FLOAT, &dt);
            MPI_Type_commit(&dt);

            MPI_Comm_rank(comm, &newrank);
            MPI_Cart_coords(comm, newrank, Ndims, coord);

            printf("I am process %d, new rank %d, new coords (%d, %d), comm=%p\n",
                   myrank, newrank, coord[0], coord[1], (void*)comm);
            MPI_Barrier(comm);

            if (newrank != 0)
            {
                for (j = 0; j < Ndims; j++)
                {
                    fcoord[j] = (float)(coord[j]);
                }
                /* send our coords to 0 */
                MPI_Send(fcoord, 1, dt, 0, 1000, comm);
            }
            else
            {
                int Nc;
                MPI_Comm_size(comm, &Nc);
                for (j = 1; j < Nc; j++)
                {
                    MPI_Recv(fcoord, 1, dt, j, 1000, comm, NULL);

                    printf("Received from %d: (%g %g)\n",
                           j, fcoord[0], fcoord[1]);

                }
            }

            /*printf("I am process %d, dt=%d\n", myrank, (int)dt);*/

            {
                int      remain_dims[Ndims];
                int      nnewrank;
                int      ncoord[Ndims];
                MPI_Comm newcomm;

                for (i = 0; i < Ndims; i++)
                {
                    remain_dims[i] = (i == 1);
                }
                MPI_Cart_sub(comm, remain_dims, &newcomm);
                MPI_Comm_rank(newcomm, &nnewrank);
                MPI_Cart_coords(newcomm, nnewrank, Ndims-1, ncoord);

                printf("SUB: I am process %d, new rank %d, new coords (%d), comm=%p\n",
                       myrank, nnewrank, ncoord[0], (void*)newcomm);
            }
        }
        else
        {
            printf("I am process %d, and I'm not in the new communicator\n",
                   myrank);
        }
        if (comm)
        {
            MPI_Comm_free(&comm);
        }
    }
}
