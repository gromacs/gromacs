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

To help us fund development, we humbly ask that you cite
any papers on the package - you can find them in the top README file.

*/


/* this file is to be #included from collective.c */

/* broadcast */
int tMPI_Bcast(void* buffer, int count, tMPI_Datatype datatype, int root,
               tMPI_Comm comm)
{
    int synct;
    struct coll_env *cev;
    int myrank;
    int ret=TMPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    /* we increase our counter, and determine which coll_env we get */
    cev=tMPI_Get_cev(comm, myrank, &synct);

    if (myrank==root)
    {
        /* first set up the data */
        tMPI_Post_multi(cev, myrank, 0, TMPI_BCAST_TAG, datatype, 
                        count*datatype->size, buffer, comm->grp.N-1, synct);
        /* and wait until everybody is done copying */
        tMPI_Wait_for_others(cev, myrank);
    }
    else
    {
        size_t bufsize=count*datatype->size;
        /* wait until root becomes available */
        tMPI_Wait_for_data(cev, root, synct);
        tMPI_Mult_recv(comm, cev, root, 0, TMPI_BCAST_TAG, datatype, bufsize, 
                       buffer, &ret);
    }
    return ret;
}





