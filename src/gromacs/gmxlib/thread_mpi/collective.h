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


/* get a pointer the next coll_env once it's ready */
struct coll_env *tMPI_Get_cev(tMPI_Comm comm, int myrank, int *synct);

/* post the availability of data in a cev.
   cev         = the collective comm environment
   myrank      = my rank
   index       = the buffer index
   tag         = the tag
   datatype    = the datatype
   busize      = the buffer size
   buf         = the buffer to xfer
   n_remaining = the number of remaining threads that need to transfer
   synct       = the multicast sync number
   dest        = -1 for all theads, or a specific rank number.
 */
int tMPI_Post_multi(struct coll_env *cev, int myrank, int index,
                    int tag, tMPI_Datatype datatype,
                    size_t bufsize, void *buf, int n_remaining,
                    int synct, int dest);

/* transfer data from cev->met[rank] to recvbuf */
void tMPI_Mult_recv(tMPI_Comm comm, struct coll_env *cev, int rank,
                    int index, int expected_tag, tMPI_Datatype recvtype,
                    size_t recvsize, void *recvbuf, int *ret);

/* do a root transfer (from root send buffer to root recv buffer) */
void tMPI_Coll_root_xfer(tMPI_Comm comm,
                         tMPI_Datatype sendtype, tMPI_Datatype recvtype,
                         size_t sendsize, size_t recvsize,
                         void* sendbuf, void* recvbuf, int *ret);

/* wait for other processes to copy data from my cev */
void tMPI_Wait_for_others(struct coll_env *cev, int myrank);
/* wait for data to become available from a specific rank */
void tMPI_Wait_for_data(struct tmpi_thread *cur, struct coll_env *cev,
                        int myrank);
/*int rank, int myrank, int synct);*/
