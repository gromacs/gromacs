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


/* request list: */
/* get a request from the thread's pre-allocated request list */
struct tmpi_req_ *tMPI_Get_req(struct req_list *rl);
/* return a request to the thread's pre-allocated request list */
void tMPI_Return_req(struct req_list *rl, struct tmpi_req_ *req);

/* initialize a request with sensible values */
void tMPI_Req_init(struct tmpi_req_ *rq, struct envelope *ev);

/* wait for incoming connections (and the completion of outgoing connections
   if spin locks are disabled), and handle them. */
void tMPI_Wait_process_incoming(struct tmpi_thread *th);





/* check for the completion of a single request */
tmpi_bool tMPI_Test_single(struct tmpi_thread *cur,
                           struct tmpi_req_   *rq);
/* check and wait for the completion of a single request */
void tMPI_Wait_single(struct tmpi_thread *cur, struct tmpi_req_ *rq);

/* check for the completion of a NULL-delimited doubly linked list of
   requests */
tmpi_bool tMPI_Test_multi(struct tmpi_thread *cur, struct tmpi_req_ *rqs,
                          tmpi_bool *any_done);



/* set a request status */
void tMPI_Set_status(struct tmpi_req_ *req, tMPI_Status *st);


/* post a send envelope */
struct envelope *tMPI_Post_send(struct tmpi_thread *cur,
                                tMPI_Comm comm,
                                struct tmpi_thread *dest,
                                void *send_buf, int send_count,
                                tMPI_Datatype datatype, int tag,
                                tmpi_bool nonblock);

/* post and match a receive envelope */
struct envelope* tMPI_Post_match_recv(struct tmpi_thread *cur,
                                      tMPI_Comm comm,
                                      struct tmpi_thread *src,
                                      void *recv_buf, int recv_count,
                                      tMPI_Datatype datatype,
                                      int tag, tmpi_bool nonblock);
