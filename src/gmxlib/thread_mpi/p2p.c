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

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


#include "impl.h"



/* free envelopes: */
static struct envelope *tMPI_Free_env_list_fetch_recv(struct free_envelope_list 
                                                      *evl);

/* return an envelope to the free envelopes list */
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct envelope *rev);



/* send envelopes: */
/* get a new envelope from the send list's free envelope list */
static struct envelope*
tMPI_Send_env_list_fetch_new(struct send_envelope_list *evl);

/* return a send envelope to the send list's free envelope list, 
    (to be used by the sending thread, who owns the send_envelope_list) */
static void tMPI_Send_env_list_return(struct envelope *ev);
#ifdef USE_SEND_RECV_COPY_BUFFER
/* return a send envelope to the sender's send list. 
    (to be used by the receiving thread). */
static void tMPI_Send_env_list_rts(struct envelope *sev);
#endif

/* remove a send envelope from the old list. Does not lock */
static void tMPI_Send_env_list_remove_old(struct envelope *sev);



/* remove a send envelope from its head_old list. Does not lock */
static void tMPI_Send_env_list_remove_old(struct envelope *sev);

/* add a send envelope to the new envelopes queue in a list */
static void tMPI_Send_env_list_add_new(struct tmpi_thread *cur, 
                                       struct send_envelope_list *evl,
                                       struct envelope *sev);
/* move a send envelope to the old envelopes queue in a list. 
   Assumes that this is safe to do without interference
   from other threads, i.e. the list it's in must have been
   detached. */
static void tMPI_Send_env_list_move_to_old(struct envelope *sev);


/* receive envelopes: */
/* add a receive envelope to a list */
static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl,
                                   struct envelope *ev);
/* remove a receive envelope from its list */
static void tMPI_Recv_env_list_remove(struct envelope *ev);




/* request list: */
/* get a request from the thread's pre-allocated request list */
static struct tmpi_req_ *tMPI_Get_req(struct req_list *rl);
/* return a request to the thread's pre-allocated request list */
static void tMPI_Return_req(struct req_list *rl, struct tmpi_req_ *req);

/* initialize a request with sensible values */
static void tMPI_Req_init(struct tmpi_req_ *rq, struct envelope *ev);

/* wait for incoming connections (and the completion of outgoing connections
   if spin locks are disabled), and handle them. */
static void tMPI_Wait_process_incoming(struct tmpi_thread *th);

/* do the actual point-to-point transfer */
static void tMPI_Xfer(struct tmpi_thread *cur, struct envelope *sev, 
                      struct envelope *rev);


/* check for the completion of a single request */
static bool tMPI_Test_single(struct tmpi_thread *cur, struct tmpi_req_ *rq);
/* check and wait for the completion of a single request */
static void tMPI_Wait_single(struct tmpi_thread *cur, struct tmpi_req_ *rq);

/* check for the completion of a NULL-delimited doubly linked list of 
   requests */
static bool tMPI_Test_multi(struct tmpi_thread *cur, struct tmpi_req_ *rqs,
                            bool *any_done);



#include "p2p_protocol.c"
#include "p2p_send_recv.c"
#include "p2p_wait.c"

