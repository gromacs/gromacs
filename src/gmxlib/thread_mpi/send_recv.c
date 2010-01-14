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


#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "tmpi_impl.h"



/* free envelopes: */
static struct recv_envelope *tMPI_Free_env_list_fetch_recv(
                                              struct free_envelope_list *evl);

/* return an envelope to the free envelopes list */
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct recv_envelope *ev);



/* send envelopes: */
/* get a new envelope from the send list's free envelope list */
static struct send_envelope*
tMPI_Send_env_list_fetch_new(struct send_envelope_list *evl);

/* return a send envelope to the send list's free envelope list, 
    (to be used by the sending thread, who owns the send_envelope_list) */
static void tMPI_Send_env_list_return(struct send_envelope *ev);
#ifdef USE_SEND_RECV_COPY_BUFFER
/* return a send envelope to the sender's send list. 
    (to be used by the receiving thread). */
static void tMPI_Send_env_list_rts(struct send_envelope *ev);
#endif

/* remove a send envelope from the old list. Does not lock */
static void tMPI_Send_env_list_remove_old(struct send_envelope *ev);



/* remove a send envelope from its head_old list. Does not lock */
static void tMPI_Send_env_list_remove_old(struct send_envelope *ev);

/* add a send envelope to the new envelopes queue in a list */
static void tMPI_Send_env_list_add_new(struct send_envelope_list *evl,
                                       struct send_envelope *ev);
/* move a send envelope to the old envelopes queue in a list. 
   Assumes that this is safe to do without interference
   from other threads, i.e. the list it's in must have been
   detached. */
static void tMPI_Send_env_list_move_to_old(struct send_envelope *ev);


/* receive envelopes: */
/* add a receive envelope to a list */
static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl,
                                   struct recv_envelope *ev);
/* remove a receive envelope from its list */
static void tMPI_Recv_env_list_remove(struct recv_envelope *ev);




/* request list: */
/* get a request from the thread's pre-allocated request list */
static struct tmpi_req_ *tMPI_Get_req(struct req_list *rl);
/* return a request to the thread's pre-allocated request list */
static void tMPI_Return_req(struct req_list *rl, struct tmpi_req_ *req);

/* wait for incoming connections (and the completion of outgoing connections
   if spin locks are disabled), and handle them. */
static void tMPI_Test_incoming(struct tmpi_thread *th);

/* do the actual point-to-point transfer */
static void tMPI_Xfer(struct send_envelope *evs, struct recv_envelope *evr);






/* these are the internal versions of tMPI_Wait/Test, etc. that use 
   pre-allocated tmpi_reqs. */
static int tMPI_Wait_r(struct tmpi_req_ *rq, tMPI_Status *status);
/* wait all on an array of pointers to requests, that progressively get 
   assigned NULL once they're done. If may_free, the pointers also
   get returned to the free request lists */
static int tMPI_Waitall_r(int count, struct tmpi_req_ *array_of_requests[],
                          tMPI_Status *array_of_statuses, bool may_free);




/* Point-to-point communication protocol functions */


void tMPI_Free_env_list_init(struct free_envelope_list *evl, int N)
{
    int i;

    /* allocate the head element */
    evl->recv_alloc_head=(struct recv_envelope*)
                                tMPI_Malloc(sizeof(struct recv_envelope)*N );
    evl->head_recv=evl->recv_alloc_head;

    for(i=0;i<N;i++)
    {
        if (i < N-1)
        {
            evl->head_recv[i].next=&(evl->head_recv[i+1]);
        }
        else
        {
            evl->head_recv[i].next=NULL;
        }
        evl->head_recv[i].list=NULL;
    }
}

void tMPI_Free_env_list_destroy(struct free_envelope_list *evl)
{
    free(evl->recv_alloc_head);
    evl->head_recv=NULL;
    evl->recv_alloc_head=NULL;
}

static struct recv_envelope* 
tMPI_Free_env_list_fetch_recv(struct free_envelope_list *evl)
{
    struct recv_envelope *ret;
    if (! evl->head_recv )
    {
        /* TODO: make this do something better than crash */
        fprintf(stderr, "Ran out of recv envelopes!!!!\n");
        abort();
    }

    ret=evl->head_recv;
    evl->head_recv=ret->next;
    ret->next=NULL;
    ret->prev=NULL;
    /*evl->N--;*/

    return ret;
}

static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct recv_envelope *ev)
{
    ev->list=NULL;
    ev->prev=NULL;
    ev->next=evl->head_recv;
    evl->head_recv=ev;
}










/* tmpi_send_envelope_list functions */

void tMPI_Send_env_list_init(struct send_envelope_list *evl, int N)
{
    int i;
#ifndef TMPI_LOCK_FREE_LISTS
    tMPI_Spinlock_init( &(evl->lock_rts) );
    tMPI_Spinlock_init( &(evl->lock_new) );
#endif
    evl->Nalloc=N;

    evl->alloc_head=(struct send_envelope*) 
                        tMPI_Malloc(sizeof(struct send_envelope)*N );
    for(i=0;i<N;i++) 
    { 
        evl->alloc_head[i].next=(i<(N-1)) ? &(evl->alloc_head[i+1]) : NULL;
        evl->alloc_head[i].prev=NULL;
        evl->alloc_head[i].list=evl;
#ifdef USE_SEND_RECV_COPY_BUFFER
        evl->alloc_head[i].cb=(void*)tMPI_Malloc(sizeof(char)*COPY_BUFFER_SIZE);
#endif
    }
 
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(evl->head_new), NULL);
    tMPI_Atomic_ptr_set(&(evl->head_rts), NULL);
#else
    evl->head_new = NULL;
    evl->head_rts = NULL;
#endif
    evl->head_free = &(evl->alloc_head[1]);
    /* initialize the head_old circular list with dummy element */
    evl->head_old = evl->alloc_head; /* the first element is a dummy */
    evl->head_old->next = evl->head_old;
    evl->head_old->prev = evl->head_old;
}

void tMPI_Send_env_list_destroy(struct send_envelope_list *evl)
{
#ifdef USE_SEND_RECV_COPY_BUFFER
    int i;
    for(i=0;i<evl->Nalloc;i++)
    {
        free(evl->alloc_head[i].cb);
    }
#endif
    free(evl->alloc_head);
    evl->alloc_head=NULL;
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new=NULL; 
#endif
    evl->head_old=NULL; /* make it crash if used after tMPI_Finalize */
}


static struct send_envelope* 
tMPI_Send_env_list_fetch_new(struct send_envelope_list *evl)
{
    struct send_envelope *ret;

    do
    {
        /* first check whether any envelopes were returned to sender */
#ifdef TMPI_LOCK_FREE_LISTS
        if ( (ret=(struct send_envelope*)tMPI_Atomic_ptr_get(&(evl->head_rts))) )
#else
        if (evl->head_rts)
#endif
        {
            /* detach the list */
#ifdef TMPI_LOCK_FREE_LISTS
            /* we detach by swapping what we expect the pointer value to be,
               with NULL. If there were a cross-platform way to atomically 
               swap  without checking, we could do that, too. */
            while(tMPI_Atomic_ptr_cas( &(evl->head_rts), ret, NULL ) !=
                  (void*)ret)
            {
                ret=(struct send_envelope*)
                          tMPI_Atomic_ptr_get(&(evl->head_rts));
            }
#else
            tMPI_Spinlock_lock( &(evl->lock_rts) );
            ret=evl->head_rts;
            evl->head_rts=NULL;
            tMPI_Spinlock_unlock( &(evl->lock_rts) );
#endif
            /* now add the items to head_free */
            while(ret)
            {
                struct send_envelope *next=ret->next;
                ret->next=evl->head_free;
                evl->head_free=ret;
                ret=next;
            }
        }

        /* get the last free one off the list */
        ret=evl->head_free;
        if (!ret)
#ifdef USE_SEND_RECV_COPY_BUFFER
        {
            /* There are no free send envelopes, so all we can do is handle
               incoming requests until we get a free send envelope. */
            tMPI_Test_incoming(tMPI_Get_current());
        }
#else
        {
            /* If this happens, it most likely indicates a bug in the 
               calling program. We could fix the situation by waiting,
               but that would most likely lead to deadlocks - even
               more difficult to debug than this. */
            fprintf(stderr, "Ran out of send envelopes!!!!\n");
            abort();
        }
#endif
    } while(!ret);

    evl->head_free=ret->next;

    ret->next=NULL;
    ret->prev=NULL;
    ret->list=evl;

    /* and return it */
    return ret;
}

static void tMPI_Send_env_list_return(struct send_envelope *ev)
{
    struct send_envelope_list *evl=ev->list;

    ev->next=evl->head_free;
    evl->head_free=ev;
}


#ifdef USE_SEND_RECV_COPY_BUFFER
static void tMPI_Send_env_list_rts(struct send_envelope *ev)
{
    struct send_envelope_list *evl=ev->list;
#ifdef TMPI_LOCK_FREE_LISTS
    struct send_envelope *evn;

    do
    {
        evn=(struct send_envelope*)tMPI_Atomic_ptr_get(&evl->head_rts);
        ev->next=evn;
        /* the cmpxchg operation is a memory fence, so we shouldn't need
           to worry about out-of-order evaluation */
    }
    while (tMPI_Atomic_ptr_cas( &(evl->head_rts), evn, ev ) != (void*)evn);
#else
    tMPI_Spinlock_lock( &(evl->lock_rts) );
    ev->next=(struct send_envelope*)evl->head_rts;
    evl->head_rts=ev;
    tMPI_Spinlock_unlock( &(evl->lock_rts) );
#endif
}
#endif



static void tMPI_Send_env_list_remove_old(struct send_envelope *ev)
{
    /* pretty straighforward because it isn't a shared list */
    if (ev->next)
        ev->next->prev=ev->prev; 
    if (ev->prev)
        ev->prev->next=ev->next; 
    ev->prev=NULL;
    ev->next=NULL;
}


static void tMPI_Send_env_list_add_new(struct send_envelope_list *evl, 
                                       struct send_envelope *ev)
{
#ifdef TMPI_LOCK_FREE_LISTS
    struct send_envelope *evl_head_new_orig;
    struct send_envelope *evl_cas;
#endif
    ev->prev=NULL;

#ifdef TMPI_LOCK_FREE_LISTS
    /* behold our lock-free shared linked list: 
       (it's actually quite simple because we only do operations at the head 
        of the list, either adding them - such as here - or detaching the whole
        list) */
    do 
    {
        /* read the old head atomically */
        evl_head_new_orig=(struct send_envelope*)
                          tMPI_Atomic_ptr_get( &(evl->head_new) );
        /* set our envelope to have that as its next */
        ev->next=evl_head_new_orig;
        /* do the compare-and-swap. 
           this operation is a memory fence, so we shouldn't need
           to worry about out-of-order stores */
        evl_cas=(struct send_envelope*)tMPI_Atomic_ptr_cas(&(evl->head_new), 
                                                           evl_head_new_orig, 
                                                           ev);
        /* and compare the results: if they aren't the same,
           somebody else got there before us: */
    } while (evl_cas != evl_head_new_orig); 
#else
    tMPI_Spinlock_lock( &(evl->lock_new) );
    /* we add to the start of the list */
    ev->next=(struct send_envelope*)evl->head_new;
    /* actually attach it to the list */
    evl->head_new=ev;
    tMPI_Spinlock_unlock( &(evl->lock_new) );
#endif

    /* signal to the thread that there is a new envelope */
#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Atomic_fetch_add( &(ev->dest->evs_new_incoming) ,1);
#else
    tMPI_Thread_mutex_lock( &(ev->dest->ev_check_lock) );
    ev->dest->evs_new_incoming++;
    tMPI_Thread_cond_signal( &(ev->dest->ev_check_cond) );
    tMPI_Thread_mutex_unlock( &(ev->dest->ev_check_lock) );
#endif
}

static void tMPI_Send_env_list_move_to_old(struct send_envelope *ev)
{
    struct send_envelope_list *evl=ev->list;

    /* remove from old list. We assume the list has been detached! */
    if (ev->next)
        ev->next->prev=ev->prev; 
    if (ev->prev)
        ev->prev->next=ev->next; 

    /* we add to the end of the list */
    ev->next=evl->head_old;
    ev->prev=evl->head_old->prev;

    ev->next->prev=ev;
    ev->prev->next=ev;
}








/* tmpi_recv_envelope_list functions */

void tMPI_Recv_env_list_init(struct recv_envelope_list *evl)
{
    evl->head = &(evl->dummy);
    evl->head->prev=evl->head;
    evl->head->next=evl->head;
}

void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl)
{
    evl->head=NULL;
}

static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl, 
                                   struct recv_envelope *ev)
{
    ev->list=evl;
    /* we add to the end of the list */
    ev->next=evl->head;
    ev->prev=evl->head->prev;

    ev->next->prev=ev;
    ev->prev->next=ev;
}

static void tMPI_Recv_env_list_remove(struct recv_envelope *ev)
{
    if (ev->next)
        ev->next->prev=ev->prev; 
    if (ev->prev)
        ev->prev->next=ev->next; 
    ev->prev=NULL;
    ev->next=NULL;
    ev->list=NULL;
}








/* tmpi_req functions */

void tMPI_Req_list_init(struct req_list *rl, int N_reqs)
{
    int i;

    rl->alloc_head=(struct tmpi_req_*)tMPI_Malloc(
                                        sizeof(struct tmpi_req_)*N_reqs);
    rl->head=rl->alloc_head;
    for(i=0;i<N_reqs;i++)
    {
        if (i==0)
            rl->head[i].prev=NULL;
        else
            rl->head[i].prev=&(rl->head[i-1]);

        if (i >= (N_reqs-1))
            rl->head[i].next=NULL;
        else
            rl->head[i].next=&(rl->head[i+1]);
    }
}

void tMPI_Req_list_destroy(struct req_list *rl)
{
    free(rl->alloc_head);
    rl->head=NULL;
    rl->alloc_head=NULL;
}

static struct tmpi_req_ *tMPI_Get_req(struct req_list *rl)
{
    struct tmpi_req_ *req=rl->head;
    

    /* we don't need locks here because requests are a per-thread property */
    if (!req)
    {
        /* this could be fixed */
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_REQUESTS);
        return NULL;
    }
    rl->head=req->next;
    req->next=NULL;

    return req;
}

static void tMPI_Return_req(struct req_list *rl, struct tmpi_req_ *req)
{
    req->next=rl->head;
    req->prev=NULL;
    rl->head=req;
}






/* Point-to-point communication protocol functions */



static void tMPI_Set_recv_status(struct recv_envelope *ev, tMPI_Status *status)
{
    if (status)
    {
        status->TMPI_SOURCE = tMPI_Comm_seek_rank(ev->comm, ev->src);
        status->TMPI_TAG = ev->tag;
        status->TMPI_ERROR = ev->error;
        if (tMPI_Atomic_get(&(ev->state))==env_finished)
            status->transferred = ev->bufsize;
        else
            status->transferred = 0;
    }
}

static void tMPI_Set_send_status(struct send_envelope *ev, tMPI_Status *status)
{
    if (status)
    {
        status->TMPI_SOURCE = tMPI_Comm_seek_rank(ev->comm, ev->src);
        status->TMPI_TAG = ev->tag;
        status->TMPI_ERROR = ev->error;
        if (tMPI_Atomic_get(&(ev->state))>env_unmatched)
            status->transferred = ev->bufsize;
        else
            status->transferred = 0;
    }
}



static bool tMPI_Envelope_matches(const struct send_envelope *send,
                                  const struct recv_envelope *recv)
{
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Envelope_matches (%d->%d)==(%d->%d),  tag=(%d==%d),       \n       datatype=(%ld==%ld), comm=(%ld,%ld),\n              finished=(%d==%d)\n",
            tMPI_This_threadnr(),
            tMPI_Threadnr(send->src), tMPI_Threadnr(send->dest),
            tMPI_Threadnr(recv->src), tMPI_Threadnr(recv->dest),
            (int)(send->tag), (int)(recv->tag),
            (long int)send->datatype, (long int)recv->datatype,
            (long int)send->comm, (long int)recv->comm,
            (int)send->state.value, (int)recv->state.value);
    fflush(stdout);
#endif
    if ( ( (recv->tag == TMPI_ANY_TAG) || (recv->tag == send->tag) ) &&
            ( send->comm == recv->comm ) &&
            ( (!recv->src)  || (recv->src == send->src) ) &&
            ( send->dest == recv->dest ) &&
            ( send->datatype == recv->datatype ) &&
            ( send->state.value < env_finished  &&
              recv->state.value == env_unmatched ) )
    {
#ifdef TMPI_DEBUG
        printf("%5d: (%d->%d) tag=%d found match\n",
                tMPI_This_threadnr(),
                tMPI_Threadnr(send->src), tMPI_Threadnr(send->dest),
                (int)(send->tag));
        fflush(stdout);
#endif
        return TRUE;
    }
    return FALSE;
}




static struct send_envelope* 
tMPI_Send_env_list_search_old(struct send_envelope_list *evl,
                              struct recv_envelope *evr)
{
    struct send_envelope *evs;

    evs=(struct send_envelope*)evl->head_old->next;
    while(evs != evl->head_old)
    {
        if (tMPI_Envelope_matches(evs, evr))
        {
            /* remove the envelope */
            tMPI_Send_env_list_remove_old(evs);
            return evs;
        }
        evs=(struct send_envelope*)evs->next;
    }
    return NULL;
}


static struct recv_envelope* 
tMPI_Recv_env_list_search_new(struct recv_envelope_list *evl,
                              struct send_envelope *evs)
{
    struct recv_envelope *evr;

    evr=evl->head->next;
    while(evr != evl->head)
    {
        if (tMPI_Envelope_matches(evs, evr))
        {
            return evr;
        }
        evr=evr->next;
    }
    return NULL;
}


#ifdef USE_SEND_RECV_COPY_BUFFER
static void tMPI_Send_copy_buffer(struct send_envelope *evs)
{
    /* Fill copy buffer, after having anounced its possible use */

    /* in the special case of a zero buffer size, we don't do anything and
       always let the receiver handle it */
    if (evs->bufsize==0) 
        return;

    /* first check whether the other side hasn't started yet */
    tMPI_Atomic_memory_barrier();
    if ((tMPI_Atomic_get( &(evs->state) ) == env_unmatched )) 
    {
        /* first copy */
        memcpy(evs->cb, evs->buf, evs->bufsize);
        /* now set state, if other side hasn't started copying yet. */
        if (tMPI_Atomic_cas( &(evs->state), env_unmatched, env_cb_available)
            == env_unmatched)
        {
            /* if it was originally unmatched, the receiver wasn't 
               copying the old buffer. We can don't need to wait,
               and the receiver is going to clean up this envelope. */
#ifdef TMPI_DEBUG
            printf("%5d: tMPI_Send_copy_buffer(%d->%d, tag=%d) completed\n", 
                   tMPI_This_threadnr(), 
                   tMPI_Threadnr(evs->src), tMPI_Threadnr(evs->dest), 
                   (int)(evs->tag));
            fflush(stdout);
#endif
            return;
        }
    }
    /* we need to wait until the receiver is finished copying */
    while(tMPI_Atomic_get( &(evs->state) ) < env_cb_available)
    {
        tMPI_Atomic_memory_barrier();
    }
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Send_copy_buffer(%d->%d, tag=%d) waiting-completed\n", 
           tMPI_This_threadnr(), 
           tMPI_Threadnr(evs->src), tMPI_Threadnr(evs->dest), (int)(evs->tag));
    fflush(stdout);
#endif
    /* and if we reached this point, we need to clean up the
       envelope ourselves */
    tMPI_Send_env_list_return( evs);
}
#endif


static struct send_envelope*
tMPI_Prep_send_envelope(struct send_envelope_list *evl, tMPI_Comm comm, 
                        struct tmpi_thread *src, struct tmpi_thread *dest, 
                        void *buf, int count, tMPI_Datatype datatype, int tag,
                        bool nonblock)
{
    /* get an envelope from the send-envelope stack */
    struct send_envelope *ev=tMPI_Send_env_list_fetch_new( evl );

    ev->tag=tag;
    ev->nonblock=nonblock;

    ev->comm=comm;

    ev->src=src;
    ev->dest=dest;

    ev->buf=buf;
    ev->bufsize=count*datatype->size;
    ev->datatype=datatype;


#ifdef USE_SEND_RECV_COPY_BUFFER
    /* check whether we'll be double buffering */
    ev->using_cb=(ev->bufsize < COPY_BUFFER_SIZE);
    /* but don't do anything yet */
#endif

    tMPI_Atomic_set(&(ev->state), env_unmatched);
    ev->error=TMPI_SUCCESS;

    return ev;
}

static struct recv_envelope* 
tMPI_Prep_recv_envelope(struct tmpi_thread *cur, tMPI_Comm comm, 
                        struct tmpi_thread *src, struct tmpi_thread *dest, 
                        void *buf, int count, tMPI_Datatype datatype, int tag,
                        bool nonblock)
{
    /* get an envelope from the stack */
    struct recv_envelope *ev=tMPI_Free_env_list_fetch_recv( &(cur->envelopes) );

    ev->tag=tag;
    ev->nonblock=nonblock;

    ev->comm=comm;

    ev->src=src;
    ev->dest=dest;

    ev->buf=buf;
    ev->bufsize=count*datatype->size;
    ev->datatype=datatype;
    
    ev->list=NULL;

    tMPI_Atomic_set(&(ev->state), env_unmatched);
    ev->error=TMPI_SUCCESS;

    return ev;
}









static struct recv_envelope* tMPI_Post_match_recv(tMPI_Comm comm, 
                                                  struct tmpi_thread *src, 
                                                  void *recv_buf, 
                                                  int recv_count, 
                                                  tMPI_Datatype datatype, 
                                                  int tag, bool nonblock)
{
    struct tmpi_thread *cur=tMPI_Get_current();
    struct tmpi_thread *dest=cur;
    struct recv_envelope *evr;
    struct send_envelope *evs=NULL;
    int src_threadnr=src ? tMPI_Threadnr(src) : Nthreads;
    int i;

    /* reserve an envelope to post */
    evr=tMPI_Prep_recv_envelope(cur, comm, src, dest, recv_buf, recv_count, 
                                datatype, tag, nonblock);

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) started\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(evr->src), tMPI_Threadnr(evr->dest), (int)(evr->tag));
    fflush(stdout);
#endif
    /* we now check the entire exisiting send queue */
    if (src)
    {
        evs=tMPI_Send_env_list_search_old( &(dest->evs[src_threadnr]), evr);
    }
    else
    {
        /* if we don't know the source, we look at all possible sources */
        for(i=0;i<Nthreads;i++)
        {
            evs=tMPI_Send_env_list_search_old(&(dest->evs[i]), evr);
            if (evs)
                break;
        } 
    }

    if (evs)
    {
#ifdef TMPI_DEBUG
        printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) found match\n", 
                tMPI_This_threadnr(), 
                tMPI_Threadnr(evr->src), tMPI_Threadnr(evr->dest), 
                (int)(evr->tag));
        fflush(stdout);
#endif
        /* we found a matching send */
        tMPI_Xfer(evs, evr);
    }
    else
    {
#ifdef TMPI_DEBUG
        printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) no match\n", 
                tMPI_This_threadnr(), 
                tMPI_Threadnr(evr->src), tMPI_Threadnr(evr->dest), 
                (int)(evr->tag));
        fflush(stdout);
#endif
        /* we post the envelope in the right list */
        tMPI_Recv_env_list_add( &(dest->evr), evr);
    }
    return evr;
}




static struct send_envelope *tMPI_Post_send(tMPI_Comm comm, 
                                            struct tmpi_thread *dest, 
                                            void *send_buf, int send_count,
                                            tMPI_Datatype datatype, int tag, 
                                            bool nonblock)
{
    struct tmpi_thread *cur=tMPI_Get_current();
    struct tmpi_thread *src=cur;
    struct send_envelope *evs;
    int src_threadnr=tMPI_Threadnr(src);
    struct send_envelope_list *evsl=&(dest->evs[src_threadnr]);

    /* reserve an envelope to post */
    evs=tMPI_Prep_send_envelope(evsl, comm, src, dest, send_buf, send_count, 
                                datatype, tag, nonblock);

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Post_send (%d->%d, tag=%d)\n", 
           tMPI_This_threadnr(), 
           tMPI_Threadnr(evs->src), tMPI_Threadnr(evs->dest), 
           (int)(evs->tag));
    fflush(stdout);
#endif
    /* we post the envelope in the right list */
    tMPI_Send_env_list_add_new( &(dest->evs[src_threadnr]), evs);
   
    return evs;
}




static void tMPI_Test_incoming(struct tmpi_thread *th)
{
    int i;
    int check_id;

    /* we check for newly arrived send envelopes */
#ifndef TMPI_NO_BUSY_WAIT
    tMPI_Atomic_memory_barrier();
    check_id=tMPI_Atomic_get( &(th->evs_new_incoming));
#else
    tMPI_Thread_mutex_lock( &(th->ev_check_lock) );
    if (th->evs_new_incoming == 0 && th->ev_received==0)
    {
        tMPI_Thread_cond_wait( &(th->ev_check_cond), &(th->ev_check_lock) );
        /* we don't need to check for spurious wakeups here, because our 
           conditional predicate is checked later on, and we'll be back here 
           if it isn't satisfied. */
    }
    check_id=th->evs_new_incoming;
    tMPI_Thread_mutex_unlock( &(th->ev_check_lock) );
#endif
    while( check_id > 0)
    {
        /*int repl=check_id;*/
        int n=0;
        /* there were new send envelopes. Let's check them all */
        for(i=0;i<Nthreads;i++)
        {
            struct send_envelope *evs_head;

#ifdef TMPI_LOCK_FREE_LISTS
            /* Behold our lock-free shared linked list:
               (see tMPI_Send_env_list_add_new for more info) */
            struct send_envelope *evl_cas;

            do
            {
                /* read old head atomically */
                evs_head=(struct send_envelope*)
                          tMPI_Atomic_ptr_get( &(th->evs[i].head_new) );
                /* do the compare-and-swap to detach the list */
                evl_cas=(struct send_envelope*)
                          tMPI_Atomic_ptr_cas(&(th->evs[i].head_new),
                                                  evs_head,
                                                  NULL);
            } while (evl_cas != evs_head);
#else
            tMPI_Spinlock_lock( &(th->evs[i].lock_new) );
            evs_head=(struct send_envelope*)th->evs[i].head_new;
            th->evs[i].head_new=NULL; /* detach the list */
            tMPI_Spinlock_unlock( &(th->evs[i].lock_new) );
#endif

            if (evs_head) /* there's a newly arrived send envelope from this 
                             thread*/
            {
                struct send_envelope *evs=evs_head;
                struct send_envelope *prev=NULL;
                struct recv_envelope *evr;

                /* first enable reversing order by creating a regular 
                   doubly-linked list from the singly-linked shared
                   linked list */
                while(evs) 
                {
                    evs->prev=prev;
                    prev=evs;
                    evs=evs->next;
                }
                /* now walk through it backwards (in order of addition) */ 
                evs=prev;
                while(evs)
                {
                    struct send_envelope *evsp=(struct send_envelope*)evs->prev;
                    n++;
                    evr=tMPI_Recv_env_list_search_new(&(th->evr), evs);
                    if (evr)
                    {
                        tMPI_Xfer(evs, evr);
                    }
                    else
                    {
                        tMPI_Send_env_list_move_to_old( evs );
                    }
                    evs=evsp;
                }
            }
        }
        /* we count down with the number of send envelopes we thought we had
           in the beginning */
#ifndef TMPI_NO_BUSY_WAIT
        /*check_id=tMPI_Atomic_fetch_add( &(th->evs_new_incoming), -n);*/
        check_id=tMPI_Atomic_add_return( &(th->evs_new_incoming), -n);
#else
        tMPI_Thread_mutex_lock( &(th->ev_check_lock) );
        check_id = (--th->evs_new_incoming);
        tMPI_Thread_mutex_unlock( &(th->ev_check_lock) );
#endif
    }
}








static bool tMPI_Test_recv(struct recv_envelope *ev, bool blocking, 
                           tMPI_Status *status)
{
    struct tmpi_thread *cur=tMPI_Get_current();

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Test_recv (%d->%d, tag=%d)\n", 
                tMPI_This_threadnr(), 
                tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), 
                (int)(ev->tag));
    fflush(stdout);
#endif

    /* we check for new send envelopes */
    /*tMPI_Test_incoming(cur);*/
    
    /* and check whether the envelope we're waiting for has finished */
    if (blocking)
    {
        while( tMPI_Atomic_get( &(ev->state) ) <  env_finished ) 
        {
            /* while blocking, we wait for incoming send envelopes */
            tMPI_Test_incoming(cur);
        }
    }
    else
    {
        if ( tMPI_Atomic_get( &(ev->state) ) <  env_finished ) 
        {
            tMPI_Test_incoming(cur);
            return FALSE;
        }
    }

    tMPI_Set_recv_status(ev, status);
    tMPI_Free_env_list_return_recv( &(cur->envelopes), ev);

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Test_recv (%d->%d, tag=%d) env_finished\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), (int)(ev->tag));
    fflush(stdout);
#endif
    return TRUE;
}


static bool tMPI_Test_send(struct send_envelope *ev, bool blocking, 
                           tMPI_Status *status)
{
    struct tmpi_thread *cur=tMPI_Get_current();

#ifdef USE_SEND_RECV_COPY_BUFFER
    if (ev->using_cb)
    {
        /* We buffer. Just do the transfer to the buffer and return saying 
           that we're done. It's now up to the receiver to return our 
           envelope.*/
        /* do our transfer. */
        tMPI_Send_copy_buffer(ev);
        tMPI_Set_send_status(ev, status);
        /* do one of these for good measure */
        /*tMPI_Test_incoming(cur);*/
        /* and exit */
        return TRUE;
    }
#endif

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Test_send (%d->%d, tag=%d)\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), 
            (int)(ev->tag));
    fflush(stdout);
#endif
 
    /* we do a check to service all incoming xmissions; this won't affect the
       current wait, of course */
    tMPI_Test_incoming(cur);
    if (blocking)
    {
        while( tMPI_Atomic_get( &(ev->state) ) <  env_finished ) 
        {
            /* while blocking, we wait for incoming send envelopes. That's
               all we can do at this moment. */
            tMPI_Test_incoming(cur);
        }
    }
    else
    {
        if ( tMPI_Atomic_get( &(ev->state) ) <  env_finished ) 
        {
            tMPI_Test_incoming(cur);
            return FALSE;
        }
    }
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_lock( &(cur->ev_check_lock) );
    /* remove one from the received list */
    cur->ev_received--;
    tMPI_Thread_mutex_unlock( &(cur->ev_check_lock) );
#endif

    tMPI_Set_send_status(ev, status);
    tMPI_Send_env_list_return( ev);
 
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Test_send (%d->%d, tag=%d) env_finished\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), (int)(ev->tag));
    fflush(stdout);
#endif
    return TRUE;
}









static void tMPI_Xfer(struct send_envelope *evs, struct recv_envelope *evr)
{
#ifdef USE_SEND_RECV_COPY_BUFFER
    /* we remove the sender's envelope only if we do the transfer, which 
       we always do if the buffer size = 0 */ 
    bool remove_sender = (evs->bufsize==0);
#endif
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Xfer (%d->%d, tag=%d) started\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(evs->src), tMPI_Threadnr(evr->dest), (int)(evs->tag));
    fflush(stdout);
#endif
    /* first set data on the receiving end so status can be updated */
    evr->src = evs->src;
    evr->tag = evs->tag;

    if (evs->bufsize) /* do the actual transfer */
    {
        void *sbuf=evs->buf; /* source buffer */
        if (evs->bufsize > evr->bufsize)
        {
            tMPI_Error((evr->comm), TMPI_ERR_XFER_BUFSIZE);
            tMPI_Atomic_set(&(evr->state), env_finished);
            tMPI_Atomic_set(&(evs->state), env_finished);
            evr->error = TMPI_ERR_XFER_BUFSIZE;
            evs->error = TMPI_ERR_XFER_BUFSIZE;
            return;
        }

#ifdef USE_SEND_RECV_COPY_BUFFER
        if (evs->using_cb)
        {
            /* check if the other side has already finished copying */
            if (tMPI_Atomic_cas( &(evs->state), env_unmatched, env_copying)
                != env_unmatched)
            {
                /* it has, and we're copying from the new buffer. 
                   We're now also tasked with removing the envelope */
                sbuf=evs->cb;
                remove_sender=TRUE;
            }
        }
#endif


        if (!evr->buf || !evs->buf)
        {
            tMPI_Error((evr->comm), TMPI_ERR_BUF);
            tMPI_Atomic_set(&(evr->state), env_finished);
            tMPI_Atomic_set(&(evs->state), env_finished);
            evr->error = TMPI_ERR_BUF;
            evs->error = TMPI_ERR_BUF;
            tMPI_Atomic_memory_barrier();
            return;
        }
        memcpy(evr->buf, sbuf, evs->bufsize);
        /* for status update */
    }
    evr->bufsize=evs->bufsize;
    /* and mark that we're finished */
#ifdef TMPI_NO_BUSY_WAIT
    tMPI_Thread_mutex_lock( &(evr->src->ev_check_lock) );
#endif
 
    tMPI_Atomic_set( &(evr->state), env_finished);
    tMPI_Atomic_set( &(evs->state), env_finished);
    tMPI_Atomic_memory_barrier();
    /* remove the receiving envelope if it's in a list */
    tMPI_Recv_env_list_remove(evr);
#ifdef USE_SEND_RECV_COPY_BUFFER
    if (remove_sender)
    {
        tMPI_Send_env_list_rts(evs);
    }
#endif
#ifdef TMPI_NO_BUSY_WAIT
    evr->src->ev_received++;
    /* wake a potentially sleeping source thread. */
    tMPI_Thread_cond_signal( &(evr->src->ev_check_cond) );
    tMPI_Thread_mutex_unlock( &(evr->src->ev_check_lock) );
#endif

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Xfer (%d->%d, tag=%d) done\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(evs->src), tMPI_Threadnr(evr->dest), (int)(evs->tag));
    fflush(stdout);
#endif
    return;
}











/* point-to-point communication */

int tMPI_Send(void* buf, int count, tMPI_Datatype datatype, int dest,
              int tag, tMPI_Comm comm)
{
    struct send_envelope *sd;
    struct tmpi_thread *send_dst=tMPI_Get_thread(comm, dest);

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Send(%p, %d, %p, %d, %d, %p)", buf, count, 
                       datatype, dest, tag, comm);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!send_dst)
    {
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST);
    }

    sd=tMPI_Post_send(comm, send_dst, buf, count, datatype, tag, FALSE);
    tMPI_Test_send(sd, TRUE, NULL);

    return sd->error;
}

int tMPI_Recv(void* buf, int count, tMPI_Datatype datatype, int source,
             int tag, tMPI_Comm comm, tMPI_Status *status)
{
    struct recv_envelope *rc;
    struct tmpi_thread *recv_src=0;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Recv(%p, %d, %p, %d, %d, %p, %p)", buf, count, 
                       datatype, source, tag, comm, status);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (source!=TMPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC); 
        }
    }

    rc=tMPI_Post_match_recv(comm, recv_src, buf, count, datatype, tag, FALSE);
    tMPI_Test_recv(rc, TRUE, status);

    return rc->error;
}


int tMPI_Sendrecv(void *sendbuf, int sendcount, tMPI_Datatype sendtype,
                  int dest, int sendtag, void *recvbuf, int recvcount,
                  tMPI_Datatype recvtype, int source, int recvtag, 
                  tMPI_Comm comm, tMPI_Status *status)
{
    struct recv_envelope *rc;
    struct send_envelope *sd;

    struct tmpi_thread *recv_src=0;
    struct tmpi_thread *send_dst=tMPI_Get_thread(comm, dest);

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Sendrecv(%p, %d, %p, %d, %d, %p, %d, %p, %d, %d, %p, %p)", 
                       sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                       recvcount, recvtype, source, recvtag, comm, status);
#endif
    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!send_dst)
    {
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST); 
    }
    if (source!=TMPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC);
        }
    }
    /* we first prepare to send */
    sd=tMPI_Post_send(comm, send_dst, sendbuf, sendcount, 
                            sendtype, sendtag, FALSE);
    /* then we prepare to receive */
    rc=tMPI_Post_match_recv(comm, recv_src, recvbuf, recvcount, 
                            recvtype, recvtag, FALSE);

    /* now we wait to receive */
    tMPI_Test_recv(rc, TRUE, status);
    if (rc->error != TMPI_SUCCESS)
        return rc->error;
    /* we wait until the send completes */
    tMPI_Test_send(sd, TRUE, NULL);
    if (rc->error != TMPI_SUCCESS)
        return rc->error;

    return TMPI_SUCCESS;
}


/* async */

int tMPI_Isend(void* buf, int count, tMPI_Datatype datatype, int dest,
               int tag, tMPI_Comm comm, tMPI_Request *request)
{
    struct req_list *rql=&(tMPI_Get_current()->rql);
    struct tmpi_req_ *rq=tMPI_Get_req(rql);
    struct tmpi_thread *send_dst=tMPI_Get_thread(comm, dest);

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Isend(%p, %d, %p, %d, %d, %p, %p)", buf, count, 
                       datatype, dest, tag, comm, request);
#endif
    if (!comm)
    {
        tMPI_Return_req(rql,rq);
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    if (!send_dst)
    {
        tMPI_Return_req(rql,rq);
        return tMPI_Error(comm, TMPI_ERR_SEND_DEST);
    }
    rq->recv=FALSE;
    rq->finished=FALSE;
    rq->evs=tMPI_Post_send(comm, send_dst, buf, count, datatype, tag, TRUE);
    *request=rq;

    return TMPI_SUCCESS;    
}


int tMPI_Irecv(void* buf, int count, tMPI_Datatype datatype, int source,
               int tag, tMPI_Comm comm, tMPI_Request *request)
{
    struct req_list *rql=&(tMPI_Get_current()->rql);
    struct tmpi_req_ *rq=tMPI_Get_req(rql);
    struct tmpi_thread *recv_src=0;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Irecv(%p, %d, %p, %d, %d, %p, %p)", buf, count, 
                       datatype, source, tag, comm, request);
#endif
    if (!comm)
    {
        tMPI_Return_req(rql,rq);
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }

    if (source!=TMPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            tMPI_Return_req(rql,rq);
            return tMPI_Error(comm, TMPI_ERR_RECV_SRC);
        }
    }
    rq->recv=TRUE;
    rq->finished=FALSE;
    rq->evr=tMPI_Post_match_recv(comm, recv_src, buf, count, datatype, tag,
                                 TRUE);

    *request=rq;
    return TMPI_SUCCESS;
}



static int tMPI_Wait_r(struct tmpi_req_ *rq, tMPI_Status *status)
{
    int ret=TMPI_SUCCESS;
    if (!rq)
        return TMPI_SUCCESS;

    if ( ! ( rq->finished ) )
    {
        if (rq->recv)
        {
            /* the receiving end just waits until the data is ready */
            tMPI_Test_recv(rq->evr, TRUE, status);
        }
        else
        {
            /* the sending end also just waits until the data is ready */
            tMPI_Test_send(rq->evs, TRUE, status);
        }
        rq->finished=TRUE;
    }
    rq->evr=NULL; /* we won't be using that envelope any more */
    rq->evs=NULL; /* we won't be using that envelope any more */

    return ret;
}

static int tMPI_Test_r(struct tmpi_req_ *rq, int *flag, tMPI_Status *status)
{    
    int ret=FALSE;

    if (!rq)
        return TMPI_SUCCESS;

    if ( !rq->finished ) 
    {
        if (rq->recv)
        {
            ret=tMPI_Test_recv(rq->evr, FALSE, status);
        }
        else
        {
            ret=tMPI_Test_send(rq->evs, FALSE, status);
        }
        if (ret)
            rq->finished=TRUE;
    }

    if (rq->finished)
    {
        *flag=TRUE;
        /* get rid of the envelope link */
        rq->evr=NULL;
        rq->evs=NULL;
    }
    else
        *flag=FALSE;

    return ret;
}

static int tMPI_Waitall_r(int count, struct tmpi_req_ *array_of_requests[],
                          tMPI_Status *array_of_statuses, bool may_free)
{
    int done;
    int flags_[MAX_PREALLOC_THREADS];
    int *flags=flags_;
    int i;
    struct req_list *rql=&(tMPI_Get_current()->rql);

    if (count > MAX_PREALLOC_THREADS)
    {
#ifdef TMPI_WARN_MALLOC
        fprintf(stderr, "Warning: malloc during tMPI_Waitall_r\n");
#endif
        flags=(int*)tMPI_Malloc(sizeof(int)*count);
    }

    for(i=0;i<count;i++)
        flags[i]=FALSE;
    /* Waitall polls all the requests by calling tMPI_Test. This
       ensures that incoming receives are handled in the order that they
       come in, but of course is busy-wait . */
    do
    {
        /* first take care of duds */
        for(i=0;i<count;i++)
        {
            if (!array_of_requests[i] || array_of_requests[i]->finished)
                flags[i]=TRUE;
        }
        /* do receives */
        for(i=0;i<count;i++)
        {
            if (!flags[i] && array_of_requests[i]->recv)
            {
                struct tmpi_req_ *rq=array_of_requests[i];
                tMPI_Status *status=NULL;
                if (array_of_statuses)
                    status=&(array_of_statuses[i]);

                flags[i]=tMPI_Test_recv(rq->evr, FALSE, status);

                if (rq->evr->error!=TMPI_SUCCESS)
                    return rq->evr->error;
                if (flags[i])
                {
                    rq->evr=NULL;
                    if (may_free)
                    {
                        tMPI_Return_req(rql,rq);
                        array_of_requests[i]=TMPI_REQUEST_NULL;
                    }
                }
            }
        }
        /* then do sends */
        for(i=0;i<count;i++)
        {
            if (!flags[i] && !(array_of_requests[i]->recv))
            {
                struct tmpi_req_ *rq=array_of_requests[i];
                tMPI_Status *status=NULL;
                if (array_of_statuses)
                    status=&(array_of_statuses[i]);

                flags[i]=tMPI_Test_send(rq->evs, FALSE, status);
                                    
                if (rq->evs->error!=TMPI_SUCCESS)
                    return rq->evs->error;
                if (flags[i])
                {
                    rq->evs=NULL;
                    if (may_free)
                    {
                        tMPI_Return_req(rql,rq);
                        array_of_requests[i]=TMPI_REQUEST_NULL;
                    }
                }
            }
        }
        /* count done flags */
        done=0;
        for(i=0;i<count;i++)
        {
            if (flags[i])
                done++;
        }
    }
    while (done<count);

    if (count > MAX_PREALLOC_THREADS)
        free(flags);

    return TMPI_SUCCESS;
}




int tMPI_Wait(tMPI_Request *request, tMPI_Status *status)
{
    int ret=TMPI_SUCCESS;
    struct req_list *rql=&(tMPI_Get_current()->rql);

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Wait(%p, %p)", request, status);
#endif
    if (!request || !(*request))
        return TMPI_SUCCESS;

    ret=tMPI_Wait_r(*request, status);

    /* deallocate if needed */
    tMPI_Return_req(rql, *request);
    *request=TMPI_REQUEST_NULL;
    return ret;
}

int tMPI_Test(tMPI_Request *request, int *flag, tMPI_Status *status)
{
    int ret=TMPI_SUCCESS;
    struct req_list *rql=&(tMPI_Get_current()->rql);

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Test(%p, %p, %p)", request, flag, status);
#endif
    if (!request || !(*request))
        return TMPI_SUCCESS;

    ret=tMPI_Test_r(*request, flag, status);

    if ((*request)->finished)
    {
        /* deallocate if needed */
        tMPI_Return_req(rql, *request);
        *request=TMPI_REQUEST_NULL;
    }
    return ret;
}


int tMPI_Waitall(int count, tMPI_Request *array_of_requests,
                 tMPI_Status *array_of_statuses)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Waitall(%d, %p, %p)", count, array_of_requests, 
                       array_of_statuses);
#endif
    return tMPI_Waitall_r(count, array_of_requests, array_of_statuses, TRUE);
}

