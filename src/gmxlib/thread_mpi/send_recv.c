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

/* Include the defines that determine which thread library to use. 
 * Note that this could also be controlled using preprocessor flags,
 * which is the method used for cmake */
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


#include "thread_mpi.h"
#include "tmpi_impl.h"



/* free envelopes: */
static struct send_envelope *tMPI_Free_env_list_fetch_send(
                                              struct free_envelope_list *evl);
static struct recv_envelope *tMPI_Free_env_list_fetch_recv(
                                              struct free_envelope_list *evl);

/* return an envelope to the free envelopes list */
static void tMPI_Free_env_list_return_send(struct free_envelope_list *evl,
                                           struct send_envelope *ev);
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct recv_envelope *ev);




/* send envelopes: */
/* remove a send envelope from its list. Does not lock */
static void tMPI_Send_env_list_remove(struct send_envelope *ev);

/* add a send envelopelope to the free envelopes list */
static void tMPI_Free_env_list_return_send(struct free_envelope_list *evl,
                                           struct send_envelope *ev);
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct recv_envelope *ev);



/* remove a send envelope from its list. Does not lock */
static void tMPI_Send_env_list_remove(struct send_envelope *ev);

/* add a send envelope to the new envelopes queue in a list */
static void tMPI_Send_env_list_add_new(struct send_envelope_list *evl,
                                       struct send_envelope *ev);
/* move a send envelope to the old envelopes queue in a list. Does not lock. */
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
    evl->send_alloc_head=tMPI_Malloc(sizeof(struct send_envelope)*N );
    evl->recv_alloc_head=tMPI_Malloc(sizeof(struct recv_envelope)*N );
    evl->head_send=evl->send_alloc_head;
    evl->head_recv=evl->recv_alloc_head;

    for(i=0;i<N;i++)
    {
        if (i < N-1)
        {
            evl->head_send[i].next=&(evl->head_send[i+1]);
            evl->head_recv[i].next=&(evl->head_recv[i+1]);
#ifdef TMPI_LOCK_FREE_LISTS
            tMPI_Atomic_ptr_set(&(evl->head_send[i].next_a),NULL);
#else
            evl->head_send[i].next_a=NULL;
#endif
        }
        else
        {
            evl->head_send[i].next=NULL;
            evl->head_recv[i].next=NULL;
#ifdef TMPI_LOCK_FREE_LISTS
            tMPI_Atomic_ptr_set(&(evl->head_send[i].next_a),NULL);
#else
            evl->head_send[i].next_a=NULL;
#endif
        }
        evl->head_send[i].list=NULL;
        evl->head_recv[i].list=NULL;
    }
}

void tMPI_Free_env_list_destroy(struct free_envelope_list *evl)
{
    free(evl->send_alloc_head);
    free(evl->recv_alloc_head);
    evl->head_send=NULL;
    evl->head_recv=NULL;
    evl->send_alloc_head=NULL;
    evl->recv_alloc_head=NULL;
}

static struct send_envelope *tMPI_Free_env_list_fetch_send(
                                            struct free_envelope_list *evl)
{
    struct send_envelope *ret;
    if (! evl->head_send )
    {
        /* TODO: make this do something better than crash */
        fprintf(stderr, "Ran out of envelopes!!!!\n");
        abort();
    }

    ret=evl->head_send;
    evl->head_send=ret->next;
    ret->next=NULL;
    ret->prev=NULL;
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(ret->next_a),NULL);
#else
    ret->next_a=NULL;
#endif
    return ret;
}

static struct recv_envelope *tMPI_Free_env_list_fetch_recv(
                                            struct free_envelope_list *evl)
{
    struct recv_envelope *ret;
    if (! evl->head_recv )
    {
        /* TODO: make this do something better than crash */
        fprintf(stderr, "Ran out of envelopes!!!!\n");
        abort();
    }

    ret=evl->head_recv;
    evl->head_recv=ret->next;
    ret->next=NULL;
    ret->prev=NULL;
    /*evl->N--;*/

    return ret;
}


static void tMPI_Free_env_list_return_send(struct free_envelope_list *evl,
                                           struct send_envelope *ev)
{
    /* just set to NULL to be sure */
    ev->list=NULL;
    ev->prev=NULL;
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(ev->next_a),NULL);
#else
    ev->next_a=NULL;
#endif

    ev->next=evl->head_send;
    evl->head_send=ev;
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

void tMPI_Send_env_list_init(struct send_envelope_list *evl)
{
    tMPI_Spinlock_init( &(evl->lock) );

#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new = NULL;
#endif
    evl->head_old = &(evl->old_dummy);
    evl->head_old->next = evl->head_old;
    evl->head_old->prev = evl->head_old;
}

void tMPI_Send_env_list_destroy(struct send_envelope_list *evl)
{
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new=NULL; 
#endif
    evl->head_old=NULL; /* make it crash if used after tMPI_Finalize */
}

static void tMPI_Send_env_list_remove(struct send_envelope *ev)
{
    if (ev->next)
        ev->next->prev=ev->prev; 
    if (ev->prev)
        ev->prev->next=ev->next; 
    ev->prev=NULL;
    ev->next=NULL;

#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(ev->next_a), NULL);
#else
    ev->next_a = NULL;
#endif

    ev->list=NULL;
}


static void tMPI_Send_env_list_add_new(struct send_envelope_list *evl, 
                                       struct send_envelope *ev)
{
#ifdef TMPI_LOCK_FREE_LISTS
    struct send_envelope *evl_head_new_orig;
    struct send_envelope *evl_cas_result;
#endif
    ev->prev=NULL;
    ev->list=evl;

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
        tMPI_Atomic_ptr_set( &(ev->next_a), evl_head_new_orig );
        /* do the compare-and-swap */
        evl_cas_result=(struct send_envelope*)tMPI_Atomic_ptr_cmpxchg( 
                                                &(evl->head_new), 
                                                evl_head_new_orig,
                                                ev);
        /* and compare the results: if they aren't the same,
           somebody else got there before us: */
    } while (evl_cas_result != evl_head_new_orig); 
#else
    tMPI_Spinlock_lock( &(evl->lock) );
    /* we add to the start of the list */
    ev->next_a=evl->head_new;
    /* actually attach it to the list */
    evl->head_new=ev;
    tMPI_Spinlock_unlock( &(evl->lock) );
#endif

    /* signal to the thread that there is a new envelope */
    tMPI_Atomic_fetch_add( &(ev->dest->evs_check_id) ,1);
}

static void tMPI_Send_env_list_move_to_old(struct send_envelope *ev)
{
    struct send_envelope_list *evl=ev->list;

    /* remove from old list. We assume next_a has been dealt with. */
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

    rl->alloc_head=tMPI_Malloc(sizeof(struct tmpi_req_)*N_reqs);
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
    

    /* we don't need locks here because requests are a per-thread
       property */
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
            status->transferred = (int)(ev->bufsize/ev->datatype->size);
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
        if (tMPI_Atomic_get(&(ev->state))==env_finished)
            status->transferred = (int)(ev->bufsize/ev->datatype->size);
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
            ( send->state.value == env_unmatched &&
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




static struct send_envelope *tMPI_Send_env_list_search_old(
                                    struct send_envelope_list *evl,
                                    struct recv_envelope *evr)
{
    struct send_envelope *evs;

    evs=(struct send_envelope*)evl->head_old->next;
    while(evs != evl->head_old)
    {
        if (tMPI_Envelope_matches(evs, evr))
        {
            /* remove the envelope */
            tMPI_Send_env_list_remove(evs);
            return evs;
        }
        evs=(struct send_envelope*)evs->next;
    }
    return NULL;
}


static struct recv_envelope *tMPI_Recv_env_list_search_new(
                                    struct recv_envelope_list *evl,
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




static struct send_envelope *tMPI_Prep_send_envelope(struct tmpi_thread *cur,
                                                     tMPI_Comm comm, 
                                                     struct tmpi_thread *src, 
                                                     struct tmpi_thread *dest, 
                                                     void *buf, int count, 
                                                     tMPI_Datatype datatype, 
                                                     int tag, bool nonblock)
{
    /* get an envelope from the stack */
    struct send_envelope *ev=tMPI_Free_env_list_fetch_send( &(cur->envelopes) );

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

static struct recv_envelope *tMPI_Prep_recv_envelope(struct tmpi_thread *cur,
                                                     tMPI_Comm comm, 
                                                     struct tmpi_thread *src, 
                                                     struct tmpi_thread *dest, 
                                                     void *buf, int count, 
                                                     tMPI_Datatype datatype, 
                                                     int tag, bool nonblock)
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









static struct recv_envelope *tMPI_Post_match_recv(tMPI_Comm comm, 
                                           struct tmpi_thread *src, 
                                           void *recv_buf, int recv_count,
                                           tMPI_Datatype datatype, int tag, 
                                           bool nonblock)
{
    struct tmpi_thread *cur=tMPI_Get_current();
    struct tmpi_thread *dest=cur;
    struct recv_envelope *evr;
    struct send_envelope *evs=NULL;
    threadnr_t src_threadnr=src ? tMPI_Threadnr(src) : Nthreads;
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
    threadnr_t src_threadnr=tMPI_Threadnr(src);

    /* reserve an envelope to post */
    evs=tMPI_Prep_send_envelope(cur, comm, src, dest, send_buf, send_count, 
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
    check_id=tMPI_Atomic_get( &(th->evs_check_id));
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

            struct send_envelope *evl_cas_result;
            do
            {
                /* read old head atomically */
                evs_head=(struct send_envelope*)
                        tMPI_Atomic_ptr_get( &(th->evs[i].head_new) );
                /* do the compare-and-swap to detach the list */
                evl_cas_result=(struct send_envelope*)tMPI_Atomic_ptr_cmpxchg( 
                                                        &(th->evs[i].head_new),
                                                        evs_head,
                                                        NULL);
            } while (evl_cas_result != evs_head);
#else
            tMPI_Spinlock_lock( &(th->evs[i].lock) );
            evs_head=(struct send_envelope*)th->evs[i].head_new;
            th->evs[i].head_new=NULL; /* detach the list */
            tMPI_Spinlock_unlock( &(th->evs[i].lock) );
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
#ifdef TMPI_LOCK_FREE_LISTS
                    evs->next=(struct send_envelope*)
                               tMPI_Atomic_ptr_get( &(evs->next_a) );
#else
                    evs->next=(struct send_envelope*)evs->next_a;
#endif
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
        check_id=tMPI_Atomic_fetch_add( &(th->evs_check_id), -n);
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
    tMPI_Test_incoming(cur);
    
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

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Test_send (%d->%d, tag=%d)\n", 
            tMPI_This_threadnr(), 
            tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), 
            (int)(ev->tag));
    fflush(stdout);
#endif
 
    /* we do a check to service all incoming sends; this won't affect the
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
            return FALSE;
        }
    }

    tMPI_Set_send_status(ev, status);
    tMPI_Free_env_list_return_send( &(cur->envelopes), ev);
 
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
        if (evs->bufsize > evr->bufsize)
        {
            tMPI_Error((evr->comm), TMPI_ERR_XFER_BUFSIZE);
            tMPI_Atomic_set(&(evr->state), env_finished);
            tMPI_Atomic_set(&(evs->state), env_finished);
            evr->error = TMPI_ERR_XFER_BUFSIZE;
            evs->error = TMPI_ERR_XFER_BUFSIZE;
            return;
        }

        if (!evr->buf || !evs->buf)
        {
            tMPI_Error((evr->comm), TMPI_ERR_BUF);
            tMPI_Atomic_set(&(evr->state), env_finished);
            tMPI_Atomic_set(&(evs->state), env_finished);
            evr->error = TMPI_ERR_BUF;
            evs->error = TMPI_ERR_BUF;
            return;
        }
        memcpy(evr->buf, evs->buf, evs->bufsize);
        /* for status update */
    }
    evr->bufsize=evs->bufsize;
    /* and mark that we're finished */
    tMPI_Atomic_set( &(evr->state), env_finished);
    tMPI_Atomic_set( &(evs->state), env_finished);
    /* remove the receiving envelope if it's in a list */
    tMPI_Recv_env_list_remove(evr);

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
    bool send_finished=FALSE; 
    bool recv_finished=FALSE;


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

    do 
    {
        /* now we wait to receive */
        if (!recv_finished)
        {
            recv_finished=tMPI_Test_recv(rc, FALSE, status);
            if (rc->error != TMPI_SUCCESS)
                return rc->error;
        }   
        /* we wait until the send completes */
        if (!send_finished)
        {
            send_finished=tMPI_Test_send(sd, FALSE, NULL);
            if (rc->error != TMPI_SUCCESS)
                return rc->error;
        }
    }
    while (! (send_finished && recv_finished) );

    return TMPI_SUCCESS;
}


/* async */

int tMPI_Isend(void* buf, int count, tMPI_Datatype datatype, int dest,
              int tag, tMPI_Comm comm, tMPI_Request *request)
{
    struct req_list *rql=&(tMPI_Get_current()->rql);
    struct tmpi_req_ *rq=tMPI_Get_req(rql);
    struct tmpi_thread *send_dst=tMPI_Get_thread(comm, dest);

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

    if ( ! ( rq->finished || 
           ( 
            ( rq->recv && tMPI_Atomic_get( &(rq->evr->state)) == env_finished )
                    ||
            ( !rq->recv && tMPI_Atomic_get( &(rq->evs->state)) == env_finished )
            )
           )
       )
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

    if (rq->recv)
        tMPI_Set_recv_status(rq->evr, status);
    else
        tMPI_Set_send_status(rq->evs, status);

    rq->evr=NULL; /* we won't be using that envelope any more */
    rq->evs=NULL; /* we won't be using that envelope any more */

    return ret;
}

static int tMPI_Test_r(struct tmpi_req_ *rq, int *flag, tMPI_Status *status)
{    
    int ret=TMPI_SUCCESS;
    bool finished=FALSE;

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
        if (finished)
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
        /* then  do sends */
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
    return tMPI_Waitall_r(count, array_of_requests, array_of_statuses, TRUE);
}

