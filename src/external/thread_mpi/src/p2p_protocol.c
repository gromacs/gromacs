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
#include "p2p.h"
#include "unused.h"

/* free envelopes: */
static struct envelope *tMPI_Free_env_list_fetch_recv(struct free_envelope_list
                                                      *evl);

/* return an envelope to the free envelopes list */
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct envelope           *rev);



/* send envelope lists: */
/* send envelopes: */
/* get a new envelope from the send list's free envelope list */
static struct envelope* tMPI_Send_env_list_fetch_new(struct
                                                     send_envelope_list *evl);

/* return a send envelope to the send list's free envelope list,
   (to be used by the sending thread, who owns the send_envelope_list) */
static void tMPI_Send_env_list_return(struct envelope *ev);
#ifdef USE_SEND_RECV_COPY_BUFFER
/* return a send envelope to the sender's send list.
   (to be used by the receiving thread). */
static void tMPI_Send_env_list_rts(struct envelope *sev);
#endif




/* send envelopes: */
/* remove a send envelope from its head_old list. Does not lock */
static void tMPI_Send_env_list_remove_old(struct envelope *sev);

/* add a send envelope to the new envelopes queue in a list */
static void tMPI_Send_env_list_add_new(struct tmpi_thread        *cur,
                                       struct send_envelope_list *evl,
                                       struct envelope           *sev);
/* move a send envelope to the old envelopes queue in a list.
   Assumes that this is safe to do without interference
   from other threads, i.e. the list it's in must have been
   detached. */
static void tMPI_Send_env_list_move_to_old(struct envelope *sev);




/* receive envelopes: */
/* add a receive envelope to a list */
static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl,
                                   struct envelope           *ev);
/* remove a receive envelope from its list */
static void tMPI_Recv_env_list_remove(struct envelope *ev);




/* do the actual point-to-point transfer */
static void tMPI_Xfer(struct tmpi_thread *cur, struct envelope *sev,
                      struct envelope *rev);




/* Point-to-point communication protocol functions */
int tMPI_Free_env_list_init(struct free_envelope_list *evl, int N)
{
    int i;

    /* allocate the head element */
    evl->recv_alloc_head = (struct envelope*)tMPI_Malloc(sizeof(struct envelope)
                                                         *N);
    if (evl->recv_alloc_head == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    evl->head_recv = evl->recv_alloc_head;

    for (i = 0; i < N; i++)
    {
        if (i < N-1)
        {
            evl->head_recv[i].next = &(evl->head_recv[i+1]);
        }
        else
        {
            evl->head_recv[i].next = NULL;
        }
        evl->head_recv[i].rlist = NULL;
        evl->head_recv[i].slist = NULL;
    }
    return TMPI_SUCCESS;
}

void tMPI_Free_env_list_destroy(struct free_envelope_list *evl)
{
    free(evl->recv_alloc_head);
    evl->head_recv       = NULL;
    evl->recv_alloc_head = NULL;
}

static struct envelope* tMPI_Free_env_list_fetch_recv(struct
                                                      free_envelope_list *evl)
{
    struct envelope *ret;
    if (!evl->head_recv)
    {
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_ENVELOPES);
        return NULL;
    }

    ret            = evl->head_recv;
    evl->head_recv = ret->next;
    ret->next      = NULL;
    ret->prev      = NULL;
    /*evl->N--;*/

    return ret;
}

static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct envelope           *rev)
{
    rev->rlist     = NULL;
    rev->slist     = NULL;
    rev->prev      = NULL;
    rev->next      = evl->head_recv;
    evl->head_recv = rev;
}

/* tmpi_send_envelope_list functions */

int tMPI_Send_env_list_init(struct send_envelope_list *evl, int N)
{
    int i;
#ifndef TMPI_LOCK_FREE_LISTS
    tMPI_Spinlock_init( &(evl->lock_rts) );
    tMPI_Spinlock_init( &(evl->lock_new) );
#endif
    evl->Nalloc = N;

    evl->alloc_head = (struct envelope*)tMPI_Malloc(sizeof(struct envelope)*N);
    if (evl->alloc_head == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    for (i = 0; i < N; i++)
    {
        evl->alloc_head[i].next  = (i < (N-1)) ? &(evl->alloc_head[i+1]) : NULL;
        evl->alloc_head[i].prev  = NULL;
        evl->alloc_head[i].slist = evl;
        evl->alloc_head[i].rlist = NULL;
#ifdef USE_SEND_RECV_COPY_BUFFER
        evl->alloc_head[i].cb = (void*)tMPI_Malloc(sizeof(char)*
                                                   COPY_BUFFER_SIZE);
        if (evl->alloc_head[i].cb == NULL)
        {
            return TMPI_ERR_NO_MEM;
        }
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
    evl->head_old       = evl->alloc_head; /* the first element is a dummy */
    evl->head_old->next = evl->head_old;
    evl->head_old->prev = evl->head_old;
    return TMPI_SUCCESS;
}

void tMPI_Send_env_list_destroy(struct send_envelope_list *evl)
{
#ifdef USE_SEND_RECV_COPY_BUFFER
    size_t i;
    for (i = 0; i < evl->Nalloc; i++)
    {
        free(evl->alloc_head[i].cb);
    }
#endif
    free(evl->alloc_head);
    evl->alloc_head = NULL;
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new = NULL;
#endif
    evl->head_old = NULL; /* make it crash if used after tMPI_Finalize */
}


static struct envelope* tMPI_Send_env_list_fetch_new(struct
                                                     send_envelope_list *evl)
{
    struct envelope *ret;

    do
    {
        /* first check whether any envelopes were returned to sender */
#ifdef TMPI_LOCK_FREE_LISTS
        if ((ret = (struct envelope*)tMPI_Atomic_ptr_get(&(evl->head_rts))))
#else
        if (evl->head_rts)
#endif
        {
            /* detach the list */
#ifdef TMPI_LOCK_FREE_LISTS
            /* we detach by swapping what we expect the pointer value to be,
               with NULL. If there were a cross-platform way to atomically
               swap  without checking, we could do that, too. */
            while (!tMPI_Atomic_ptr_cas( &(evl->head_rts), ret, NULL ))
            {
                ret = (struct envelope*)tMPI_Atomic_ptr_get(&(evl->head_rts));
            }
#else
            tMPI_Spinlock_lock( &(evl->lock_rts) );
            ret           = evl->head_rts;
            evl->head_rts = NULL;
            tMPI_Spinlock_unlock( &(evl->lock_rts) );
#endif
            /* now add the items to head_free */
            while (ret)
            {
                struct envelope *next = ret->next;
                ret->next      = evl->head_free;
                evl->head_free = ret;
                ret            = next;
            }
        }

        /* get the last free one off the list */
        ret = evl->head_free;
        if (!ret)
#ifdef USE_SEND_RECV_COPY_BUFFER
        {
            /* There are no free send envelopes, so all we can do is handle
               incoming requests until we get a free send envelope. */
#if defined(TMPI_DEBUG)  || defined(TMPI_WARNINGS)
            printf("Ran out of send envelopes!!\n");
            fflush(stdout);
#endif
            tMPI_Wait_process_incoming(tMPI_Get_current());
        }
#else
        {
            /* If this happens, it most likely indicates a bug in the
               calling program. We could fix the situation by waiting,
               but that would most likely lead to deadlocks - even
               more difficult to debug than this. */
            tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_ENVELOPES);
            return NULL;
        }
#endif
    }
    while (!ret);

    evl->head_free = ret->next;

    ret->next  = NULL;
    ret->prev  = NULL;
    ret->slist = evl;
    ret->rlist = NULL;

    /* and return it */
    return ret;
}

static void tMPI_Send_env_list_return(struct envelope *sev)
{
    struct send_envelope_list *evl = sev->slist;

    sev->next      = evl->head_free;
    evl->head_free = sev;
}


#ifdef USE_SEND_RECV_COPY_BUFFER
static void tMPI_Send_env_list_rts(struct envelope *sev)
{
    struct send_envelope_list *evl = sev->slist;
#ifdef TMPI_LOCK_FREE_LISTS
    struct envelope           *sevn;

    do
    {
        sevn      = (struct envelope*)tMPI_Atomic_ptr_get(&evl->head_rts);
        sev->next = sevn;
        /* the cmpxchg operation is a memory fence, so we shouldn't need
           to worry about out-of-order evaluation */
    }
    while (!tMPI_Atomic_ptr_cas( &(evl->head_rts), sevn, sev ));
#else
    tMPI_Spinlock_lock( &(evl->lock_rts) );
    ev->next      = (struct envelope*)evl->head_rts;
    evl->head_rts = sev;
    tMPI_Spinlock_unlock( &(evl->lock_rts) );
#endif
}
#endif



static void tMPI_Send_env_list_remove_old(struct envelope *sev)
{
    /* pretty straighforward because it isn't a shared list */
    if (sev->next)
    {
        sev->next->prev = sev->prev;
    }
    if (sev->prev)
    {
        sev->prev->next = sev->next;
    }
    sev->prev = NULL;
    sev->next = NULL;
}


static void tMPI_Send_env_list_add_new(struct tmpi_thread        tmpi_unused *cur,
                                       struct send_envelope_list             *evl,
                                       struct envelope                       *sev)
{
#ifdef TMPI_LOCK_FREE_LISTS
    struct envelope *evl_head_new_orig;
#endif
    sev->prev = NULL;

#ifdef TMPI_LOCK_FREE_LISTS
    /* behold our lock-free shared linked list:
       (it's actually quite simple because we only do operations at the head
        of the list, either adding them - such as here - or detaching the whole
        list) */
    do
    {
        /* read the old head atomically */
        evl_head_new_orig = (struct envelope*) tMPI_Atomic_ptr_get(
                    &(evl->head_new) );
        /* set our envelope to have that as its next */
        sev->next = evl_head_new_orig;
        /* do the compare-and-swap.
           this operation is a memory fence, so we shouldn't need
           to worry about out-of-order stores. If it returns false,
           somebody else got there before us: */
    }
    while (!tMPI_Atomic_ptr_cas(&(evl->head_new), evl_head_new_orig, sev));

#else
    tMPI_Spinlock_lock( &(evl->lock_new) );
    /* we add to the start of the list */
    sev->next = (struct send_envelope*)evl->head_new;
    /* actually attach it to the list */
    evl->head_new = sev;
    tMPI_Spinlock_unlock( &(evl->lock_new) );
#endif

#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_start(cur);
#endif
    /* signal to the thread that there is a new envelope */
    tMPI_Event_signal( &(sev->dest->p2p_event) );
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_P2p_signal);
#endif
}

static void tMPI_Send_env_list_move_to_old(struct envelope *sev)
{
    struct send_envelope_list *evl = sev->slist;

    /* remove from old list. We assume the list has been detached! */
    if (sev->next)
    {
        sev->next->prev = sev->prev;
    }
    if (sev->prev)
    {
        sev->prev->next = sev->next;
    }

    /* we add to the end of the list */
    sev->next = evl->head_old;
    sev->prev = evl->head_old->prev;

    sev->next->prev = sev;
    sev->prev->next = sev;
}

/* tmpi_recv_envelope_list functions */

int tMPI_Recv_env_list_init(struct recv_envelope_list *evl)
{
    evl->head       = &(evl->dummy);
    evl->head->prev = evl->head;
    evl->head->next = evl->head;

    return TMPI_SUCCESS;
}

void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl)
{
    evl->head = NULL;
}

static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl,
                                   struct envelope           *rev)
{
    rev->rlist = evl;
    /* we add to the end of the list */
    rev->next = evl->head;
    rev->prev = evl->head->prev;

    rev->next->prev = rev;
    rev->prev->next = rev;
}

static void tMPI_Recv_env_list_remove(struct envelope *rev)
{
    if (rev->next)
    {
        rev->next->prev = rev->prev;
    }
    if (rev->prev)
    {
        rev->prev->next = rev->next;
    }
    rev->prev  = NULL;
    rev->next  = NULL;
    rev->rlist = NULL;
}

/* tmpi_req functions */

int tMPI_Req_list_init(struct req_list *rl, int N_reqs)
{
    int i;

    rl->alloc_head = (struct tmpi_req_*)tMPI_Malloc(
                sizeof(struct tmpi_req_)*N_reqs);
    if (rl->alloc_head == 0)
    {
        return TMPI_ERR_NO_MEM;
    }
    rl->head = rl->alloc_head;
    for (i = 0; i < N_reqs; i++)
    {
        if (i == 0)
        {
            rl->head[i].prev = NULL;
        }
        else
        {
            rl->head[i].prev = &(rl->head[i-1]);
        }

        if (i >= (N_reqs-1))
        {
            rl->head[i].next = NULL;
        }
        else
        {
            rl->head[i].next = &(rl->head[i+1]);
        }
    }
    return TMPI_SUCCESS;
}

void tMPI_Req_list_destroy(struct req_list *rl)
{
    free(rl->alloc_head);
    rl->head       = NULL;
    rl->alloc_head = NULL;
}



struct tmpi_req_ *tMPI_Get_req(struct req_list *rl)
{
    struct tmpi_req_ *req = rl->head;


    /* we don't need locks here because requests are a per-thread property */
    if (!req)
    {
        /* this could be fixed */
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_REQUESTS);
        return NULL;
    }
    rl->head  = req->next;
    req->next = NULL;

    return req;
}

void tMPI_Return_req(struct req_list *rl, struct tmpi_req_ *req)
{
    req->next = rl->head;
    req->prev = NULL;
    rl->head  = req;
}

void tMPI_Req_init(struct tmpi_req_ *rq, struct envelope *ev)
{
    rq->ev       = ev;
    rq->finished = FALSE;
    rq->next     = rq;
    rq->prev     = rq;

    rq->source      = ev->src;
    rq->comm        = ev->comm;
    rq->tag         = TMPI_ANY_TAG;
    rq->error       = TMPI_SUCCESS;
    rq->transferred = 0;
    rq->cancelled   = FALSE;
}

/* Point-to-point communication protocol functions */
void tMPI_Set_req(struct envelope *ev, struct tmpi_req_ *req)
{
    req->source = ev->src;
    req->comm   = ev->comm;
    req->tag    = ev->tag;
    req->error  = ev->error;
    if (ev->send)
    {
        if (tMPI_Atomic_get(&(ev->state)) > env_unmatched)
        {
            req->transferred = ev->bufsize;
        }
        else
        {
            req->transferred = 0;
        }
    }
    else
    {
        if (tMPI_Atomic_get(&(ev->state)) == env_finished)
        {
            req->transferred = ev->bufsize;
        }
        else
        {
            req->transferred = 0;
        }
    }
}

void tMPI_Set_status(struct tmpi_req_ *req, tMPI_Status *st)
{
    if (st)
    {
        st->TMPI_SOURCE = tMPI_Comm_seek_rank(req->comm, req->source);
        st->TMPI_TAG    = req->tag;
        st->TMPI_ERROR  = req->error;
        st->transferred = req->transferred;
        st->cancelled   = req->cancelled;
    }
}

tmpi_bool tMPI_Envelope_matches(const struct envelope *sev,
                                const struct envelope *rev)
{
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Envelope_matches (%d->%d)==(%d->%d),  tag=(%d==%d),       \n       datatype=(%ld==%ld), comm=(%ld,%ld),\n              finished=(%d==%d)\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(sev->src), tMPI_Threadnr(sev->dest),
           tMPI_Threadnr(rev->src), tMPI_Threadnr(rev->dest),
           (int)(sev->tag), (int)(rev->tag),
           (long int)sev->datatype, (long int)rev->datatype,
           (long int)sev->comm, (long int)rev->comm,
           (int)sev->state.value, (int)rev->state.value);
    fflush(stdout);
#endif
    if ( ( (rev->tag == TMPI_ANY_TAG) || (rev->tag == sev->tag) ) &&
         ( sev->comm == rev->comm ) &&
         ( (!rev->src)  || (rev->src == sev->src) ) &&
         ( sev->dest == rev->dest ) &&
         ( sev->datatype == rev->datatype ) &&
         ( tMPI_Atomic_get(&(sev->state)) < env_finished  &&
           tMPI_Atomic_get(&(rev->state)) == env_unmatched ) )
    {
#ifdef TMPI_DEBUG
        printf("%5d: (%d->%d) tag=%d found match\n",
               tMPI_This_threadnr(),
               tMPI_Threadnr(sev->src), tMPI_Threadnr(sev->dest),
               (int)(sev->tag));
        fflush(stdout);
#endif
        return TRUE;
    }
    return FALSE;
}

struct envelope* tMPI_Send_env_list_search_old(struct send_envelope_list *evl,
                                               struct envelope           *rev)
{
    struct envelope *sev;

    sev = (struct envelope*)evl->head_old->next;
    while (sev != evl->head_old)
    {
        if (tMPI_Envelope_matches(sev, rev))
        {
            /* remove the envelope */
            tMPI_Send_env_list_remove_old(sev);
            return sev;
        }
        sev = sev->next;
    }
    return NULL;
}

struct envelope* tMPI_Recv_env_list_search_new(struct recv_envelope_list *evl,
                                               struct envelope           *sev)
{
    struct envelope *rev;

    rev = evl->head->next;
    while (rev != evl->head)
    {
        if (tMPI_Envelope_matches(sev, rev))
        {
            return rev;
        }
        rev = rev->next;
    }
    return NULL;
}

#ifdef USE_SEND_RECV_COPY_BUFFER
void tMPI_Send_copy_buffer(struct envelope *sev, struct tmpi_req_ *req)
{
    int state;
    /* Fill copy buffer, after having anounced its possible use */

    /* in the special case of a zero buffer size, we don't do anything and
       always let the receiver handle it */
    if (sev->bufsize == 0)
    {
        return;
    }

    /* first check whether the other side hasn't started yet */
    state = tMPI_Atomic_get( &(sev->state) );
    tMPI_Atomic_memory_barrier_acq();
    if (state == env_unmatched)
    {
        /* first copy */
        memcpy(sev->cb, sev->buf, sev->bufsize);
        /* now set state, if other side hasn't started copying yet. */
        tMPI_Atomic_memory_barrier_rel();
        if (tMPI_Atomic_cas( &(sev->state), env_unmatched, env_cb_available))
        {
            /* if it was originally unmatched, the receiver wasn't
               copying the old buffer. We can don't need to wait,
               and the receiver is going to clean up this envelope. */
#ifdef TMPI_DEBUG
            printf("%5d: tMPI_Send_copy_buffer(%d->%d, tag=%d) completed\n",
                   tMPI_This_threadnr(),
                   tMPI_Threadnr(sev->src), tMPI_Threadnr(sev->dest),
                   (int)(sev->tag));
            fflush(stdout);
#endif
            return;
        }
    }
    /* and if we reached this point, the receiver had already started
       copying, and we need to clean up the envelope ourselves.

       we first need to wait until the receiver is finished copying. We
       know this is a short wait (since the buffer was small enough to be
       buffered in the first place), so we just spin-wait.  */
    tMPI_Atomic_memory_barrier(); /* a full barrier to make sure that the
                                     sending doesn't interfere with the
                                     waiting */
    while (tMPI_Atomic_get( &(sev->state) ) < env_cb_available)
    {
        tMPI_Atomic_memory_barrier_acq();
    }
    tMPI_Atomic_memory_barrier_acq();
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Send_copy_buffer(%d->%d, tag=%d) waiting-completed\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(sev->src), tMPI_Threadnr(sev->dest), (int)(sev->tag));
    fflush(stdout);
#endif
    tMPI_Set_req(sev, req);
    /* and now we clean up */
    tMPI_Send_env_list_return(sev);
}
#endif

struct envelope* tMPI_Prep_send_envelope(struct send_envelope_list *evl,
                                         tMPI_Comm comm,
                                         struct tmpi_thread *src,
                                         struct tmpi_thread *dest,
                                         void *buf, int count,
                                         tMPI_Datatype datatype,
                                         int tag, tmpi_bool nonblock)
{
    /* get an envelope from the send-envelope stack */
    struct envelope *ev = tMPI_Send_env_list_fetch_new( evl );
    if (ev == NULL)
    {
        return NULL;
    }

    ev->tag      = tag;
    ev->nonblock = nonblock;

    ev->comm = comm;

    ev->src  = src;
    ev->dest = dest;

    ev->buf      = buf;
    ev->bufsize  = count*datatype->size;
    ev->datatype = datatype;

    ev->send = TRUE;

    ev->rlist = NULL;

#ifdef USE_SEND_RECV_COPY_BUFFER
    /* check whether we'll be double buffering */
    ev->using_cb = (ev->bufsize < COPY_BUFFER_SIZE);
    /* but don't do anything yet */
#endif

    tMPI_Atomic_set(&(ev->state), env_unmatched);

    ev->error = TMPI_SUCCESS;
    if (count < 0)
    {
        tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
        ev->error = TMPI_ERR_XFER_BUFSIZE;
    }

    return ev;
}

struct envelope* tMPI_Prep_recv_envelope(struct tmpi_thread *cur,
                                         tMPI_Comm comm,
                                         struct tmpi_thread *src,
                                         struct tmpi_thread *dest,
                                         void *buf, int count,
                                         tMPI_Datatype datatype, int tag,
                                         tmpi_bool nonblock)
{
    /* get an envelope from the stack */
    struct envelope *ev = tMPI_Free_env_list_fetch_recv( &(cur->envelopes) );
    if (ev == NULL)
    {
        return NULL;
    }

    ev->tag      = tag;
    ev->nonblock = nonblock;

    ev->comm = comm;

    ev->src  = src;
    ev->dest = dest;

    ev->buf      = buf;
    ev->bufsize  = count*datatype->size;
    ev->datatype = datatype;

    ev->send = FALSE;

    ev->slist = NULL;
    ev->rlist = NULL;

    tMPI_Atomic_set(&(ev->state), env_unmatched);

    ev->error = TMPI_SUCCESS;
    if (count < 0)
    {
        tMPI_Error(comm, TMPI_ERR_XFER_BUFSIZE);
        ev->error = TMPI_ERR_XFER_BUFSIZE;
    }

    return ev;
}

static void tMPI_Xfer(struct tmpi_thread tmpi_unused *cur, struct envelope *sev,
                      struct envelope *rev)
{
#ifdef USE_SEND_RECV_COPY_BUFFER
    /* we remove the sender's envelope only if we do the transfer, which
       we always do if the buffer size = 0 */
    tmpi_bool remove_sender = (sev->bufsize == 0);
#endif
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Xfer (%d->%d, tag=%d) started\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(sev->src), tMPI_Threadnr(rev->dest), (int)(sev->tag));
    fflush(stdout);
#endif
    /* first set data on the receiving end so status can be updated */
    rev->src = sev->src;
    rev->tag = sev->tag;

    if (sev->bufsize)          /* do the actual transfer */
    {
        void *sbuf = sev->buf; /* source buffer */
        if (sev->bufsize > rev->bufsize)
        {
            tMPI_Error((rev->comm), TMPI_ERR_XFER_BUFSIZE);
            tMPI_Atomic_set(&(rev->state), env_finished);
            tMPI_Atomic_set(&(sev->state), env_finished);
            rev->error = TMPI_ERR_XFER_BUFSIZE;
            sev->error = TMPI_ERR_XFER_BUFSIZE;
            return;
        }

#ifdef USE_SEND_RECV_COPY_BUFFER
        if (sev->using_cb)
        {
            /* check if the other side has already finished copying */
            if (!tMPI_Atomic_cas( &(sev->state), env_unmatched, env_copying))
            {
                /* it has, and we're copying from the new buffer.
                   We're now also tasked with removing the envelope */
                sbuf          = sev->cb;
                remove_sender = TRUE;
#ifdef TMPI_PROFILE
                tMPI_Profile_count_buffered_p2p_xfer(cur);
#endif
            }
        }
#endif

        if (!rev->buf || !sev->buf)
        {
            tMPI_Error((rev->comm), TMPI_ERR_BUF);
            tMPI_Atomic_set(&(rev->state), env_finished);
            tMPI_Atomic_set(&(sev->state), env_finished);
            rev->error = TMPI_ERR_BUF;
            sev->error = TMPI_ERR_BUF;
            return;
        }
        memcpy(rev->buf, sbuf, sev->bufsize);
#ifdef TMPI_PROFILE
        tMPI_Profile_count_p2p_xfer(cur);
#endif
        /* for status update */
    }
    rev->bufsize = sev->bufsize;
    /* and mark that we're finished */
#if defined(TMPI_PROFILE)
    {
        tMPI_Profile_wait_start(cur);
#endif
    tMPI_Atomic_set( &(rev->state), env_finished);
    tMPI_Atomic_set( &(sev->state), env_finished);

    /* signal to a potentially waiting thread that we're done. */
    tMPI_Atomic_fetch_add( &(rev->src->ev_outgoing_received), 1);
    tMPI_Event_signal(&(rev->src->p2p_event));

    /* remove the receiving envelope if it's in a list */
    tMPI_Recv_env_list_remove(rev);
#ifdef USE_SEND_RECV_COPY_BUFFER
    if (remove_sender)
    {
        tMPI_Send_env_list_rts(sev);
    }
#endif
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_P2p_signal);
}
#endif


#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Xfer (%d->%d, tag=%d) done\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(sev->src), tMPI_Threadnr(rev->dest), (int)(sev->tag));
    fflush(stdout);
#endif
    return;
}

struct envelope* tMPI_Post_match_recv(struct tmpi_thread *cur,
                                      tMPI_Comm comm,
                                      struct tmpi_thread *src,
                                      void *recv_buf, int recv_count,
                                      tMPI_Datatype datatype,
                                      int tag, tmpi_bool nonblock)
{
    struct tmpi_thread *dest = cur;
    struct envelope    *rev;
    struct envelope    *sev          = NULL;
    int                 src_threadnr = src ? tMPI_Threadnr(src) : Nthreads;
    int                 i;

    /* reserve an envelope to post */
    rev = tMPI_Prep_recv_envelope(cur, comm, src, dest, recv_buf, recv_count,
                                  datatype, tag, nonblock);
    if (rev == NULL)
    {
        return NULL;
    }

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) started\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(rev->src), tMPI_Threadnr(rev->dest), (int)(rev->tag));
    fflush(stdout);
#endif
    /* we now check the entire exisiting send queue */
    if (src)
    {
        sev = tMPI_Send_env_list_search_old( &(dest->evs[src_threadnr]), rev);
    }
    else
    {
        /* if we don't know the source, we look at all possible sources */
        for (i = 0; i < Nthreads; i++)
        {
            sev = tMPI_Send_env_list_search_old(&(dest->evs[i]), rev);
            if (sev)
            {
                break;
            }
        }
    }

    if (sev)
    {
#ifdef TMPI_DEBUG
        printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) found match\n",
               tMPI_This_threadnr(),
               tMPI_Threadnr(rev->src), tMPI_Threadnr(rev->dest),
               (int)(rev->tag));
        fflush(stdout);
#endif
        /* we found a matching send */
        tMPI_Xfer(cur, sev, rev);
    }
    else
    {
#ifdef TMPI_DEBUG
        printf("%5d: tMPI_Post_match_recv (%d->%d, tag=%d) no match\n",
               tMPI_This_threadnr(),
               tMPI_Threadnr(rev->src), tMPI_Threadnr(rev->dest),
               (int)(rev->tag));
        fflush(stdout);
#endif
        /* we post the envelope in the right list */
        tMPI_Recv_env_list_add( &(dest->evr), rev);
    }
    return rev;
}

struct envelope *tMPI_Post_send(struct tmpi_thread *cur,
                                tMPI_Comm comm,
                                struct tmpi_thread *dest,
                                void *send_buf, int send_count,
                                tMPI_Datatype datatype, int tag,
                                tmpi_bool nonblock)
{
    struct tmpi_thread        *src = cur;
    struct envelope           *sev;
    int                        src_threadnr = tMPI_Threadnr(src);
    struct send_envelope_list *sevl         = &(dest->evs[src_threadnr]);

    /* reserve an envelope to post */
    sev = tMPI_Prep_send_envelope(sevl, comm, src, dest, send_buf, send_count,
                                  datatype, tag, nonblock);
    if (sev == NULL)
    {
        return NULL;
    }

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Post_send (%d->%d, tag=%d)\n",
           tMPI_This_threadnr(),
           tMPI_Threadnr(sev->src), tMPI_Threadnr(sev->dest),
           (int)(sev->tag));
    fflush(stdout);
#endif
    /* we post the envelope in the right list */
    tMPI_Send_env_list_add_new(cur, &(dest->evs[src_threadnr]), sev);

    return sev;
}

void tMPI_Wait_process_incoming(struct tmpi_thread *cur)
{
    int i;
    int check_id;
    int n_handled = 0;


#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_start(cur);
#endif
    /* we check for newly arrived send envelopes and finished send
       envelopes */
    check_id = tMPI_Event_wait( &(cur->p2p_event));
    /* the outgoing_received items are handled 'automatically'
       by the function calling this function */
#if defined(TMPI_PROFILE)
    tMPI_Profile_wait_stop(cur, TMPIWAIT_P2p);
#endif
    n_handled = tMPI_Atomic_get(&(cur->ev_outgoing_received));
    tMPI_Atomic_fetch_add( &(cur->ev_outgoing_received), -n_handled);
    check_id -= n_handled;

    if (check_id > 0)
    {
        /*int repl=check_id;*/
        /*int n=0;*/
        /* there were new send envelopes. Let's check them all */
        for (i = 0; i < Nthreads; i++)
        {
            struct envelope *sev_head;

#ifdef TMPI_LOCK_FREE_LISTS
            /* Behold our lock-free shared linked list:
               (see tMPI_Send_env_list_add_new for more info) */
            do
            {
                /* read old head atomically */
                sev_head = (struct envelope*)
                    tMPI_Atomic_ptr_get( &(cur->evs[i].head_new) );
                /* do the compare-and-swap to detach the list */
            }
            while (!tMPI_Atomic_ptr_cas(&(cur->evs[i].head_new), sev_head,
                                        NULL));
#else
            tMPI_Spinlock_lock( &(cur->evs[i].lock_new) );
            sev_head             = (struct send_envelope*)cur->evs[i].head_new;
            cur->evs[i].head_new = NULL; /* detach the list */
            tMPI_Spinlock_unlock( &(cur->evs[i].lock_new) );
#endif

            if (sev_head) /* there's a newly arrived send envelope from this
                             thread*/
            {
                struct envelope *sev    = sev_head;
                struct envelope *prev_s = NULL;
                struct envelope *rev;

                /* first enable reversing order by creating a regular
                   doubly-linked list from the singly-linked shared
                   linked list */
                while (sev)
                {
                    sev->prev = prev_s;
                    prev_s    = sev;
                    sev       = sev->next;
                }
                /* now walk through it backwards (in order of addition) */
                sev = prev_s;
                while (sev)
                {
                    struct envelope *sevp = sev->prev;
                    n_handled++;
                    rev = tMPI_Recv_env_list_search_new(&(cur->evr), sev);
                    if (rev)
                    {
                        tMPI_Xfer(cur, sev, rev);
                    }
                    else
                    {
                        tMPI_Send_env_list_move_to_old( sev );
                    }
                    sev = sevp;
                }
            }
        }
    }
    tMPI_Event_process( &(cur->p2p_event), n_handled);
}

tmpi_bool tMPI_Test_single(struct tmpi_thread *cur, struct tmpi_req_ *rq)
{
    struct envelope *ev = rq->ev;

    if (ev && !(rq->finished) )
    {
#ifdef USE_SEND_RECV_COPY_BUFFER
        if (ev->send && ev->using_cb)
        {
            /* We buffer-copy. Just do the transfer to the buffer and
               return saying that we're done. It's now up to the
               receiver to return our envelope.*/
            /* do our transfer and are guaranteed a finished
               envelope. */
            tMPI_Send_copy_buffer(ev, rq);
            /* get the results */
            rq->error    = rq->ev->error;
            rq->finished = TRUE;
        }
        else
#endif
        {
            if (tMPI_Atomic_get( &(ev->state) ) >= env_finished)
            {
                rq->finished = TRUE;
                /* get the results */
                rq->error = rq->ev->error;
                tMPI_Set_req(ev, rq);
                /* and release the envelope. After this point, the envelope
                   may be reused, so its contents shouldn't be relied on. */
                if (ev->send)
                {
                    tMPI_Send_env_list_return(ev);
                }
                else
                {
                    tMPI_Free_env_list_return_recv( &(cur->envelopes), ev);
                }
            }
        }
    }
    return rq->finished;
}

void tMPI_Wait_single(struct tmpi_thread *cur, struct tmpi_req_ *rq)
{
    do
    {
        if (tMPI_Test_single(cur, rq))
        {
            return;
        }
        tMPI_Wait_process_incoming(cur);
    }
    while (TRUE);
}

tmpi_bool tMPI_Test_multi(struct tmpi_thread *cur, struct tmpi_req_ *rqs,
                          tmpi_bool *any_done)
{
    tmpi_bool         all_done = TRUE;
    struct tmpi_req_ *creq     = rqs;

    int               i = 0;
    if (any_done)
    {
        *any_done = FALSE;
    }

    while (creq)
    {
        tmpi_bool finished = tMPI_Test_single(cur, creq);
        i++;

        /* now do the check */
        if (!finished)
        {
            all_done = FALSE;
        }
        else
        {
            /* remove the request from the list we've been given. */
            if (creq->prev)
            {
                creq->prev->next = creq->next;
            }
            if (creq->next)
            {
                creq->next->prev = creq->prev;
            }
            if (any_done)
            {
                *any_done = TRUE;
            }
        }

        creq = creq->next;
    }

    return all_done;
}
