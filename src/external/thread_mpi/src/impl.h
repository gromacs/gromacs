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


/* this is the header file for the implementation side of the thread_mpi
   library. It contains the definitions for all the internal data structures
   and the prototypes for all the internal functions that aren't static.  */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#endif

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "settings.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/threads.h"
#include "thread_mpi/event.h"
#include "thread_mpi/tmpi.h"
#include "thread_mpi/collective.h"
#include "thread_mpi/barrier.h"
#include "thread_mpi/lock.h"
#ifdef TMPI_PROFILE
#include "profile.h"
#endif



/**************************************************************************

   BASIC DEFINITIONS

 **************************************************************************/


typedef int tmpi_bool;
#define TRUE 1
#define FALSE 0



#ifdef USE_COLLECTIVE_COPY_BUFFER
/**************************************************************************

   PRE-ALLOCATED COMMUNICATION BUFFERS

 **************************************************************************/


/* Buffer structure for collective communications. Every thread structure
   has several of these ready to be used when the collective data
   transmission is small enough for double copying to occur (i.e. the size
   of the transmission is less than N*MAX_COPY_BUFFER_SIZE, where N is the
   number of receiving threads).  */
struct copy_buffer
{
    void               *buf;  /* the actual buffer */
    struct copy_buffer *next; /* pointer to next free buffer in buffer_list */
    size_t              size; /* allocated size of buffer */
};

/* a list of copy_buffers of a specific size. */
struct copy_buffer_list
{
    struct copy_buffer *cb;       /* pointer to the first copy_buffer */
    size_t              size;     /* allocated size of buffers in this list */
    struct copy_buffer *cb_alloc; /* list as allocated */
    int                 Nbufs;    /* number of allocated buffers */
};
#endif













/**************************************************************************

   POINT-TO-POINT COMMUNICATION DATA STRUCTURES

 **************************************************************************/

/* the message envelopes (as described in the MPI standard).
   These fully describes the message, and make each message unique (enough).

   Transmitting data works by having the sender put a pointer to an envelope
   onto the receiver's new envelope list corresponding to the originating
   thread.
   The sender then waits until the receiver finishes the transmission, while
   matching all incoming new envelopes against its own list of receive
   envelopes.

   The receiver either directly matches its receiving envelope against
   all previously un-matched sending envelopes, or, if no suitable envelope
   is found, it puts the receive envelope on a receive list.
   Once waiting for completion, the receiver matches against all incoming
   new envelopes.  */

/* the state of an individual point-to-point transmission */
enum envelope_state
{
    env_unmatched       = 0, /* the envelope has not had a match yet */
    env_copying         = 1, /* busy copying (only used for send envelope
                                by receiver if using_cpbuf is true,
                                but cb was still NULL).  */
    env_cb_available    = 2, /* the copy buffer is available. Set by
                                the sender on a send_buffer. */
    env_finished        = 3  /* the transmission has finished */
};


/* the envelope. Held in tmpi_thread->evs[src_thread] for send envelopes,
   or in tmpi_thread->evl for receive envelopes */
struct envelope
{
    int       tag;                  /* the tag */
    tMPI_Comm comm;                 /* this is a structure shared across threads, so we
                                       can test easily whether two threads are talking
                                       about the same comm. */

    struct tmpi_thread *src, *dest; /* these are pretty obvious */

    void               *buf;        /* buffer to be sent  */
    size_t              bufsize;    /* the size of the data to be transmitted */
    tMPI_Datatype       datatype;   /* the data type */

    tmpi_bool           nonblock;   /* whether the receiver is non-blocking */

    /* state, values from enum_envelope_state .
       (there's a few busy-waits relying on this flag).
       status=env_unmatched  is the initial state.*/
    tMPI_Atomic_t state;

    /* the error condition */
    int error;

    /* the message status */
    /*tMPI_Status *status;*/

    /* prev and next envelopes in the send/recv_envelope_list linked list  */
    struct envelope *prev, *next;

    tmpi_bool        send; /* whether this is a send envelope (if TRUE), or a receive
                              envelope (if FALSE) */
#ifdef USE_SEND_RECV_COPY_BUFFER
    tmpi_bool using_cb;    /* whether a copy buffer is (going to be) used */
    void    * cb;          /* the allocated copy buffer pointer */
#endif
    /* the next and previous envelopes in the request list */
    struct envelope *prev_req, *next_req;

    /* the list I'm in */
    struct recv_envelope_list *rlist;
    struct send_envelope_list *slist;
};



/* singly linked lists of free send & receive envelopes belonging to a
   thread. */
struct free_envelope_list
{
    struct envelope *head_recv;       /* the first element in the linked list */
    struct envelope *recv_alloc_head; /* the allocated recv list */
};

/* collection of send envelopes to a specific thread */
struct send_envelope_list
{
    struct envelope *head_free;  /* singly linked list with free send
                                    envelopes. A single-thread LIFO.*/
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_t head_new;  /* singly linked list with the new send
                                    envelopes (i.e. those that are put there by
                                    the sending thread, but not yet checked by
                                    the receiving thread). This is a lock-free
                                    shared detachable list.*/
    tMPI_Atomic_ptr_t head_rts;  /* singly linked list with free send
                                    envelopes returned by the other thread.
                                    This is a lock-free shared LIFO.*/
#else
    struct envelope *head_new;   /* singly linked list with the new send
                                    envelopes (i.e. those that are put there by
                                    the sending thread, but not yet checked by
                                    the receiving thread). */
    struct envelope *head_rts;   /* singly linked list with free send envelopes */
    tMPI_Spinlock_t  lock_new;   /* this locks head_new */
    tMPI_Spinlock_t  lock_rts;   /* this locks head_rts */
#endif
    struct envelope *head_old;   /* the old send envelopes, in a circular doubly
                                    linked list. These have been checked by the
                                    receiving thread against the existing
                                    recv_envelope_list. */

    struct envelope *alloc_head; /* the allocated send list */
    size_t           Nalloc;     /* number of allocted sends */
};

struct recv_envelope_list
{
    struct envelope *head;  /* first envelope in this list */
    struct envelope  dummy; /* the dummy element for the list */
};


/* the request object for asynchronious operations. */
struct tmpi_req_
{
    tmpi_bool           finished;    /* whether it's finished */
    struct envelope    *ev;          /* the envelope */

    struct tmpi_thread *source;      /* the message source (for receives) */
    tMPI_Comm           comm;        /* the comm */
    int                 tag;         /* the tag */
    int                 error;       /* error code */
    size_t              transferred; /* the number of transferred bytes */
    tmpi_bool           cancelled;   /* whether the transmission was canceled */

    struct tmpi_req_   *next, *prev; /* next,prev request in linked list,
                                        used in the req_list, but also in
                                        tMPI_Test_mult().  */
};

/* pre-allocated  request object list */
struct req_list
{
    struct tmpi_req_ *head;       /* pre-allocated singly linked list of requests.
                                     (i.e. reqs->prev is undefined). */
    struct tmpi_req_ *alloc_head; /* the allocated block */
};
















/**************************************************************************

   MULTICAST COMMUNICATION DATA STRUCTURES

 **************************************************************************/

/* these are data structures meant for keeping track of multicast operations
   (tMPI_Bcast, tMPI_Gather, etc.). Because these operations are all collective
   across the comm, and are always blocking, the protocol can be much simpler
   than that for point-to-point communication through tMPI_Send/Recv, etc. */

/* unique tags for multicast & collective operations */
#define TMPI_BCAST_TAG      1
#define TMPI_GATHER_TAG     2
#define TMPI_GATHERV_TAG    3
#define TMPI_SCATTER_TAG    4
#define TMPI_SCATTERV_TAG   5
#define TMPI_REDUCE_TAG     6
#define TMPI_ALLTOALL_TAG   7
#define TMPI_ALLTOALLV_TAG  8


/* thread-specific part of the coll_env */
struct coll_env_thread
{
    tMPI_Atomic_t current_sync; /* sync counter value for the current
                                   communication */
    tMPI_Atomic_t n_remaining;  /* remaining threads count for each thread */

    int           tag;          /* collective communication type */
    tMPI_Datatype datatype;     /* datatype */

    void        **buf;          /* array of send/recv buffer values */
    size_t       *bufsize;      /* array of number of bytes to send/recv */

#ifdef USE_COLLECTIVE_COPY_BUFFER
    tmpi_bool     using_cb;      /* whether a copy buffer is (going to be) used */
    tMPI_Atomic_t buf_readcount; /* Number of threads reading from buf
                                    while using_cpbuf is true, but cpbuf
                                    is still NULL.  */
    tMPI_Atomic_ptr_t  *cpbuf;   /* copy_buffer pointers. */
    struct copy_buffer *cb;      /* the copy buffer cpbuf points to */
#endif

    tMPI_Event send_ev;   /* event associated with being the sending thread.
                             Triggered when last receiving thread is ready,
                             and the coll_env_thread is ready for re-use. */
    tMPI_Event recv_ev;   /* event associated with being a receiving thread. */

    tmpi_bool *read_data; /* whether we read data from a specific thread. */
};

/* Collective communications once sync. These run in parallel with
   the collection of coll_env_threads*/
struct coll_env_coll
{
    /* collective sync data */
    tMPI_Atomic_t current_sync; /* sync counter value for the current
                                   communication */
    tMPI_Atomic_t n_remaining;  /* remaining threads count */

    void         *res;          /* result data for once calls. */
};

/* the collective communication envelope. There's a few of these per
   comm, and each one stands for one collective communication call.  */
struct coll_env
{
    struct coll_env_thread *met; /* thread-specific collective envelope data.*/

    struct coll_env_coll    coll;
    int                     N;
};

/* multicast synchronization data structure. There's one of these for
   each thread in each tMPI_Comm structure */
struct coll_sync
{
    int         synct;  /* sync counter for coll_env_thread.  */
    int         syncs;  /* sync counter for coll_env_coll.  */

    tMPI_Event *events; /* One event for each other thread */
    int         N;      /* the number of threads */
};











/**************************************************************************

   THREAD DATA STRUCTURES

 **************************************************************************/

/* information about a running thread. This structure is put in a
   globally available array; the envelope exchange, etc. are all done through
   the elements of this array.*/
struct tmpi_thread
{
    tMPI_Thread_t thread_id; /* this thread's id */

    /* p2p communication structures: */

    /* the receive envelopes posted for other threads to check */
    struct recv_envelope_list  evr;
    /* the send envelopes posted by other threadas */
    struct send_envelope_list *evs;
    /* free send and receive envelopes */
    struct free_envelope_list  envelopes;
    /* number of finished send envelopes */
    tMPI_Atomic_t              ev_outgoing_received;
    /* the p2p communication events (incoming envelopes + finished send
       envelopes generate events) */
    tMPI_Event      p2p_event;
    TMPI_YIELD_WAIT_DATA /* data associated with waiting */
    struct req_list rql; /* list of pre-allocated requests */

    /* collective communication structures: */
#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* copy buffer list for multicast communications */
    struct copy_buffer_list cbl_multi;
#endif

    /* miscellaneous data: */

    tMPI_Comm self_comm; /* comms for MPI_COMM_SELF */
#ifdef TMPI_PROFILE
    /* the per-thread profile structure that keeps call counts & wait times. */
    struct tmpi_profile profile;
#endif
    /* The start function (or NULL, if a main()-style start function is to
       be called) */
    void  (*start_fn)(void*);
    /* The main()-style start function */
    int   (*start_fn_main)(int, char**);
    /* the argument to the start function, if it's not main()*/
    void *start_arg;

    /* we copy these for each thread (providing these to main() is not
       required by the MPI standard, but it's convenient). Note that we copy,
       because some programs (like Gromacs) like to manipulate these. */
    int    argc;
    char **argv;
};






/**************************************************************************

   ERROR HANDLER DATA STRUCTURES

 **************************************************************************/


/* the error handler  */
struct tmpi_errhandler_
{
    int                err;
    tMPI_Errhandler_fn fn;
};

/* standard error handler functions */
void tmpi_errors_are_fatal_fn(tMPI_Comm *comm, int *err);
void tmpi_errors_return_fn(tMPI_Comm *comm, int *err);





/**************************************************************************

   GLOBAL DATA STRUCTURE

 **************************************************************************/

/* global MPI information */
struct tmpi_global
{
    /* list of pointers to all user-defined types */
    struct tmpi_datatype_ **usertypes;
    int                     N_usertypes;
    int                     Nalloc_usertypes;

    /* spinlock/mutex for manipulating tmpi_user_types */
    tMPI_Spinlock_t datatype_lock;

    /* Lock to prevent multiple threads manipulating the linked list of comm
       structures.*/
    tMPI_Thread_mutex_t comm_link_lock;

    /* barrier for tMPI_Finalize(), etc. */
    tMPI_Thread_barrier_t barrier;

    /* the timer for tMPI_Wtime() */
    tMPI_Thread_mutex_t timer_mutex;
#if !(defined( _WIN32 ) || defined( _WIN64 ) )
    /* the time at initialization. */
    struct timeval timer_init;
#else
    /* the time at initialization. */
    DWORD timer_init;
#endif
};















/**************************************************************************

   COMMUNICATOR DATA STRUCTURES

 **************************************************************************/


struct tmpi_group_
{
    int                  N;     /* the number of threads */
    struct tmpi_thread **peers; /* the list of peers to communicate with */
#if 0
    int                  Nrefs; /* the number of references to this structure */
#endif
};


/* the communicator objects are globally shared. */
struct tmpi_comm_
{
    struct tmpi_group_ grp; /* the communicator group */

    /* the barrier for tMPI_Barrier() */
    tMPI_Barrier_t barrier;


    /* List of barriers for reduce operations.
       reduce_barrier[0] contains a list of N/2 barriers for N threads
       reduce_barrier[1] contains a list of N/4 barriers for N/2 threads
       reduce_barrier[2] contains a list of N/8 barriers for N/4 threads
       and so on. (until N/x reaches 1)
       This is to facilitate tree-based algorithms for tMPI_Reduce, etc.  */
    tMPI_Barrier_t **reduce_barrier;
    int             *N_reduce;      /* the number of barriers in each iteration */
    int              N_reduce_iter; /* the number of iterations */


    struct coll_env  *cev;   /* list of multicast envelope objecs */
    struct coll_sync *csync; /* list of multicast sync objecs */

    /* lists of globally shared send/receive buffers for tMPI_Reduce. */
    tMPI_Atomic_ptr_t *reduce_sendbuf, *reduce_recvbuf;

    /* mutex for communication object creation. Traditional mutexes are
       better here because communicator creation should not be done in
       time-critical sections of code.   */
    tMPI_Thread_mutex_t comm_create_lock;
    tMPI_Thread_cond_t  comm_create_prep;
    tMPI_Thread_cond_t  comm_create_finish;

    tMPI_Comm          *new_comm; /* newly created communicators */

    /* the split structure is shared among the comm threads and is
       allocated & deallocated during tMPI_Comm_split */
    struct tmpi_split *split;

    /* the topologies (only cartesian topology is currently implemented */
    struct cart_topol *cart;
    /*struct tmpi_graph_topol_ *graph;*/

    tMPI_Errhandler erh;

    /* links for a global circular list of all comms that starts at
       TMPI_COMM_WORLD. Used to de-allocate the comm structures after
       tMPI_Finalize(). */
    struct tmpi_comm_ *next, *prev;

    /* A counter that counts up to N before the comm is freed. */
    tMPI_Atomic_t destroy_counter;
};



/* specific for tMPI_Split: */
struct tmpi_split
{
    volatile int       Ncol_init;
    volatile int       Ncol_destroy;
    volatile tmpi_bool can_finish;
    volatile int      *colors;
    volatile int      *keys;
};

/* cartesian topology */
struct cart_topol
{
    int  ndims;   /* number of dimensions */
    int *dims;    /* procs per coordinate */
    int *periods; /* whether the grid is periodic, per dimension */
};

#if 0
/* graph topology */
struct tmpi_graph_topol_
{
};
#endif





/**************************************************************************

   DATA TYPE DATA STRUCTURES

 **************************************************************************/

/* tMPI_Reduce Op functions */
typedef void (*tMPI_Op_fn)(void*, void*, void*, int);


struct tmpi_datatype_component
{
    struct tmpi_datatype_ *type;
    unsigned int           count;
};

/* we don't support datatypes with holes (yet)  */
struct tmpi_datatype_
{
    size_t                          size;         /* full extent of type. */
    tMPI_Op_fn                     *op_functions; /* array of op functions for this datatype */
    int                             N_comp;       /* number of components */
    struct tmpi_datatype_component *comps;        /* the components */
    tmpi_bool                       committed;    /* whether the data type is committed */
};
/* just as a shorthand:  */
typedef struct tmpi_datatype_ tmpi_dt;








/**************************************************************************

   GLOBAL VARIABLES

 **************************************************************************/


/* the threads themselves (tmpi_comm only contains lists of pointers to this
         structure */
extern struct tmpi_thread *threads;
extern int                 Nthreads;

/* thread info */
extern tMPI_Thread_key_t id_key; /* the key to get the thread id */

/* misc. global information about MPI */
extern struct tmpi_global *tmpi_global;








/**************************************************************************

   FUNCTION PROTOTYPES & MACROS

 **************************************************************************/

#ifdef TMPI_TRACE
void tMPI_Trace_print(const char *fmt, ...);
#endif

/* error-checking malloc/realloc: */
void *tMPI_Malloc(size_t size);
void *tMPI_Realloc(void *p, size_t size);
void tMPI_Free(void *p);


/* get the current thread structure pointer */
#define tMPI_Get_current() ((struct tmpi_thread*) \
                            tMPI_Thread_getspecific(id_key))

/* get the number of this thread */
/*#define tMPI_This_threadnr() (tMPI_Get_current() - threads)*/

/* get the number of a specific thread. We convert to the resulting size_t to
   int, which is unlikely to cause problems in the foreseeable future. */
#define tMPI_Threadnr(th) (int)(th - threads)

/* get thread associated with rank  */
#define tMPI_Get_thread(comm, rank) (comm->grp.peers[rank])


#if 0
/* get the current thread structure pointer */
struct tmpi_thread *tMPI_Get_current(void);
/* get the thread belonging to comm with rank rank */
struct tmpi_thread *tMPI_Get_thread(tMPI_Comm comm, int rank);

#endif

/* handle an error, returning the errorcode */
int tMPI_Error(tMPI_Comm comm, int tmpi_errno);



/* check whether we're the main thread */
tmpi_bool tMPI_Is_master(void);
/* check whether the current process is in a group */
tmpi_bool tMPI_In_group(tMPI_Group group);

/* find the rank of a thread in a comm */
int tMPI_Comm_seek_rank(tMPI_Comm comm, struct tmpi_thread *th);
/* find the size of a comm */
int tMPI_Comm_N(tMPI_Comm comm);

/* allocate a comm object, making space for N threads */
int tMPI_Comm_alloc(tMPI_Comm *newcomm, tMPI_Comm parent, int N);
/* de-allocate a comm object */
int tMPI_Comm_destroy(tMPI_Comm comm, tmpi_bool do_link_lock);
/* allocate a group object */
tMPI_Group tMPI_Group_alloc(void);

/* topology functions */
/* de-allocate a cartesian topology structure. (it is allocated with
   the internal function tMPI_Cart_init()) */
void tMPI_Cart_destroy(struct cart_topol *top);






/* initialize a free envelope list with N envelopes */
int tMPI_Free_env_list_init(struct free_envelope_list *evl, int N);
/* destroy a free envelope list */
void tMPI_Free_env_list_destroy(struct free_envelope_list *evl);


/* initialize a send envelope list */
int tMPI_Send_env_list_init(struct send_envelope_list *evl, int N);
/* destroy a send envelope list */
void tMPI_Send_env_list_destroy(struct send_envelope_list *evl);






/* initialize a recv envelope list */
int tMPI_Recv_env_list_init(struct recv_envelope_list *evl);
/* destroy a recv envelope list */
void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl);




/* initialize request list */
int tMPI_Req_list_init(struct req_list *rl, int N_reqs);
/* destroy request list */
void tMPI_Req_list_destroy(struct req_list *rl);



/* collective data structure ops */


/* initialize a coll env structure */
int tMPI_Coll_env_init(struct coll_env *mev, int N);
/* destroy a coll env structure */
void tMPI_Coll_env_destroy(struct coll_env *mev);

/* initialize a coll sync structure */
int tMPI_Coll_sync_init(struct coll_sync *msc, int N);
/* destroy a coll sync structure */
void tMPI_Coll_sync_destroy(struct coll_sync *msc);

#ifdef USE_COLLECTIVE_COPY_BUFFER
/* initialize a copy_buffer_list */
int tMPI_Copy_buffer_list_init(struct copy_buffer_list *cbl, int Nbufs,
                               size_t size);
/* initialize a copy_buffer_list */
void tMPI_Copy_buffer_list_destroy(struct copy_buffer_list *cbl);
/* get a copy buffer from a list */
struct copy_buffer *tMPI_Copy_buffer_list_get(struct copy_buffer_list *cbl);
/* return a copy buffer to a list */
void tMPI_Copy_buffer_list_return(struct copy_buffer_list *cbl,
                                  struct copy_buffer      *cb);
/* initialize a copy buffer */
int tMPI_Copy_buffer_init(struct copy_buffer *cb, size_t size);
void tMPI_Copy_buffer_destroy(struct copy_buffer *cb);
#endif


/* reduce ops: run a single iteration of a reduce operation on a, b -> dest */
int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b,
                       tMPI_Datatype datatype, int count, tMPI_Op op,
                       tMPI_Comm comm);


/* and we need this prototype */
int main(int argc, char **argv);
