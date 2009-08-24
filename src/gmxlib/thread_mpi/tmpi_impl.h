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




/*#define TMPI_DEBUG*/

/* if this is defined, MPI will warn/hang/crash on practices that don't conform
   to the MPI standard (such as not calling tMPI_Comm_free on all threads that
   are part of the comm being freed). */
#define TMPI_STRICT 



/* whether to warn if there are mallocs at performance-critical sections 
   (due to preallocations being too small) */
#define TMPI_WARN_MALLOC




/* the max. number of envelopes waiting to be read at thread */
#define MAX_ENVELOPES 128

/* the normal maximum number of threads for pre-defined arrays
   (if the actual number of threads is bigger than this, it'll
    allocate/deallocate arrays, so no problems will arise).*/
#define MAX_PREALLOC_THREADS 32


/* whether to use lock&wait-free lists using compare-and-swap (cmpxchg on x86)
   pointer functions. Message passing using blocking Send/Recv is still
   blocking, of course. */
#define TMPI_LOCK_FREE_LISTS



/* BASIC DEFINITIONS */


#ifndef __cplusplus
typedef int bool;
#define TRUE 1
#define FALSE 0
#else
#ifndef TRUE
#define TRUE true
#endif
#ifndef FALSE
#define FALSE false
#endif
#endif



/* POINT-TO-POINT COMMUNICATION DATA STRUCTURES */

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
    env_finished        = 1  /* the transmission has finished */
};


/* the send envelope. Held in tmpi_thread->evs[src_thread]  */
struct send_envelope
{
    int tag; /* the tag */
    tMPI_Comm comm; /* this is a structure shared across threads, so we
                      can test easily whether two threads are talking
                      about the same comm. */

    struct tmpi_thread *src, *dest; /* these are pretty obvious */

    void *buf; /* buffer to be sent  */
    size_t bufsize; /* the size of the data to be transmitted */
    tMPI_Datatype datatype; /* the data type */

    bool nonblock; /* whether the receiver is non-blocking */

     /* state, values from enum_envelope_state .  
        this is volatile because it's probably read by a thread 
        while another one is writing to them (there's a few busy-waits 
        relying on these flags). 
        status=env_unmatched  is the initial state.*/
    tMPI_Atomic_t state;

    /* the error condition */
    int error;

    /* the message status */
    tMPI_Status *status;
    /* prev and next envelopes in the linked list  */
    struct send_envelope *prev,*next;
#ifdef TMPI_LOCK_FREE_LISTS
    /* next pointer for lock-free lists */
    tMPI_Atomic_ptr_t next_a; 
#else
    /* next pointer for shared lists */
    volatile struct send_envelope *next_a; 
#endif
    /* the list I'm in */
    struct send_envelope_list *list;
};

/* the receive envelope. Held in tmpi_thread->evl  */
struct recv_envelope
{
    int tag; /* the tag */
    /* transmission type flags */


    tMPI_Comm comm; /* this is a structure shared across threads, so we
                      can test easily whether two threads are talking
                      about the same comm. */

    struct tmpi_thread *src, *dest; /* these are pretty obvious */

    void *buf; /* buffer to be sent  */
    size_t bufsize; /* the size of the data to be transmitted, or the
                       size of the actually transmitted data after Xfer. */
    tMPI_Datatype datatype; /* the data type */

    bool nonblock; /* whether the receiver is non-blocking */

     /* state, values from enum_envelope_state .  
        this is volatile because it's probably read by a thread 
        while another one is writing to them (there's a few busy-waits 
        relying on these flags). 
        status=env_unmatched  is the initial state.*/
    tMPI_Atomic_t state;

    /* the error condition */
    int error;

    /* prev and next envelopes in the linked list  */
    struct recv_envelope *prev,*next;
    /* the list I'm in */
    struct recv_envelope_list *list;
};


/* singly linked lists of free send & receive envelopes belonging to 
   a thread. */
struct free_envelope_list 
{
    struct send_envelope *head_send; /* the first element in the 
                                             linked list */
    struct recv_envelope *head_recv; /* the first element in the 
                                             linked list */

    struct send_envelope *send_alloc_head;  /* the allocated send list */
    struct recv_envelope *recv_alloc_head;  /* the allocated recv list */
};

/* collection of send envelopes to a specific thread */
struct send_envelope_list
{
#ifdef TMPI_LOCK_FREE_LISTS
    tMPI_Atomic_ptr_t head_new; /* singly linked list with the new send 
                                  envelopes (i.e. those that are put there 
                                  by the sending thread, but not yet checked 
                                  by the receiving thread). */
#else
    volatile struct send_envelope *head_new; /* singly linked list with 
                                                the new send envelopes 
                                                (i.e. those that are put 
                                                there by the sending thread,
                                                but not yet checked by the
                                                receiving thread). */
#endif
    struct send_envelope *head_old; /* the old send envelopes,
                                       in a circular doubly linked 
                                       list. These have been checked
                                       by the receiving thread 
                                       against the existing 
                                       recv_envelope_list. */

    tMPI_Spinlock_t lock; /* this locks head_new */

    struct send_envelope old_dummy; /* the dummy element for the head_old
                                       list */
};

struct recv_envelope_list
{
    struct recv_envelope *head; /* first envelope in this list */
    struct recv_envelope dummy; /* the dummy element for the list */
};


/* the request object for asynchronious operations. */
struct tmpi_req_
{
    bool recv; /* whether it's a receive request. */
    bool finished; /* whether it's finished */

    struct send_envelope *evs; /* the envelope */
    struct recv_envelope *evr; /* the envelope */

    struct tmpi_status_ st;

    struct tmpi_req_ *next,*prev; /* next,prev request in linked list */
};

/* pre-allocated  request object list */
struct req_list
{
    struct tmpi_req_ *head; /* pre-allocated singly linked list of requests. 
                              (i.e. reqs->prev is undefined). */
    struct tmpi_req_ *alloc_head; /* the allocated block */
};








/* MULTICAST COMMUNICATION DATA STRUCTURES */

/* these are data structures meant for keeping track of multicast operations
   (tMPI_Bcast, tMPI_Gather, etc.). Because these operations are all collective
   across the comm, and are always blocking, the protocol can be much simpler
   than that for point-to-point communication through tMPI_Send/Recv, etc. */

/* unique tags for multicast */
#define TMPI_BCAST_TAG      1
#define TMPI_GATHER_TAG     2
#define TMPI_GATHERV_TAG    3
#define TMPI_SCATTER_TAG    4
#define TMPI_SCATTERV_TAG   5
#define TMPI_REDUCE_TAG     6
#define TMPI_ALLTOALL_TAG   7
#define TMPI_ALLTOALLV_TAG  8



/* the multicast envelope. There's a few of these for each 
   multi_sync structure */
struct multi_env
{
    tMPI_Atomic_t current_counter; /* counter value for the current
                                     communication */
    tMPI_Atomic_t n_remaining; /* remaining count for root thread */


    int tag; /* multicast communication type */
    volatile void **buf; /* array of send/recv buffer values */
    volatile size_t *bufsize; /* array of number of bytes to send/recv */
    bool *read_data; /* whether we read data from a specific thread*/
    tMPI_Datatype datatype;
};

/* multicast synchronization data structure. There's one of these for 
   each thread in each tMPI_Comm structure */
struct multi_sync
{
    int counter; /* sync counter for list in mev */
#define N_MULTI_SYNC 2
    struct multi_env mev[N_MULTI_SYNC];
};










/* GLOBALLY AVAILABLE DATA STRUCTURES */

/* information about a running thread. This structure is put in a 
   globally available array; the envelope exchange, etc. are all done through
   the elements of this array.*/
struct tmpi_thread
{
    tMPI_Thread_t thread_id;

    /* the receive envelopes posted for other threads to check */
    struct recv_envelope_list evr;

    /* the send envelopes posted by other threadas */
    struct send_envelope_list *evs;

    /* change indicator */
    tMPI_Atomic_t evs_check_id;

    struct req_list rql;  /* the requests */

    /* the free envelopes */
    struct free_envelope_list envelopes; 

    tMPI_Comm self_comm;

    void (*start_fn)(void*); /* The start function (or NULL, if main() is to be
                                called) */
    void *start_arg; /* the argument to the start function (if it's not main()*/

    /* we copy these for each thread (this is not required by the 
       MPI standard, but it's convenient). Note that we copy, because
       some programs (like Gromacs) like to manipulate these. */    
    int argc;
    char **argv;
};


/* the error handler  */
struct tmpi_errhandler_
{
    int err;
    tMPI_Errhandler_fn fn;
};

/* standard error handler functions */
void tmpi_errors_are_fatal_fn(tMPI_Comm *comm, int *err);
void tmpi_errors_return_fn(tMPI_Comm *comm, int *err);


/* global MPI information */
struct tmpi_global
{
    /* list of pointers to all user-defined types */
    struct tmpi_datatype_ **usertypes;
    int N_usertypes;
    int Nalloc_usertypes;

    /* spinlock/mutex for manipulating tmpi_user_types */
    tMPI_Spinlock_t  datatype_lock;
};


















/* COMMUNICATOR DATA STRUCTURES */



struct tmpi_group_
{
    int N; /* the number of threads */
    struct tmpi_thread **peers; /* the list of peers to communicate with */
#if 0
    int Nrefs; /* the number of references to this structure */
#endif
};


/* the communicator objects are globally shared. */
struct tmpi_comm_
{
    struct tmpi_group_ grp; /* the communicator group */

    /* list of barrier_t's. 
       multicast_barrier[0] contains a barrier for N threads 
       (N=the number of threads in the communicator)
       multicast_barrier[1] contains a barrier for N/2 threads
       multicast_barrier[2] contains a barrier for N/4 threads
       multicast_barrier[3] contains a barrier for N/8 threads
       and so on. (until N/x reaches 1)
       This is to facilitate tree-based algorithms for tMPI_Reduce, etc.  */
    tMPI_Spinlock_barrier_t *multicast_barrier;   
    int *N_multicast_barrier;
    int Nbarriers;


    struct multi_sync *msc; /* list of multicast sync objecs */

    /* lists of globally shared send/receive buffers for tMPI_Reduce, etc. */
    volatile void **sendbuf, **recvbuf; 
    
    /* mutex for communication object creation. Traditional mutexes are 
       better here because communicator creation should not be done in 
       time-critical sections of code.   */ 
    tMPI_Thread_mutex_t comm_create_lock;
    tMPI_Thread_cond_t comm_create_prep; 
    tMPI_Thread_cond_t comm_create_finish;

    tMPI_Comm *new_comm; /* newly created communicators */
    struct tmpi_split *split;

    /* the topologies (only cartesian topology is currently implemented */
    struct cart_topol *cart;
    /*struct tmpi_graph_topol_ *graph;*/

    tMPI_Errhandler erh;
};


/* specific for tMPI_Split: */
struct tmpi_split
{ 
    volatile int Ncol_init;
    volatile int Ncol_destroy;
    volatile bool can_finish;
    volatile int *colors;
    volatile int *keys;
};

/* cartesian topology */
struct cart_topol
{
    int ndims; /* number of dimensions */
    int *dims; /* procs per coordinate */
    int *periods; /* whether the grid is periodic, per dimension */
};

#if 0
/* graph topology */
struct tmpi_graph_topol_
{
};
#endif







/* tMPI_Reduce Op functions */
typedef void (*tMPI_Op_fn)(void*, void*, void*, int);


struct tmpi_datatype_component
{
    struct tmpi_datatype_ *type;
    unsigned int count;
};

/* we don't support datatypes with holes (yet)  */
struct tmpi_datatype_
{
    size_t size; /* full extent of type. */   
    tMPI_Op_fn *op_functions; /* array of op functions for this datatype */
    int N_comp; /* number of components */
    struct tmpi_datatype_component *comps; /* the components */
    bool committed; /* whether the data type is committed */
};
/* just as a shorthand:  */
typedef struct tmpi_datatype_ tmpi_dt;



/* global variables */

/* the threads themselves (tmpi_comm only contains lists of pointers to this
         structure */
extern struct tmpi_thread *threads;
extern int Nthreads;

/* thread info */
extern tMPI_Thread_key_t id_key; /* the key to get the thread id */

/* misc. global information about MPI */
extern struct tmpi_global *tmpi_global;






/* error-checking malloc/realloc: */
void *tMPI_Malloc(size_t size);
void *tMPI_Realloc(void *p, size_t size);


/* get the current thread structure pointer */
#define tMPI_Get_current() ((struct tmpi_thread*)tMPI_Thread_getspecific(id_key))

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
bool tMPI_Is_master(void);
/* check whether the current process is in a group */
bool tMPI_In_group(tMPI_Group group);

/* find the rank of a thread in a comm */
int tMPI_Comm_seek_rank(tMPI_Comm comm, struct tmpi_thread *th);
/* find the size of a comm */
int tMPI_Comm_N(tMPI_Comm comm);

/* allocate a comm object, making space for N threads */
tMPI_Comm tMPI_Comm_alloc(tMPI_Comm parent, int N);
/* allocate a group object */
tMPI_Group tMPI_Group_alloc(void);





/* initialize a free envelope list with N envelopes */
void tMPI_Free_env_list_init(struct free_envelope_list *evl, int N);
/* destroy a free envelope list */
void tMPI_Free_env_list_destroy(struct free_envelope_list *evl);


/* initialize a send envelope list */
void tMPI_Send_env_list_init(struct send_envelope_list *evl);
/* destroy a send envelope list */
void tMPI_Send_env_list_destroy(struct send_envelope_list *evl);






/* initialize a recv envelope list */
void tMPI_Recv_env_list_init(struct recv_envelope_list *evl);
/* destroy a recv envelope list */
void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl);




/* initialize request list */
void tMPI_Req_list_init(struct req_list *rl, int N_reqs);
/* destroy request list */
void tMPI_Req_list_destroy(struct req_list *rl);




/* multicast functions */
/* initialize a multi sync structure */
void tMPI_Multi_sync_init(struct multi_sync *msc, int N);

/* destroy a multi sync structure */
void tMPI_Multi_sync_destroy(struct multi_sync *msc);



/* and we need this prototype */
int main(int argc, char **argv);






