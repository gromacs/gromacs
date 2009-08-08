/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-s
tyle: "stroustrup"; -*- 
*
* $Id: 
* 
* This file is part of Gromacs        Copyright (c) 1991-2009
* David van der Spoel, Erik Lindahl, University of Groningen.
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* To help us fund GROMACS development, we humbly ask that you cite
* the research papers on the package. Check out http://www.gromacs.org
* 
* And Hey:
* Gnomes, ROck Monsters And Chili Sauce
*/

/* Include the defines that determine which thread library to use.  
      */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#ifndef THREAD_MPI_STANDALONE
#include "main.h"
#include "statutil.h"
#include "ctype.h"
/*#include "gmx_fatal.h"*/

#else

#define GMX_THREAD_MPI  1
#ifndef cplusplus__
typedef int bool;
#define TRUE 1
#define FALSE 0
#else
#define TRUE true
#define FALSE false
#endif
#endif

#ifdef GMX_THREAD_MPI

#include "smalloc.h"
#include "gmx_thread.h"
#include "gmx_atomic.h"
#include "thread_mpi.h"





/*#define TMPI_DEBUG*/

/* if this is defined, MPI will warn/hang/crash on practices that don't conform
   to the MPI standard (such as not calling MPI_Comm_free on all threads that
   are part of the comm being freed). */
#define TMPI_STRICT 



/* whether to warn if there are mallocs at critical sections (due to 
    preallocations being too small) */
#define TMPI_WARN_MALLOC




/* the max. number of envelopes waiting to be read at thread */
#define MAX_ENVELOPES 128

/* the normal maximum number of threads for pre-defined arrays
   (if the actual number of threads is bigger than this, it'll
    allocate/deallocate arrays, so no problems will arise).*/
#define MAX_PREALLOC_THREADS 32


/* whether to use lock/wait free lists using compare-and-swap (cmpxchg on x86)
   pointer functions. */
#define TMPI_LOCK_FREE_LISTS


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


/* the send envelope. Held in mpi_thread->evs[src_thread]  */
struct send_envelope
{
    int tag; /* the tag */
    /* transmission type flags */
    struct 
    {
        unsigned int nonblock : 1; /* whether the sender is non-blocking */ 
        unsigned int multicast :1; /* whether the message is multicast
                                      (instead of point-to-point) */
    } type;
    MPI_Comm comm; /* this is a structure shared across threads, so we
                      can test easily whether two threads are talking
                      about the same comm. */

    struct mpi_thread *src, *dest; /* these are pretty obvious */

    void *buf; /* buffer to be sent  */
    size_t bufsize; /* the size of the data to be transmitted */
    MPI_Datatype datatype; /* the data type */

     /* state, values from enum_envelope_state .  
        this is volatile because it's probably read by a thread 
        while another one is writing to them (there's a few busy-waits 
        relying on these flags). 
        status=env_unmatched  is the initial state.*/
    gmx_atomic_t state;

    /* the error condition */
    int error;

    /* the message status */
    MPI_Status *status;
    /* prev and next envelopes in the linked list  */
    struct send_envelope *prev,*next;
#ifdef TMPI_LOCK_FREE_LISTS
    /* next pointer for lock-free lists */
    gmx_atomic_ptr_t next_a; 
#else
    /* next pointer for shared lists */
    volatile struct send_envelope *next_a; 
#endif
    /* the list I'm in */
    struct send_envelope_list *list;
};

/* the receive envelope. Held in mpi_thread->evl  */
struct recv_envelope
{
    int tag; /* the tag */
    /* transmission type flags */
    struct 
    {
        unsigned int nonblock : 1; /* whether the receiver is non-blocking */ 
        unsigned int multicast :1; /* whether the message is multicast
                                      (instead of point-to-point)*/
    } type;
    MPI_Comm comm; /* this is a structure shared across threads, so we
                      can test easily whether two threads are talking
                      about the same comm. */

    struct mpi_thread *src, *dest; /* these are pretty obvious */

    void *buf; /* buffer to be sent  */
    size_t bufsize; /* the size of the data to be transmitted, or the
                       size of the actually transmitted data after Xfer. */
    MPI_Datatype datatype; /* the data type */

     /* state, values from enum_envelope_state .  
        this is volatile because it's probably read by a thread 
        while another one is writing to them (there's a few busy-waits 
        relying on these flags). 
        status=env_unmatched  is the initial state.*/
    gmx_atomic_t state;

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
    gmx_atomic_ptr_t head_new; /* singly linked list with the new send 
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

    gmx_spinlock_t lock; /* this locks head_new */

    struct send_envelope old_dummy; /* the dummy element for the head_old
                                       list */
};

struct recv_envelope_list
{
    struct recv_envelope *head; /* first envelope in this list */
    struct recv_envelope dummy; /* the dummy element for the list */
};


/* the request object for asynchronious operations. */
struct mpi_req_
{
    bool recv; /* whether it's a receive request. */
    bool finished; /* whether it's finished */

    struct send_envelope *evs; /* the envelope */
    struct recv_envelope *evr; /* the envelope */

    struct mpi_status_ st;

    struct mpi_req_ *next,*prev; /* next,prev request in linked list */
};

/* pre-allocated  request object list */
struct req_list
{
    struct mpi_req_ *head; /* pre-allocated singly linked list of requests. 
                              (i.e. reqs->prev is undefined). */
    struct mpi_req_ *alloc_head; /* the allocated block */
};








/* MULTICAST COMMUNICATION DATA STRUCTURES */


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
    gmx_atomic_t current_counter; /* counter value for the current
                                     communication */
    gmx_atomic_t n_remaining; /* remaining count for root thread */


    int tag; /* multicast communication type */
    volatile void **buf; /* array of send/recv buffer values */
    volatile size_t *bufsize; /* array of number of bytes to send/recv */
    bool *read_data; /* whether we read data from a specific thread*/
    MPI_Datatype datatype;
};

/* multicast synchronization data structure. There's one of these for 
   each thread in each MPI_Comm structure */
struct multi_sync
{
    int counter;
#define N_MULTI_SYNC 2
    struct multi_env mev[N_MULTI_SYNC];
};








/* GLOBALLY AVAILABLE DATA STRUCTURES */

/* information about a running thread. This structure is put in a 
   globally available array; the envelope exchange, etc. are all done through
   the elements of this array.*/
struct mpi_thread
{
    gmx_thread_t thread_id;

    /* the receive envelopes posted for other threads to check */
    struct recv_envelope_list evr;

    /* the send envelopes posted by other threadas */
    struct send_envelope_list *evs;

    /* change indicator */
    gmx_atomic_t evs_check_id;

    struct req_list rql;  /* the requests */

    /* the free envelopes */
    struct free_envelope_list envelopes; 

    MPI_Comm self_comm;

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
struct mpi_errhandler_
{
    int err;
    MPI_Errhandler_fn fn;
};

/* standard error handler functions */
static void mpi_errors_are_fatal_fn(MPI_Comm *comm, int *err);
static void mpi_errors_return_fn(MPI_Comm *comm, int *err);


/* global MPI information */
struct mpi_global
{
    /* list of pointers to all user-defined types */
    struct mpi_datatype_ **usertypes;
    int N_usertypes;
    int Nalloc_usertypes;

    /* spinlock/mutex for manipulating mpi_user_types */
    gmx_spinlock_t  datatype_lock;
};



















/* COMMUNICATOR DATA STRUCTURES */



/* TODO: make this reference counting */
struct mpi_group_
{
    int N; /* the number of threads */
    struct mpi_thread **peers; /* the list of peers to communicate with */
#if 0
    int Nrefs; /* the number of references to this structure */
#endif
};


/* the communicator objects are globally shared. */
struct mpi_comm_
{
    struct mpi_group_ grp; /* the communicator group */

    /* list of barrier_t's. 
       multicast_barrier[0] contains a barrier for N threads 
       (N=the number of threads in the communicator)
       multicast_barrier[1] contains a barrier for N/2 threads
       multicast_barrier[2] contains a barrier for N/4 threads
       multicast_barrier[3] contains a barrier for N/8 threads
       and so on. (until N/x reaches 1)
       This is to facilitate tree-based algorithms for MPI_Reduce, etc.  */
    gmx_spinlock_barrier_t *multicast_barrier;   
    int *N_multicast_barrier;
    int Nbarriers;


    struct multi_sync *msc; /* list of multicast sync objecs */

    /* lists of globally shared send/receive buffers for MPI_Reduce, etc. */
    volatile void **sendbuf, **recvbuf; 
    
    /* mutex for communication object creation. Traditional mutexes are 
       better here because communicator creation should not be done in 
       time-critical sections of code.   */ 
    gmx_thread_mutex_t comm_create_lock;
    gmx_thread_cond_t comm_create_prep; 
    gmx_thread_cond_t comm_create_finish;

    MPI_Comm *new_comm; /* newly created communicators */
    struct mpi_split *split;

    /* the topologies (only cartesian topology is currently implemented */
    struct cart_topol *cart;
    /*struct mpi_graph_topol_ *graph;*/

    MPI_Errhandler erh;
};


/* specific for MPI_Split: */
struct mpi_split
{ 
    int Ncol_init;
    int Ncol_destroy;
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
struct mpi_graph_topol_
{
};
#endif







/* MPI_Reduce Op functions */
typedef void (*tMPI_Op_fn)(void*, void*, void*, int);


struct mpi_datatype_component
{
    struct mpi_datatype_ *type;
    unsigned int count;
};

/* we don't support datatypes with holes (yet)  */
struct mpi_datatype_
{
    size_t size; /* full extent of type. */   
    tMPI_Op_fn *op_functions; /* array of op functions for this datatype */
    int N_comp; /* number of components */
    struct mpi_datatype_component *comps; /* the components */
    bool committed; /* whether the data type is committed */
};
/* just as a shorthand:  */
typedef struct mpi_datatype_ mpi_dt_;



/* error messages */
static const char *mpi_errmsg[] = 
{
    "No error",
    "MPI Initialization error",
    "MPI Finalize error",
    "Invalid MPI_Group",
    "Invalid MPI_Comm",
    "Invalid MPI_Status",
    "Invalid MPI_Group rank",
    "Invalid Cartesian topology dimensions",
    "Invalid Cartesian topology coordinates",
    "Insufficient number processes for Cartesian topology",
    "Invalid counterpart for MPI transfer",
    "Receive buffer size too small for transmission",
    "Overlapping send/receive buffers: probably due to thread-unsafe code.",
    "Invalid send destination",
    "Invalid receive source",
    "Invalid buffer (null pointer in send or receive buffer)",
    "Multicast operation mismatch (multicast not collective across comm)",
    "Invalid reduce operator",
    "Out of receive envelopes: this shouldn't happen (probably a bug).",
    "Out of receive requests: this shouldn't happen (probably a bug).",
    "Transmission failure",
    "Unknown MPI error"
};



/* get the current thread structure pointer */
static struct mpi_thread *tMPI_Get_current(void);

/* handle an error, returning the errorcode */
static int tMPI_Error(MPI_Comm comm, int mpi_errno);


/* check whether we're the main thread */
static bool tMPI_Is_master(void);
/* check whether the current process is in a group */
static bool tMPI_In_group(MPI_Group group);

/* find the rank of a thread in a comm */
static int tMPI_Comm_seek_rank(MPI_Comm comm, struct mpi_thread *th);
/* find the size of a comm */
static int tMPI_Comm_N(MPI_Comm comm);

/* start N threads with argc, argv (used by MPI_Init)*/
static void tMPI_Start_threads(int N, int *argc, char ***argv, 
                               void (*start_fn)(void*), void *start_arg);
/* starter function for threads; takes a void pointer to a
   struct mpi_starter_, which calls main() if mpi_start_.fn == NULL */
static void* tMPI_Thread_starter(void *arg);

/* allocate a comm object, making space for N threads */
static MPI_Comm tMPI_Comm_alloc(MPI_Comm parent, int N);
/* allocate a group object */
static MPI_Group tMPI_Group_alloc(void);

/* helper function for MPI_Comm_split. Splits N entities with color and key
   out so that the output contains Ngroups groups each with elements
   of the same color. The group array contains the entities in each group. */
static void tMPI_split_colors(int N, const int *color, const int *key, 
                              int *Ngroups, int *grp_N, int *grp_color, 
                              int *group);


/* get the thread belonging to comm with rank rank */
static struct mpi_thread *tMPI_Get_thread(MPI_Comm comm, int rank);


/* initialize a free envelope list with N envelopes */
static void tMPI_Free_env_list_init(struct free_envelope_list *evl, int N);
/* destroy a free envelope list */
static void tMPI_Free_env_list_destroy(struct free_envelope_list *evl);


/* get a free send envelope from the top of a list */
static struct send_envelope *tMPI_Free_env_list_fetch_send(
                                      struct free_envelope_list *evl);
static struct recv_envelope *tMPI_Free_env_list_fetch_recv(
                                      struct free_envelope_list *evl);

/* return an envelope to the free envelopes list */
static void tMPI_Free_env_list_return_send(struct free_envelope_list *evl,
                                           struct send_envelope *ev);
static void tMPI_Free_env_list_return_recv(struct free_envelope_list *evl,
                                           struct recv_envelope *ev);




/* initialize a send envelope list */
static void tMPI_Send_env_list_init(struct send_envelope_list *evl);
/* destroy a send envelope list */
static void tMPI_Send_env_list_destroy(struct send_envelope_list *evl);

/* remove a send envelope from its list. Does not lock */
static void tMPI_Send_env_list_remove(struct send_envelope *ev);

/* add a send envelope to the new envelopes queue in a list */
static void tMPI_Send_env_list_add_new(struct send_envelope_list *evl,
                                       struct send_envelope *ev);
/* move a send envelope to the old envelopes queue in a list. Does not
   lock. */
static void tMPI_Send_env_list_move_to_old(struct send_envelope *ev);







/* initialize a recv envelope list */
static void tMPI_Recv_env_list_init(struct recv_envelope_list *evl);
/* destroy a recv envelope list */
static void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl);

/* add a receive envelope to a list */
static void tMPI_Recv_env_list_add(struct recv_envelope_list *evl,
                                           struct recv_envelope *ev);
/* remove a receive envelope from its list */
static void tMPI_Recv_env_list_remove(struct recv_envelope *ev);





/* do the actual transfer */
static void tMPI_Xfer(struct send_envelope *evs, struct recv_envelope *evr);


/* initialize request list */
static void tMPI_Req_list_init(struct req_list *rl, int N_reqs);
/* destroy request list */
static void tMPI_Req_list_destroy(struct req_list *rl);
/* get a request from the thread's pre-allocated request list */
static struct mpi_req_ *tMPI_Get_req(struct req_list *rl);
/* return a request to the thread's pre-allocated request list */
static void tMPI_Return_req(struct req_list *rl, struct mpi_req_ *req);







/* these are the internal versions of MPI_Isend/Irecv, etc that use 
   pre-allocated mpi_reqs, and can be used for multicast messages. 
   They are used by the regular MPI_Isend/Irecv, etc.  */
static int tMPI_Recv_r(void* buf, int count, MPI_Datatype datatype, 
                      int source, int tag, MPI_Comm comm, 
                      MPI_Status *status, bool multicast);
static int tMPI_Isend_r(void* buf, int count, MPI_Datatype datatype, 
                        int dest, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, bool multicast);
static int tMPI_Irecv_r(void* buf, int count, MPI_Datatype datatype, 
                        int source, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, bool multicast);
static int tMPI_Test_r(struct mpi_req_ *rq, int *flag, MPI_Status *status);
static int tMPI_Wait_r(struct mpi_req_ *rq, MPI_Status *status);
/* wait all on an array of pointers to requests, that progressively get 
   assigned NULL once they're done. If may_free, the pointers also
   get free()d */
static int tMPI_Waitall_r(int count, struct mpi_req_ *array_of_requests[],
                          MPI_Status *array_of_statuses, bool may_free);


/* multicast functions */
/* initialize a multi sync structure */
void tMPI_Multi_sync_init(struct multi_sync *msc, int N);

/* destroy a multi sync structure */
void tMPI_Multi_sync_destroy(struct multi_sync *msc);

/* do a mulitcast transfer, with checks */
int tMPI_Mult_xfer(MPI_Comm comm, int rank, struct multi_env *rmev,
                   void *recvbuf, size_t recvsize, MPI_Datatype recvtype,
                   int expected_tag, int *ret);




/* run a single binary reduce operation on src_a and src_b, producing dest. 
   dest and src_a may be identical */
static int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
                        MPI_Datatype datatype, int count, MPI_Op op, 
                        MPI_Comm comm);


/* and we need this prototype */
int main(int argc, char **argv);







struct mpi_errhandler_ mpi_errors_are_fatal = { 0, mpi_errors_are_fatal_fn };
struct mpi_errhandler_ mpi_errors_return = { 0, mpi_errors_return_fn };


MPI_Errhandler MPI_ERRORS_ARE_FATAL=&mpi_errors_are_fatal;
MPI_Errhandler MPI_ERRORS_RETURN=&mpi_errors_return;




/* there are a few global variables that maintain information about the
   running threads. Some are defined by the MPI standard: */
MPI_Comm MPI_COMM_WORLD=NULL;
MPI_Group MPI_GROUP_EMPTY=NULL;

/* the threads themselves (mpi_comm only contains lists of pointers to this
      structure */
static struct mpi_thread *threads=NULL;
static int Nthreads=0;

/* thread info */
static gmx_thread_key_t id_key; /* the key to get the thread id */

/* whether MPI has finalized (we need this to distinguish pre-inited from
       post-finalized states */
static bool mpi_finalized=FALSE;

/* misc. global information about MPI */
static struct mpi_global *mpi_global=NULL;

/* this is where all the MPI_Reduce ops are included from thread_mpi_ops.c */
#define THREAD_MPI_OPS 1

#define TYPE char
#define TYPENM CHAR
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE short
#define TYPENM SHORT
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE int
#define TYPENM INT
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE long
#define TYPENM LONG
#define INTTYPE 1
#include "thread_mpi_ops.c"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE long long
#define TYPENM L_LONG
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE long long int
#define TYPENM L_L_INT
#define INTTYPE 1
#include "thread_mpi_ops.c"

#endif

#define TYPE signed char
#define TYPENM S_CHAR
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE unsigned char
#define TYPENM U_CHAR
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE unsigned short
#define TYPENM U_SHORT
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE unsigned 
#define TYPENM UNSIGNED
#define INTTYPE 1
#include "thread_mpi_ops.c"

#define TYPE unsigned long
#define TYPENM U_LONG
#define INTTYPE 1
#include "thread_mpi_ops.c"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE unsigned long long
#define TYPENM U_L_LONG
#define INTTYPE 1
#include "thread_mpi_ops.c"

#endif

#define TYPE float
#define TYPENM FLOAT
#define INTTYPE 0
#include "thread_mpi_ops.c"

#define TYPE double
#define TYPENM DOUBLE
#define INTTYPE 0
#include "thread_mpi_ops.c"

#define TYPE long double
#define TYPENM L_DOUBLE
#define INTTYPE 0
#include "thread_mpi_ops.c"

#define TYPE char
#define TYPENM BYTE
#define INTTYPE 1
#include "thread_mpi_ops.c"


mpi_dt_ mpi_char_    ={sizeof(char),              oplist_CHAR,     0,NULL,TRUE};
mpi_dt_ mpi_short_   ={sizeof(short),             oplist_SHORT,    0,NULL,TRUE};
mpi_dt_ mpi_int_     ={sizeof(int),               oplist_INT,      0,NULL,TRUE};
mpi_dt_ mpi_long_    ={sizeof(long),              oplist_LONG,     0,NULL,TRUE};
#ifdef SIZEOF_LONG_LONG_INT
mpi_dt_ mpi_l_long_  ={sizeof(long long),         oplist_L_LONG,   0,NULL,TRUE};
mpi_dt_ mpi_l_l_int_ ={sizeof(long long int),     oplist_L_L_INT,  0,NULL,TRUE};
#endif
mpi_dt_ mpi_s_char_  ={sizeof(signed char),       oplist_S_CHAR,   0,NULL,TRUE};
mpi_dt_ mpi_u_char_  ={sizeof(unsigned char),     oplist_U_CHAR,   0,NULL,TRUE};
mpi_dt_ mpi_u_short_ ={sizeof(unsigned short),    oplist_U_SHORT,  0,NULL,TRUE};
mpi_dt_ mpi_unsigned_={sizeof(unsigned),          oplist_UNSIGNED, 0,NULL,TRUE};
mpi_dt_ mpi_u_long_  ={sizeof(unsigned long),     oplist_U_LONG,   0,NULL,TRUE};
#ifdef SIZEOF_LONG_LONG_INT
mpi_dt_ mpi_u_l_long_={sizeof(unsigned long long),oplist_U_L_LONG, 0,NULL,TRUE};
#endif
mpi_dt_ mpi_float_   ={sizeof(float),             oplist_FLOAT,    0,NULL,TRUE};
mpi_dt_ mpi_double_  ={sizeof(double),            oplist_DOUBLE,   0,NULL,TRUE};
mpi_dt_ mpi_l_double_={sizeof(long double),       oplist_L_DOUBLE, 0,NULL,TRUE};
mpi_dt_ mpi_byte_    ={sizeof(char),              oplist_CHAR,     0,NULL,TRUE};




MPI_Datatype MPI_CHAR               = &mpi_char_;
MPI_Datatype MPI_SHORT              = &mpi_short_;
MPI_Datatype MPI_INT                = &mpi_int_;
MPI_Datatype MPI_LONG               = &mpi_long_;
#ifdef SIZEOF_LONG_LONG_INT
MPI_Datatype MPI_LONG_LONG          = &mpi_l_long_;
MPI_Datatype MPI_LONG_LONG_INT      = &mpi_l_l_int_;
#endif
MPI_Datatype MPI_SIGNED_CHAR        = &mpi_s_char_;
MPI_Datatype MPI_UNSIGNED_CHAR      = &mpi_u_char_;
MPI_Datatype MPI_UNSIGNED_SHORT     = &mpi_u_short_;
MPI_Datatype MPI_UNSIGNED           = &mpi_unsigned_;
MPI_Datatype MPI_UNSIGNED_LONG      = &mpi_u_long_;
#ifdef SIZEOF_LONG_LONG_INT
MPI_Datatype MPI_UNSIGNED_LONG_LONG = &mpi_u_l_long_;
#endif

MPI_Datatype MPI_FLOAT              = &mpi_float_;
MPI_Datatype MPI_DOUBLE             = &mpi_double_;
MPI_Datatype MPI_LONG_DOUBLE        = &mpi_l_double_;

/*extern MPI_Datatype MPI_UNSIGNED_WCHAR*/
MPI_Datatype MPI_BYTE               = &mpi_byte_;



struct mpi_thread *tMPI_Get_current(void)
{
    if (!threads)
        return NULL;

    return (struct mpi_thread*)gmx_thread_getspecific(id_key);
} 


unsigned int tMPI_Threadnr(struct mpi_thread *thr)
{
    return thr-threads;
}

unsigned int tMPI_This_threadnr(void)
{
    return tMPI_Get_current()-threads;
}

static struct mpi_thread *tMPI_Get_thread(MPI_Comm comm, int rank)
{
    /* check destination */
    if ( (rank < 0) || (rank > comm->grp.N) )
    {
        tMPI_Error(comm, MPI_ERR_GROUP_RANK);
        return NULL;
    }
    return comm->grp.peers[rank];
}


static bool tMPI_Is_master(void)
{
    /* if there are no other threads, we're the main thread */
    if ( (!MPI_COMM_WORLD) || MPI_COMM_WORLD->grp.N==0)
        return TRUE;

    /* otherwise we know this through thread specific data: */
    /* whether the thread pointer points to the head of the threads array */
    return (bool)(tMPI_Get_current() == threads); 
}

MPI_Comm tMPI_Get_comm_self(void)
{
    struct mpi_thread* th=tMPI_Get_current();
    return th->self_comm;
}

static int tMPI_Error(MPI_Comm comm, int mpi_errno)
{
    if (comm)
    {
        comm->erh->err=mpi_errno;
        comm->erh->fn(&comm, &mpi_errno);
    }
    else
    {
        /* initialization errors have no comm */
        mpi_errors_are_fatal_fn(NULL, &mpi_errno);
    }
    return mpi_errno;
}

int MPI_Error_string(int errorcode, char *strn, int *resultlen)
{
    if (errorcode<0 || errorcode>=N_MPI_ERR)
        errorcode=MPI_ERR_UNKNOWN;

    strncpy(strn, mpi_errmsg[errorcode], MPI_MAX_ERROR_STRING);
    *resultlen=strlen(strn);
    return MPI_SUCCESS;
}

int MPI_Create_errhandler(MPI_Errhandler_fn *function, 
                          MPI_Errhandler *errhandler) 
{
    smalloc(*errhandler, sizeof(struct mpi_errhandler_));
    (*errhandler)->err=0;
    (*errhandler)->fn=*function;
    return MPI_SUCCESS;
}

int MPI_Errhandler_free(MPI_Errhandler *errhandler)
{
    sfree(*errhandler);
    return MPI_SUCCESS;
}


int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler)
{
    comm->erh = errhandler;
    return MPI_SUCCESS;
}

int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler)
{
    *errhandler=comm->erh;
    return MPI_SUCCESS;
}

static void mpi_errors_are_fatal_fn(MPI_Comm *comm, int *err)
{
    char errstr[MPI_MAX_ERROR_STRING];
    int len;

    MPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "MPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "MPI error: %s\n", errstr);
    }
    abort();
    /*exit(0);*/
}

static void mpi_errors_return_fn(MPI_Comm *comm, int *err)
{
    char errstr[MPI_MAX_ERROR_STRING];
    int len;

    MPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "MPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "MPI error: %s\n", errstr);
    }
    return;
}

int tMPI_Get_N(int *argc, char ***argv)
{
    int i=0;
    int np=0;

    for(i=1;i<*argc;i++)
    {
        if (strcmp("-np", (*argv)[i]) == 0)
        {
            if (i+1 < (*argc))
            {
                /* the number of processes is an argument */
                char *end;
                np=strtol((*argv)[i+1], &end, 0);
                if ( !end || (*end != 0) )
                    np=0;
            }
            break;
        }
    }
    if (np<1)
        np=1;
    return np;
}


static void* tMPI_Thread_starter(void *arg)
{
    struct mpi_thread *th=(struct mpi_thread*)arg;
    int N_envelopes=(Nthreads+1)*(Nthreads+1)*8;  /*AARGH arbitrary number*/
    int N_reqs=(Nthreads+1)*(Nthreads+1)*2;  /*AARGH arbitrary number*/
    int i;

    gmx_thread_setspecific(id_key, arg);


    /* allocate comm.self */
    th->self_comm=tMPI_Comm_alloc(MPI_COMM_WORLD, 1);
    th->self_comm->grp.peers[0]=th;

    /* allocate envelopes */
    tMPI_Free_env_list_init( &(th->envelopes), N_envelopes );
    /* recv list */
    tMPI_Recv_env_list_init( &(th->evr));
    /* send lists */
    smalloc(th->evs, sizeof(struct send_envelope_list)*Nthreads);
    for(i=0;i<Nthreads;i++)
    {
        tMPI_Send_env_list_init( &(th->evs[i]));
    }
    gmx_atomic_set( &(th->evs_check_id), 0);

    /* allocate requests */
    tMPI_Req_list_init(&(th->rql), N_reqs);

    /* now wait for all other threads to come on line, before we
       start the MPI program */
    MPI_Barrier(MPI_COMM_WORLD);

    if (! th->start_fn )
        main(th->argc, th->argv);
    else
    {
        th->start_fn(th->start_arg);
        if (!mpi_finalized)
            MPI_Finalize();
    }

    return 0;
}


void tMPI_Start_threads(int N, int *argc, char ***argv, 
                        void (*start_fn)(void*), void *start_arg)
{
    if (N>0) 
    {
        int i;

        mpi_finalized=FALSE;
        Nthreads=N;

        /* allocate global data */
        smalloc(mpi_global, sizeof(struct mpi_global));
        mpi_global->usertypes=NULL;
        mpi_global->N_usertypes=0;
        mpi_global->Nalloc_usertypes=0;
        gmx_spinlock_init(&(mpi_global->datatype_lock));

        /* allocate world and thread data */
        smalloc(threads, sizeof(struct mpi_thread)*N);
        MPI_COMM_WORLD=tMPI_Comm_alloc(NULL, N);
        MPI_GROUP_EMPTY=tMPI_Group_alloc();

        MPI_COMM_WORLD->grp.N=N;

        if (gmx_thread_key_create(&id_key, NULL))
        {
            tMPI_Error(MPI_COMM_WORLD, MPI_ERR_INIT);
        }
        /*printf("thread keys created\n"); fflush(NULL);*/
        for(i=0;i<N;i++)
        {
            MPI_COMM_WORLD->grp.peers[i]=&(threads[i]);

            /* copy argc, argv */
            if (argc && argv)
            {
                int j;
                threads[i].argc=*argc;
                smalloc(threads[i].argv, threads[i].argc*sizeof(char*));
                for(j=0;j<threads[i].argc;j++)
                {
                    threads[i].argv[j]=strdup( (*argv)[j] );
                }
            }
            else
            {
                threads[i].argc=0;
                threads[i].argv=NULL;
            }
            threads[i].start_fn=start_fn;
            threads[i].start_arg=start_arg;
        }
        /* set the main thread specific item */
        gmx_thread_setspecific(id_key, (void*)&(threads[0]));
        threads[0].thread_id=NULL;

        for(i=1;i<N;i++) /* zero is the main thread */
        {
            if (gmx_thread_create(&(threads[i].thread_id), 
                                  tMPI_Thread_starter,
                                  (void*)&(threads[i]) ) )
            {
                tMPI_Error(MPI_COMM_WORLD, MPI_ERR_INIT);
            }
        }
        /* the main thread now also runs start_fn */
        threads[0].thread_id=NULL;
        tMPI_Thread_starter((void*)&(threads[0]));
    }
}


int MPI_Init(int *argc, char ***argv)
{
    if (MPI_COMM_WORLD==0) /* we're the main process */
    {
        int N=tMPI_Get_N(argc, argv);
        tMPI_Start_threads(N, argc, argv, NULL, NULL);
    }
    else
    {
        /* if we're a sub-thread we need don't need to do anyhing, because 
           everything has already been set up by either the main thread, 
           or the thread runner function.*/
    }
    return MPI_SUCCESS;
}

int tMPI_Init_fn(int N, void (*start_function)(void*), void *arg)
{
    if (MPI_COMM_WORLD==0 && N>=1) /* we're the main process */
    {
        tMPI_Start_threads(N, 0, 0, start_function, arg);
    }
    return MPI_SUCCESS;
}

int MPI_Initialized(int *flag)
{
    *flag=(MPI_COMM_WORLD && !mpi_finalized);

    return MPI_SUCCESS;
}

int MPI_Finalize(void)
{
    int i;
    struct mpi_thread *th=tMPI_Get_current();

#ifdef TMPI_DEBUG
    printf("%5d: MPI_Finalize called\n", tMPI_This_threadnr());
    fflush(stdout);
#endif

    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0;i<(Nthreads+1);i++)
    {
        tMPI_Recv_env_list_destroy( &(th->evr));
    }
    for(i=0;i<Nthreads;i++)
    {
        tMPI_Send_env_list_destroy( &(th->evs[i]));
    }
    tMPI_Free_env_list_destroy( &(th->envelopes) );

    tMPI_Req_list_destroy( &(th->rql) );


    if (tMPI_Is_master())
    {
        /* we just wait for all threads to finish; the order isn't very 
           relevant, as all threads should arrive at their endpoints soon. */
        for(i=1;i<Nthreads;i++)
        {
            if (gmx_thread_join(threads[i].thread_id, NULL))
            {
                tMPI_Error(MPI_COMM_WORLD, MPI_ERR_FINALIZE);
            }
        }
        sfree(threads);
        sfree(MPI_COMM_WORLD);
        sfree(MPI_GROUP_EMPTY);
        threads=0;
        MPI_COMM_WORLD=NULL;
        MPI_GROUP_EMPTY=NULL;
        Nthreads=0;
        mpi_finalized=TRUE;
    }
    else
    {
        gmx_thread_exit(0);
    }
    return MPI_SUCCESS;
}

int MPI_Finalized(int *flag)
{
    *flag=mpi_finalized;

    return MPI_SUCCESS;
}



int MPI_Abort(MPI_Comm comm, int errorcode)
{
    /* we abort(). This way we can run a debugger on it */
    fprintf(stderr, "MPI_Abort called with error code %d",errorcode);
    if (comm==MPI_COMM_WORLD)
        fprintf(stderr, " on MPI_COMM_WORLD");
    fprintf(stderr,"\n");
    fflush(0);

    abort();
#if 0
    /* we just kill all threads, but not the main process */
    int i;
    struct mpi_thread *me=tMPI_Get_current();
    /* kill all threads */
    for(i=0;i<comm->grp.N;i++)
    {
        if (comm->grp.peers[i] != me && threads[i].thread_id)
            gmx_thread_cancel(threads[i].thread_id);
    }
    /* kill myself */
    if (me->thread_id)
        gmx_thread_cancel(me->thread_id);
#endif
    return MPI_SUCCESS;
}


int MPI_Get_processor_name(char *name, int *resultlen)
{
    unsigned int nr=tMPI_Threadnr(tMPI_Get_current());
    unsigned int digits=0;
    const unsigned int base=10;

    /* we don't want to call sprintf here (it turns out to be not entirely
       thread-safe on Mac OS X, for example), so we do it our own way: */

    /* first determine number of digits */
    {
        unsigned int rest=nr;
        while(rest > 0)
        {
            rest /= base;
            digits++;
        }
        if (digits==0)
            digits=1;
    }
    strcpy(name, "thread #");
    /* now construct the number */
    {
        unsigned int len=strlen(name);
        unsigned int i;
        unsigned int rest=nr;

        for(i=0;i<digits;i++)
        {
            int pos=len + (digits-i-1);
            if (pos < (MPI_MAX_PROCESSOR_NAME -1) )
                name[ pos ]=(char)('0' + rest%base);
            rest /= base;
        }
        if ( (digits+len) < MPI_MAX_PROCESSOR_NAME)
            name[digits + len]='\0';
        else
            name[MPI_MAX_PROCESSOR_NAME]='\0';

    }
    if (resultlen)
        *resultlen=strlen(name);
    return MPI_SUCCESS;
}

/* TODO: there must be better ways to do this */
double MPI_Wtime(void)
{
    struct timeval tv;
    double ret;
    gettimeofday(&tv, NULL);
    ret=tv.tv_sec + 1000000.*tv.tv_usec;
    return ret;
}

#if 0
double MPI_Wtick(void)
{
    /* this is just lying: */
    return 1./10000.;
}
#endif



int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count)
{
    if (!status)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_STATUS);
    }
    *count = status->transferred/datatype->size;
    return MPI_SUCCESS;
}



int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype)
{
    struct mpi_datatype_ *ntp;

    smalloc(ntp, sizeof(struct mpi_datatype_));
    ntp->size=count*oldtype->size;
    ntp->op_functions=NULL;

    /* establish components */
    ntp->N_comp=1;
    smalloc(ntp->comps, sizeof(struct mpi_datatype_component)*1);
    ntp->comps[0].type=oldtype;
    ntp->comps[0].count=1;
    ntp->committed=FALSE;

    /* now add it to the list.  */
    gmx_spinlock_lock(&(mpi_global->datatype_lock));
    /* check whether there's space */
    if (mpi_global->N_usertypes + 1 >= mpi_global->Nalloc_usertypes)
    {
        /* make space */
        mpi_global->Nalloc_usertypes=Nthreads*(mpi_global->N_usertypes) + 1;
        srealloc(mpi_global->usertypes, 
                 (sizeof(struct mpi_datatype_ *)*mpi_global->Nalloc_usertypes));
    }
    /* add to the list */
    mpi_global->usertypes[mpi_global->N_usertypes]=ntp;
    mpi_global->N_usertypes++;
    *newtype=ntp;
    gmx_spinlock_unlock(&(mpi_global->datatype_lock));

    return MPI_SUCCESS;
}


int MPI_Type_commit(MPI_Datatype *datatype)
{
    int i,j;
    struct mpi_datatype_ *dt=*datatype;

    if (dt->committed)
        return MPI_SUCCESS;

    /* search the list for a matching committed type, because if there's
       already a committed type that has the same composition, we just 
       make the datatype pointer point to it, ensuring we share datatype 
       information across threads. */
    gmx_spinlock_lock(&(mpi_global->datatype_lock));
    for(i=0;i<mpi_global->N_usertypes;i++)
    {
        struct mpi_datatype_ *lt=mpi_global->usertypes[i];
        if (lt->committed && lt->N_comp==dt->N_comp)
        {
            bool found=TRUE;
            for(j=0;j<lt->N_comp;j++)
            {
                if ( (lt->comps[j].type  != dt->comps[j].type) ||
                     (lt->comps[j].count != dt->comps[j].count) )
                {
                    found=FALSE;
                    break;
                }
            }
            if (found)
            {
                dt=lt;
            }
        }
    }
    if (dt != *datatype)
    {
        bool found=FALSE;
        /* we remove the old one from the list */
        for(i=0;i<mpi_global->N_usertypes;i++)
        {
            if (mpi_global->usertypes[i]==*datatype)
            {
                found=TRUE;
                break;
            }
        }
        if (found)
        {
            /* we put the last one in the list in our slot */
            mpi_global->usertypes[i]=mpi_global->
                usertypes[mpi_global->N_usertypes-1];
            mpi_global->N_usertypes--;
        }
        sfree( (*datatype)->comps);
        sfree(  *datatype );

        /* and overwrite the pointer with the new data type */
        *datatype=dt;
    }
    else
    {
        /* it was the first one of its type */
        dt->committed=TRUE;
    }
    gmx_spinlock_unlock(&(mpi_global->datatype_lock));
    return MPI_SUCCESS;
}


/* Group query & manipulation functions */

static bool tMPI_In_group(MPI_Group group)
{
    int i;
    struct mpi_thread *cur;

    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
            return TRUE;
    }
    return FALSE;
}

int MPI_Group_size(MPI_Group group, int *size)
{
    if (group)
        *size = group->N;
    else
        *size = 0;
    return MPI_SUCCESS;
}

int MPI_Group_rank(MPI_Group group, int *rank)
{    
    int i;
    struct mpi_thread *cur;

    if (!group)
        return MPI_UNDEFINED;

    /* search for my id in the list of peers */
    cur=tMPI_Get_current();
    for(i=0;i<group->N;i++)
    {
        if (group->peers[i] == cur)
        {
            *rank=i;
            return MPI_SUCCESS;
        }
    }
    return MPI_UNDEFINED;
}



static MPI_Group tMPI_Group_alloc(void)
{
    struct mpi_group_ *ret;

    smalloc(ret, sizeof(struct mpi_group_));
    smalloc(ret->peers, sizeof(struct mpi_thread*)*Nthreads);
    ret->N=0;
#if 0
    ret->Nrefs=1;
#endif

    return ret;
}

int MPI_Group_free(MPI_Group *group)
{
    if (group)
    {
        sfree((*group)->peers);
        sfree(*group);
    }
    return MPI_SUCCESS;
}

int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)
{
    int i;
    struct mpi_group_ *ret=tMPI_Group_alloc();

    ret->N=comm->grp.N;
    for(i=0;i<comm->grp.N;i++)
    {
        ret->peers[i]=comm->grp.peers[i];
    }
    *group=ret;
#if 0
    if (comm)
    {
        *group=&(comm->grp);
    }
    else
    {
        *group=NULL;
    }
#endif

    return MPI_SUCCESS;
}


int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup)
{
    int i;
    MPI_Group ng;

    /* just allocate and copy */
    ng=tMPI_Group_alloc();
    ng->N=n;
    for(i=0;i<n;i++)
    {
        if (ranks[i] < 0 || !group || ranks[i] >= group->N)
        {
            return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_GROUP_RANK);
        }
        ng->peers[i]=group->peers[ranks[i]];
    }
    *newgroup=ng;
    return MPI_SUCCESS;
}




/* communicator query&manipulation functions */
static int tMPI_Comm_N(MPI_Comm comm)
{
    if (!comm)
        return 0;
    return comm->grp.N;
}

int MPI_Comm_size(MPI_Comm comm, int *size)
{
    return MPI_Group_size(&(comm->grp), size);
}

int MPI_Comm_rank(MPI_Comm comm, int *rank)
{
    return MPI_Group_rank(&(comm->grp), rank);
}

static MPI_Comm tMPI_Comm_alloc(MPI_Comm parent, int N)
{
    struct mpi_comm_ *ret;
    int i,Nbarriers,Nred;

    smalloc(ret, sizeof(struct mpi_comm_));
    smalloc(ret->grp.peers, sizeof(struct mpi_thread*)*Nthreads);
    ret->grp.N=N;

    gmx_thread_mutex_init( &(ret->comm_create_lock) );
    gmx_thread_cond_init( &(ret->comm_create_prep) );
    gmx_thread_cond_init( &(ret->comm_create_finish) );

    ret->split = NULL;
    ret->new_comm = NULL;
    /* we have no topology to start out with */
    ret->cart=NULL;
    /*ret->graph=NULL;*/

    /* calculate the number of multicast barriers */
    Nbarriers=0;
    Nred=N;
    while(Nred>1) {
        Nbarriers+=1;
        Nred = Nred/2 + Nred%2;
    } 

    ret->Nbarriers=Nbarriers;
    smalloc(ret->multicast_barrier, 
            sizeof(gmx_spinlock_barrier_t)*(Nbarriers+1));
    smalloc(ret->N_multicast_barrier, sizeof(int)*(Nbarriers+1));
    Nred=N;
    for(i=0;i<Nbarriers;i++)
    {
        gmx_spinlock_barrier_init( &(ret->multicast_barrier[i]), Nred);
        ret->N_multicast_barrier[i]=Nred;
        /* Nred is now Nred/2 + a rest term because solitary 
           process at the end of the list must still be accounter for */
        Nred = Nred/2 + Nred%2;
    }
    smalloc(ret->sendbuf, sizeof(void*)*Nthreads);
    smalloc(ret->recvbuf, sizeof(void*)*Nthreads);


    if (parent)
    {
        ret->erh=parent->erh;
    }
    else
    {
        ret->erh=MPI_ERRORS_ARE_FATAL;
    }

    /* multi_sync objects */
    smalloc( ret->msc, sizeof(struct multi_sync)*N);
    for(i=0;i<N;i++)
        tMPI_Multi_sync_init( &(ret->msc[i]), N);

    return ret;
}

int MPI_Comm_free(MPI_Comm *comm)
{
#ifndef TMPI_STRICT
    int myrank=tMPI_Comm_seek_rank(*comm, tMPI_Get_current());
    if (! *comm)
        return MPI_SUCCESS;

    if ((*comm)->grp.N > 1)
    {
        /* we remove ourselves from the comm. */
        gmx_thread_mutex_lock(&((*comm)->comm_create_lock));
        (*comm)->grp.peers[myrank] = (*comm)->grp.peers[(*comm)->grp.N-1];
        (*comm)->grp.N--;
        gmx_thread_mutex_unlock(&((*comm)->comm_create_lock));
    }
    else
    {
        /* we're the last one so we can safely destroy it */
        sfree((*comm)->grp.peers);
        sfree((*comm)->multicast_barrier);
        sfree((*comm)->sendbuf);
        sfree((*comm)->recvbuf);
        if ( (*comm)->cart)
        {
            sfree((*comm)->cart->dims);
            sfree((*comm)->cart->periods);
            sfree((*comm)->cart);
        }
        sfree(*comm);
    }
#else
    int myrank=tMPI_Comm_seek_rank(*comm, tMPI_Get_current());
    /* This is correct if programs actually treat Comm_free as a 
       collective call */
    /* we need to barrier because the comm is a shared structure and
       we have to be sure that nobody else is using it 
       (for example, to get its rank, like above) before destroying it*/
    MPI_Barrier(*comm);
    /* this is a collective call on a shared data structure, so only 
       one process (rank[0] in this case) should do anything */
    if (myrank==0)
    {
        sfree((*comm)->grp.peers);
        sfree((*comm)->multicast_barrier);
        sfree((*comm)->sendbuf);
        sfree((*comm)->recvbuf);
        if ( (*comm)->cart)
        {
            sfree((*comm)->cart->dims);
            sfree((*comm)->cart->periods);
            sfree((*comm)->cart);
        }
        sfree(*comm);
    }
#endif
    return MPI_SUCCESS;
}

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    /* we just call Comm_split because it already contains all the
       neccesary synchronization constructs. */
    return MPI_Comm_split(comm, 0, tMPI_Comm_seek_rank(comm, 
                            tMPI_Get_current()), newcomm);
}


int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    int color=MPI_UNDEFINED;
    int key=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    if (tMPI_In_group(group))
    {
        color=1;
    }
    /* the MPI specs specifically say that this is equivalent */
    return MPI_Comm_split(comm, color, key, newcomm);
}

static void tMPI_split_colors(int N, const int *color, const int *key, 
                              int *Ngroups, int *grp_N, int *grp_color, 
                              int *group)
{
    int i,j;
    bool found;

    /* reset groups */
    for(i=0;i<N;i++)
        grp_N[i]=0;
    for(i=0;i<N;i++)
    {
        if (color[i] != MPI_UNDEFINED)
        {
            found=FALSE;
            for(j=0;j<(*Ngroups);j++)
            {
                if (grp_color[j] == color[i])
                {
                    /* we insert where we need to, by counting back */
                    int k=grp_N[j];

                    while (k>0 && ( key[group[N*j + k-1]]>key[i]) )
                    {
                        /* shift up */
                        group[N*j + k]=group[N*j + k-1];
                        k--;
                    }
                    group[N*j+k]=i;
                    grp_N[j]++;
                    found=TRUE;
                }
            }
            if (!found)
            {
                /* not found. just add a new color */
                grp_N[(*Ngroups)]=1;
                grp_color[(*Ngroups)]=color[i];
                group[N*(*Ngroups) + 0]=i;
                (*Ngroups)++;
            }
        }
    }
}

/* this is the main comm creation function. All other functions that create
    comms use this*/
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    int i,j;
    int N=tMPI_Comm_N(comm);
    volatile MPI_Comm *newcomm_list;
    volatile int colors[MAX_PREALLOC_THREADS]; /* array with the colors 
                                                  of each thread */
    volatile int keys[MAX_PREALLOC_THREADS]; /* same for keys (only one of 
                                                the threads actually suplies 
                                                these arrays to the comm 
                                                structure) */
    bool i_am_first=FALSE;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    struct mpi_split *spl;

    if (!comm)
    {
        *newcomm=NULL;
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    /*printf("**Calling MPI_Comm_split with color %d, key %d\n",color, key);*/

    gmx_thread_mutex_lock(&(comm->comm_create_lock));    
    /* first get the colors */
    if (!comm->new_comm)
    {
        smalloc(comm->split, sizeof(struct mpi_split));
        smalloc(comm->new_comm, N*sizeof(MPI_Comm));
        if (N<=MAX_PREALLOC_THREADS)
        {
            comm->split->colors=colors;
            comm->split->keys=keys;
        }
        else
        {
            smalloc(comm->split->colors, N*sizeof(int));
            smalloc(comm->split->keys, N*sizeof(int));
        }
        comm->split->Ncol_init=tMPI_Comm_N(comm); 
        i_am_first=TRUE;
        /* the main communicator contains a list the size of grp.N */
    }
    newcomm_list=comm->new_comm; /* we copy it to the local stacks because
                                    we can later erase comm->new_comm safely */
    spl=comm->split; /* we do the same for spl */
    spl->colors[myrank] = color;
    spl->keys[myrank] = key;
    spl->Ncol_init--;

    if (spl->Ncol_init == 0)
        gmx_thread_cond_signal(&(comm->comm_create_prep));

    if (!i_am_first)
    {
        /* all other threads can just wait until the creator thread is 
           finished */
        gmx_thread_cond_wait( &(comm->comm_create_finish) ,
                              &(comm->comm_create_lock) );
    }
    else
    {
        int Ncomms=0;
        int comm_color_[MAX_PREALLOC_THREADS]; 
        int comm_N_[MAX_PREALLOC_THREADS]; 
        int *comm_color=comm_color_; /* there can't be more comms than N*/
        int *comm_N=comm_N_; /* the number of procs in a group */

        int *comm_groups; /* the groups */
        MPI_Comm *comms; /* the communicators */

        /* wait for the colors to be done */
        if (N>1)
        {
            gmx_thread_cond_wait( &(comm->comm_create_prep), 
                                  &(comm->comm_create_lock));
        }

        /* reset the state so that a new comm creating function can run */
        spl->Ncol_destroy=N;
        comm->new_comm=0;
        comm->split=0;

        smalloc(comm_groups, N*N*sizeof(int));
        if (N>MAX_PREALLOC_THREADS)
        {
            smalloc(comm_color, N*sizeof(int));
            smalloc(comm_N, N*sizeof(int));
        }

        /* count colors, allocate and split up communicators */
        tMPI_split_colors(N, (int*)spl->colors, 
                             (int*)spl->keys, 
                             &Ncomms, 
                             comm_N, comm_color, comm_groups);


        /* allocate a bunch of communicators */
        smalloc(comms, Ncomms*sizeof(MPI_Comm));
        for(i=0;i<Ncomms;i++)
            comms[i]=tMPI_Comm_alloc(comm, comm_N[i]);

        /* now distribute the comms */
        for(i=0;i<Ncomms;i++)
        {
            comms[i]->grp.N=comm_N[i];
            for(j=0;j<comm_N[i];j++)
                comms[i]->grp.peers[j]=
                    comm->grp.peers[comm_groups[i*comm->grp.N + j]];
        }
        /* and put them into the newcomm_list */
        for(i=0;i<N;i++)
        {
            newcomm_list[i]=MPI_COMM_NULL;
            for(j=0;j<Ncomms;j++)
            {
                if (spl->colors[i] == comm_color[j])
                {
                    newcomm_list[i] = comms[j];
                    break;
                }
            }
        }

#ifdef TMPI_DEBUG
        /* output */
        for(i=0;i<Ncomms;i++)
        {
            printf("Group %d (color %d) has %d members: ",
                    i, comm_color[i], comm_N[i]);
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",comm_groups[comm->grp.N*i + j]);

            printf(" rank: ");
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",spl->keys[comm_groups[N*i + j]]);
            printf(" color: ");
            for(j=0;j<comm_N[i];j++)
                printf(" %d ",spl->colors[comm_groups[N*i + j]]);
            printf("\n");
        }
#endif
        if (N>MAX_PREALLOC_THREADS)
        {
            sfree((int*)spl->colors);
            sfree((int*)spl->keys);
            sfree(comm_color);
            sfree(comm_N);
        }
        sfree(comm_groups);

        /* tell the waiting threads that there's a comm ready */
        gmx_thread_cond_broadcast(&(comm->comm_create_finish));
    }
    /* here the individual threads get their comm object */
    *newcomm=newcomm_list[myrank];

    /* free when we have assigned them all, so we can reuse the object*/
    spl->Ncol_destroy--;
    if (spl->Ncol_destroy==0)
    {
        sfree((void*)newcomm_list);
        sfree(spl);
    }

    gmx_thread_mutex_unlock(&(comm->comm_create_lock));

    return MPI_SUCCESS;    
}

static int tMPI_Comm_seek_rank(MPI_Comm comm, struct mpi_thread *th)
{
    int i;
    if (!comm)
        return -1;

    for(i=0;i<comm->grp.N;i++)
    {
        if (comm->grp.peers[i] == th)
            return i;
    }
    return -1;
}



/* topology functions */
int MPI_Topo_test(MPI_Comm comm, int *status)
{
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }

    if (comm->cart)
        *status=MPI_CART;
    /*else if (comm->graph)
        status=MPI_GRAPH;*/
    else 
        *status=MPI_UNDEFINED;

    return MPI_SUCCESS;
}

int MPI_Cartdim_get(MPI_Comm comm, int *ndims)
{
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
    {
        return MPI_SUCCESS;
    }
    *ndims=comm->cart->ndims;
    return MPI_SUCCESS;
}


int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods,
                 int *coords)
{
    int i;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return MPI_SUCCESS;

    MPI_Cart_coords(comm, myrank, maxdims, coords);

    for(i=0;i<comm->cart->ndims;i++)
    {
        if (i>=maxdims)
        {
            return tMPI_Error(comm, MPI_ERR_DIMS);
        }
        dims[i]=comm->cart->dims[i];
        periods[i]=comm->cart->periods[i];
    }

    return MPI_SUCCESS;
}

int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank)
{
    int i,mul=1,ret=0;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return MPI_SUCCESS;

    /* because of row-major ordering, we count the dimensions down */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        int rcoord=coords[i];
        if (comm->cart->periods[i])
        {
            /* apply periodic boundary conditions */
            rcoord = rcoord % comm->cart->dims[i];
            if (rcoord < 0)
                rcoord += comm->cart->dims[i];
        }
        else
        {
            if (rcoord < 0 || rcoord >= comm->cart->dims[i])
            {
                return tMPI_Error(comm, MPI_ERR_DIMS);
            }
        }
        ret += mul*rcoord;
        mul *= comm->cart->dims[i];
    }
    *rank=ret;
    return MPI_SUCCESS;
}

int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords)
{
    int i;
    int rank_left=rank;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!comm->cart || comm->cart->ndims==0)
        return MPI_SUCCESS;
    if (maxdims < comm->cart->ndims)
    {
        return tMPI_Error(comm, MPI_ERR_DIMS);
    }

    /* again, row-major ordering */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        coords[i]=rank_left%comm->cart->dims[i];
        rank_left /= comm->cart->dims[i];
    }   

    return 0;
}



int MPI_Cart_map(MPI_Comm comm, int ndims, int *dims, int *periods, 
                 int *newrank)
{
    /* this function doesn't actually do anything beyond returning the current 
       rank (or MPI_UNDEFINED if it doesn't fit in the new topology */
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int Ntot=1;
    int i;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    /* calculate the total number of procs in cartesian comm */
    for(i=0;i<ndims;i++)
    {
        Ntot *= dims[i];
    }

    if (myrank >= Ntot)
    {
        *newrank=MPI_UNDEFINED;
    }
    else
    {
        *newrank=myrank;
    }

    return MPI_SUCCESS;
}



int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
                    int reorder, MPI_Comm *comm_cart)
{
    int myrank=tMPI_Comm_seek_rank(comm_old, tMPI_Get_current());
    int key=myrank;
    int newrank=-1;
    int color=0;
    int Ntot=1;
    int i;
    

    if (!comm_old)
    {
        return tMPI_Error(comm_old, MPI_ERR_COMM);
    }
    /* calculate the total number of procs in cartesian comm */
    for(i=0;i<ndims;i++)
    {
        Ntot *= dims[i];
    }
    /* refuse to create if there's not enough procs */
    if (comm_old->grp.N < Ntot)
    {
        *comm_cart=MPI_COMM_NULL;
#if 1
        return tMPI_Error(comm_old, MPI_ERR_CART_CREATE_NPROCS);
#endif
    }

    if (key >= Ntot)
        key=MPI_UNDEFINED;

    if (reorder)
    {
        MPI_Cart_map(comm_old, ndims, dims, periods, &key);
    }

    if (key==MPI_UNDEFINED)
    {
        color=MPI_UNDEFINED;
    }

    MPI_Comm_split(comm_old, color, key, comm_cart);

    if (*comm_cart)
    {
        MPI_Comm_rank(*comm_cart, &newrank);
    }

    if (newrank==0)
    {
        smalloc((*comm_cart)->cart, sizeof(struct cart_topol));
        smalloc((*comm_cart)->cart->dims, ndims*sizeof(int));
        smalloc((*comm_cart)->cart->periods, ndims*sizeof(int));
        (*comm_cart)->cart->ndims=ndims;
        for(i=0;i<ndims;i++)
        {
            (*comm_cart)->cart->dims[i]=dims[i];
            (*comm_cart)->cart->periods[i]=periods[i];
        }
    }

    /* and we add a barrier to make sure the cart object is seen by 
       every thread that is part of the new communicator */
    if (*comm_cart)
    {
        gmx_spinlock_barrier_wait( &( (*comm_cart)->multicast_barrier[0]) );
    }


    return MPI_SUCCESS;
}

/* Point-to-point communication protocol functions */


static void tMPI_Free_env_list_init(struct free_envelope_list *evl, int N)
{
    int i;

    /* allocate the head element */
    smalloc(evl->send_alloc_head, sizeof(struct send_envelope)*N );
    smalloc(evl->recv_alloc_head, sizeof(struct recv_envelope)*N );
    evl->head_send=evl->send_alloc_head;
    evl->head_recv=evl->recv_alloc_head;

    for(i=0;i<N;i++)
    {
        if (i < N-1)
        {
            evl->head_send[i].next=&(evl->head_send[i+1]);
            evl->head_recv[i].next=&(evl->head_recv[i+1]);
#ifdef TMPI_LOCK_FREE_LISTS
            gmx_atomic_ptr_set(&(evl->head_send[i].next_a),NULL);
#else
            evl->head_send[i].next_a=NULL;
#endif
        }
        else
        {
            evl->head_send[i].next=NULL;
            evl->head_recv[i].next=NULL;
#ifdef TMPI_LOCK_FREE_LISTS
            gmx_atomic_ptr_set(&(evl->head_send[i].next_a),NULL);
#else
            evl->head_send[i].next_a=NULL;
#endif
        }
        evl->head_send[i].list=NULL;
        evl->head_recv[i].list=NULL;
    }
}

static void tMPI_Free_env_list_destroy(struct free_envelope_list *evl)
{
    sfree(evl->send_alloc_head);
    sfree(evl->recv_alloc_head);
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
    gmx_atomic_ptr_set(&(ret->next_a),NULL);
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
    gmx_atomic_ptr_set(&(ev->next_a),NULL);
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






/* mpi_send_envelope_list functions */

static void tMPI_Send_env_list_init(struct send_envelope_list *evl)
{
    gmx_spinlock_init( &(evl->lock) );

#ifdef TMPI_LOCK_FREE_LISTS
    gmx_atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new = NULL;
#endif
    evl->head_old = &(evl->old_dummy);
    evl->head_old->next = evl->head_old;
    evl->head_old->prev = evl->head_old;
}

static void tMPI_Send_env_list_destroy(struct send_envelope_list *evl)
{
#ifdef TMPI_LOCK_FREE_LISTS
    gmx_atomic_ptr_set(&(evl->head_new), NULL);
#else
    evl->head_new=NULL; 
#endif
    evl->head_old=NULL; /* make it crash if used after MPI_Finalize */
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
    gmx_atomic_ptr_set(&(ev->next_a), NULL);
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
                          gmx_atomic_ptr_get( &(evl->head_new) );
        /* set our envelope to have that as its next */
        gmx_atomic_ptr_set( &(ev->next_a), evl_head_new_orig );
        /* do the compare-and-swap */
        evl_cas_result=(struct send_envelope*)gmx_atomic_ptr_cmpxchg( 
                                                &(evl->head_new), 
                                                evl_head_new_orig,
                                                ev);
        /* and compare the results: if they aren't the same,
           somebody else got there before us: */
    } while (evl_cas_result != evl_head_new_orig); 
#else
    gmx_spinlock_lock( &(evl->lock) );
    /* we add to the start of the list */
    ev->next_a=evl->head_new;
    /* actually attach it to the list */
    evl->head_new=ev;
    gmx_spinlock_unlock( &(evl->lock) );
#endif

    /* signal to the thread that there is a new envelope */
    gmx_atomic_fetch_add( &(ev->dest->evs_check_id) ,1);
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








/* mpi_recv_envelope_list functions */

static void tMPI_Recv_env_list_init(struct recv_envelope_list *evl)
{
    evl->head = &(evl->dummy);
    evl->head->prev=evl->head;
    evl->head->next=evl->head;
}

static void tMPI_Recv_env_list_destroy(struct recv_envelope_list *evl)
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








/* mpi_req functions */

static void tMPI_Req_list_init(struct req_list *rl, int N_reqs)
{
    int i;

    smalloc(rl->alloc_head, sizeof(struct mpi_req_)*N_reqs);
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

static void tMPI_Req_list_destroy(struct req_list *rl)
{
    sfree(rl->alloc_head);
    rl->head=NULL;
    rl->alloc_head=NULL;
}

static struct mpi_req_ *tMPI_Get_req(struct req_list *rl)
{
    struct mpi_req_ *req=rl->head;
    

    /* we don't need locks here because requests are a per-thread
       property */
    if (!req)
    {
        /* this could be fixed */
        tMPI_Error(MPI_COMM_WORLD, MPI_ERR_REQUESTS);
        return NULL;
    }
    rl->head=req->next;
    req->next=NULL;

    return req;
}

static void tMPI_Return_req(struct req_list *rl, struct mpi_req_ *req)
{
    req->next=rl->head;
    req->prev=NULL;
    rl->head=req;
}






/* Point-to-point communication protocol functions */



static void tMPI_Set_recv_status(struct recv_envelope *ev, MPI_Status *status)
{
    if (status)
    {
        status->MPI_SOURCE = tMPI_Comm_seek_rank(ev->comm, ev->src);
        status->MPI_TAG = ev->tag;
        status->MPI_ERROR = ev->error;
        if (gmx_atomic_get(&(ev->state))==env_finished)
            status->transferred = ev->bufsize/ev->datatype->size;
        else
            status->transferred = 0;
    }
}

static void tMPI_Set_send_status(struct send_envelope *ev, MPI_Status *status)
{
    if (status)
    {
        status->MPI_SOURCE = tMPI_Comm_seek_rank(ev->comm, ev->src);
        status->MPI_TAG = ev->tag;
        status->MPI_ERROR = ev->error;
        if (gmx_atomic_get(&(ev->state))==env_finished)
            status->transferred = ev->bufsize/ev->datatype->size;
        else
            status->transferred = 0;
    }
}



static bool tMPI_Envelope_matches(const struct send_envelope *send,
                                  const struct recv_envelope *recv)
{
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Envelope_matches (%d->%d)==(%d->%d),  tag=(%d==%d),       \n       xfertype=(%d==%d), datatype=(%ld==%ld), comm=(%ld,%ld),\n              finished=(%d==%d)\n",
            tMPI_This_threadnr(),
            tMPI_Threadnr(send->src), tMPI_Threadnr(send->dest),
            tMPI_Threadnr(recv->src), tMPI_Threadnr(recv->dest),
            (int)(send->tag), (int)(recv->tag),
            (send->type.multicast), (recv->type.multicast),
            (long int)send->datatype, (long int)recv->datatype,
            (long int)send->comm, (long int)recv->comm,
            (int)send->state.value, (int)recv->state.value);
    fflush(stdout);
#endif
    if ( ( (recv->tag == MPI_ANY_TAG) || (recv->tag == send->tag) ) &&
            ( send->type.multicast == recv->type.multicast ) &&
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




static struct send_envelope *tMPI_Prep_send_envelope(struct mpi_thread *cur,
                                                     MPI_Comm comm, 
                                                     struct mpi_thread *src, 
                                                     struct mpi_thread *dest, 
                                                     void *buf, int count, 
                                                     MPI_Datatype datatype, 
                                                     int tag, bool multicast, 
                                                     bool nonblock)
{
    /* get an envelope from the stack */
    struct send_envelope *ev=tMPI_Free_env_list_fetch_send( &(cur->envelopes) );

    ev->tag=tag;
    ev->type.multicast=multicast;
    ev->type.nonblock=nonblock;

    ev->comm=comm;

    ev->src=src;
    ev->dest=dest;

    ev->buf=buf;
    ev->bufsize=count*datatype->size;
    ev->datatype=datatype;

    ev->list=NULL;

    gmx_atomic_set(&(ev->state), env_unmatched);
    ev->error=MPI_SUCCESS;

    return ev;
}

static struct recv_envelope *tMPI_Prep_recv_envelope(struct mpi_thread *cur,
                                                     MPI_Comm comm, 
                                                     struct mpi_thread *src, 
                                                     struct mpi_thread *dest, 
                                                     void *buf, int count, 
                                                     MPI_Datatype datatype, 
                                                     int tag, bool multicast, 
                                                     bool nonblock)
{
    /* get an envelope from the stack */
    struct recv_envelope *ev=tMPI_Free_env_list_fetch_recv( &(cur->envelopes) );

    ev->tag=tag;
    ev->type.multicast=multicast;
    ev->type.nonblock=nonblock;

    ev->comm=comm;

    ev->src=src;
    ev->dest=dest;

    ev->buf=buf;
    ev->bufsize=count*datatype->size;
    ev->datatype=datatype;
    
    ev->list=NULL;

    gmx_atomic_set(&(ev->state), env_unmatched);
    ev->error=MPI_SUCCESS;

    return ev;
}









static struct recv_envelope *tMPI_Post_match_recv(MPI_Comm comm, 
                                           struct mpi_thread *src, 
                                           void *recv_buf, int recv_count,
                                           MPI_Datatype datatype, int tag, 
                                           bool multicast, bool nonblock)
{
    struct mpi_thread *cur=tMPI_Get_current();
    struct mpi_thread *dest=cur;
    struct recv_envelope *evr;
    struct send_envelope *evs=NULL;
    int src_threadnr=src ? tMPI_Threadnr(src) : Nthreads;
    int i;

    /* reserve an envelope to post */
    evr=tMPI_Prep_recv_envelope(cur, comm, src, dest, recv_buf, recv_count, 
                                datatype, tag, multicast, nonblock);

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




static struct send_envelope *tMPI_Post_send(MPI_Comm comm, 
                                            struct mpi_thread *dest, 
                                            void *send_buf, int send_count,
                                            MPI_Datatype datatype, int tag, 
                                            bool multicast, bool nonblock)
{
    struct mpi_thread *cur=tMPI_Get_current();
    struct mpi_thread *src=cur;
    struct send_envelope *evs;
    int src_threadnr=tMPI_Threadnr(src);

    /* reserve an envelope to post */
    evs=tMPI_Prep_send_envelope(cur, comm, src, dest, send_buf, send_count, 
                                datatype, tag, multicast, nonblock);

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




static void tMPI_Test_incoming(struct mpi_thread *th)
{
    int i;
    int check_id;

    /* we check for newly arrived send envelopes */
    check_id=gmx_atomic_get( &(th->evs_check_id));
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
                        gmx_atomic_ptr_get( &(th->evs[i].head_new) );
                /* do the compare-and-swap to detach the list */
                evl_cas_result=(struct send_envelope*)gmx_atomic_ptr_cmpxchg( 
                                                        &(th->evs[i].head_new),
                                                        evs_head,
                                                        NULL);
            } while (evl_cas_result != evs_head);
#else
            gmx_spinlock_lock( &(th->evs[i].lock) );
            evs_head=(struct send_envelope*)th->evs[i].head_new;
            th->evs[i].head_new=NULL; /* detach the list */
            gmx_spinlock_unlock( &(th->evs[i].lock) );
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
                               gmx_atomic_ptr_get( &(evs->next_a) );
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
        check_id=gmx_atomic_fetch_add( &(th->evs_check_id), -n);
    }
}








static bool tMPI_Test_recv(struct recv_envelope *ev, bool blocking, 
                           MPI_Status *status)
{
    struct mpi_thread *cur=tMPI_Get_current();

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
        while( gmx_atomic_get( &(ev->state) ) <  env_finished ) 
        {
            /* while blocking, we wait for incoming send envelopes */
            tMPI_Test_incoming(cur);
        }
    }
    else
    {
        if ( gmx_atomic_get( &(ev->state) ) <  env_finished ) 
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
                           MPI_Status *status)
{
    struct mpi_thread *cur=tMPI_Get_current();

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
        while( gmx_atomic_get( &(ev->state) ) <  env_finished ) 
        {
            /* while blocking, we wait for incoming send envelopes. That's
               all we can do at this moment. */
            tMPI_Test_incoming(cur);
        }
    }
    else
    {
        if ( gmx_atomic_get( &(ev->state) ) <  env_finished ) 
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
            tMPI_Error((evr->comm), MPI_ERR_XFER_BUFSIZE);
            gmx_atomic_set(&(evr->state), env_finished);
            gmx_atomic_set(&(evs->state), env_finished);
            evr->error = MPI_ERR_XFER_BUFSIZE;
            evs->error = MPI_ERR_XFER_BUFSIZE;
            return;
        }

        if (!evr->buf || !evs->buf)
        {
            tMPI_Error((evr->comm), MPI_ERR_BUF);
            gmx_atomic_set(&(evr->state), env_finished);
            gmx_atomic_set(&(evs->state), env_finished);
            evr->error = MPI_ERR_BUF;
            evs->error = MPI_ERR_BUF;
            return;
        }
        memcpy(evr->buf, evs->buf, evs->bufsize);
        /* for status update */
    }
    evr->bufsize=evs->bufsize;
    /* and mark that we're finished */
    gmx_atomic_set( &(evr->state), env_finished);
    gmx_atomic_set( &(evs->state), env_finished);
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


static int tMPI_Send_r(void* buf, int count, MPI_Datatype datatype, int dest,
                      int tag, MPI_Comm comm, bool multicast)
{
    struct send_envelope *sd;
    struct mpi_thread *send_dst=tMPI_Get_thread(comm, dest);

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!send_dst)
    {
        return tMPI_Error(comm, MPI_ERR_SEND_DEST);
    }

    sd=tMPI_Post_send(comm, send_dst, buf, count, datatype, tag, 
                      multicast, FALSE);
    tMPI_Test_send(sd, TRUE, NULL);

    return sd->error;
}

static int tMPI_Recv_r(void* buf, int count, MPI_Datatype datatype, 
                      int source, int tag, MPI_Comm comm, 
                      MPI_Status *status, bool multicast)
{
    struct recv_envelope *rc;
    struct mpi_thread *recv_src=0;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }

    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, MPI_ERR_RECV_SRC); 
        }
    }

    rc=tMPI_Post_match_recv(comm, recv_src, buf, count, datatype, tag, 
                            multicast, FALSE);
    tMPI_Test_recv(rc, TRUE, status);

    return rc->error;
}


int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 int dest, int sendtag, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype, int source, int recvtag, 
                 MPI_Comm comm, MPI_Status *status)
{
    struct recv_envelope *rc;
    struct send_envelope *sd;

    struct mpi_thread *recv_src=0;
    struct mpi_thread *send_dst=tMPI_Get_thread(comm, dest);
    bool send_finished=FALSE; 
    bool recv_finished=FALSE;


    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!send_dst)
    {
        return tMPI_Error(comm, MPI_ERR_SEND_DEST); 
    }
    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, MPI_ERR_RECV_SRC);
        }
    }
    /* we first prepare to send */
    sd=tMPI_Post_send(comm, send_dst, sendbuf, sendcount, 
                            sendtype, sendtag, FALSE, FALSE);
    /* then we prepare to receive */
    rc=tMPI_Post_match_recv(comm, recv_src, recvbuf, recvcount, 
                            recvtype, recvtag, FALSE, FALSE);

    do 
    {
        /* now we wait to receive */
        if (!recv_finished)
        {
            recv_finished=tMPI_Test_recv(rc, FALSE, status);
            if (rc->error != MPI_SUCCESS)
                return rc->error;
        }   
        /* we wait until the send completes */
        if (!send_finished)
        {
            send_finished=tMPI_Test_send(sd, FALSE, NULL);
            if (rc->error != MPI_SUCCESS)
                return rc->error;
        }
    }
    while (! (send_finished && recv_finished) );

    return MPI_SUCCESS;
}


/* async */
static int tMPI_Isend_r(void* buf, int count, MPI_Datatype datatype, int dest,
                        int tag, MPI_Comm comm, struct mpi_req_ *rq, 
                        bool multicast)
{
    struct mpi_thread *send_dst=tMPI_Get_thread(comm, dest);

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!send_dst)
    {
        return tMPI_Error(comm, MPI_ERR_SEND_DEST);
    }
    if (!buf && count!=0)
    {
        return tMPI_Error(comm, MPI_ERR_BUF);
    }
    rq->recv=FALSE;
    rq->finished=FALSE;
    rq->evs=tMPI_Post_send(comm, send_dst, buf, count, datatype, tag,
                           multicast, TRUE);
    return MPI_SUCCESS;    
}

static int tMPI_Irecv_r(void* buf, int count, MPI_Datatype datatype, 
                        int source, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, bool multicast)
{
    struct mpi_thread *recv_src=0;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }

    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            return tMPI_Error(comm, MPI_ERR_RECV_SRC);
        }
    }
    rq->recv=TRUE;
    rq->finished=FALSE;
    rq->evr=tMPI_Post_match_recv(comm, recv_src, buf, count, datatype, tag,
                                 multicast, TRUE);
    return MPI_SUCCESS;    
}

static int tMPI_Wait_r(struct mpi_req_ *rq, MPI_Status *status)
{
    int ret=MPI_SUCCESS;
    if (!rq)
        return MPI_SUCCESS;

    if ( ! ( rq->finished || 
           ( 
            ( rq->recv && gmx_atomic_get( &(rq->evr->state)) == env_finished )
                    ||
            ( !rq->recv && gmx_atomic_get( &(rq->evs->state)) == env_finished )
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

static int tMPI_Test_r(struct mpi_req_ *rq, int *flag, MPI_Status *status)
{    
    int ret=MPI_SUCCESS;
    bool finished=FALSE;

    if (!rq)
        return MPI_SUCCESS;

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

static int tMPI_Waitall_r(int count, struct mpi_req_ *array_of_requests[],
                          MPI_Status *array_of_statuses, bool may_free)
{
    int done;
    int flags_[MAX_PREALLOC_THREADS];
    int *flags=flags_;
    int i;
    struct req_list *rql=&(tMPI_Get_current()->rql);

    if (count > MAX_PREALLOC_THREADS)
    {
#ifdef TMPI_WARN_MALLOC
        fprintf(stderr, "Warning: malloc during MPI_Waitall_r\n");
#endif
        smalloc(flags, sizeof(int)*count);
    }

    for(i=0;i<count;i++)
        flags[i]=FALSE;
    /* Waitall polls all the requests by calling MPI_Test. This
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
                struct mpi_req_ *rq=array_of_requests[i];
                MPI_Status *status=NULL;
                if (array_of_statuses)
                    status=&(array_of_statuses[i]);

                flags[i]=tMPI_Test_recv(rq->evr, FALSE, status);

                if (rq->evr->error!=MPI_SUCCESS)
                    return rq->evr->error;
                if (flags[i])
                {
                    rq->evr=NULL;
                    if (may_free)
                    {
                        tMPI_Return_req(rql,rq);
                        array_of_requests[i]=MPI_REQUEST_NULL;
                    }
                }
            }
        }
        /* then  do sends */
        for(i=0;i<count;i++)
        {
            if (!flags[i] && !(array_of_requests[i]->recv))
            {
                struct mpi_req_ *rq=array_of_requests[i];
                MPI_Status *status=NULL;
                if (array_of_statuses)
                    status=&(array_of_statuses[i]);

                flags[i]=tMPI_Test_send(rq->evs, FALSE, status);
                                    
                if (rq->evs->error!=MPI_SUCCESS)
                    return rq->evs->error;
                if (flags[i])
                {
                    rq->evs=NULL;
                    if (may_free)
                    {
                        tMPI_Return_req(rql,rq);
                        array_of_requests[i]=MPI_REQUEST_NULL;
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
        sfree(flags);

    return MPI_SUCCESS;
}



/* the real MPI functions just call their tMPI_*r counterparts */

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm)
{
    return tMPI_Send_r(buf, count, datatype, dest, tag, comm, FALSE);
}


int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
             int tag, MPI_Comm comm, MPI_Status *status)
{
    return tMPI_Recv_r(buf, count, datatype, source, tag, comm, status, FALSE);
}


int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    struct req_list *rql=&(tMPI_Get_current()->rql);
    struct mpi_req_ *rq=tMPI_Get_req(rql);

    if ( (ret=tMPI_Isend_r(buf, count, datatype, dest, tag, comm, rq, FALSE)) 
                != MPI_SUCCESS)
    {
#ifdef TMPI_DEBUG
        printf("%5d Freeing isend request after error\n", 
               tMPI_This_threadnr()); 
        fflush(0);
#endif
        tMPI_Return_req(rql,rq);
        *request=NULL;
    }
    else
    {
        *request=rq;
    }
    return ret;
}

int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request *request)
{
    int ret;
    struct req_list *rql=&(tMPI_Get_current()->rql);
    struct mpi_req_ *rq=tMPI_Get_req(rql);

    if ( (ret=tMPI_Irecv_r(buf, count, datatype, source, tag, comm, rq, FALSE))
            != MPI_SUCCESS)
    {
#ifdef TMPI_DEBUG
        printf("%5d Freeing irecv request after error\n", 
                tMPI_This_threadnr()); 
        fflush(0);
#endif
        tMPI_Return_req(rql,rq);
        *request=NULL;
    }
    else
    {
        *request=rq;
    }
    return ret;
}

int MPI_Wait(MPI_Request *request, MPI_Status *status)
{
    int ret=MPI_SUCCESS;
    struct req_list *rql=&(tMPI_Get_current()->rql);
    if (!request || !(*request))
        return MPI_SUCCESS;

    ret=tMPI_Wait_r(*request, status);

    /* deallocate if needed */
    tMPI_Return_req(rql, *request);
    *request=MPI_REQUEST_NULL;
    return ret;
}

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{
    int ret=MPI_SUCCESS;
    struct req_list *rql=&(tMPI_Get_current()->rql);
    if (!request || !(*request))
        return MPI_SUCCESS;

    ret=tMPI_Test_r(*request, flag, status);

    if ((*request)->finished)
    {
        /* deallocate if needed */
        tMPI_Return_req(rql, *request);
        *request=MPI_REQUEST_NULL;
    }
    return ret;
}


int MPI_Waitall(int count, MPI_Request *array_of_requests,
                MPI_Status *array_of_statuses)
{
    return tMPI_Waitall_r(count, array_of_requests, array_of_statuses, TRUE);
}


void tMPI_Multi_sync_init(struct multi_sync *msc, int N)
{
    int i;

    msc->counter=0;
    for(i=0;i<N_MULTI_SYNC;i++)
    {
        gmx_atomic_set( &(msc->mev[i].current_counter), 0);
        gmx_atomic_set( &(msc->mev[i].n_remaining), 0);
        smalloc(msc->mev[i].buf, sizeof(void*)*N);
        smalloc(msc->mev[i].bufsize, sizeof(size_t)*N);
        smalloc(msc->mev[i].read_data, sizeof(bool)*N);
    }
}

void tMPI_Multi_sync_destroy(struct multi_sync *msc)
{
    int i;

    for(i=0;i<N_MULTI_SYNC;i++)
    {
        sfree((void*)msc->mev[i].buf);
        sfree((void*)msc->mev[i].bufsize);
        sfree((void*)msc->mev[i].read_data);
    }
}

int tMPI_Mult_xfer(MPI_Comm comm, int rank, struct multi_env *rmev, 
                   void *recvbuf, size_t recvsize, MPI_Datatype recvtype, 
                   int expected_tag, int *ret)
{
    int sendsize=rmev->bufsize[rank];
    /* check tags, types */
    if ( (rmev->datatype != recvtype ) || (rmev->tag != expected_tag) )
    {
        return tMPI_Error(comm, MPI_ERR_MULTI_MISMATCH);
    }
  
    if (sendsize) /* we allow NULL ptrs if there's nothing to xmit */
    {
        if ( sendsize > recvsize ) 
        {
            return tMPI_Error(comm, MPI_ERR_XFER_BUFSIZE);
        }

        if ( rmev->buf[rank] == recvbuf )
        {
            return tMPI_Error(MPI_COMM_WORLD,MPI_ERR_XFER_BUF_OVERLAP);
        }
        /* copy data */
        memcpy((char*)recvbuf, (void*)(rmev->buf[rank]), sendsize);
    }
    /* signal one thread ready */
    gmx_atomic_fetch_add( &(rmev->n_remaining), -1);
    return MPI_SUCCESS;
}


int MPI_Barrier(MPI_Comm comm) 
{
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }

    if (comm->grp.N>1)
    {
        gmx_spinlock_barrier_wait( &(comm->multicast_barrier[0]));
    }
    return MPI_SUCCESS;
}


/* multi */
int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root,
              MPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=MPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        struct multi_env *mev=&(msc->mev[mevi]);

        /* first set up the data */
        mev->tag=TMPI_BCAST_TAG;
        mev->datatype=datatype;
        mev->buf[0]=buffer;
        mev->bufsize[0]=count*datatype->size;
        gmx_atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);
        /* and wait until everybody is done copying */
        while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=count*datatype->size;
        /* wait until root becomes available */
        while( gmx_atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }

        tMPI_Mult_xfer(comm, 0, rmev, buffer, bufsize, datatype, 
                       TMPI_BCAST_TAG, &ret);
    }
    return ret;
}





int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=MPI_SUCCESS;
    struct multi_env *mev;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);
    if (myrank==root)
    {
        int i;
        int n_remaining=comm->grp.N-1;
        /* do root transfer */
        {
            int recvsize=recvtype->size*recvcount;
            int sendsize=sendtype->size*sendcount;
            if (recvsize < sendsize)
            {
                return tMPI_Error(comm, MPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, MPI_ERR_MULTI_MISMATCH);
            }
            if ( sendbuf == (char*)recvbuf+myrank*recvcount*recvtype->size )
            {
                return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy((char*)recvbuf+myrank*recvcount*recvtype->size, sendbuf, 
                   sendsize);
        }
        /* poll data availability */
        while(n_remaining>0)
        {
            for(i=0;i<comm->grp.N;i++)
            {
                struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
                if ((i!=myrank) && 
                    (gmx_atomic_get(&(rmev->current_counter))==msc->counter)&& 
                    (gmx_atomic_get(&(rmev->n_remaining)) > 0) )
                {
                    ret=tMPI_Mult_xfer(comm, 0, rmev,
                                       (char*)recvbuf+
                                            i*recvcount*recvtype->size,
                                       recvcount*recvtype->size,
                                       recvtype, TMPI_GATHERV_TAG,&ret);
                    if (ret!=MPI_SUCCESS)
                        return ret;
                    n_remaining--;
                }
            }
        }
    }
    else
    {
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, MPI_ERR_BUF);
        }

        /* first set up the data */
        mev->tag=TMPI_GATHERV_TAG;
        mev->datatype=sendtype;
        mev->buf[0]=sendbuf;
        mev->bufsize[0]=bufsize;

        gmx_atomic_set(&(mev->n_remaining), 1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);

        /* and wait until root is done copying */
        while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }

    return ret;
}


int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int *recvcounts, int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=MPI_SUCCESS;
    struct multi_env *mev;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);
    if (myrank==root)
    {
        int i;
        int n_remaining=comm->grp.N-1;
        /* do root transfer */
        {
            int recvsize=recvtype->size*recvcounts[myrank];
            int sendsize=sendtype->size*sendcount;
            if (recvsize < sendsize)
            {
                return tMPI_Error(comm, MPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, MPI_ERR_MULTI_MISMATCH);
            }
            if ( sendbuf == (char*)recvbuf+displs[myrank]*recvtype->size )
            {
                return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy((char*)recvbuf+displs[myrank]*recvtype->size, 
                   sendbuf, sendsize);
        }
        /* poll data availability */
        while(n_remaining>0)
        {
            for(i=0;i<comm->grp.N;i++)
            {
                struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
                if ((i!=myrank) && 
                    (gmx_atomic_get(&(rmev->current_counter))==msc->counter)&& 
                    (gmx_atomic_get(&(rmev->n_remaining)) > 0) )
                {
                    ret=tMPI_Mult_xfer(comm, 0, rmev,
                                       (char*)recvbuf+displs[i]*recvtype->size,
                                       recvcounts[i]*recvtype->size,
                                       recvtype, TMPI_GATHERV_TAG,&ret);
                    if (ret!=MPI_SUCCESS)
                        return ret;
                    n_remaining--;
                }
            }
        }
    }
    else
    {
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, MPI_ERR_BUF);
        }

        /* first set up the data */
        mev->tag=TMPI_GATHERV_TAG;
        mev->datatype=sendtype;
        mev->buf[0]=sendbuf;
        mev->bufsize[0]=bufsize;

        gmx_atomic_set(&(mev->n_remaining), 1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);

        /* and wait until root is done copying */
        while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }

    return ret;
}


int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, 
                MPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=MPI_SUCCESS;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        int i;
        struct multi_env *mev=&(msc->mev[mevi]);
        size_t bufsize = sendcount*sendtype->size;

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, MPI_ERR_BUF);
        }
        /* first set up the data */
        mev->tag=TMPI_SCATTER_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf + bufsize*i;
            mev->bufsize[i]=bufsize;
        }
        gmx_atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);

        /* do the root transfer */
        {
            int recvsize=recvtype->size*recvcount;
            if (recvsize < bufsize)
            {
                return tMPI_Error(comm, MPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, MPI_ERR_MULTI_MISMATCH);
            }
            if ( recvbuf == mev->buf[myrank] )
            {
                return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy(recvbuf, (void*)(mev->buf[myrank]), bufsize);
        }

        /* and wait until everybody is done copying */
        while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=recvcount*recvtype->size;
        /* wait until root becomes available */
        while( gmx_atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }
        tMPI_Mult_xfer(comm, myrank, rmev, recvbuf, bufsize, recvtype,
                       TMPI_SCATTER_TAG, &ret);
    }
    return ret;
}



int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype sendtype, void* recvbuf, int recvcount,
                 MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    struct multi_sync *msc;
    int mevi;
    int myrank;
    int ret=MPI_SUCCESS;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());


    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    if (myrank==root)
    {
        int i;
        struct multi_env *mev=&(msc->mev[mevi]);

        if (!sendbuf) /* don't do pointer arithmetic on a NULL ptr */
        {
            return tMPI_Error(comm, MPI_ERR_BUF);
        }
        /* first set up the data */
        mev->tag=TMPI_SCATTERV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf + displs[i]*sendtype->size;
            mev->bufsize[i]=sendtype->size*sendcounts[i];
        }
        gmx_atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);

        /* do the root transfer */
        {
            int recvsize=recvtype->size*recvcount;
            if (recvsize < mev->bufsize[myrank])
            {
                return tMPI_Error(comm, MPI_ERR_XFER_BUFSIZE);
            }
            if (recvtype != sendtype)
            {
                return tMPI_Error(comm, MPI_ERR_MULTI_MISMATCH);
            }
            if ( recvbuf == mev->buf[myrank] )
            {
                return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_XFER_BUF_OVERLAP);
            }

            memcpy(recvbuf, (void*)(mev->buf[myrank]), mev->bufsize[myrank]);
        }

        /* and wait until everybody is done copying */
        while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
        {
        }
    }
    else
    {
        /* get the root mev */
        struct multi_env *rmev = &(comm->msc[root].mev[mevi]);
        size_t bufsize=recvcount*recvtype->size;
        /* wait until root becomes available */
        while( gmx_atomic_get( &(rmev->current_counter)) != msc->counter)
        {
        }
        tMPI_Mult_xfer(comm, myrank, rmev, recvbuf, bufsize, recvtype,
                       TMPI_SCATTERV_TAG, &ret);
    }
    return ret;

}



int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm)
{
    struct multi_sync *msc;
    struct multi_env *mev;
    int mevi;
    int i;
    int myrank;
    int n_remaining;
    int ret=MPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, MPI_ERR_BUF);
    }

    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);

    /* post our pointers */
    {
        /* first set up the data */
        mev->tag=TMPI_ALLTOALLV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf+i*sendcount*sendtype->size;
            mev->bufsize[i]=sendcount*sendtype->size;
            mev->read_data[i]=FALSE;
        }

        gmx_atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);
    }
    /* do root transfer */
    {
        int recvsize=recvtype->size*recvcount;
        tMPI_Mult_xfer(comm, myrank, mev,
                       (char*)recvbuf+myrank*recvcount*recvtype->size, 
                       recvsize, recvtype, TMPI_ALLTOALLV_TAG, &ret);
        if (ret!=MPI_SUCCESS)
            return ret;
                    
        mev->read_data[myrank]=TRUE;
    }

    n_remaining=comm->grp.N-1;
    /* poll data availability */
    while(n_remaining>0)
    {
        for(i=0;i<comm->grp.N;i++)
        {
            struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
            if ((!mev->read_data[i]) &&
                (gmx_atomic_get(&(rmev->current_counter))==msc->counter) )
            {
                ret=tMPI_Mult_xfer(comm, myrank, rmev, 
                                  (char*)recvbuf+i*recvcount*recvtype->size,
                                  recvcount*recvtype->size, 
                                  recvtype, TMPI_ALLTOALLV_TAG,&ret);
                if (ret!=MPI_SUCCESS)
                    return ret;

                mev->read_data[i]=TRUE;
                n_remaining--;
            }
        }
    }
    /* and wait until everybody is done copying our data */
    while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
    {
    }

    return ret;
}


int MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype sendtype, 
                  void* recvbuf, int *recvcounts, int *rdispls, 
                  MPI_Datatype recvtype, 
                  MPI_Comm comm)

{
    struct multi_sync *msc;
    struct multi_env *mev;
    int mevi;
    int i;
    int myrank;
    int n_remaining;
    int ret=MPI_SUCCESS;

    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!sendbuf || !recvbuf) /* don't do pointer arithmetic on a NULL ptr */
    {
        return tMPI_Error(comm, MPI_ERR_BUF);
    }

    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    msc=&(comm->msc[myrank]);
    /* we increase our counter. */
    msc->counter++;
    mevi=msc->counter % N_MULTI_SYNC;
    mev=&(msc->mev[mevi]);

    /* post our pointers */
    {
        /* first set up the data */
        mev->tag=TMPI_ALLTOALLV_TAG;
        mev->datatype=sendtype;
        for(i=0;i<comm->grp.N;i++)
        {
            mev->buf[i]=(char*)sendbuf+sdispls[i]*sendtype->size;
            mev->bufsize[i]=sendcounts[i]*sendtype->size;
            mev->read_data[i]=FALSE;
        }

        gmx_atomic_set(&(mev->n_remaining), comm->grp.N-1);
        /* we now publish our counter */
        gmx_atomic_set(&(mev->current_counter), msc->counter);
    }
    /* do root transfer */
    {
        int recvsize=recvtype->size*recvcounts[myrank];
        tMPI_Mult_xfer(comm, myrank, mev,
                       (char*)recvbuf+rdispls[myrank]*recvtype->size, 
                       recvsize, recvtype, TMPI_ALLTOALLV_TAG, &ret);
        if (ret!=MPI_SUCCESS)
            return ret;
                    
        mev->read_data[myrank]=TRUE;
    }

    n_remaining=comm->grp.N-1;
    /* poll data availability */
    while(n_remaining>0)
    {
        for(i=0;i<comm->grp.N;i++)
        {
            struct multi_env *rmev=&(comm->msc[i].mev[mevi]);
            if ((!mev->read_data[i]) &&
                (gmx_atomic_get(&(rmev->current_counter))==msc->counter) )
            {
                ret=tMPI_Mult_xfer(comm, myrank, rmev, 
                                  (char*)recvbuf+rdispls[i]*recvtype->size,
                                  recvcounts[i]*recvtype->size, 
                                  recvtype, TMPI_ALLTOALLV_TAG,&ret);
                if (ret!=MPI_SUCCESS)
                    return ret;

                mev->read_data[i]=TRUE;
                n_remaining--;
            }
        }
    }
    /* and wait until everybody is done copying our data */
    while (gmx_atomic_get( &(mev->n_remaining) ) > 0)
    {
    }

    return ret;
}






int tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
                        MPI_Datatype datatype, int count, MPI_Op op, 
                        MPI_Comm comm)
{
    tMPI_Op_fn fn=datatype->op_functions[op];

    if (src_a==src_b)
    {
        return tMPI_Error(comm, MPI_ERR_XFER_BUF_OVERLAP);
    }
    fn(dest, src_a, src_b, count);
    return MPI_SUCCESS;
}

int tMPI_Reduce_fast(void* sendbuf, void* recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int i;

    /* this function uses a binary tree-like reduction algorithm: */
    int N=tMPI_Comm_N(comm);
    int myrank_rtr=(N+myrank-root)%N; /* my rank relative to root */
    int Nnbrs=N; /* number of neighbours to communicate with 
                    (decreases exponentially) */
    int nbr_dist=1; /* distance between communicating neighbours 
                       (increases exponentially) */
    int stepping=2; /* distance between non-communicating neighbours
                       (increases exponentially) */
    int iteration=0;
  
    if (count==0)
        return MPI_SUCCESS;
    if (!comm)
    {
        return tMPI_Error(MPI_COMM_WORLD, MPI_ERR_COMM);
    }
    if (!recvbuf)
    {
        return tMPI_Error(comm, MPI_ERR_BUF);
    }
    if ( (!datatype->op_functions) || (!datatype->op_functions[op]) )
    {
        return tMPI_Error(comm, MPI_ERR_OP_FN);
    }

    if (!sendbuf)/* i.e. sendbuf == MPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }
    comm->sendbuf[myrank]=sendbuf;
    comm->recvbuf[myrank]=recvbuf;
    /* there's a barrier to wait for all the processes to put their 
       send/recvbuf in the global list */
    gmx_spinlock_barrier_wait( &(comm->multicast_barrier[0]));

    /* check the buffers */
    for(i=0;i<N;i++)
    {
        if ( (i!=myrank) && ( (comm->recvbuf[i]==recvbuf) || 
                              (comm->sendbuf[i]==sendbuf) ) )
        {
            return tMPI_Error(comm, MPI_ERR_XFER_BUF_OVERLAP);
        }
    }

    while(Nnbrs>1)
    {
        int nbr;
        /* reduce between myrank and myrank+nbr_dist, if there is such
            a neighbour. */
#ifdef TMPI_DEBUG
        printf("%d: iteration=%d, Nnbrs=%d, barrier size=%d and I'm still in here\n", 
                myrank, iteration, comm->N_multicast_barrier[iteration], Nnbrs);
        fflush(0);
#endif
        if (myrank_rtr%stepping == 0 )
        {
            if (myrank_rtr+nbr_dist<N) 
            {
                void *a,*b;
                int ret;
#ifdef TMPI_DEBUG
                printf("%d: reducing with %d, iteration=%d\n", 
                        myrank, myrank+nbr_dist, iteration);
                fflush(0);
#endif
               /* we reduce with our neighbour*/
                nbr=(N+myrank+nbr_dist)%N; /* here we use the real rank */
                if (iteration==0)
                {
                    a=sendbuf;
                    b=(void*)(comm->sendbuf[nbr]);
                }
                else
                {
                    a=recvbuf;
                    b=(void*)(comm->recvbuf[nbr]);
                }
                if ((ret=tMPI_Reduce_run_op(recvbuf, a, b, datatype, 
                                            count, op, comm)) != MPI_SUCCESS)
                    return ret;
            }
            else
            {
                /* we still need to put things in the right buffer for the next
                   iteration */
                if (iteration==0)
                    memcpy(recvbuf, sendbuf, datatype->size*count);
            }
            /* split barrier */
            gmx_spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
        }
        else 
        {
            /* this barrier is split because the threads not actually 
               calculating still need to be waiting for the thread using its
               data to finish calculating */
#ifdef TMPI_DEBUG
            printf("%d: barrier waiting, iteration=%d\n", myrank,  iteration);
            fflush(0);
#endif
            gmx_spinlock_barrier_wait( &(comm->multicast_barrier[iteration]));
            break;
        }
#ifdef TMPI_DEBUG
        printf("%d: barrier over, iteration=%d\n", myrank,  iteration);
        fflush(0);
#endif
        Nnbrs = Nnbrs/2 + Nnbrs%2;
        nbr_dist*=2;
        stepping*=2;
        iteration++;
    }
    return MPI_SUCCESS;
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int ret;

    if (myrank==root)
    {
        if (!sendbuf) /* i.e. sendbuf == MPI_IN_PLACE */
        {
            sendbuf=recvbuf;
        }
    }
    else
    {
#ifdef TMPI_WARN_MALLOC
        fprintf(stderr, "Warning: malloc during MPI_Reduce\n");
#endif
        smalloc(recvbuf, datatype->size*count);
    }
    ret=tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, root, comm);
    if (myrank!=root)
    {
        sfree(recvbuf);
    }
    return ret;
}

int MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    void *rootbuf=NULL; /* root process' receive buffer */
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    int ret;

    if (count==0)
        return MPI_SUCCESS;
    if (!recvbuf)
    {
        return tMPI_Error(comm, MPI_ERR_BUF);
    }
    if (!sendbuf) /* i.e. sendbuf == MPI_IN_PLACE */
    {
        sendbuf=recvbuf;
    }

    ret=tMPI_Reduce_fast(sendbuf, recvbuf, count, datatype, op, 0, comm);
    /* distribute rootbuf */
    rootbuf=(void*)comm->recvbuf[0];

    gmx_spinlock_barrier_wait( &(comm->multicast_barrier[0]));
    /* and now we just copy things back inefficiently. We should make
       a better MPI_Scatter, and use that. */
    if (myrank != 0)
    {
        if (rootbuf==recvbuf)
        {
            return tMPI_Error(comm, MPI_ERR_XFER_BUF_OVERLAP);
        }
        memcpy(recvbuf, rootbuf, datatype->size*count );
    }
    gmx_spinlock_barrier_wait( &(comm->multicast_barrier[0]));

    return MPI_SUCCESS;
}



#endif /* GMX_THREAD_MPI */

