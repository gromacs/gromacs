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

#ifndef THREAD_MPI_TESTING
#include "main.h"
#include "statutil.h"
#include "ctype.h"
#endif

#ifdef GMX_THREAD_MPI

#include "gmx_fatal.h"
#include "smalloc.h"
#include "gmx_thread.h"
#include "thread_mpi.h"

/*#define TMPI_DEBUG*/


/* the max. number of envelopes waiting to be read at thread */
#define MAX_ENVELOPES 128

/* the normal maximum number of threads for pre-defined arrays
   (if the actual number of threads is bigger than this, it'll
    allocate/deallocate arrays, so no problems will arise
    if the number of threads is bigger than this) */
#define MAX_THREADS 32


/* information about a running thread. This structure is put in a 
   globally available array.*/
struct mpi_thread_
{
    gmx_thread_t thread_id;

    int N_evr; /* the number of receive envelopes currently waiting */
    /* the receive envelopes pending for this thread: */
    struct mpi_recv_envelope_ *evr[MAX_ENVELOPES]; 

    int N_evs; /* the number of send envelopes currently waiting */
    /* the send envelopes pending for this thread: */
    struct mpi_send_envelope_ *evs[MAX_ENVELOPES]; 


    /* mutex for exchanging envelopes. Either the sending or the receiving
       process put an envelope pointer on the envelope list of the receiving 
       thread. This list is protected by mutexes; if there already is
       a matching envelope, the other process will complete the data
       structure, allowing the receiving end (in nonblocking communication
       also the sending end) to transfer the data. */
    gmx_thread_mutex_t envelope_mutex;
    /* mutex for transmitting data in the case of non-blocking receives.
       If receives are non-blocking, the sender may need to initiate transfer
       so this mutex is used for locking out one or the other party. */
    gmx_thread_mutex_t xmit_mutex;

    MPI_Comm self_comm;

    /* we copy these for each thread (this is not required by the 
       MPI standard, but it's convenient). Note that we copy, because
       some programs (like Gromacs) like to manipulate these. */    
    int argc;
    char **argv;
};


/* TODO: make this reference counting */
struct mpi_group_
{
    int N; /* the number of threads */
    struct mpi_thread_ **peers; /* the list of peers to communicate with */
#if 0
    int Nrefs; /* the number of references to this structure */
#endif
};


/* the communicator objects are globally shared. */
struct mpi_comm_
{
    struct mpi_group_ grp; /* the communicator group */

    /* mutex and condition types for comm syncs */
    gmx_thread_mutex_t multicast_mutex;
    gmx_thread_cond_t multicast_cond;
    /* list of barrier_t's. 
       multicast_barrier[0] contains a barrier for N threads 
       (N=the number of threads in the communicator)
       multicast_barrier[1] contains a barrier for N/2 threads
       multicast_barrier[2] contains a barrier for N/4 threads
       multicast_barrier[3] contains a barrier for N/8 threads
       and so on. (until N/x reaches 1)
       This is to facilitate tree-based algorithms for MPI_Reduce, etc.  */
    gmx_thread_barrier_t *multicast_barrier;   
    int *N_multicast_barrier;
    int Nbarriers;

    /* lists of globally shared send/receive buffers for MPI_Reduce, etc. */
    volatile void **sendbuf, **recvbuf; 
    
    /* mutex for communcation object creation */ 
    gmx_thread_mutex_t comm_create_mutex;
    gmx_thread_cond_t comm_create_prep; 
    gmx_thread_cond_t comm_create_finish;

    MPI_Comm *new_comm; /* newly created communicators */
    struct mpi_split_ *split;

    /* the topologies (only cartesian topology is currently implemented */
    struct mpi_cart_topol_ *cart;
    /*struct mpi_graph_topol_ *graph;*/

    MPI_Errhandler erh;
};


/* specific for MPI_Split: */
struct mpi_split_
{ 
    int Ncol_init;
    int Ncol_destroy;
    volatile int *colors;
    volatile int *keys;
};

/* cartesian topology */
struct mpi_cart_topol_
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

enum xfertype_
{
    point_to_point,
    multicast
};

/* the message envelopes (as described in the MPI standard). 
   These fully describes the message, and make each message unique (enough).
  
   Transmitting data works by either the transmitting or receiving process
   putting a pointer to an envelope on the envelope list (the envelope
   list is protected by mutexes). The other end then puts a pointer in
   the 'counterpart' field. If it's the sending side, it then waits 
   for the receiver to complete the transfer. The receiving side starts
   the transfer as soon as the data is available from the sender. 

   If the receiver is non-blocking, the sender may also initiate transfers,
   and transfers become excluded by mutexes.  */
struct mpi_send_envelope_
{
    MPI_Comm comm;  /* this is a structure shared across threads, so we
                       can test easily whether two threads are talking
                       about the same comm. */

    struct mpi_thread_ *src, *dest; /* these are pretty obvious */
    int tag; /* the tag */
    enum xfertype_ xt; /* the type of transfer (to make sure we're not 
                          mixing up point-to-point messages with 
                          multicast messages)*/

    void *buf; /* buffer to be sent  */
    size_t bufsize; /* the size of the data to be sent  */
    MPI_Datatype datatype; /* the data type */

    /* these are volatile because they're probably read by a thread while
       another one is writing to them (there's a few busy-waits relying on
       these flags). There are also flags for the receiver (which are 
       set by the receiver), that make waiting for conditions a lot easier */
    volatile bool ready; /* whether the data is ready to be sent */
    volatile bool nonblock; /* whether the sender is non-blocking */
    volatile bool finished; /* whether the data transmission is finished */
    volatile bool recv_ready; /* whether the data is ready to be received */
    volatile bool recv_nonblock; /* whether the receiver is non-blocking */

    /* A link to the other envelope (the sending envelope for the receiver,
       or the receiving envelope for the sender). Status updates are
       done for both the envelope itself, and for its counterpart. 
       This way, we can throw away the local data structure if we're
       finished with it, but the other side still has a valid status.

       The counterpart can be relied on to exist until finished is true.  */
    volatile struct mpi_recv_envelope_ *counterpart; 
};
/* the corresponding receive envelope */
struct mpi_recv_envelope_
{
    MPI_Comm comm;  /* this is a structure shared across threads, so we
                       can test easily whether two threads are talking
                       about the same comm. */

    struct mpi_thread_ *src, *dest; /* these are pretty obvious */
    int tag; /* the tag */
    enum xfertype_ xt; /* the type of transfer (to make sure we're not 
                          mixing up point-to-point messages with 
                          multicast messages)*/

    void *buf; /* buffer to be received  */
    size_t bufsize; /* the size of the data to be received  */
    MPI_Datatype datatype; /* the data type */
    size_t transferred; /* the size of data actually transferred */

    /* these are volatile because they're probably read by a thread while
       another one is writing to them (there's a few busy-waits relying on
       these flags). There are also flags for the sender (which are 
       set by the sender), that make waiting for conditions a lot easier. */
    volatile bool ready; /* whether the data is ready to be received */
    volatile bool nonblock; /* whether the receiver is non-blocking */
    volatile bool finished; /* whether the data transmission is finished */
    volatile bool send_ready; /* whether the data is ready to be sent  */
    volatile bool send_nonblock; /* whether the sender is non-blocking */

    /* A link to the other envelope. Status updates are
       done for both the envelope itself, and for its counterpart. 
       This way, we can throw away the local data structure if we're
       finished with it, but the other side still has a valid status.

       The counterpart can be relied on to exist until finished is true.  */
    volatile struct mpi_send_envelope_ *counterpart; 
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

struct mpi_errhandler_ mpi_errors_are_fatal = { 0, mpi_errors_are_fatal_fn };
struct mpi_errhandler_ mpi_errors_return = { 0, mpi_errors_return_fn };

MPI_Errhandler MPI_ERRORS_ARE_FATAL=&mpi_errors_are_fatal;;
MPI_Errhandler MPI_ERRORS_RETURN=&mpi_errors_return;;



/* global MPI information */
struct mpi_global_
{
    /* list of pointers to all user-defined types */
    struct mpi_datatype_ **usertypes;
    int N_usertypes;
    int Nalloc_usertypes;

    /* mutex for manipulating mpi_user_types */
    gmx_thread_mutex_t datatype_mutex;
};


/* MPI_Reduce Op functions */
typedef void (*tMPI_Op_fn)(void*, void*, void*, int);


/* the request object for asynchronious operations. */
struct mpi_req_
{
    bool recv; /* whether it's a receive request. */
    bool finished; /* whether it's finished */

    struct mpi_recv_envelope_ evr; /* receive envelope */
    struct mpi_send_envelope_ evs; /* send envelope */
};


struct mpi_datatype_component_
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
    struct mpi_datatype_component_ *comps; /* the components */
    bool committed; /* whether the data type is committed */
};
/* just as a shorthand:  */
typedef struct mpi_datatype_ mpi_dt_;


/* there are a few global variables that maintain information about the
   running threads. Some are defined by the MPI standard: */
MPI_Comm MPI_COMM_WORLD=NULL;
MPI_Group MPI_GROUP_EMPTY=NULL;

/* the threads themselves (mpi_comm only contains lists of pointers to this
   structure */
struct mpi_thread_ *threads=NULL;
int Nthreads=0;

/* thread info */
gmx_thread_key_t id_key; /* the key to get the thread id */

/* whether MPI has finalized (we need this to distinguish pre-inited from
    post-finalized states */
static bool mpi_finalized=FALSE;

/* misc. global information about MPI */
struct mpi_global_ *mpi_global=NULL;




/* unique tags for multicast */
#define TMPI_BCAST_TAG      1
#define TMPI_GATHER_TAG     2
#define TMPI_SCATTER_TAG    3
#define TMPI_REDUCE_TAG     4
#define TMPI_ALLTOALL_TAG   5
#define TMPI_ALLTOALLV_TAG  6



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


/* error messages */
char *mpi_errmsg[] = 
{
    "No error",
    "Invalid MPI_Group",
    "Invalid MPI_Comm",
    "Invalid MPI_Status",
    "Invalid Cartesian topology dimensions",
    "Invalid Cartesian topology coordinates",
    "Insufficient number processes for Cartesian topology",
    "Invalid counterpart for MPI transfer",
    "Receive buffer size too small for transmission",
    "Invalid send destination",
    "Invalid receive source",
    "Invalid buffer",
    "Invalid reduce operator",
    "Transmission failure"
    "Unknown MPI error"
};



/* get the current thread structure pointer */
static struct mpi_thread_ *tMPI_Get_current(void);
/* mostly for debugging: gives current thread if thr==0*/
static int tMPI_Threadnr(struct mpi_thread_ *thr);


/* check whether we're the main thread */
static bool tMPI_Is_master(void);
/* check whether the current process is in a group */
static bool tMPI_In_group(MPI_Group group);

/* find the rank of a thread in a comm */
static int tMPI_Comm_seek_rank(MPI_Comm comm, struct mpi_thread_ *th);
/* find the size of a comm */
static int tMPI_Comm_N(MPI_Comm comm);

/* start N threads with argc, argv (used by MPI_Init) */
static void tMPI_Start_threads(int N, int *argc, char ***argv);
/* starter function for threads; calls main() */
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
static struct mpi_thread_ *tMPI_Get_thread(MPI_Comm comm, int rank);

/* fill the pre-allocated envelope sd with the right values */
static void tMPI_Prep_send_envelope(struct mpi_send_envelope_ *ev, 
                                    MPI_Comm comm, struct mpi_thread_ *dest,
                                    void *send_buf, int send_count, 
                                    MPI_Datatype datatype, int tag,
                                    enum xfertype_ xt);

/* fill the pre-allocated envelope sd with the right values */
static void tMPI_Prep_recv_envelope(struct mpi_recv_envelope_ *ev, 
                                    MPI_Comm comm, struct mpi_thread_ *src, 
                                    void *recv_buf, int recv_count,
                                    MPI_Datatype datatype, int tag,
                                    enum xfertype_ xt);



/* check whether a send and receive envelope match (same dest, tag, etc.) */
static bool tMPI_Envelope_matches(const struct mpi_send_envelope_ *send, 
                                  const struct mpi_recv_envelope_ *recv);

/* try to post an envelope as a sender */
static void tMPI_Try_put_send_envelope(struct mpi_send_envelope_ *ev);

/* try to post an envelope as a receiver (indicating readiness to receive)*/
static void tMPI_Try_put_recv_envelope(struct mpi_recv_envelope_ *ev);



/* wait for a specific envelope's data to be come ready */
static void tMPI_Wait_recv(struct mpi_recv_envelope_ *ev);
/* wait for a specific envelope transfer finish */
static void tMPI_Wait_send(struct mpi_send_envelope_ *ev);


/* do a transfer from an envelope send to recv buffer, updating status. 
   This function will lock a mutex if the envelope says that there's a
   non-blocking receiver (in that case, the sender may have to initiate
   transer). If try=true, it will just try to lock (and return if it doesn't
   get a lock), if try=false, it will wait until it gets a lock. 
   
   There are two versions of Xfer: one for recv-initiated transfers (Xfer_recv)
   and one for sender-initiated transfers (Xfer-send). */
static int tMPI_Xfer_recv(struct mpi_recv_envelope_ *evr, MPI_Status *status, 
                          bool try);
static int tMPI_Xfer_send(struct mpi_send_envelope_ *evs, MPI_Status *status, 
                          bool try);



/* these are the internal versions of MPI_Isend/Irecv, etc that use 
   pre-allocated mpi_reqs so as not to incur the overhead of allocations. 
   They are used by the regular MPI_Isend/Irecv, etc.  */
static int MPI_Recv_r(void* buf, int count, MPI_Datatype datatype, 
                      int source, int tag, MPI_Comm comm, 
                      MPI_Status *status, enum xfertype_ xt);
static int tMPI_Isend_r(void* buf, int count, MPI_Datatype datatype, 
                        int dest, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, enum xfertype_ xt);
static int tMPI_Irecv_r(void* buf, int count, MPI_Datatype datatype, 
                        int source, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, enum xfertype_ xt);
static int tMPI_Test_r(struct mpi_req_ *rq, int *flag, MPI_Status *status);
static int tMPI_Wait_r(struct mpi_req_ *rq, MPI_Status *status);
/* wait all on an array of pointers to requests, that progressively get 
   assigned NULL once they're done. If may_free, the pointers also
   get free()d */
static int tMPI_Waitall_r(int count, struct mpi_req_ *array_of_requests[],
                          MPI_Status *array_of_statuses, bool may_free);

/* run a single binary reduce operation on src_a and src_b, producing dest. 
   dest and src_a may be identical */
static void tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
                               MPI_Datatype datatype, int count, MPI_Op op);


/* and we need this prototype */
int main(int argc, char **argv);





static struct mpi_thread_ *tMPI_Get_current(void)
{
    if (!threads)
        return NULL;

    return (struct mpi_thread_*)gmx_thread_getspecific(id_key);
} 
static int tMPI_Threadnr(struct mpi_thread_ *thr)
{
    if (thr)
        return tMPI_Comm_seek_rank(MPI_COMM_WORLD, thr);

    return tMPI_Comm_seek_rank(MPI_COMM_WORLD, tMPI_Get_current());
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
    struct mpi_thread_* th=tMPI_Get_current();
    return th->self_comm;
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
    fprintf(stderr, "MPI error: %s\n", errstr);
    abort();
    /*exit(0);*/
}

static void mpi_errors_return_fn(MPI_Comm *comm, int *err)
{
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
    struct mpi_thread_ *th=(struct mpi_thread_*)arg;
    gmx_thread_setspecific(id_key, arg);

    main(th->argc, th->argv);
    return 0;
}


void tMPI_Start_threads(int N, int *argc, char ***argv)
{
    if (N>1) 
    {
        int i;

        mpi_finalized=FALSE;
        Nthreads=N;

        /* allocate global data */
        smalloc(mpi_global, sizeof(struct mpi_global_));
        mpi_global->usertypes=NULL;
        mpi_global->N_usertypes=0;
        mpi_global->Nalloc_usertypes=0;

        /* allocate world and thread data */
        smalloc(threads, sizeof(struct mpi_thread_)*N);
        MPI_COMM_WORLD=tMPI_Comm_alloc(NULL, N);
        MPI_GROUP_EMPTY=tMPI_Group_alloc();

        /* allocate mutex for communicator creation */
        gmx_thread_mutex_init( &(MPI_COMM_WORLD->comm_create_mutex) );
        gmx_thread_cond_init( &(MPI_COMM_WORLD->comm_create_prep) );
        MPI_COMM_WORLD->split = NULL;
        MPI_COMM_WORLD->new_comm = NULL;


        MPI_COMM_WORLD->grp.N=N;

        if (gmx_thread_key_create(&id_key, NULL))
        {
            gmx_fatal(errno, __FILE__, __LINE__, 
                    "thread_key_create failed");
        }
        /*printf("thread keys created\n"); fflush(NULL);*/
        for(i=0;i<N;i++)
        {
            threads[i].N_evr=0;
            threads[i].N_evs=0;
            MPI_COMM_WORLD->grp.peers[i]=&(threads[i]);

            /* create the mutexes and conditions */
            gmx_thread_mutex_init(&(threads[i].envelope_mutex));
            gmx_thread_mutex_init(&(threads[i].xmit_mutex));

            /* allocate comm.self */
            threads[i].self_comm=tMPI_Comm_alloc(MPI_COMM_WORLD, 1);
            threads[i].self_comm->grp.peers[0]=&(threads[i]);

            /* copy argc, argv */
            if (argc && argv)
            {
                int j;
                threads[i].argc=*argc;
                smalloc(threads[i].argv, threads[i].argc*sizeof(char*));
                for(j=0;j<threads[i].argc;j++)
                {
                    size_t len=strlen((*argv)[j]);
                    smalloc(threads[i].argv[j], len+1);
                    strncpy(threads[i].argv[j], (*argv)[j], len);
                }
            }
            else
            {
                threads[i].argc=0;
                threads[i].argv=0;
            }
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
                gmx_fatal(errno, __FILE__, __LINE__, "thread_create failed");
            }
        }
    }
}


int MPI_Init(int *argc, char ***argv)
{
    if (MPI_COMM_WORLD==0) /* we're the main process */
    {
        int N=tMPI_Get_N(argc, argv);
        tMPI_Start_threads(N, argc, argv);
    }
    else
    {
        /* if we're a sub-thread we need don't need to do anyhing, because 
           everything has already been set up by either the main thread, 
           or the thread runner function.*/
    }
    return MPI_SUCCESS;
}

int MPI_Init_N(int N)
{
    if (MPI_COMM_WORLD==0 && N>1) /* we're the main process */
    {
        tMPI_Start_threads(N, 0, 0);
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

#ifdef TMPI_DEBUG
    printf("%d: MPI_Finalize called\n", tMPI_Threadnr(0));
    fflush(0);
#endif
    if (tMPI_Is_master())
    {
#ifdef TMPI_DEBUG
        for(i=0;i<Nthreads;i++)
        {
            if (threads[i].N_evs)
            {
                printf(
                    "MPI_Finalize: %d open send envelopes in thread %d\n",
                    threads[i].N_evs, i);
            }
            if (threads[i].N_evr)
            {
                printf(
                    "MPI_Finalize: %d open recv envelopes in thread %d\n",
                    threads[i].N_evr, i);
            }
        }
#endif
        /* we just wait for all threads to finish; the order isn't very 
           relevant, as all threads should arrive at their endpoints soon. */
        for(i=1;i<Nthreads;i++)
        {
            if (gmx_thread_join(threads[i].thread_id, NULL))
                gmx_fatal(errno, __FILE__, __LINE__, "thread_join failed");
        }
        for(i=0;i<Nthreads;i++)
        {
            gmx_thread_mutex_destroy(&(threads[i].envelope_mutex));
            gmx_thread_mutex_destroy(&(threads[i].xmit_mutex));
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
    fprintf(stderr, "MPI_Abort called");
    if (comm==MPI_COMM_WORLD)
        fprintf(stderr, " on MPI_COMM_WORLD");
    fprintf(stderr,"\n");
    fflush(0);

    abort();
#if 0
    /* we just kill all threads, but not the main process */
    int i;
    struct mpi_thread_ *me=tMPI_Get_current();
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
    /*int *ret;

    ret=gmx_thread_getspecific(id_key);*/

    *resultlen=snprintf(name, MPI_MAX_PROCESSOR_NAME, "thread #%d", 
                        tMPI_Threadnr(0));
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
        MPI_COMM_WORLD->erh->err=MPI_ERR_STATUS;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD,
                                       &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
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
    smalloc(ntp->comps, sizeof(struct mpi_datatype_component_)*1);
    ntp->comps[0].type=oldtype;
    ntp->comps[0].count=1;
    ntp->committed=FALSE;

    /* now add it to the list.  */
    gmx_thread_mutex_lock(&(mpi_global->datatype_mutex));
    /* check whether there's space */
    if (mpi_global->N_usertypes + 1 >= mpi_global->Nalloc_usertypes)
    {
        /* make space */
        mpi_global->Nalloc_usertypes=2*(mpi_global->N_usertypes + 1);
        mpi_global->usertypes=srealloc(mpi_global->usertypes, 
                                        (sizeof(struct mpi_datatype_ *)*
                                          mpi_global->Nalloc_usertypes));
    }
    /* add to the list */
    mpi_global->usertypes[mpi_global->N_usertypes]=ntp;
    mpi_global->N_usertypes++;
    gmx_thread_mutex_unlock(&(mpi_global->datatype_mutex));

    *newtype=ntp;

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
    gmx_thread_mutex_lock(&(mpi_global->datatype_mutex));
    for(i=0;i<mpi_global->N_usertypes;i++)
    {
        struct mpi_datatype_ *lt=mpi_global->usertypes[i];
        if (lt->committed && lt->N_comp==dt->N_comp)
        {
            bool found=TRUE;
            for(j=0;j<lt->N_comp;j++)
            {
                if ( (lt->comps[i].type  != dt->comps[i].type) ||
                     (lt->comps[i].count != dt->comps[i].count) )
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
    gmx_thread_mutex_unlock(&(mpi_global->datatype_mutex));
    return MPI_SUCCESS;
}



/* Group query & manipulation functions */

static bool tMPI_In_group(MPI_Group group)
{
    int i;
    struct mpi_thread_ *cur;

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
    struct mpi_thread_ *cur;

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
    smalloc(ret->peers, sizeof(struct mpi_thread_*)*Nthreads);
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
    if (comm)
    {
        *group=&(comm->grp);
    }
    else
    {
        *group=NULL;
    }

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
            gmx_fatal(FARGS, "MPI_Group_incl: group rank illegal");
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
    smalloc(ret->grp.peers, sizeof(struct mpi_thread_*)*Nthreads);
    ret->grp.N=N;

    gmx_thread_mutex_init( &(ret->multicast_mutex) );
    gmx_thread_cond_init( &(ret->multicast_cond) );

    /* calculate the number of multicast barriers */
    Nbarriers=0;
    Nred=N;
    while(Nred>1) {
        Nbarriers+=1;
        Nred = Nred/2 + Nred%2;
    } 

    ret->Nbarriers=Nbarriers;
    smalloc(ret->multicast_barrier, 
            sizeof(gmx_thread_barrier_t)*(Nbarriers+1));
    smalloc(ret->N_multicast_barrier, sizeof(int)*(Nbarriers+1));
    Nred=N;
    for(i=0;i<Nbarriers;i++)
    {
        gmx_thread_barrier_init( &(ret->multicast_barrier[i]), Nred);
        ret->N_multicast_barrier[i]=Nred;
        /* Nred is now Nred/2 + a rest term because solitary 
           process at the end of the list must still be accounter for */
        Nred = Nred/2 + Nred%2;
    }
    smalloc(ret->sendbuf, sizeof(void*)*Nthreads);
    smalloc(ret->recvbuf, sizeof(void*)*Nthreads);

    gmx_thread_mutex_init( &(ret->comm_create_mutex) );
    gmx_thread_cond_init( &(ret->comm_create_prep) );
    gmx_thread_cond_init( &(ret->comm_create_finish) );

    if (parent)
    {
        ret->erh=parent->erh;
    }
    else
    {
        ret->erh=MPI_ERRORS_ARE_FATAL;
    }


    /* we have no topology to start out with */
    ret->cart=NULL;
    /*ret->graph=NULL;*/

#if 0
    ret->grp.Nref=1;
#endif
    /*ret->context_tag=mpi_context_tags++;*/
    return ret;
}

int MPI_Comm_free(MPI_Comm *comm)
{
    if (! *comm)
        return MPI_SUCCESS;

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
    volatile int colors[MAX_THREADS]; /* array with the colors of each thread */
    volatile int keys[MAX_THREADS]; /* same for keys (only one of the 
                                       threads actually suplies these arrays 
                                       to the comm structure) */
    bool i_am_first=FALSE;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());
    struct mpi_split_ *spl;

    if (!comm)
    {
        *newcomm=NULL;
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    /*printf("**Calling MPI_Comm_split with color %d, key %d\n",color, key);*/

    gmx_thread_mutex_lock(&(comm->comm_create_mutex));    
    /* first get the colors */
    if (!comm->new_comm)
    {
        smalloc(comm->split, sizeof(struct mpi_split_));
        smalloc(comm->new_comm, N*sizeof(MPI_Comm));
        if (N<=MAX_THREADS)
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
                              &(comm->comm_create_mutex) );
    }
    else
    {
        int Ncomms=0;
        int comm_color_[MAX_THREADS]; 
        int comm_N_[MAX_THREADS]; 
        int *comm_color=comm_color_; /* there can't be more comms than N*/
        int *comm_N=comm_N_; /* the number of procs in a group */

        int *comm_groups; /* the groups */
        MPI_Comm *comms; /* the communicators */

        /* wait for the colors to be done */
        gmx_thread_cond_wait( &(comm->comm_create_prep), 
                              &(comm->comm_create_mutex));

        /* reset the state so that a new comm creating function can run */
        spl->Ncol_destroy=N;
        comm->new_comm=0;
        comm->split=0;

        smalloc(comm_groups, N*N*sizeof(int));
        if (N>MAX_THREADS)
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
        if (N>MAX_THREADS)
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

    gmx_thread_mutex_unlock(&(comm->comm_create_mutex));    

    return MPI_SUCCESS;    
}

static int tMPI_Comm_seek_rank(MPI_Comm comm, struct mpi_thread_ *th)
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
int MPI_Topo_test(MPI_Comm comm, int status)
{
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }

    if (comm->cart)
        status=MPI_CART;
    /*else if (comm->graph)
        status=MPI_GRAPH;*/
    else 
        status=MPI_UNDEFINED;

    return MPI_SUCCESS;
}

int MPI_Cartdim_get(MPI_Comm comm, int *ndims)
{
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
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
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!comm->cart || comm->cart->ndims==0)
        return MPI_SUCCESS;

    MPI_Cart_coords(comm, myrank, maxdims, coords);

    for(i=0;i<comm->cart->ndims;i++)
    {
        if (i>=maxdims)
        {
            comm->erh->err=MPI_ERR_DIMS;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
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
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
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
                comm->erh->err=MPI_ERR_DIMS;
                comm->erh->fn(&comm, &(comm->erh->err));
                return comm->erh->err;
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
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!comm->cart || comm->cart->ndims==0)
        return MPI_SUCCESS;
    if (maxdims < comm->cart->ndims)
    {
        comm->erh->err=MPI_ERR_DIMS;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    /* again, row-major ordering */
    for(i=comm->cart->ndims-1;i>=0;i--)
    {
        /*printf("dims[%d]=%d\n",i,dims[i]);*/
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
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
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
        comm_old->erh->err=MPI_ERR_COMM;
        comm_old->erh->fn(&comm_old, &(comm_old->erh->err));
        return comm_old->erh->err;
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
        comm_old->erh->err=MPI_ERR_CART_CREATE_NPROCS;
        comm_old->erh->fn(&comm_old, &(comm_old->erh->err));
        return comm_old->erh->err;
#endif
        return MPI_FAILURE;
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
        smalloc((*comm_cart)->cart, sizeof(struct mpi_cart_topol_));
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
        gmx_thread_barrier_wait( &( (*comm_cart)->multicast_barrier[0]) );
    }


    return MPI_SUCCESS;
}


/* Point-to-point communication protocol functions */

static struct mpi_thread_ *tMPI_Get_thread(MPI_Comm comm, int rank)
{
    /* check destination */
    if ( (rank < 0) || (rank > comm->grp.N) )
    {
        gmx_fatal(FARGS, "illegal rank");
        return 0;
    }
    return comm->grp.peers[rank];
}


static void tMPI_Prep_send_envelope(struct mpi_send_envelope_ *ev, 
                                    MPI_Comm comm, struct mpi_thread_ *dest,
                                    void *send_buf, int send_count,
                                    MPI_Datatype datatype, int tag, 
                                    enum xfertype_ xt)
{
    ev->comm=comm;
    ev->src=tMPI_Get_current();;
    ev->dest=dest;

    ev->tag=tag;
    ev->xt=xt;

    ev->buf=send_buf;
    ev->bufsize=send_count*datatype->size;
    ev->datatype=datatype;

    ev->ready=FALSE; 
    ev->nonblock=FALSE;
    ev->finished=FALSE;
    ev->recv_ready=FALSE;
    ev->recv_nonblock=FALSE;

    ev->counterpart=NULL;
}


static void tMPI_Prep_recv_envelope(struct mpi_recv_envelope_ *ev, 
                                    MPI_Comm comm, struct mpi_thread_ *src, 
                                    void *recv_buf, int recv_count,
                                    MPI_Datatype datatype, int tag, 
                                    enum xfertype_ xt)
{
    ev->comm=comm;
    ev->src=src;
    ev->dest=tMPI_Get_current();

    ev->tag=tag;
    ev->xt=xt;

    ev->buf=recv_buf;
    ev->bufsize=recv_count*datatype->size; 
    ev->datatype=datatype;

    ev->ready=FALSE; 
    ev->nonblock=FALSE;
    ev->finished=FALSE;
    ev->send_ready=FALSE;
    ev->send_nonblock=FALSE;

    ev->counterpart=NULL;
}



static bool tMPI_Envelope_matches(const struct mpi_send_envelope_ *send, 
                                  const struct mpi_recv_envelope_ *recv)
{
#ifdef TMPI_DEBUG
    printf("%d: (%d->%d)==(%d->%d),  tag=(%d==%d), xfertype=(%d==%d), datatype=(%d==%d), finished=(%d==%d)\n",
            tMPI_Threadnr(0),
            tMPI_Threadnr(send->src), tMPI_Threadnr(send->dest),
            tMPI_Threadnr(recv->src), tMPI_Threadnr(recv->dest),
            (int)(send->tag), (int)(recv->tag), 
            (int)(send->xt), (int)(recv->xt), 
            (int)recv->datatype, (int)send->datatype,
            recv->finished, send->finished);
#endif

    if ( (send->comm == recv->comm) &&
         (send->dest == recv->dest) &&
         ( recv->xt == send->xt ) &&
         ( (!recv->src)  || (recv->src == send->src) ) &&
         ( (recv->tag == MPI_ANY_TAG) || (recv->tag == send->tag) ) &&
         ( recv->datatype == send->datatype ) &&
         ( ! (recv->finished || send->finished) ) ) 
    {
#ifdef TMPI_DEBUG
        printf("%d: (%d->%d) tag=%d found match\n",
            tMPI_Threadnr(0),
            tMPI_Threadnr(send->src), tMPI_Threadnr(send->dest),
            (int)(send->tag));
#endif
        return TRUE;
    }
    return FALSE;
}


static void tMPI_Try_put_send_envelope(struct mpi_send_envelope_ *ev)
{
    int i;
    bool add=TRUE;
#ifdef TMPI_DEBUG
    printf("%d: Entering try_put_send_envelope: (%d->%d), tag=%d\n",
            tMPI_Threadnr(0), tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest),
            ev->tag); 
    fflush(0);
#endif

    gmx_thread_mutex_lock(&(ev->dest->envelope_mutex));
    /* search for matching receive envelopes */
    for(i=0;i<ev->dest->N_evr;i++)
    {
        struct mpi_recv_envelope_ *evr=ev->dest->evr[i];
        if (tMPI_Envelope_matches(ev, evr))
        {
#ifdef TMPI_DEBUG
            printf("%d: Send: found envelope\n", tMPI_Threadnr(0)); fflush(0);
#endif
            add=FALSE;
            /* establish reciprocity */
            ev->counterpart=evr;
            evr->counterpart=ev;

            /* the receiving end doesn't neccesarily know these */
            evr->src=ev->src; 
            evr->tag=ev->tag; 

            /* get flags from the other end*/
            ev->recv_ready=evr->ready;
            ev->recv_nonblock=evr->nonblock;
            /* set flags on the other end */
            evr->send_nonblock=ev->nonblock;
            evr->send_ready=ev->ready;

            /* now remove the envelope from the list, it's not intended
               to be found anymore. 
               We just put the last one in our slot:*/
            ev->dest->evr[i] = ev->dest->evr[ev->dest->N_evr-1];
            ev->dest->N_evr--;
            break;
        }
    }

    if (add)
    {
#ifdef TMPI_DEBUG
        printf("%d: Send: adding envelope\n", tMPI_Threadnr(0)); fflush(0);
#endif
        /* we add the envelope to the receiver's list */
        if (ev->dest->N_evs >= MAX_ENVELOPES)
            gmx_fatal(FARGS, "A receive thread has reached the maximum number of send envelopes. This is probably a bug\n");
        ev->dest->evs[ev->dest->N_evs]=ev;
        ev->dest->N_evs++;
    }
    ev->ready=TRUE;
    if (ev->counterpart)
    {
#ifdef TMPI_DEBUG
        printf("%d: Send: signaling send_ready to %d\n", tMPI_Threadnr(0),
                    tMPI_Threadnr(ev->dest));
#endif
        ev->counterpart->send_ready=TRUE;
    }
    /* now we let other threads do their thing */
    gmx_thread_mutex_unlock(&(ev->dest->envelope_mutex));
}

static void tMPI_Try_put_recv_envelope(struct mpi_recv_envelope_ *ev)
{
    int i;
    bool add=TRUE;
#ifdef TMPI_DEBUG
    printf("%d: Entering try_put_recv_envelope: %d->%d, tag=%d\n",
            tMPI_Threadnr(0),tMPI_Threadnr(ev->src),tMPI_Threadnr(ev->dest),
            ev->tag); 
    fflush(0);
#endif

    gmx_thread_mutex_lock(&(ev->dest->envelope_mutex));
    /* search for matching envelopes */
    for(i=0;i<ev->dest->N_evs;i++)
    {
        struct mpi_send_envelope_ *evs=ev->dest->evs[i];
        if (tMPI_Envelope_matches(evs, ev))
        {
#ifdef TMPI_DEBUG
            printf("%d: Recv: found envelope\n", tMPI_Threadnr(0)); fflush(0);
#endif            
            add=FALSE;
            /* establish reciprocity */
            ev->counterpart=evs;
            evs->counterpart=ev;

            /* the receiving end doesn't neccesarily know these */
            ev->src=evs->src; 
            ev->tag=evs->tag; 

            /* get flags from the other end */
            ev->send_ready=evs->ready;
            ev->send_nonblock=evs->nonblock;
            /* set flags on the other end */
            evs->recv_nonblock=ev->nonblock;
            evs->recv_ready=ev->ready;

            /* now remove the envelope from the list, it's not intended
               to be found anymore. 
               We just put the last one in our slot:*/
            ev->dest->evs[i] = ev->dest->evs[ev->dest->N_evs-1];
            ev->dest->N_evs--;

            break;
        }
    }

    if (add)
    {
        /* we add the envelope to the receiver's list */
#ifdef TMPI_DEBUG
        printf("%d: Recv: adding envelope\n", tMPI_Threadnr(0)); fflush(0);
#endif
        if (ev->dest->N_evr >= MAX_ENVELOPES)
        {
            gmx_fatal(FARGS, "A receive thread has reached the maximum number of recv envelopes. This is probably a bug\n");
        }
        ev->dest->evr[ev->dest->N_evr]=ev;
        ev->dest->N_evr++;
    }
    ev->ready=TRUE;
    if (ev->counterpart)
    {
#ifdef TMPI_DEBUG
        printf("%d: Send: signaling recv_ready to %d\n", tMPI_Threadnr(0),
                tMPI_Threadnr(ev->src));
#endif
        ev->counterpart->recv_ready=TRUE;
    }
    /* we let other threads do their thing */
    gmx_thread_mutex_unlock(&(ev->dest->envelope_mutex));
}

static void tMPI_Wait_recv(struct mpi_recv_envelope_ *ev)
{
#ifdef TMPI_DEBUG
    printf("%d: Recv: waiting (%d->%d, tag=%d)\n", tMPI_Threadnr(0),
           tMPI_Threadnr(ev->src), tMPI_Threadnr(ev->dest), ev->tag); 
    fflush(0);
#endif
    /* busy-wait until the data is ready. Busy waiting is not
       so bad in this case because we need low latency, and 
       all threads scheduled. */
    while( ! (ev->send_ready || ev->finished ) )
    {
    }
}

static void tMPI_Wait_send(struct mpi_send_envelope_ *ev)
{
#ifdef TMPI_DEBUG
    printf("%d: Send: waiting (%d->%d, tag=%d)\n", tMPI_Threadnr(0),
                                                   tMPI_Threadnr(ev->src), 
                                                   tMPI_Threadnr(ev->dest), 
                                                   ev->tag);
    fflush(0);
#endif
    /* we busy wait until the data is either picked up, or we might
       need ot send it ourselves.  */
    while( ! (ev->finished || ( ev->recv_nonblock && ev->recv_ready) ) )
    {
    }
}

static int tMPI_Xfer_recv(struct mpi_recv_envelope_ *evr, MPI_Status *status, 
                         bool try)
{
    bool locked=FALSE;
    size_t transferred;
    int ret=MPI_SUCCESS;
    int rets=MPI_SUCCESS;
#ifdef TMPI_DEBUG
    printf("%d: Xfering data from %d -> %d (tag=%d), recv\n", tMPI_Threadnr(0), 
            tMPI_Threadnr(evr->src), tMPI_Threadnr(evr->dest), evr->tag );
    fflush(0);
#endif
    if (evr->nonblock)
    {
        if (try)
        {
            int ret=gmx_thread_mutex_trylock(&(evr->dest->xmit_mutex));
            if (ret==EBUSY)
                return MPI_SUCCESS;
        }
        else
        {
            gmx_thread_mutex_lock(&(evr->dest->xmit_mutex));
        }
        locked=TRUE;
    }

    /* we wait with this check until we either have a lock, or we know
       we don't need one */
    if (evr->finished ) 
    {
        goto end;
    }

    if (!evr->counterpart)
    {
        evr->comm->erh->err=MPI_ERR_XFER_COUNTERPART;
        evr->comm->erh->fn(&(evr->comm), &(evr->comm->erh->err));
        rets=MPI_FAILURE;
        ret=evr->comm->erh->err;
        goto end;
    }

    transferred=evr->counterpart->bufsize;
    if (evr->bufsize < transferred)
    {
        evr->comm->erh->err=MPI_ERR_XFER_BUFSIZE;
        evr->comm->erh->fn(&evr->comm, &(evr->comm->erh->err));
        rets=MPI_FAILURE;
        ret=evr->comm->erh->err;
        evr->finished=TRUE;
        evr->counterpart->finished=TRUE;
        goto end;
    }
    /* we do the actual transfer */
    memcpy(evr->buf, evr->counterpart->buf, transferred);
    evr->transferred=transferred;
    if(status)
    {
        status->MPI_SOURCE=tMPI_Comm_seek_rank(evr->comm, evr->src);
        status->MPI_TAG=evr->tag;
        status->MPI_ERROR=MPI_SUCCESS;
        status->transferred=transferred;
    }

    /* and we set the finished flag. This makes the blocking MPI_Send/MPI_Wait,
       etc. exit., and the counterpart structure invalid.*/
    evr->finished=TRUE;
    evr->counterpart->finished=TRUE;
end:
    if (status)
        status->MPI_ERROR=rets;
    if (locked)
        gmx_thread_mutex_unlock(&(evr->dest->xmit_mutex));
    return ret;
}

static int tMPI_Xfer_send(struct mpi_send_envelope_ *evs, 
                          MPI_Status *status, bool try)
{
    bool locked=FALSE;
    size_t transferred;
    int ret=MPI_SUCCESS;
    int rets=MPI_SUCCESS;
#ifdef TMPI_DEBUG
    printf("%d: Xfering data from %d -> %d (tag=%d), send (wrong direction)\n", 
            tMPI_Threadnr(0), 
            tMPI_Threadnr(evs->src), tMPI_Threadnr(evs->dest), evs->tag );
    fflush(0);
#endif
    if (evs->recv_nonblock)
    {
        if (try)
        {
            int ret=gmx_thread_mutex_trylock(&(evs->dest->xmit_mutex));
            if (ret==EBUSY)
                return MPI_SUCCESS;
        }
        else
        {
            gmx_thread_mutex_lock(&(evs->dest->xmit_mutex));
        }
        locked=TRUE;
    }
    /* we wait with this check until we either have a lock, or we know
       we don't need one */
    if (evs->finished) 
    {
        ret=MPI_SUCCESS;
        goto end;
    }
    if (!evs->counterpart)
    {
        evs->comm->erh->err=MPI_ERR_XFER_COUNTERPART;
        evs->comm->erh->fn(&(evs->comm), &(evs->comm->erh->err));
        rets=MPI_FAILURE;
        ret=evs->comm->erh->err;
        goto end;
    }
    transferred=evs->bufsize;
    if (evs->counterpart->bufsize < transferred)
    {
        evs->comm->erh->err=MPI_ERR_XFER_BUFSIZE;
        evs->comm->erh->fn(&evs->comm, &(evs->comm->erh->err));
        rets=MPI_FAILURE;
        ret=evs->comm->erh->err;
        evs->finished=TRUE;
        evs->counterpart->finished=TRUE;
        goto end;
    }
    /* we do the actual transfer */
    memcpy(evs->counterpart->buf, evs->buf, transferred);
    evs->counterpart->transferred=transferred;
    if(status)
    {
        status->MPI_SOURCE=tMPI_Comm_seek_rank(evs->comm, evs->src);
        status->MPI_TAG=evs->tag;
        status->MPI_ERROR=MPI_SUCCESS;
        status->transferred=transferred;
    }

    /* and we set the finished flag. This makes the blocking MPI_Send/MPI_Wait,
       etc. exit., and the counterpart structure invalid.*/
    evs->finished=TRUE;
    evs->counterpart->finished=TRUE;
end:
    if (status)
        status->MPI_ERROR=ret;
    if (locked)
        gmx_thread_mutex_unlock(&(evs->dest->xmit_mutex));
    return ret;
}

/* point-to-point communication */


static int MPI_Send_r(void* buf, int count, MPI_Datatype datatype, int dest,
                      int tag, MPI_Comm comm, enum xfertype_ xt)
{
    struct mpi_send_envelope_ sd;
    struct mpi_thread_ *send_dst=tMPI_Get_thread(comm, dest);
    int ret=MPI_SUCCESS;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!send_dst)
    {
        comm->erh->err=MPI_ERR_SEND_DEST;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
    if (!buf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    tMPI_Prep_send_envelope(&sd, comm, send_dst, buf, count, datatype, tag, xt);
    tMPI_Try_put_send_envelope(&sd);
    tMPI_Wait_send(&sd);

    /* this xfer will only happen if there is a non-blocking receive on the
       other end. */
    if (!sd.finished)
    {
        ret=tMPI_Xfer_send(&sd, 0, FALSE);
    }
    return ret;    
}

static int MPI_Recv_r(void* buf, int count, MPI_Datatype datatype, 
                      int source, int tag, MPI_Comm comm, 
                      MPI_Status *status, enum xfertype_ xt)
{
    struct mpi_recv_envelope_ rc;
    struct mpi_thread_ *recv_src=0;
    int ret;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!buf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    if(status)
    {
        status->MPI_ERROR=MPI_FAILURE; /* assume we fail */
    }

    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            comm->erh->err=MPI_ERR_RECV_SRC;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
    }

    tMPI_Prep_recv_envelope(&rc, comm, recv_src, buf, count, datatype, tag, xt);
    tMPI_Try_put_recv_envelope(&rc);
    tMPI_Wait_recv(&rc);

    /* now we copy the data to our side */
    ret=tMPI_Xfer_recv(&rc, status, FALSE);
    return ret;
}


int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                 int dest, int sendtag, void *recvbuf, int recvcount,
                 MPI_Datatype recvtype, int source, int recvtag, 
                 MPI_Comm comm, MPI_Status *status)
{
    struct mpi_recv_envelope_ rc;
    struct mpi_send_envelope_ sd;

    struct mpi_thread_ *recv_src=0;
    struct mpi_thread_ *send_dst=tMPI_Get_thread(comm, dest);
    int ret=MPI_SUCCESS;
    int ret2=MPI_SUCCESS;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!send_dst)
    {
        comm->erh->err=MPI_ERR_SEND_DEST;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
    if (!sendbuf || !recvbuf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
    if(status)
    {
        status->MPI_ERROR=MPI_FAILURE; /* assume we fail */
    }
    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            comm->erh->err=MPI_ERR_RECV_SRC;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
    }
    /* we first prepare to send */
    tMPI_Prep_send_envelope(&sd, comm, send_dst, sendbuf, sendcount, 
                            sendtype, sendtag, point_to_point);
    tMPI_Try_put_send_envelope(&sd);

    /* then we prepare to receive */
    tMPI_Prep_recv_envelope(&rc, comm, recv_src, recvbuf, recvcount, 
                            recvtype, recvtag, point_to_point);
    tMPI_Try_put_recv_envelope(&rc);

    /* now we wait to receive */
    tMPI_Wait_recv(&rc);
    /* we actually receive */
    ret=tMPI_Xfer_recv(&rc, status, FALSE);

    /* we wait until the send completes */
    tMPI_Wait_send(&sd);
    /* this xfer will only happen if there is a non-blocking receive on the
       other end. */
    if (!sd.finished)
    {
        ret2=tMPI_Xfer_send(&sd, 0, FALSE);
    }

    if (ret2!=MPI_SUCCESS)
        ret=ret2;

    return ret;
}


/* async */
static int tMPI_Isend_r(void* buf, int count, MPI_Datatype datatype, int dest,
                        int tag, MPI_Comm comm, struct mpi_req_ *rq, 
                        enum xfertype_ xt)
{
    struct mpi_thread_ *send_dst=tMPI_Get_thread(comm, dest);

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!send_dst)
    {
        comm->erh->err=MPI_ERR_SEND_DEST;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
    if (!buf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    rq->recv=FALSE;
    rq->finished=FALSE;
    tMPI_Prep_send_envelope(&(rq->evs), comm, send_dst, buf, count, 
                            datatype, tag, xt);
    rq->evs.nonblock=TRUE;
    tMPI_Try_put_send_envelope(&(rq->evs));
    /* now we don't wait and assume that data isn't touched until MPI_WAIT */
    return MPI_SUCCESS;    
}

static int tMPI_Irecv_r(void* buf, int count, MPI_Datatype datatype, 
                        int source, int tag, MPI_Comm comm, 
                        struct mpi_req_ *rq, enum xfertype_ xt)
{
    struct mpi_thread_ *recv_src=0;
    int ret=MPI_SUCCESS;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!buf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    if (source!=MPI_ANY_SOURCE)
    {
        recv_src=tMPI_Get_thread(comm, source);
        if (!recv_src)
        {
            comm->erh->err=MPI_ERR_RECV_SRC;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
    }
    rq->recv=TRUE;
    rq->finished=FALSE;

    tMPI_Prep_recv_envelope(&(rq->evr), comm, recv_src, buf, count, 
                            datatype, tag, xt);
    rq->evr.nonblock=TRUE;
    tMPI_Try_put_recv_envelope(&(rq->evr));
    /* if there's already data waiting, we just xfer */
    if (rq->evr.send_ready)
    {
        ret=tMPI_Xfer_recv(&rq->evr, 0, TRUE); 
    }

    /* we're probably going to do the actual transfers during Wait or Waitall, 
        so this one is non-blocking. */
    return MPI_SUCCESS;    
}

static int tMPI_Wait_r(struct mpi_req_ *rq, MPI_Status *status)
{
    int ret=MPI_SUCCESS;
    if (!rq)
        return MPI_SUCCESS;

#ifdef TMPI_DEBUG
    printf("%d: MPI_Wait on (%d -> %d, tag=%d)\n", 
           tMPI_Threadnr(NULL),
           tMPI_Threadnr( rq->recv ? rq->evr.src : rq->evs.src),
           tMPI_Threadnr( rq->recv ? rq->evr.dest : rq->evs.dest ),
           rq->recv ? rq->evr.tag : rq->evs.tag);
    fflush(0);
#endif


    if ( ! ( rq->finished || ( (( rq->recv) && rq->evr.finished)  ||
                               ((!rq->recv) && rq->evs.finished) ) ) )
    {
#ifdef TMPI_DEBUG
        printf("%d: MPI_Wait request was not finished\n", tMPI_Threadnr(0));
        fflush(0);
#endif
        if (rq->recv)
        {
            /* the receiving end just waits until the data is ready */
            tMPI_Wait_recv(&(rq->evr));
            ret=tMPI_Xfer_recv(&(rq->evr), status, FALSE);
        }
        else
        {
            /* we are now the waiting party: we need to initiate sending data.
               We do need to wait until sending is OK, though */
#ifdef TMPI_DEBUG
            printf(
            "%d: Active send waiting for recv_waiting to (%d->%d, tag %d)\n",
                    tMPI_Threadnr(0), 
                    tMPI_Threadnr(rq->recv ? rq->evr.src : rq->evs.src),
                    tMPI_Threadnr(rq->recv ? rq->evr.dest : rq->evs.dest),
                    rq->recv ? rq->evr.tag : rq->evs.tag);
            fflush(0);
#endif
            tMPI_Wait_send( &(rq->evs) );
#ifdef TMPI_DEBUG
            printf("%d: Sending data in the wrong direction (%d->%d) %d\n",
                    tMPI_Threadnr(0),
                    tMPI_Threadnr(rq->recv ? rq->evr.src : rq->evs.src),
                    tMPI_Threadnr(rq->recv ? rq->evr.dest : rq->evs.dest), 
                    rq->recv ? rq->evr.tag : rq->evs.tag);
#endif
            if (! rq->evs.finished)
            {
                ret=tMPI_Xfer_send(&(rq->evs), status, FALSE);
            }
        }
        rq->finished=TRUE;
    }

    if (status && rq->recv)
    {
        status->MPI_SOURCE=tMPI_Comm_seek_rank(rq->evr.comm, rq->evr.src);
        status->MPI_TAG=rq->evr.tag;
        status->MPI_ERROR=MPI_SUCCESS;
        status->transferred=rq->evr.transferred;
    }

    return ret;
}

static int tMPI_Test_r(struct mpi_req_ *rq, int *flag, MPI_Status *status)
{    
    int ret=MPI_SUCCESS;
    if (!rq)
        return MPI_SUCCESS;

    if ( ! ( rq->finished || ( (( rq->recv) && rq->evr.finished)  ||
                               ((!rq->recv) && rq->evs.finished) ) ) )
    {
       if (rq->recv)
       {
           /* the receiving end just checks  until the data is ready */
           if ( rq->evr.send_ready ) 
           {
               ret=tMPI_Xfer_recv(&(rq->evr), status, TRUE);
           }
       }
       else
       {
           /* only transfer if we really need to */
           if ( rq->evs.recv_nonblock && rq->evs.recv_ready )
           {
               ret=tMPI_Xfer_send(&(rq->evs), status, TRUE);
           }

       }
    }
    if (( rq->recv && rq->evr.finished) || ((!rq->recv) && rq->evs.finished) )
    {
        rq->finished=TRUE;
    }

    
    if  ( rq->finished )     
    {
        *flag=TRUE;
        if (status && rq->recv)
        {
            status->MPI_SOURCE=tMPI_Comm_seek_rank(rq->evr.comm, rq->evr.src);
            status->MPI_TAG=rq->evr.tag;
            status->MPI_ERROR=MPI_SUCCESS;
            status->transferred=rq->evr.transferred;
        }
    }
    else
        *flag=FALSE;
    return ret;
}

static int tMPI_Waitall_r(int count, struct mpi_req_ *array_of_requests[],
                          MPI_Status *array_of_statuses, bool may_free)
{
    int done;
    int flags_[MAX_THREADS];
    int *flags=flags_;
    int i;

    if (count > MAX_THREADS)
    {
        smalloc(flags, sizeof(int)*count);
    }

    for(i=0;i<count;i++)
        flags[i]=FALSE;
    /* Waitall polls all the requests by calling MPI_Test. This
       ensures that incoming receives are handled in the order that they
       come in, but of course is another busy-wait type function. */
    do
    {
        /* first take care of duds */
        for(i=0;i<count;i++)
        {
            if (!array_of_requests[i])
                flags[i]=TRUE;
        }
        /* do receives */
        for(i=0;i<count;i++)
        {
            MPI_Status *st=0; if (array_of_statuses) st=&(array_of_statuses[i]);
            if (!flags[i] && array_of_requests[i]->recv)
            {
                tMPI_Test_r(array_of_requests[i], &(flags[i]), st);
                if (flags[i])
                {
                    if (may_free)
                        sfree(array_of_requests[i]);
                    array_of_requests[i]=MPI_REQUEST_NULL;
                }
            }
        }
        /* then  do sends */
        for(i=0;i<count;i++)
        {
            MPI_Status *st=0; if (array_of_statuses) st=&(array_of_statuses[i]);
            if (!flags[i] && !(array_of_requests[i]->recv))
            {
                tMPI_Test_r(array_of_requests[i], &(flags[i]), st);
                if (flags[i])
                {
                    if (may_free)
                        sfree(array_of_requests[i]);
                    array_of_requests[i]=MPI_REQUEST_NULL;
                }
            }
        }
        /* count done flags */
        done=0;
        for(i=0;i<count;i++)
            if (flags[i])
                done++;
    }
    while (done<count);

    if (count > MAX_THREADS)
        sfree(flags);

    return MPI_SUCCESS;
}



/* the real MPI functions just call their tMPI_*r counterparts */

int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
             int tag, MPI_Comm comm)
{
    return MPI_Send_r(buf, count, datatype, dest, tag, comm, point_to_point);
}


int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
             int tag, MPI_Comm comm, MPI_Status *status)
{
    return MPI_Recv_r(buf, count, datatype, source, tag, comm, status, 
                      point_to_point);
}


int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, MPI_Request *request)
{
    struct mpi_req_ *rq;
    int ret;
    smalloc(rq, sizeof(struct mpi_req_));
    if ( (ret=tMPI_Isend_r(buf, count, datatype, dest, tag, comm, rq,
                           point_to_point)) != 
            MPI_SUCCESS)
    {
#ifdef TMPI_DEBUG
        printf("%d Freeing isend request after error\n", tMPI_Threadnr(0)); 
        fflush(0);
#endif
        sfree(rq);
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
    struct mpi_req_ *rq;
    int ret;
    smalloc(rq, sizeof(struct mpi_req_));
    if ( (ret=tMPI_Irecv_r(buf, count, datatype, source, tag, comm, rq,
                           point_to_point)) != 
            MPI_SUCCESS)
    {
#ifdef TMPI_DEBUG
        printf("%d Freeing irecv request after error\n", tMPI_Threadnr(0)); 
        fflush(0);
#endif
        sfree(rq);
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
    if (!request || !(*request))
        return MPI_SUCCESS;

    ret=tMPI_Wait_r(*request, status);

    /* deallocate if needed */
    sfree(*request);
    *request=MPI_REQUEST_NULL;
    return ret;
}

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{
    int ret=MPI_SUCCESS;
    if (!request || !(*request))
        return MPI_SUCCESS;

    ret=tMPI_Test_r(*request, flag, status);

    if ((*request)->finished)
    {
        /* deallocate if needed */
        sfree(*request);
        *request=MPI_REQUEST_NULL;
    }

    return ret;
}


int MPI_Waitall(int count, MPI_Request *array_of_requests,
                MPI_Status *array_of_statuses)
{
    return tMPI_Waitall_r(count, array_of_requests, array_of_statuses, TRUE);
}




int MPI_Barrier(MPI_Comm comm) 
{
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }

    gmx_thread_barrier_wait( &(comm->multicast_barrier[0]));
    return MPI_SUCCESS;
}


/* multi */
int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root,
              MPI_Comm comm)
{
    int myrank,ret=MPI_SUCCESS;
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!buffer)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (myrank==root)
    {
        int i;
        int N=tMPI_Comm_N(comm);
        struct mpi_req_ rqr_[MAX_THREADS]; 
        MPI_Request rq_[MAX_THREADS]; 
        struct mpi_req_ *rqr=rqr_; /* this is where we allocate the requests */
        MPI_Request *rq=rq_; /* pointers to the requests for tMPI_Waitall */

        if (N>MAX_THREADS)
        {
            smalloc(rqr, sizeof(struct mpi_req_)*N);
            smalloc(rq, sizeof(MPI_Request)*N);
        }

        for(i=0;i<N;i++)
        {
            if (i != myrank)
            {
                int ret2;
                if ((ret2=tMPI_Isend_r(buffer, count, datatype, i, 
                                TMPI_BCAST_TAG, comm, &(rqr[i]), multicast)) !=
                        MPI_SUCCESS)
                {
                    if (N>MAX_THREADS)
                    {
                        sfree(rqr);
                        sfree(rq);
                    }
                    return ret2;
                }
                rq[i]=&(rqr[i]);
            }
            else
            {
                rq[i]=0;
            }
        }
        ret=tMPI_Waitall_r(N, rq, 0, FALSE);

        if (N>MAX_THREADS)
        {
            sfree(rqr);
            sfree(rq);
        }
    }
    else
    {
        ret=MPI_Recv_r(buffer, count, datatype, root, TMPI_BCAST_TAG, comm, 0,
                       multicast);
    }
    return ret;
}



int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
               void* recvbuf, int recvcount, MPI_Datatype recvtype, int root,
               MPI_Comm comm)
{
    int i;
    int N=tMPI_Comm_N(comm);
    int displs_[MAX_THREADS];
    int recvcounts_[MAX_THREADS];
    int *displs=displs_;
    int *recvcounts=recvcounts_;
    int ret;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!sendbuf || !recvbuf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }


    if (N>MAX_THREADS)
    {
        smalloc(displs, sizeof(int)*N);
        smalloc(recvcounts, sizeof(int)*N);
    }

    for(i=0;i<N;i++)
    {
        recvcounts[i]=recvcount;
        displs[i]=recvcount*i;
    }
    ret=MPI_Gatherv(sendbuf, sendcount, sendtype,
                       recvbuf, recvcounts, displs, recvtype, root, comm);

    if (N>MAX_THREADS)
    {
        sfree(displs);
        sfree(recvcounts);
    }
    return ret;
}



int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int *recvcounts, int *displs,
                MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int myrank,ret=MPI_SUCCESS;
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (myrank==root)
    {
        int i;
        int N=tMPI_Comm_N(comm);
        struct mpi_req_ rqr_[MAX_THREADS]; 
        MPI_Request rq_[MAX_THREADS]; 
        struct mpi_req_ *rqr=rqr_; /* this is where we allocate the requests */
        MPI_Request *rq=rq_; /* pointers to the requests for tMPI_Waitall */

        if (!recvbuf)
        {
            comm->erh->err=MPI_ERR_BUF;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }

        if (N>MAX_THREADS)
        {
            smalloc(rqr, sizeof(struct mpi_req_)*N);
            smalloc(rq, sizeof(MPI_Request)*N);
        }

        for(i=0;i<N;i++)
        {
            if (i!=myrank)
            {
                int ret2;
                if ((ret2=tMPI_Irecv_r((char*)recvbuf+ 
                                       displs[i]*recvtype->size, 
                                       recvcounts[i], recvtype, i, 
                                       TMPI_GATHER_TAG, comm, &(rqr[i]),
                                       multicast)) !=
                        MPI_SUCCESS)
                {
                    if (N>MAX_THREADS)
                    {
                        sfree(rqr);
                        sfree(rq);
                    }
                    return ret2;
                }
                rq[i]=&(rqr[i]);
            }
            else
            {
                rq[i]=0;
            }
        }
        if (sendbuf) /* i.e.  ( sendbuf!=MPI_IN_PLACE) */
        {
            /* do the root transfer */
            int recvsize=recvtype->size*recvcounts[myrank];
            int sendsize=sendtype->size*sendcount;
            if (recvsize < sendsize)
            {
                if (N>MAX_THREADS)
                {
                    sfree(rqr);
                    sfree(rq);
                }
                comm->erh->err=MPI_ERR_XFER_BUFSIZE;
                comm->erh->fn(&comm, &(comm->erh->err));
                return comm->erh->err;
            }
            memcpy((char*)recvbuf + displs[myrank]*recvtype->size,
                   sendbuf, sendsize);
        }
        ret=tMPI_Waitall_r(N, rq, 0, FALSE);
        if (N>MAX_THREADS)
        {
            sfree(rqr);
            sfree(rq);
        }
    }
    else
    {
        if (!sendbuf)
        {
            comm->erh->err=MPI_ERR_BUF;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
        ret=MPI_Send_r(sendbuf, sendcount, sendtype, root, TMPI_GATHER_TAG, 
                       comm, multicast);
    }
 
    return ret;
}

int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, 
                MPI_Comm comm)
{
    int i;
    int N=tMPI_Comm_N(comm);
    int displs_[MAX_THREADS];
    int sendcounts_[MAX_THREADS];
    int *displs=displs_;
    int *sendcounts=sendcounts_;
    int ret;

    if (N>MAX_THREADS)
    {
        smalloc(displs, sizeof(struct mpi_req_)*N);
        smalloc(sendcounts, sizeof(MPI_Request)*N);
    }

    for(i=0;i<N;i++)
    {
        sendcounts[i]=sendcount;
        displs[i]=sendcount*i;
    }
    ret=MPI_Scatterv(sendbuf, sendcounts, displs, sendtype,
                        recvbuf, recvcount, recvtype, root, comm);

    if (N>MAX_THREADS)
    {
        sfree(displs);
        sfree(sendcounts);
    }


    return ret;
}



int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs,
                 MPI_Datatype sendtype, void* recvbuf, int recvcount,
                 MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    /* TODO: make a better version of this function. This 
       should probably be binary tree-based */
    int myrank,ret=MPI_SUCCESS;
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (myrank==root)
    {
        int i;
        int N=tMPI_Comm_N(comm);
        struct mpi_req_ rqr_[MAX_THREADS]; 
        MPI_Request rq_[MAX_THREADS]; 
        struct mpi_req_ *rqr=rqr_; /* this is where we allocate the requests */
        MPI_Request *rq=rq_; /* pointers to the requests for tMPI_Waitall */

        if (!sendbuf)
        {
            comm->erh->err=MPI_ERR_BUF;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }

        if (N>MAX_THREADS)
        {
            smalloc(rqr, sizeof(struct mpi_req_)*N);
            smalloc(rq, sizeof(MPI_Request)*N);
        }
        for(i=0;i<N;i++)
        {
            if (i!=myrank)
            {
                int ret2;
                if ((ret2=tMPI_Isend_r((char*)sendbuf+ displs[i]*sendtype->size,
                                sendcounts[i], sendtype, i, 
                                TMPI_SCATTER_TAG, comm, &(rqr[i]), multicast))
                        != MPI_SUCCESS)
                {
                    if (N>MAX_THREADS)
                    {
                        smalloc(rqr, sizeof(struct mpi_req_)*N);
                        smalloc(rq, sizeof(MPI_Request)*N);
                    }
                    return ret2;
                }
                rq[i]=&(rqr[i]);
            }
            else
            {
                rq[i]=0;
            }
        }
        if (recvbuf) /* i.e.  ( recvbuf!=MPI_IN_PLACE) */
        {
            /* do the root transfer */
            int recvsize=recvtype->size*recvcount;
            int sendsize=sendtype->size*sendcounts[myrank];
            if (recvsize < sendsize)
            {
                comm->erh->err=MPI_ERR_XFER_BUFSIZE;
                comm->erh->fn(&comm, &(comm->erh->err));
                return comm->erh->err;
            }
            memcpy(recvbuf, (char*)sendbuf+ displs[myrank]*sendtype->size,
                   sendsize);
        }
        ret=tMPI_Waitall_r(N, rq, 0, FALSE);
        if (N>MAX_THREADS)
        {
            smalloc(rqr, sizeof(struct mpi_req_)*N);
            smalloc(rq, sizeof(MPI_Request)*N);
        }
    }
    else
    {
        if (!recvbuf)
        {
            comm->erh->err=MPI_ERR_BUF;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
        ret=MPI_Recv_r(recvbuf, recvcount, recvtype, root, TMPI_SCATTER_TAG, 
                       comm, 0, multicast);
    }
 
    return ret;
}




int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm)

{
    int i;
    int N=tMPI_Comm_N(comm);
    int sdispls_[MAX_THREADS];
    int rdispls_[MAX_THREADS];
    int sendcounts_[MAX_THREADS];
    int recvcounts_[MAX_THREADS];
    int *sdispls=sdispls_;
    int *rdispls=rdispls_;
    int *sendcounts=sendcounts_;
    int *recvcounts=recvcounts_;
    int ret;

    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (N>MAX_THREADS)
    {
        smalloc(sdispls, N*sizeof(int));
        smalloc(rdispls, N*sizeof(int));
        smalloc(sendcounts, N*sizeof(int));
        smalloc(recvcounts, N*sizeof(int));
    }
 

    for(i=0;i<N;i++)
    {
        sendcounts[i]=sendcount;
        recvcounts[i]=recvcount;
        sdispls[i]=sendcount*i;
        rdispls[i]=recvcount*i;
    }
    ret= MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype,
                       recvbuf, recvcounts, rdispls, recvtype, 
                       comm);
    if (N>MAX_THREADS)
    {
        sfree(sdispls);
        sfree(rdispls);
        sfree(sendcounts);
        sfree(recvcounts);
    }
    return ret;
}


int MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype sendtype, void* recvbuf, int *recvcounts,
                  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)

{
    int i,myrank,ret=MPI_SUCCESS;
    int N=tMPI_Comm_N(comm);

    struct mpi_req_ rqs_[MAX_THREADS]; /* send requests */
    struct mpi_req_ rqr_[MAX_THREADS]; /* receive requests */
    /* pointers to the requests for tMPI_Waitall */
    MPI_Request rq_[2*MAX_THREADS]; 

    struct mpi_req_ *rqs=rqs_; /* send requests */
    struct mpi_req_ *rqr=rqr_; /* receive requests */
    MPI_Request *rq=rq_; /* pointers to the requests for tMPI_Waitall */
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
    if (!sendbuf || !recvbuf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
    myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    if (N>MAX_THREADS)
    {
        smalloc(rqs, sizeof(struct mpi_req_)*N);
        smalloc(rqr, sizeof(struct mpi_req_)*N);
        smalloc(rq, sizeof(MPI_Request)*N*2);
    }

    /* first post async sends */
    for(i=0;i<N;i++)
    {
        if (i!=myrank ) /*&& sendcounts[i]>0)*/
        {
            int ret2;
            if ( (ret2=tMPI_Isend_r((char*)sendbuf+ sdispls[i]*sendtype->size, 
                            sendcounts[i], sendtype, i, 
                            TMPI_ALLTOALL_TAG, comm, &(rqs[i]), multicast))
                    != MPI_SUCCESS)
            {
                return ret2;
            }
            rq[0*N+i]=&(rqs[i]);
        }
        else
        {
            rq[0*N+i]=0;
        }
    }

    /* do local transfer */
    {
        int recvsize=recvtype->size*recvcounts[myrank];
        int sendsize=sendtype->size*sendcounts[myrank];
        if (recvsize < sendsize)
        {
            comm->erh->err=MPI_ERR_XFER_BUFSIZE;
            comm->erh->fn(&comm, &(comm->erh->err));
            return comm->erh->err;
        }
        memcpy((char*)recvbuf + rdispls[myrank]*recvtype->size, 
               (char*)sendbuf + sdispls[myrank]*sendtype->size,
               sendsize);
    }

    /* now post async recvs */
    for(i=0;i<N;i++)
    {
        if (i!=myrank ) /*&& recvcounts[i]>0 )*/
        {
            int ret2;
            if ((ret2=tMPI_Irecv_r((char*)recvbuf+ rdispls[i]*recvtype->size, 
                            recvcounts[i], recvtype, i, TMPI_ALLTOALL_TAG, 
                            comm, &(rqr[i]), multicast)) !=
                    MPI_SUCCESS)
            {
                return ret2;
            }
            rq[1*N + i]=&(rqr[i]);
        }
        else
        {
            rq[1*N + i]=0;
        }
    }
    /* and wait for everything to complete */
    ret=tMPI_Waitall_r(2*N, rq, 0, FALSE);
    return ret;
}






void tMPI_Reduce_run_op(void *dest, void *src_a, void *src_b, 
                        MPI_Datatype datatype, int count, MPI_Op op)
{
    tMPI_Op_fn fn=datatype->op_functions[op];

    fn(dest, src_a, src_b, count);
}

int MPI_Reduce(void* sendbuf, void* recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
{
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

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
   
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
#if 0
    if (!sendbuf || recvbuf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
#endif
    if ( (!datatype->op_functions) || (!datatype->op_functions[op]) )
    {
        comm->erh->err=MPI_ERR_OP_FN;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }


    if (myrank==root)
    {
        if (!recvbuf) /* i.e. recvbuf == MPI_IN_PLACE */
        {
            recvbuf=sendbuf;
        }
    }
    comm->sendbuf[myrank]=sendbuf;
    comm->recvbuf[myrank]=recvbuf;
    /* there's a barrier to wait for all the processes to put their 
       send/recvbuf in the global list */
    gmx_thread_barrier_wait( &(comm->multicast_barrier[0]));

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
#ifdef TMPI_DEBUG
                printf("%d: reducing with %d, iteration=%d\n", 
                        myrank, myrank+nbr_dist, iteration);
                fflush(0);
#endif
               /* we reduce with our neighbour*/
                nbr=(N+myrank+nbr_dist)%N; /* here we use the real rank */
                void *a,*b;
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
                tMPI_Reduce_run_op(recvbuf, a, b, datatype, count, op);
            }
            else
            {
                /* we still need to put things in the right buffer for the next
                   iteration */
                if (iteration==0)
                    memcpy(recvbuf, sendbuf, datatype->size*count);
            }
            /* split barrier */
            gmx_thread_barrier_wait( &(comm->multicast_barrier[iteration]));
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
            gmx_thread_barrier_wait( &(comm->multicast_barrier[iteration]));
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


int MPI_Allreduce(void* sendbuf, void* recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    void *rootbuf=0;
    int myrank=tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    /* this function uses a binary tree-like reduction algorithm: */
    int N=tMPI_Comm_N(comm);
    int Nnbrs=N; /* number of neighbours to communicate with 
                    (decreases exponentially) */
    int nbr_dist=1; /* distance between communicating neighbours 
                       (increases exponentially) */
    int stepping=2; /* distance between non-communicating neighbours
                       (increases exponentially) */
    int iteration=0;
   
    if (!comm)
    {
        MPI_COMM_WORLD->erh->err=MPI_ERR_COMM;
        MPI_COMM_WORLD->erh->fn(&MPI_COMM_WORLD, &(MPI_COMM_WORLD->erh->err));
        return MPI_COMM_WORLD->erh->err;
    }
#if 0
    if (!sendbuf || recvbuf)
    {
        comm->erh->err=MPI_ERR_BUF;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }
#endif
    if ( (!datatype->op_functions) || (!datatype->op_functions[op]) )
    {
        comm->erh->err=MPI_ERR_OP_FN;
        comm->erh->fn(&comm, &(comm->erh->err));
        return comm->erh->err;
    }

    if (myrank==0)
    {
        if (!recvbuf) /* i.e. recvbuf == MPI_IN_PLACE */
        {
            recvbuf=sendbuf;
        }
        rootbuf=recvbuf;
    }
    comm->sendbuf[myrank]=sendbuf;
    comm->recvbuf[myrank]=recvbuf;
    /* there's a barrier to wait for all the processes to put their 
       send/recvbuf in the global list */
    gmx_thread_barrier_wait( &(comm->multicast_barrier[0]));

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
        if (myrank%stepping == 0 )
        {
            if (myrank+nbr_dist<N) 
            {
#ifdef TMPI_DEBUG
                printf("%d: reducing with %d, iteration=%d\n", 
                        myrank, myrank+nbr_dist, iteration);
                fflush(0);
#endif
               /* we reduce with our neighbour*/
                nbr=(N+myrank+nbr_dist)%N; /* here we use the real rank */
                void *a,*b;
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
                tMPI_Reduce_run_op(recvbuf, a, b, datatype, count, op);
            }
            else
            {
                /* we still need to put things in the right buffer for the next
                   iteration */
                if (iteration==0)
                    memcpy(recvbuf, sendbuf, datatype->size*count);
            }
            /* split barrier */
            gmx_thread_barrier_wait( &(comm->multicast_barrier[iteration]));
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
            gmx_thread_barrier_wait( &(comm->multicast_barrier[iteration]));
            break;
        }
#ifdef TMPI_DEBUG
        printf("%d: barrier over, iteration=%d\n", myrank,  iteration);
        fflush(0);
#endif
        Nnbrs = Nnbrs/2 + Nnbrs%2;
        nbr_dist*=2;
        stepping*=2;
        /* there's a barrier to wait for all the processes to put their 
           send/recvbuf in the global list */
        iteration++;
    }

    gmx_thread_barrier_wait( &(comm->multicast_barrier[0]));
    /* and now we just copy things back inefficiently. We should make
       a better MPI_Scatter, and use that. */
    if (myrank != 0)
    {
        memcpy(recvbuf, rootbuf, datatype->size*count );
    }
    gmx_thread_barrier_wait( &(comm->multicast_barrier[0]));

    return MPI_SUCCESS;
}

#endif /* GMX_THREAD_MPI */
