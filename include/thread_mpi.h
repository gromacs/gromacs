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
#ifndef _GMX_THREAD_MPI_H_
#define _GMX_THREAD_MPI_H_

/** \file gmx_thread_mpi.h
 *
 * \brief Partial implementation of MPI using only threads. 
 * 
 * See the MPI specification at
 * http://www.mpi-forum.org/docs/docs.html
 * for an explanation of what these functions do.
 *
 * Because this is a thread-based library, be very careful with global
 * variables and static variables in functions: they will be shared across
 * all threads an lead to conflicts if not properly mutex-ed or barrier-ed
 * out.
 *
 * This library supports all of MPI that is being used in Gromacs, but can 
 * still use some improvement:
 * - the gmx_mutexes should be replaced by busy-waits on atomic operations
 *   for performance reasons (the aim of a pthreads mutex: scheduling out
 *   waiting threads, is antithetical to the requirements of Gromacs: low 
 *   latency and high throughput).
 * - Some of the global communication functions (bcast, scatter, alltoall) 
 *   could perhaps use a binary tree-like distribution method rather than 
 *   simply letting each receiver thread read from one distributor.
 * 
 * Right now, this library can only be enabled using cmake (although some
 * work has been done on autoconf). The relevant option is GMX_THREADED.
 * 
*/

#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif


/* The MPI_Comm structure contains the group of processes to communicate
   with (defines the scope for global operations such as broadcast) */
typedef struct mpi_comm_ *MPI_Comm;
/* The group part of the MPI-Comm structure */
typedef struct mpi_group_ *MPI_Group;
/* Request structure for holding data about non-blocking transfers */
typedef struct mpi_req_ *MPI_Request;
/* status of receives */
typedef struct mpi_status_ MPI_Status;
/* data types */
typedef struct mpi_datatype_ *MPI_Datatype;


/** MPI data types as specified by the MPI standard. 
    Note that not all are available.  */
extern MPI_Datatype MPI_CHAR;
extern MPI_Datatype MPI_SHORT;
extern MPI_Datatype MPI_INT;
extern MPI_Datatype MPI_LONG;
#ifdef SIZEOF_LONG_LONG_INT
extern MPI_Datatype MPI_LONG_LONG;
extern MPI_Datatype MPI_LONG_LONG_INT;
#endif

extern MPI_Datatype MPI_SIGNED_CHAR;
extern MPI_Datatype MPI_UNSIGNED_CHAR;
extern MPI_Datatype MPI_UNSIGNED_SHORT;
extern MPI_Datatype MPI_UNSIGNED;
extern MPI_Datatype MPI_UNSIGNED_LONG;
#ifdef SIZEOF_LONG_LONG_INT
extern MPI_Datatype MPI_UNSIGNED_LONG_LONG;
#endif

extern MPI_Datatype MPI_FLOAT;
extern MPI_Datatype MPI_DOUBLE;
extern MPI_Datatype MPI_LONG_DOUBLE;

/*extern MPI_Datatype MPI_UNSIGNED_WCHAR */
extern MPI_Datatype MPI_BYTE;


/* error codes */
#define MPI_SUCCESS                     0
#define MPI_ERR_INIT                    1
#define MPI_ERR_FINALIZE                2
#define MPI_ERR_GROUP                   3
#define MPI_ERR_COMM                    4
#define MPI_ERR_STATUS                  5
#define MPI_ERR_GROUP_RANK              6
#define MPI_ERR_DIMS                    7
#define MPI_ERR_COORDS                  8
#define MPI_ERR_CART_CREATE_NPROCS      9
#define MPI_ERR_XFER_COUNTERPART        10
#define MPI_ERR_XFER_BUFSIZE            11
#define MPI_ERR_XFER_BUF_OVERLAP        12
#define MPI_ERR_SEND_DEST               13
#define MPI_ERR_RECV_SRC                14
#define MPI_ERR_BUF                     15
#define MPI_ERR_MULTI_MISMATCH          16
#define MPI_ERR_OP_FN                   17
#define MPI_ERR_ENVELOPES               18
#define MPI_ERR_REQUESTS                19
#define MPI_FAILURE                     20
#define MPI_ERR_UNKNOWN                 21

#define N_MPI_ERR                       22

#define MPI_MAX_ERROR_STRING            256

#define MPI_UNDEFINED -1

/* error handling */
typedef void (*MPI_Errhandler_fn)(MPI_Comm*, int*);
typedef struct mpi_errhandler_ *MPI_Errhandler;

extern MPI_Errhandler MPI_ERRORS_ARE_FATAL;
extern MPI_Errhandler MPI_ERRORS_RETURN;


/* miscelaneous defines */
#define MPI_ANY_SOURCE -1
#define MPI_ANY_TAG -1

/* topology test defines */
#define MPI_CART 1
#define MPI_GRAPH 2


/** All communicators */
extern MPI_Comm MPI_COMM_WORLD;
/* these are 0 instead of NULL so that we can compare against them */
#define MPI_COMM_NULL 0
#define MPI_GROUP_NULL 0

/** empty group */
extern MPI_Group MPI_GROUP_EMPTY;


#define MPI_MAX_PROCESSOR_NAME 128


/* MPI status */
#define MPI_STATUS_IGNORE 0
#define MPI_STATUSES_IGNORE 0

/* the status object is user-maintained. */
struct mpi_status_
{
    int MPI_SOURCE; /* the message source rank */
    int MPI_TAG; /* the message source tag */
    int MPI_ERROR; /* the message error */
    int transferred;
};


#define MPI_REQUEST_NULL 0

/* collective communication specials: */
#define MPI_IN_PLACE 0


/**  MPI_Reduce operators.
    These all work (except obviously bad combinations like bitwise
    and/or/xor on floats,etc): */
typedef enum
{
    MPI_MAX, /* maximum */
    MPI_MIN, /* minimum */
    MPI_SUM, /* sum */
    MPI_PROD, /* product */
    MPI_LAND, /* logical and */
    MPI_BAND, /* binary and */
    MPI_LOR, /* logical or */
    MPI_BOR, /* binary or */
    MPI_LXOR, /* logical xor */
    MPI_BXOR /* binary xor */
} MPI_Op;


/* function for MPI_COMM_SELF */
MPI_Comm tMPI_Get_comm_self(void);
/* this must be a funtion because it's a thread-local property: */
#define MPI_COMM_SELF (tMPI_Get_comm_self())


/** MPI initializer. Seeks the argument '-np n', where n is the number of 
    threads that will be created. These new threads then run main() again,
    with the original argc and argv. */
int MPI_Init(int *argc, char ***argv);

/** Alternate thread MPI intializer. Creates N threads (including main thread) 
    that run the function start_function, which takes a void* argument, 
    given by arg. The function start_function also gets called by the main
    thread. When the function start_function returns it, will behave 
    as if MPI_Finalize is called, and if it's a sub-thread it will
    stop running. */
int tMPI_Init_fn(int N, void (*start_function)(void*), void *arg);

/** get the number of threads that will be requested (can be called before 
    MPI_Init() ) */
int tMPI_Get_N(int *argc, char ***argv);


/** mostly for debugging */
struct mpi_thread;
unsigned int tMPI_Threadnr(struct mpi_thread *thr);
unsigned int tMPI_This_threadnr(void);


/** waits for all threads to join() */
int MPI_Finalize(void);
/** just kills all threads. Not really neccesary because exit() will do 
    that for us anyway */
int MPI_Abort(MPI_Comm comm, int errorcode);
/** whether MPI_Init, but not yet MPI_Finalize, has been run*/
int MPI_Initialized(int *flag);
/** whether MPI_Finalize has been run */
int MPI_Finalized(int *flag);


/** create an error handler object from a function */
int MPI_Create_errhandler(MPI_Errhandler_fn *function,
                          MPI_Errhandler *errhandler);
/** free the error handler object */
int MPI_Errhandler_free(MPI_Errhandler *errhandler);

/** set the error handler */
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
/** get the error handler */
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);

/** get the error string associated with an error code */
int MPI_Error_string(int errorcode, char *string, int *resultlen);




/* system query: */
/** returns string with thread # */
int MPI_Get_processor_name(char *name, int *resultlen);
/** get an elapsed time value as a double, in seconds */
double MPI_Wtime(void);
#if 0
/** get the resolution of MPI_Wtime as a double, in seconds */
double MPI_Wtick(void);
#endif




/** check the size of a group */
int MPI_Group_size(MPI_Group group, int *size);
/** check the rank of a group */
int MPI_Group_rank(MPI_Group group, int *rank);
/** create a new group as a union of an existing group and new ranks*/
int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
/** get a pointer to the group in the comm */
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
/** de-allocate a group */
int MPI_Group_free(MPI_Group *group);

/** get the comm size */
int MPI_Comm_size(MPI_Comm comm, int *size);
/** get the rank in comm of the current process */
int MPI_Comm_rank(MPI_Comm comm, int *rank);
/** de-allocate a comm */
int MPI_Comm_free(MPI_Comm *comm);
/** create a comm based on a group */
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
/** split up a group into same-colored sub-groups ordered by key */
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
/** make a duplicate of a comm*/
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);

/* topology functions */
/** check what type of topology the comm has */
int MPI_Topo_test(MPI_Comm comm, int status);
/** check which dimensionality a topology has */
int MPI_Cartdim_get(MPI_Comm comm, int *ndims);
/** check which size and pbc a Cartesian topology has */
int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods, 
                 int *coords);
/** check which rank a set of process coordinates has in a Cartesian topology */
int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);
/** check which coordinates a process rank has in a Cartesian topology */
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int *coords);
/** check which rank this process would have in a Cartesian topology */
int MPI_Cart_map(MPI_Comm comm, int ndims, int *dims, int *periods, 
                         int *newrank);
/** create a comm with a Cartesian topology */
int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods, 
                    int reorder, MPI_Comm *comm_cart);


/** create a contiguous data type (the only type possible right now */
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, 
                        MPI_Datatype *newtype);
/** make the data type ready for use */
int MPI_Type_commit(MPI_Datatype *datatype);



/** wait for all process in comm to arrive here */
int MPI_Barrier(MPI_Comm comm);



/** blocking transfers. The actual transfer (copy) is done on the receiving end 
    (so that the receiver's cache already contains the data that it presumably
     will use soon).  */
/* send message; waits until finished.  */
int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest, 
             int tag, MPI_Comm comm);
/** receive message; waits until finished.  */
int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, 
             int tag, MPI_Comm comm, MPI_Status *status);
/** send & receive message at the same time; waits until finished.  */
int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                 int dest, int sendtag, void *recvbuf, int recvcount, 
                 MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, 
                 MPI_Status *status);
/** get the number of actually transferred items from a transfer status */
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);


/** async send/recv. The actual transfer is usually done on the receiving 
    end, during MPI_Wait, MPI_Waitall or MPI_Test. For MPI_Waitall, 
    the incoming messages are processed in the order they come in. 
   
    In the case of async receives, the sender may initiate transfer,
    and there's a lock in the envelope to make sure that it doesn't
    happen on both ends simultaneously. */
/** initiate sending a message */
int MPI_Isend(void* buf, int count, MPI_Datatype datatype, int dest, 
              int tag, MPI_Comm comm, MPI_Request *request);
/** initiate receiving a message */
int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source, 
              int tag, MPI_Comm comm, MPI_Request *request);
/** test whether message is sent */
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
/** wait until message is sent */
int MPI_Wait(MPI_Request *request, MPI_Status *status);
/** wait for several message sending requests */
int MPI_Waitall(int count, MPI_Request *array_of_requests, 
                MPI_Status *array_of_statuses);





/** multicast */
/** broadcast over entire comm from root */
int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype, int root, 
              MPI_Comm comm);

/** gather data from all processes in comm to root */
int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
               void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, 
               MPI_Comm comm);
/** gather irregularly laid out data from all processes in comm to root */
int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                void* recvbuf, int *recvcounts, int *displs, 
                MPI_Datatype recvtype, int root, MPI_Comm comm);

/** spread parts of sendbuf to all processes in comm from root */
int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, 
                MPI_Comm comm);
/** spread irregularly laid out parts of sendbuf to all processes from root */
int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
                 MPI_Datatype sendtype, void* recvbuf, int recvcount, 
                 MPI_Datatype recvtype, int root, MPI_Comm comm); 




/** spread out parts of sendbuf to all processes from all processes */
int MPI_Alltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                 void* recvbuf, int recvcount, MPI_Datatype recvtype, 
                 MPI_Comm comm);
/** spread out irregularly laid out parts of sendbuf to all processes 
    from all processes */
int MPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, 
                  MPI_Datatype sendtype, void* recvbuf, int *recvcounts, 
                  int *rdispls, MPI_Datatype recvtype, MPI_Comm comm);




/** Do an operation between all locally held buffers on all items in the 
    buffers and send the results to root*/
int MPI_Reduce(void* sendbuf, void* recvbuf, int count, 
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
/** Do an MPI_Reduce, but with the following assumption:
    recvbuf points to a valid buffer in all calling threads, or 
    sendbuf has the value MPI_IN_PLACE (in which case the values of 
    sendbuf may be changed in that thread).  
    This avoids unnecesary memory allocations associated with the normal
    MPI_Reduce. */
int tMPI_Reduce_fast(void* sendbuf, void* recvbuf, int count, 
                     MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

/** Do an operation between all locally held buffers on all items in the 
    buffers and broadcast the results */
int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

#ifdef __cplusplus
} /* closing extern "C" */
#endif

#endif /* _THREAD_MPI_H_ */
