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

#ifndef _TMPI_H_
#define _TMPI_H_

/** \file thread_mpi.h
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
 * 
*/


#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif

#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"


/* The tMPI_Comm structure contains the group of processes to communicate
   with (defines the scope for global operations such as broadcast) */
typedef struct tmpi_comm_ *tMPI_Comm;
/* The group part of the MPI-Comm structure */
typedef struct tmpi_group_ *tMPI_Group;
/* Request structure for holding data about non-blocking transfers */
typedef struct tmpi_req_ *tMPI_Request;
/* status of receives */
typedef struct tmpi_status_ tMPI_Status;
/* data types */
typedef struct tmpi_datatype_ *tMPI_Datatype;


/** MPI data types as specified by the MPI standard. 
    Note that not all are available.  */
extern tMPI_Datatype TMPI_CHAR;
extern tMPI_Datatype TMPI_SHORT;
extern tMPI_Datatype TMPI_INT;
extern tMPI_Datatype TMPI_LONG;
#ifdef SIZEOF_LONG_LONG_INT
extern tMPI_Datatype TMPI_LONG_LONG;
extern tMPI_Datatype TMPI_LONG_LONG_INT;
#endif

extern tMPI_Datatype TMPI_SIGNED_CHAR;
extern tMPI_Datatype TMPI_UNSIGNED_CHAR;
extern tMPI_Datatype TMPI_UNSIGNED_SHORT;
extern tMPI_Datatype TMPI_UNSIGNED;
extern tMPI_Datatype TMPI_UNSIGNED_LONG;
#ifdef SIZEOF_LONG_LONG_INT
extern tMPI_Datatype TMPI_UNSIGNED_LONG_LONG;
#endif

extern tMPI_Datatype TMPI_FLOAT;
extern tMPI_Datatype TMPI_DOUBLE;
extern tMPI_Datatype TMPI_LONG_DOUBLE;

/*extern tMPI_Datatype tMPI_UNSIGNED_WCHAR */
extern tMPI_Datatype TMPI_BYTE;


/* error codes */
enum
{
    TMPI_SUCCESS=0,
    TMPI_ERR_MALLOC,
    TMPI_ERR_INIT,
    TMPI_ERR_FINALIZE,
    TMPI_ERR_GROUP,
    TMPI_ERR_COMM,
    TMPI_ERR_STATUS,
    TMPI_ERR_GROUP_RANK,
    TMPI_ERR_DIMS,
    TMPI_ERR_COORDS,
    TMPI_ERR_CART_CREATE_NPROCS,
    TMPI_ERR_XFER_COUNTERPART,
    TMPI_ERR_XFER_BUFSIZE,
    TMPI_ERR_XFER_BUF_OVERLAP,
    TMPI_ERR_SEND_DEST,
    TMPI_ERR_RECV_SRC,
    TMPI_ERR_BUF,
    TMPI_ERR_MULTI_MISMATCH,
    TMPI_ERR_OP_FN,
    TMPI_ERR_ENVELOPES,
    TMPI_ERR_REQUESTS,
    TMPI_FAILURE,
    TMPI_ERR_UNKNOWN,
    N_TMPI_ERR  /* this must be the last one */
};

#define TMPI_MAX_ERROR_STRING            256

#define TMPI_UNDEFINED -1

/* error handling */
typedef void (*tMPI_Errhandler_fn)(tMPI_Comm*, int*);
typedef struct tmpi_errhandler_ *tMPI_Errhandler;

extern tMPI_Errhandler TMPI_ERRORS_ARE_FATAL;
extern tMPI_Errhandler TMPI_ERRORS_RETURN;


/* miscelaneous defines */
#define TMPI_ANY_SOURCE -1
#define TMPI_ANY_TAG -1

/* topology test defines */
#define TMPI_CART 1
#define TMPI_GRAPH 2


/** All communicators */
extern tMPI_Comm TMPI_COMM_WORLD;
/* these are 0 instead of NULL so that we can compare against them */
#define TMPI_COMM_NULL 0
#define TMPI_GROUP_NULL 0

/** empty group */
extern tMPI_Group tMPI_GROUP_EMPTY;


#define TMPI_MAX_PROCESSOR_NAME 128


/* MPI status */
#define TMPI_STATUS_IGNORE 0
#define TMPI_STATUSES_IGNORE 0

/* the status object is user-maintained. */
struct tmpi_status_
{
    int TMPI_SOURCE; /* the message source rank */
    int TMPI_TAG; /* the message source tag */
    int TMPI_ERROR; /* the message error */
    int transferred;
};


#define TMPI_REQUEST_NULL 0

/* collective communication specials: */
#define TMPI_IN_PLACE 0


/**  tMPI_Reduce operators.
    These all work (except obviously bad combinations like bitwise
    and/or/xor on floats,etc): */
typedef enum
{
    TMPI_MAX, /* maximum */
    TMPI_MIN, /* minimum */
    TMPI_SUM, /* sum */
    TMPI_PROD, /* product */
    TMPI_LAND, /* logical and */
    TMPI_BAND, /* binary and */
    TMPI_LOR, /* logical or */
    TMPI_BOR, /* binary or */
    TMPI_LXOR, /* logical xor */
    TMPI_BXOR /* binary xor */
} tMPI_Op;


/* function for tMPI_COMM_SELF */
tMPI_Comm tMPI_Get_comm_self(void);
/* this must be a funtion because it's a thread-local property: */
#define TMPI_COMM_SELF (tMPI_Get_comm_self())


/** traditional MPI initializer; spawns threads that start at main(). 
  
    Seeks the argument '-np n', where n is the number of 
    threads that will be created. These new threads then run main() again,
    with the original argc and argv. 
    
    \param argc     argc of original main() invocation, or NULL
    \param argv     argv of original main() invocation, or NULL. 

    \return  MPI_SUCCESS on success, MPI_FAILURE on failure.  */
int tMPI_Init(int *argc, char ***argv);

/** Alternate thread MPI intializer and thread spawner.
  
    Creates N threads (including main thread) 
    that run the function start_function, which takes a void* argument, 
    given by arg. The function start_function also gets called by the main
    thread. When the function start_function returns it, will behave 
    as if tMPI_Finalize is called, and if it's a sub-thread it will
    stop running. 
    
    \param N                the number of threads to start.
    \param start_function   the function to start threads at (including 
                            main thread). 
    \param arg              an optional argument for start_function().

    \return  MPI_FAILURE on failure, MPI_SUCCESS on succes, after all
             threads have finished.
    */
int tMPI_Init_fn(int N, void (*start_function)(void*), void *arg);

/** get the number of threads that will be requested (can be called before 
    tMPI_Init() )
    
    \return  MPI_SUCCESS on success, MPI_FAILURE on failure.  */
int tMPI_Get_N(int *argc, char ***argv);


/* mostly for debugging */
/** get the number of this thread. 

    \return the thread number. */
#define tMPI_This_threadnr() (tMPI_Get_current() - threads)


/** waits for all other threads to finish and cleans up 

    \return  MPI_SUCCESS on success, MPI_FAILURE on failure.  */
int tMPI_Finalize(void);


/** just kills all threads. 
  
    Not really neccesary because exit() would do that for us anyway 
    
    \return Never returns. */
int tMPI_Abort(tMPI_Comm comm, int errorcode);
/** whether tMPI_Init, but not yet tMPI_Finalize, has been run

    \param flag     set to TRUE if tMPI_Init() has been called, FALSE if not.
    
    \return     always returns MPI_SUCCESS. */
int tMPI_Initialized(int *flag);
/** whether tMPI_Finalize has been run.

    \param flag     set to TRUE if tMPI_Finalize() has been called, FALSE 
                    if not.
     
    \return     always returns MPI_SUCCESS.  */
int tMPI_Finalized(int *flag);


/** Create an error handler object from a function.

    \param function     the function to make an error handler of.
    \param errhandler   the error handler.

    \return  MPI_SUCCESS on success, MPI_FAILURE on failure.  */
int tMPI_Create_errhandler(tMPI_Errhandler_fn *function,
                           tMPI_Errhandler *errhandler);


/** free the error handler object */
int tMPI_Errhandler_free(tMPI_Errhandler *errhandler);

/** set the error handler */
int tMPI_Comm_set_errhandler(tMPI_Comm comm, tMPI_Errhandler errhandler);
/** get the error handler */
int tMPI_Comm_get_errhandler(tMPI_Comm comm, tMPI_Errhandler *errhandler);

/** get the error string associated with an error code */
int tMPI_Error_string(int errorcode, char *string, size_t *resultlen);




/* system query: */
/** returns string with thread # */
int tMPI_Get_processor_name(char *name, size_t *resultlen);
/** get an elapsed time value as a double, in seconds */
double tMPI_Wtime(void);
#if 0
/** get the resolution of tMPI_Wtime as a double, in seconds */
double tMPI_Wtick(void);
#endif




/** check the size of a group */
int tMPI_Group_size(tMPI_Group group, int *size);
/** check the rank of a group */
int tMPI_Group_rank(tMPI_Group group, int *rank);
/** create a new group as a union of an existing group and new ranks*/
int tMPI_Group_incl(tMPI_Group group, int n, int *ranks, tMPI_Group *newgroup);
/** get a pointer to the group in the comm */
int tMPI_Comm_group(tMPI_Comm comm, tMPI_Group *group);
/** de-allocate a group */
int tMPI_Group_free(tMPI_Group *group);

/** get the comm size */
int tMPI_Comm_size(tMPI_Comm comm, int *size);
/** get the rank in comm of the current process */
int tMPI_Comm_rank(tMPI_Comm comm, int *rank);
/** de-allocate a comm */
int tMPI_Comm_free(tMPI_Comm *comm);
/** create a comm based on a group */
int tMPI_Comm_create(tMPI_Comm comm, tMPI_Group group, tMPI_Comm *newcomm);
/** split up a group into same-colored sub-groups ordered by key */
int tMPI_Comm_split(tMPI_Comm comm, int color, int key, tMPI_Comm *newcomm);
/** make a duplicate of a comm*/
int tMPI_Comm_dup(tMPI_Comm comm, tMPI_Comm *newcomm);

/* topology functions */
/** check what type of topology the comm has */
int tMPI_Topo_test(tMPI_Comm comm, int *status);
/** check which dimensionality a topology has */
int tMPI_Cartdim_get(tMPI_Comm comm, int *ndims);
/** check which size and pbc a Cartesian topology has */
int tMPI_Cart_get(tMPI_Comm comm, int maxdims, int *dims, int *periods, 
                 int *coords);
/** check which rank a set of process coordinates has in a Cartesian topology */
int tMPI_Cart_rank(tMPI_Comm comm, int *coords, int *rank);
/** check which coordinates a process rank has in a Cartesian topology */
int tMPI_Cart_coords(tMPI_Comm comm, int rank, int maxdims, int *coords);
/** check which rank this process would have in a Cartesian topology */
int tMPI_Cart_map(tMPI_Comm comm, int ndims, int *dims, int *periods, 
                         int *newrank);
/** create a comm with a Cartesian topology */
int tMPI_Cart_create(tMPI_Comm comm_old, int ndims, int *dims, int *periods, 
                    int reorder, tMPI_Comm *comm_cart);


/** create a contiguous data type (the only type possible right now */
int tMPI_Type_contiguous(int count, tMPI_Datatype oldtype, 
                        tMPI_Datatype *newtype);
/** make the data type ready for use */
int tMPI_Type_commit(tMPI_Datatype *datatype);



/** wait for all process in comm to arrive here */
int tMPI_Barrier(tMPI_Comm comm);



/** blocking transfers. The actual transfer (copy) is done on the receiving end 
    (so that the receiver's cache already contains the data that it presumably
     will use soon).  */
/* send message; waits until finished.  */
int tMPI_Send(void* buf, int count, tMPI_Datatype datatype, int dest, 
             int tag, tMPI_Comm comm);
/** receive message; waits until finished.  */
int tMPI_Recv(void* buf, int count, tMPI_Datatype datatype, int source, 
             int tag, tMPI_Comm comm, tMPI_Status *status);
/** send & receive message at the same time; waits until finished.  */
int tMPI_Sendrecv(void *sendbuf, int sendcount, tMPI_Datatype sendtype, 
                 int dest, int sendtag, void *recvbuf, int recvcount, 
                 tMPI_Datatype recvtype, int source, int recvtag, tMPI_Comm comm, 
                 tMPI_Status *status);
/** get the number of actually transferred items from a transfer status */
int tMPI_Get_count(tMPI_Status *status, tMPI_Datatype datatype, int *count);


/** async send/recv. The actual transfer is usually done on the receiving 
    end, during tMPI_Wait, tMPI_Waitall or tMPI_Test. For tMPI_Waitall, 
    the incoming messages are processed in the order they come in. 
   
    In the case of async receives, the sender may initiate transfer,
    and there's a lock in the envelope to make sure that it doesn't
    happen on both ends simultaneously. */
/** initiate sending a message */
int tMPI_Isend(void* buf, int count, tMPI_Datatype datatype, int dest, 
              int tag, tMPI_Comm comm, tMPI_Request *request);
/** initiate receiving a message */
int tMPI_Irecv(void* buf, int count, tMPI_Datatype datatype, int source, 
              int tag, tMPI_Comm comm, tMPI_Request *request);
/** test whether message is sent */
int tMPI_Test(tMPI_Request *request, int *flag, tMPI_Status *status);
/** wait until message is sent */
int tMPI_Wait(tMPI_Request *request, tMPI_Status *status);
/** wait for several message sending requests */
int tMPI_Waitall(int count, tMPI_Request *array_of_requests, 
                tMPI_Status *array_of_statuses);





/** multicast */
/** broadcast over entire comm from root */
int tMPI_Bcast(void* buffer, int count, tMPI_Datatype datatype, int root, 
              tMPI_Comm comm);

/** gather data from all processes in comm to root */
int tMPI_Gather(void* sendbuf, int sendcount, tMPI_Datatype sendtype, 
               void* recvbuf, int recvcount, tMPI_Datatype recvtype, int root, 
               tMPI_Comm comm);
/** gather irregularly laid out data from all processes in comm to root */
int tMPI_Gatherv(void* sendbuf, int sendcount, tMPI_Datatype sendtype, 
                void* recvbuf, int *recvcounts, int *displs, 
                tMPI_Datatype recvtype, int root, tMPI_Comm comm);

/** spread parts of sendbuf to all processes in comm from root */
int tMPI_Scatter(void* sendbuf, int sendcount, tMPI_Datatype sendtype, 
                void* recvbuf, int recvcount, tMPI_Datatype recvtype, int root, 
                tMPI_Comm comm);
/** spread irregularly laid out parts of sendbuf to all processes from root */
int tMPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
                 tMPI_Datatype sendtype, void* recvbuf, int recvcount, 
                 tMPI_Datatype recvtype, int root, tMPI_Comm comm); 




/** spread out parts of sendbuf to all processes from all processes */
int tMPI_Alltoall(void* sendbuf, int sendcount, tMPI_Datatype sendtype, 
                 void* recvbuf, int recvcount, tMPI_Datatype recvtype, 
                 tMPI_Comm comm);
/** spread out irregularly laid out parts of sendbuf to all processes 
    from all processes */
int tMPI_Alltoallv(void* sendbuf, int *sendcounts, int *sdispls, 
                  tMPI_Datatype sendtype, void* recvbuf, int *recvcounts, 
                  int *rdispls, tMPI_Datatype recvtype, tMPI_Comm comm);




/** Do an operation between all locally held buffers on all items in the 
    buffers and send the results to root*/
int tMPI_Reduce(void* sendbuf, void* recvbuf, int count, 
               tMPI_Datatype datatype, tMPI_Op op, int root, tMPI_Comm comm);
/** Do an tMPI_Reduce, but with the following assumption:
    recvbuf points to a valid buffer in all calling threads, or 
    sendbuf has the value tMPI_IN_PLACE (in which case the values of 
    sendbuf may be changed in that thread).  
    This avoids unnecesary memory allocations associated with the normal
    tMPI_Reduce. */
int tMPI_Reduce_fast(void* sendbuf, void* recvbuf, int count, 
                     tMPI_Datatype datatype, tMPI_Op op, int root, 
                     tMPI_Comm comm);

/** Do an operation between all locally held buffers on all items in the 
    buffers and broadcast the results */
int tMPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
                  tMPI_Datatype datatype, tMPI_Op op, tMPI_Comm comm);

#ifdef __cplusplus
} /* closing extern "C" */
#endif

#endif /* _TMPI_H_ */
