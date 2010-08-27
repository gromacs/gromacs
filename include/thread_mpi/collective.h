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

#ifndef _THREAD_MPI_COLLECTIVE_H_
#define _THREAD_MPI_COLLECTIVE_H_

/** \file 
 *
 * \brief Collective functions
 * 
*/


#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif

/** Execute function once over comm
    
    Executes a given function only once per collective call over comm. 

    Collective function.

    \param[in] function     the function to call
    \param[in] param        the parameter to the function
    \param[out] was_first   set to 1 if the current thread was the one to 
                            execute the function. unchanged if current thread
                            was not the one to execute the function; ignored
                            if NULL.
    \param[in] comm         The communicator.
    \returns MPI_SUCCESS on success.
*/
int tMPI_Once(tMPI_Comm comm,void (*function)(void*), void *param, 
              int *was_first);

/** Execute function once over comm, and wait for the function to return
    its value.
    
    Executes a given function only once per collective call over comm. 

    Collective function.

    \param[in] function     the function to call
    \param[in] param        the parameter to the function
    \param[out] was_first   set to 1 if the current thread was the one to 
                            execute the function. unchanged if current thread
                            was not the one to execute the function; ignored
                            if NULL.
    \param[in] comm         The communicator.
    \returns the return value of the function
*/
void* tMPI_Once_wait(tMPI_Comm comm,void* (*function)(void*), void *param, 
                     int *was_first);

/** Allocate a shared block of memory as a collective call.
    
    Collective function.

    \param[in] comm         The communicator
    \param[in] size         The size in bytes to allocate 
*/
void* tMPI_Shmalloc(tMPI_Comm comm, size_t size);






#include "atomic.h" 

typedef struct
{
    tMPI_Atomic_t n_remaining; /* number of remaining operations */
    void *res;                 /* result pointer */
    tMPI_Comm comm;
} tMPI_Reduce_req;

tMPI_Reduce_req *tMPI_Reduce_req_alloc(tMPI_Comm comm);
#if 0
/** Execute fast a asynchronious reduce over comm. 

  Reduces array input with supplied funtion. This function may return before 
  the input array is ready to be written to again; to check for its completion,
  use the tMPI_Reduce_wait functions.

  \param[in] function   The function to reduce with (takes three arguments,
                        a size argument, two input pointers, and an output
                        pointer (which may be the same as one of the input
                        arguments, if res==NULL).
  \param[in] n          A size argument to pass to the function.
  \param[inout] input   The array of input data.
  \param[out] res       A temporary results array. May be NULL, in which case
                        the input data in 'in' gets overwritten.
*/
void tMPI_Reduce_async(tMPI_Reduce_req *req, 
                       void (*function)(int, void*, void*, void *), 
                       size_t n, void *input, void *res);


void tMPI_Reduce_wait(tMPI_Reduce_req *req);

void tMPI_Reduce_wait_results(tMPI_Reduce_req *req, void *res);

#endif


#ifdef __cplusplus
} /* closing extern "C" */
#endif

#endif /* _THREAD_MPI_COLLECTIVE_H_ */
