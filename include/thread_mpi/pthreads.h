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




/*! \file 
 *
 *  \brief POSIX threads specific data structures
 *  
 *  Most data structures used in thread_mpi/threads.h map directly onto 
 *  POSIX threads structures. This file contains typedefs for these data 
 *  types for use with POSIX threads.
 *
 *  \sa thread_mpi/threads.h for documentation.
 */

#ifndef DOXYGEN

/* we need this for all the data types */
#include <pthread.h>

/*! \brief Status for one-time initialization of thread stuff.
 *
 *  \internal
 *
 *  This is used both for the static initialization, and as
 *  a safeguard to catch errors where the user
 *  is sloppy and fails to initialize various things in the 
 *  thread system. It is only defined here to enable static
 *  initialization - don't use it outside the thread code.
 */
enum tMPI_Thread_once_status
{
    TMPI_THREAD_ONCE_STATUS_NOTCALLED = 0,   /*!< Not yet initialized     */
    TMPI_THREAD_ONCE_STATUS_PROGRESS = 1,    /*!< Somebody working on it  */
    TMPI_THREAD_ONCE_STATUS_READY = 2        /*!< Everything completed    */
}; 





/*! \brief Pthread implementation of the abstract tMPI_Thread type
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef pthread_t tMPI_Thread_t;


/*! \brief Opaque mutex datatype 
 *
 *  This type is only defined in the header to enable static
 *  initialization with TMPI_THREAD_MUTEX_INITIALIZER.
 *  You should _never_ touch the contents or create a variable
 *  with automatic storage class without calling tMPI_Thread_mutex_init().
 */
typedef pthread_mutex_t tMPI_Thread_mutex_t;
/*! \brief Static initializer for tMPI_Thread_mutex_t
 *
 *  See the description of the tMPI_Thread_mutex_t datatype for instructions
 *  on how to use this. Note that any variables initialized with this value
 *  MUST have static storage allocation.
 */
#define TMPI_THREAD_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER 


/*! \brief Pthread implementation of the abstract tMPI_Thread_key type 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef pthread_key_t tMPI_Thread_key_t;


/*! \brief One-time initialization data for thread
 *
 *  This is an opaque datatype which is necessary for tMPI_Thread_once(),
 *  but since it needs to be initialized statically it must be defined
 *  in the header. You will be sorry if you touch the contents.
 *  Variables of this type should always be initialized statically to
 *  TMPI_THREAD_ONCE_INIT.
 *
 *  This type is used as control data for single-time initialization.
 *  The most common example is a mutex at file scope used when calling 
 *  a non-threadsafe function, e.g. the FFTW initialization routines.
 *
 */
typedef pthread_once_t tMPI_Thread_once_t;
/*! \brief Static initializer for tMPI_Thread_once_t
 *
 *  See the description of the tMPI_Thread_once_t datatype for instructions
 *  on how to use this. Normally, all variables of that type should be 
 *  initialized statically to this value.
 */
#define TMPI_THREAD_ONCE_INIT PTHREAD_ONCE_INIT


/*! \brief Condition variable handle for threads
 *
 *  Condition variables are useful for synchronization together
 *  with a mutex: Lock the mutex and check if our thread is the last
 *  to the barrier. If no, wait for the condition to be signaled.
 *  If yes, reset whatever data you want and then signal the condition.
 *
 *  This should be considered an opaque structure, but since it is sometimes
 *  useful to initialize it statically it must go in the header. 
 *  You will be sorry if you touch the contents.
 *  
 *  There are two alternatives: Either initialize it as a static variable
 *  with TMPI_THREAD_COND_INITIALIZER, or call tMPI_Thread_cond_init()
 *  before using it.
 */
typedef pthread_cond_t tMPI_Thread_cond_t;


/*! \brief Static initializer for tMPI_Thread_cond_t
*
*  See the description of the tMPI_Thread_cond_t datatype for instructions
*  on how to use this. Note that any variables initialized with this value
*  MUST have static storage allocation.
*/
#define TMPI_THREAD_COND_INITIALIZER PTHREAD_COND_INITIALIZER

#ifdef TMPI_RWLOCK 
/*! \brief Read-write lock for threads

    Pthread implementation of the read-write lock (a lock that allows
    multiple readers, but only a single writer). 
  */
typedef pthread_rwlock_t tMPI_Thread_rwlock_t;


#if 0
/* this, for some reason, is not part of the POSIX standard, so we're not
   going to implement this in a library that needs to conform. */
/*! \brief static initializer for tMPI_Thread_rwlock_t

  See the description of the tMPI_Thread_cond_t datatype for instructions
  on how to use this. Note that any variables initialized with this value
  MUST have static storage allocation.
*/
#define TMPI_THREAD_RWLOCK_INITIALIZER PTHREAD_RWLOCK_INITIALIZER
#endif
#endif


/*! \brief Pthread implementation of barrier type. 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef struct tMPI_Thread_pthread_barrier
{
    pthread_mutex_t   mutex;     /*!< Lock for the barrier contents          */
    pthread_cond_t    cv;        /*!< Condition to signal barrier completion */
    int               threshold; /*!< Total number of members in barrier     */
    int               count;     /*!< Remaining count before completion      */
    int               cycle;     /*!< Alternating 0/1 to indicate round      */
}
tMPI_Thread_barrier_t;
/*! \brief Static initializer for tMPI_Thread_barrier_t
 *
 *  See the description of the tMPI_Thread_barrier_t datatype for instructions
 *  on how to use this. Note that variables initialized with this value
 *  MUST have static storage allocation.
 *
 * \param count  Threshold for barrier
 */
#define TMPI_THREAD_BARRIER_INITIALIZER(count)   {\
        PTHREAD_MUTEX_INITIALIZER,\
        PTHREAD_COND_INITIALIZER,\
        count, count, 0 \
        }


#endif



