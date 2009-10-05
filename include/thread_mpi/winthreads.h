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



/*! \file winthreads.h
 *
 *  @brief Windows threads specific data structures 
 *
 *  For use of Windows threads in thread_mpi/threads.h
 * 
 *  \sa thread_mpi/threads.h for documentation.
 */


#ifndef DOXYGEN

/* we need this for all the data types. We use WIN32_LEAN_AND_MEAN to avoid 
   polluting the global namespace. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN


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





/*! \brief win32 implementation of the abstract tMPI_Thread type
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef HANDLE tMPI_Thread_t;


/*! \brief Opaque mutex datatype 
 *
 *  This type is only defined in the header to enable static
 *  initialization with TMPI_THREAD_MUTEX_INITIALIZER.
 *  You should _never_ touch the contents or create a variable
 *  with automatic storage class without calling tMPI_Thread_mutex_init().
 */
typedef struct 
{ 
    volatile enum tMPI_Thread_once_status init_state; 
    CRITICAL_SECTION cs;
} tMPI_Thread_mutex_t;
/*! \brief Static initializer for tMPI_Thread_mutex_t
 *
 *  See the description of the tMPI_Thread_mutex_t datatype for instructions
 *  on how to use this. Note that any variables initialized with this value
 *  MUST have static storage allocation.
 */
#define TMPI_THREAD_MUTEX_INITIALIZER { TMPI_THREAD_ONCE_STATUS_NOTCALLED } 


/*! \brief Pthread implementation of the abstract tMPI_Thread_key type 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.  */
typedef DWORD tMPI_Thread_key_t;


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
#if 0
/* use this once Vista is the oldest supported Windows version: */
typedef INIT_ONCE tMPI_Thread_once_t;
#else
typedef volatile int tMPI_Thread_once_t;
#endif
/*! \brief Static initializer for tMPI_Thread_once_t
 *
 *  See the description of the tMPI_Thread_once_t datatype for instructions
 *  on how to use this. Normally, all variables of that type should be 
 *  initialized statically to this value.
 */
#define TMPI_THREAD_ONCE_INIT 0


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
typedef struct 
{ 
    volatile enum tMPI_Thread_once_status init_state; 


#if 0
    /* this works since Windows Vista: */
    CONDITION_VARIABLE cv;
#else
    /* this data structure and its algorithms are based on 
       'Strategies for Implementing POSIX Condition Variables on Win32'
       by
       Douglas C. Schmidt and Irfan Pyarali
       Department of Computer Science
       Washington University, St. Louis, Missouri
       http://www.cs.wustl.edu/~schmidt/win32-cv-1.html */
    int Nwaiters; /* number of waiting threads */
    CRITICAL_SECTION wtr_lock; /* lock for Nwaiters */
    int Nrelease; /* number of threads to release in broadcast/signal */
    int cycle; /* cycle number so threads can't steal signals */
    HANDLE ev; /* the event used to trigger WaitForSingleObject.  
                  Is a manual reset event.  */
#endif
} tMPI_Thread_cond_t;

/*typedef pthread_cond_t tMPI_Thread_cond_t;*/


/*! \brief Static initializer for tMPI_Thread_cond_t
*
*  See the description of the tMPI_Thread_cond_t datatype for instructions
*  on how to use this. Note that any variables initialized with this value
*  MUST have static storage allocation.
*/
#define TMPI_THREAD_COND_INITIALIZER { TMPI_THREAD_ONCE_STATUS_NOTCALLED }



#ifdef TMPI_RWLOCK
/*! \brief Read-write lock for threads

  Windows implementation of the read-write lock (a lock that allows
  multiple readers, but only a single writer). 
*/
typedef SRWLOCK tMPI_Thread_rwlock_t;
#endif

/*! \brief Pthread implementation of barrier type. 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef struct tMPI_Thread_pthread_barrier
{
    volatile enum tMPI_Thread_once_status init_state; 
    int                            init_threshold; /* < N */
    int                threshold; /*!< Total number of members in barrier     */
    int                count;     /*!< Remaining count before completion      */
    int                cycle;     /*!< Alternating 0/1 to indicate round      */

#if 0
    /* use this once Vista is the oldest supported windows version: */
     CRITICAL_SECTION   cs;        /*!< Lock for the barrier contents          */
     CONDITION_VARIABLE cv;        /*!< Condition to signal barrier completion */
#else
    tMPI_Thread_mutex_t cs;       /*!< Lock for the barrier contents          */
    tMPI_Thread_cond_t cv;        /*!< Condition to signal barrier completion */
#endif
}tMPI_Thread_barrier_t;

/*! \brief Static initializer for tMPI_Thread_barrier_t
 *
 *  See the description of the tMPI_Thread_barrier_t datatype for instructions
 *  on how to use this. Note that variables initialized with this value
 *  MUST have static storage allocation.
 *
 * \param count  Threshold for barrier
 */
#define TMPI_THREAD_BARRIER_INITIALIZER(count)   { \
                                    TMPI_THREAD_ONCE_STATUS_NOTCALLED, \
                                    count, count, 0 }



#endif

