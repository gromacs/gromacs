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


#ifndef TMPI_THREADS_H_
#define TMPI_THREADS_H_

/*! \file threads.h
 *
 *  \brief Platform-independent multithreading support.
 *
 *  This file provides an portable thread interface very similar to POSIX
 *  threads, as a thin wrapper around the threads provided operating system
 *  (whether they be POSIX threads or something else).
 *
 *  In other words, while the naming conventions are very similar to
 *  pthreads, you should NOT assume that a thread_mpi thread type
 *  (thread,mutex,key, etc) is the same as the Pthreads equivalent,
 *  even on platforms where we are using pthreads.
 *
 *  Because the synchronization functions here are similar to the basic
 *  mutexes/conditions/barriers provided by the operating system,
 *  performance will most likely be worse than when using the atomic
 *  synchronization functions of atomic.h. On the other hand, because the
 *  operating system can schedule out waiting threads using these functions,
 *  they are the appropriate ones for I/O and initialization.
 *
 *  Since this module is merely intended to be a transparent wrapper around
 *  a system-dependent thread implementation, we simply echo errors to stderr.
 *  The user should check the return codes\] and take appropriate action
 *  when using these functions (fat chance, but errors are rare).
 */



#include <stdio.h>

#include "visibility.h"
#include "atomic.h"


#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif





/*! \brief Thread ID: abstract tMPI_Thread type
 *
 *  The contents of this structure depends on the actual threads
 *  implementation used.
 */
typedef struct tMPI_Thread* tMPI_Thread_t;



/*! \brief Opaque mutex datatype
 *
 *  This type is only defined in the header to enable static
 *  initialization with TMPI_THREAD_MUTEX_INITIALIZER.
 *  You should _never_ touch the contents or create a variable
 *  with automatic storage class without calling tMPI_Thread_mutex_init().
 */
typedef struct
{
    tMPI_Atomic_t      initialized; /*!< Whether \a mutex has been
                                       initialized. */
    struct tMPI_Mutex* mutex;       /*!< Actual mutex data structure. */
}  tMPI_Thread_mutex_t;
/*! \brief Static initializer for tMPI_Thread_mutex_t
 *
 *  See the description of the tMPI_Thread_mutex_t datatype for instructions
 *  on how to use this. Note that any variables initialized with this value
 *  MUST have static storage allocation.
 */
#define TMPI_THREAD_MUTEX_INITIALIZER { {0}, NULL }





/*! \brief Pthread implementation of the abstract tMPI_Thread_key type
 *
 *  The contents of this structure depends on the actual threads
 *  implementation used.
 */
typedef struct
{
    tMPI_Atomic_t           initialized; /*!< Whether \a key has been
                                            initialized. */
    struct tMPI_Thread_key *key;         /*!< Actual key data structure. */
} tMPI_Thread_key_t;





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
typedef struct
{
    tMPI_Atomic_t once; /*!< Whether the operation has been performed. */
} tMPI_Thread_once_t;
/*! \brief Static initializer for tMPI_Thread_once_t
 *
 *  See the description of the tMPI_Thread_once_t datatype for instructions
 *  on how to use this. Normally, all variables of that type should be
 *  initialized statically to this value.
 */
#define TMPI_THREAD_ONCE_INIT { {0} }




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
    tMPI_Atomic_t            initialized; /*!< Whether \a condp has been
                                             initialized. */
    struct tMPI_Thread_cond* condp;       /*!< Actual condition variable data
                                             structure. */
} tMPI_Thread_cond_t;
/*! \brief Static initializer for tMPI_Thread_cond_t
 *
 *  See the description of the tMPI_Thread_cond_t datatype for instructions
 *  on how to use this. Note that any variables initialized with this value
 *  MUST have static storage allocation.
 */
#define TMPI_THREAD_COND_INITIALIZER { {0}, NULL}






/*! \brief Pthread implementation of barrier type.
 *
 *  The contents of this structure depends on the actual threads
 *  implementation used.
 */
typedef struct
{
    tMPI_Atomic_t               initialized; /*!< Whether \a barrierp has been initialized. */
    struct tMPI_Thread_barrier* barrierp;    /*!< Actual barrier data structure. */
    volatile int                threshold;   /*!< Total number of members in barrier     */
    volatile int                count;       /*!< Remaining count before completion      */
    volatile int                cycle;       /*!< Alternating 0/1 to indicate round      */
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
        NULL, count, count, 0 \
}






/** Thread support status enumeration */
enum tMPI_Thread_support
{
    TMPI_THREAD_SUPPORT_NO  = 0, /*!< Starting threads will fail */
    TMPI_THREAD_SUPPORT_YES = 1  /*!< Thread support available   */
};


/** Thread setaffinity support status enumeration */
enum tMPI_Thread_setaffinity_support
{
    TMPI_SETAFFINITY_SUPPORT_NO  = 0, /*!< Setting thread affinity not
                                           supported */
    TMPI_SETAFFINITY_SUPPORT_YES = 1  /*!< Setting thread affinity supported */
};




/** handle a fatal error.

   \param file     source code file name of error.
   \param line     source code line number of error.
   \param message  format string for error message.
 */
TMPI_EXPORT
void tMPI_Fatal_error(const char *file, int line, const char *message, ...);
/** Convenience macro for the first two arguments to tMPI_Fatal_error(). */
#define TMPI_FARGS __FILE__, __LINE__



/*! \name Thread creation, destruction, and inspection
 \{ */
/** Check if threads are supported
 *
 *  This routine provides a cleaner way to check if threads are supported
 *  instead of sprinkling your code with preprocessor conditionals.
 *
 *  All thread functions are still available even without thread support,
 *  but some of them might return failure error codes, for instance if you try
 *  to start a thread.
 *
 *  \return 1 if threads are supported, 0 if not.
 */
TMPI_EXPORT
enum tMPI_Thread_support tMPI_Thread_support(void);


/** Get the number of hardware threads that can be run simultaneously.

    Returns the total number of cores and SMT threads that can run.

    \returns The maximum number of threads that can run simulataneously.
        If this number cannot be determined for the current architecture,
        0 is returned.
 */
TMPI_EXPORT
int tMPI_Thread_get_hw_number(void);


/** Create a new thread
 *
 *  The new thread will call start_routine() with the argument arg.
 *
 *  Please be careful not to change arg after calling this function.
 *
 *  \param[out] thread          Pointer to thread ID
 *  \param[in] start_routine    The function to call in the new thread
 *  \param[in] arg              Argument to call with
 *
 *  \return Status - 0 on success, or an error code.
 */
TMPI_EXPORT
int tMPI_Thread_create(tMPI_Thread_t *thread,
                       void           * (*start_routine)(void *),
                       void         * arg);




/** Wait for a specific thread to finish executing
 *
 *  If the thread has already finished the routine returns immediately.
 *
 *  \param[in] thread       Pointer to thread ID
 *  \param[out] value_ptr   Pointer to location where to store pointer to
 *                          exit value from threads that called
 *                          tMPI_Thread_exit().
 *
 *  \return 0 if the join went ok, or a non-zero error code.
 */
TMPI_EXPORT
int tMPI_Thread_join(tMPI_Thread_t thread, void **value_ptr);


/** Terminate calling thread
 *
 *  Die voluntarily.
 *
 *  \param value_ptr   Pointer to a return value. Threads waiting for us to
 *                     join them can read this value if they try.
 *  \return
 */
TMPI_EXPORT
void tMPI_Thread_exit(void *value_ptr);



/** Ask a thread to exit
 *
 *  This routine tries to end the execution of another thread, but there are
 *  no guarantees it will succeed.
 *
 *  \param thread     Handle to thread we want to see dead.
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cancel(tMPI_Thread_t thread);




/** Get a thread ID of the calling thread.
 *
 * This function also works on threads not started with tMPI_Thread_create,
 * or any other function in thread_mpi. This makes it possible to, for
 * example assign thread affinities to any thread.
 *
 * \return A thread ID of the calling thread */
TMPI_EXPORT
tMPI_Thread_t tMPI_Thread_self(void);



/** Check whether two thread pointers point to the same thread
 *
 * \param[in]  t1  Thread ID 1
 * \param[in]  t2  Thread ID 2
 * \return non-zero if the thread structs refer to the same thread,
            0 if the threads are different*/
TMPI_EXPORT
int tMPI_Thread_equal(tMPI_Thread_t t1, tMPI_Thread_t t2);


/** Check whether this platform supports setting of thread affinity
 *
 * This function returns TMPI_SETAFFINITY_SUPPORT_YES if setting thread
 * affinity is supported by the platform, and TMPI_SETAFFINITY_SUPPORT_NO
 * if not. If this function returns 0, the function
 * tMPI_Thread_setaffinity_single will simply return 0 itself, effectively
 * ignoring the request.
 *
 * \return TMPI_SETAFFINITY_SUPPORT_YES if setting affinity is supported,
            TMPI_SETAFFINITY_SUPPORT_NO otherwise */
TMPI_EXPORT
enum tMPI_Thread_setaffinity_support tMPI_Thread_setaffinity_support(void);


/** Set thread affinity to a single core
 *
 * This function sets the thread affinity of a thread to a a specific
 * numbered processor. This only works if the underlying operating system
 * supports it. The processor number must be between 0 and the number returned
 * by tMPI_Thread_get_hw_number().
 *
 * \param[in] thread   Thread ID of the thread to set affinity for
 * \param[in] nr       Processor number to set affinity to
 * \return zero on success, non-zero on error */
TMPI_EXPORT
int tMPI_Thread_setaffinity_single(tMPI_Thread_t thread, unsigned int nr);


/*! \} */
/*! \name Mutexes
 \{ */


/** Initialize a new mutex
 *
 *  This routine must be called before using any mutex not initialized
 *  with static storage class and TMPI_THREAD_MUTEX_INITIALIZER.
 *
 *  \param mtx   Pointer to a mutex opaque type.
 *  \return      0 or an error code.
 */
TMPI_EXPORT
int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx);




/** Kill a mutex you no longer need
 *
 *  Note that this call only frees resources allocated inside the mutex. It
 *  does not free the tMPI_Thread_mutex_t memory area itself if you created it
 *  with dynamic memory allocation.
 *
 *  \param mtx  Pointer to a mutex variable to get rid of.
 *  \return 0 or a non-zero error code.
 */
TMPI_EXPORT
int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx);




/** Wait for exclusive access to a mutex
 *
 *  This routine does not return until the mutex has been acquired.
 *
 *  \param mtx  Pointer to the mutex to lock
 *  \return 0 or a non-zero error code.
 */
TMPI_EXPORT
int tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx);




/** Try to lock a mutex, return if busy
 *
 *  This routine always return directly. If the mutex was available and
 *  we successfully locked it we return 0, otherwise a non-zero
 *  return code (usually meaning the mutex was already locked).
 *
 *  \param mtx  Pointer to the mutex to try and lock
 *  \return 0 if locked, non-zero if not locked or an error occurred.
 */
TMPI_EXPORT
int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx);




/** Release the exclusive access to a mutex
 *
 *  \param mtx  Pointer to the mutex to release
 *  \return 0 or a non-zero error code.
 */
TMPI_EXPORT
int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx);



/*! \} */
/*! \name Thread-specific storage
 \{ */


/** Initialize thread-specific-storage handle
 *
 *  The tMPI_Thread_key_t handle must always be initialized dynamically with
 *  this routine. If you need to initialize it statically in a file, use the
 *  tMPI_Thread_once() routine and corresponding data to initialize the
 *  thread-specific-storage key the first time you access it.
 *
 *  \param  key          Pointer to opaque thread key type.
 *  \param  destructor   Routine to call (to free memory of key) when we quit
 *
 *  \return status - 0 on sucess or a standard error code.
 *
 */
TMPI_EXPORT
int tMPI_Thread_key_create(tMPI_Thread_key_t *key, void (*destructor)(void *));




/** Delete thread-specific-storage handle
 *
 *  Calling this routine will kill the handle, and invoke the automatic
 *  destructor routine for each non-NULL value pointed to by key.
 *
 *  \param  key  Opaque key type to destroy
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_key_delete(tMPI_Thread_key_t key);




/** Get value for thread-specific-storage in this thread
 *
 *  If it has not yet been set, NULL is returned.
 *
 *  \param key   Thread-specific-storage handle.
 *  \return Pointer-to-void, the value of the data in this thread.
 */
TMPI_EXPORT
void * tMPI_Thread_getspecific(tMPI_Thread_key_t key);



/** Set value for thread-specific-storage in this thread
 *
 *  \param key     Thread-specific-storage handle.
 *  \param value   What to set the data to (pointer-to-void).
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_setspecific(tMPI_Thread_key_t key, void *value);


/*! \} */
/*! \name Run-once
 \{ */


/** Run the provided routine exactly once
 *
 *  The control data must have been initialized before calling this routine,
 *  but you can do it with the static initialzer TMPI_THREAD_ONCE_INIT.
 *
 *  tMPI_Thread_once() will not return to any of the calling routines until
 *  the initialization function has been completed.
 *
 *  \param once_data     Initialized one-time execution data
 *  \param init_routine  Function to call exactly once
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_once(tMPI_Thread_once_t *once_data,
                     void                (*init_routine)(void));


/*! \} */
/*! \name Condition variables
 \{ */

/** Initialize condition variable
 *
 *  This routine must be called before using any condition variable
 *  not initialized with static storage class and TMPI_THREAD_COND_INITIALIZER.
 *
 *  \param cond  Pointer to previously allocated condition variable
 *  \return      0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond);



/** Destroy condition variable
 *
 *  This routine should be called when you are done with a condition variable.
 *  Note that it only releases memory allocated internally, not the
 *  tMPI_Thread_cond_t structure you provide a pointer to.
 *
 *  \param cond Pointer to condition variable.
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond);



/** Wait for a condition to be signaled
 *
 *  This routine releases the mutex, and waits for the condition to be
 *  signaled by another thread before it returns.
 *
 *  Note that threads waiting for conditions with tMPI_Thread_cond_wait
 *  may be subject to spurious wakeups: use this function in a while loop
 *  and check the state of a predicate associated with the wakeup
 *  before leaving the loop.
 *
 * \param cond  condition variable
 *  \param mtx  Mutex protecting the condition variable
 *
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cond_wait(tMPI_Thread_cond_t  *cond,
                          tMPI_Thread_mutex_t *mtx);




/** Unblock one waiting thread
 *
 *  This routine signals a condition variable to one
 *  thread (if any) waiting for it after calling
 *  tMPI_Thread_cond_wait().
 *
 * \param cond  condition variable
 *
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cond_signal(tMPI_Thread_cond_t *cond);


/** Unblock all waiting threads
 *
 *  This routine signals a condition variable to all
 *  (if any) threads that are waiting for it after calling
 *  tMPI_Thread_cond_wait().
 *
 * \param cond  condition variable
 *
 *  \return 0 or a non-zero error message.
 */
TMPI_EXPORT
int tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *cond);



/*! \} */
/*! \name Barriers
 \{ */


/** Initialize a synchronization barrier type
 *
 *  You only need to initialize a barrier once. They cycle
 *  automatically, so after release it is immediately ready
 *  to accept waiting threads again.
 *
 *  \param barrier  Pointer to previously allocated barrier type
 *  \param count    Number of threads to synchronize. All threads
 *                  will be released after \a count calls to
 *                  tMPI_Thread_barrier_wait().
 */
TMPI_EXPORT
int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier, int count);



/** Release data in a barrier datatype
 *
 *  \param barrier  Pointer to previously
 *                  initialized barrier.
 */
TMPI_EXPORT
int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier);


/** Perform barrier synchronization
 *
 *  This routine blocks until it has been called N times,
 *  where N is the count value the barrier was initialized with.
 *  After N total calls all threads return. The barrier automatically
 *  cycles, and thus requires another N calls to unblock another time.
 *
 *  \param barrier  Pointer to previously create barrier.
 *
 *  \return The last thread returns -1, all the others 0.
 */
TMPI_EXPORT
int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *barrier);

/*! \} */


#ifdef __cplusplus
}
#endif

#endif /* TMPI_THREADS_H_ */
