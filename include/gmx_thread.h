/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
*
* $Id$
* 
* This file is part of Gromacs        Copyright (c) 1991-2004
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
#ifndef _GMX_THREAD_H_
#define _GMX_THREAD_H_

/*! \file gmx_thread.h
 *
 *  @brief Platform-independent multithreading support.
 *
 *  This file provides an portable thread interface very similar to POSIX 
 *  threads, but the actual implementation is sometimes done with inline atomic
 *  operations in assembly for performance. 
 *
 *  In other words, while the naming conventions are very similar to 
 *  pthreads, you should NOT assume that a Gromacs thread type
 *  (thread,mutex,key, etc) is the same as the Pthreads equivalent,
 *  even on platforms where we are using pthreads. For instance, on 
 *  all platforms where it is possible we try to use lower-latency spinlocks
 *  instead of real posix style mutexes.
 *
 *  IMPORTANT:
 *  The Gromacs thread implementation is used to guarantee threadsafe 
 *  operation for the gmx_message.h functions, which in turn is used by
 *  the gmx_memory.h allocation stuff.
 *  This means we cannot use gmx_message() or gmx_new() memory calls 
 *  in the implementation.
 *
 *  Since this module is merely intended to be a transparent wrapper around
 *  a system-dependent thread implementation, we simply echo errors to stderr.
 *  The user should check the return codes\] and take appropriate action 
 *  when using these functions (fat chance, but errors are rare).
 */



#include <stdio.h>



#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif




/*! \brief Opaque datatype for gromacs thread
 *
 *  A new thread will be initiated when you start a thread
 *  with gmx_thread_create(), and it will be destroyed when you
 *  join the thread with gmx_thread_join(). Remember to join your threads!
 */
typedef struct gmx_thread *            gmx_thread_t;





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
enum gmx_thread_once_status
{
    GMX_THREAD_ONCE_STATUS_NOTCALLED = 0,   /*!< Not yet initialized     */
    GMX_THREAD_ONCE_STATUS_PROGRESS = 1,    /*!< Somebody working on it  */
    GMX_THREAD_ONCE_STATUS_READY = 2        /*!< Everything completed    */
}; 




/*! \brief Opaque mutex datatype for Gromacs threads
 *
 *  This type is only defined in the header to enable static
 *  initialization with GMX_THREAD_MUTEX_INITIALIZER.
 *  You should _never_ touch the contents or create a variable
 *  with automatic storage class without calling gmx_thread_mutex_init().
 */
typedef struct gmx_mutex  
{
    enum gmx_thread_once_status       status;       /*!< Indicates completed init */
    void *                            actual_mutex; /*!< Implementation pointer   */
} 
gmx_thread_mutex_t;





/*! \brief Statical initializer for gmx_thread_mutex_t
 *
 *  See the description of the gmx_thread_mutex_t datatype for instructions
 *  on how to use this. Note that any variables initialized with this value
 *  MUST have static storage allocation.
 */
#define GMX_THREAD_MUTEX_INITIALIZER     { GMX_THREAD_ONCE_STATUS_NOTCALLED, NULL }





/*! \brief One-time initialization data for Gromacs thread
 *
 *  This is an opaque datatype which is necessary for gmx_thread_once(),
 *  but since it needs to be initialized statically it must be defined
 *  in the header. You will be sorry if you touch the contents.
 *  Variables of this type should always be initialized statically to
 *  GMX_THREAD_ONCE_INIT.
 *
 *  This type is used as control data for single-time initialization.
 *  The most common example is a mutex at file scope used when calling 
 *  a non-threadsafe function, e.g. the FFTW initialization routines.
 *
 */
typedef struct
{
    enum gmx_thread_once_status     status; /*!< not called, progress, or ready */
}
gmx_thread_once_t;



/*! \brief Statical initializer for gmx_thread_once_t
 *
 *  See the description of the gmx_thread_once_t datatype for instructions
 *  on how to use this. Normally, all variables of that type should be 
 *  initialized statically to this value.
 */
#define GMX_THREAD_ONCE_INIT       GMX_THREAD_ONCE_STATUS_NOTCALLED



/*! \brief Thread-specific-data handle for Gromacs threads
 *
 *  This handle is used to keep track of thread-local storage.
 *  After it has been initialized with gmx_thread_key_create()
 *  you can use gmx_thread_setspecific() and gmx_thread_getspecific()
 *  to access a pointer which will be private to each thread.
 */
typedef struct gmx_thread_key *         gmx_thread_key_t;



/*! \brief Condition variable handle for Gromacs threads
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
 *  with GMX_THREAD_COND_INITIALIZER, or call gmx_thread_cond_init()
 *  before using it.
 */
typedef struct gmx_thread_cond
{
    enum gmx_thread_once_status         status;      /*!< Initialized or not */
    void *                              actual_cond; /*!< Implementation ptr */
}
gmx_thread_cond_t;



/*! \brief Statical initializer for gmx_thread_cond_t
*
*  See the description of the gmx_thread_cond_t datatype for instructions
*  on how to use this. Note that any variables initialized with this value
*  MUST have static storage allocation.
*/
#define GMX_THREAD_COND_INITIALIZER     { GMX_THREAD_ONCE_STATUS_NOTCALLED, NULL }



/*! \brief Thread synchronization barrier datatype
 *
 *  Barrier variables are used to guarantee that all threads
 *  reach an execution point before they continue. This is the general
 *  Gromacs barrier type which uses a mutex and condition variable to
 *  suspend waiting threads. In performance-critical sections you might
 *  consider using the busy-waiting spinlock barrier instead.
 *
 *  Note that barrier variables must be initialized with
 *  gmx_thread_barrier_init().
 */
typedef struct gmx_thread_barrier 
{
    enum gmx_thread_once_status         status;         /*!< Initialized or not */
    void *                              actual_barrier; /*!< Implementation ptr */
    int                                 init_threshold; /*!< For static init    */
}
gmx_thread_barrier_t;



/*! \brief Statical initializer for gmx_thread_barrier_t
 *
 *  See the description of the gmx_thread_barrier_t datatype for instructions
 *  on how to use this. Note that variables initialized with this value
 *  MUST have static storage allocation.
 *
 * \param cnt  Threshold for barrier
 */
#define GMX_THREAD_BARRIER_INITIALIZER(cnt)   { GMX_THREAD_ONCE_STATUS_NOTCALLED, NULL, (cnt) }



/*! \brief Thread support status enumeration */
enum gmx_thread_support
{
    GMX_THREAD_SUPPORT_NO = 0,  /*!< Starting threads will fail */
    GMX_THREAD_SUPPORT_YES = 1  /*!< Thread support available   */
};



/*! \brief Check if threads are supported
 *
 *  This routine provides a cleaner way to check if threads are supported
 *  instead of sprinkling your code with preprocessor conditionals.
 *
 *  All Gromacs thread functions are still available even without thread support,
 *  but some of them might return failure error codes, for instance if you try
 *  to start a thread.
 * 
 *  \return 1 if threads are supported, 0 if not.
 */
enum gmx_thread_support
gmx_thread_support(void);



/*! \brief Create a new thread
 *
 *  The new thread will call start_routine() with the argument arg.
 *  Please be careful not to change arg after calling this function.
 * 
 *  \param thread          Pointer to opaque thread datatype
 *  \param start_routine   The function to call in the new thread
 *  \param arg             Argument to call with
 *  
 *  \return Status - 0 on success, or an error code.
 */
int
gmx_thread_create   (gmx_thread_t *      thread,
                     void *            (*start_routine)(void *),
                     void *              arg);




/*! \brief Wait for a specific thread to finish executing
 *
 *  If the thread has already finished the routine returns immediately.
 *
 *  \param thread      Opaque thread datatype to wait for.
 *  \param value_ptr   Pointer to location where to store pointer to exit value
 *                     from threads that called gmx_thread_exit().
 *  
 *  \return 0 if the join went ok, or a non-zero error code.
 */
int
gmx_thread_join     (gmx_thread_t      thread,
                     void **           value_ptr);



/*! \brief Initialize a new mutex
 *
 *  This routine must be called before using any mutex not initialized
 *  with static storage class and GMX_THREAD_MUTEX_INITIALIZER.
 *
 *  \param mtx   Pointer to a gromacs mutex opaque type.
 *  \return      0 or an error code.
 */
int
gmx_thread_mutex_init(gmx_thread_mutex_t *    mtx);




/*! \brief Kill a mutex you no longer need
*
*  Note that this call only frees resources allocated inside the mutex. It
*  does not free the gmx_thread_mutex_t memory area itself if you created it
*  with dynamic memory allocation.
* 
*  \param mtx  Pointer to a mutex variable to get rid of.
*  \return 0 or a non-zero error code.
*/
int
gmx_thread_mutex_destroy(gmx_thread_mutex_t *    mtx);




/*! \brief Wait for exclusive access to a mutex
*
*  This routine does not return until the mutex has been acquired.
*
*  \param mtx  Pointer to the mutex to lock
*  \return 0 or a non-zero error code.
*/
int
gmx_thread_mutex_lock(gmx_thread_mutex_t *    mtx);




/*! \brief Try to lock a mutex, return if busy
 *
 *  This routine always return directly. If the mutex was available and
 *  we successfully locked it we return 0, otherwise a non-zero
 *  error code (usually meaning the mutex was already locked).
 *
 *  \param mtx  Pointer to the mutex to try and lock
 *  \return 0 or a non-zero error code.
 */
int
gmx_thread_mutex_trylock(gmx_thread_mutex_t *   mtx);




/*! \brief Release the exclusive access to a mutex
 *
 *  \param mtx  Pointer to the mutex to release
 *  \return 0 or a non-zero error code.
 */
int
gmx_thread_mutex_unlock(gmx_thread_mutex_t *   mtx);




/*! \brief Initialize thread-specific-storage handle
 *
 *  The gmx_thread_key_t handle must always be initialized dynamically with
 *  this routine. If you need to initialize it statically in a file, use the
 *  gmx_thread_once() routine and corresponding data to initialize the 
 *  thread-specific-storage key the first time you access it.
 *
 *  \param  key          Pointer to opaque Gromacs thread key type.
 *  \param  destructor   Routine to call (to free memory of key) when we quit
 *
 *  \return status - 0 on sucess or a standard error code.
 *
 */
int
gmx_thread_key_create(gmx_thread_key_t *       key,
                      void                   (*destructor)(void *));




/*! \brief Delete thread-specific-storage handle
 *
 *  Calling this routine will kill the handle, and invoke the automatic 
 *  destructor routine for each non-NULL value pointed to by key.
 *
 *  \param  key  Opaque Gromacs key type to destroy
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_key_delete(gmx_thread_key_t       key);




/*! \brief Get value for thread-specific-storage in this thread
 *
 *  If it has not yet been set, NULL is returned.
 *  
 *  \param key   Thread-specific-storage handle.
 *  \return Pointer-to-void, the value of the data in this thread.
*/
void *
gmx_thread_getspecific(gmx_thread_key_t    key);



/*! \brief Set value for thread-specific-storage in this thread
 *
 *  \param key     Thread-specific-storage handle.
 *  \param value   What to set the data to (pointer-to-void).
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_setspecific(gmx_thread_key_t       key, 
                       void *                 value);



/*! \brief Run the provided routine exactly once
 *
 *  The control data must have been initialized before calling this routine,
 *  but you can do it with the static initialzer GMX_THREAD_ONCE_INIT.
 *
 *  gmx_thread_once() will not return to any of the calling routines until
 *  the initialization function has been completed.
 *
 *  \param once_data     Initialized one-time execution data
 *  \param init_routine  Function to call exactly once
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_once(gmx_thread_once_t *      once_data,
                void                    (*init_routine)(void));    



/*! \brief Initialize condition variable
 *
 *  This routine must be called before using any condition variable
 *  not initialized with static storage class and GMX_THREAD_COND_INITIALIZER.
 *
 *  \param cond  Pointer to previously allocated condition variable
 *  \return      0 or a non-zero error message.
 */
int
gmx_thread_cond_init(gmx_thread_cond_t *     cond);



/*! \brief Destroy condition variable
 *
 *  This routine should be called when you are done with a condition variable.
 *  Note that it only releases memory allocated internally, not the 
 *  gmx_thread_cond_t structure you provide a pointer to.
 *
 *  \param cond Pointer to gromacs condition variable.
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_cond_destroy(gmx_thread_cond_t *    cond);



/*! \brief Wait for a condition to be signaled
 *
 *  This routine releases the mutex, and waits for the
 *  condition to be signaled by another thread before
 *  it returns.
 *
 * \param cond  Gromacs condition variable
 *  \param mtx  Mutex protecting the condition variable
 *
 *  \return 0 or a non-zero error message.
 */
int 
gmx_thread_cond_wait(gmx_thread_cond_t *    cond,
                     gmx_thread_mutex_t *   mtx);




/*! \brief Unblock one waiting thread
 *
 *  This routine signals a condition variable to one
 *  thread (if any) waiting for it after calling
 *  gmx_thread_cond_wait().
 * 
 * \param cond  Gromacs condition variable
 *
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_cond_signal(gmx_thread_cond_t *  cond);


/*! \brief Unblock all waiting threads
*
*  This routine signals a condition variable to all
*  (if any) threads that are waiting for it after calling
*  gmx_thread_cond_wait().
* 
* \param cond  Gromacs condition variable
*
*  \return 0 or a non-zero error message.
*/
int
gmx_thread_cond_broadcast(gmx_thread_cond_t *  cond);




/*! \brief Terminate calling thread
 *
 *  Die voluntarily.
 *
 *  \param value_ptr   Pointer to a return value. Threads waiting for us to
 *                     join them can read this value if they try.
 *  \return 
 */
void
gmx_thread_exit(void *      value_ptr);



/*! \brief Ask a thread to exit
 *
 *  This routine tries to end the execution of another thread, but there are
 *  no guarantees it will succeed.
 *
 *  \param thread     Handle to thread we want to see dead.
 *  \return 0 or a non-zero error message.
 */
int
gmx_thread_cancel(gmx_thread_t      thread);


/*! \brief Initialize a synchronization barrier type
 *
 *  You only need to initialize a barrier once. They cycle 
 *  automatically, so after release it is immediately ready
 *  to accept waiting threads again.
 *
 *  \param barrier  Pointer to previously allocated barrier type
 *  \param count    Number of threads to synchronize. All threads
 *                  will be released after \a count calls to 
 *                  gmx_thread_barrier_wait(). 
 */
int
gmx_thread_barrier_init(gmx_thread_barrier_t *      barrier,
                        int                         count);



/*! \brief Release data in a barrier datatype
 *
 *  \param barrier  Pointer to previously 
 *                  initialized barrier.
 */
int
gmx_thread_barrier_destroy(gmx_thread_barrier_t *   barrier);


/*! \brief Perform barrier synchronization
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
int
gmx_thread_barrier_wait(gmx_thread_barrier_t *   barrier);



/*! \brief Lock a file so only one thread can use it
 *
 *  Call this routine before writing to logfiles or standard out, in order
 *  to avoid mixing output from multiple threads.
 */
void
gmx_lockfile(FILE *   stream);


/*! \brief Unlock a file (allow other threads to use it)
 *
 *  Call this routine when you finish a write statement to a file, so other
 *  threads can use it again.
 */
void
gmx_unlockfile(FILE *   stream);


#endif /* _GMX_THREAD_H_ */

