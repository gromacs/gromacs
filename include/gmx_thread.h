/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- 
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
 *  A key difference compared to pthreads is that you must use dynamic
 *  (instead of static) initialization of all variables by calling a routine.
 *  The reason for this is that static initialization doesn't work for many
 *  other thread libraries, not that we are lazy...
 *  The only type which can be initialized statically is gmx_thread_once_t.
 *  Note that you can only declare variable as a pointer to a mutex 
 *  and not a mutex variable directly. Declare a pointer, and call 
 *  gmx_thread_mutex_create() to get a pointer to a new mutex. Don't forget to
 *  embed it in a function called by gmx_thread_once() if your routine can be
 *  executed by multiple threads...
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
 *  Note that in contrast to pthreads you must always work with pointers
 *  to Gromacs threads. This is unfortunately necessary in order to hide the
 *  underlying implementation - we can't let this header file to depend
 *  on config.h since it might be used with a different compiler, or on 
 *  a different system from where Gromacs was compiled.
 */
typedef struct gmx_thread              gmx_thread_t;





/*! \brief Opaque mutex datatype for Gromacs threads
 *
 *  Note that in contrast to pthreads you must always work with pointers
 *  to Gromacs mutexes. This is unfortunately necessary in order to hide the
 *  underlying implementation - we can't let this header file to depend
 *  on config.h since it might be used with a different compiler, or on 
 *  a different system from where Gromacs was compiled.
 */
typedef struct gmx_thread_mutex        gmx_thread_mutex_t;



/*! \brief Opaque mutex datatype for Gromacs threads
 *
 *  Note that in contrast to pthreads you must always work with pointers
 *  to Gromacs threads. This is unfortunately necessary in order to hide the
 *  underlying implementation - we can't let this header file to depend
 *  on config.h since it might be used with a different compiler, or on 
 *  a different system from where Gromacs was compiled.
 */


/*! \brief One-time initialization data for Gromacs thread
 *
 *  This is an opaque datatype which is necessary for gmx_thread_once(),
 *  but since it needs to be initialized statically it must be defined
 *  in the header. You will be sorry if you touch the contents.
 *  Variables of this type should always be initialized statically to
 *  GMX_THREAD_ONCE_INITIALIZER.
 *
 *  This type is used as control data for single-time initialization.
 *  The most common example is a mutex at file scope used when calling 
 *  a non-threadsafe function, e.g. the FFTW initialization routines.
 *  In contrast to pthreads, many thread libraries (think windows) do not
 *  support static initialization of the mutex. This means it needs to be
 *  initialized by the first thread trying to use it, but we must also 
 *  make sure it is only initialized by exactly one thread.
 *
 *  The situation is solved by declaring the mutex statically in the file,
 *  but not initialized. However, you should also declare a static variable
 *  of the gmx_thread_once_t, and initialize it to GMX_THREAD_ONCE_INITIALIZER.
 *
 *  Write a small static routine (no arguments) that initializes your mutex,
 *  and in the routine where you need to use the mutex you first make a call to
 *  gmx_thread_once() with the once control data type and your initialization
 *  routines as arguments.
 */
typedef struct
 {
    int  started; /*!< -1 if nobody has touched it, >=0 when init. started */
    int  done;    /*!<  1 once the initialization is complete, 0 otherwise */
} gmx_thread_once_t;



/*! \brief Statical initializer for gmx_thread_once_t
 *
 *  See the description of the gmx_thread_once_t datatype for instructions
 *  on how to use this. Normally, all variables of that type should be 
 *  initialized statically to this value.
 */
#define GMX_THREAD_ONCE_INITIALIZER     { -1, 0 }




/*! \brief Thread-specific-data handle for Gromacs threads
 *
 *  This handle is used to keep track of thread-local storage.
 *  After it has been initialized with gmx_thread_key_create()
 *  you can use gmx_thread_setspecific() and gmx_thread_getspecific()
 *  to access a pointer which will be private to each thread.
 */
typedef struct gmx_thread_key          gmx_thread_key_t;





/*! \brief Create a new thread
 *
 *  The new thread will call start_routine() with the argument arg.
 *  Please be careful not to change arg after calling this function.
 * 
 *  \param start_routine   The function to call in the new thread
 *  \param arg             Argument to call with
 *  
 *  \return Pointer to a thread identifier, or NULL if an error occured.
 */
gmx_thread_t *
gmx_thread_create   (void *            (*start_routine)(void *),
                     void *              arg);




/*! \brief Wait for a specific thread to finish executing
 *
 *  If the thread has already finished the routine returns immediately.
 *
 *  \param thread      Thread to wait for.
 *  \param value_ptr   Pointer to location where to store pointer to exit value
 *                     from threads that called gmx_thread_exit().
 *  
 *  \return GMX_SUCCESS if the join went ok, or a non-zero error code.
 */
int
gmx_thread_join     (gmx_thread_t *    thread,
                     void **           value_ptr);



/*! \brief Create a new mutex
 *
 *  \return Pointer to the new mutex, NULL if creation failed.
 */
gmx_thread_mutex_t *
gmx_thread_mutex_create(void);




/*! \brief Kill a mutex you no longer need
*
*  Note that this call only frees resources allocated inside the mutex. It
*  does not free the gmx_thread_mutex_t memory area itself if you created it
*  with dynamic memory allocation.
* 
*  \param mtx  Pointer to a mutex variable to get rid of.
*  \return GMX_SUCCESS or a non-zero error code.
*/
int
gmx_thread_mutex_destroy(gmx_thread_mutex_t *mtx);




/*! \brief Wait for exclusive access to a mutex
*
*  This routine does not return until the mutex has been acquired.
*
*  \param mtx  Pointer to the mutex to lock
*  \return GMX_SUCCESS or a non-zero error code.
*/
int
gmx_thread_mutex_lock(gmx_thread_mutex_t *mtx);




/*! \brief Try to lock a mutex, return if busy
 *
 *  This routine always return directly. If the mutex was available and
 *  we successfully locked it we return GMX_SUCCESS, otherwise a non-zero
 *  error code (usually meaning the mutex was already locked).
 *
 *  \param mtx  Pointer to the mutex to try and lock
 *  \return GMX_SUCCESS or a non-zero error code.
 */
int
gmx_thread_mutex_trylock(gmx_thread_mutex_t *mtx);




/*! \brief Release the exclusive access to a mutex
 *
 *  \param mtx  Pointer to the mutex to release
 *  \return GMX_SUCCESS or a non-zero error code.
 */
int
gmx_thread_mutex_unlock(gmx_thread_mutex_t *mtx);




/*! \brief Initialize thread-specific-storage handle
 *
 *  The gmx_thread_key_t handle must always be initialized dynamically with
 *  this routine. If you need to initialize it statically in a file, use the
 *  gmx_thread_once() routine and corresponding data to initialize the 
 *  thread-specific-storage key the first time you access it.
 *
 *  \param  destructor   Routine to call (to free memory of key) when we quit
 *  \return Pointer to the new key, or NULL if an error occured.
 */
gmx_thread_key_t *     
gmx_thread_key_create(void                   (*destructor)(void *));




/*! \brief Delete thread-specific-storage handle
 *
 *  Calling this routine will kill the handle, and invoke the automatic 
 *  destructor routine for each non-NULL value pointed to by key.
 *
 *  \param  key  Handle to the key to get rid off.
 *  \return GMX_SUCCESS or a non-zero error message.
 */
int
gmx_thread_key_delete(gmx_thread_key_t *      key);




/*! \brief Get value for thread-specific-storage in this thread
 *
 *  If it has not yet been set, NULL is returned.
 *  
 *  \param key   Thread-specific-storage handle.
 *  \return Pointer-to-void, the value of the data in this thread.
*/
void *
gmx_thread_getspecific(gmx_thread_key_t *key);



/*! \brief Set value for thread-specific-storage in this thread
 *
 *  \param key     Thread-specific-storage handle.
 *  \param value   What to set the data to (pointer-to-void).
 *  \return GMX_SUCCESS or a non-zero error message.
 */
int
gmx_thread_setspecific(gmx_thread_key_t *key, void *value);



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
 *  \return GMX_SUCCESS or a non-zero error message.
 */
int
gmx_thread_once(gmx_thread_once_t *     once_data,
                void                    (*init_routine)(void));    




/*! \brief Terminate calling thread
 *
 *  Die voluntarily.
 *
 *  \param value_ptr   Pointer to a return value. Threads waiting for us to
 *                     join them can read this value if they try.
 */
int
gmx_thread_exit(void *      value_ptr);



/*! \brief Ask a thread to exit
 *
 *  This routine tries to end the execution of another thread, but there are
 *  no guarantees it will succeed.
 *
 *  \param thread     Handle to thread we want to see dead.
*  \return GMX_SUCCESS or a non-zero error message.
 */
int
gmx_thread_cancel(gmx_thread_t *    thread);


/*! \brief Lock a file so only one thread can use it
 *
 *  Call this routine before writing to logfiles or standard out, in order
 *  to avoid mixing output from multiple threads.
 */
int
gmx_lockfile(FILE *stream);


/*! \brief Unlock a file (allow other threads to use it)
 *
 *  Call this routine when you finish a write statement to a file, so other
 *  threads can use it again.
 */
int
gmx_unlockfile(FILE *stream);



#endif
