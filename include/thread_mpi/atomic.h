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

#ifndef _TMPI_ATOMIC_H_
#define _TMPI_ATOMIC_H_

/*! \file atomic.h
 *
 *  \brief Atomic operations for fast SMP synchronization
 *
 *  This file defines atomic integer operations and spinlocks for 
 *  fast synchronization in performance-critical regions.
 *
 *  In general, the best option is to use functions without explicit 
 *  locking, e.g. tMPI_Atomic_fetch_add() or tMPI_Atomic_cmpxchg().
 *  
 *  Depending on the architecture/compiler, these operations may either
 *  be provided as functions or macros; be aware that those macros may
 *  reference their arguments repeatedly, possibly leading to multiply 
 *  evaluated code with side effects: be careful with what you use as 
 *  arguments.
 *
 *  Not all architectures support atomic operations though inline assembly,
 *  and even if they do it might not be implemented here. In that case
 *  we use a fallback mutex implementation, so you can always count on
 *  the function interfaces working.
 *
 *  Don't use spinlocks in non-performance-critical regions like file I/O.
 *  Since they always spin busy they would waste CPU cycles instead of 
 *  properly yielding to a computation thread while waiting for the disk.
 *
 *  Finally, note that all our spinlock operations are defined to return
 *  0 if initialization or locking completes successfully.
 *  This is the opposite of some other implementations, but the same standard
 *  as used for pthread mutexes. So, if e.g. are trying to lock a spinlock,
 *  you will have gotten the lock if the return value is 0.
 * 
 *  tMPI_Spinlock_islocked(x) obviously still returns 1 if the lock is locked,
 *  and 0 if it is available, though...
 */
/* Se the comments on the non-atomic versions for explanations */

#include <stdio.h>



#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif

#ifdef __GNUC__
#define TMPI_GCC_VERSION (__GNUC__ * 10000 \
                          + __GNUC_MINOR__ * 100 \
                          + __GNUC_PATCHLEVEL__)
#endif



/* first check for gcc/icc platforms */
#if ( defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__PATHSCALE__) )

/* now check specifically for several architectures: */
#if (defined(i386) || defined(__x86_64__)) 
/* first x86: */
#include "atomic/gcc_x86.h"
/*#include "atomic/gcc.h"*/

#elif (defined(__ia64__))
/* then ia64: */
#include "atomic/gcc_ia64.h"

#elif (defined(__powerpc__) || (defined(__ppc__)) )
/* and powerpc: */
#include "atomic/gcc_ppc.h"

#elif
/* otherwise, there's a generic gcc intrinsics version: */
#include "atomic/gcc.h"
#endif /* end of check for gcc specific architectures */

/* not gcc: */
#elif (defined(_MSC_VER) && (_MSC_VER >= 1200))
/* Microsoft Visual C on x86, define taken from FFTW who got it from 
   Morten Nissov */
#include "atomic/msvc.h"

#elif ( (defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM))  && \
        (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using xlC inline assembly. 
 * Recent versions of xlC (>=7.0) _partially_ support GCC inline assembly
 * if you use the option -qasm=gcc but we have had to hack things a bit, in 
 * particular when it comes to clobbered variables. Since this implementation
 * _could_ be buggy, we have separated it from the known-to-be-working gcc
 * one above.
 */
#include "atomic/xlc_ppc.h"

#elif defined(__xlC__) && defined (_AIX)
/* IBM xlC compiler on AIX */
#include "atomic/xlc_aix.h"

#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
/* HP compiler on ia64 */
#include "atomic/hpux.h"




#else
/* No atomic operations, use mutex fallback. Documentation is in x86 section */

#ifdef TMPI_CHECK_ATOMICS
#error No atomic operations implemented for this cpu/compiler combination. 
#endif


/** Memory barrier operation

 Modern CPUs rely heavily on out-of-order execution, and one common feature
 is that load/stores might be reordered. Also, when using inline assembly
 the compiler might already have loaded the variable we are changing into
 a register, so any update to memory won't be visible.
 
 This command creates a memory barrier, i.e. all memory results before
 it in the code should be visible to all memory operations after it - the
 CPU cannot propagate load/stores across it.

 \hideinitializer
 */
#define tMPI_Atomic_memory_barrier()

/* System mutex used for locking to guarantee atomicity */
static tMPI_Thread_mutex_t tMPI_Atomic_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

/** Atomic operations datatype
 *
 *  Portable synchronization primitives like mutexes are effective for
 *  many purposes, but usually not very high performance.
 *  One of the problem is that you have the overhead of a function call,
 *  and another is that Mutexes often have extra overhead to make the
 *  scheduling fair. Finally, if performance is important we don't want
 *  to suspend the thread if we cannot lock a mutex, but spin-lock at 100%
 *  CPU usage until the resources is available (e.g. increment a counter).
 *
 *  These things can often be implemented with inline-assembly or other
 *  system-dependent functions, and we provide such functionality for the
 *  most common platforms. For portability we also have a fallback 
 *  implementation using a mutex for locking.
 *
 *  Performance-wise, the fastest solution is always to avoid locking 
 *  completely (obvious, but remember it!). If you cannot do that, the
 *  next best thing is to use atomic operations that e.g. increment a
 *  counter without explicit locking. Spinlocks are useful to lock an
 *  entire region, but leads to more overhead and can be difficult to
 *  debug - it is up to you to make sure that only the thread owning the
 *  lock unlocks it!
 *
 *  You should normally NOT use atomic operations for things like 
 *  I/O threads. These should yield to other threads while waiting for 
 *  the disk instead of spinning at 100% CPU usage.
 *
 *  It is imperative that you use the provided routines for reading
 *  and writing, since some implementations require memory barriers before
 *  the CPU or memory sees an updated result. The structure contents is
 *  only visible here so it can be inlined for performance - it might
 *  change without further notice.
 *
 *  \note No initialization is required for atomic variables.
 *
 *  Currently, we have (real) atomic operations for:
 *
 *  - gcc version 4.1 and later (all platforms)
 *  - x86 or x86_64, using GNU compilers
 *  - x86 or x86_64, using Intel compilers 
 *  - x86 or x86_64, using Pathscale compilers
 *  - Itanium, using GNU compilers 
 *  - Itanium, using Intel compilers
 *  - Itanium, using HP compilers
 *  - PowerPC, using GNU compilers 
 *  - PowerPC, using IBM AIX compilers 
 *  - PowerPC, using IBM compilers >=7.0 under Linux or Mac OS X.
 *
 * \see
 * - tMPI_Atomic_get
 * - tMPI_Atomic_set
 * - tMPI_Atomic_cmpxchg
 * - tMPI_Atomic_add_return
 * - tMPI_Atomic_fetch_add
 */
typedef struct tMPI_Atomic
{
        int value;  
}
tMPI_Atomic_t;


/** Atomic pointer type equivalent to tMPI_Atomic_t
 *
 * Useful for lock-free and wait-free data structures.
 * The only operations available for this type are:
 * \see
 * - tMPI_Atomic_ptr_get
 * - tMPI_Atomic_ptr_set
 * - tMPI_Atomic_ptr_cmpxchg
*/
typedef struct tMPI_Atomic_ptr
{
        void* value;  
}
tMPI_Atomic_ptr_t;


/** Spinlock
 *
 *  Spinlocks provide a faster synchronization than mutexes,
 *  although they consume CPU-cycles while waiting. They are implemented
 *  with atomic operations and inline assembly whenever possible, and
 *  otherwise we use a fallback implementation where a spinlock is identical
 *  to a mutex (this is one of the reasons why you have to initialize them).
 *
 *  There are no guarantees whatsoever about fair scheduling or
 *  debugging if you make a mistake and unlock a variable somebody
 *  else has locked - performance is the primary goal of spinlocks.
 *
 * \see
 * - tMPI_Spinlock_init
 * - tMPI_Spinlock_lock
 * - tMPI_Spinlock_unlock
 * - tMPI_Spinlock_trylock
 * - tMPI_Spinlock_wait
 */
typedef struct 
{
#ifndef DOXYGEN
    tMPI_Thread_mutex_t lock; /* we don't want this documented */
#endif
} tMPI_Spinlock_t;
/*#define tMPI_Spinlock_t     tMPI_Thread_mutex_t*/

 /*! \def TMPI_SPINLOCK_INITIALIZER
 * \brief Spinlock static initializer
 *
 *  This is used for static spinlock initialization, and has the same
 *  properties as TMPI_THREAD_MUTEX_INITIALIZER has for mutexes.
 *  This is only for inlining in the tMPI_Thread.h header file. Whether
 *  it is 0, 1, or something else when unlocked depends on the platform.
 *  Don't assume anything about it. It might even be a mutex when using the
 *  fallback implementation!
 *
 *  \hideinitializer
 */
#  define TMPI_SPINLOCK_INITIALIZER   { TMPI_THREAD_MUTEX_INITIALIZER }

/* Since mutexes guarantee memory barriers this works fine */
/** Return value of an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable to read
 *  \return     Integer value of the atomic variable
 *
 *  \hideinitializer
 */
#ifdef DOXYGEN
static inline int tMPI_Atomic_get( tMPI_Atomic_t &a);
#else
#define tMPI_Atomic_get(a)   ((a)->value)
#endif

/** Write value to an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable
 *  \param  i   Integer to set the atomic variable to.
 *
 *  \hideinitializer
 */
static inline void tMPI_Atomic_set(tMPI_Atomic_t *   a, 
                                  int              i)
{
    /* Mutexes here are necessary to guarantee memory visibility */
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = i;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}


/** Return value of an atomic pointer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable to read
 *  \return     Pointer value of the atomic variable
 *
 *  \hideinitializer
 */
#ifdef DOXYGEN
static inline void* tMPI_Atomic_ptr_get( tMPI_Atomic_ptr_t &a);
#else
#define tMPI_Atomic_ptr_get(a)   ((a)->value)
#endif




/** Write value to an atomic pointer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable
 *  \param  p   Pointer value to set the atomic variable to.
 *
 *  \hideinitializer
 */
static inline void tMPI_Atomic_ptr_set(tMPI_Atomic_t *   a, 
                                      void*            p)
{
    /* Mutexes here are necessary to guarantee memory visibility */
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = (void*)p;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}


/** Add integer to atomic variable
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param a   atomic datatype to modify
 *  \param i   integer to increment with. Use i<0 to subtract atomically.
 *
 *  \return The new value (after summation).
 */
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *   a, 
                                        int              i)
{
    int t;
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    t = a->value + i;
    a->value = t;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}



/** Add to variable, return the old value.
 *
 *  This operation is quite useful for synchronization counters.
 *  By performing a fetchadd with N, a thread can e.g. reserve a chunk 
 *  with the next N iterations, and the return value is the index
 *  of the first element to treat.
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param a   atomic datatype to modify
 *  \param i   integer to increment with. Use i<0 to subtract atomically.
 *
 *  \return    The value of the atomic variable before addition.
 */
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *   a,
                                       int              i)
{
    int old_value;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    old_value  = a->value;
    a->value   = old_value + i;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return old_value;
}



/** Atomic compare-exchange operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is updated to the new value, 
 *   and otherwise left unchanged. 
 *  
 *   This is a very useful synchronization primitive: You can start by reading
 *   a value (without locking anything), perform some calculations, and then
 *   atomically try to update it in memory unless it has changed. If it has
 *   changed you will get an error return code - reread the new value
 *   an repeat the calculations in that case.
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param old_val  Integer value read from the atomic type at an earlier point
 *   \param new_val  New value to write to the atomic type if it currently is
 *                   identical to the old value.
 *
 *   \return The value of the atomic memory variable in memory when this 
 *           instruction was executed. This, if the operation succeeded the
 *           return value was identical to the \a old parameter, and if not
 *           it returns the updated value in memory so you can repeat your
 *           operations on it. 
 *
 *   \note   The exchange occured if the return value is identical to \a old.
 */
static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *           a, 
                                     int                      old_val,
                                     int                      new_val)
{
    int t;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    t=old_val;
    if (a->value == old_val)
    {
        a->value = new_val;
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}




/** Atomic pointer compare-exchange operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is updated to the new value, 
 *   and otherwise left unchanged. 
 *  
 *   This is essential for implementing wait-free lists and other data
 *   structures. 
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param old_val  Pointer value read from the atomic type at an earlier point
 *   \param new_val  New value to write to the atomic type if it currently is
 *                   identical to the old value.
 *
 *   \return The value of the atomic pointer in memory when this 
 *           instruction was executed. This, if the operation succeeded the
 *           return value was identical to the \a old parameter, and if not
 *           it returns the updated value in memory so you can repeat your
 *           operations on it. 
 *
 *   \note   The exchange occured if the return value is identical to \a old.
 */
static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t * a, 
                                           void*              old_val,
                                           void*              new_val)
{
    void *t;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    t=old_val;
    if (a->value == old_val)
    {
        a->value = new_val;
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}


/** Initialize spinlock
 *
 *  In theory you can call this from multiple threads, but remember
 *  that we don't check for errors. If the first thread proceeded to
 *  lock the spinlock after initialization, the second will happily
 *  overwrite the contents and unlock it without warning you.
 *
 *  \param x      Spinlock pointer.
 *
 *  \hideinitializer
 */
#ifdef DOXYGEN
void tMPI_Spinlock_init( tMPI_Spinlock_t &x);
#else
#define tMPI_Spinlock_init(x)       tMPI_Thread_mutex_init((x)->lock)
#endif

/** Acquire spinlock
 *
 *  This routine blocks until the spinlock is available, and
 *  the locks it again before returning.
 *
 *  \param x     Spinlock pointer
 */
#ifdef DOXYGEN
void tMPI_Spinlock_lock( tMPI_Spinlock_t &x);
#else
#define tMPI_Spinlock_lock(x)       tMPI_Thread_mutex_lock((x)->lock)
#endif


/** Attempt to acquire spinlock
 *
 * This routine acquires the spinlock if possible, but if 
 * already locked it return an error code immediately.
 *
 *  \param x     Spinlock pointer
 *
 * \return 0 if the mutex was available so we could lock it,
 *         otherwise a non-zero integer (1) if the lock is busy.
 */
#ifdef DOXYGEN
void tMPI_Spinlock_trylock( tMPI_Spinlock_t &x);
#else
#define tMPI_Spinlock_trylock(x)    tMPI_Thread_mutex_trylock((x)->lock)
#endif

/** Release spinlock
 *
 *  \param x     Spinlock pointer
 *
 *  Unlocks the spinlock, regardless if which thread locked it.
 */
#ifdef DOXYGEN
void tMPI_Spinlock_unlock( tMPI_Spinlock_t &x);
#else
#define tMPI_Spinlock_unlock(x)     tMPI_Thread_mutex_unlock((x)->lock)
#endif



/** Check if spinlock is locked
 *
 *  This routine returns immediately with the lock status.
 *
 *  \param x  Spinlock pointer
 *
 *  \return 1 if the spinlock is locked, 0 otherwise.
 */
static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    int rc;
    
    if(tMPI_Spinlock_trylock(x) != 0)
    {
        /* It was locked */
        return 1;
    }
    else
    {
        /* We just locked it */
        tMPI_Spinlock_unlock(x);
        return 0;
    }
}

/** Wait for a spinlock to become available
 *
 *  This routine blocks until the spinlock is unlocked, 
 *  but in contrast to tMPI_Spinlock_lock() it returns without 
 *  trying to lock the spinlock.
 *
 *  \param x  Spinlock pointer
 */
static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    int rc;
    
    tMPI_Spinlock_lock(x);
    /* Got the lock now, so the waiting is over */
    tMPI_Spinlock_unlock(x);
}


#endif




/** Spinlock-based barrier type
 *
 *  This barrier has the same functionality as the standard
 *  tMPI_Thread_barrier_t, but since it is based on spinlocks
 *  it provides faster synchronization at the cost of busy-waiting.
 *
 *  Variables of this type should be initialized by calling
 *  tMPI_Spinlock_barrier_init() to set the number of threads
 *  that should be synchronized.
 * 
 * \see
 * - tMPI_Spinlock_barrier_init
 * - tMPI_Spinlock_barrier_wait
 */
typedef struct tMPI_Spinlock_barrier
{
        tMPI_Atomic_t      count;     /*!< Number of threads remaining     */
        int               threshold; /*!< Total number of threads         */
        volatile int      cycle;     /*!< Current cycle (alternating 0/1) */
}
tMPI_Spinlock_barrier_t;
 



/** Initialize spinlock-based barrier
 *
 *  \param barrier  Pointer to _spinlock_ barrier. Note that this is not
 *                  the same datatype as the full, thread based, barrier.
 *  \param count    Number of threads to synchronize. All threads
 *                  will be released after \a count calls to 
 *                  tMPI_Spinlock_barrier_wait().  
 */
static inline void tMPI_Spinlock_barrier_init(
                                    tMPI_Spinlock_barrier_t *       barrier,
                                    int                              count)
{
        barrier->threshold = count;
        barrier->cycle     = 0;
        tMPI_Atomic_set(&(barrier->count),count);
}




/** Perform busy-waiting barrier synchronization
*
*  This routine blocks until it has been called N times,
*  where N is the count value the barrier was initialized with.
*  After N total calls all threads return. The barrier automatically
*  cycles, and thus requires another N calls to unblock another time.
*
*  Note that spinlock-based barriers are completely different from
*  standard ones (using mutexes and condition variables), only the 
*  functionality and names are similar.
*
*  \param barrier  Pointer to previously create barrier.
*
*  \return The last thread returns -1, all the others 0.
*/
static inline int tMPI_Spinlock_barrier_wait(
                                tMPI_Spinlock_barrier_t *barrier)
{
    int    cycle;
    int    status;
    /*int    i;*/

    /* We don't need to lock or use atomic ops here, since the cycle index 
     * cannot change until after the last thread has performed the check
     * further down. Further, they cannot reach this point in the next 
     * barrier iteration until all of them have been released, and that 
     * happens after the cycle value has been updated.
     *
     * No synchronization == fast synchronization.
     */
    cycle = barrier->cycle;

    /* Decrement the count atomically and check if it is zero.
     * This will only be true for the last thread calling us.
     */
    if( tMPI_Atomic_add_return( &(barrier->count), -1 ) <= 0)
    { 
        tMPI_Atomic_set(&(barrier->count), barrier->threshold);
        barrier->cycle = !barrier->cycle;

        status = -1;
    }
    else
    {
        /* Wait until the last thread changes the cycle index.
         * We are both using a memory barrier, and explicit
         * volatile pointer cast to make sure the compiler
         * doesn't try to be smart and cache the contents.
         */
        do
        { 
            tMPI_Atomic_memory_barrier();
        } 
        while( *(volatile int *)(&(barrier->cycle)) == cycle);

        status = 0;
    }
    return status;
}

#ifdef inline_defined_in_atomic
#undef inline
#endif


#ifdef __cplusplus
}
#endif


#endif /* _TMPI_ATOMIC_H_ */
