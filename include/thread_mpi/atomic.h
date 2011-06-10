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
 *  locking, e.g. tMPI_Atomic_fetch_add() or tMPI_Atomic_cas().
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


/* first check for gcc/icc platforms. 
   Some compatible compilers, like icc on linux+mac will take this path, 
   too */
#if ( (defined(__GNUC__) || defined(__PATHSCALE__) || defined(__PGI)) && (!defined(__xlc__)) )

/* now check specifically for several architectures: */
#if (defined(i386) || defined(__x86_64__)) 
/* first x86: */
#include "atomic/gcc_x86.h"
/*#include "atomic/gcc.h"*/

#elif (defined(__ia64__))
/* then ia64: */
#include "atomic/gcc_ia64.h"

/* for now we use gcc intrinsics on gcc: */
/*#elif (defined(__powerpc__) || (defined(__ppc__)) )*/
/*#include "atomic/gcc_ppc.h"*/

#else
/* otherwise, there's a generic gcc intrinsics version: */
#include "atomic/gcc.h"

#endif /* end of check for gcc specific architectures */

/* not gcc: */
#elif (defined(_MSC_VER) && (_MSC_VER >= 1200))
/* Microsoft Visual C on x86, define taken from FFTW who got it from 
   Morten Nissov. icc on windows will take this path.  */
#include "atomic/msvc.h"

#elif ( (defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM))  && \
        (defined(__powerpc__) || defined(__ppc__)))

/* PowerPC using xlC intrinsics.  */

#include "atomic/xlc_ppc.h"

#elif defined(__xlC__)  || defined(__xlc__)
/* IBM xlC compiler */
#include "atomic/xlc_ppc.h"


#elif defined (__sun) && (defined(__sparcv9) || defined(__sparc))
/* Solaris on SPARC (Sun C Compiler, Solaris Studio) */
#include "atomic/suncc-sparc.h"




#else
/* No atomic operations, use mutex fallback. Documentation is in x86 section */

#ifdef TMPI_CHECK_ATOMICS
#error No atomic operations implemented for this cpu/compiler combination. 
#endif

#define TMPI_NO_ATOMICS


/** Memory barrier operation

 Modern CPUs rely heavily on out-of-order execution, and one common feature
 is that load/stores might be reordered. Also, when using inline assembly
 the compiler might already have loaded the variable we are changing into
 a register, so any update to memory won't be visible.
 
 This command creates a memory barrier, i.e. all memory results before
 it in the code should be visible to all memory operations after it - the
 CPU cannot propagate load/stores across it.

 This barrier is a full barrier: all load and store operations of
 instructions before it are completed, while all load and store operations
 that are in instructions after it won't be done before this barrier.

 \hideinitializer
 */
#define tMPI_Atomic_memory_barrier()

/** Memory barrier operation with acquire semantics

 This barrier is a barrier with acquire semantics: the terminology comes
 from its common use after acquiring a lock: all load/store instructions 
 after this barrier may not be re-ordered to happen before this barrier.

 \hideinitializer
 */
#define tMPI_Atomic_memory_barrier_acq()

/** Memory barrier operation with release semantics

 This barrier is a barrier with release semantics: the terminology comes
 from its common use before releasing a lock: all load/store instructions 
 before this barrier may not be re-ordered to happen after this barrier.

 \hideinitializer
 */
#define tMPI_Atomic_memory_barrier_rel()




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
 * - tMPI_Atomic_cas
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
 * - tMPI_Atomic_ptr_cas
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
static inline int tMPI_Atomic_get(tMPI_Atomic_t &a);
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
static inline void tMPI_Atomic_set(tMPI_Atomic_t *a, int i)
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
static inline void* tMPI_Atomic_ptr_get(tMPI_Atomic_ptr_t &a);
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
static inline void tMPI_Atomic_ptr_set(tMPI_Atomic_t *a, void *p)
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
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
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
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
{
    int old_value;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    old_value  = a->value;
    a->value   = old_value + i;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return old_value;
}



/** Atomic compare-and-swap operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is swapped with the new value, 
 *   and otherwise left unchanged. 
 * 
 *   This is *the* synchronization primitive: it has a consensus number of
 *   infinity, and is available in some form on all modern CPU architectures.
 *   In the words of Herlihy&Shavit (The art of multiprocessor programming),
 *   it is the 'king of all wild things'. 
 *  
 *   In practice, use it as follows: You can start by reading a value 
 *   (without locking anything), perform some calculations, and then
 *   atomically try to update it in memory unless it has changed. If it has
 *   changed you will get an error return code - reread the new value
 *   an repeat the calculations in that case.
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param old_val  Integer value read from the atomic type at an earlier point
 *   \param new_val  New value to write to the atomic type if it currently is
 *                   identical to the old value.
 *
 *   \return    True (1) if the swap occurred: i.e. if the value in a was equal
 *              to old_val. False (0) if the swap didn't occur and the value
 *              was not equal to old_val.
 * 
 *   \note   The exchange occured if the return value is identical to \a old.
 */
static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int old_val, int new_val)
{
    int t=0;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    if (a->value == old_val)
    {
        a->value = new_val;
        t=1;
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}




/** Atomic pointer compare-and-swap operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is swapped with the new value, 
 *   and otherwise left unchanged. 
 *  
 *   This is essential for implementing wait-free lists and other data
 *   structures. See 'tMPI_Atomic_cas()'.
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param old_val  Pointer value read from the atomic type at an earlier point
 *   \param new_val  New value to write to the atomic type if it currently is
 *                   identical to the old value.
 *
 *   \return    True (1) if the swap occurred: i.e. if the value in a was equal
 *              to old_val. False (0) if the swap didn't occur and the value
 *              was not equal to old_val.
 * 
 *   \note   The exchange occured if the return value is identical to \a old.
 */
static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t * a, void *old_val,
                                      void *new_val)
{
    int t=0;
    
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    if (a->value == old_val)
    {
        a->value = new_val;
        t=1;
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
static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
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
static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    int rc;
    
    tMPI_Spinlock_lock(x);
    /* Got the lock now, so the waiting is over */
    tMPI_Spinlock_unlock(x);
}


#endif



/* only do this if there was no better solution */
#ifndef TMPI_HAVE_SWAP
/** Atomic swap operation.

  Atomically swaps the data in the tMPI_Atomic_t operand with the value of b.
  NOTE: DON'T USE YET! (This has no good asm counterparts on many architectures).

  \param a  Pointer to atomic type
  \param b  Value to swap 
  \return the original value of a
*/
static inline int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    int oldval;
    do
    {
        oldval=(int)(a->value);
    } while(!tMPI_Atomic_cas(a, oldval, b));
    return oldval;
}
/** Atomic swap pointer operation.

  Atomically swaps the pointer in the tMPI_Atomic_ptr_t operand with the 
  value of b.
  NOTE: DON'T USE YET! (This has no good asm counterparts on many architectures).

  \param a  Pointer to atomic type
  \param b  Value to swap 
  \return the original value of a
*/
static inline void *tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    void *oldval;
    do
    {
        oldval=(void*)(a->value);
    } while(!tMPI_Atomic_ptr_cas(a, oldval, b));
    return oldval;
}
#endif

/* only define this if there were no separate acquire and release barriers */
#ifndef TMPI_HAVE_ACQ_REL_BARRIERS

/* if they're not defined explicitly, we just make full barriers out of both */
#define tMPI_Atomic_memory_barrier_acq tMPI_Atomic_memory_barrier
#define tMPI_Atomic_memory_barrier_rel tMPI_Atomic_memory_barrier

#endif

/* this allows us to use the inline keyword without breaking support for 
   some compilers that don't support it: */
#ifdef inline_defined_in_atomic
#undef inline
#endif


#ifdef __cplusplus
}
#endif


#endif /* _TMPI_ATOMIC_H_ */
