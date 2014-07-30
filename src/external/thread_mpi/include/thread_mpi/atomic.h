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

#ifndef TMPI_ATOMIC_H_
#define TMPI_ATOMIC_H_

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

#include "visibility.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif

/* Setting TMPI_ATOMICS_DISABLED permits the build to enforce that no
 * atomic operations are used. This is used when building to run
 * ThreadSanitzer.
 *
 * It could also be useful as a temporary measure on some
 * compiler+hardware for which the detection below fails to produce a
 * correct result. Performance will be greatly improved by using
 * whatever atomic operations are available, so make sure such a
 * measure is only temporary! */
#ifdef TMPI_ATOMICS_DISABLED

#ifndef DOXYGEN
#define TMPI_NO_ATOMICS
#endif

#else

/* first check for gcc/icc platforms.
   Some compatible compilers, like icc on linux+mac will take this path,
   too */
#if ( (defined(__GNUC__) || defined(__PATHSCALE__) || defined(__PGI)) && \
    (!defined(__xlc__)) && (!defined(_CRAYC)) && (!defined(TMPI_TEST_NO_ATOMICS)) )

#ifdef __GNUC__
#define TMPI_GCC_VERSION (__GNUC__ * 10000 \
                          + __GNUC_MINOR__ * 100 \
                          + __GNUC_PATCHLEVEL__)
#endif

/* now check specifically for several architectures: */
#if ((defined(__i386__) || defined(__x86_64__)) && !defined(__OPEN64__))
/* first x86: */
#include "atomic/gcc_x86.h"

#elif (defined(__ia64__))
/* then ia64: */
#include "atomic/gcc_ia64.h"

/* for now we use gcc intrinsics on gcc: */
/*#elif (defined(__powerpc__) || (defined(__ppc__)) )*/
/*#include "atomic/gcc_ppc.h"*/

#elif defined(__FUJITSU) && ( defined(__sparc_v9__) || defined (__sparcv9) )

/* Fujitsu FX10 SPARC compiler */
#include "atomic/fujitsu_sparc64.h"

#else
/* otherwise, there's a generic gcc intrinsics version: */
#include "atomic/gcc.h"

#endif /* end of check for gcc specific architectures */

/* not gcc: */
#elif (defined(_MSC_VER) && (_MSC_VER >= 1200) && \
    (!defined(TMPI_TEST_NO_ATOMICS)) )

/* Microsoft Visual C on x86, define taken from FFTW who got it from
   Morten Nissov. icc on windows will take this path.  */
#include "atomic/msvc.h"

#elif ( (defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM))  && \
    (defined(__powerpc__) || defined(__ppc__)) && \
    (!defined(TMPI_TEST_NO_ATOMICS)) )

/* PowerPC using xlC intrinsics.  */

#include "atomic/xlc_ppc.h"

#elif ( ( defined(__xlC__)  || defined(__xlc__) ) && \
    (!defined(TMPI_TEST_NO_ATOMICS)) )
/* IBM xlC compiler */
#include "atomic/xlc_ppc.h"


#elif (defined (__sun) && (defined(__sparcv9) || defined(__sparc)) && \
    (!defined(TMPI_TEST_NO_ATOMICS)) )
/* Solaris on SPARC (Sun C Compiler, Solaris Studio) */
#include "atomic/suncc-sparc.h"

#elif defined(__FUJITSU) && defined(__sparc__)

/* Fujitsu FX10 SPARC compiler requires gcc compatibility with -Xg */
#warning Atomics support for Fujitsu FX10 compiler requires -Xg (gcc compatibility)
#define TMPI_NO_ATOMICS

#elif defined(_CRAYC)

/* Cray compiler */
#include "atomic/cce.h"
#else

#ifndef DOXYGEN
/** Indicates that no support for atomic operations is present. */
#define TMPI_NO_ATOMICS
#endif

#endif /* platform-specific checks */

#endif /* TMPI_NO_ATOMICS */

#ifdef TMPI_NO_ATOMICS

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

#ifndef DOXYGEN
/* signal that they exist */
#define TMPI_HAVE_ACQ_REL_BARRIERS
#endif

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
 *  - Sparc64, using Fujitsu compilers.
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
    int value; /**< The atomic value.*/
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
    void *value; /**< The atomic pointer. */
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
typedef struct tMPI_Spinlock *tMPI_Spinlock_t;

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
#define TMPI_SPINLOCK_INITIALIZER   { NULL }

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
TMPI_EXPORT
int tMPI_Atomic_get(const tMPI_Atomic_t *a);

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
TMPI_EXPORT
void tMPI_Atomic_set(tMPI_Atomic_t *a, int i);


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
TMPI_EXPORT
void* tMPI_Atomic_ptr_get(const tMPI_Atomic_ptr_t *a);




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
TMPI_EXPORT
void tMPI_Atomic_ptr_set(tMPI_Atomic_ptr_t *a, void *p);

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
TMPI_EXPORT
int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i);
#ifndef DOXYGEN
#define TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN
#endif



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
TMPI_EXPORT
int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i);
#ifndef DOXYGEN
#define TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD
#endif



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
TMPI_EXPORT
int tMPI_Atomic_cas(tMPI_Atomic_t *a, int old_val, int new_val);




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
TMPI_EXPORT
int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t * a, void *old_val,
                        void *new_val);

/** Atomic swap operation.

   Atomically swaps the data in the tMPI_Atomic_t operand with the value of b.
   Note: This has no good assembly counterparts on many architectures, so
         it might not be faster than a repreated CAS.

   \param a  Pointer to atomic type
   \param b  Value to swap
   \return the original value of a
 */
TMPI_EXPORT
int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b);

/** Atomic swap pointer operation.

   Atomically swaps the pointer in the tMPI_Atomic_ptr_t operand with the
   value of b.
   Note: This has no good assembly counterparts on many architectures, so
         it might not be faster than a repreated CAS.

   \param a  Pointer to atomic type
   \param b  Value to swap
   \return the original value of a
 */
TMPI_EXPORT
void *tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b);
#ifndef DOXYGEN
#define TMPI_ATOMIC_HAVE_NATIVE_SWAP
#endif


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
TMPI_EXPORT
void tMPI_Spinlock_init( tMPI_Spinlock_t *x);
#ifndef DOXYGEN
#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK
#endif

/** Acquire spinlock
 *
 *  This routine blocks until the spinlock is available, and
 *  the locks it again before returning.
 *
 *  \param x     Spinlock pointer
 */
TMPI_EXPORT
void tMPI_Spinlock_lock( tMPI_Spinlock_t *x);


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
TMPI_EXPORT
int tMPI_Spinlock_trylock( tMPI_Spinlock_t *x);

/** Release spinlock
 *
 *  \param x     Spinlock pointer
 *
 *  Unlocks the spinlock, regardless if which thread locked it.
 */
TMPI_EXPORT
void tMPI_Spinlock_unlock( tMPI_Spinlock_t *x);



/** Check if spinlock is locked
 *
 *  This routine returns immediately with the lock status.
 *
 *  \param x  Spinlock pointer
 *
 *  \return 1 if the spinlock is locked, 0 otherwise.
 */
TMPI_EXPORT
int tMPI_Spinlock_islocked( tMPI_Spinlock_t *x);

/** Wait for a spinlock to become available
 *
 *  This routine blocks until the spinlock is unlocked,
 *  but in contrast to tMPI_Spinlock_lock() it returns without
 *  trying to lock the spinlock.
 *
 *  \param x  Spinlock pointer
 */
TMPI_EXPORT
void tMPI_Spinlock_wait(tMPI_Spinlock_t *x);


#endif /* TMPI_NO_ATOMICS */

/* now define all the atomics that are not avaible natively. These
   are done on the assumption that a native CAS does exist. */
#include "atomic/derived.h"

/* this allows us to use the inline keyword without breaking support for
   some compilers that don't support it: */
#ifdef inline_defined_in_atomic
#undef inline
#endif

#if !defined(TMPI_NO_ATOMICS) && !defined(TMPI_ATOMICS)
/* Set it here to make sure the user code can check this without having to have
   a config.h */
/** Indicates that support for atomic operations is present. */
#define TMPI_ATOMICS
#endif


#ifdef __cplusplus
}
#endif


#endif /* TMPI_ATOMIC_H_ */
