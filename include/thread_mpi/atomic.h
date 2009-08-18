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
 *  @brief Atomic operations for fast SMP synchronization
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


#include <stdio.h>



#ifdef __cplusplus
extern "C" 
{  
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)
#endif


#if ( ( (defined(__GNUC__) || defined(__INTEL_COMPILER) ||  \
       defined(__PATHSCALE__)) && (defined(i386) || defined(__x86_64__)) ) \
      || defined (DOXYGEN) )


#include <limits.h>
#include <stdint.h>
/* This code is executed for x86 and x86-64, with these compilers:
 * GNU
 * Intel 
 * Pathscale
 * All these support GCC-style inline assembly. 
 * We also use this section for the documentation.
 */

/*! \brief Memory barrier operation

 Modern CPUs rely heavily on out-of-order execution, and one common feature
 is that load/stores might be reordered. Also, when using inline assembly
 the compiler might already have loaded the variable we are changing into
 a register, so any update to memory won't be visible.
 
 This command creates a memory barrier, i.e. all memory results before
 it in the code should be visible to all memory operations after it - the
 CPU cannot propagate load/stores across it.

 \hideinitializer
 */
#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("": : :"memory")

/* Only gcc and Intel support this check, otherwise set it to true (skip doc) */
#if (!defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined DOXYGEN)
#define __builtin_constant_p(i) (1)
#endif


/*! \brief Atomic operations datatype
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
 *  - x86 or x86_64, using GNU compilers
 *  - x86 or x86_64, using Intel compilers 
 *  - x86 or x86_64, using Pathscale compilers
 *  - Itanium, using GNU compilers 
 *  - Itanium, using Intel compilers
 *  - Itanium, using HP compilers
 *  - PowerPC, using GNU compilers 
 *  - PowerPC, using IBM AIX compilers 
 *  - PowerPC, using IBM compilers >=7.0 under Linux or Mac OS X.
 */
typedef struct tMPI_Atomic
{
        volatile int       value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

/*! \brief Atomic pointer type equivalent to tMPI_Atomic_t
 *
 * Useful for lock-free and wait-free data structures.
 * The only operations available for this type are
 * tMPI_Atomic_ptr_get
 * tMPI_Atomic_ptr_set
 * tMPI_Atomic_ptr_cmpxch
*/
typedef struct tMPI_Atomic_ptr
{
        void* volatile*    value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


/*! \brief Spinlock
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
 */
typedef struct tMPI_Spinlock
{
    volatile unsigned int  lock;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;





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
#define TMPI_SPINLOCK_INITIALIZER   { 1 }



/*! \brief Return value of an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable to read
 *  \return     Integer value of the atomic variable
 *
 *  \hideinitializer
 */
#define tMPI_Atomic_get(a)  ((a)->value) 

 
/*! \brief Write value to an atomic integer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable
 *  \param  i   Integer to set the atomic variable to.
 *
 *  \hideinitializer
 */
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))


/*! \brief Return value of an atomic pointer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable to read
 *  \return     Pointer value of the atomic variable
 *
 *  \hideinitializer
 */
#define tMPI_Atomic_ptr_get(a)  ((a)->value) 

 
/*! \brief Write value to an atomic pointer 
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param  a   Atomic variable
 *  \param  i   Pointer value to set the atomic variable to.
 *
 *  \hideinitializer
 */
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))


 
/*! \brief Add integer to atomic variable
 *
 *  Also implements proper memory barriers when necessary.
 *  The actual implementation is system-dependent.
 *
 *  \param a   atomic datatype to modify
 *  \param i   integer to increment with. Use i<0 to subtract atomically.
 *
 *  \return The new value (after summation).
 */
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *     a, 
                                        volatile int       i)
{
    int __i;
    
    __i = i;
    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(i) :"m"(a->value), "0"(i));
    return i + __i;
}  
  

/*! \brief Add to variable, return the old value.
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
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *     a,
                                       volatile int       i)
{
#if 0
    int __i;

    __i = i;
#endif
    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(i) :"m"(a->value), "0"(i));
    return i;
}


/*! \brief Atomic compare-exchange operation
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
 *   \param oldval   Integer value read from the atomic type at an earlier point
 *   \param newval   New value to write to the atomic type if it currently is
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
static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *    a, 
                                     int               oldval,
                                     int               newval)
{
    volatile unsigned long prev;
    
    __asm__ __volatile__("lock ; cmpxchgl %1,%2"
                         : "=a"(prev)
                         : "q"(newval), "m"(a->value), "0"(oldval)
                         : "memory");
    
    return prev;
}


/*! \brief Atomic pointer compare-exchange operation
 *
 *   The \a old value is compared with the memory value in the atomic datatype.
 *   If the are identical, the atomic type is updated to the new value, 
 *   and otherwise left unchanged. 
 *  
 *   This is essential for implementing wait-free lists and other data
 *   structures. 
 *
 *   \param a        Atomic datatype ('memory' value)
 *   \param oldval   Pointer value read from the atomic type at an earlier point
 *   \param newval   New value to write to the atomic type if it currently is
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
static inline void* volatile* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t* a, 
                                                    void*             oldval,
                                                    void*             newval)
{
    void* volatile *prev;
#ifndef __x86_64__ 
    __asm__ __volatile__("lock ; cmpxchgl %1,%2"
                         : "=a"(prev)
                         : "q"(newval), "m"(a->value), "0"(oldval)
                         : "memory");
#else 
    __asm__ __volatile__("lock ; cmpxchgq %1,%2"
                         : "=a"(prev)
                         : "q"(newval), "m"(a->value), "0"(oldval)
                         : "memory");
#endif
    return prev;
}



/*! \brief Initialize spinlock
 *
 *  In theory you can call this from multiple threads, but remember
 *  that we don't check for errors. If the first thread proceeded to
 *  lock the spinlock after initialization, the second will happily
 *  overwrite the contents and unlock it without warning you.
 *
 *  \param x      Spinlock pointer.
 */
static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *   x)
{
    x->lock = 1;
}



/*! \brief Acquire spinlock
 *
 *  This routine blocks until the spinlock is available, and
 *  the locks it again before returning.
 *
 *  \param x     Spinlock pointer
 */
static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *  x)
{
        __asm__ __volatile__("\n1:\t" 
                             "lock ; decb %0\n\t" 
                             "jns 3f\n" 
                             "2:\t" 
                             "rep;nop\n\t" 
                             "cmpb $0,%0\n\t" 
                             "jle 2b\n\t" 
                             "jmp 1b\n" 
                             "3:\n\t" 
                             :"=m" (x->lock) : : "memory"); 
}


/*! \brief Attempt to acquire spinlock
 *
 * This routine acquires the spinlock if possible, but if 
 * already locked it return an error code immediately.
 *
 *  \param x     Spinlock pointer
 *
 * \return 0 if the mutex was available so we could lock it,
 *         otherwise a non-zero integer (1) if the lock is busy.
 */
static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *  x)
{
        char old_value;
        
    __asm__ __volatile__("xchgb %b0,%1"
                         :"=q" (old_value), "=m" (x->lock)
                         :"0" (0) : "memory");
    return (old_value <= 0);
}


/*! \brief Release spinlock
 *
 *  \param x     Spinlock pointer
 *
 *  Unlocks the spinlock, regardless if which thread locked it.
 */
static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *  x)
{
        char old_value = 1;
        
        __asm__ __volatile__(
                         "xchgb %b0, %1" 
                         :"=q" (old_value), "=m" (x->lock) 
                         :"0" (old_value) : "memory"
                         );
}
 

/*! \brief Check if spinlock is locked
 *
 *  This routine returns immediately with the lock status.
 *
 *  \param x  Spinlock pointer
 *
 *  \return 1 if the spinlock is locked, 0 otherwise.
 */
static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *  x)
{
    return (*(volatile signed char *)(&(x)->lock) <= 0);
}


/*! \brief Wait for a spinlock to become available
 *
 *  This routine blocks until the spinlock is unlocked, 
 *  but in contrast to tMPI_Spinlock_lock() it returns without 
 *  trying to lock the spinlock.
 *
 *  \param x  Spinlock pointer
 */
static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    do 
    {
        tMPI_Atomic_memory_barrier(); 
    } 
    while(tMPI_Spinlock_islocked(x));
}


#elif ( defined(__GNUC__) && (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using proper GCC inline assembly. 
 * Recent versions of xlC (>=7.0) _partially_ support this, but since it is
 * not 100% compatible we provide a separate implementation for xlC in
 * the next section.
 */

/* Compiler-dependent stuff: GCC memory barrier */
#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("": : :"memory")



typedef struct tMPI_Atomic
{
        volatile int       value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile*     value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


typedef struct tMPI_Spinlock
{
    volatile unsigned int   lock;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define tMPI_Atomic_get(a)        ((a)->value) 
#define tMPI_Atomic_set(a,i)     (((a)->value) = (i))

#define tMPI_Atomic_ptr_get(a)    ((a)->value) 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))

static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *    a, 
                                        int               i)
{
    int t;
    
    __asm__ __volatile__("1:     lwarx   %0,0,%2\n"
                         "\tadd     %0,%1,%0\n"
                         "\tstwcx.  %0,0,%2 \n"
                         "\tbne-    1b\n"
                         "\tisync\n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value)
                         : "cc" , "memory");
    return t;
}



static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *     a,
                                       int                i)
{
    int t;
    
    __asm__ __volatile__("\teieio\n"
                         "1:     lwarx   %0,0,%2\n"                         
                         "\tadd     %0,%1,%0\n"
                         "\tstwcx.  %0,0,%2 \n"
                         "\tbne-    1b\n"
                         "\tisync\n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value)
                         : "cc", "memory");
    
    return (t - i);    
}


static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *       a,
                                     int                  oldval,
                                     int                  newval)
{
    int prev;
    
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\tcmpw    0,%0,%3 \n"
                          "\tbne     2f \n"
                          "\tstwcx.  %4,0,%2 \n"
                          "bne-    1b\n"
                          "\tsync\n"
                          "2:\n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value)
                          : "cc", "memory");
    
    return prev;
}


static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t *   a,
                                           void *               oldval,
                                           void *               newval)
{
    void *prev;
   
#if (!defined(__PPC64__)) && (!defined(__ppc64))
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\tcmpw    0,%0,%3 \n"
                          "\tbne     2f \n"
                          "\tstwcx.  %4,0,%2 \n"
                          "bne-    1b\n"
                          "\tsync\n"
                          "2:\n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value)
                          : "cc", "memory");
#else
    __asm__ __volatile__ ("1:    ldarx   %0,0,%2 \n"
                          "\tcmpd    0,%0,%3 \n"
                          "\tbne     2f \n"
                          "\tstdcx.  %4,0,%2 \n"
                          "bne-    1b\n"
                          "\tsync\n"
                          "2:\n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value)
                          : "cc", "memory");
#endif
    return prev;
}




static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}



static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *  x)
{
    unsigned int tmp;
    
    __asm__ __volatile__("\tb      1f\n"
                         "2:      lwzx    %0,0,%1\n"
                         "\tcmpwi   0,%0,0\n"
                         "\tbne+    2b\n"
                         "1:      lwarx   %0,0,%1\n"
                         "\tcmpwi   0,%0,0\n"
                         "\tbne-    2b\n"
                         "\tstwcx.  %2,0,%1\n"
                         "bne-    2b\n"
                         "\tisync\n"
                         : "=&r"(tmp)
                         : "r"(&x->lock), "r"(1)
                         : "cr0", "memory");
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *  x)
{
    unsigned int old, t;
    unsigned int mask = 1;
    volatile unsigned int *p = &x->lock;
    
    __asm__ __volatile__("\teieio\n"
                         "1:      lwarx   %0,0,%4 \n"
                         "\tor      %1,%0,%3 \n"
                         "\tstwcx.  %1,0,%4 \n"
                         "\tbne     1b\n"
                         "\tsync\n"
                         : "=&r" (old), "=&r" (t), "=m" (*p)
                         : "r" (mask), "r" (p), "m" (*p)
                         : "cc", "memory");
    
    return ((old & mask) != 0);    
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *  x)
{
    __asm__ __volatile__("\teieio\n": : :"memory");
    x->lock = 0;
}


static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return ( x->lock != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    do 
    {
        tMPI_Atomic_memory_barrier(); 
    }
    while(tMPI_Spinlock_islocked(x));
}



#elif ( (defined(__IBM_GCC_ASM) || defined(__IBM_STDCPP_ASM))  && \
        (defined(__powerpc__) || defined(__ppc__)))
/* PowerPC using xlC inline assembly. 
 * Recent versions of xlC (>=7.0) _partially_ support GCC inline assembly
 * if you use the option -qasm=gcc but we have had to hack things a bit, in 
 * particular when it comes to clobbered variables. Since this implementation
 * _could_ be buggy, we have separated it from the known-to-be-working gcc
 * one above.
 */

/* memory barrier - no idea how to create one with xlc! */
#define tMPI_Atomic_memory_barrier()



typedef struct tMPI_Atomic
{
        volatile int       value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;


typedef struct tMPI_Atomic_ptr
{
        void* volatile*     value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;



typedef struct tMPI_Spinlock
{
    volatile unsigned int   lock;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define tMPI_Atomic_get(a)   ((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))
#define tMPI_Atomic_ptr_get(a)   ((a)->value) 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))


static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *    a, 
                                        int               i)
{
    int t;
    
    __asm__ __volatile__("1:     lwarx   %0,0,%2 \n"
                         "\t add     %0,%1,%0 \n"
                         "\t stwcx.  %0,0,%2 \n"
                         "\t bne-    1b \n"
                         "\t isync \n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value) );
    return t;
}



static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *     a,
                                       int                i)
{
    int t;
    
    __asm__ __volatile__("\t eieio\n"
                         "1:     lwarx   %0,0,%2 \n"                         
                         "\t add     %0,%1,%0 \n"
                         "\t stwcx.  %0,0,%2 \n"
                         "\t bne-    1b \n"
                         "\t isync \n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value));
    
    return (t - i);    
}


static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *       a,
                                     int                  oldval,
                                     int                  newval)
{
    int prev;
    
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\t cmpw    0,%0,%3 \n"
                          "\t bne     2f \n"
                          "\t stwcx.  %4,0,%2 \n"
                          "\t bne-    1b \n"
                          "\t sync \n"
                          "2: \n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value));
    
    return prev;
}

static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t *   a,
                                           void*                oldval,
                                           void*                newval)
{
    void* prev;
   

#if (!defined(__PPC64__)) && (!defined(__ppc64))
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\t cmpw    0,%0,%3 \n"
                          "\t bne     2f \n"
                          "\t stwcx.  %4,0,%2 \n"
                          "\t bne-    1b \n"
                          "\t sync \n"
                          "2: \n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value));
    
#else
    __asm__ __volatile__ ("1:    ldarx   %0,0,%2 \n"
                          "\t cmpd    0,%0,%3 \n"
                          "\t bne     2f \n"
                          "\t stdcx.  %4,0,%2 \n"
                          "\t bne-    1b \n"
                          "\t sync \n"
                          "2: \n"
                          : "=&r" (prev), "=m" (a->value)
                          : "r" (&a->value), "r" (oldval), "r" (newval), 
                            "m" (a->value));
#endif
    return prev;
}


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}



static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *  x)
{
    unsigned int tmp;
    
    __asm__ __volatile__("\t b      1f \n"
                         "2:      lwzx    %0,0,%1 \n"
                         "\t cmpwi   0,%0,0 \n"
                         "\t bne+    2b \n"
                         "1:      lwarx   %0,0,%1 \n"
                         "\t cmpwi   0,%0,0 \n"
                         "\t bne-    2b \n"
                         "\t stwcx.  %2,0,%1 \n"
                         "\t bne-    2b \n"
                         "\t isync\n"
                         : "=&r"(tmp)
                         : "r"(&x->lock), "r"(1));
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *  x)
{
    unsigned int old, t;
    unsigned int mask = 1;
    volatile unsigned int *p = &x->lock;
    
    __asm__ __volatile__("\t eieio\n"
                         "1:      lwarx   %0,0,%4 \n"
                         "\t or      %1,%0,%3 \n"
                         "\t stwcx.  %1,0,%4 \n"
                         "\t bne     1b \n"
                         "\t sync \n"
                         : "=&r" (old), "=&r" (t), "=m" (*p)
                         : "r" (mask), "r" (p), "m" (*p));
    
    return ((old & mask) != 0);    
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *  x)
{
    __asm__ __volatile__("\t eieio \n");
    x->lock = 0;
}


static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return ( x->lock != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    
    do 
    {
        tMPI_Atomic_memory_barrier();
    }
    while(spin_islocked(x));
}




#elif (defined(__ia64__) && (defined(__GNUC__) || defined(__INTEL_COMPILER)))
/* ia64 with GCC or Intel compilers. Since we need to define everything through
* cmpxchg and fetchadd on ia64, we merge the different compilers and only 
* provide different implementations for that single function. 
* Documentation? Check the gcc/x86 section.
*/


typedef struct tMPI_Atomic
{
    volatile int       value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
    void* volatile    value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


typedef struct tMPI_Spinlock
{
    volatile unsigned int   lock; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define tMPI_Atomic_get(a)   ((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))

#define tMPI_Atomic_ptr_get(a)   ((a)->value) 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (i))


/* Compiler thingies */
#ifdef __INTEL_COMPILER
/* prototypes are neccessary for these intrisics: */
#include <ia64intrin.h>
void __memory_barrier(void);
int _InterlockedCompareExchange(volatile int *dest, int xchg, int comp);
/*void* _InterlockedCompareExchangePointer(void* volatile **dest, void* xchg, 
                                         void* comp);*/
unsigned __int64 __fetchadd4_rel(unsigned int *addend, const int increment);
/* ia64 memory barrier */
#define tMPI_Atomic_memory_barrier() __memory_barrier()
/* ia64 cmpxchg */
#define tMPI_Atomic_cmpxchg(a, oldval, newval) _InterlockedCompareExchange(&((a)->value),newval,oldval)
/* ia64 pointer cmpxchg */
#define tMPI_Atomic_ptr_cmpxchg(a, oldval, newval) _InterlockedCompareExchangePointer(&((a)->value),newval,oldval)

/*#define tMPI_Atomic_ptr_cmpxchg(a, oldval, newval) __sync_val_compare_and_swap(&((a)->value),newval,oldval)*/


/* ia64 fetchadd, but it only works with increments +/- 1,4,8,16 */
#define tMPI_ia64_fetchadd(a, inc)  __fetchadd4_rel(a, inc)

#elif defined __GNUC__  
/* ia64 memory barrier */
#  define tMPI_Atomic_memory_barrier() asm volatile ("":::"memory")
/* ia64 cmpxchg */
static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *   a,
                                     int              oldval,
                                     int              newval)
{
#if GCC_VERSION < 40200
    volatile int res;
    asm volatile ("mov ar.ccv=%0;;" :: "rO"(oldval));
    asm volatile ("cmpxchg4.acq %0=[%1],%2,ar.ccv":                    
                  "=r"(res) : "r"(&a->value), "r"(newval) : "memory"); 
                          
    return res;
#else
    return __sync_val_compare_and_swap( &(a->value), oldval, newval);
#endif
}

/* ia64 ptr cmpxchg */
static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t * a,
                                           void*              oldval,
                                           void*              newval)
{
#if GCC_VERSION < 40200
    void* volatile* res;
    asm volatile ("mov ar.ccv=%0;;" :: "rO"(oldval));
    asm volatile ("cmpxchg8.acq %0=[%1],%2,ar.ccv":                    
                  "=r"(res) : "r"(&a->value), "r"(newval) : "memory"); 
                          
    return (void*)res;
#else
    return (void*)__sync_val_compare_and_swap( &(a->value), oldval, newval);
#endif
}


/* fetchadd, but on ia64 it only works with increments +/- 1,4,8,16 */
#define tMPI_ia64_fetchadd(a, inc)                                             \
({  unsigned long res;                                                        \
    asm volatile ("fetchadd4.rel %0=[%1],%2"                                  \
                  : "=r"(res) : "r"(a), "i" (inc) : "memory");                \
                  res;                                                        \
})


#else /* Unknown compiler */
#  error Unknown ia64 compiler (not GCC or ICC) - modify tMPI_Thread.h!
#endif



static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *       a, 
                                        volatile int         i)
{
    volatile int oldval,newval;    
    volatile int __i = i;

    /* Use fetchadd if, and only if, the increment value can be determined
     * at compile time (otherwise this check is optimized away) and it is
     * a value supported by fetchadd (1,4,8,16,-1,-4,-8,-16).
     */                         
    if (__builtin_constant_p(i) &&
        ( (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) || 
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) ) )
    {
        oldval = tMPI_ia64_fetchadd((unsigned int*)&(a->value),__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = tMPI_Atomic_get(a);
            newval = oldval + i;
        }
        while(tMPI_Atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return newval;
}



static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *     a,
                                       volatile int       i)
{
    volatile int oldval,newval;    
    volatile int __i = i;
    
    /* Use ia64 fetchadd if, and only if, the increment value can be determined
     * at compile time (otherwise this check is optimized away) and it is
     * a value supported by fetchadd (1,4,8,16,-1,-4,-8,-16).
     */                         
    if (__builtin_constant_p(i) &&
        ( (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) || 
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) ) )
    {
        oldval = tMPI_ia64_fetchadd((unsigned int*)&(a->value),__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = tMPI_Atomic_get(a);
            newval = oldval + i;
        }
        while(tMPI_Atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return oldval;
}


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *   x)
{
    tMPI_Atomic_t *a = (tMPI_Atomic_t *) x;
    unsigned long value;                                                 
    value = tMPI_Atomic_cmpxchg(a, 0, 1);                             
    if (value)                                                           
    {                                                                    
        do                                                               
        {                                                                
            while (a->value != 0)   
            {                                                            
                tMPI_Atomic_memory_barrier();                             
            }                                                            
            value = tMPI_Atomic_cmpxchg(a, 0, 1);                       
        }                                                                
        while (value);                                                   
    }                                                                    
} 


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *   x)
{
    return (tMPI_Atomic_cmpxchg( ((tMPI_Atomic_t *)x), 0, 1) != 0);
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *   x)
{
    do
    {
        tMPI_Atomic_memory_barrier(); 
        x->lock = 0;
    } 
    while (0);
}


static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return (x->lock != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    
    do 
    {
        tMPI_Atomic_memory_barrier();
    }
    while(tMPI_Spinlock_islocked(x));
}


#undef tMPI_ia64_fetchadd



#elif (defined(__hpux) || defined(__HP_cc)) && defined(__ia64)
/* HP compiler on ia64 */
#include <machine/sys/inline.h>

#define tMPI_Atomic_memory_barrier() _Asm_mf()

#define tMPI_hpia64_fetchadd(a, i)                           \
    _Asm_fetchadd((_Asm_fasz)_FASZ_W,(_Asm_sem)_SEM_REL,    \
                  (UInt32*)a,(unsigned int) i,              \
                  (_Asm_ldhint)LDHINT_NONE)
 

typedef struct tMPI_Atomic
{
        volatile int       value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile*     value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;



typedef struct tMPI_Spinlock
{
    volatile unsigned int   lock; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *   a,
                                     int              oldval,
                                     int              newval)
{
    int ret;
    
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,(Uint32)oldval,                  
                   (_Asm_fence)(_UP_CALL_FENCE | _UP_SYS_FENCE |         
                                _DOWN_CALL_FENCE | _DOWN_SYS_FENCE));
                   
    ret = _Asm_cmpxchg((_Asm_sz)SZ_W,(_Asm_sem)_SEM_ACQ,(Uint32*)a,    
                       (Uint32)newval,(_Asm_ldhint)_LDHINT_NONE);
                   
    return ret;
}



static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t *  a,
                                           void*               oldval,
                                           void*               newval)
{
    void *ret;

    /* todo: fix this */
    
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,(Uint64)oldval,                  
                   (_Asm_fence)(_UP_CALL_FENCE | _UP_SYS_FENCE |         
                                _DOWN_CALL_FENCE | _DOWN_SYS_FENCE));
                   
    ret = _Asm_cmpxchg((_Asm_sz)SZ_W,(_Asm_sem)_SEM_ACQ,(Uint64)a,    
                       (Uint64)newval,(_Asm_ldhint)_LDHINT_NONE);
                   
    return ret;
}




#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define tMPI_Atomic_get(a)   ((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))


static inline void tMPI_Atomic_add_return(tMPI_Atomic_t *       a, 
                                         int                  i)
{
    int old,new;    
    int __i = i;
    
    /* On HP-UX we don't know any macro to determine whether the increment
     * is known at compile time, but hopefully the call uses something simple
     * like a constant, and then the optimizer should be able to do the job.
     */                         
    if (  (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||  
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) )
    {
        oldval = tMPI_hpia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = tMPI_Atomic_get(a);
            newval = oldval + i;
        }
        while(tMPI_Atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return newval;
}



static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *     a,
                                       int                i)
{
    int oldval,newval;    
    int __i = i;
    
    /* On HP-UX we don't know any macro to determine whether the increment
     * is known at compile time, but hopefully the call uses something simple
     * like a constant, and then the optimizer should be able to do the job.
     */                         
    if (  (__i ==   1) || (__i ==   4)  || (__i ==   8) || (__i ==  16) ||
          (__i ==  -1) || (__i ==  -4)  || (__i ==  -8) || (__i == -16) )
    {
        oldval = tMPI_hpia64_fetchadd(a,__i);
        newval = oldval + i;
    }
    else
    {
        /* Use compare-exchange addition that works with any value */
        do
        {
            oldval = tMPI_Atomic_get(a);
            newval = oldval + i;
        }
        while(tMPI_Atomic_cmpxchg(a,oldval,newval) != oldval);
    }
    return oldval;
}


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}





static inline void tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    int rc;

    rc = _Asm_xchg((_Asm_sz)_SZ_W, (unsigned int *)x, 1        
                    (_Asm_ldhit)_LDHINT_NONE);
    
    return ( (rc>0) ? 1 : 0);
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
    int      status = 1;
    
    do
    {
        if( *((unsigned int *)x->lock) == 0 ) 
        {
            status = tMPI_Spinlock_trylock(x);
        }
    } while( status != 0);
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *   x)
{
    _Asm_fetchadd((_Asm_fasz)_SZ_W,(_Asm_sem)_SEM_REL,                  
                  (unsigned int *)x,-1,(_Asm_ldhint)_LDHINT_NONE);
}



static inline void tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return ( x->lock != 0 );
}



static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    do
    {
        tMPI_Atomic_memory_barrier(); 
    } 
    while(spin_islocked(x));
}


#undef tMPI_hpia64_fetchadd



#elif (defined(_MSC_VER) && (_MSC_VER >= 1200))
/* Microsoft Visual C on x86, define taken from FFTW who got it from Morten Nissov */

#include <windows.h>

#if (!defined(inline)) && (!defined(__cplusplus))
#define inline_defined_in_atomic 1
#define inline __inline
#endif

#define tMPI_Atomic_memory_barrier()


typedef struct tMPI_Atomic
{
        LONG volatile      value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile      value; /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;



typedef struct tMPI_Spinlock
{
    LONG volatile      lock;      /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }




#define tMPI_Atomic_get(a)  ((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))


#define tMPI_Atomic_ptr_get(a)    ((a)->value) 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))


#define tMPI_Atomic_fetch_add(a, i)  \
    InterlockedExchangeAdd((LONG volatile *)(a), (LONG) (i))

#define tMPI_Atomic_add_return(a, i)  \
    ( (i) + InterlockedExchangeAdd((LONG volatile *)(a), (LONG) (i)) )

#define tMPI_Atomic_cmpxchg(a, oldval, newval) \
    InterlockedCompareExchange((LONG volatile *)(a), (LONG) (newval), (LONG) (oldval))

#define tMPI_Atomic_ptr_cmpxchg(a, oldval, newval) \
    InterlockedCompareExchangePointer(&((a)->value), (PVOID) (newval),  \
                                      (PVOID) (oldval))


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *   x)
{
    x->lock = 0;
}

# define tMPI_Spinlock_lock(x)   \
    while((InterlockedCompareExchange((LONG volatile *)(x), 1, 0))!=0)


#define tMPI_Spinlock_trylock(x)   \
    InterlockedCompareExchange((LONG volatile *)(x), 1, 0)


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *   x)
{
    x->lock = 0;
}


static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return (*(volatile signed char *)(&(x)->lock) != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    while(tMPI_Spinlock_islocked(x))
    {
        /*Sleep(0);*/
    }
}



#elif defined(__xlC__) && defined (_AIX)
/* IBM xlC compiler on AIX */
#include <sys/atomic_op.h>


#define tMPI_Atomic_memory_barrier()


typedef struct tMPI_Atomic
{
        volatile int       value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;


typedef struct tMPI_Atomic_ptr
{
        void* volatile*     value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;




typedef struct tMPI_Spinlock
{
    volatile unsigned int   lock;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


static inline int tMPI_Atomic_cmpxchg(tMPI_Atomic_t *    a,
                                     int               oldval,
                                     int               newval)
{
    int t;
    
    if(__check_lock((atomic_p)&a->value, oldval, newval))
    {
        /* Not successful - value had changed in memory. Reload value. */
        t = a->value;
    }
    else
    {
        /* replacement suceeded */
        t = oldval;
    }
    return t;        
}


static inline void* tMPI_Atomic_ptr_cmpxchg(tMPI_Atomic_ptr_t *a,
                                           void*             oldval,
                                           void*             newval)
{
    void *t;
    
    if(__check_lock((atomic_p)&a->value, oldval, newval))
    {
        /* Not successful - value had changed in memory. Reload value. */
        t = a->value;
    }
    else
    {
        /* replacement suceeded */
        t = oldval;
    }
    return t;        
}




static inline void tMPI_Atomic_add_return(tMPI_Atomic_t *       a, 
                                         int                  i)
{
    int oldval,newval;    
    
    do
    {
        oldval = tMPI_Atomic_get(a);
        newval = oldval + i;
    }
    while(__check_lock((atomic_p)&a->value, oldval, newval));

    return newval;
}



static inline void tMPI_Atomic_fetch_add(tMPI_Atomic_t *       a, 
                                        int                  i)
{
    int oldval,newval;    
    
    do
    {
        oldval = tMPI_Atomic_get(a);
        newval = oldval + i;
    }
    while(__check_lock((atomic_p)&a->value, oldval, newval));
    
    return oldval;
}


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *   x)
{
    __clear_lock((atomic_p)x,0);
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *   x)
{
    do
    {
        ;
    }
    while(__check_lock((atomic_p)x, 0, 1));
}


static inline void tMPI_Spinlock_trylock(tMPI_Spinlock_t *   x)
{
    /* Return 0 if we got the lock */
    return (__check_lock((atomic_p)x, 0, 1) != 0)
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *   x)
{
    __clear_lock((atomic_p)x,0);
}


static inline void tMPI_Spinlock_islocked(tMPI_Spinlock_t *   x)
{
    return (*((atomic_p)x) != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *    x)
{
    while(spin_islocked(x)) { ; } 
}


#else
/* No atomic operations, use mutex fallback. Documentation is in x86 section */


#define tMPI_Atomic_memory_barrier()

/* System mutex used for locking to guarantee atomicity */
static tMPI_Thread_mutex_t tMPI_Atomic_mutex = TMPI_THREAD_MUTEX_INITIALIZER;


typedef struct tMPI_Atomic
{
        volatile int value;
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile value;
}
tMPI_Atomic_ptr_t;



#define tMPI_Spinlock_t     tMPI_Thread_mutex_t

 
#  define TMPI_SPINLOCK_INITIALIZER   TMPI_THREAD_MUTEX_INITIALIZER

/* Since mutexes guarantee memory barriers this works fine */
#define tMPI_Atomic_get(a)   ((a)->value)
#define tMPI_Atomic_ptr_get(a)   ((a)->value)


static inline void tMPI_Atomic_set(tMPI_Atomic_t *   a, 
                                  int              i)
{
    /* Mutexes here are necessary to guarantee memory visibility */
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = i;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}

static inline void tMPI_Atomic_ptr_set(tMPI_Atomic_t *   a, 
                                      void*            p)
{
    /* Mutexes here are necessary to guarantee memory visibility */
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = (void*)p;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}



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



#define tMPI_Spinlock_init(lock)       tMPI_Thread_mutex_init(lock)
#define tMPI_Spinlock_lock(lock)       tMPI_Thread_mutex_lock(lock)
#define tMPI_Spinlock_trylock(lock)    tMPI_Thread_mutex_trylock(lock)
#define tMPI_Spinlock_unlock(lock)     tMPI_Thread_mutex_unlock(lock)

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


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    int rc;
    
    tMPI_Spinlock_lock(x);
    /* Got the lock now, so the waiting is over */
    tMPI_Spinlock_unlock(x);
}


#endif




/*! \brief Spinlock-based barrier type
 *
 *  This barrier has the same functionality as the standard
 *  tMPI_Thread_barrier_t, but since it is based on spinlocks
 *  it provides faster synchronization at the cost of busy-waiting.
 *
 *  Variables of this type should be initialized by calling
 *  tMPI_Spinlock_barrier_init() to set the number of threads
 *  that should be synchronized.
 */
typedef struct tMPI_Spinlock_barrier
{
        tMPI_Atomic_t      count;     /*!< Number of threads remaining     */
        int               threshold; /*!< Total number of threads         */
        volatile int      cycle;     /*!< Current cycle (alternating 0/1) */
}
tMPI_Spinlock_barrier_t;
 



/*! \brief Initialize spinlock-based barrier
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




/*! \brief Perform busy-waiting barrier synchronization
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
