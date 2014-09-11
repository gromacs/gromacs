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



#include <limits.h>
#include <stdint.h>
/* This code is executed for x86 and x86-64, with these compilers:
 * GNU
 * Intel
 * Pathscale
 * All these support GCC-style inline assembly.
 * We also use this section for the documentation.
 */


#if 0
/* Only gcc and Intel support this check, otherwise set it to true (skip doc) */
#if (!defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined DOXYGEN)
#define __builtin_constant_p(i) (1)
#endif
#endif

/* we put all of these on their own cache line by padding the data structure
   to the size of a cache line on x86 (64 bytes): */
#define TMPI_SIZEOF_X86_CACHE_LINE 64
typedef struct tMPI_Atomic
{
    int  value;
    char padding[TMPI_SIZEOF_X86_CACHE_LINE-sizeof(int)];
} tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
    void* value;
    char  padding[TMPI_SIZEOF_X86_CACHE_LINE-sizeof(void*)];
} tMPI_Atomic_ptr_t;

typedef struct tMPI_Spinlock
{
    unsigned int lock;
    char         padding[TMPI_SIZEOF_X86_CACHE_LINE-sizeof(unsigned int)];
} tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }

#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK



/* these are guaranteed to be  atomic on x86 and x86_64 */
#define tMPI_Atomic_get(a)  ((a)->value)
#define tMPI_Atomic_set(a, i)  (((a)->value) = (i))

#define tMPI_Atomic_ptr_get(a)  ((a)->value)
#define tMPI_Atomic_ptr_set(a, i)  (((a)->value) = (void*)(i))


/* do the intrinsics.

   We disable this for 32-bit builds because the target may be 80386,
   which didn't have cmpxchg, etc (they were introduced as only as 'recently'
   as the 486, and gcc on some Linux versions still target 80386 by default).

   We also specifically check for icc, because intrinsics are not always
   supported there.

   llvm has issues with inline assembly and also in 32 bits has support for
   the gcc intrinsics */
#if ( ( (TMPI_GCC_VERSION >= 40100) && defined(__x86_64__) &&  \
    !defined(__INTEL_COMPILER) )  || defined(__llvm__) )
#include "gcc_intrinsics.h"

#else
/* older versions of gcc don't support atomic intrinsics */

#ifndef __MIC__
#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("sfence;" : : : "memory")
#else
/* MIC is in-order and does not need nor support sfense */
#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("" ::: "memory")
#endif

#define TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
{
    volatile int res = i;
    /* volatile because we read and write back to the same variable in the
       asm section.  some compilers requires this to be volatile */
    __asm__ __volatile__("lock ; xaddl %0, %1;"      /* swap-add */
                         : "=r" (res)                /* with register as
                                                        output*/
                         : "m" (a->value), "0" (res) /* and memory as input */
                         : "memory");
    return res;
}

#define TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
{
    int          orig = i;
    volatile int res  = i;

    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         : "=r" (res)
                         : "m" (a->value), "0" (res)
                         :  "memory");
    return res + orig; /* then add again from the right value */
}



static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
    int prev;

    __asm__ __volatile__("lock ; cmpxchgl %1,%2"
                         : "=a" (prev)
                         : "q" (newval), "m" (a->value), "0" (oldval)
                         : "memory");

    return prev == oldval;
}

static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t *a,
                                      void              *oldval,
                                      void              *newval)
{
    void* prev;
#if (defined(__x86_64__) && !defined(__ILP32__))
    __asm__ __volatile__("lock ; cmpxchgq %1,%2"
                         : "=a" (prev)
                         : "q" (newval), "m" (a->value), "0" (oldval)
                         : "memory");
#elif (defined(__x86_64__) && defined(__ILP32__)) || defined(__i386__)
    __asm__ __volatile__("lock ; cmpxchgl %1,%2"
                         : "=a" (prev)
                         : "q" (newval), "m" (a->value), "0" (oldval)
                         : "memory");
#else
#    error Cannot detect whether this is a 32-bit or 64-bit x86 build.
#endif
    return prev == oldval;
}

#endif /* end of check for gcc intrinsics */


#define TMPI_ATOMIC_HAVE_NATIVE_SWAP
/* do the swap fns; we told the intrinsics that we have them. */
static inline int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    volatile int ret = b;
    __asm__ __volatile__("\txchgl %0, %1;"
                         : "+r" (ret), "+m" (a->value)
                         :
                         : "memory");
    return (int)ret;
}

static inline void *tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    void *volatile *ret = (void* volatile*)b;
#if (defined(__x86_64__) && !defined(__ILP32__))
    __asm__ __volatile__("\txchgq %0, %1;"
                         : "+r" (ret), "+m" (a->value)
                         :
                         : "memory");
#elif (defined(__x86_64__) && defined(__ILP32__)) || defined(__i386__)
    __asm__ __volatile__("\txchgl %0, %1;"
                         : "+r" (ret), "+m" (a->value)
                         :
                         : "memory");
#else
#    error Cannot detect whether this is a 32-bit or 64-bit x86 build.
#endif
    return (void*)ret;
}



/* spinlocks : */

static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}



static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
    /* this is a spinlock with a double loop, as recommended by Intel
       it pauses in the outer loop (the one that just checks for the
       availability of the lock), and thereby reduces bus contention and
       prevents the pipeline from flushing. */
    __asm__ __volatile__("1:\tcmpl $0, %0\n"    /* check the lock */
                         "\tje 2f\n"            /* try to lock if it is
                                                   free by jumping forward */
                         "\tpause\n"            /* otherwise: small pause
                                                   as recommended by Intel */
                         "\tjmp 1b\n"           /* and jump back */

                         "2:\tmovl $1, %%eax\n" /* set eax to 1, the locked
                                                   value of the lock */
                         "\txchgl %%eax, %0\n"  /* atomically exchange
                                                   eax with the lock value */
                         "\tcmpl $0, %%eax\n"   /* compare the exchanged
                                                   value with 0 */
                         "\tjne 1b"             /* jump backward if we didn't
                                                   just lock */
                         : "+m" (x->lock)       /* input & output var */
                         :
                         : "%eax", "memory"     /* we changed memory */
                         );
}



static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x)
{
    /* this is apparently all that is needed for unlocking a lock */
    __asm__ __volatile__(
        "\n\tmovl $0, %0\n"
        : "=m" (x->lock) : : "memory" );
}



static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    int old_value = 1;

    __asm__ __volatile__("\tmovl %2, %0\n"     /* set eax to 1, the locked
                                                  value of the lock */
                         "\txchgl %0, %1\n"    /* atomically exchange
                                                  eax with the address in
                                                  rdx. */
                         : "+r" (old_value), "+m" (x->lock)
                         : "i" (1)
                         : "memory");
    return (old_value);
}



static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
{
    return ( (*((volatile int*)(&(x->lock)))) != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    /* this is the spinlock without the xchg.  */
    __asm__ __volatile__("1:\tcmpl $0, %0\n" /* check the lock */
                         "\tje 2f\n"         /* try to lock if it is
                                                free by jumping forward */
                         "\tpause\n"         /* otherwise: small pause
                                                as recommended by Intel */
                         "\tjmp 1b\n"        /* and jump back */
                         "2:\tnop\n"         /* jump target for end
                                                of wait */
                         : "+m" (x->lock)    /* input & output var */
                         :
                         : "memory"          /* we changed memory */
                         );
#if 0
    do
    {
        tMPI_Atomic_memory_barrier();
    }
    while (tMPI_Spinlock_islocked(x));
#endif
}
