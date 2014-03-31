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


/* NOTE:

 ***************************************************************************
   this file is not used any more. gcc intrinsics take care of the atomics
 ***************************************************************************

 */

#error included gcc_ppc.h. This file is outdated



typedef struct tMPI_Atomic
{
    volatile int value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
    void* volatile* value;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


typedef struct tMPI_Spinlock
{
    volatile unsigned int lock;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;
#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK


#define TMPI_SPINLOCK_INITIALIZER   { 0 }

#define tMPI_Atomic_get(a)        ((a)->value)
#define tMPI_Atomic_set(a, i)     (((a)->value) = (i))

#define tMPI_Atomic_ptr_get(a)    (void*)((a)->value)
#define tMPI_Atomic_ptr_set(a, i)  (((a)->value) = (void*)(i))


#if (TMPI_GCC_VERSION >= 40100)

#include "gcc_intrinsics.h"

#else

/* Compiler-dependent stuff: GCC memory barrier */
#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("isync" : : : "memory")



#define TMPI_ATOMIC_HAVE_NATIVE_SWAP
static inline int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    int ret;

    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\tstwcx.  %3,0,%2 \n"
                          "\tbne-    1b\n"
                          : "=&r" (ret), "=m" (a->value)
                          : "r" (&(a->value)), "r" (b)
                          : "cc", "memory");

    return ret;
}

static inline void* tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    int ret;

#if (!defined(__PPC64__)) && (!defined(__ppc64))
    __asm__ __volatile__ ("1:    lwarx   %0,0,%2 \n"
                          "\tstwcx.  %3,0,%2 \n"
                          "\tbne-    1b\n"
                          : "=&r" (ret), "=m" (a->value)
                          : "r" (&(a->value)), "r" (b)
                          : "cc", "memory");
#else
    __asm__ __volatile__ ("1:    ldarx   %0,0,%2 \n"
                          "\tstdcx.  %3,0,%2 \n"
                          "\tbne-    1b\n"
                          : "=&r" (ret), "=m" (a->value)
                          : "r" (&(a->value)), "r" (b)
                          : "cc", "memory");
#endif

    return ret;
}



static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
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

    return prev == oldval;
}


static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t *a, void *oldval,
                                      void *newval)
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
    return prev == oldval;
}

#define TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
{
    int t;

    __asm__ __volatile__("1:     lwarx   %0,0,%2\n"
                         "\tadd     %0,%1,%0\n"
                         "\tstwcx.  %0,0,%2 \n"
                         "\tbne-    1b\n"
                         "\tisync\n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value)
                         : "cc", "memory");
    return t;
}



#define TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
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


#endif


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}



static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
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
                         : "=&r" (tmp)
                         : "r" (&x->lock), "r" (1)
                         : "cr0", "memory");
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    unsigned int           old, t;
    unsigned int           mask = 1;
    volatile unsigned int *p    = &x->lock;

    __asm__ __volatile__("\teieio\n"
                         "1:      lwarx   %0,0,%4 \n"
                         "\tor      %1,%0,%3 \n"
                         "\tstwcx.  %1,0,%4 \n"
                         "\tbne     1b\n"
                         "\tsync\n"
                         : "=&r" (old), "=&r" (t), "=m" (*p)
                         : "r" (mask), "r" (p), "m" (*p)
                         : "cc", "memory");

    return (old & mask);
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x)
{
    __asm__ __volatile__("\teieio\n" : : : "memory");
    x->lock = 0;
}


static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
{
    return ( x->lock != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    do
    {
        tMPI_Atomic_memory_barrier();
    }
    while (tMPI_Spinlock_islocked(x));
}
