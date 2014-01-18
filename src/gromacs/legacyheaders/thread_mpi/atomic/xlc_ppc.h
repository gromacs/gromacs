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

/* PowerPC using xlC inline assembly.
 * Recent versions of xlC (>=7.0) _partially_ support GCC inline assembly
 * if you use the option -qasm=gcc but we have had to hack things a bit, in
 * particular when it comes to clobbered variables. Since this implementation
 * _could_ be buggy, we have separated it from the known-to-be-working gcc
 * one above.
 *
 * For now, we just disable the inline keyword if we're compiling C code:
 */
#if 1
#if (!defined(__cplusplus)) && (!defined(inline))
#define inline_defined_in_atomic 1
#define inline
#endif
#endif


/* IBM xlC compiler */
#ifdef __cplusplus
#include <builtins.h>
#endif


#define TMPI_XLC_INTRINSICS

/* ppc has many memory synchronization instructions */
/*#define tMPI_Atomic_memory_barrier() { __fence(); __sync(); __fence();}*/
/*#define tMPI_Atomic_memory_barrier() __isync();*/
/*#define tMPI_Atomic_memory_barrier() __lwsync();*/

/* for normal memory, this should be enough: */
#define tMPI_Atomic_memory_barrier() { __fence(); __eieio(); __fence(); }
#define tMPI_Atomic_memory_barrier_acq() { __eieio(); __fence(); }
#define tMPI_Atomic_memory_barrier_rel() { __fence(); __eieio(); }
#define TMPI_HAVE_ACQ_REL_BARRIERS

/*#define tMPI_Atomic_memory_barrier() __eieio();*/


typedef struct tMPI_Atomic
{
    volatile int value __attribute__ ((aligned(64)));
}
tMPI_Atomic_t;


typedef struct tMPI_Atomic_ptr
{
    /* volatile char* volatile is not a bug, but means a volatile pointer
       to a volatile value. This is needed for older versions of
       xlc. */
    volatile char* volatile value __attribute__ ((aligned(64)));  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


typedef struct tMPI_Spinlock
{
    volatile int lock __attribute__ ((aligned(64)));
}
tMPI_Spinlock_t;
#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK




#define tMPI_Atomic_get(a)   (int)((a)->value)
#define tMPI_Atomic_set(a, i)  (((a)->value) = (i))
#define tMPI_Atomic_ptr_get(a)   ((a)->value)
#define tMPI_Atomic_ptr_set(a, i)  (((a)->value) = (i))

#define TMPI_SPINLOCK_INITIALIZER   { 0 }


static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
#ifdef TMPI_XLC_INTRINSICS
    int ret;

    __fence(); /* this one needs to be here to avoid ptr. aliasing issues */
    __eieio();
    ret = (__compare_and_swap(&(a->value), &oldval, newval));
    __isync();
    __fence(); /* and this one needs to be here to avoid aliasing issues */
    return ret;
#else
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

    return prev == oldval;
#endif
}


static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t *a, void* oldval,
                                      void* newval)
{
    int                     ret;
    volatile char* volatile oldv = (char*)oldval;
    volatile char* volatile newv = (char*)newval;

    __fence(); /* this one needs to be here to avoid ptr. aliasing issues */
    __eieio();
#if (!defined (__LP64__) ) && (!defined(__powerpc64__) )
    ret = __compare_and_swap((int *)&(a->value), (int*)&oldv, (int)newv);
#else
    ret = __compare_and_swaplp((long *)&(a->value), (long*)&oldv, (long)newv);
#endif
    __isync();
    __fence();

    return ret;
}




static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
{
#ifdef TMPI_XLC_INTRINSICS
    int oldval, newval;

    do
    {
        __fence();
        __eieio(); /* these memory barriers are neccesary */
        oldval = tMPI_Atomic_get(a);
        newval = oldval + i;
    }
    /*while(!__compare_and_swap( &(a->value), &oldval, newval));*/
    while (__check_lock_mp( (int*)&(a->value), oldval, newval));

    /*__isync();*/
    __fence();

    return newval;
#else
    int t;

    __asm__ __volatile__("1:     lwarx   %0,0,%2 \n"
                         "\t add     %0,%1,%0 \n"
                         "\t stwcx.  %0,0,%2 \n"
                         "\t bne-    1b \n"
                         "\t isync \n"
                         : "=&r" (t)
                         : "r" (i), "r" (&a->value) );
    return t;
#endif
}
#define TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN



static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
{
#ifdef TMPI_XLC_INTRINSICS
    int oldval, newval;

    do
    {
        __fence();
        __eieio(); /* these memory barriers are neccesary */
        oldval = tMPI_Atomic_get(a);
        newval = oldval + i;
    }
    /*while(__check_lock_mp((const int*)&(a->value), oldval, newval));*/
    while (__check_lock_mp( (int*)&(a->value), oldval, newval));
    /*while(!__compare_and_swap( &(a->value), &oldval, newval));*/
    /*__isync();*/
    __fence();

    return oldval;
#else
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
#endif
}
#define TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    __fence();
    __clear_lock_mp((const int*)x, 0);
    __fence();
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
    __fence();
    do
    {
    }
    while (__check_lock_mp((int*)&(x->lock), 0, 1));
    tMPI_Atomic_memory_barrier_acq();
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    int ret;
    /* Return 0 if we got the lock */
    __fence();
    ret = __check_lock_mp((int*)&(x->lock), 0, 1);
    tMPI_Atomic_memory_barrier_acq();
    return ret;
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x)
{
    tMPI_Atomic_memory_barrier_rel();
    __clear_lock_mp((int*)&(x->lock), 0);
}


static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
{
    int ret;
    __fence();
    ret = ((x->lock) != 0);
    tMPI_Atomic_memory_barrier_acq();
    return ret;
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    do
    {
    }
    while (tMPI_Spinlock_islocked(x));
}
