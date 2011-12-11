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


/* File contributed by Sergey Klyaus */

#include <atomic.h>

/* this is for newer versions of gcc that have built-in intrinsics,
   on platforms not explicitly supported with inline assembly. */

#define tMPI_Atomic_memory_barrier()  do { membar_consumer(); membar_producer(); } while (0)

/* Only gcc and Intel support this check, otherwise set it to true (skip doc) */
#if (!defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined DOXYGEN)
#define __builtin_constant_p(i) (1)
#endif


typedef struct tMPI_Atomic
{
    volatile uint_t value;   
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
    void* volatile value;   
}
tMPI_Atomic_ptr_t;



#define TMPI_SPINLOCK_INITIALIZER   { 0 }


/* for now we simply assume that int and void* assignments are atomic */
#define tMPI_Atomic_get(a)  ((int)( (a)->value) )
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))


#define tMPI_Atomic_ptr_get(a)  ((void*)((a)->value) )
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))

static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, volatile int i)
{
    return (int) atomic_add_int_nv(&a->value, i);
}

static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, volatile int i)
{
    return (int) atomic_add_int_nv(&a->value, i) - i;
}


static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
    return (int)atomic_cas_uint(&a->value, (uint_t)oldval, (uint_t)newval);
}


static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t* a, void *oldval, 
                                      void *newval)
{
    /*atomic_cas_ptr always returns value stored in a, so*/
    return atomic_cas_ptr(&(a->value), oldval, newval) == oldval;
}



typedef struct tMPI_Spinlock
{
    volatile unsigned long  lock;
} tMPI_Spinlock_t;

#define TMPI_SPINLOCK_INITIALIZER   { 0 }

static inline unsigned long tas(volatile unsigned long *ptr)
{
    unsigned long result;
    __asm__ __volatile__("          \
            ldstub [%1], %0         "
				: "=r"(result)
                : "r"(ptr)
                : "memory");
    return result;
}


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
    do
    {
    } while (tas(&(x->lock)) == 1);
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    return tas(&(x->lock));
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *  x)
{
    x->lock = 0;
}
 
static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *  x)
{
    tMPI_Atomic_memory_barrier();
    return ( x->lock == 1 );
}

static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    do
    {
    } while (x->lock == 1);
    tMPI_Atomic_memory_barrier();
}

