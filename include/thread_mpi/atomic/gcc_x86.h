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



#include <limits.h>
#include <stdint.h>
/* This code is executed for x86 and x86-64, with these compilers:
 * GNU
 * Intel 
 * Pathscale
 * All these support GCC-style inline assembly. 
 * We also use this section for the documentation.
 */

#define tMPI_Atomic_memory_barrier() __asm__ __volatile__("": : :"memory")

/* Only gcc and Intel support this check, otherwise set it to true (skip doc) */
#if (!defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined DOXYGEN)
#define __builtin_constant_p(i) (1)
#endif

typedef struct tMPI_Atomic
{
        volatile int       value;   /* Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile*    value;   /* Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;

typedef struct tMPI_Spinlock
{
    volatile unsigned int  lock;   /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;



#define TMPI_SPINLOCK_INITIALIZER   { 1 }



#define tMPI_Atomic_get(a)  ((a)->value) 

#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))

#define tMPI_Atomic_ptr_get(a)  ((a)->value) 

 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))


 
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *     a, 
                                        volatile int       i)
{
    int __i;
    
    __i = i;
    __asm__ __volatile__("lock ; xaddl %0, %1;"
                         :"=r"(i) :"m"(a->value), "0"(i));
    return i + __i;
}  
  

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


static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *   x)
{
    x->lock = 1;
}



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


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *  x)
{
        char old_value;
        
    __asm__ __volatile__("xchgb %b0,%1"
                         :"=q" (old_value), "=m" (x->lock)
                         :"0" (0) : "memory");
    return (old_value <= 0);
}

static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *  x)
{
        char old_value = 1;
        
        __asm__ __volatile__(
                         "xchgb %b0, %1" 
                         :"=q" (old_value), "=m" (x->lock) 
                         :"0" (old_value) : "memory"
                         );
}
 

static inline int tMPI_Spinlock_islocked(tMPI_Spinlock_t *  x)
{
    return (*(volatile signed char *)(&(x)->lock) <= 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    do 
    {
        tMPI_Atomic_memory_barrier(); 
    } 
    while(tMPI_Spinlock_islocked(x));
}

