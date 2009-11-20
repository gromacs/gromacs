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
 */

#define tMPI_Atomic_memory_barrier()  { __asm__ __volatile__("\t isync\n"\
                                                             : : : ); }




typedef struct tMPI_Atomic
{
    volatile int value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_t;


typedef struct tMPI_Atomic_ptr
{
    void* volatile *value;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;



typedef struct tMPI_Spinlock
{
    volatile unsigned int lock;  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Spinlock_t;


#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define tMPI_Atomic_get(a)   (int)((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))
#define tMPI_Atomic_ptr_get(a)   (void*)((a)->value) 
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


static inline int tMPI_Atomic_cas(tMPI_Atomic_t *       a,
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

static inline void* tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t *   a,
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




