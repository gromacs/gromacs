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



