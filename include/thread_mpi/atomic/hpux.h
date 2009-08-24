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



