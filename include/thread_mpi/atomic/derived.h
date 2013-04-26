/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2013, Sander Pronk, Erik Lindahl.
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


/* These functions are fallback definitions for when there are no native
   variants for fetch-add, spinlock, etc., but there is a native
   compare-and-swap. */


/* only define this if there were no separate acquire and release barriers */
#ifndef TMPI_HAVE_ACQ_REL_BARRIERS

/* if they're not defined explicitly, we just make full barriers out of both */
#define tMPI_Atomic_memory_barrier_acq tMPI_Atomic_memory_barrier
#define tMPI_Atomic_memory_barrier_rel tMPI_Atomic_memory_barrier

#endif

#ifndef TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD
TMPI_EXPORT
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
{
    int newval, oldval;
    do
    {
        tMPI_Atomic_memory_barrier_acq();
        oldval = tMPI_Atomic_get(a);
        newval = oldval + i;
    }
    while (!tMPI_Atomic_cas(a, oldval, newval));
    tMPI_Atomic_memory_barrier_rel();
    return oldval;
}
#endif /* TMPI_HAVE_FETCH_ADD */


#ifndef TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN
TMPI_EXPORT
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
{
    /* implement in terms of fetch-add */
    return tMPI_Atomic_fetch_add(a, i) + i;
}
#endif /* TMPI_HAVE_ADD_RETURN */





/* only do this if there was no better solution */
#ifndef TMPI_ATOMIC_HAVE_NATIVE_SWAP
TMPI_EXPORT
static inline int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    int oldval;
    do
    {
        oldval = (int)(a->value);
    }
    while (!tMPI_Atomic_cas(a, oldval, b));
    return oldval;
}


TMPI_EXPORT
static inline void *tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    void *oldval;
    do
    {
        oldval = (void*)(a->value);
    }
    while (!tMPI_Atomic_ptr_cas(a, oldval, b));
    return oldval;
}
#endif



#ifndef TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK

typedef struct tMPI_Spinlock
{
    tMPI_Atomic_t a;
}
tMPI_Spinlock_t;

#define TMPI_SPINLOCK_INITIALIZER   { 0 }



TMPI_EXPORT
static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    tMPI_Atomic_set(&(x->a), 0);
}


TMPI_EXPORT
static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
    tMPI_Atomic_memory_barrier_acq();
    do
    {
        while (tMPI_Atomic_get(&(x->a)) == 1)
        {
            tMPI_Atomic_memory_barrier_acq();
        }
    }
    while (!tMPI_Atomic_cas(&(x->a), 0, 1));
    tMPI_Atomic_memory_barrier_acq();
}


TMPI_EXPORT
static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    int ret;
    tMPI_Atomic_memory_barrier_acq();
    ret = !tMPI_Atomic_cas(&(x->a), 0, 1);
    return ret;
}


TMPI_EXPORT
static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x)
{
    tMPI_Atomic_memory_barrier_rel();
    tMPI_Atomic_set(&(x->a), 0);
    tMPI_Atomic_memory_barrier_rel();
}


TMPI_EXPORT
static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
{
    int ret;
    tMPI_Atomic_memory_barrier_rel();
    ret = (tMPI_Atomic_get(&(x->a)) != 0);
    return ret;
}


TMPI_EXPORT
static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    do
    {
    }
    while (tMPI_Spinlock_islocked(x));
}
#endif /* TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK */
