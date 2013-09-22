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


/* this is for newer versions of gcc that have built-in intrinsics.
   These are the generic spinlocks:*/


typedef struct tMPI_Spinlock
{
    volatile unsigned int lock /*__attribute__ ((aligned(64)))*/;
} tMPI_Spinlock_t;

#define TMPI_SPINLOCK_INITIALIZER   { 0 }

#define TMPI_ATOMIC_HAVE_NATIVE_SPINLOCK



static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *x)
{
    x->lock = 0;
}


static inline void tMPI_Spinlock_lock(tMPI_Spinlock_t *x)
{
#if 1
    while (__sync_lock_test_and_set(&(x->lock), 1) == 1)
    {
        /* this is nicer on the system bus: */
        while (x->lock == 1)
        {
        }
    }
#else
    do
    {
    }
    while (__sync_lock_test_and_set(&(x->lock), 1) == 1);
#endif
}


static inline int tMPI_Spinlock_trylock(tMPI_Spinlock_t *x)
{
    return __sync_lock_test_and_set(&(x->lock), 1);
}


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *x)
{
    __sync_lock_release(&(x->lock));
}

static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *x)
{
    __sync_synchronize();
    return ( x->lock == 1 );
}

static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    do
    {
    }
    while (x->lock == 1);
    __sync_synchronize();
}
