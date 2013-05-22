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

#include "impl.h"

/* This file is only needed when no intrinsic atomic operations are present. */
#ifdef TMPI_NO_ATOMICS

/** System mutex used for locking to guarantee atomicity */
static tMPI_Thread_mutex_t tMPI_Atomic_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

struct tMPI_Spinlock
{
    tMPI_Thread_mutex_t *lock;
};

int tMPI_Atomic_get(const tMPI_Atomic_t *a)
{
    int ret;
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    ret = a->value;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return ret;
}

void tMPI_Atomic_set(tMPI_Atomic_t *a, int value)
{
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = value;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}


void* tMPI_Atomic_ptr_get(const tMPI_Atomic_ptr_t *a)
{
    void* ret;
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    ret = a->value;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return ret;
}

void tMPI_Atomic_ptr_set(tMPI_Atomic_ptr_t *a, void *value)
{
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    a->value = value;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}

int tMPI_Atomic_add_return(tMPI_Atomic_t *a, int i)
{
    int t;
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    t        = a->value + i;
    a->value = t;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}

int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, int i)
{
    int old_value;

    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    old_value = a->value;
    a->value  = old_value + i;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return old_value;
}

int tMPI_Atomic_cas(tMPI_Atomic_t *a, int old_val, int new_val)
{
    int t = 0;

    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    if (a->value == old_val)
    {
        a->value = new_val;
        t        = 1;
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}


int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t * a, void *old_val, void *new_val)
{
    int t = 0;

    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    if (a->value == old_val)
    {
        a->value = new_val;
        t        = 1;
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
    return t;
}

int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    int ret;
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    ret      = a->value;
    a->value = b;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);

    return ret;
}

void *tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    void *ret;

    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    ret      = a->value;
    a->value = b;
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);

    return ret;
}


void tMPI_Spinlock_init( tMPI_Spinlock_t *x)
{
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    *x         = (tMPI_Spinlock_t)malloc(sizeof(tMPI_Spinlock_t));
    (*x)->lock = (tMPI_Thread_mutex_t*)malloc(sizeof(tMPI_Thread_mutex_t));
    tMPI_Thread_mutex_init((*x)->lock);
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}

/* NOTE: assumes atomic mutex is locked */
static void tMPI_Spinlock_init_once(tMPI_Spinlock_t *x)
{
    tMPI_Thread_mutex_lock(&tMPI_Atomic_mutex);
    if (!*x)
    {
        *x         = (tMPI_Spinlock_t)malloc(sizeof(tMPI_Spinlock_t));
        (*x)->lock = (tMPI_Thread_mutex_t*)malloc(sizeof(tMPI_Thread_mutex_t));
        tMPI_Thread_mutex_init((*x)->lock);
    }
    tMPI_Thread_mutex_unlock(&tMPI_Atomic_mutex);
}


void tMPI_Spinlock_lock( tMPI_Spinlock_t *x)
{
    tMPI_Spinlock_init_once(x);
    tMPI_Thread_mutex_lock((*x)->lock);
}

void tMPI_Spinlock_unlock( tMPI_Spinlock_t *x)
{
    tMPI_Spinlock_init_once(x);
    tMPI_Thread_mutex_unlock((*x)->lock);
}

int tMPI_Spinlock_trylock( tMPI_Spinlock_t *x)
{
    int ret;
    tMPI_Spinlock_init_once(x);
    ret = tMPI_Thread_mutex_trylock((*x)->lock);
    return ret;
}

int tMPI_Spinlock_islocked(tMPI_Spinlock_t *x)
{
    int ret;
    tMPI_Spinlock_init_once(x);
    ret = tMPI_Thread_mutex_trylock((*x)->lock);
    if (ret == 0)
    {
        tMPI_Thread_mutex_unlock((*x)->lock);
        ret = 0;
    }
    else
    {
        ret = 1;
    }

    return ret;
}


void tMPI_Spinlock_wait(tMPI_Spinlock_t *x)
{
    tMPI_Spinlock_init_once(x);

    tMPI_Spinlock_lock(x);
    /* Got the lock now, so the waiting is over */
    tMPI_Spinlock_unlock(x);
}

#else

/* just to have some symbols */
int _tMPI_Atomics = 1;

#endif
