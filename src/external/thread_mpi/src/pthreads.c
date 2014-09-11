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



/* Include the defines that determine which thread library to use.
 * We do not use HAVE_PTHREAD_H directly, since we might want to
 * turn off thread support explicity (e.g. for debugging).
 */
#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#ifdef THREAD_PTHREADS

#ifdef HAVE_PTHREAD_SETAFFINITY
#define _GNU_SOURCE
#endif

/* pthread.h must be the first header, apart from the defines in config.h */
#include <pthread.h>


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "thread_mpi/atomic.h"
#include "thread_mpi/threads.h"
#include "impl.h"
#include "unused.h"

#include "pthreads.h"

/* mutex for initializing mutexes */
static pthread_mutex_t mutex_init = PTHREAD_MUTEX_INITIALIZER;
/* mutex for initializing barriers */
static pthread_mutex_t once_init = PTHREAD_MUTEX_INITIALIZER;
/* mutex for initializing thread_conds */
static pthread_mutex_t cond_init = PTHREAD_MUTEX_INITIALIZER;
/* mutex for initializing barriers */
static pthread_mutex_t barrier_init = PTHREAD_MUTEX_INITIALIZER;

/* mutex for managing  thread IDs */
static pthread_mutex_t thread_id_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_key_t   thread_id_key;
static int             thread_id_key_initialized = 0;




enum tMPI_Thread_support tMPI_Thread_support(void)
{
    return TMPI_THREAD_SUPPORT_YES;
}


int tMPI_Thread_get_hw_number(void)
{
    int ret = 0;
#ifdef HAVE_SYSCONF
#if defined(_SC_NPROCESSORS_ONLN)
    ret = sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    ret = sysconf(_SC_NPROC_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    ret = sysconf(_SC_NPROCESSORS_CONF);
#elif defined(_SC_NPROC_CONF)
    ret = sysconf(_SC_NPROC_CONF);
#endif
#endif

    return ret;
}

/* destructor for thread ids */
static void tMPI_Destroy_thread_id(void* thread_id)
{
    struct tMPI_Thread *thread = (struct tMPI_Thread*)thread_id;
    if (!thread->started_by_tmpi)
    {
        /* if the thread is started by tMPI, it must be freed in the join()
           call. */
        free(thread_id);
    }
}


/* Set the thread id var for this thread
    Returns a pointer to the thread object if succesful, NULL if ENOMEM */
static struct tMPI_Thread* tMPI_Set_thread_id_key(int started_by_tmpi)
{
    struct tMPI_Thread *th;

    th = (struct tMPI_Thread*)malloc(sizeof(struct tMPI_Thread)*1);
    if (th == NULL)
    {
        return NULL;
    }
    th->th              = pthread_self();
    th->started_by_tmpi = started_by_tmpi;
    /* we ignore errors because any thread that needs this value will
       re-generate it in the next iteration. */
    pthread_setspecific(thread_id_key, th);
    return th;
}

/* initialize the thread id vars if not already initialized */
static int tMPI_Init_thread_ids(void)
{
    int ret;
    ret = pthread_mutex_lock( &thread_id_mutex );
    if (ret != 0)
    {
        return ret;
    }

    if (!thread_id_key_initialized)
    {
        /* initialize and set the thread id thread-specific variable */
        struct tMPI_Thread *th;

        thread_id_key_initialized = 1;
        ret = pthread_key_create(&thread_id_key, tMPI_Destroy_thread_id);
        if (ret != 0)
        {
            goto err;
        }
        th = tMPI_Set_thread_id_key(0);
        if (th == NULL)
        {
            ret = ENOMEM;
            goto err;
        }
    }

    ret = pthread_mutex_unlock( &thread_id_mutex );
    return ret;
err:
    pthread_mutex_unlock( &thread_id_mutex );
    return ret;
}

/* structure to hold the arguments for the thread_starter function */
struct tMPI_Thread_starter
{
    struct tMPI_Thread *thread;
    void               *(*start_routine)(void*);
    void               *arg;
    pthread_mutex_t     cond_lock; /* lock for initialization of thread
                                      structure */
};

/* the thread_starter function that sets the thread id */
#ifdef __MINGW32__
__attribute__((force_align_arg_pointer))
#endif
static void *tMPI_Thread_starter(void *arg)
{
    struct tMPI_Thread_starter *starter = (struct tMPI_Thread_starter *)arg;
    void *(*start_routine)(void*);
    void *parg;
    int   ret;

    /* first wait for the parent thread to signal that the starter->thread
       structure is ready. That's done by unlocking the starter->cond_lock */
    ret = pthread_mutex_lock(&(starter->cond_lock));
    if (ret != 0)
    {
        return NULL;
    }
    ret = pthread_mutex_unlock(&(starter->cond_lock));
    if (ret != 0)
    {
        return NULL;
    }

    /* now remember the tMPI_thread_t structure for this thread */
    ret = pthread_setspecific(thread_id_key, starter->thread);
    if (ret != 0)
    {
        return NULL;
    }
    start_routine = starter->start_routine;
    parg          = starter->arg;

    /* deallocate the starter structure. Errors here are non-fatal. */
    pthread_mutex_destroy(&(starter->cond_lock));
    free(starter);
    return (*start_routine)(parg);
}

int tMPI_Thread_create(tMPI_Thread_t *thread, void *(*start_routine)(void *),
                       void *arg)
{
    int ret;
    struct tMPI_Thread_starter *starter;

    if (thread == NULL)
    {
        return EINVAL;
    }
    tMPI_Init_thread_ids();

    *thread = (struct tMPI_Thread*)malloc(sizeof(struct tMPI_Thread)*1);
    if (*thread == NULL)
    {
        return ENOMEM;
    }
    (*thread)->started_by_tmpi = 1;
    starter                    = (struct tMPI_Thread_starter*)
        malloc(sizeof(struct tMPI_Thread_starter)*1);
    if (starter == NULL)
    {
        return ENOMEM;
    }
    /* fill the starter structure */
    starter->thread        = *thread;
    starter->start_routine = start_routine;
    starter->arg           = arg;

    ret = pthread_mutex_init(&(starter->cond_lock), NULL);
    if (ret != 0)
    {
        return ret;
    }
    /* now lock the mutex so we can unlock it once we know the data in
       thread->th is safe. */
    ret = pthread_mutex_lock(&(starter->cond_lock));
    if (ret != 0)
    {
        return ret;
    }

    ret = pthread_create(&((*thread)->th), NULL, tMPI_Thread_starter,
                         (void*)starter);
    if (ret != 0)
    {
        return ret;
    }

    /* Here we know thread->th is safe. */
    ret = pthread_mutex_unlock(&(starter->cond_lock));

    return ret;
}



int tMPI_Thread_join(tMPI_Thread_t thread, void **value_ptr)
{
    int       ret;
    pthread_t th = thread->th;

    ret = pthread_join( th, value_ptr );
    if (ret != 0)
    {
        return ret;
    }
    free(thread);
    return 0;
}


tMPI_Thread_t tMPI_Thread_self(void)
{
    tMPI_Thread_t th;
    int           ret;

    /* make sure the key var is set */
    ret = tMPI_Init_thread_ids();
    if (ret != 0)
    {
        return NULL;
    }

    th = pthread_getspecific(thread_id_key);

    /* check if it is already in our list */
    if (th == NULL)
    {
        th = tMPI_Set_thread_id_key(0);
    }
    return th;
}

int tMPI_Thread_equal(tMPI_Thread_t t1, tMPI_Thread_t t2)
{
    return pthread_equal(t1->th, t2->th);
}


enum tMPI_Thread_setaffinity_support tMPI_Thread_setaffinity_support(void)
{
#ifdef HAVE_PTHREAD_SETAFFINITY
    cpu_set_t set;
    int       ret;

    /* run getaffinity to check whether we get back ENOSYS */
    ret = pthread_getaffinity_np(pthread_self(), sizeof(set), &set);
    if (ret == 0)
    {
        return TMPI_SETAFFINITY_SUPPORT_YES;
    }
    else
    {
        return TMPI_SETAFFINITY_SUPPORT_NO;
    }
#else
    return TMPI_SETAFFINITY_SUPPORT_NO;
#endif
}


/* set thread's own affinity to a processor number n */
int tMPI_Thread_setaffinity_single(tMPI_Thread_t tmpi_unused thread,
                                   unsigned int  tmpi_unused nr)
{
#ifdef HAVE_PTHREAD_SETAFFINITY
    int       nt = tMPI_Thread_get_hw_number();
    cpu_set_t set;

    if (nt < nr)
    {
        return TMPI_ERR_PROCNR;
    }

    CPU_ZERO(&set);
    CPU_SET(nr, &set);
    return pthread_setaffinity_np(thread->th, sizeof(set), &set);
#else
    return 0;
#endif
}




int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    if (mtx == NULL)
    {
        return EINVAL;
    }

    mtx->mutex = (struct tMPI_Mutex*)malloc(sizeof(struct tMPI_Mutex)*1);
    if (mtx->mutex == NULL)
    {
        return ENOMEM;
    }
    ret = pthread_mutex_init(&(mtx->mutex->mtx), NULL);
    if (ret != 0)
    {
        return ret;
    }

#ifndef TMPI_NO_ATOMICS
    tMPI_Atomic_set(&(mtx->initialized), 1);
#else
    mtx->initialized.value = 1;
#endif
    return 0;
}

static inline int tMPI_Thread_mutex_init_once(tMPI_Thread_mutex_t *mtx)
{
    int ret = 0;

#ifndef TMPI_NO_ATOMICS
    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
#endif
    {
        /* we're relying on the memory barrier semantics of mutex_lock/unlock
           for the check preceding this function call to have worked */
        ret = pthread_mutex_lock( &(mutex_init) );
        if (ret != 0)
        {
            return ret;
        }

        if (mtx->mutex == NULL)
        {
            mtx->mutex = (struct tMPI_Mutex*)malloc(sizeof(struct tMPI_Mutex));
            if (mtx->mutex == NULL)
            {
                ret = ENOMEM;
                goto err;
            }
            ret = pthread_mutex_init( &(mtx->mutex->mtx), NULL);
            if (ret != 0)
            {
                goto err;
            }
        }
        ret = pthread_mutex_unlock( &(mutex_init) );
    }
    return ret;
err:
    pthread_mutex_unlock( &(mutex_init) );
    return ret;
}


int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    if (mtx == NULL)
    {
        return EINVAL;
    }

    ret = pthread_mutex_destroy( &(mtx->mutex->mtx) );
    if (ret != 0)
    {
        return ret;
    }
    free(mtx->mutex);
    return ret;
}



int tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the mutex is initialized */
    ret = tMPI_Thread_mutex_init_once(mtx);
    if (ret != 0)
    {
        return ret;
    }

    ret = pthread_mutex_lock(&(mtx->mutex->mtx));
    return ret;
}




int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the mutex is initialized */
    ret = tMPI_Thread_mutex_init_once(mtx);
    if (ret != 0)
    {
        return ret;
    }

    ret = pthread_mutex_trylock(&(mtx->mutex->mtx));
    return ret;
}



int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the mutex is initialized */
    ret = tMPI_Thread_mutex_init_once(mtx);
    if (ret != 0)
    {
        return ret;
    }

    ret = pthread_mutex_unlock(&(mtx->mutex->mtx));
    return ret;
}



int tMPI_Thread_key_create(tMPI_Thread_key_t *key, void (*destructor)(void *))
{
    int ret;

    if (key == NULL)
    {
        return EINVAL;
    }


    key->key = (struct tMPI_Thread_key*)malloc(sizeof(struct
                                                      tMPI_Thread_key)*1);
    if (key->key == NULL)
    {
        return ENOMEM;
    }
    ret = pthread_key_create(&((key)->key->pkey), destructor);
    if (ret != 0)
    {
        return ret;
    }

    tMPI_Atomic_set(&(key->initialized), 1);
    return 0;
}


int tMPI_Thread_key_delete(tMPI_Thread_key_t key)
{
    int ret;

    ret = pthread_key_delete((key.key->pkey));
    if (ret != 0)
    {
        return ret;
    }
    free(key.key);

    return 0;
}



void * tMPI_Thread_getspecific(tMPI_Thread_key_t key)
{
    void *p = NULL;

    p = pthread_getspecific((key.key->pkey));

    return p;
}


int tMPI_Thread_setspecific(tMPI_Thread_key_t key, void *value)
{
    int ret;

    ret = pthread_setspecific((key.key->pkey), value);

    return ret;
}


int tMPI_Thread_once(tMPI_Thread_once_t *once_control,
                     void                (*init_routine)(void))
{
    int ret;
    if (!once_control || !init_routine)
    {
        return EINVAL;
    }

    /* really ugly hack - and it's slow... */
    ret = pthread_mutex_lock( &once_init );
    if (ret != 0)
    {
        return ret;
    }
    if (tMPI_Atomic_get(&(once_control->once)) == 0)
    {
        (*init_routine)();
        tMPI_Atomic_set(&(once_control->once), 1);
    }
    ret = pthread_mutex_unlock( &once_init );

    return ret;
}




int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond)
{
    int ret;

    if (cond == NULL)
    {
        return EINVAL;
    }

    cond->condp = (struct tMPI_Thread_cond*)malloc(
                sizeof(struct tMPI_Thread_cond));
    if (cond->condp == NULL)
    {
        return ENOMEM;
    }

    ret = pthread_cond_init(&(cond->condp->cond), NULL);
    if (ret != 0)
    {
        return ret;
    }
    tMPI_Atomic_set(&(cond->initialized), 1);
    tMPI_Atomic_memory_barrier();

    return 0;
}


static int tMPI_Thread_cond_init_once(tMPI_Thread_cond_t *cond)
{
    int ret = 0;

    /* we're relying on the memory barrier semantics of mutex_lock/unlock
       for the check preceding this function call to have worked */
    ret = pthread_mutex_lock( &(cond_init) );
    if (ret != 0)
    {
        return ret;
    }
    if (cond->condp == NULL)
    {
        cond->condp = (struct tMPI_Thread_cond*)
            malloc(sizeof(struct tMPI_Thread_cond)*1);
        if (cond->condp == NULL)
        {
            ret = ENOMEM;
            goto err;
        }
        ret = pthread_cond_init( &(cond->condp->cond), NULL);
        if (ret != 0)
        {
            goto err;
        }
    }
    ret = pthread_mutex_unlock( &(cond_init) );
    return ret;
err:
    /* try to unlock anyway */
    pthread_mutex_unlock( &(cond_init) );
    return ret;
}



int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond)
{
    int ret;

    if (cond == NULL)
    {
        return EINVAL;
    }

    ret = pthread_cond_destroy(&(cond->condp->cond));
    if (ret != 0)
    {
        return ret;
    }
    free(cond->condp);

    return 0;
}


int tMPI_Thread_cond_wait(tMPI_Thread_cond_t *cond, tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        ret = tMPI_Thread_cond_init_once(cond);
        if (ret != 0)
        {
            return ret;
        }
    }
    /* the mutex must have been initialized because it should be locked here */

    ret = pthread_cond_wait( &(cond->condp->cond), &(mtx->mutex->mtx) );

    return ret;
}




int tMPI_Thread_cond_signal(tMPI_Thread_cond_t *cond)
{
    int ret;

    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        ret = tMPI_Thread_cond_init_once(cond);
        if (ret != 0)
        {
            return ret;
        }
    }

    ret = pthread_cond_signal( &(cond->condp->cond) );

    return ret;
}



int tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *cond)
{
    int ret;

    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        ret = tMPI_Thread_cond_init_once(cond);
        if (ret != 0)
        {
            return ret;
        }
    }

    ret = pthread_cond_broadcast( &(cond->condp->cond) );

    return ret;
}




void tMPI_Thread_exit(void *value_ptr)
{
    pthread_exit(value_ptr);
}


int tMPI_Thread_cancel(tMPI_Thread_t thread)
{
    #ifdef __native_client__
    return ENOSYS;
    #endif
    return pthread_cancel(thread->th);
}




int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier, int n)
{
    int ret;
    /*tMPI_Thread_pthread_barrier_t *p;*/

    if (barrier == NULL)
    {
        return EINVAL;
    }

    barrier->barrierp = (struct tMPI_Thread_barrier*)
        malloc(sizeof(struct tMPI_Thread_barrier)*1);
    if (barrier->barrierp == NULL)
    {
        return ENOMEM;
    }

    ret = pthread_mutex_init(&(barrier->barrierp->mutex), NULL);
    if (ret != 0)
    {
        return ret;
    }

    ret = pthread_cond_init(&(barrier->barrierp->cv), NULL);
    if (ret != 0)
    {
        return ret;
    }

    barrier->threshold = n;
    barrier->count     = n;
    barrier->cycle     = 0;

    tMPI_Atomic_set(&(barrier->initialized), 1);
    return 0;
}

static int tMPI_Thread_barrier_init_once(tMPI_Thread_barrier_t *barrier)
{
    int ret = 0;

    /* we're relying on the memory barrier semantics of mutex_lock/unlock
       for the check preceding this function call to have worked */
    ret = pthread_mutex_lock( &(barrier_init) );
    if (ret != 0)
    {
        return ret;
    }

    if (barrier->barrierp == NULL)
    {
        barrier->barrierp = (struct tMPI_Thread_barrier*)
            malloc(sizeof(struct tMPI_Thread_barrier)*1);
        if (barrier->barrierp == NULL)
        {
            ret = ENOMEM;
            goto err;
        }

        ret = pthread_mutex_init(&(barrier->barrierp->mutex), NULL);

        if (ret != 0)
        {
            goto err;
        }

        ret = pthread_cond_init(&(barrier->barrierp->cv), NULL);

        if (ret != 0)
        {
            goto err;
        }
    }
    ret = pthread_mutex_unlock( &(barrier_init) );
    return ret;
err:
    pthread_mutex_unlock( &(barrier_init) );
    return ret;
}




int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier)
{
    int ret;

    if (barrier == NULL)
    {
        return EINVAL;
    }

    ret = pthread_mutex_destroy(&(barrier->barrierp->mutex));
    if (ret != 0)
    {
        return ret;
    }
    ret = pthread_cond_destroy(&(barrier->barrierp->cv));
    if (ret != 0)
    {
        return ret;
    }

    free(barrier->barrierp);

    return 0;
}


int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t * barrier)
{
    int cycle;
    int ret;

    /* check whether the barrier is initialized */
    if (tMPI_Atomic_get( &(barrier->initialized)  ) == 0)
    {
        tMPI_Thread_barrier_init_once(barrier);
    }


    ret = pthread_mutex_lock(&barrier->barrierp->mutex);
    if (ret != 0)
    {
        return ret;
    }

    cycle = barrier->cycle;

    /* Decrement the count atomically and check if it is zero.
     * This will only be true for the last thread calling us.
     */
    if (--barrier->count <= 0)
    {
        barrier->cycle = !barrier->cycle;
        barrier->count = barrier->threshold;
        ret            = pthread_cond_broadcast(&barrier->barrierp->cv);

        if (ret == 0)
        {
            goto err;
        }
    }
    else
    {
        while (cycle == barrier->cycle)
        {
            ret = pthread_cond_wait(&barrier->barrierp->cv,
                                    &barrier->barrierp->mutex);
            if (ret != 0)
            {
                goto err;
            }
        }
    }

    ret = pthread_mutex_unlock(&barrier->barrierp->mutex);
    return ret;
err:
    pthread_mutex_unlock(&barrier->barrierp->mutex);
    return ret;

}

#else

/* just to have some symbols */
int tMPI_Thread_pthreads = 0;

#endif /* THREAD_PTHREADS */
