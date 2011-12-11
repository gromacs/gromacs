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


#include "pthreads.h"


/* mutex for initializing mutexes */
static pthread_mutex_t mutex_init=PTHREAD_MUTEX_INITIALIZER;   
/* mutex for initializing barriers */
static pthread_mutex_t once_init=PTHREAD_MUTEX_INITIALIZER;    
/* mutex for initializing thread_conds */
static pthread_mutex_t cond_init=PTHREAD_MUTEX_INITIALIZER;    
/* mutex for initializing barriers */
static pthread_mutex_t barrier_init=PTHREAD_MUTEX_INITIALIZER; 


/* lock & variable for setting main thread affinity */
static tMPI_Spinlock_t main_thread_aff_lock=TMPI_SPINLOCK_INITIALIZER;
static tMPI_Atomic_t main_thread_aff_set={ 0 };
static tMPI_Atomic_t aff_thread_number;




/* TODO: this needs to go away!  (there's another one in winthreads.c)
   fatal errors are thankfully really rare*/
void tMPI_Fatal_error(const char *file, int line, const char *message, ...)
{
    va_list ap;

    fprintf(stderr, "tMPI Fatal error in %s, line %d: ", file, line);
    va_start(ap, message);
    vfprintf(stderr, message, ap);
    va_end(ap);
    fprintf(stderr,"\n");
    abort();
}


enum tMPI_Thread_support tMPI_Thread_support(void)
{
    return TMPI_THREAD_SUPPORT_YES;
}


int tMPI_Thread_get_hw_number(void)
{
    int ret=0;
#ifdef HAVE_SYSCONF
#if defined(_SC_NPROCESSORS_ONLN)
    ret=sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
    ret=sysconf(_SC_NPROC_ONLN);
#elif defined(_SC_NPROCESSORS_CONF)
    ret=sysconf(_SC_NPROCESSORS_CONF);
#elif defined(_SC_NPROC_CONF)
    ret=sysconf(_SC_NPROC_CONF);
#endif
#endif

    return ret;
}

int tMPI_Thread_create(tMPI_Thread_t *thread, void *(*start_routine)(void *),
                       void *arg)
{
    int ret;

    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }

    *thread=(struct tMPI_Thread*)malloc(sizeof(struct tMPI_Thread)*1);
    ret=pthread_create(&((*thread)->th),NULL,start_routine,arg);

    if(ret!=0)
    {
        /* Cannot use tMPI_error() since messages use threads for locking */
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create POSIX thread, rc=%d",ret);
        /* Use system memory allocation routines */
        return -1;
    }

    return 0;
}


/* set thread's own affinity to a processor number n */
static int tMPI_Set_affinity(int n)
{
#ifdef HAVE_PTHREAD_SETAFFINITY
    cpu_set_t set;

    CPU_ZERO(&set);
    CPU_SET(n, &set);
    return pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
#endif
    return 0;
}

int tMPI_Thread_create_aff(tMPI_Thread_t *thread, 
                           void *(*start_routine)(void *),
                           void *arg)
{
    int ret;

    /* set the calling thread's affinity mask */
    if (tMPI_Atomic_get(&main_thread_aff_set) == 0)
    {
#ifdef HAVE_PTHREAD_SETAFFINITY
        cpu_set_t set;
#endif
        /* this can be a spinlock because the chances of collision are low. */
        tMPI_Spinlock_lock( &main_thread_aff_lock );
        tMPI_Atomic_set( &aff_thread_number, 0);
#ifdef HAVE_PTHREAD_SETAFFINITY
        CPU_ZERO(&set);
        CPU_SET(0, &set);
        pthread_setaffinity_np(pthread_self(), sizeof(set), &set);
        /*fprintf(stderr, "Setting affinity.\n");*/
#endif
        tMPI_Atomic_set( &main_thread_aff_set, 1);
        tMPI_Spinlock_unlock( &main_thread_aff_lock );
    }


    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }

    *thread=(struct tMPI_Thread*)malloc(sizeof(struct tMPI_Thread)*1);
    ret=pthread_create(&((*thread)->th),NULL,start_routine,arg);

    if(ret!=0)
    {
        /* Cannot use tMPI_error() since messages use threads for locking */
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create POSIX thread, rc=%d",ret);
        /* Use system memory allocation routines */
        return -1;
    }
    else
    {
#ifdef HAVE_PTHREAD_SETAFFINITY
        int n;
        cpu_set_t set;

        n=tMPI_Atomic_add_return(&aff_thread_number, 1);
        CPU_ZERO(&set);
        CPU_SET(n, &set);
        return pthread_setaffinity_np((*thread)->th, sizeof(set), &set);
#else
        return 0;
#endif
    }
}





int tMPI_Thread_join(tMPI_Thread_t thread, void **value_ptr)
{
    int ret;
    pthread_t th=thread->th;

    free(thread);
    
    ret = pthread_join( th, value_ptr );

    if(ret != 0 )
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to join POSIX thread. rc=%d",ret);
    }
    return ret;
}



int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx) 
{
    int ret;
  
    if (mtx == NULL)
    {
        return EINVAL;
    }

    mtx->mutex=(struct tMPI_Mutex*)tMPI_Malloc(sizeof(struct tMPI_Mutex)*1);
    ret = pthread_mutex_init(&(mtx->mutex->mtx),NULL);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX mutex. rc=%d");
        /* Use system memory allocation routines */
        return ret;
    }

    tMPI_Atomic_set(&(mtx->initialized), 1);
    return 0;
}

static int tMPI_Thread_mutex_init_once(tMPI_Thread_mutex_t *mtx)
{
    int ret=0;

    /* we're relying on the memory barrier semantics of mutex_lock/unlock
       for the check preceding this function call to have worked */
    pthread_mutex_lock( &(mutex_init) );    
    if(mtx->mutex==NULL)
    {
        mtx->mutex=(struct tMPI_Mutex*)tMPI_Malloc(sizeof(struct tMPI_Mutex)*1);
        ret=pthread_mutex_init( &(mtx->mutex->mtx), NULL);
    }
    pthread_mutex_unlock( &(mutex_init) );    
    return ret;
}


int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx) 
{
    int ret;

    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_mutex_destroy( &(mtx->mutex->mtx) );
    free(mtx->mutex);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error destroying POSIX mutex. rc=%d",ret);
        /* Use system memory allocation routines */
    }
    return ret;
}



int tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
    {
        ret=tMPI_Thread_mutex_init_once(mtx);
        if (ret)
            return ret;
    }
   
    ret=pthread_mutex_lock(&(mtx->mutex->mtx));

    return ret;
}

 


int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
    {
        ret=tMPI_Thread_mutex_init_once(mtx);
        if (ret)
            return ret;
    }

    ret=pthread_mutex_trylock(&(mtx->mutex->mtx));
    
    return ret;
}



int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{
    int ret;
 
    /* check whether the mutex is initialized */
    if (tMPI_Atomic_get( &(mtx->initialized)  ) == 0)
    {
        ret=tMPI_Thread_mutex_init_once(mtx);
        if (ret)
            return ret;
    }
 
    ret = pthread_mutex_unlock(&(mtx->mutex->mtx));
    
    return ret;
}



int tMPI_Thread_key_create(tMPI_Thread_key_t *key, void (*destructor)(void *))
{
    int ret;

    if(key==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid key pointer.");
        return EINVAL;
    }


    key->key=(struct tMPI_Thread_key*)tMPI_Malloc(sizeof(struct 
                                                         tMPI_Thread_key)*1);
    ret = pthread_key_create(&((key)->key->pkey), destructor);
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create thread key, rc=%d.",ret);
        fflush(stderr);
        return -1;
    }

    tMPI_Atomic_set(&(key->initialized), 1);
    return 0;
}


int tMPI_Thread_key_delete(tMPI_Thread_key_t key)
{
    int ret;

    ret=pthread_key_delete((key.key->pkey));
    free(key.key);

    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to delete thread key, rc=%d.",ret);
        fflush(stderr);
    }
    
    return ret;
}



void * tMPI_Thread_getspecific(tMPI_Thread_key_t key)
{
    void *p = NULL;

    p=pthread_getspecific((key.key->pkey));

    return p;
}


int tMPI_Thread_setspecific(tMPI_Thread_key_t key, void *value)
{
    int ret;
    
    ret=pthread_setspecific((key.key->pkey),value);
    
    return ret;
}


int tMPI_Thread_once(tMPI_Thread_once_t *once_control,
                     void (*init_routine)(void))
{
    int ret;
    if (!once_control || !init_routine)
    {
        return EINVAL;
    }

    /* really ugly hack - and it's slow... */
    if ( (ret=pthread_mutex_lock( &once_init )) )
        return ret;
    if (tMPI_Atomic_get(&(once_control->once)) == 0)
    {
        (*init_routine)();
        tMPI_Atomic_set(&(once_control->once), 1);
    }
    pthread_mutex_unlock( &once_init );

    return 0;
}




int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond) 
{
    int ret;
    
    if(cond==NULL)
    {
        return EINVAL;
    }
   
    cond->condp=(struct tMPI_Thread_cond*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_cond)*1);
    ret = pthread_cond_init(&(cond->condp->cond), NULL);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX condition variable. rc=%d",ret);
        fflush(stderr);
    }
    tMPI_Atomic_set(&(cond->initialized),1);
    return ret;
}


static int tMPI_Thread_cond_init_once(tMPI_Thread_cond_t *cond)
{
    int ret=0;

    /* we're relying on the memory barrier semantics of mutex_lock/unlock
       for the check preceding this function call to have worked */
    pthread_mutex_lock( &(cond_init) );    
    if(cond->condp==NULL)
    {
        cond->condp=(struct tMPI_Thread_cond*)
                  tMPI_Malloc(sizeof(struct tMPI_Thread_cond)*1);
        ret=pthread_cond_init( &(cond->condp->cond), NULL);
    }
    pthread_mutex_unlock( &(cond_init) );    
    return ret;
}



int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond) 
{
    int ret;
    
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_cond_destroy(&(cond->condp->cond));
    free(cond->condp);
   
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,
                         "Error destroying POSIX condition variable. rc=%d",
                         ret);
        fflush(stderr);
    }
    return ret;
}


int tMPI_Thread_cond_wait(tMPI_Thread_cond_t *cond, tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* check whether the condition is initialized */
    if (tMPI_Atomic_get( &(cond->initialized)  ) == 0)
    {
        tMPI_Thread_cond_init_once(cond);
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
        tMPI_Thread_cond_init_once(cond);
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
        tMPI_Thread_cond_init_once(cond);
    }
   
    ret = pthread_cond_broadcast( &(cond->condp->cond) );
    
    return ret;
}




void tMPI_Thread_exit(void *      value_ptr)
{
    pthread_exit(value_ptr);
}


int tMPI_Thread_cancel(tMPI_Thread_t     thread)
{
    return pthread_cancel(thread->th);
}




int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier, int n)
{
    int ret;
    /*tMPI_Thread_pthread_barrier_t *p;*/
    
    if(barrier==NULL)
    {
        return EINVAL;
    }
    
    barrier->barrierp=(struct tMPI_Thread_barrier*)
              tMPI_Malloc(sizeof(struct tMPI_Thread_barrier)*1);
    ret = pthread_mutex_init(&(barrier->barrierp->mutex),NULL);
        
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX mutex. rc=%d",
                         ret);
        return ret;
    }
    
    ret = pthread_cond_init(&(barrier->barrierp->cv),NULL);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,
                         "Error initializing POSIX condition variable. rc=%d",
                         ret);
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
    int ret=0;

    /* we're relying on the memory barrier semantics of mutex_lock/unlock
       for the check preceding this function call to have worked */
    pthread_mutex_lock( &(barrier_init) );    
    if(barrier->barrierp==NULL)
    {
        barrier->barrierp=(struct tMPI_Thread_barrier*)
                  tMPI_Malloc(sizeof(struct tMPI_Thread_barrier)*1);
        ret = pthread_mutex_init(&(barrier->barrierp->mutex),NULL);

        if(ret!=0)
        {
            tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX mutex. rc=%d",
                             ret);
            return ret;
        }

        ret = pthread_cond_init(&(barrier->barrierp->cv),NULL);

        if(ret!=0)
        {
            tMPI_Fatal_error(TMPI_FARGS,
                             "Error initializing POSIX condition variable. rc=%d",
                             ret);
            return ret;
        }
    }
    pthread_mutex_unlock( &(barrier_init) );    
    return ret;
}




int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier)
{
    if(barrier==NULL)
    {
        return EINVAL;
    }

    pthread_mutex_destroy(&(barrier->barrierp->mutex));
    pthread_cond_destroy(&(barrier->barrierp->cv));

    free(barrier->barrierp);
    
    return 0;
}
 

int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *   barrier)
{
    int    cycle;
    int    rc;
    
    /* check whether the barrier is initialized */
    if (tMPI_Atomic_get( &(barrier->initialized)  ) == 0)
    {
        tMPI_Thread_barrier_init_once(barrier);
    }


    rc = pthread_mutex_lock(&barrier->barrierp->mutex);

    
    if(rc != 0)
        return EBUSY;

    cycle = barrier->cycle;
    
    /* Decrement the count atomically and check if it is zero.
        * This will only be true for the last thread calling us.
        */
    if( --barrier->count <= 0 )
    { 
        barrier->cycle = !barrier->cycle;
        barrier->count = barrier->threshold;
        rc = pthread_cond_broadcast(&barrier->barrierp->cv);
        
        if(rc == 0)
            rc = -1;
    }
    else
    {
        while(cycle == barrier->cycle)
        {
            rc = pthread_cond_wait(&barrier->barrierp->cv,
                                   &barrier->barrierp->mutex);
            if(rc != 0) break;
        }
    }
    
    pthread_mutex_unlock(&barrier->barrierp->mutex);
    return rc;
}

#else

/* just to have some symbols */
int tMPI_Thread_pthreads=0;

#endif /* THREAD_PTHREADS */

