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



/* Include the defines that determine which thread library to use.
* We do not use HAVE_PTHREAD_H directly, since we might want to
* turn off thread support explicity (e.g. for debugging).
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef THREAD_PTHREADS 

/* pthread.h must be the first header, apart from the defines in config.h */
#include <pthread.h>


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "thread_mpi/threads.h"
/*#include "thread_mpi_impl.h"*/




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



int tMPI_Thread_create(tMPI_Thread_t *   thread,
                       void *            (*start_routine)(void *),
                       void *            arg)
{
    int ret;

    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }

    ret=pthread_create(thread,NULL,start_routine,arg);

    if(ret!=0)
    {
        /* Cannot use tMPI_error() since messages use threads for locking */
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create POSIX thread, rc=%d",ret);
        /* Use system memory allocation routines */
        return -1;
    }

    return 0;
}



int tMPI_Thread_join     (tMPI_Thread_t     thread,
                          void **          value_ptr)
{
    int ret;

    
    ret = pthread_join( thread, value_ptr );

    if(ret != 0 )
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to join POSIX thread. rc=%d",ret);
    }
    return ret;
}




int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx) 
{
    int ret;
    
    if(mtx==NULL)
    {
        return EINVAL;
    }
   
    ret = pthread_mutex_init(mtx,NULL);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX mutex. rc=%d");
        /* Use system memory allocation routines */
        return ret;
    }

    return 0;
}


int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx) 
{
    int ret;

    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_mutex_destroy( mtx );
    
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
   
    ret=pthread_mutex_lock(mtx);

    return ret;
}

 


int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    int ret;
    ret=pthread_mutex_trylock(mtx);
    
    return ret;
}



int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{
    int ret;
    
    ret = pthread_mutex_unlock(mtx);
    
    return ret;
}



int tMPI_Thread_key_create(tMPI_Thread_key_t *       key, 
        void                      (*destructor)(void *))
{
    int ret;

    if(key==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid key pointer.");
        return EINVAL;
    }


    ret = pthread_key_create(key, destructor);
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create thread key, rc=%d.",ret);
        fflush(stderr);
        return -1;
    }

    return 0;
}


int tMPI_Thread_key_delete(tMPI_Thread_key_t       key)
{
    int ret;

    ret=pthread_key_delete(key);

    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to delete thread key, rc=%d.",ret);
        fflush(stderr);
    }
    
    return ret;
}



void * tMPI_Thread_getspecific(tMPI_Thread_key_t   key)
{
    void *p = NULL;

    p=pthread_getspecific(key);

    return p;
}


int tMPI_Thread_setspecific(tMPI_Thread_key_t    key, 
                            void *              value)
{
    int ret;
    
    ret=pthread_setspecific(key,value);
    
    return ret;
}



int tMPI_Thread_once(tMPI_Thread_once_t *     once_control,
                     void                    (*init_routine)(void))
{
    int ret;
    ret=pthread_once(once_control, init_routine);
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed run thread_once, rc=%d.",ret);
        fflush(stderr);
        
    }
    return ret;
}
    




int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond) 
{
    int ret;
    
    if(cond==NULL)
    {
        return EINVAL;
    }
   
    ret = pthread_cond_init(cond, NULL);
    
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX condition variable. rc=%d",ret);
        fflush(stderr);
    }
    return ret;
}


int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond) 
{
    int ret;
    
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_cond_destroy(cond);
   
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
   
    ret = pthread_cond_wait( cond, mtx );
    
    return ret;
}




int tMPI_Thread_cond_signal(tMPI_Thread_cond_t *cond)
{
    int ret;
    
    ret = pthread_cond_signal( cond );
    
    return ret;
}



int tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *cond)
{
    int ret;
   
    ret = pthread_cond_broadcast( cond );
    
    return ret;
}




void tMPI_Thread_exit(void *      value_ptr)
{
    pthread_exit(value_ptr);
}




int tMPI_Thread_cancel(tMPI_Thread_t     thread)
{
    return pthread_cancel(thread);
}


int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier,
                             int                   n)
{
    int ret;
    /*tMPI_Thread_pthread_barrier_t *p;*/
    
    if(barrier==NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_mutex_init(&(barrier->mutex),NULL);
        
    if(ret!=0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Error initializing POSIX mutex. rc=%d",
                         ret);
        return ret;
    }
    
    ret = pthread_cond_init(&(barrier->cv),NULL);
    
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

    return 0;
}



int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier)
{
    if(barrier==NULL)
    {
        return EINVAL;
    }

    pthread_mutex_destroy(&(barrier->mutex));
    pthread_cond_destroy(&(barrier->cv));
    
    return 0;
}
 

int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *   barrier)
{
    int    cycle;
    int    rc;
    
    rc = pthread_mutex_lock(&barrier->mutex);

    
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
        rc = pthread_cond_broadcast(&barrier->cv);
        
        if(rc == 0)
            rc = -1;
    }
    else
    {
        while(cycle == barrier->cycle)
        {
            rc = pthread_cond_wait(&barrier->cv,&barrier->mutex);
            if(rc != 0) break;
        }
    }
    
    pthread_mutex_unlock(&barrier->mutex);
    return rc;
}



void tMPI_Lockfile(FILE *stream)
{
    flockfile(stream);
}


void tMPI_Unlockfile(FILE *stream)
{
    funlockfile(stream);
}

#endif /* THREAD_PTHREADS */
