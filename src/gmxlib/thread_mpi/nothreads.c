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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Hide this implementation from Doxygen, since it conflicts with pthreads */
#ifndef DOXYGEN

#if ! (defined(THREAD_PTHREADS) || defined(THREAD_WINDOWS))

#include "thread_mpi/thread.h"


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "tMPI_Fatal_error.h"

/* Number of thread-specific storage elements.
 *
 * It might seem stupid to make the length constant instead of dynamic,
 * but e.g. the Pthreads standard only guarantees 128 key entries, so for
 * portability we shouldn't use more.
 */ 
#define TMPI_THREAD_NOTHREADS_NSPECIFICDATA 128


/*
 * Static variables used to keep track of thread-private data
 * when we are compiling without thread support. This is actually quite
 * simple - we just keep a list of structures with the value (pointer)
 * and destructor.
 *
 */
static struct
{
    int            used;                 /*!< 1 if this element is used, 0 if it is free */
    void *         value;                /*!< Pointer to the specific data for key n     */
    void          (*destructor)(void *); /*!< Pointer to destructor function for key n   */
} 
tMPI_Thread_nothreads_specificdata[TMPI_THREAD_NOTHREADS_NSPECIFICDATA];


/*! \brief flag to determine if static key storage has been initialized.
*/
static int            
tMPI_Thread_nothreads_specific_init_done = 0;



/* Dummy implementation of abstract thread datatype */
struct tMPI_Thread
{
    int      locked; /* Just useful for debugging */
};




/* Dummy implementation of abstract thread-specific data key type */
struct tMPI_Thread_key
{
    int      index;  /*!< Index into our static list of private data */
};




enum tMPI_Thread_support
tMPI_Thread_support(void)
{
    return TMPI_THREAD_SUPPORT_NO;
}



int
tMPI_Thread_create   (tMPI_Thread_t *    thread,
                     void *            (*start_routine)(void *),
                     void *            arg)
{
    tMPI_Fatal_error(FARGS,"Cannot start threads without thread support.\n");
    
    return 0;
}



int
tMPI_Thread_join     (tMPI_Thread_t     thread,
                     void **          value_ptr)
{
    tMPI_Fatal_error(FARGS,"Cannot join threads without thread support.\n");
    
    return 0;
}




int
tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx==NULL)
    {
        return EINVAL;
    }
    
    /* We use an integer as a mutex for debugging. */
    mtx->actual_mutex = malloc(sizeof(int));
    
    if(mtx->actual_mutex==NULL)
    {
        fprintf(stderr,
                "error [%s, line %d]: Failed to allocate mutex memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    * ( (int *) (mtx->actual_mutex) ) = 0;  /* Unlocked */
        
    mtx->status = TMPI_THREAD_ONCE_STATUS_READY;

    return 0;
}



int
tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    /* Must use system free() since memory allocation depends on
     * messages and threads working.
     */
    free(mtx->actual_mutex);
    
    return 0;
}




static int
tMPI_Thread_mutex_init_once(tMPI_Thread_mutex_t *mtx)
{
    int rc;
    /* No threads = nothing to worry about */
    if(mtx->status != TMPI_THREAD_ONCE_STATUS_READY)
        rc = tMPI_Thread_mutex_init(mtx);
    else
        rc = 0;

    return rc;
}



int
tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx)
{
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }
    
    /* The mutex is now guaranteed to be valid. */
    * ( (int *) (mtx->actual_mutex) ) = 1;

    return 0;
}

 


int
tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    int ret;
    int *p;

    /* Ccheck whether this mutex is initialized */
    if(mtx->status != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }
    
    p = (int *) (mtx->actual_mutex);

    if( *p == 0 )
    {
        *p = 1;
        ret = 0;
    }
    else
    {
        ret = EBUSY;
    }
    
    return ret;
}



int
tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{    
    * ( (int *) (mtx->actual_mutex) ) = 0;  /* Unlocked */
    
    return 0;
}



int
tMPI_Thread_key_create(tMPI_Thread_key_t *       key,
                      void                   (*destructor)(void *))
{
    int i;

    if(key==NULL)
    {
        fprintf(stderr,
                "error [%s, line %d]: Invalid key pointer.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }
    
    /*
     * Allocate memory for the pthread key. We must use the system malloc
     * here since the other memory allocation depends on tMPI_Thread.h.
     */
    *key = malloc(sizeof(struct tMPI_Thread_key));
    
    if(*key==NULL)
    {
        fprintf(stderr,
                "error [%s, line %d]: Failed to allocate thread key memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }

    if( tMPI_Thread_nothreads_specific_init_done == 0)
    {

        for(i=0;i<TMPI_THREAD_NOTHREADS_NSPECIFICDATA;i++)
        {
            tMPI_Thread_nothreads_specificdata[i].used       = 0;
            tMPI_Thread_nothreads_specificdata[i].value      = NULL;
            tMPI_Thread_nothreads_specificdata[i].destructor = NULL;
        }
        
        tMPI_Thread_nothreads_specific_init_done = 1;   
    }
        
    /* Try to find an empty spot */
    for(i=0;i<TMPI_THREAD_NOTHREADS_NSPECIFICDATA;i++)
    {
        if(tMPI_Thread_nothreads_specificdata[i].used == 0)
            break;
    }

    if(i==TMPI_THREAD_NOTHREADS_NSPECIFICDATA)
    {
        tMPI_Fatal_error(FARGS,"Already used all %d private data keys.\n");
        return -1;
    }
    
    tMPI_Thread_nothreads_specificdata[i].used = 1;
    (*key)->index = i;    
    
	return 0;
}


int
tMPI_Thread_key_delete(tMPI_Thread_key_t        key)
{
    int      i;
    void *   value;             
    void    (*destructor)(void *);    

    if ( tMPI_Thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return EINVAL;
    }
    
    i = key->index;
    
    if(tMPI_Thread_nothreads_specificdata[i].value != NULL)
    {
        destructor = tMPI_Thread_nothreads_specificdata[i].destructor;
        value      = tMPI_Thread_nothreads_specificdata[i].value;
        (*destructor)(value);
    }
    
    free(key);
    
    return 0;
}



void *
tMPI_Thread_getspecific(tMPI_Thread_key_t  key)
{
    if ( tMPI_Thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return NULL;
    }

    return tMPI_Thread_nothreads_specificdata[key->index].value;
}


int
tMPI_Thread_setspecific(tMPI_Thread_key_t    key, 
                       void *              value)
{
    
    if ( tMPI_Thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return EINVAL;
    }
    
    tMPI_Thread_nothreads_specificdata[key->index].value = value;

    return 0;
}



int
tMPI_Thread_once(tMPI_Thread_once_t *     once_control,
                void                    (*init_routine)(void))
{
    if(once_control->status != TMPI_THREAD_ONCE_STATUS_READY)
    {
        (*init_routine)();
        once_control->status = TMPI_THREAD_ONCE_STATUS_READY;
    }
    
    return 0;
}
    




int
tMPI_Thread_cond_init(tMPI_Thread_cond_t *   cond) 
{
    int ret;
    
    /* Condition variables are completely useless without threads, since
     * we cannot wait for some signal to happen if we are the only thread
     * executing.
     *
     * Still, as long as we don't try to call the wait routine it won't hurt,
     * and we want the contents to be clean for debugging.
     */
    cond->actual_cond = NULL;
    cond->status = TMPI_THREAD_ONCE_STATUS_READY;
    
    return 0;
}


int
tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *   cond) 
{
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    return 0;
}




static int
tMPI_Thread_cond_init_once(tMPI_Thread_cond_t *     cond)
{
    /* No threads = nothing to worry about */
    if(cond->status != TMPI_THREAD_ONCE_STATUS_READY)
        return tMPI_Thread_cond_init(cond);
}



int
tMPI_Thread_cond_wait(tMPI_Thread_cond_t *   cond,
                     tMPI_Thread_mutex_t *  mtx)
{
    tMPI_Fatal_error(FARGS,"Called tMPI_Thread_cond_wait() without thread support. This is a major\n"
              "error, since nobody else can signal the condition variable - exiting.\n");
    
    return 0;
}




int
tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *   cond)
{
    tMPI_Fatal_error(FARGS,"Called tMPI_Thread_broadcast() without thread support.\n");
    
    return 0;
}




void
tMPI_Thread_exit(void *      value_ptr)
{
    tMPI_Fatal_error(FARGS,"Called tMPI_Thread_exit() without thread support.\n");
}




int
tMPI_Thread_cancel(tMPI_Thread_t    thread)
{
    tMPI_Fatal_error(FARGS,"Called tMPI_Thread_cancel() without thread support.\n");
    
    return 0;
}



int
tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *    barrier,
                        int                       n)
{
    /* Barriers don't do anything without multiple threads */
    return 0;
}


int
tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *    barrier)
{
    /* Barriers don't do anything without multiple threads */
    return 0;
}


int
tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *    barrier)
{
    /* Barriers don't do anything without multiple threads.
     * Since we are the only thread we're the master!
     */
    return -1;
}




void
tMPI_lockfile(FILE *    stream)
{
    /* Nothing to worry about without threads */
    return;
}


void
tMPI_unlockfile(FILE *   stream)
{
    /* Nothing to worry about without threads */
    return;
}



 
#endif /* ifndef DOXYGEN */
 
#endif /* no threads */
