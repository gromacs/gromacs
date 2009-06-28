/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
*
* 
* This file is part of Gromacs        Copyright (c) 1991-2004
* David van der Spoel, Erik Lindahl, University of Groningen.
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* To help us fund GROMACS development, we humbly ask that you cite
* the research papers on the package. Check out http://www.gromacs.org
* 
* And Hey:
* Gnomes, ROck Monsters And Chili Sauce
*/

/* Include the defines that determine which thread library to use.
* We do not use HAVE_PTHREAD_H directly, since we might want to
* turn off thread support explicity (e.g. for debugging).
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if ( !defined(GMX_THREAD_PTHREADS) && !defined(GMX_THREADS_WINDOWS) )
/* Hide this implementation from Doxygen, since it conflicts with pthreads */
#ifndef DOXYGEN




/*  IMPORTANT:
 *  The Gromacs thread implementation is used to guarantee threadsafe 
 *  operation for the gmx_message.h functions, which in turn is used by
 *  the gmx_memory.h allocation stuff.
 *  This means we cannot use gmx_message() or gmx_new() memory calls 
 *  in the implementation.
 */

#include "gmx_thread.h"


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

#include "gmx_fatal.h"

/* Number of thread-specific storage elements.
 *
 * It might seem stupid to make the length constant instead of dynamic,
 * but e.g. the Pthreads standard only guarantees 128 key entries, so for
 * portability we shouldn't use more.
 */ 
#define GMX_THREAD_NOTHREADS_NSPECIFICDATA 128


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
gmx_thread_nothreads_specificdata[GMX_THREAD_NOTHREADS_NSPECIFICDATA];


/*! \brief flag to determine if static key storage has been initialized.
*/
static int            
gmx_thread_nothreads_specific_init_done = 0;



/* Dummy implementation of abstract Gromacs thread datatype */
struct gmx_thread
{
    int      locked; /* Just useful for debugging */
};




/* Dummy implementation of abstract Gromacs thread-specific data key type */
struct gmx_thread_key
{
    int      index;  /*!< Index into our static list of private data */
};




enum gmx_thread_support
gmx_thread_support(void)
{
    return GMX_THREAD_SUPPORT_NO;
}



int
gmx_thread_create   (gmx_thread_t *    thread,
                     void *            (*start_routine)(void *),
                     void *            arg)
{
    gmx_fatal(FARGS,"Cannot start threads without thread support.\n");
    
    return 0;
}



int
gmx_thread_join     (gmx_thread_t     thread,
                     void **          value_ptr)
{
    gmx_fatal(FARGS,"Cannot join threads without thread support.\n");
    
    return 0;
}




int
gmx_thread_mutex_init(gmx_thread_mutex_t *mtx) 
{
    if(mtx==NULL)
    {
        return EINVAL;
    }
    
    /* We use an integer as a mutex for debugging. */
    mtx->actual_mutex = malloc(sizeof(int));
    
    if(mtx->actual_mutex==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Gromacs error [%s, line %d]: Failed to allocate mutex memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }
    
    * ( (int *) (mtx->actual_mutex) ) = 0;  /* Unlocked */
        
    mtx->status = GMX_THREAD_ONCE_STATUS_READY;

    return 0;
}



int
gmx_thread_mutex_destroy(gmx_thread_mutex_t *mtx) 
{
    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    /* Must use system free() since Gromacs memory allocation depends on
     * messages and threads working.
     */
    free(mtx->actual_mutex);
    
    return 0;
}




static int
gmx_thread_mutex_init_once(gmx_thread_mutex_t *mtx)
{
    int rc;
    /* No threads = nothing to worry about */
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
        rc = gmx_thread_mutex_init(mtx);
    else
        rc = 0;

    return rc;
}



int
gmx_thread_mutex_lock(gmx_thread_mutex_t *mtx)
{
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_mutex_init_once(mtx);
    }
    
    /* The mutex is now guaranteed to be valid. */
    * ( (int *) (mtx->actual_mutex) ) = 1;

    return 0;
}

 


int
gmx_thread_mutex_trylock(gmx_thread_mutex_t *mtx)
{
    int ret;
    int *p;

    /* Ccheck whether this mutex is initialized */
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_mutex_init_once(mtx);
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
gmx_thread_mutex_unlock(gmx_thread_mutex_t *mtx)
{    
    * ( (int *) (mtx->actual_mutex) ) = 0;  /* Unlocked */
    
    return 0;
}



int
gmx_thread_key_create(gmx_thread_key_t *       key,
                      void                   (*destructor)(void *))
{
    int i;

    if(key==NULL)
    {
        fprintf(stderr,
                "Gromacs error [%s, line %d]: Invalid key pointer.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return EINVAL;
    }
    
    /*
     * Allocate memory for the pthread key. We must use the system malloc
     * here since the gromacs memory allocation depends on gmx_message.h and gmx_thread.h.
     */
    *key = (gmx_thread_key_t)malloc(sizeof(struct gmx_thread_key));
    
    if(*key==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        fprintf(stderr,
                "Gromacs error [%s, line %d]: Failed to allocate thread key memory.\n",
                __FILE__,__LINE__);
        fflush(stderr);
        return ENOMEM;
    }

    if( gmx_thread_nothreads_specific_init_done == 0)
    {

        for(i=0;i<GMX_THREAD_NOTHREADS_NSPECIFICDATA;i++)
        {
            gmx_thread_nothreads_specificdata[i].used       = 0;
            gmx_thread_nothreads_specificdata[i].value      = NULL;
            gmx_thread_nothreads_specificdata[i].destructor = NULL;
        }
        
        gmx_thread_nothreads_specific_init_done = 1;   
    }
        
    /* Try to find an empty spot */
    for(i=0;i<GMX_THREAD_NOTHREADS_NSPECIFICDATA;i++)
    {
        if(gmx_thread_nothreads_specificdata[i].used == 0)
            break;
    }

    if(i==GMX_THREAD_NOTHREADS_NSPECIFICDATA)
    {
        gmx_fatal(FARGS,"Already used all %d private data keys.\n");
        return -1;
    }
    
    gmx_thread_nothreads_specificdata[i].used = 1;
    (*key)->index = i;    
    
	return 0;
}


int
gmx_thread_key_delete(gmx_thread_key_t        key)
{
    int      i;
    void *   value;             
    void    (*destructor)(void *);    

    if ( gmx_thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return EINVAL;
    }
    
    i = key->index;
    
    if(gmx_thread_nothreads_specificdata[i].value != NULL)
    {
        destructor = gmx_thread_nothreads_specificdata[i].destructor;
        value      = gmx_thread_nothreads_specificdata[i].value;
        (*destructor)(value);
    }
    
    free(key);
    
    return 0;
}



void *
gmx_thread_getspecific(gmx_thread_key_t  key)
{
    if ( gmx_thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return NULL;
    }

    return gmx_thread_nothreads_specificdata[key->index].value;
}


int
gmx_thread_setspecific(gmx_thread_key_t    key, 
                       void *              value)
{
    
    if ( gmx_thread_nothreads_specific_init_done == 0 || key == NULL)
    {
        return EINVAL;
    }
    
    gmx_thread_nothreads_specificdata[key->index].value = value;

    return 0;
}



int
gmx_thread_once(gmx_thread_once_t *     once_control,
                void                    (*init_routine)(void))
{
    if(once_control->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        (*init_routine)();
        once_control->status = GMX_THREAD_ONCE_STATUS_READY;
    }
    
    return 0;
}
    




int
gmx_thread_cond_init(gmx_thread_cond_t *   cond) 
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
    cond->status = GMX_THREAD_ONCE_STATUS_READY;
    
    return 0;
}


int
gmx_thread_cond_destroy(gmx_thread_cond_t *   cond) 
{
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    return 0;
}




static int
gmx_thread_cond_init_once(gmx_thread_cond_t *     cond)
{
    /* No threads = nothing to worry about */
    if(cond->status != GMX_THREAD_ONCE_STATUS_READY)
        return gmx_thread_cond_init(cond);
    return 0;
}



int
gmx_thread_cond_wait(gmx_thread_cond_t *   cond,
                     gmx_thread_mutex_t *  mtx)
{
    gmx_fatal(FARGS,"Called gmx_thread_cond_wait() without thread support. This is a major\n"
              "error, since nobody else can signal the condition variable - exiting.\n");
    
    return 0;
}




int
gmx_thread_cond_broadcast(gmx_thread_cond_t *   cond)
{
    gmx_fatal(FARGS,"Called gmx_thread_broadcast() without thread support.\n");
    
    return 0;
}




void
gmx_thread_exit(void *      value_ptr)
{
    gmx_fatal(FARGS,"Called gmx_thread_exit() without thread support.\n");
}




int
gmx_thread_cancel(gmx_thread_t    thread)
{
    gmx_fatal(FARGS,"Called gmx_thread_cancel() without thread support.\n");
    
    return 0;
}



int
gmx_thread_barrier_init(gmx_thread_barrier_t *    barrier,
                        int                       n)
{
    /* Barriers don't do anything without multiple threads */
    return 0;
}


int
gmx_thread_barrier_destroy(gmx_thread_barrier_t *    barrier)
{
    /* Barriers don't do anything without multiple threads */
    return 0;
}


int
gmx_thread_barrier_wait(gmx_thread_barrier_t *    barrier)
{
    /* Barriers don't do anything without multiple threads.
     * Since we are the only thread we're the master!
     */
    return -1;
}




void
gmx_lockfile(FILE *    stream)
{
    /* Nothing to worry about without threads */
    return;
}


void
gmx_unlockfile(FILE *   stream)
{
    /* Nothing to worry about without threads */
    return;
}



 
#endif /* ifndef DOXYGEN */
 

#endif /* PTHREADS OR WINDOWS */
