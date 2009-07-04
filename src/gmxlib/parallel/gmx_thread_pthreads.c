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

#ifdef GMX_THREAD_PTHREADS 

/* pthread.h must be the first header, apart from the defines in config.h */
#include <pthread.h>


#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

/*  IMPORTANT:
 *  The Gromacs thread implementation is used to guarantee threadsafe 
 *  operation for the gmx_message.h functions, which in turn is used by
 *  the gmx_memory.h allocation stuff.
 *  This means we cannot use gmx_message() or gmx_new() memory calls 
 *  in the implementation.
 */


#include "gmx_thread.h"
#include "gmx_fatal.h"

/*! \brief System mutex for all one-time initialization 
 *
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 */
static pthread_mutex_t 
gmx_thread_pthreads_system_mtx = PTHREAD_MUTEX_INITIALIZER;


/*! \brief System condition variable for one-time initialization
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 *
 * However, remember that the mutex/condition variables are static, and thus common
 * for all one-time initialization calls. This means that e.g. a thread
 * waiting for the condition variable to signal should check whether the
 * condition signaled came from the one-time initialization we were waiting for, 
 * and if not continue to wait.
 */
static pthread_cond_t 
gmx_thread_pthreads_system_cond = PTHREAD_COND_INITIALIZER;





/*! \brief Pthread implementation of the abstract gmx_thread type
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
struct gmx_thread
{
    pthread_t        pthread; /*!< Pointer to POSIX thread datatype */
};




/*! \brief Pthread implementation of the abstract gmx_thread_key type 
*
*  The contents of this structure depends on the actual threads 
*  implementation used.
*/
struct gmx_thread_key
{
    pthread_key_t    pthread_key; /*!< Pointer to POSIX thread key datatype */
};


/*! \brief Pthread implementation of barrier type. 
 *
 *  The contents of this structure depends on the actual threads 
 *  implementation used.
 */
typedef struct gmx_thread_pthread_barrier
{
    pthread_mutex_t   mutex;     /*!< Lock for the barrier contents          */
    pthread_cond_t    cv;        /*!< Condition to signal barrier completion */
    int               threshold; /*!< Total number of members in barrier     */
    int               count;     /*!< Remaining count before completion      */
    int               cycle;     /*!< Alternating 0/1 to indicate round      */
}
gmx_thread_pthread_barrier_t;



enum gmx_thread_support
gmx_thread_support(void)
{
    return GMX_THREAD_SUPPORT_YES;
}



int
gmx_thread_create   (gmx_thread_t *    thread,
                     void *            (*start_routine)(void *),
                     void *            arg)
{
    int ret;

    if(thread==NULL)
    {
        gmx_fatal(FARGS,"Invalid thread pointer.");
        return EINVAL;
    }
        
    /* We cannot use gromacs memory operations since they
     * are dependent on messages, which in turn require thread support.
     */
    *thread = (gmx_thread_t)malloc(sizeof(struct gmx_thread));

    if(*thread==NULL)
    {
        gmx_fatal(FARGS,"Failed to allocate thread memory.");
        return ENOMEM;
    }
    
	ret=pthread_create(&((*thread)->pthread),NULL,start_routine,arg);
    
    if(ret!=0)
    {
        /* Cannot use gmx_error() since messages use threads for locking */
        gmx_fatal(FARGS,"Failed to create POSIX thread, rc=%d",ret);
        /* Use system memory allocation routines */
        free(*thread);
        return -1;
    }

	return 0;
}



int
gmx_thread_join     (gmx_thread_t     thread,
                     void **          value_ptr)
{
    int ret;

    
    ret = pthread_join( thread->pthread , value_ptr );

    if(ret == 0 )
    {
        /* Free (with system implementation) memory resources used by thread structure */
        free(thread);
    }
    else
    {
        gmx_fatal(FARGS,"Failed to join POSIX thread. rc=%d",ret);
    }
    
    /* Gromacs error numbers are compatible with UNIX, so 
     * we can just pass the return value through.
     */
    
    return ret;
}




int
gmx_thread_mutex_init(gmx_thread_mutex_t *mtx) 
{
    int ret;
    
    if(mtx==NULL)
    {
        return EINVAL;
    }
    
    /*
     * Allocate memory for the pthread mutex. We must use the system malloc
     * here since the gromacs memory allocation depends on gmx_message.h and gmx_thread.h.
     */
    mtx->actual_mutex = malloc(sizeof(pthread_mutex_t));
    
    if(mtx->actual_mutex==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        gmx_fatal(FARGS,"Failed to allocate mutex memory.");
        return ENOMEM;
    }
    
    ret = pthread_mutex_init((pthread_mutex_t *)(mtx->actual_mutex),NULL);
    
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error initializing POSIX mutex. rc=%d");
        /* Use system memory allocation routines */
        free(mtx->actual_mutex);
        return ret;
    }

    mtx->status = GMX_THREAD_ONCE_STATUS_READY;
    
    return 0;
}


int
gmx_thread_mutex_destroy(gmx_thread_mutex_t *mtx) 
{
    int ret;

    if(mtx == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_mutex_destroy( (pthread_mutex_t *)(mtx->actual_mutex) );
    
    /* Must use system free() since Gromacs memory allocation depends on
     * messages and threads working.
     */
    free(mtx->actual_mutex);
    
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error destroying POSIX mutex. rc=%d",ret);
        /* Use system memory allocation routines */
        free(mtx->actual_mutex);
    }
    return ret;
}




static int
gmx_thread_mutex_init_once(gmx_thread_mutex_t *mtx)
{
    int ret;
    
    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the mutex init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 

    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
    * In that case, the mutex will also be unlocked.
    */
    while (mtx->status == GMX_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&gmx_thread_pthreads_system_cond,&gmx_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        mtx->status = GMX_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
        * Only one thread can do it anyway.
        */
        pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
        ret=gmx_thread_mutex_init(mtx);
        pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
        
        /* Status will be marked as ready by gmx_thread_mutex_init(). */ 
        pthread_cond_broadcast (&gmx_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);

    return ret;
}



int
gmx_thread_mutex_lock(gmx_thread_mutex_t *mtx)
{
    int ret;
    
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_mutex_init_once(mtx);
    }
    
    /* The mutex is now guaranteed to be valid. */
    ret=pthread_mutex_lock((pthread_mutex_t *)(mtx->actual_mutex));

    return ret;
}

 


int
gmx_thread_mutex_trylock(gmx_thread_mutex_t *mtx)
{
    int ret;
    
    /* Ccheck whether this mutex is initialized */
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_mutex_init_once(mtx);
    }
    
    /* The mutex is now guaranteed to be valid. */
    ret=pthread_mutex_trylock((pthread_mutex_t *)(mtx->actual_mutex));
    
    return ret;
}



int
gmx_thread_mutex_unlock(gmx_thread_mutex_t *mtx)
{
    int ret;
    
    ret = pthread_mutex_unlock((pthread_mutex_t *)(mtx->actual_mutex));
    
    return ret;
}



int     
gmx_thread_key_create(gmx_thread_key_t *       key,
                      void                   (*destructor)(void *))
{
    int ret;

    if(key==NULL)
    {
        gmx_fatal(FARGS,"Invalid key pointer.");
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
        gmx_fatal(FARGS,"Failed to allocate thread key memory.");
        return ENOMEM;
    }
    
	ret = pthread_key_create(&((*key)->pthread_key),destructor);
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Failed to create thread key, rc=%d.",ret);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(*key);
        return -1;
    }

	return 0;
}


int
gmx_thread_key_delete(gmx_thread_key_t       key)
{
    int ret;

	ret=pthread_key_delete(key->pthread_key);

    if(ret!=0)
    {
        gmx_fatal(FARGS,"Failed to delete thread key, rc=%d.",ret);
        fflush(stderr);
    }
    
    return ret;
}



void *
gmx_thread_getspecific(gmx_thread_key_t   key)
{
    void *p = NULL;

	p=pthread_getspecific(key->pthread_key);

	return p;
}


int
gmx_thread_setspecific(gmx_thread_key_t    key, 
                       void *              value)
{
    int ret;
    
    ret = pthread_setspecific(key->pthread_key,value);
    
    return ret;
}



int
gmx_thread_once(gmx_thread_once_t *     once_control,
                void                    (*init_routine)(void))
{
    /* Do a preliminary check without locking the mutex so we can return 
     * immediately if it is already completed. 
     */
    if(once_control->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        /* Lock the common one-time init mutex so we can check carefully */
        pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
        
        /* If somebody is already working, wait until he is finished.
         * In that case, the mutex will also be unlocked.
        */
        while (once_control->status == GMX_THREAD_ONCE_STATUS_PROGRESS)
            pthread_cond_wait (&gmx_thread_pthreads_system_cond,&gmx_thread_pthreads_system_mtx);
        
        /* Do the actual (locked) check - mutex is locked if we get here */
        if (once_control->status != GMX_THREAD_ONCE_STATUS_READY)
        {
            once_control->status = GMX_THREAD_ONCE_STATUS_PROGRESS;

            /* No need to keep the lock during execution -
             * Only one thread can do it anyway.
             */
            pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
            (*init_routine)();
            pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);

            once_control->status = GMX_THREAD_ONCE_STATUS_READY;
            pthread_cond_broadcast (&gmx_thread_pthreads_system_cond);
        }
        
        pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
    }
    
    return 0;
}
    




int
gmx_thread_cond_init(gmx_thread_cond_t *cond) 
{
    int ret;
    
    if(cond==NULL)
    {
        return EINVAL;
    }
    
    /*
     * Allocate memory for the pthread condition variable. We must use the system malloc
     * here since the gromacs memory allocation depends on gmx_message.h and gmx_thread.h.
     */
    cond->actual_cond = malloc(sizeof(pthread_cond_t));
    
    if(cond->actual_cond==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        gmx_fatal(FARGS,"Failed to allocate condition variable memory.");
        fflush(stderr);
        return ENOMEM;
    }
    
    ret = pthread_cond_init((pthread_cond_t *)(cond->actual_cond),NULL);
    
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error initializing POSIX condition variable. rc=%d",ret);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(cond->actual_cond);
    }
 
    cond->status = GMX_THREAD_ONCE_STATUS_READY;

    return ret;
}


int
gmx_thread_cond_destroy(gmx_thread_cond_t *cond) 
{
    int ret;
    
    if(cond == NULL)
    {
        return EINVAL;
    }
    
    ret = pthread_cond_destroy((pthread_cond_t *)(cond->actual_cond));
    
    /* Must use system free() since Gromacs memory allocation depends on
        * messages and threads working.
        */
    free(cond->actual_cond);
    
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error destroying POSIX condition variable. rc=%d",ret);
        fflush(stderr);
        /* Use system memory allocation routines */
        free(cond->actual_cond);
    }
    return ret;
}




/*! \brief Static init routine for pthread barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for gmx_thread.h
 * 
 * \param cond  Condition variable, must be statically initialized
 *  
 * \return status - 0 on success, or a standard error code.
 */
static int
gmx_thread_cond_init_once(gmx_thread_cond_t *cond)
{
    int ret;
    
    /* This is essentially a copy of the code from the one-time
    * initialization, but with a call to the cond init routine instead.
    * It might seem like overkill, but it will only be executed the first
    * time you call a static condition variable, and it is important to get 
    * the memory barriers right. Trust me, you don't want a deadlock here...
    */ 
    
    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
     * In that case, the mutex will also be unlocked.
     */
    while (cond->status == GMX_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&gmx_thread_pthreads_system_cond,&gmx_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (cond->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        cond->status = GMX_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
         * Only one thread can reach this code!
         */
        pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
        ret=gmx_thread_cond_init(cond);
        pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
        
        /* Status will be marked as GMX_THREAD_ONCE_STATUS_READY by gmx_thread_mutex_init(). */ 
        pthread_cond_broadcast (&gmx_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
    
    return ret;
}



int
gmx_thread_cond_wait(gmx_thread_cond_t *cond, gmx_thread_mutex_t *mtx)
{
    int ret;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_cond_init_once(cond);
    }
    if(mtx->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_mutex_init_once(mtx);
    }
    
    /* The condition variable and mutex are now guaranteed to be valid. */
    ret = pthread_cond_wait( (pthread_cond_t *) (cond->actual_cond) , 
                             (pthread_mutex_t *) (mtx->actual_mutex) );
    
    return ret;
}




int
gmx_thread_cond_signal(gmx_thread_cond_t *cond)
{
    int ret;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_cond_init_once(cond);
    }
    
    /* The condition variable is now guaranteed to be valid. */
    ret = pthread_cond_signal( (pthread_cond_t *) (cond->actual_cond) );
    
    return ret;
}



int
gmx_thread_cond_broadcast(gmx_thread_cond_t *cond)
{
    int ret;
    
    /* Ccheck whether this condition variable is initialized */
    if(cond->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_cond_init_once(cond);
    }
    
    /* The condition variable is now guaranteed to be valid. */
    ret = pthread_cond_broadcast( (pthread_cond_t *) (cond->actual_cond) );
    
    return ret;
}




void
gmx_thread_exit(void *      value_ptr)
{
    pthread_exit(value_ptr);
}




int
gmx_thread_cancel(gmx_thread_t     thread)
{
    return pthread_cancel(thread->pthread);
}


int
gmx_thread_barrier_init(gmx_thread_barrier_t *    barrier,
                        int                       n)
{
    int ret;
    gmx_thread_pthread_barrier_t *p;
    
    if(barrier==NULL)
    {
        return EINVAL;
    }
    
    /*
     * Allocate memory for the barrier variable. We must use the system malloc
     * here since the gromacs memory allocation depends on gmx_message.h and gmx_thread.h.
     */
    barrier->actual_barrier = malloc(sizeof(gmx_thread_pthread_barrier_t));
    
    if(barrier->actual_barrier==NULL)
    {
        /* Write to stderr since we cannot use gromacs messages. */
        gmx_fatal(FARGS,"Failed to allocate barrier memory.");
        return ENOMEM;
    }
    
    p = (gmx_thread_pthread_barrier_t*)barrier->actual_barrier;
    
    ret = pthread_mutex_init(&p->mutex,NULL);
        
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error initializing POSIX mutex. rc=%d",ret);
        /* Use system memory allocation routines */
        free(barrier->actual_barrier);
        return ret;
    }
    
    ret = pthread_cond_init(&p->cv,NULL);
    
    if(ret!=0)
    {
        gmx_fatal(FARGS,"Error initializing POSIX condition variable. rc=%d",ret);
        /* Use system memory allocation routines */
        free(barrier->actual_barrier);
        return ret;
    }
        
    p->threshold = n;
    p->count     = n;
    p->cycle     = 0;

    barrier->status = GMX_THREAD_ONCE_STATUS_READY;

    return 0;
}



int
gmx_thread_barrier_destroy(gmx_thread_barrier_t *barrier)
{
    gmx_thread_pthread_barrier_t *p;
    
    if(barrier==NULL)
    {
        return EINVAL;
    }

    p = (gmx_thread_pthread_barrier_t*)barrier->actual_barrier;
    
    if(barrier->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_fatal(FARGS,"Cannot destroy uninitialized barrier.");
        return EINVAL;
    }
    
    pthread_mutex_destroy(&p->mutex);
    pthread_cond_destroy(&p->cv);
    
    return 0;
}
 


/*! \brief Static init routine for pthread barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for gmx_thread.h
 *
 * \param barrier Statically initialized barrier type
 * \param n       Number of members in barrier
 * 
 * \return status - 0 on success, or a standard error code.
 */
static int
gmx_thread_barrier_init_once(gmx_thread_barrier_t *    barrier,
                             int                       n)
{
    int ret;
    
    /* This is essentially a copy of the code from the general one-time
     * initialization, but with a call to the barrier init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 
    
    /* Lock the common one-time init mutex so we can check carefully */
    pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
    
    /* If somebody is already initializing, wait until he is finished.
     * In that case, the mutex will also be unlocked.
     */
    while (barrier->status == GMX_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&gmx_thread_pthreads_system_cond,&gmx_thread_pthreads_system_mtx);
    
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (barrier->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        barrier->status = GMX_THREAD_ONCE_STATUS_PROGRESS;
        
        /* No need to keep the lock during execution -
        * Only one thread can do it anyway.
        */
        pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
        ret=gmx_thread_barrier_init(barrier,n);
        pthread_mutex_lock (&gmx_thread_pthreads_system_mtx);
        
        /* Status will be marked as ready by gmx_thread_barrier_init(). */ 
        pthread_cond_broadcast (&gmx_thread_pthreads_system_cond);
    }
    else
    {
        ret = 0;
    }
    
    pthread_mutex_unlock (&gmx_thread_pthreads_system_mtx);
    
    return ret;
}



int
gmx_thread_barrier_wait(gmx_thread_barrier_t *   barrier)
{
    int    cycle;
    int    status;
    int    i;
    int    rc;
    gmx_thread_pthread_barrier_t *p;

    if(barrier->status != GMX_THREAD_ONCE_STATUS_READY)
    {
        gmx_thread_barrier_init_once(barrier,barrier->init_threshold);        
    }

    p = (gmx_thread_pthread_barrier_t*)barrier->actual_barrier;
    
    rc = pthread_mutex_lock(&p->mutex);

    
    if(rc != 0)
        return EBUSY;

    cycle = p->cycle;
    
    /* Decrement the count atomically and check if it is zero.
        * This will only be true for the last thread calling us.
        */
    if( --p->count == 0 )
    { 
        p->cycle = !p->cycle;
        p->count = p->threshold;
        rc = pthread_cond_broadcast(&p->cv);
        
        if(rc == 0)
            rc = -1;
    }
    else
    {
        while(cycle == p->cycle)
        {
            rc = pthread_cond_wait(&p->cv,&p->mutex);
            if(rc != 0) break;
        }
    }
    
    pthread_mutex_unlock(&p->mutex);
    return rc;
}



void
gmx_lockfile(FILE *stream)
{
    flockfile(stream);
}


void
gmx_unlockfile(FILE *stream)
{
    funlockfile(stream);
}

#endif /* GMX_THREAD_PTHREADS */
