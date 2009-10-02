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

#ifdef THREAD_WINDOWS

/* the win32 header */
#include <windows.h>


#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"

/*! \brief System mutex for all one-time initialization 
 *
 *  This static variable is necessary in order to make the header file 
 *  independent of the thread library implementation. Anyway, it
 *  will only be locked a handful of times at the start of program execution.
 */
/*
   enum tMPI_Thread_once_status tMPI_Thread_system_lock_state=
   TMPI_THREAD_ONCE_STATUS_NOTCALLED;
   static CRITICAL_SECTION tMPI_Thread_system_lock;
   */
tMPI_Spinlock_t tMPI_Thread_system_lock=TMPI_SPINLOCK_INITIALIZER;


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

struct tMPI_Thread_starter_param
{
    void *(*start_routine)(void*); /* the function */
    void *param; /* its parameter */
};

static DWORD WINAPI tMPI_Win32_thread_starter( LPVOID lpParam ) 
{
    struct tMPI_Thread_starter_param *prm=
              (struct tMPI_Thread_starter_param*)lpParam;

    (prm->start_routine)(prm->param);
    return 0;
}


int tMPI_Thread_create(tMPI_Thread_t *thread,
                       void *(*start_routine)(void *), void *arg)
{
    DWORD thread_id;
    struct tMPI_Thread_starter_param *prm;

    /* a small memory leak to be sure that it doesn't get deallocated 
       once this function ends */
    prm=(struct tMPI_Thread_starter_param*)
              malloc(sizeof(struct tMPI_Thread_starter_param));
    prm->start_routine= start_routine;
    prm->param=arg;

    if(thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid thread pointer.");
        return EINVAL;
    }

    *thread = CreateThread(NULL, 0, tMPI_Win32_thread_starter, prm, 0, 
                           &thread_id);

    if(*thread==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to create thread, error code=%d",
                         GetLastError());
        return -1;
    }

    return 0;
}



int tMPI_Thread_join(tMPI_Thread_t thread, void **value_ptr)
{
    DWORD ret,retval;

    ret = WaitForSingleObject(thread, INFINITE);

    if (ret != 0)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to join thread. error code=%d",
                         GetLastError());
        return -1;
    }

    if (value_ptr)
    {
        if (!GetExitCodeThread(thread, &retval))
        {
            /* TODO: somehow assign value_ptr */
            tMPI_Fatal_error(TMPI_FARGS,
                             "Failed to get thread exit code: error=%d",
                             GetLastError());
            return -1;
        }
    }
    CloseHandle(thread);

    return 0;
}


void tMPI_Thread_exit(void *value_ptr)
{
    /* TODO: fix exit code */
    /* TODO: call destructors for thread-local storage */
    ExitThread( 0 );
}




int tMPI_Thread_cancel(tMPI_Thread_t thread)
{
    if (!TerminateThread( thread, -1) )
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed thread_cancel, error code=%d",
                         GetLastError());
        return -1;
    }
    return 0;
}




int tMPI_Thread_mutex_init(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx==NULL)
    {
        return EINVAL;
    }

    InitializeCriticalSection(&(mtx->cs));
    mtx->init_state = TMPI_THREAD_ONCE_STATUS_READY;

    return 0;
}


int tMPI_Thread_mutex_destroy(tMPI_Thread_mutex_t *mtx) 
{
    if(mtx == NULL)
    {
        return EINVAL;
    }

    DeleteCriticalSection(&(mtx->cs));

    return 0;
}




static int tMPI_Thread_mutex_init_once(tMPI_Thread_mutex_t *mtx)
{
    int ret;

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the mutex init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static mutex, and it is important to get all the 
     * memory barriers right. Trust me, you don't want a deadlock here...
     */ 

    /* Lock the common one-time init mutex so we can check carefully */
    /*EnterCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_lock( &tMPI_Thread_system_lock );


#if 0
    /* If somebody is already initializing, wait until he is finished.
     * In that case, the mutex will also be unlocked.
     */
    while (mtx->status == TMPI_THREAD_ONCE_STATUS_PROGRESS)
        pthread_cond_wait (&tMPI_Thread_pthreads_system_cond,
                           &tMPI_Thread_pthreads_system_mtx);
#endif

    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (mtx->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        /*mtx->status = TMPI_THREAD_ONCE_STATUS_PROGRESS;*/

        /* No need to keep the lock during execution -
         * Only one thread can do it anyway.
         */
        /*pthread_mutex_unlock (&tMPI_Thread_pthreads_system_mtx);*/
        ret=tMPI_Thread_mutex_init(mtx);
        /*pthread_mutex_lock (&tMPI_Thread_pthreads_system_mtx);*/

        /* Status will be marked as ready by tMPI_Thread_mutex_init(). */ 
        /*pthread_cond_broadcast (&tMPI_Thread_pthreads_system_cond);*/
    }
    else
    {
        ret = 0;
    }

    /*LeaveCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_unlock( &tMPI_Thread_system_lock );

    return ret;
}



int tMPI_Thread_mutex_lock(tMPI_Thread_mutex_t *mtx)
{
    /* Ccheck whether this mutex is initialized */
    if(mtx->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }

    /* The mutex is now guaranteed to be valid. */
    EnterCriticalSection( &(mtx->cs) );

    return 0;
}




int tMPI_Thread_mutex_trylock(tMPI_Thread_mutex_t *mtx)
{
    BOOL ret;

    /* Ccheck whether this mutex is initialized */
    if(mtx->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }

    /* The mutex is now guaranteed to be valid. */
    ret=TryEnterCriticalSection( &(mtx->cs) );

    return (ret != 0);
}



int tMPI_Thread_mutex_unlock(tMPI_Thread_mutex_t *mtx)
{
    LeaveCriticalSection( &(mtx->cs) );

    return 0;
}



int tMPI_Thread_key_create(tMPI_Thread_key_t *key, void (*destructor)(void *))
{
    if(key==NULL)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Invalid key pointer.");
        return EINVAL;
    }


    /* TODO: make list of destructors for thread-local storage */
    *key=TlsAlloc();

    if ( *key == TLS_OUT_OF_INDEXES ) 
    {
        tMPI_Fatal_error(TMPI_FARGS,
                         "Failed to create thread key, error code=%d.",
                         GetLastError());
        return -1;
    }

    return 0;
}


int tMPI_Thread_key_delete(tMPI_Thread_key_t key)
{
    TlsFree(key);

    return 0;
}



void * tMPI_Thread_getspecific(tMPI_Thread_key_t key)
{
    void *p = NULL;

    p=TlsGetValue(key);

    return p;
}


int tMPI_Thread_setspecific(tMPI_Thread_key_t key, void *value)
{
    BOOL ret;

    ret = TlsSetValue(key, value);

    return ret==0;
}


static BOOL CALLBACK InitHandleWrapperFunction(PINIT_ONCE InitOnce,
                                               PVOID Parameter,
                                               PVOID *lpContext)
{
    void (*fn)(void)=(void (*)(void))Parameter;

    fn();

    return TRUE;
}

CRITICAL_SECTION tMPI_Once_cs;
tMPI_Spinlock_t tMPI_Once_cs_lock=TMPI_SPINLOCK_INITIALIZER;
volatile int tMPI_Once_init=0;


int tMPI_Thread_once(tMPI_Thread_once_t *once_control, 
                     void (*init_routine)(void))
{
#if 0
    BOOL bStatus;
    bStatus = InitOnceExecuteOnce(once_control, InitHandleWrapperFunction, 
                                  init_routine, NULL);

    if (!bStatus)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed to run thread_once routine");
        return -1;
    }
#else
   /* ugly hack to initialize the critical section */
    if (!tMPI_Once_init)
    {
        tMPI_Spinlock_lock( &tMPI_Once_cs_lock);
        InitializeCriticalSection(&(tMPI_Once_cs));
        tMPI_Once_init=1;
        tMPI_Spinlock_unlock( &tMPI_Once_cs_lock);
    }
    EnterCriticalSection(&tMPI_Once_cs);
    if (*once_control == 0)
    {
        (*init_routine)(); /* call the function */
        *once_control=1; /* flag that we're done */
    }
    LeaveCriticalSection(&tMPI_Once_cs);
#endif
    return 0;
}





int tMPI_Thread_cond_init(tMPI_Thread_cond_t *cond) 
{
    if(cond==NULL)
    {
        return EINVAL;
    }

#if 0
    /* use this code once Vista is the minimum version required */
    InitializeConditionVariable( &(cond->cv) );
#else
    cond->Nwaiters=0;
    InitializeCriticalSection(&(cond->wtr_lock));
    cond->Nrelease=0;
    cond->cycle=0;
    /* a manual reset, unsignalled event */
    cond->ev = CreateEvent(NULL, TRUE, FALSE, NULL); 
#endif

    cond->init_state=TMPI_THREAD_ONCE_STATUS_READY;
    return 0;
}


int tMPI_Thread_cond_destroy(tMPI_Thread_cond_t *cond) 
{
#if 0
    /* use this code once Vista is the minimum version required */
    /* windows doesnt have this function */
#else
    DeleteCriticalSection(&(cond->wtr_lock));
#endif
    return 0;
}



/*! \brief Static init routine for pthread barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for tMPI_Thread.h
 * 
 * \param cond  Condition variable, must be statically initialized
 *  
 * \return status - 0 on success, or a standard error code.
 */
static int tMPI_Thread_cond_init_once(tMPI_Thread_cond_t *cond)
{
    int ret;

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the cond init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static condition variable, and it is important to get 
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */ 
    /* Lock the common one-time init mutex so we can check carefully */
    /*EnterCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_lock( &tMPI_Thread_system_lock );
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (cond->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        ret=tMPI_Thread_cond_init(cond);
    }
    else
    {
        ret = 0;
    }
    /*LeaveCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_unlock( &tMPI_Thread_system_lock );

    return ret;
}



int tMPI_Thread_cond_wait(tMPI_Thread_cond_t *cond, tMPI_Thread_mutex_t *mtx)
{
    BOOL wait_done=FALSE;
    BOOL last_waiter=FALSE;
    int my_cycle;

    /* Ccheck whether this condition variable is initialized */
    if(cond->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    if(mtx->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_mutex_init_once(mtx);
    }
#if 0
    /* use this code once Vista is the minimum version required */
    ret=SleepConditionVariableCS (&(cond->cv), &(mtx->cs), INFINITE);

    if (!ret)
    {
        tMPI_Fatal_error(TMPI_FARGS,"Failed wait for condition, error code=%d",
                         GetLastError());
        return -1;
    }
#else
    /* serially increase waiter count */
    EnterCriticalSection(&(cond->wtr_lock));
    cond->Nwaiters++;
    my_cycle = cond->cycle;
    LeaveCriticalSection(&(cond->wtr_lock));

    /* now it's safe to release the mutex from the fn call */
    LeaveCriticalSection(&(mtx->cs));

    /* Loop a wait until we found out we've waited for the right event.
       Note that this loop is potentially a busy-wait loop in bad
       circumstances (higher priority threads, for example). */
    do
    {
        /* do the actual waiting */
        if (WaitForSingleObject( cond->ev, INFINITE )== WAIT_FAILED)
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed event reset, error code=%d",
                             GetLastError());
            return -1;
        }

        /* serially check whether we got the right event.  */
        EnterCriticalSection(&(cond->wtr_lock));
        wait_done = (cond->Nrelease > 0) && (cond->cycle!=my_cycle);
        LeaveCriticalSection(&(cond->wtr_lock));
    }
    while(!wait_done);

    /* We obtain the mutex from the function call */
    EnterCriticalSection(&(mtx->cs));

    /* we serially decrease the waiter count and release count */
    EnterCriticalSection(&(cond->wtr_lock));
    cond->Nwaiters--;
    cond->Nrelease--;
    last_waiter=(cond->Nrelease==0);
    LeaveCriticalSection(&(cond->wtr_lock));

    /* manually release the event if everybody's done with it */
    if (last_waiter)
    {
        if (!ResetEvent( cond->ev ))
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed event reset, error code=%d",
                             GetLastError());
            return -1;
        }
    }
#endif

    return 0;
}




int tMPI_Thread_cond_signal(tMPI_Thread_cond_t *cond)
{
    /* Ccheck whether this condition variable is initialized */
    if(cond->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    /* The condition variable is now guaranteed to be valid. */
#if 0
    /* use this code once Vista is the minimum version required */
    WakeConditionVariable( &(cond->cv) );
#else
    EnterCriticalSection(&(cond->wtr_lock));
    /* check if we're not still busy with a release. If we are, do nothing. */
    if (cond->Nwaiters > cond->Nrelease)
    {
        cond->Nrelease++;
        cond->cycle++;
        if (!SetEvent(cond->ev)) /* actually release the waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed SetEvent, error code=%d",
                             GetLastError());
            return -1;
        }
    }
    LeaveCriticalSection(&(cond->wtr_lock));
#endif

    return 0;
}



int tMPI_Thread_cond_broadcast(tMPI_Thread_cond_t *cond)
{
    /* Ccheck whether this condition variable is initialized */
    if(cond->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_cond_init_once(cond);
    }
    /* The condition variable is now guaranteed to be valid. */
#if 0
    /* use this code once Vista is the minimum version required */
    WakeAllConditionVariable( &(cond->cv) );
#else
    EnterCriticalSection(&(cond->wtr_lock));
    /* check whether there are any waiters */
    if (cond->Nwaiters > 0)
    {
        cond->Nrelease=cond->Nwaiters;
        cond->cycle++;
        if (!SetEvent(cond->ev)) /* actually release the waiting threads */
        {
            tMPI_Fatal_error(TMPI_FARGS,"Failed SetEvent, error code=%d",
                             GetLastError());
            return -1;
        }
    }
    LeaveCriticalSection(&(cond->wtr_lock));
#endif
    return 0;
}

#ifdef TMPI_RWLOCK
int tMPI_Thread_rwlock_init(tMPI_Thread_rwlock_t *rwlock)
{
    InitializeSRWLock(rwlock);
    return 0; /* no error value returned by the above function in Windows */
}

int tMPI_Thread_rwlock_destroy(tMPI_Thread_rwlock_t *rwlock)
{
    /* Windows doesn't define a function for this. Presumably this means
       that the lock's size is too small to bother with */
    return 0;
}   

int tMPI_Thread_rwlock_rdlock(tMPI_Thread_rwlock_t *rwlock)
{
    AcquireSRWLockShared(rwlock);    
    return 0; /* no error value returned by the above function in Windows */
}   

int tMPI_Thread_rwlock_tryrdlock(tMPI_Thread_rwlock_t *rwlock)
{
    /*TryAcquireSRWLockShared(rwlock);   */ 
    return 0; /* no error value returned by the above function in Windows */
}

int tMPI_Thread_rwlock_wrlock(tMPI_Thread_rwlock_t *rwlock)
{
    AcquireSRWLockExclusive(rwlock);    
    return 0; /* no error value returned by the above function in Windows */
}

int tMPI_Thread_rwlock_trywrlock(tMPI_Thread_rwlock_t *rwlock)
{
    /*TryAcquireSRWLockExclusive(rwlock); */   
    return 0; /* no error value returned by the above function in Windows */
}

int tMPI_Thread_rwlock_rdunlock(tMPI_Thread_rwlock_t *rwlock)
{
    ReleaseSRWLockShared(rwlock);    
    return 0; /* no error value returned by the above function in Windows */
}

int tMPI_Thread_rwlock_wrunlock(tMPI_Thread_rwlock_t *rwlock)
{
    ReleaseSRWLockExclusive(rwlock);    
    return 0; /* no error value returned by the above function in Windows */
}
#endif





int tMPI_Thread_barrier_init(tMPI_Thread_barrier_t *barrier, int n)
{
    if(barrier==NULL)
    {
        return EINVAL;
    }


#if 0
 /* use this once Vista is the oldest supported windows version: */
    InitializeCriticalSection(&(barrier->cs));
    InitializeConditionVariable(&(barrier->cv));
#else
    tMPI_Thread_mutex_init(&(barrier->cs));
    tMPI_Thread_cond_init(&(barrier->cv));
#endif

    barrier->threshold = n;
    barrier->count     = n;
    barrier->cycle     = 0;

    barrier->init_state = TMPI_THREAD_ONCE_STATUS_READY;

    return 0;
}



int tMPI_Thread_barrier_destroy(tMPI_Thread_barrier_t *barrier)
{   
    if(barrier==NULL)
    {
        return EINVAL;
    }

#if 0
    DeleteCriticalSection(&(barrier->cs));
#else
    tMPI_Thread_mutex_destroy(&(barrier->cs));
#endif

    tMPI_Thread_cond_destroy(&(barrier->cv));

    return 0;
}



/*! \brief Static init routine for pthread barrier 
 *
 * \internal
 *
 * This is only used as a wrapper to enable static initialization
 * of posix thread types together with out abstraction layer for tMPI_Thread.h
 *
 * \param barrier Statically initialized barrier type
 * \param n       Number of members in barrier
 * 
 * \return status - 0 on success, or a standard error code.
 */
static int tMPI_Thread_barrier_init_once(tMPI_Thread_barrier_t *barrier, int n)
{
    int ret;

    /* This is essentially a copy of the code from the one-time
     * initialization, but with a call to the cond init routine instead.
     * It might seem like overkill, but it will only be executed the first
     * time you call a static condition variable, and it is important to get 
     * the memory barriers right. Trust me, you don't want a deadlock here...
     */ 
    /* Lock the common one-time init mutex so we can check carefully */
    /*EnterCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_lock( &tMPI_Thread_system_lock );
    /* Do the actual (locked) check - system mutex is locked if we get here */
    if (barrier->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        ret=tMPI_Thread_barrier_init(barrier, n);
    }
    else
    {
        ret = 0;
    }
    /*LeaveCriticalSection( &tMPI_Thread_system_lock );*/
    tMPI_Spinlock_lock( &tMPI_Thread_system_lock );

    return ret;
}



int tMPI_Thread_barrier_wait(tMPI_Thread_barrier_t *barrier)
{
    int    cycle;
    BOOL    rc=FALSE;
    int     ret=0;
    /*tMPI_Thread_pthread_barrier_t *p;*/

    if(barrier->init_state != TMPI_THREAD_ONCE_STATUS_READY)
    {
        tMPI_Thread_barrier_init_once(barrier,barrier->threshold);        
    }

    /*p = (tMPI_Thread_pthread_barrier_t*)barrier->actual_barrier;*/
#if 0
    EnterCriticalSection( &(barrier->cs)  );
#else
    tMPI_Thread_mutex_lock( &(barrier->cs) );
#endif



    cycle = barrier->cycle;

    /* Decrement the count atomically and check if it is zero.
     * This will only be true for the last thread calling us.
     */
    if( --(barrier->count) <= 0 )
    { 
        barrier->cycle = !barrier->cycle;
        barrier->count = barrier->threshold;
#if 0
        WakeAllConditionVariable( &(barrier->cv) );
#else
        tMPI_Thread_cond_broadcast( &(barrier->cv) );
#endif
    }
    else
    {
        while(cycle == barrier->cycle)
        {
#if 0
            rc=SleepConditionVariableCS (&(barrier->cv), &(barrier->cs), 
                                         INFINITE);
            if(!rc) 
            {
                ret=-1;
                break;
            }
#else
            rc = tMPI_Thread_cond_wait(&barrier->cv,&barrier->cs);
            if(rc != 0) break;
#endif
        }
    }
#if 0
    LeaveCriticalSection( &(barrier->cs)  );
#else
    tMPI_Thread_mutex_unlock( &(barrier->cs) );
#endif
    return ret;
}



void tMPI_lockfile(FILE *stream)
{
    /*    flockfile(stream);*/
    /* TODO: implement this */
}


void tMPI_unlockfile(FILE *stream)
{
    /*    funlockfile(stream);*/
    /* TODO: implement this */
}

#endif /* THREAD_WINDOWS  */
