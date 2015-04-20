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


#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if !(defined( _WIN32 ) || defined( _WIN64 ) )
#include <sys/time.h>

#endif


#include "impl.h"

#ifdef TMPI_TRACE
#include <stdarg.h>
#endif






/* there are a few global variables that maintain information about the
   running threads. Some are defined by the MPI standard: */
/* TMPI_COMM_WORLD is in tmpi_malloc.c due to technical reasons */
tMPI_Group TMPI_GROUP_EMPTY = NULL;


/* the threads themselves (tmpi_comm only contains lists of pointers to this
      structure */
struct tmpi_thread *threads  = NULL;
int                 Nthreads = 0;

/* thread info */
tMPI_Thread_key_t id_key; /* the key to get the thread id */



/* whether MPI has finalized (we need this to distinguish pre-inited from
       post-finalized states */
static tmpi_bool tmpi_finalized = FALSE;

/* misc. global information about MPI */
struct tmpi_global *tmpi_global = NULL;






/* start N threads with argc, argv (used by tMPI_Init)*/
int tMPI_Start_threads(tmpi_bool main_returns, int N,
                       tMPI_Affinity_strategy aff_strategy,
                       int *argc, char ***argv,
                       void (*start_fn)(void*), void *start_arg,
                       int (*start_fn_main)(int, char**));

/* starter function for threads; takes a void pointer to a
      struct tmpi_starter_, which calls main() if tmpi_start_.fn == NULL */
static void* tMPI_Thread_starter(void *arg);

/* allocate and initialize the data associated with a thread structure */
static int tMPI_Thread_init(struct tmpi_thread *th);
/* deallocate the data associated with a thread structure */
static void tMPI_Thread_destroy(struct tmpi_thread *th);




#ifdef TMPI_TRACE
void tMPI_Trace_print(const char *fmt, ...)
{
    va_list                    argp;
    struct tmpi_thread       * th  = NULL;
    static tMPI_Thread_mutex_t mtx = TMPI_THREAD_MUTEX_INITIALIZER;

    /* don't check for errors during trace */
    tMPI_Thread_mutex_lock(&mtx);
    if (threads)
    {
        th = tMPI_Get_current();
        printf("THREAD %02d: ", (int)(th-threads));
    }
    else
    {
        printf("THREAD main: ");
    }
    va_start(argp, fmt);
    vprintf(fmt, argp);
    printf("\n");
    fflush(stdout);
    va_end(argp);
    tMPI_Thread_mutex_unlock(&mtx);
}
#endif


tmpi_bool tMPI_Is_master(void)
{
    /* if there are no other threads, we're the main thread */
    if ( (!TMPI_COMM_WORLD) || TMPI_COMM_WORLD->grp.N == 0)
    {
        return TRUE;
    }

    /* otherwise we know this through thread specific data: */
    /* whether the thread pointer points to the head of the threads array */
    return (tmpi_bool)(tMPI_Get_current() == threads);
}

tMPI_Comm tMPI_Get_comm_self(void)
{
    struct tmpi_thread* th = tMPI_Get_current();
    return th->self_comm;
}


int tMPI_Get_N(int *argc, char ***argv, const char *optname, int *nthreads)
{
    int i;
    int ret = TMPI_SUCCESS;

    *nthreads = 0;
    if (!optname)
    {
        i = 0;
    }
    else
    {
        for (i = 1; i < *argc; i++)
        {
            if (strcmp(optname, (*argv)[i]) == 0)
            {
                break;
            }
        }
    }
    if (i+1 < (*argc))
    {
        /* the number of processes is an argument */
        char *end;
        *nthreads = strtol((*argv)[i+1], &end, 10);
        if (!end || (*end != 0) )
        {
            *nthreads = 0;
            ret       = TMPI_FAILURE;
        }
    }
    if (*nthreads < 1)
    {
        int nth = tMPI_Thread_get_hw_number();

        if (nth < 1)
        {
            nth = 1;      /* make sure it's at least 1 */
        }
        *nthreads = nth;
    }

    return ret;
}

static int tMPI_Thread_init(struct tmpi_thread *th)
{
    int ret;
    int N_envelopes      = (Nthreads+1)*N_EV_ALLOC;
    int N_send_envelopes = N_EV_ALLOC;
    int N_reqs           = (Nthreads+1)*N_EV_ALLOC;
    int i;

    /* we set our thread id, as a thread-specific piece of global data. */
    ret = tMPI_Thread_setspecific(id_key, th);
    if (ret != 0)
    {
        return ret;
    }

    /* allocate comm.self */
    ret = tMPI_Comm_alloc( &(th->self_comm), TMPI_COMM_WORLD, 1);
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }
    th->self_comm->grp.peers[0] = th;

    /* allocate envelopes */
    ret = tMPI_Free_env_list_init( &(th->envelopes), N_envelopes );
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }
    /* recv list */
    ret = tMPI_Recv_env_list_init( &(th->evr));
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }
    /* send lists */
    th->evs = (struct send_envelope_list*)tMPI_Malloc(
                sizeof(struct send_envelope_list)*Nthreads);
    if (th->evs == NULL)
    {
        return TMPI_ERR_NO_MEM;
    }
    for (i = 0; i < Nthreads; i++)
    {
        ret = tMPI_Send_env_list_init( &(th->evs[i]), N_send_envelopes);
        if (ret != TMPI_SUCCESS)
        {
            return ret;
        }
    }

    tMPI_Atomic_set( &(th->ev_outgoing_received), 0);

    tMPI_Event_init( &(th->p2p_event) );

    /* allocate requests */
    ret = tMPI_Req_list_init(&(th->rql), N_reqs);
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }


#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* allcate copy_buffer list */
    ret = tMPI_Copy_buffer_list_init(&(th->cbl_multi),
                                     (Nthreads+1)*(N_COLL_ENV+1),
                                     Nthreads*COPY_BUFFER_SIZE);
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }
#endif

#ifdef TMPI_PROFILE
    ret = tMPI_Profile_init(&(th->profile));
    if (ret != TMPI_SUCCESS)
    {
        return ret;
    }
#endif
    /* now wait for all other threads to come on line, before we
       start the MPI program */
    ret = tMPI_Thread_barrier_wait( &(tmpi_global->barrier) );
    if (ret != 0)
    {
        return ret;;
    }
    return ret;
}


static void tMPI_Thread_destroy(struct tmpi_thread *th)
{
    int i;

    tMPI_Recv_env_list_destroy( &(th->evr));
    for (i = 0; i < Nthreads; i++)
    {
        tMPI_Send_env_list_destroy( &(th->evs[i]));
    }
    free(th->evs);
    tMPI_Free_env_list_destroy( &(th->envelopes) );
    tMPI_Event_destroy( &(th->p2p_event) );
    tMPI_Req_list_destroy( &(th->rql) );

#ifdef USE_COLLECTIVE_COPY_BUFFER
    tMPI_Copy_buffer_list_destroy(&(th->cbl_multi));
#endif

    for (i = 0; i < th->argc; i++)
    {
        free(th->argv[i]);
    }
}

static int tMPI_Global_init(struct tmpi_global *g, int Nthreads)
{
    int ret;

    g->usertypes        = NULL;
    g->N_usertypes      = 0;
    g->Nalloc_usertypes = 0;
    ret                 = tMPI_Thread_mutex_init(&(g->timer_mutex));
    if (ret != 0)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
    }
    tMPI_Spinlock_init(&(g->datatype_lock));

    ret = tMPI_Thread_barrier_init( &(g->barrier), Nthreads);
    if (ret != 0)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
    }

    ret = tMPI_Thread_mutex_init(&(g->comm_link_lock));
    if (ret != 0)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
    }


#if !(defined( _WIN32 ) || defined( _WIN64 ) )
    /* the time at initialization. */
    gettimeofday( &(g->timer_init), NULL);
#else
    /* the time at initialization. */
    g->timer_init = GetTickCount();
#endif
    return TMPI_SUCCESS;
}

static void tMPI_Global_destroy(struct tmpi_global *g)
{
    tMPI_Thread_barrier_destroy(&(g->barrier));
    tMPI_Thread_mutex_destroy(&(g->timer_mutex));
    tMPI_Thread_mutex_destroy(&(g->comm_link_lock));
}




static void* tMPI_Thread_starter(void *arg)
{
    int                 ret;
    struct tmpi_thread *th = (struct tmpi_thread*)arg;

#ifdef TMPI_TRACE
    tMPI_Trace_print("Created thread nr. %d", (int)(th-threads));
#endif

    ret = tMPI_Thread_init(th);
    if (ret != TMPI_SUCCESS)
    {
        return NULL;
    }

    /* start_fn, start_arg, argc and argv were set by the calling function */
    if (!th->start_fn)
    {
        th->start_fn_main(th->argc, th->argv);
    }
    else
    {
        th->start_fn(th->start_arg);
        if (!tmpi_finalized)
        {
            tMPI_Finalize();
        }
    }

    return NULL;
}


int tMPI_Start_threads(tmpi_bool main_returns, int N,
                       tMPI_Affinity_strategy aff_strategy,
                       int *argc, char ***argv,
                       void (*start_fn)(void*), void *start_arg,
                       int (*start_fn_main)(int, char**))
{
    int ret;
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Start_threads(%d, %d, %d, %d, %d, %p, %p, %p, %p)",
                     main_returns, N, aff_strategy, argc, argv, start_fn,
                     start_arg);
#endif
    if (N > 0)
    {
        int i;
        int set_affinity = FALSE;

        tmpi_finalized = FALSE;
        Nthreads       = N;

        /* allocate global data */
        tmpi_global = (struct tmpi_global*)
            tMPI_Malloc(sizeof(struct tmpi_global));
        if (tmpi_global == 0)
        {
            return TMPI_ERR_NO_MEM;
        }
        ret = tMPI_Global_init(tmpi_global, N);
        if (ret != TMPI_SUCCESS)
        {
            return ret;
        }

        /* allocate world and thread data */
        threads = (struct tmpi_thread*)
            tMPI_Malloc(sizeof(struct tmpi_thread)*N);
        if (threads == NULL)
        {
            return TMPI_ERR_NO_MEM;
        }
        ret = tMPI_Comm_alloc(&TMPI_COMM_WORLD, NULL, N);
        if (ret != TMPI_SUCCESS)
        {
            return ret;
        }
        TMPI_GROUP_EMPTY = tMPI_Group_alloc();

        if (tMPI_Thread_key_create(&id_key, NULL))
        {
            return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_INIT);
        }
        for (i = 0; i < N; i++)
        {
            TMPI_COMM_WORLD->grp.peers[i] = &(threads[i]);

            /* copy argc, argv */
            if (argc && argv)
            {
                int j;
                threads[i].argc = *argc;
                threads[i].argv = (char**)tMPI_Malloc(threads[i].argc*
                                                      sizeof(char*));
                for (j = 0; j < threads[i].argc; j++)
                {
#if !(defined( _WIN32 ) || defined( _WIN64 ) )
                    threads[i].argv[j] = strdup( (*argv)[j] );
#else
                    threads[i].argv[j] = _strdup( (*argv)[j] );
#endif
                }
            }
            else
            {
                threads[i].argc = 0;
                threads[i].argv = NULL;
            }
            threads[i].start_fn      = start_fn;
            threads[i].start_fn_main = start_fn_main;
            threads[i].start_arg     = start_arg;
        }

        /* now check whether to set affinity */
        if (aff_strategy == TMPI_AFFINITY_ALL_CORES)
        {
            int nhw = tMPI_Thread_get_hw_number();
            if ((nhw > 1) && (nhw == N))
            {
                set_affinity = TRUE;
            }
        }

        /* set thread 0's properties */
        threads[0].thread_id = tMPI_Thread_self();
        if (set_affinity)
        {
            /* set the main thread's affinity */
            tMPI_Thread_setaffinity_single(threads[0].thread_id, 0);
        }

        for (i = 1; i < N; i++) /* zero is the main thread */
        {
            ret = tMPI_Thread_create(&(threads[i].thread_id),
                                     tMPI_Thread_starter,
                                     (void*)&(threads[i]) );

            if (set_affinity)
            {
                tMPI_Thread_setaffinity_single(threads[i].thread_id, i);
            }
            if (ret != TMPI_SUCCESS)
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_INIT);
            }
        }
        /* the main thread also runs start_fn if we don't want
           it to return */
        if (!main_returns)
        {
            tMPI_Thread_starter((void*)&(threads[0]));

        }
        else
        {
            ret = tMPI_Thread_init(&(threads[0]));
            if (ret != 0)
            {
                return ret;
            }
        }
    }
    return TMPI_SUCCESS;
}


int tMPI_Init(int *argc, char ***argv,
              int (*start_function)(int, char**))
{
    int ret;
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Init(%p, %p, %p)", argc, argv, start_function);
#endif

    if (TMPI_COMM_WORLD == 0) /* we're the main process */
    {
        int N = 0;
        tMPI_Get_N(argc, argv, "-nt", &N);
        ret = tMPI_Start_threads(TRUE, N, TMPI_AFFINITY_ALL_CORES, argc, argv,
                                 NULL, NULL, start_function) != 0;
        if (ret != 0)
        {
            return ret;
        }
    }
    else
    {
        /* if we're a sub-thread we need don't need to do anyhing, because
           everything has already been set up by either the main thread,
           or the thread runner function.*/
    }
    return TMPI_SUCCESS;
}




int tMPI_Init_fn(int main_thread_returns, int N,
                 tMPI_Affinity_strategy aff_strategy,
                 void (*start_function)(void*), void *arg)
{
    int ret;
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Init_fn(%d, %p, %p)", N, start_function, arg);
#endif

    if (N < 1)
    {
        N = tMPI_Thread_get_hw_number();
        if (N < 1)
        {
            N = 1;    /*because that's what the fn returns if it doesn't know*/
        }
    }

    if (TMPI_COMM_WORLD == 0 && N >= 1) /* we're the main process */
    {
        ret = tMPI_Start_threads(main_thread_returns, N, aff_strategy,
                                 0, 0, start_function, arg, NULL);
        if (ret != 0)
        {
            return ret;
        }
    }
    return TMPI_SUCCESS;
}

int tMPI_Initialized(int *flag)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Initialized(%p)", flag);
#endif

    *flag = (TMPI_COMM_WORLD && !tmpi_finalized);

    return TMPI_SUCCESS;
}

int tMPI_Finalize(void)
{
    int i;
    int ret;
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Finalize()");
#endif
#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Finalize called\n", tMPI_This_threadnr());
    fflush(stdout);
#endif

#ifdef TMPI_PROFILE
    {
        struct tmpi_thread *cur = tMPI_Get_current();

        tMPI_Profile_stop( &(cur->profile) );
        ret = tMPI_Thread_barrier_wait( &(tmpi_global->barrier) );
        if (ret != 0)
        {
            return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
        }

        if (tMPI_Is_master())
        {
            tMPI_Profiles_summarize(Nthreads, threads);
        }
    }
#endif
    ret = tMPI_Thread_barrier_wait( &(tmpi_global->barrier) );
    if (ret != 0)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
    }



    if (tMPI_Is_master())
    {

        /* we just wait for all threads to finish; the order isn't very
           relevant, as all threads should arrive at their endpoints soon. */
        for (i = 1; i < Nthreads; i++)
        {
            if (tMPI_Thread_join(threads[i].thread_id, NULL))
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_FINALIZE);
            }
            tMPI_Thread_destroy(&(threads[i]));
        }
        /* at this point, we are the only thread left, so we can
           destroy the global structures with impunity. */
        tMPI_Thread_destroy(&(threads[0]));
        free(threads);

        tMPI_Thread_key_delete(id_key);
        /* de-allocate all the comm stuctures. */
        {
            tMPI_Comm cur;

            ret = tMPI_Thread_mutex_lock(&(tmpi_global->comm_link_lock));
            if (ret != 0)
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
            }
            cur = TMPI_COMM_WORLD->next;
            while (cur && (cur != TMPI_COMM_WORLD) )
            {
                tMPI_Comm next = cur->next;
                ret = tMPI_Comm_destroy(cur, FALSE);
                if (ret != 0)
                {
                    tMPI_Thread_mutex_unlock(&(tmpi_global->comm_link_lock));
                    return ret;
                }
                cur = next;
            }
            ret = tMPI_Comm_destroy(TMPI_COMM_WORLD, FALSE);
            if (ret != 0)
            {
                tMPI_Thread_mutex_unlock(&(tmpi_global->comm_link_lock));
                return ret;
            }
            ret = tMPI_Thread_mutex_unlock(&(tmpi_global->comm_link_lock));
            if (ret != 0)
            {
                return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_IO);
            }

        }

        tMPI_Group_free(&TMPI_GROUP_EMPTY);
        threads          = 0;
        TMPI_COMM_WORLD  = NULL;
        TMPI_GROUP_EMPTY = NULL;
        Nthreads         = 0;

        /* deallocate the 'global' structure */
        tMPI_Global_destroy(tmpi_global);
        free(tmpi_global);

        tmpi_finalized = TRUE;
    }
    else
    {
        tMPI_Thread_exit(0);
    }
    return TMPI_SUCCESS;
}


int tMPI_Finalized(int *flag)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Finalized(%p)", flag);
#endif
    *flag = tmpi_finalized;

    return TMPI_SUCCESS;
}



int tMPI_Abort(tMPI_Comm comm, int errorcode)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Abort(%p, %d)", comm, errorcode);
#endif
#if 0
    /* we abort(). This way we can run a debugger on it */
    fprintf(stderr, "tMPI_Abort called with error code %d", errorcode);
    if (comm == TMPI_COMM_WORLD)
    {
        fprintf(stderr, " on TMPI_COMM_WORLD");
    }
    fprintf(stderr, "\n");
    fflush(stdout);

    abort();
#else
    /* we just kill all threads, but not the main process */

    if (tMPI_Is_master())
    {
        if (comm == TMPI_COMM_WORLD)
        {
            fprintf(stderr,
                    "tMPI_Abort called on TMPI_COMM_WORLD main with errorcode=%d\n",
                    errorcode);
        }
        else
        {
            fprintf(stderr,
                    "tMPI_Abort called on main thread with errorcode=%d\n",
                    errorcode);
        }
        fflush(stderr);
        exit(errorcode);
    }
    else
    {
        int *ret;
        /* kill myself */
        fprintf(stderr, "tMPI_Abort called with error code %d on thread %d\n",
                errorcode, tMPI_This_threadnr());
        fflush(stderr);
        ret = (int*)malloc(sizeof(int));
        tMPI_Thread_exit(ret);
    }
#endif
    return TMPI_SUCCESS;
}


int tMPI_Get_processor_name(char *name, int *resultlen)
{
    int                nr     = tMPI_Threadnr(tMPI_Get_current());
    unsigned int       digits = 0;
    const unsigned int base   = 10;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Get_processor_name(%p, %p)", name, resultlen);
#endif
    /* we don't want to call sprintf here (it turns out to be not entirely
       thread-safe on Mac OS X, for example), so we do it our own way: */

    /* first determine number of digits */
    {
        int rest = nr;
        while (rest > 0)
        {
            rest /= base;
            digits++;
        }
        if (digits == 0)
        {
            digits = 1;
        }
    }
#ifndef _MSC_VER
    strcpy(name, "thread #");
#else
    strncpy_s(name, TMPI_MAX_PROCESSOR_NAME, "thread #", TMPI_MAX_PROCESSOR_NAME);
#endif
    /* now construct the number */
    {
        size_t       len = strlen(name);
        unsigned int i;
        int          rest = nr;

        for (i = 0; i < digits; i++)
        {
            size_t pos = len + (digits-i-1);
            if (pos < (TMPI_MAX_PROCESSOR_NAME -1) )
            {
                name[ pos ] = (char)('0' + rest%base);
            }
            rest /= base;
        }
        if ( (digits+len) < TMPI_MAX_PROCESSOR_NAME)
        {
            name[digits + len] = '\0';
        }
        else
        {
            name[TMPI_MAX_PROCESSOR_NAME] = '\0';
        }

    }
    if (resultlen)
    {
        *resultlen = (int)strlen(name); /* For some reason the MPI standard
                                           uses ints instead of size_ts for
                                           sizes. */
    }
    return TMPI_SUCCESS;
}





/* TODO: there must be better ways to do this */
double tMPI_Wtime(void)
{
    double ret = 0;

#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Wtime()");
#endif

#if !(defined( _WIN32 ) || defined( _WIN64 ) )
    {
        struct timeval tv;
        long int       secdiff;
        int            usecdiff;

        gettimeofday(&tv, NULL);
        secdiff  = tv.tv_sec - tmpi_global->timer_init.tv_sec;
        usecdiff = tv.tv_usec - tmpi_global->timer_init.tv_usec;

        ret = (double)secdiff + 1e-6*usecdiff;
    }
#else
    {
        DWORD tv = GetTickCount();

        /* the windows absolute time GetTickCount() wraps around in ~49 days,
           so it's safer to always use differences, and assume that our
           program doesn't run that long.. */
        ret = 1e-3*((unsigned int)(tv - tmpi_global->timer_init));
    }
#endif
    return ret;
}

double tMPI_Wtick(void)
{
#if !(defined( _WIN32 ) || defined( _WIN64 ) )
    /* In Unix, we don't really know. Any modern OS should be at least
       this precise, though */
    return 1./100.;
#else
    /* According to the Windows documentation, this is about right: */
    return 1./100.;
#endif
}

int tMPI_Get_count(tMPI_Status *status, tMPI_Datatype datatype, int *count)
{
#ifdef TMPI_TRACE
    tMPI_Trace_print("tMPI_Get_count(%p, %p, %p)", status, datatype, count);
#endif
    if (!status)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_STATUS);
    }
    *count = (int)(status->transferred/datatype->size);
    return TMPI_SUCCESS;
}
