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
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
#include <sys/time.h>
#endif

#include "thread_mpi/threads.h"
#include "thread_mpi/atomic.h"
#include "thread_mpi/tmpi.h"
#include "tmpi_impl.h"







/* there are a few global variables that maintain information about the
   running threads. Some are defined by the MPI standard: */
tMPI_Comm TMPI_COMM_WORLD=NULL;
tMPI_Group tMPI_GROUP_EMPTY=NULL;


/* the threads themselves (tmpi_comm only contains lists of pointers to this
      structure */
struct tmpi_thread *threads=NULL;
int Nthreads=0;

/* thread info */
tMPI_Thread_key_t id_key; /* the key to get the thread id */



/* whether MPI has finalized (we need this to distinguish pre-inited from
       post-finalized states */
static bool tmpi_finalized=FALSE;

/* misc. global information about MPI */
struct tmpi_global *tmpi_global=NULL;








/* start N threads with argc, argv (used by tMPI_Init)*/
static void tMPI_Start_threads(int N, int *argc, char ***argv,
                               void (*start_fn)(void*), void *start_arg);
/* starter function for threads; takes a void pointer to a
      struct tmpi_starter_, which calls main() if tmpi_start_.fn == NULL */
static void* tMPI_Thread_starter(void *arg);








void *tMPI_Malloc(size_t size)
{
    void *ret=(void*)malloc(size);
    if (!ret)
    {
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_MALLOC);
    }
    return ret;
}

void *tMPI_Realloc(void *p, size_t size)
{
    void *ret=(void*)realloc(p, size);
    if (!ret)
    {
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_MALLOC);
    }
    return ret;
}


#if 0
struct tmpi_thread *tMPI_Get_current(void)
{
    if (!threads)
        return NULL;

    return (struct tmpi_thread*)tMPI_thread_getspecific(id_key);
} 


unsigned int tMPI_Threadnr(struct tmpi_thread *thr)
{
    return thr-threads;
}
#endif
#if 0
unsigned int tMPI_This_threadnr(void)
{
    return tMPI_Get_current()-threads;
}

struct tmpi_thread *tMPI_Get_thread(tMPI_Comm comm, int rank)
{
    /* check destination */
    if ( (rank < 0) || (rank > comm->grp.N) )
    {
        tMPI_Error(comm, TMPI_ERR_GROUP_RANK);
        return NULL;
    }
    return comm->grp.peers[rank];
}
#endif

bool tMPI_Is_master(void)
{
    /* if there are no other threads, we're the main thread */
    if ( (!TMPI_COMM_WORLD) || TMPI_COMM_WORLD->grp.N==0)
        return TRUE;

    /* otherwise we know this through thread specific data: */
    /* whether the thread pointer points to the head of the threads array */
    return (bool)(tMPI_Get_current() == threads); 
}

tMPI_Comm tMPI_Get_comm_self(void)
{
    struct tmpi_thread* th=tMPI_Get_current();
    return th->self_comm;
}


int tMPI_Get_N(int *argc, char ***argv, const char *optname, int *nthreads)
{
    int i;
    int ret=TMPI_FAILURE;

    *nthreads=1;
    if (!optname)
    {
        i=0;
    }
    else
    {
        for(i=1;i<*argc;i++)
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
        *nthreads=strtol((*argv)[i+1], &end, 10);
        if ( !end || (*end != 0) )
        {
            *nthreads=1;
        }
        else if (*nthreads > 0)
        {
            ret=TMPI_SUCCESS;
        }
    }
    if (*nthreads<1)
        *nthreads=1;
    return ret;
}


static void* tMPI_Thread_starter(void *arg)
{
    struct tmpi_thread *th=(struct tmpi_thread*)arg;
    int N_envelopes=(Nthreads+1)*(Nthreads+1)*8;  /*AARGH arbitrary number*/
    int N_send_envelopes=(Nthreads+1)*8;  /*AARGH arbitrary number*/
    int N_reqs=(Nthreads+1)*(Nthreads+1)*2;  /*AARGH arbitrary number*/
    int i;

    tMPI_Thread_setspecific(id_key, arg);


    /* allocate comm.self */
    th->self_comm=tMPI_Comm_alloc(TMPI_COMM_WORLD, 1);
    th->self_comm->grp.peers[0]=th;

    /* allocate envelopes */
    tMPI_Free_env_list_init( &(th->envelopes), N_envelopes );
    /* recv list */
    tMPI_Recv_env_list_init( &(th->evr));
    /* send lists */
    th->evs=(struct send_envelope_list*)tMPI_Malloc(
                        sizeof(struct send_envelope_list)*Nthreads);
    for(i=0;i<Nthreads;i++)
    {
        tMPI_Send_env_list_init( &(th->evs[i]), N_send_envelopes);
    }
    tMPI_Atomic_set( &(th->evs_check_id), 0);

    /* allocate requests */
    tMPI_Req_list_init(&(th->rql), N_reqs);

#ifdef USE_COLLECTIVE_COPY_BUFFER
    /* allcate copy_buffer list */
    tMPI_Copy_buffer_list_init(&(th->cbl_multi), 
                               (Nthreads+1)*(N_COLL_ENV+1),
                               Nthreads*COPY_BUFFER_SIZE);
#endif

    /* now wait for all other threads to come on line, before we
       start the MPI program */
    tMPI_Barrier(TMPI_COMM_WORLD);

    if (! th->start_fn )
        main(th->argc, th->argv);
    else
    {
        th->start_fn(th->start_arg);
        if (!tmpi_finalized)
            tMPI_Finalize();
    }

    /* TODO: freeing the thread structures and the mpi_comms still needs to 
       be implemented. */

    return 0;
}


void tMPI_Start_threads(int N, int *argc, char ***argv, 
                        void (*start_fn)(void*), void *start_arg)
{
    if (N>0) 
    {
        int i;

        tmpi_finalized=FALSE;
        Nthreads=N;

        /* allocate global data */
        tmpi_global=(struct tmpi_global*)
                        tMPI_Malloc(sizeof(struct tmpi_global));
        tmpi_global->usertypes=NULL;
        tmpi_global->N_usertypes=0;
        tmpi_global->Nalloc_usertypes=0;
        tMPI_Spinlock_init(&(tmpi_global->datatype_lock));

        /* allocate world and thread data */
        threads=(struct tmpi_thread*)tMPI_Malloc(sizeof(struct tmpi_thread)*N);
        TMPI_COMM_WORLD=tMPI_Comm_alloc(NULL, N);
        tMPI_GROUP_EMPTY=tMPI_Group_alloc();

        TMPI_COMM_WORLD->grp.N=N;

        if (tMPI_Thread_key_create(&id_key, NULL))
        {
            tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_INIT);
        }
        /*printf("thread keys created\n"); fflush(NULL);*/
        for(i=0;i<N;i++)
        {
            TMPI_COMM_WORLD->grp.peers[i]=&(threads[i]);

            /* copy argc, argv */
            if (argc && argv)
            {
                int j;
                threads[i].argc=*argc;
                threads[i].argv=(char**)tMPI_Malloc(threads[i].argc*
                                                   sizeof(char*));
                for(j=0;j<threads[i].argc;j++)
                {
#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
                    threads[i].argv[j]=strdup( (*argv)[j] );
#else
                    threads[i].argv[j]=_strdup( (*argv)[j] );
#endif
                }
            }
            else
            {
                threads[i].argc=0;
                threads[i].argv=NULL;
            }
            threads[i].start_fn=start_fn;
            threads[i].start_arg=start_arg;
        }
        for(i=1;i<N;i++) /* zero is the main thread */
        {
            if (tMPI_Thread_create(&(threads[i].thread_id), 
                                  tMPI_Thread_starter,
                                  (void*)&(threads[i]) ) )
            {
                tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_INIT);
            }
        }
        /* the main thread now also runs start_fn */
        /*threads[0].thread_id=NULL; we can't count on this being a pointer*/
        tMPI_Thread_starter((void*)&(threads[0]));
    }
}


int tMPI_Init(int *argc, char ***argv)
{
    if (TMPI_COMM_WORLD==0) /* we're the main process */
    {
        int N;
        tMPI_Get_N(argc, argv, "-np", &N);
        tMPI_Start_threads(N, argc, argv, NULL, NULL);
    }
    else
    {
        /* if we're a sub-thread we need don't need to do anyhing, because 
           everything has already been set up by either the main thread, 
           or the thread runner function.*/
    }
    return TMPI_SUCCESS;
}

int tMPI_Init_fn(int N, void (*start_function)(void*), void *arg)
{
    if (TMPI_COMM_WORLD==0 && N>=1) /* we're the main process */
    {
        tMPI_Start_threads(N, 0, 0, start_function, arg);
    }
    return TMPI_SUCCESS;
}

int tMPI_Initialized(int *flag)
{
    *flag=(TMPI_COMM_WORLD && !tmpi_finalized);

    return TMPI_SUCCESS;
}

int tMPI_Finalize(void)
{
    int i;
    struct tmpi_thread *th=tMPI_Get_current();

#ifdef TMPI_DEBUG
    printf("%5d: tMPI_Finalize called\n", tMPI_This_threadnr());
    fflush(stdout);
#endif

    tMPI_Barrier(TMPI_COMM_WORLD);

    for(i=0;i<(Nthreads+1);i++)
    {
        tMPI_Recv_env_list_destroy( &(th->evr));
    }
    for(i=0;i<Nthreads;i++)
    {
        tMPI_Send_env_list_destroy( &(th->evs[i]));
    }
    tMPI_Free_env_list_destroy( &(th->envelopes) );

    tMPI_Req_list_destroy( &(th->rql) );


    if (tMPI_Is_master())
    {
        /* we just wait for all threads to finish; the order isn't very 
           relevant, as all threads should arrive at their endpoints soon. */
        for(i=1;i<Nthreads;i++)
        {
            if (tMPI_Thread_join(threads[i].thread_id, NULL))
            {
                tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_FINALIZE);
            }
        }
        free(threads);
        free(TMPI_COMM_WORLD);
        free(tMPI_GROUP_EMPTY);
        threads=0;
        TMPI_COMM_WORLD=NULL;
        tMPI_GROUP_EMPTY=NULL;
        Nthreads=0;
        tmpi_finalized=TRUE;
    }
    else
    {
        tMPI_Thread_exit(0);
    }
    return TMPI_SUCCESS;
}

int tMPI_Finalized(int *flag)
{
    *flag=tmpi_finalized;

    return TMPI_SUCCESS;
}



int tMPI_Abort(tMPI_Comm comm, int errorcode)
{
#if 0
    /* we abort(). This way we can run a debugger on it */
    fprintf(stderr, "tMPI_Abort called with error code %d",errorcode);
    if (comm==TMPI_COMM_WORLD)
        fprintf(stderr, " on TMPI_COMM_WORLD");
    fprintf(stderr,"\n");
    fflush(0);

    abort();
#else
    /* we just kill all threads, but not the main process */
    
    if (tMPI_Is_master())
    {
        if (comm==TMPI_COMM_WORLD)
            fprintf(stderr, 
               "tMPI_Abort called on TMPI_COMM_WORLD main with errorcode=%d\n",
               errorcode);
        else
        fprintf(stderr, "tMPI_Abort called on main thread with errorcode=%d\n",
                errorcode);
        fflush(stderr);
        /*sleep(1);*/
        abort();
    }
    else
    {
        int *ret;
        /* kill myself */
        fprintf(stderr, "tMPI_Abort called wiht error code %d on thread %d\n", 
                        errorcode, tMPI_This_threadnr());
        fflush(stderr);
        ret=(int*)malloc(sizeof(int));
        tMPI_Thread_exit(ret);
    }
#endif
    return TMPI_SUCCESS;
}


int tMPI_Get_processor_name(char *name, int *resultlen)
{
    int nr=tMPI_Threadnr(tMPI_Get_current());
    unsigned int digits=0;
    const unsigned int base=10;

    /* we don't want to call sprintf here (it turns out to be not entirely
       thread-safe on Mac OS X, for example), so we do it our own way: */

    /* first determine number of digits */
    {
        int rest=nr;
        while(rest > 0)
        {
            rest /= base;
            digits++;
        }
        if (digits==0)
            digits=1;
    }
#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
    strcpy(name, "thread #");
#else
    strncpy_s(name, TMPI_MAX_PROCESSOR_NAME, "thread #", TMPI_MAX_PROCESSOR_NAME);
#endif
    /* now construct the number */
    {
        size_t len=strlen(name);
        unsigned int i;
        int rest=nr;

        for(i=0;i<digits;i++)
        {
            size_t pos=len + (digits-i-1);
            if (pos < (TMPI_MAX_PROCESSOR_NAME -1) )
                name[ pos ]=(char)('0' + rest%base);
            rest /= base;
        }
        if ( (digits+len) < TMPI_MAX_PROCESSOR_NAME)
            name[digits + len]='\0';
        else
            name[TMPI_MAX_PROCESSOR_NAME]='\0';

    }
    if (resultlen)
        *resultlen=(int)strlen(name); /* For some reason the MPI standard
                                         uses ints instead of size_ts for
                                         sizes. */
    return TMPI_SUCCESS;
}

/* TODO: there must be better ways to do this */
double tMPI_Wtime(void)
{
    double ret=0;
#if ! (defined( _WIN32 ) || defined( _WIN64 ) )
    struct timeval tv;
        
    gettimeofday(&tv, NULL);
    ret=tv.tv_sec + 1000000.*tv.tv_usec;
#else
    /* TODO: make timing function work */
#endif
    return ret;
}

#if 0
double tMPI_Wtick(void)
{
    /* this is just lying: */
    return 1./10000.;
}
#endif



int tMPI_Get_count(tMPI_Status *status, tMPI_Datatype datatype, int *count)
{
    if (!status)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_STATUS);
    }
    *count = (int)(status->transferred/datatype->size);
    return TMPI_SUCCESS;
}




