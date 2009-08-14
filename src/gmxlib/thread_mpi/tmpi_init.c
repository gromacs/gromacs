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
 * Note that this could also be controlled using preprocessor flags,
 * which is the method used for cmake */
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
#ifndef _WIN32
#include <sys/time.h>
#endif

#include "thread_mpi.h"
#include "tmpi_impl.h"






struct tmpi_errhandler_ tmpi_errors_are_fatal = { 0, tmpi_errors_are_fatal_fn };
struct tmpi_errhandler_ tmpi_errors_return = { 0, tmpi_errors_return_fn };


tMPI_Errhandler TMPI_ERRORS_ARE_FATAL=&tmpi_errors_are_fatal;
tMPI_Errhandler TMPI_ERRORS_RETURN=&tmpi_errors_return;




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

/* this is where all the tMPI_Reduce ops are included from thread_tmpi_ops.c */
#define THREAD_MPI_OPS 1

#define TYPE char
#define TYPENM CHAR
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE short
#define TYPENM SHORT
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE int
#define TYPENM INT
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE long
#define TYPENM LONG
#define INTTYPE 1
#include "tmpi_ops.c"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE long long
#define TYPENM L_LONG
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE long long int
#define TYPENM L_L_INT
#define INTTYPE 1
#include "tmpi_ops.c"

#endif

#define TYPE signed char
#define TYPENM S_CHAR
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE unsigned char
#define TYPENM U_CHAR
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE unsigned short
#define TYPENM U_SHORT
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE unsigned 
#define TYPENM UNSIGNED
#define INTTYPE 1
#include "tmpi_ops.c"

#define TYPE unsigned long
#define TYPENM U_LONG
#define INTTYPE 1
#include "tmpi_ops.c"

#ifdef SIZEOF_LONG_LONG_INT

#define TYPE unsigned long long
#define TYPENM U_L_LONG
#define INTTYPE 1
#include "tmpi_ops.c"

#endif

#define TYPE float
#define TYPENM FLOAT
#define INTTYPE 0
#include "tmpi_ops.c"

#define TYPE double
#define TYPENM DOUBLE
#define INTTYPE 0
#include "tmpi_ops.c"

#define TYPE long double
#define TYPENM L_DOUBLE
#define INTTYPE 0
#include "tmpi_ops.c"

#define TYPE char
#define TYPENM BYTE
#define INTTYPE 1
#include "tmpi_ops.c"


tmpi_dt tmpi_char    ={sizeof(char),              oplist_CHAR,     0,NULL,TRUE};
tmpi_dt tmpi_short   ={sizeof(short),             oplist_SHORT,    0,NULL,TRUE};
tmpi_dt tmpi_int     ={sizeof(int),               oplist_INT,      0,NULL,TRUE};
tmpi_dt tmpi_long    ={sizeof(long),              oplist_LONG,     0,NULL,TRUE};
#ifdef SIZEOF_LONG_LONG_INT
tmpi_dt tmpi_l_long  ={sizeof(long long),         oplist_L_LONG,   0,NULL,TRUE};
tmpi_dt tmpi_l_l_int ={sizeof(long long int),     oplist_L_L_INT,  0,NULL,TRUE};
#endif
tmpi_dt tmpi_s_char  ={sizeof(signed char),       oplist_S_CHAR,   0,NULL,TRUE};
tmpi_dt tmpi_u_char  ={sizeof(unsigned char),     oplist_U_CHAR,   0,NULL,TRUE};
tmpi_dt tmpi_u_short ={sizeof(unsigned short),    oplist_U_SHORT,  0,NULL,TRUE};
tmpi_dt tmpi_unsigned={sizeof(unsigned),          oplist_UNSIGNED, 0,NULL,TRUE};
tmpi_dt tmpi_u_long  ={sizeof(unsigned long),     oplist_U_LONG,   0,NULL,TRUE};
#ifdef SIZEOF_LONG_LONG_INT
tmpi_dt tmpi_u_l_long={sizeof(unsigned long long),oplist_U_L_LONG, 0,NULL,TRUE};
#endif
tmpi_dt tmpi_float   ={sizeof(float),             oplist_FLOAT,    0,NULL,TRUE};
tmpi_dt tmpi_double  ={sizeof(double),            oplist_DOUBLE,   0,NULL,TRUE};
tmpi_dt tmpi_l_double={sizeof(long double),       oplist_L_DOUBLE, 0,NULL,TRUE};
tmpi_dt tmpi_byte    ={sizeof(char),              oplist_CHAR,     0,NULL,TRUE};




tMPI_Datatype TMPI_CHAR               = &tmpi_char;
tMPI_Datatype TMPI_SHORT              = &tmpi_short;
tMPI_Datatype TMPI_INT                = &tmpi_int;
tMPI_Datatype TMPI_LONG               = &tmpi_long;
#ifdef SIZEOF_LONG_LONG_INT
tMPI_Datatype TMPI_LONG_LONG          = &tmpi_l_long;
tMPI_Datatype TMPI_LONG_LONG_INT      = &tmpi_l_l_int;
#endif
tMPI_Datatype TMPI_SIGNED_CHAR        = &tmpi_s_char;
tMPI_Datatype TMPI_UNSIGNED_CHAR      = &tmpi_u_char;
tMPI_Datatype TMPI_UNSIGNED_SHORT     = &tmpi_u_short;
tMPI_Datatype TMPI_UNSIGNED           = &tmpi_unsigned;
tMPI_Datatype TMPI_UNSIGNED_LONG      = &tmpi_u_long;
#ifdef SIZEOF_LONG_LONG_INT
tMPI_Datatype TMPI_UNSIGNED_LONG_LONG = &tmpi_u_l_long;
#endif

tMPI_Datatype TMPI_FLOAT              = &tmpi_float;
tMPI_Datatype TMPI_DOUBLE             = &tmpi_double;
tMPI_Datatype TMPI_LONG_DOUBLE        = &tmpi_l_double;

/*extern tMPI_Datatype tMPI_UNSIGNED_WCHAR*/
tMPI_Datatype TMPI_BYTE               = &tmpi_byte;






/* error messages. Must match error codes in thread_mpi.h */
static const char *tmpi_errmsg[] =
{
    "No error",
    "malloc failure in tMPI (out of memory)",
    "tMPI Initialization error",
    "tMPI Finalize error",
    "Invalid tMPI_Group",
    "Invalid tMPI_Comm",
    "Invalid tMPI_Status",
    "Invalid tMPI_Group rank",
    "Invalid Cartesian topology dimensions",
    "Invalid Cartesian topology coordinates",
    "Insufficient number processes for Cartesian topology",
    "Invalid counterpart for MPI transfer",
    "Receive buffer size too small for transmission",
    "Overlapping send/receive buffers: probably due to thread-unsafe code.",
    "Invalid send destination",
    "Invalid receive source",
    "Invalid buffer (null pointer in send or receive buffer)",
    "Multicast operation mismatch (multicast not collective across comm)",
    "Invalid reduce operator",
    "Out of receive envelopes: this shouldn't happen (probably a bug).",
    "Out of receive requests: this shouldn't happen (probably a bug).",
    "Transmission failure",
    "Unknown tMPI error"
};






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

int tMPI_Error(tMPI_Comm comm, int tmpi_errno)
{
    if (comm)
    {
        comm->erh->err=tmpi_errno;
        comm->erh->fn(&comm, &tmpi_errno);
    }
    else
    {
        /* initialization errors have no comm */
        tmpi_errors_are_fatal_fn(NULL, &tmpi_errno);
    }
    return tmpi_errno;
}


int tMPI_Error_string(int errorcode, char *strn, size_t *resultlen)
{
    if (errorcode<0 || errorcode>=N_TMPI_ERR)
        errorcode=TMPI_ERR_UNKNOWN;

#ifndef _WIN32
    strncpy(strn, tmpi_errmsg[errorcode], TMPI_MAX_ERROR_STRING);
#else
    strncpy_s(strn, TMPI_MAX_ERROR_STRING, tmpi_errmsg[errorcode], TMPI_MAX_ERROR_STRING);
#endif
    *resultlen=strlen(strn);
    return TMPI_SUCCESS;
}

int tMPI_Create_errhandler(tMPI_Errhandler_fn *function, 
                           tMPI_Errhandler *errhandler) 
{
    /* we don't use a special malloc here because this is the error handler
       creation function. */
    *errhandler=(tMPI_Errhandler)malloc(sizeof(struct tmpi_errhandler_));
    if (!*errhandler)
    {
        fprintf(stderr, "tMPI fatal error (%s), bailing out\n", 
                tmpi_errmsg[TMPI_ERR_MALLOC]);
        abort();
    }
    (*errhandler)->err=0;
    (*errhandler)->fn=*function;
    return TMPI_SUCCESS;
}

int tMPI_Errhandler_free(tMPI_Errhandler *errhandler)
{
    free(*errhandler);
    return TMPI_SUCCESS;
}


int tMPI_Comm_set_errhandler(tMPI_Comm comm, tMPI_Errhandler errhandler)
{
    comm->erh = errhandler;
    return TMPI_SUCCESS;
}

int tMPI_Comm_get_errhandler(tMPI_Comm comm, tMPI_Errhandler *errhandler)
{
    *errhandler=comm->erh;
    return TMPI_SUCCESS;
}

void tmpi_errors_are_fatal_fn(tMPI_Comm *comm, int *err)
{
    char errstr[TMPI_MAX_ERROR_STRING];
    size_t len;

    tMPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "MPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "MPI error: %s\n", errstr);
    }
    abort();
    /*exit(0);*/
}

void tmpi_errors_return_fn(tMPI_Comm *comm, int *err)
{
    char errstr[TMPI_MAX_ERROR_STRING];
    size_t len;

    tMPI_Error_string(*err, errstr, &len);
    if (comm)
    {
        fprintf(stderr, "MPI error: %s (in valid comm)\n", errstr);
    }
    else
    {
        fprintf(stderr, "MPI error: %s\n", errstr);
    }
    return;
}

int tMPI_Get_N(int *argc, char ***argv)
{
    int i=0;
    int np=1;

    for(i=1;i<*argc;i++)
    {
        if (strcmp("-np", (*argv)[i]) == 0)
        {
            if (i+1 < (*argc))
            {
                /* the number of processes is an argument */
                char *end;
                np=strtol((*argv)[i+1], &end, 0);
                if ( !end || (*end != 0) )
                    np=1;
            }
            break;
        }
    }
    if (np<1)
        np=1;
    return np;
}


static void* tMPI_Thread_starter(void *arg)
{
    struct tmpi_thread *th=(struct tmpi_thread*)arg;
    int N_envelopes=(Nthreads+1)*(Nthreads+1)*8;  /*AARGH arbitrary number*/
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
        tMPI_Send_env_list_init( &(th->evs[i]));
    }
    tMPI_Atomic_set( &(th->evs_check_id), 0);

    /* allocate requests */
    tMPI_Req_list_init(&(th->rql), N_reqs);

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
        tmpi_global=(struct tmpi_global*)tMPI_Malloc(sizeof(struct tmpi_global));
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
#ifndef _WIN32
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
        int N=tMPI_Get_N(argc, argv);
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
    /* we abort(). This way we can run a debugger on it */
    fprintf(stderr, "tMPI_Abort called with error code %d",errorcode);
    if (comm==TMPI_COMM_WORLD)
        fprintf(stderr, " on TMPI_COMM_WORLD");
    fprintf(stderr,"\n");
    fflush(0);

    abort();
#if 0
    /* we just kill all threads, but not the main process */
    int i;
    struct tmpi_thread *me=tMPI_Get_current();
    /* kill all threads */
    for(i=0;i<comm->grp.N;i++)
    {
        if (comm->grp.peers[i] != me && threads[i].thread_id)
            tMPI_thread_cancel(threads[i].thread_id);
    }
    /* kill myself */
    if (me->thread_id)
        tMPI_thread_cancel(me->thread_id);
#endif
    return TMPI_SUCCESS;
}


int tMPI_Get_processor_name(char *name, size_t *resultlen)
{
    threadnr_t nr=tMPI_Threadnr(tMPI_Get_current());
    unsigned int digits=0;
    const unsigned int base=10;

    /* we don't want to call sprintf here (it turns out to be not entirely
       thread-safe on Mac OS X, for example), so we do it our own way: */

    /* first determine number of digits */
    {
        threadnr_t rest=nr;
        while(rest > 0)
        {
            rest /= base;
            digits++;
        }
        if (digits==0)
            digits=1;
    }
#ifndef _WIN32
    strcpy(name, "thread #");
#else
    strncpy_s(name, TMPI_MAX_PROCESSOR_NAME, "thread #", TMPI_MAX_PROCESSOR_NAME);
#endif
    /* now construct the number */
    {
        size_t len=strlen(name);
        unsigned int i;
        threadnr_t rest=nr;

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
        *resultlen=strlen(name);
    return TMPI_SUCCESS;
}

/* TODO: there must be better ways to do this */
double tMPI_Wtime(void)
{
    double ret=0;
#ifndef _WIN32
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
    *count = status->transferred;
    return TMPI_SUCCESS;
}



int tMPI_Type_contiguous(int count, tMPI_Datatype oldtype, 
                         tMPI_Datatype *newtype)
{
    struct tmpi_datatype_ *ntp;

    ntp=(struct tmpi_datatype_*)tMPI_Malloc(sizeof(struct tmpi_datatype_));
    ntp->size=count*oldtype->size;
    ntp->op_functions=NULL;

    /* establish components */
    ntp->N_comp=1;
    ntp->comps=(struct tmpi_datatype_component*)tMPI_Malloc(
                        sizeof(struct tmpi_datatype_component)*1);
    ntp->comps[0].type=oldtype;
    ntp->comps[0].count=1;
    ntp->committed=FALSE;

    /* now add it to the list.  */
    tMPI_Spinlock_lock(&(tmpi_global->datatype_lock));
    /* check whether there's space */
    if (tmpi_global->N_usertypes + 1 >= tmpi_global->Nalloc_usertypes)
    {
        /* make space */
        tmpi_global->Nalloc_usertypes=Nthreads*(tmpi_global->N_usertypes) + 1;
        tmpi_global->usertypes=tMPI_Realloc(tmpi_global->usertypes, 
                 (sizeof(struct tmpi_datatype_ *)*tmpi_global->Nalloc_usertypes));

    }
    /* add to the list */
    tmpi_global->usertypes[tmpi_global->N_usertypes]=ntp;
    tmpi_global->N_usertypes++;
    *newtype=ntp;
    tMPI_Spinlock_unlock(&(tmpi_global->datatype_lock));

    return TMPI_SUCCESS;
}


int tMPI_Type_commit(tMPI_Datatype *datatype)
{
    int i,j;
    struct tmpi_datatype_ *dt=*datatype;

    if (dt->committed)
        return TMPI_SUCCESS;

    /* search the list for a matching committed type, because if there's
       already a committed type that has the same composition, we just 
       make the datatype pointer point to it, ensuring we share datatype 
       information across threads. */
    tMPI_Spinlock_lock(&(tmpi_global->datatype_lock));
    for(i=0;i<tmpi_global->N_usertypes;i++)
    {
        struct tmpi_datatype_ *lt=tmpi_global->usertypes[i];
        if (lt->committed && lt->N_comp==dt->N_comp)
        {
            bool found=TRUE;
            for(j=0;j<lt->N_comp;j++)
            {
                if ( (lt->comps[j].type  != dt->comps[j].type) ||
                     (lt->comps[j].count != dt->comps[j].count) )
                {
                    found=FALSE;
                    break;
                }
            }
            if (found)
            {
                dt=lt;
            }
        }
    }
    if (dt != *datatype)
    {
        bool found=FALSE;
        /* we remove the old one from the list */
        for(i=0;i<tmpi_global->N_usertypes;i++)
        {
            if (tmpi_global->usertypes[i]==*datatype)
            {
                found=TRUE;
                break;
            }
        }
        if (found)
        {
            /* we put the last one in the list in our slot */
            tmpi_global->usertypes[i]=tmpi_global->
                usertypes[tmpi_global->N_usertypes-1];
            tmpi_global->N_usertypes--;
        }
        free( (*datatype)->comps);
        free(  *datatype );

        /* and overwrite the pointer with the new data type */
        *datatype=dt;
    }
    else
    {
        /* it was the first one of its type */
        dt->committed=TRUE;
    }
    tMPI_Spinlock_unlock(&(tmpi_global->datatype_lock));
    return TMPI_SUCCESS;
}

