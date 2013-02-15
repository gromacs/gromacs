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

/* the types that were defined in include/thread_mpi/threads.h */


struct tMPI_Thread
{
    HANDLE th;              /* the thread handle */
    DWORD  id;              /* the thread ID */
    int    started_by_tmpi; /* whether this thread was started by tmpi */
};

struct tMPI_Thread_key
{
    DWORD wkey;
};

struct tMPI_Mutex
{
    CRITICAL_SECTION cs;
};

struct tMPI_Thread_once
{
    int dum;
};

struct tMPI_Thread_cond
{
#if 0
    /* this works since Windows Vista: */
    CONDITION_VARIABLE cv;
#else
    /* this data structure and its algorithms are based on
       'Strategies for Implementing POSIX Condition Variables on Win32'
       by
       Douglas C. Schmidt and Irfan Pyarali
       Department of Computer Science
       Washington University, St. Louis, Missouri
       http://www.cs.wustl.edu/~schmidt/win32-cv-1.html */
    int              Nwaiters; /* number of waiting threads */
    CRITICAL_SECTION wtr_lock; /* lock for Nwaiters */
    int              Nrelease; /* number of threads to release in broadcast/signal */
    int              cycle;    /* cycle number so threads can't steal signals */
    HANDLE           ev;       /* the event used to trigger WaitForSingleObject.
                                  Is a manual reset event.  */
#endif
};

struct tMPI_Thread_barrier
{
#if 0
    /* use this once Vista is the oldest supported windows version: */
    CRITICAL_SECTION   cs;    /*!< Lock for the barrier
                                 contents          */
    CONDITION_VARIABLE cv;    /*!< Condition to signal barrier
                                 completion */
#else
    tMPI_Thread_mutex_t cs;   /*!< Lock for the barrier contents          */
    tMPI_Thread_cond_t  cv;   /*!< Condition to signal barrier completion */
#endif
};
