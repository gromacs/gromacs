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


#ifndef TMPI_WAIT_H_
#define TMPI_WAIT_H_

#ifndef TMPI_WAIT_FOR_NO_ONE

#if !(defined( _WIN32 ) || defined( _WIN64 ) )
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SCHED_H
#include <sched.h>
#endif

/* for now we just do sched_yield(). It's in POSIX. */
/* the data associated with waiting. */
#define TMPI_YIELD_WAIT_DATA
/* the initialization  associated with waiting. */
#define TMPI_YIELD_WAIT_DATA_INIT(data)

/* the waiting macro */
#define TMPI_YIELD_WAIT(data)  sched_yield()

#else
/* and in Windows, we do SwitchToThread() alternated with Sleep(0). This
   is apparently recommende practice (SwitchToThread() alone just gives
   up the slice for threads on the current core, and Sleep(0) alone could
   lead to starvation. This mixed approach actually gives better real-world
   performance in the test program.*/
/* the data associated with waiting. */
#define TMPI_YIELD_WAIT_DATA  int yield_wait_counter;
/* the initialization  associated with waiting. */
#define TMPI_YIELD_WAIT_DATA_INIT(data) { (data)->yield_wait_counter = 0; }

/* the waiting macro is so complicated because using SwitchToThread only schedules */
#define TMPI_YIELD_WAIT(data)  { \
        if ( ((data)->yield_wait_counter++)%100 == 0) \
        { \
            SwitchToThread(); \
        } \
        else \
        { \
            Sleep(0); \
        } \
}

#endif


#else /* !TMPI_WAIT_FOR_NO_ONE */

/* the data associated with waiting. */
#define TMPI_YIELD_WAIT_DATA
/* the initialization  associated with waiting. */
#define TMPI_YIELD_WAIT_DATA_INIT(data)

/* the waiting macro */
#define TMPI_YIELD_WAIT(data)  tMPI_Atomic_memory_barrier()


#endif /* !TMPI_WAIT_FOR_NO_ONE */

#endif /* TMPI_WAIT_H_ */
