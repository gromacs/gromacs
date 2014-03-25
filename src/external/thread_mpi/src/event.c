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

#include "impl.h"

#ifdef TMPI_TRACE
#include <stdarg.h>
#endif


void tMPI_Event_init(tMPI_Event *ev)
{
    tMPI_Atomic_set(&(ev->sync), 0);
    ev->last_sync = 0;
}

void tMPI_Event_destroy(tMPI_Event *ev)
{
    tMPI_Atomic_set(&(ev->sync), 0);
    ev->last_sync = 0;
}

int tMPI_Event_wait(tMPI_Event *ev)
{
    int ret;
    /* for most OSes yielding waits result in much better performance
       (by an order of magnitude) than using the OS-provided wait functions
       such as pthread_cond_wait(). That's why we do a busy-wait loop here.*/
    while ( (tMPI_Atomic_get(&(ev->sync)) - (ev->last_sync)) <= 0)
    {
        TMPI_YIELD_WAIT(ev);
    }
    tMPI_Atomic_memory_barrier_acq();
    ret = tMPI_Atomic_get(&(ev->sync)) - (ev->last_sync);
    return ret;
}
