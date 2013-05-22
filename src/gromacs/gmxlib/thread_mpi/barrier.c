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
#include <stdarg.h>
#include <string.h>

#include "impl.h"



void tMPI_Barrier_init(tMPI_Barrier_t *barrier, int count)
{
    barrier->threshold = count;
    tMPI_Atomic_set(&(barrier->cycle), 0);
    tMPI_Atomic_set(&(barrier->count), count);
    TMPI_YIELD_WAIT_DATA_INIT(barrier);
}


int tMPI_Barrier_wait(tMPI_Barrier_t *barrier)
{
    int cycle;
    int status;

    /* We don't need to lock or use atomic ops here, since the cycle index
     * cannot change until after the last thread has performed the check
     * further down. Further, they cannot reach this point in the next
     * barrier iteration until all of them have been released, and that
     * happens after the cycle value has been updated.
     *
     * No synchronization == fast synchronization.
     */
    cycle = tMPI_Atomic_get( &(barrier->cycle) );

    /* Decrement the count atomically and check if it is zero.
     * This will only be true for the last thread calling us.
     */
    if (tMPI_Atomic_fetch_add( &(barrier->count), -1 ) <= 1)
    {
        tMPI_Atomic_memory_barrier();
        tMPI_Atomic_set(&(barrier->count), barrier->threshold);
        tMPI_Atomic_fetch_add(&(barrier->cycle), 1);

        status = -1;
    }
    else
    {
        /* Wait until the last thread changes the cycle index.
         * We are both using a memory barrier, and explicit
         * volatile pointer cast to make sure the compiler
         * doesn't try to be smart and cache the contents.
         */
        do
        {
            /*tMPI_Atomic_memory_barrier();*/
            TMPI_YIELD_WAIT(barrier);
        }
        while (tMPI_Atomic_get( &(barrier->cycle) ) == cycle);
        tMPI_Atomic_memory_barrier();

        status = 0;
    }
    return status;
}
