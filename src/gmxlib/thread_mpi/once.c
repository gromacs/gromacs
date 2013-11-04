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


/* once */
int tMPI_Once(tMPI_Comm comm, void (*function)(void*), void *param,
              int *was_first)
{
    int               myrank;
    int               ret = TMPI_SUCCESS;
    struct coll_sync *csync;
    struct coll_env  *cev;
    int               syncs;


    if (!comm)
    {
        return tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
    }
    myrank = tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    /* we increase our counter, and determine which coll_env we get */
    csync = &(comm->csync[myrank]);
    csync->syncs++;
    cev = &(comm->cev[csync->syncs % N_COLL_ENV]);

    /* now do a compare-and-swap on the current_syncc */
    syncs = tMPI_Atomic_get( &(cev->coll.current_sync));
    if ((csync->syncs - syncs > 0) && /* check if sync was an earlier number.
                                         If it is a later number, we can't
                                         have been the first to arrive here. */
        tMPI_Atomic_cas(&(cev->coll.current_sync), syncs, csync->syncs))
    {
        /* we're the first! */
        function(param);
        if (was_first)
        {
            *was_first = TRUE;
        }
    }
    return ret;
}

void* tMPI_Once_wait(tMPI_Comm comm, void* (*function)(void*), void *param,
                     int *was_first)
{
    int               myrank;
    struct coll_sync *csync;
    struct coll_env  *cev;
    int               syncs;
    void             *ret;


    if (!comm)
    {
        tMPI_Error(TMPI_COMM_WORLD, TMPI_ERR_COMM);
        return NULL;
    }
    myrank = tMPI_Comm_seek_rank(comm, tMPI_Get_current());

    /* we increase our counter, and determine which coll_env we get */
    csync = &(comm->csync[myrank]);
    csync->syncs++;
    cev = &(comm->cev[csync->syncs % N_COLL_ENV]);

    /* now do a compare-and-swap on the current_syncc */
    syncs = tMPI_Atomic_get( &(cev->coll.current_sync));
    tMPI_Atomic_memory_barrier_acq();
    if ((csync->syncs - syncs > 0) && /* check if sync was an earlier number.
                                         If it is a later number, we can't
                                         have been the first to arrive here.
                                         Calculating the difference instead
                                         of comparing directly avoids ABA
                                         problems. */
        tMPI_Atomic_cas(&(cev->coll.current_sync), syncs, csync->syncs))
    {
        /* we're the first! */
        ret = function(param);
        if (was_first)
        {
            *was_first = TRUE;
        }

        /* broadcast the output data */
        cev->coll.res = ret;

        tMPI_Atomic_memory_barrier_rel();
        /* signal that we're done */
        tMPI_Atomic_fetch_add(&(cev->coll.current_sync), 1);
        /* we need to keep being in sync */
        csync->syncs++;
    }
    else
    {
        /* we need to wait until the current_syncc gets increased again */
        csync->syncs++;
        do
        {
            /*tMPI_Atomic_memory_barrier();*/
            syncs = tMPI_Atomic_get( &(cev->coll.current_sync) );
        }
        while (csync->syncs - syncs > 0);   /* difference again due to ABA
                                               problems */
        tMPI_Atomic_memory_barrier_acq();
        ret = cev->coll.res;
    }
    return ret;
}


static void *tMPI_Shmallocator(void *prm)
{
    size_t sz = *((size_t*)prm);
    return malloc(sz);
}

void *tMPI_Shmalloc(tMPI_Comm comm, size_t size)
{
    return tMPI_Once_wait(comm, tMPI_Shmallocator, &size, NULL);
}
