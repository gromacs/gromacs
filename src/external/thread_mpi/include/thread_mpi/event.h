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

#ifndef TMPI_EVENT_H_
#define TMPI_EVENT_H_

#include "visibility.h"
#include "wait.h"

/*! \file

   \brief Event notification wait and signaling functions.

   The event structure offers lightweight signaling and scheduler-yielding
   (but still spinning) waiting for waits of intermediate durations (i.e.
   longer than appropriate for a spin lock but shorter than mutex_cond_wait().
   These functions only take care of the waiting and signaling; the user is
   responsible for the handling of the actual event structures and condition
   variables.
 */

/*! \brief Event notification structure.

   Structure for notifying one thread that something has happened, allowing
   that thread to wait, and possibly give up its timeslice. This structure
   only takes care of the notification itself, and will not handle data
   about incoming events - that is left up to the user.

   This structure allows notification of a single thread by any number of
   threads*/
typedef struct tMPI_Event_t tMPI_Event;
struct tMPI_Event_t
{
    tMPI_Atomic_t sync;      /* the event sync counter */
    int           last_sync; /* the last sync event looked at */
    TMPI_YIELD_WAIT_DATA     /* data associated with yielding */
};



/*! \brief Initialize the event object.

    \param ev The event structure to be intialized. */
TMPI_EXPORT
void tMPI_Event_init(tMPI_Event *ev);

/*! \brief Deallocate the internals of the contents of event object.

    \param ev The event structure to be destroyed. */
TMPI_EXPORT
void tMPI_Event_destroy(tMPI_Event *ev);

/*! \brief Wait for an event to occur.

   Sets the number of events that had occurred during the wait in N.
   \param ev The event structure to wait on.
   \returns  The number of events that have occurred at function
             return time. */
TMPI_EXPORT
int tMPI_Event_wait(tMPI_Event *ev);

#ifdef DOXYGEN
/*! \brief Signal an event, possibly waking an tMPI_Event_wait().

    \param ev  The event to signal. */
TMPI_EXPORT
void tMPI_Event_signal(tMPI_Event *ev);
#else
#define tMPI_Event_signal(ev) \
    { \
        tMPI_Atomic_memory_barrier_rel(); \
        tMPI_Atomic_fetch_add( &((ev)->sync), 1); \
    }
#endif

#ifdef DOXYGEN
/*! \brief Signal processing of an event.

   Each event that is handled by the receiver, must be processed through
   tMPI_Event_process(). Unprocessed events will make tMPI_Event_wait() return
   immediately.

   \param ev  The event object.
   \param N   The number of processed events. */
TMPI_EXPORT
void tMPI_Event_process(tMPI_Event *ev, int N);
#else
#define tMPI_Event_process(ev, N) \
    { \
        (ev)->last_sync += N; \
    }
#endif

#endif /* TMPI_EVENT_H_ */
