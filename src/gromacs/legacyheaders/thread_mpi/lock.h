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

#ifndef TMPI_FASTLOCK_H_
#define TMPI_FASTLOCK_H_

#include "visibility.h"
#include "wait.h"
#include "atomic.h"

/** Fast (possibly busy-wait-based) lock type
 *
 *  This lock type forms an intermediate between the spinlocks and mutexes:
 *  it is based on a busy-wait loop, but yields to the scheduler if the lock
 *  is locked.  This is therefore the preferred type of  lock for when waits
 *  are expected to be reasonably short.
 *
 *  Variables of this type should be initialized by calling
 *  tMPI_Lock_init().
 *
 * \see
 * - tMPI_Lock_init
 * - tMPI_Lock_lock
 */
typedef struct tMPI_Lock tMPI_Lock_t;
struct tMPI_Lock
{
    tMPI_Spinlock_t   lock;      /*!< The underlying spin lock */
    TMPI_YIELD_WAIT_DATA
};


/** Initialize lock
 *
 *  \param lock     Pointer to the new lock.
 */
TMPI_EXPORT
void tMPI_Lock_init(tMPI_Lock_t *lock);


/** Perform yielding, busy-waiting locking
 *
 *  This function blocks until the lock is locked.
 *
 *  \param lock  Pointer to previously created lock.
 */
TMPI_EXPORT
void tMPI_Lock_lock(tMPI_Lock_t *lock);

/** Unlock the lock
 *
 *  \param lock  Pointer to previously created lock.
 */
TMPI_EXPORT
void tMPI_Lock_unlock(tMPI_Lock_t *lock);

/** Try to lock the lock but don't block if it is locked.
 *
 *  \param lock  Pointer to previously created lock.
 */
TMPI_EXPORT
int tMPI_Lock_trylock(tMPI_Lock_t *lock);

/** Check the status of the lock without affecting its state
 *
 *  \param lock  Pointer to previously created lock.
 */
TMPI_EXPORT
int tMPI_Lock_islocked(tMPI_Lock_t *lock);



#endif /* TMPI_FASTLOCK_H_ */
