/*
 * This file is part of the GROMACS molecular simulation package,
 * version 4.6
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
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

#ifndef TMPI_BARRIER_H_
#define TMPI_BARRIER_H_

#include "visibility.h"
#include "wait.h"

/** Fast (possibly busy-wait-based) barrier type
 *
 *  This barrier has the same functionality as the standard
 *  tMPI_Thread_barrier_t, but since it is based on spinlocks that yield 
 *  to the scheduler in case of waiting, it provides faster synchronization 
 *  at the cost of busy-waiting, while still behaving relatively nicely
 *  to other processes/threads. This is therefore the preferred type of
 *  barrier for when waits are expected to be reasonably short.
 *
 *  Variables of this type should be initialized by calling
 *  tMPI_Barrier_init() to set the number of threads
 *  that should be synchronized.
 * 
 * \see
 * - tMPI_Barrier_init
 * - tMPI_Barrier_wait
 */
typedef struct tMPI_Barrier_t tMPI_Barrier_t;
struct tMPI_Barrier_t
{
    tMPI_Atomic_t     count;     /*!< Number of threads remaining     */
    int               threshold; /*!< Total number of threads         */
    tMPI_Atomic_t     cycle;     /*!< Current cycle (alternating 0/1) */
    TMPI_YIELD_WAIT_DATA
};



/** Initialize barrier
 *
 *  \param barrier  Pointer to _spinlock_ barrier. Note that this is not
 *                  the same datatype as the full, thread based, barrier.
 *  \param count    Number of threads to synchronize. All threads
 *                  will be released after \a count calls to 
 *                  tMPI_Barrier_wait().  
 */
TMPI_EXPORT
void tMPI_Barrier_init(tMPI_Barrier_t *barrier, int count);


/** Perform yielding, busy-waiting barrier synchronization
  *
  *  This function blocks until it has been called N times,
  *  where N is the count value the barrier was initialized with.
  *  After N total calls all threads return. The barrier automatically
  *  cycles, and thus requires another N calls to unblock another time.
  *
  *  \param barrier  Pointer to previously created barrier.
  *
  *  \return The last thread returns -1, all the others 0.
  */
TMPI_EXPORT
int tMPI_Barrier_wait(tMPI_Barrier_t *barrier);


#ifdef DOXYGEN
/** Get the number of threads to synchronize for a barrier
  * 
  *  This function returns the total number of threads the barrier 
  *  synchronizes.
  *    
  *  \param barrier  Pointer to barrier.
  *
  *  \return the number of threads to synchronize
  */
TMPI_EXPORT
int tMPI_Barrier_N(tMPI_Barrier_t *barrier);
#else
#define tMPI_Barrier_N(barrier)  ((barrier)->threshold)
#endif

#endif /* TMPI_BARRIER_H_ */
