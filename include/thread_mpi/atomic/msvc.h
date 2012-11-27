/*
 * This file is part of the GROMACS molecular simulation package.
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


/* Microsoft Visual C on x86, define taken from FFTW who got it from Morten Nissov */

/* we need this for all the data types. We use WIN32_LEAN_AND_MEAN to avoid 
      polluting the global namespace. */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN

#if (!defined(inline)) && (!defined(__cplusplus))
#define inline_defined_in_atomic 1
#define inline __inline
#endif

#define tMPI_Atomic_memory_barrier()


typedef struct tMPI_Atomic
{
        LONG volatile      value; /*!< Volatile, to avoid compiler aliasing */
} tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
        void* volatile      value; /*!< Volatile, to avoid compiler aliasing */
} tMPI_Atomic_ptr_t;

typedef struct tMPI_Spinlock
{
    LONG volatile      lock;      /*!< Volatile, to avoid compiler aliasing */
} tMPI_Spinlock_t;

#define TMPI_SPINLOCK_INITIALIZER   { 0 }


#define TMPI_HAVE_SWAP


#define tMPI_Atomic_get(a)  ((a)->value) 
#define tMPI_Atomic_set(a,i)  (((a)->value) = (i))


#define tMPI_Atomic_ptr_get(a)    ((a)->value) 
#define tMPI_Atomic_ptr_set(a,i)  (((a)->value) = (void*)(i))


#define tMPI_Atomic_fetch_add(a, i)  \
    InterlockedExchangeAdd((LONG volatile *)(a), (LONG) (i))

#define tMPI_Atomic_add_return(a, i)  \
    ( (i) + InterlockedExchangeAdd((LONG volatile *)(a), (LONG) (i)) )

#define tMPI_Atomic_cas(a, oldval, newval) \
    (InterlockedCompareExchange((LONG volatile *)(a), (LONG) (newval), (LONG) (oldval)) == (LONG)oldval)

#define tMPI_Atomic_ptr_cas(a, oldval, newval) \
    (InterlockedCompareExchangePointer(&((a)->value), (PVOID) (newval),  \
                                      (PVOID) (oldval)) == (PVOID)oldval)

#define tMPI_Atomic_swap(a, b) \
    InterlockedExchange((LONG volatile *)(a), (LONG) (b))

#define tMPI_Atomic_ptr_swap(a, b) \
    InterlockedExchangePointer(&((a)->value), (PVOID) (b))



static inline void tMPI_Spinlock_init(tMPI_Spinlock_t *   x)
{
    x->lock = 0;
}

# define tMPI_Spinlock_lock(x)   \
    while((InterlockedCompareExchange((LONG volatile *)(x), 1, 0))!=0)


#define tMPI_Spinlock_trylock(x)   \
    InterlockedCompareExchange((LONG volatile *)(x), 1, 0)


static inline void tMPI_Spinlock_unlock(tMPI_Spinlock_t *   x)
{
    x->lock = 0;
}


static inline int tMPI_Spinlock_islocked(const tMPI_Spinlock_t *   x)
{
    return (*(volatile signed char *)(&(x)->lock) != 0);
}


static inline void tMPI_Spinlock_wait(tMPI_Spinlock_t *   x)
{
    while(tMPI_Spinlock_islocked(x))
    {
        /*Sleep(0);*/
    }
}



