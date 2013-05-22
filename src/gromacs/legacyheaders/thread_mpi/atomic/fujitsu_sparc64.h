/*
   This source code file is part of thread_mpi.
   Written by Sander Pronk, Erik Lindahl, and possibly others.

   Copyright (c) 2013, Sander Pronk, Erik Lindahl.
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

#define tMPI_Atomic_memory_barrier() { asm ("membar   #StoreStore | #LoadStore | #LoadLoad | #StoreLoad "); }
#define tMPI_Atomic_memory_barrier_acq() { asm ("membar   #StoreStore | #StoreLoad ");  }
#define tMPI_Atomic_memory_barrier_rel() { asm ("membar   #LoadStore | #StoreStore ");  }
#define TMPI_HAVE_ACQ_REL_BARRIERS


typedef struct tMPI_Atomic
{
    volatile int value __attribute__ ((aligned(64)));
}
tMPI_Atomic_t;


typedef struct tMPI_Atomic_ptr
{
    volatile char* volatile* value __attribute__ ((aligned(64)));  /*!< Volatile, to avoid compiler aliasing */
}
tMPI_Atomic_ptr_t;


/* On sparc64, aligned 32-bit and 64-bit memory accesses are atomic */
#define tMPI_Atomic_get(a)   (int)((a)->value)
#define tMPI_Atomic_set(a, i)  (((a)->value) = (i))
#define tMPI_Atomic_ptr_get(a)   ((a)->value)
#define tMPI_Atomic_ptr_set(a, i)  (((a)->value) = (i))

#define TMPI_SPINLOCK_INITIALIZER   { 0 }

/* we just define the CAS operation. Fetch-and-add and spinlocks are
   implemented through derived.h; this follows the recommendations of the
   Sparc v9 programming specs. */

static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
    asm ("cas [%2], %1, %0"
         : "=&r" (newval)
         : "r" (oldval), "r" (&(a->value)), "0" (newval)
         : "memory");
    return newval == oldval;
}


static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t *a, void* oldval,
                                      void* newval)
{
    asm ("casx [%2], %1, %0         "
         : "=&r" (newval)
         : "r" (oldval), "r" (&(a->value)), "0" (newval)
         : "memory");
    return newval == oldval;
}
