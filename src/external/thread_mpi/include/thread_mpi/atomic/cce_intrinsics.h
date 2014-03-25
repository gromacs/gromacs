/*
   This source code file is part of thread_mpi.
   Original for gcc written by Sander Pronk, Erik Lindahl, and possibly
   others. Modified for the Cray compiler by Daniel Landau.


   Copyright (c) 2009, Sander Pronk, Erik Lindahl.
   Copyright 2014, Cray Inc.
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

#include <intrinsics.h>

#define tMPI_Atomic_memory_barrier() __builtin_ia32_mfence()

TMPI_EXPORT
static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
    return __sync_val_compare_and_swap(&(a->value), oldval, newval) == oldval;
}

TMPI_EXPORT
static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t* a, void *oldval,
                                      void *newval)
{
    return __sync_val_compare_and_swap((size_t*)&(a->value), (size_t)oldval, (size_t)newval) == (size_t)oldval;
}

TMPI_EXPORT
static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, volatile int i)
{
    return __sync_add_and_fetch( &(a->value), i);
}
#define TMPI_ATOMIC_HAVE_NATIVE_ADD_RETURN


TMPI_EXPORT
static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, volatile int i)
{
    return __sync_fetch_and_add( &(a->value), i);
}
#define TMPI_ATOMIC_HAVE_NATIVE_FETCH_ADD
