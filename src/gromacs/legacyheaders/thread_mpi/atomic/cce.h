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


typedef struct tMPI_Atomic
{
    volatile long value;
}
tMPI_Atomic_t;

typedef struct tMPI_Atomic_ptr
{
    volatile void* value;
}
tMPI_Atomic_ptr_t;


/* these are guaranteed to be  atomic on x86 and x86_64 */
#define tMPI_Atomic_get(a)  ((int)( (a)->value) )
#define tMPI_Atomic_set(a, i)  (((a)->value) = (i))


#define tMPI_Atomic_ptr_get(a)  ((void*)((a)->value) )
#define tMPI_Atomic_ptr_set(a, i)  (((a)->value) = (void*)(i))


#include "cce_intrinsics.h"

#include "cce_spinlock.h"
