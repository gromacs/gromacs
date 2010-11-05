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


/* this is for newer versions of gcc that have built-in intrinsics */

#define tMPI_Atomic_memory_barrier()  __sync_synchronize()


static inline int tMPI_Atomic_add_return(tMPI_Atomic_t *a, volatile int i)
{
    return __sync_add_and_fetch( &(a->value), i);
}

static inline int tMPI_Atomic_fetch_add(tMPI_Atomic_t *a, volatile int i)
{
    return __sync_fetch_and_add( &(a->value), i);
}


static inline int tMPI_Atomic_cas(tMPI_Atomic_t *a, int oldval, int newval)
{
    return __sync_bool_compare_and_swap( &(a->value), oldval, newval);
}


#if 0 
/* these definitions are only used if there's no assembly versions for them:
   they're inefficient because they use compare-and-swap instead of just
   swap. */
static inline int tMPI_Atomic_swap(tMPI_Atomic_t *a, int b)
{
    int oldval;
    do
    {
        oldval=a->value;
    } while(__sync_val_compare_and_swap( &(a->value), oldval, b) != oldval);

    return oldval;
}

static inline void* tMPI_Atomic_ptr_swap(tMPI_Atomic_ptr_t *a, void *b)
{
    void *oldval;
    do
    {
        oldval=a->value;
    } while(__sync_val_compare_and_swap( &(a->value), oldval, b) != oldval);

    return oldval;
}
#endif



static inline int tMPI_Atomic_ptr_cas(tMPI_Atomic_ptr_t* a, void *oldval, 
                                      void *newval)
{
#if !defined(__INTEL_COMPILER)
    return __sync_bool_compare_and_swap( &(a->value), oldval, newval);
#else
    /* the intel compilers need integer type arguments for compare_and_swap.
        on the platforms supported by icc, size_t is always the size of
        a pointer. */
    return (__sync_bool_compare_and_swap( (size_t*)&(a->value), 
                                          (size_t)oldval, 
                                          (size_t)newval) );
#endif
}



