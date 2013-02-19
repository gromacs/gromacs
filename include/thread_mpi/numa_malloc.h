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


#ifndef TMPI_NUMA_MALLOC_H_
#define TMPI_NUMA_MALLOC_H_

/*! \file

    \brief NUMA aware memory allocators.

    This file provides a NUMA aware version of malloc() and friends to
    force local allocation, preventing expensive 'remote' memory access.

    Note that memory allocated with tMPI_Malloc_local(), etc. MUST be
    freed with tMPI_Free_numa().

    Currently this is only implemented on Windows. Check for the presence
    of these functions with
    \code
   #ifdef TMPI_NUMA_MALLOC
    ....
   #endif
    \endcode
 */


#include "visibility.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* Avoids screwing up auto-indentation */
#endif


#if (defined(WIN32) || defined( _WIN32 ) || defined(WIN64) || defined( _WIN64 )) && !defined (__CYGWIN__) && !defined (__CYGWIN32__)

#define TMPI_NUMA_MALLOC

#endif


/*! \brief Allocate local memory by size in bytes.

    \see malloc()

    \param[in] size  = size in units of sizeof() of the memory to be allocated

    \return Pointer to allocated memory, or NULL in case of failure
 */
TMPI_EXPORT
void *tMPI_Malloc_local(size_t size);


/*! \brief Allocate local memory by array and element size.

    \see calloc()

    \param[in] nmemb  = number of array members
    \param[in] size  = size in units of sizeof() of a single array element

    \return Pointer to allocated memory, or NULL in case of failure
 */
TMPI_EXPORT
void *tMPI_Calloc_local(size_t nmemb, size_t size);


/*! \brief Re-allocate to local memory.

    \see realloc()

    \param[in] ptr   = input pointer of originally allocated memory (can
                        be allocated with NUMA or non-NUMA malloc())
    \param[in] size  = new size in units of sizeof().

    \return Pointer to allocated memory, or NULL in case of failure
 */
TMPI_EXPORT
void *tMPI_Realloc_local(void *ptr, size_t size);


/*! \brief Free memory allocate with any of the NUMA-aware allocators

    \see calloc()

    \param[in] ptr = pointer to block of memory allocated with a NUMA allocator

    \return Returns debug info: 0 if the memory was allocated with a NUMA
            allocator, 1 if it was allocated by another allocator.
 */
TMPI_EXPORT
int tMPI_Free_numa(void *ptr);


#ifdef __cplusplus
}
#endif

#endif /* TMPI_NUMA_MALLOC_H_ */
