/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

#ifndef _smalloc_h
#define _smalloc_h
#include "visibility.h"
#include <stdlib.h>

/*
 * Memory allocation routines in gromacs:
 *
 * If an allocation fails, the program is halted by means of the
 * fatal_error routine, which outputs source file and line number
 * and the name of the variable involved.
 *
 * Macro's which can be used:
 *
 * snew(ptr,nelem)
 *    Allocates memory for nelem elements and returns this in ptr.
 *    The allocated memory is initialized to zeros.
 *
 * srenew(ptr,nelem)
 *    Reallocates memory for nelem elements and returns this in ptr.
 *
 * smalloc(ptr,size)
 *    Allocates memory for size bytes and returns this in ptr.
 *
 * scalloc(ptr,nelem,elsize)
 *    Allocates memory for nelem elements of size elsize and returns
 *    this in ptr.
 *
 * srealloc(ptr,size)
 *    Reallocates memory for size bytes and returns this in ptr.
 *
 * sfree(ptr)
 *    Frees memory referenced by ptr.
 *
 * snew_aligned(ptr,nelem,alignment)
 *    Allocates memory for nelem elements and returns this in ptr.
 *    The allocated memory is initialized to zeroes.
 *    alignment=n will constrain ptr to be n-byte aligned.
 *    This pointer should only be freed with sfree_aligned, since
 *    it may not be the value returned by the underlying malloc.
 *
 * sfree_aligned(ptr)
 *    Frees aligned memory referenced by ptr.
 *
 ****************************************************************************
 *
 * Functions which are used by the macro's:
 *
 * extern void *save_malloc(char *name,char *file,int line,int size);
 *    Like alloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_calloc(char *name,char *file,int line,
 *                          size_t nelem,size_t elsize);
 *    Like calloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *
 * extern void *save_realloc(char *name,char *file,int line,
 *                           void *ptr,size_t size);
 *    Like realloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *    If ptr equals NULL, malloc is called in stead of realloc, in this way
 *    it is possible to combine first and later allocations.
 *
 * extern void save_free(char *name,char *file,int line, void *ptr);
 *    Like free, uses name, file and line to generate an error message when
 *    the free failed.
 *
 * extern size_t maxavail();
 *    Returns the maximum available allocation unit, by applying a binary
 *    search on the largest block of memory available. After allocation
 *    it invokes free to restore the original state. So it is important
 *    that free can undo the effect of a malloc.
 *
 * extern size_t memavail();
 *    Returns the total of available allocation unit, by applying maxavail
 *    until no space is left, it then frees all allocated space and returns
 *    the sum of the previously allocated space. As mentioned with maxavail,
 *    it is important that free can undo the effect of a malloc.
 *
 * extern void *save_malloc_aligned(char *name,char *file,int line,size_t size,size_t alignment);
 *    Like alloc, returns a pointer to the allocated space, uses name, file
 *    and line to generate an error message when allocation failed.
 *    The returned pointer will be n-byte aligned, where n=alignment.
 *    The pointer should only be freed with a call to save_free.
 *
 * extern void save_free_aligned(char *name,char *file,int line, void *ptr);
 *    Like free, uses name, file and line to generate an error message when
 *    the free failed. This function is intended to be called for
 *    pointers allocated with save_malloc_aligned, and may not work
 *    on normal pointers.
 */

#ifdef __cplusplus
extern "C" {
#endif

GMX_LIBGMX_EXPORT
void *save_malloc(const char *name, const char *file, int line, size_t size);
GMX_LIBGMX_EXPORT
void *save_calloc(const char *name, const char *file, int line,
                  size_t nelem, size_t elsize);
GMX_LIBGMX_EXPORT
void *save_realloc(const char *name, const char *file, int line,
                   void *ptr, size_t nelem, size_t elsize);
GMX_LIBGMX_EXPORT
void save_free(const char *name, const char *file, int line, void *ptr);
size_t maxavail(void);
size_t memavail(void);

/* Aligned-memory counterparts */

GMX_LIBGMX_EXPORT
void *save_malloc_aligned(const char *name, const char *file, int line,
                          unsigned nelem, size_t elsize, size_t alignment);
GMX_LIBGMX_EXPORT
void *save_calloc_aligned(const char *name, const char *file, int line,
                          unsigned nelem, size_t elsize, size_t alignment);
GMX_LIBGMX_EXPORT
void save_free_aligned(const char *name, const char *file, int line, void *ptr);

#ifdef __cplusplus
}

/* Use of sizeof(T) in _snew() and _srenew() can cause obscure bugs if
 * several files define distinct data structures with identical names and
 * allocate memory for them using the macros below.
 * For this reason, the size of an element is passed as a parameter.
 *
 * The C versions work fine in such cases, but when compiled with a C++
 * compiler (and if the compiler does not inline the calls), the linker cannot
 * tell that data structures with identical names are actually different and
 * links calls to these template functions incorrectly, which can result in
 * allocation of an incorrect amount of memory if the element size is computed
 * within the function. Even with the size passed as a parameter, incorrect
 * linkage will occur, but as the type is now only present in the cast, it
 * should not cause problems.
 */
template <typename T>
void _snew(const char *name, const char *file, int line,
           T * &ptr, size_t nelem, size_t elsize)
{
    ptr = (T *)save_calloc(name, file, line, nelem, elsize);
}
template <typename T>
void _srenew(const char *name, const char *file, int line,
             T * &ptr, size_t nelem, size_t elsize)
{
    ptr = (T *)save_realloc(name, file, line, ptr, nelem, elsize);
}
template <typename T>
void _smalloc(const char *name, const char *file, int line, T * &ptr, size_t size)
{
    ptr = (T *)save_malloc(name, file, line, size);
}
template <typename T>
void _srealloc(const char *name, const char *file, int line, T * &ptr, size_t size)
{
    ptr = (T *)save_realloc(name, file, line, ptr, size, sizeof(char));
}
template <typename T>
void _snew_aligned(const char *name, const char *file, int line,
                   T * &ptr, size_t nelem, size_t elsize, size_t alignment)
{
    ptr = (T *)save_calloc_aligned(name, file, line, nelem, elsize, alignment);
}

#define snew(ptr, nelem) _snew(#ptr, __FILE__, __LINE__, (ptr), (nelem), sizeof(*(ptr)))
#define srenew(ptr, nelem) _srenew(#ptr, __FILE__, __LINE__, (ptr), (nelem), sizeof(*(ptr)))
#define smalloc(ptr, size) _smalloc(#ptr, __FILE__, __LINE__, (ptr), (size))
#define srealloc(ptr, size) _srealloc(#ptr, __FILE__, __LINE__, (ptr), (size))
#define snew_aligned(ptr, nelem, alignment) _snew_aligned(#ptr, __FILE__, __LINE__, (ptr), (nelem), sizeof(*(ptr)), alignment)

#else

/* These macros work in C, not in C++ */
#define snew(ptr, nelem) (ptr) = save_calloc(#ptr, __FILE__, __LINE__, \
                                             (nelem), sizeof(*(ptr)))
#define srenew(ptr, nelem) (ptr) = save_realloc(#ptr, __FILE__, __LINE__, \
                                                (ptr), (nelem), sizeof(*(ptr)))
#define smalloc(ptr, size) (ptr) = save_malloc(#ptr, __FILE__, __LINE__, size)
#define scalloc(ptr, nelem, elsize) \
    (ptr) = save_calloc(#ptr, __FILE__, __LINE__, nelem, elsize)
#define srealloc(ptr, size) (ptr) = save_realloc(#ptr, __FILE__, __LINE__, \
                                                 (ptr), size, 1)
#define snew_aligned(ptr, nelem, alignment) (ptr) = save_calloc_aligned(#ptr, __FILE__, __LINE__, (nelem), sizeof(*(ptr)), alignment)
#endif

#define sfree(ptr) save_free(#ptr, __FILE__, __LINE__, (ptr))

/* call this ONLY with a pointer obtained through snew_aligned or
   smalloc_aligned: */
#define sfree_aligned(ptr) save_free_aligned(#ptr, __FILE__, __LINE__, (ptr))

#endif  /* _smalloc_h */
