/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
/*! \file
 * \brief
 * C memory allocation routines for \Gromacs.
 *
 * This header provides macros snew(), srenew(), smalloc(), and sfree() for
 * C memory management.  Additionally, snew_aligned() and sfree_aligned() are
 * provided for managing memory with a specified byte alignment.
 *
 * If an allocation fails, the program is halted by calling gmx_fatal(), which
 * outputs source file and line number and the name of the variable involved.
 * This frees calling code from the trouble of checking the result of the
 * allocations everywhere.  It also provides a location for centrally logging
 * memory allocations for diagnosing memory usage (currently can only enabled
 * by changing the source code).  Additionally, sfree() works also with a
 * `NULL` parameter, which standard free() does not.
 *
 * The macros forward the calls to functions save_malloc(), save_calloc(),
 * save_realloc(), save_free(), save_calloc_aligned(), and save_free_aligned().
 * There are a few low-level locations in \Gromacs that call these directly,
 * but generally the macros should be used.
 * save_malloc_aligned() exists for this purpose, although there is no macro to
 * invoke it.
 *
 * \if internal
 * As an implementation detail, the macros need a different internal
 * implementation for C and C++ code.  This is because C accepts conversions
 * from `void *` to any pointer type, but C++ doesn't.  And in order to cast
 * the returned pointer to a correct type, a C++ template needs to be used to
 * get access to the type.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_SMALLOC_H
#define GMX_UTILITY_SMALLOC_H

#include <stddef.h>

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * \Gromacs wrapper for malloc().
 *
 * \param[in] name   Variable name identifying the allocation.
 * \param[in] file   Source code file where the allocation originates from.
 * \param[in] line   Source code line where the allocation originates from.
 * \param[in] size   Number of bytes to allocate.
 * \returns   Pointer to the allocated space.
 *
 * This should generally be called through smalloc(), not directly.
 */
void *save_malloc(const char *name, const char *file, int line, size_t size);
/*! \brief
 * \Gromacs wrapper for calloc().
 *
 * \param[in] name   Variable name identifying the allocation.
 * \param[in] file   Source code file where the allocation originates from.
 * \param[in] line   Source code line where the allocation originates from.
 * \param[in] nelem  Number of elements to allocate.
 * \param[in] elsize Number of bytes per element.
 * \returns   Pointer to the allocated space.
 *
 * This should generally be called through snew(), not directly.
 */
void *save_calloc(const char *name, const char *file, int line,
                  size_t nelem, size_t elsize);
/*! \brief
 * \Gromacs wrapper for realloc().
 *
 * \param[in] name   Variable name identifying the allocation.
 * \param[in] file   Source code file where the allocation originates from.
 * \param[in] line   Source code line where the allocation originates from.
 * \param[in] ptr    Pointer to the previously allocated memory (can be NULL).
 * \param[in] nelem  Number of elements to allocate.
 * \param[in] elsize Number of bytes per element.
 * \returns   Pointer to the allocated space.
 *
 * As with realloc(), if \p ptr is NULL, memory is allocated as if malloc() was
 * called.
 * This should generally be called through srenew(), not directly.
 *
 * Note that the allocated memory is not initialized to zero.
 */
void *save_realloc(const char *name, const char *file, int line,
                   void *ptr, size_t nelem, size_t elsize);
/*! \brief
 * \Gromacs wrapper for free().
 *
 * \param[in] name   Variable name identifying the deallocation.
 * \param[in] file   Source code file where the deallocation originates from.
 * \param[in] line   Source code line where the deallocation originates from.
 * \param[in] ptr    Pointer to the allocated memory (can be NULL).
 *
 * If \p ptr is NULL, does nothing.
 * This should generally be called through sfree(), not directly.
 * This never fails.
 */
void save_free(const char *name, const char *file, int line, void *ptr);

/*! \brief
 * \Gromacs wrapper for allocating aligned memory.
 *
 * \param[in] name   Variable name identifying the allocation.
 * \param[in] file   Source code file where the allocation originates from.
 * \param[in] line   Source code line where the allocation originates from.
 * \param[in] nelem  Number of elements to allocate.
 * \param[in] elsize Number of bytes per element.
 * \param[in] alignment Requested alignment in bytes.
 * \returns   Pointer to the allocated space, aligned at `alignment`-byte
 *     boundary.
 *
 * There is no macro that invokes this function.
 *
 * The returned pointer should only be freed with a call to save_free_aligned().
 */
void *save_malloc_aligned(const char *name, const char *file, int line,
                          size_t nelem, size_t elsize, size_t alignment);
/*! \brief
 * \Gromacs wrapper for allocating zero-initialized aligned memory.
 *
 * \param[in] name   Variable name identifying the allocation.
 * \param[in] file   Source code file where the allocation originates from.
 * \param[in] line   Source code line where the allocation originates from.
 * \param[in] nelem  Number of elements to allocate.
 * \param[in] elsize Number of bytes per element.
 * \param[in] alignment Requested alignment in bytes.
 * \returns   Pointer to the allocated space, aligned at `alignment`-byte
 *     boundary.
 *
 * This should generally be called through snew_aligned(), not directly.
 *
 * The returned pointer should only be freed with a call to save_free_aligned().
 */
void *save_calloc_aligned(const char *name, const char *file, int line,
                          size_t nelem, size_t elsize, size_t alignment);
/*! \brief
 * \Gromacs wrapper for freeing aligned memory.
 *
 * \param[in] name   Variable name identifying the deallocation.
 * \param[in] file   Source code file where the deallocation originates from.
 * \param[in] line   Source code line where the deallocation originates from.
 * \param[in] ptr    Pointer to the allocated memory (can be NULL).
 *
 * If \p ptr is NULL, does nothing.
 * \p ptr should have been allocated with save_malloc_aligned() or
 * save_calloc_aligned().
 * This should generally be called through sfree_aligned(), not directly.
 * This never fails.
 */
void save_free_aligned(const char *name, const char *file, int line, void *ptr);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/*! \cond internal */
/*! \name Implementation templates for C++ memory allocation macros
 *
 * These templates are used to implement the snew() etc. macros for C++, where
 * an explicit cast is needed from `void *` (the return value of the allocation
 * wrapper functions) to the thpe of \p ptr.
 *
 * Having these as `static` avoid some obscure bugs if several files define
 * distinct data structures with identical names and allocate memory for them
 * using snew().  By the C++ standard, such declarations cause undefined
 * behavior, but can be difficult to spot in the existing C code.
 * Without the `static` (and if the compiler does not inline the calls), the
 * linker cannot that data structures with identical names are actually
 * different and links calls to these template functions incorrectly, which can
 * result in allocation of an incorrect amount of memory if the element size is
 * computed within the function.
 *
 * The size cannot be passed as a parameter to the function either, since that
 * provokes warnings from cppcheck for some invocations, where a complex
 * expression is passed as \p ptr.
 */
/*! \{ */
/** C++ helper for snew(). */
template <typename T> static inline
void gmx_snew_impl(const char *name, const char *file, int line,
                   T * &ptr, size_t nelem)
{
    ptr = (T *)save_calloc(name, file, line, nelem, sizeof(T));
}
/** C++ helper for srenew(). */
template <typename T> static inline
void gmx_srenew_impl(const char *name, const char *file, int line,
                     T * &ptr, size_t nelem)
{
    ptr = (T *)save_realloc(name, file, line, ptr, nelem, sizeof(T));
}
/** C++ helper for smalloc(). */
template <typename T> static inline
void gmx_smalloc_impl(const char *name, const char *file, int line,
                      T * &ptr, size_t size)
{
    ptr = (T *)save_malloc(name, file, line, size);
}
/** C++ helper for snew_aligned(). */
template <typename T> static inline
void gmx_snew_aligned_impl(const char *name, const char *file, int line,
                           T * &ptr, size_t nelem, size_t alignment)
{
    ptr = (T *)save_calloc_aligned(name, file, line, nelem, sizeof(T), alignment);
}
/*! \] */
/*! \endcond */
#endif /* __cplusplus */

/*! \def snew
 * \brief
 * Allocates memory for a given number of elements.
 *
 * \param[out] ptr   Pointer to allocate.
 * \param[in]  nelem Number of elements to allocate.
 *
 * Allocates memory for \p nelem elements of type \p *ptr and sets this to
 * \p ptr.  The allocated memory is initialized to zeros.
 *
 * \hideinitializer
 */
/*! \def srenew
 * \brief
 * Reallocates memory for a given number of elements.
 *
 * \param[in,out] ptr   Pointer to allocate/reallocate.
 * \param[in]     nelem Number of elements to allocate.
 *
 * (Re)allocates memory for \p ptr such that it can hold \p nelem elements of
 * type \p *ptr, and sets the new pointer to \p ptr.
 * If \p ptr is `NULL`, memory is allocated as if it was new.
 * If \p nelem is zero, \p ptr is freed (if not `NULL`).
 * Note that the allocated memory is not initialized, unlike with snew().
 *
 * \hideinitializer
 */
/*! \def smalloc
 * \brief
 * Allocates memory for a given number of bytes.
 *
 * \param[out] ptr  Pointer to allocate.
 * \param[in]  size Number of bytes to allocate.
 *
 * Allocates memory for \p size bytes and sets this to \p ptr.
 * The allocated memory is initialized to zero.
 *
 * \hideinitializer
 */
/*! \def snew_aligned
 * \brief
 * Allocates aligned memory for a given number of elements.
 *
 * \param[out] ptr       Pointer to allocate.
 * \param[in]  nelem     Number of elements to allocate.
 * \param[in]  alignment Requested alignment in bytes.
 *
 * Allocates memory for \p nelem elements of type \p *ptr and sets this to
 * \p ptr.  The returned pointer is `alignment`-byte aligned.
 * The allocated memory is initialized to zeros.
 *
 * The returned pointer should only be freed with sfree_aligned().
 *
 * \hideinitializer
 */
#ifdef __cplusplus

/* C++ implementation */
#define snew(ptr, nelem) \
    gmx_snew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))
#define srenew(ptr, nelem) \
    gmx_srenew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))
#define smalloc(ptr, size) \
    gmx_smalloc_impl(#ptr, __FILE__, __LINE__, (ptr), (size))
#define snew_aligned(ptr, nelem, alignment) \
    gmx_snew_aligned_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem), alignment)

#else

/* C implementation */
#define snew(ptr, nelem) \
    (ptr) = save_calloc(#ptr, __FILE__, __LINE__, (nelem), sizeof(*(ptr)))
#define srenew(ptr, nelem) \
    (ptr) = save_realloc(#ptr, __FILE__, __LINE__, (ptr), (nelem), sizeof(*(ptr)))
#define smalloc(ptr, size) \
    (ptr) = save_malloc(#ptr, __FILE__, __LINE__, size)
#define snew_aligned(ptr, nelem, alignment) \
    (ptr) = save_calloc_aligned(#ptr, __FILE__, __LINE__, (nelem), sizeof(*(ptr)), alignment)

#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief
 * Frees memory referenced by \p ptr.
 *
 * \p ptr is allowed to be NULL, in which case nothing is done.
 *
 * \hideinitializer
 */
#define sfree(ptr) save_free(#ptr, __FILE__, __LINE__, (ptr))

/*! \brief
 * Frees aligned memory referenced by \p ptr.
 *
 * This must only be called with a pointer obtained through snew_aligned().
 * \p ptr is allowed to be NULL, in which case nothing is done.
 *
 * \hideinitializer
 */
#define sfree_aligned(ptr) save_free_aligned(#ptr, __FILE__, __LINE__, (ptr))

/*! \brief
 * Over allocation factor for memory allocations.
 *
 * Memory (re)allocation can be VERY slow, especially with some
 * MPI libraries that replace the standard malloc and realloc calls.
 * To avoid slow memory allocation we use over_alloc to set the memory
 * allocation size for large data blocks. Since this scales the size
 * with a factor, we use log(n) realloc calls instead of n.
 * This can reduce allocation times from minutes to seconds.
 *
 * This factor leads to 4 realloc calls to double the array size.
 */
#define OVER_ALLOC_FAC 1.19

/*! \brief
 * Turns over allocation for variable size atoms/cg/top arrays on or off,
 * default is off.
 *
 * \todo
 * This is mdrun-specific, so it might be better to put this and
 * over_alloc_dd() much higher up.
 */
void set_over_alloc_dd(gmx_bool set);

/*! \brief
 * Returns new allocation count for domain decomposition allocations.
 *
 * Returns n when domain decomposition over allocation is off.
 * Returns OVER_ALLOC_FAC*n + 100 when over allocation in on.
 * This is to avoid frequent reallocation during domain decomposition in mdrun.
 */
int over_alloc_dd(int n);

/** Over allocation for small data types: int, real etc. */
#define over_alloc_small(n) (int)(OVER_ALLOC_FAC*(n) + 8000)

/** Over allocation for large data types: complex structs */
#define over_alloc_large(n) (int)(OVER_ALLOC_FAC*(n) + 1000)

#ifdef __cplusplus
}
#endif

#endif
