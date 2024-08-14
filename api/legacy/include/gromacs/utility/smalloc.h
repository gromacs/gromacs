/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 * \brief
 * C-style memory allocation routines for \Gromacs.
 *
 * This header provides macros snew(), srenew(), smalloc(), and sfree() for
 * C-style memory management.  Additionally, snew_aligned() and sfree_aligned() are
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
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_SMALLOC_H
#define GMX_UTILITY_SMALLOC_H

#include <cstddef>

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
void* save_malloc(const char* name, const char* file, int line, size_t size);
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
void* save_calloc(const char* name, const char* file, int line, size_t nelem, size_t elsize);
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
void* save_realloc(const char* name, const char* file, int line, void* ptr, size_t nelem, size_t elsize);
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
void save_free(const char* name, const char* file, int line, void* ptr);

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
void* save_malloc_aligned(const char* name, const char* file, int line, size_t nelem, size_t elsize, size_t alignment);
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
void* save_calloc_aligned(const char* name, const char* file, int line, size_t nelem, size_t elsize, size_t alignment);
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
void save_free_aligned(const char* name, const char* file, int line, void* ptr);

#include <type_traits>

/*! \cond internal */
/*! \name Implementation templates for C++ memory allocation macros
 *
 * These templates are used to implement the snew() etc. macros for C++, where
 * an explicit cast is needed from `void *` (the return value of the allocation
 * wrapper functions) to the type of \p ptr.
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
template<typename T>
static inline void gmx_snew_impl(const char* name, const char* file, int line, T*& ptr, size_t nelem)
{
    // TODO: Use std::is_standard_layout_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value, "snew() called on C++ type");
    // NOLINTNEXTLINE bugprone-sizeof-expression
    ptr = static_cast<T*>(save_calloc(name, file, line, nelem, sizeof(T)));
}
/** C++ helper for srenew(). */
template<typename T>
static inline void gmx_srenew_impl(const char* name, const char* file, int line, T*& ptr, size_t nelem)
{
    // TODO: Use std::is_standard_layout_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value, "srenew() called on C++ type");
    // NOLINTNEXTLINE bugprone-sizeof-expression
    ptr = static_cast<T*>(save_realloc(name, file, line, ptr, nelem, sizeof(T)));
}
/** C++ helper for smalloc(). */
template<typename T>
static inline void gmx_smalloc_impl(const char* name, const char* file, int line, T*& ptr, size_t size)
{
    // TODO: Use std::is_standard_layout_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value, "smalloc() called on C++ type");
    ptr = static_cast<T*>(save_malloc(name, file, line, size));
}
/** C++ helper for snew_aligned(). */
template<typename T>
static inline void
gmx_snew_aligned_impl(const char* name, const char* file, int line, T*& ptr, size_t nelem, size_t alignment)
{
    // TODO: Use std::is_standard_layout_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value, "snew_aligned() called on C++ type");
    ptr = static_cast<T*>(save_calloc_aligned(name, file, line, nelem, sizeof(T), alignment));
}
/** C++ helper for sfree(). */
template<typename T>
static inline void gmx_sfree_impl(const char* name, const char* file, int line, T* ptr)
{
    // TODO: Use std::is_pod_v and std::is_void_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value || std::is_void<T>::value,
                  "sfree() called on C++ type");
    save_free(name, file, line, ptr);
}
/** C++ helper for sfree_aligned(). */
template<typename T>
static inline void gmx_sfree_aligned_impl(const char* name, const char* file, int line, T* ptr)
{
    // TODO: Use std::is_pod_v and std::is_void_v when CUDA 11 is a requirement.
    static_assert(std::is_standard_layout<T>::value || std::is_void<T>::value,
                  "sfree_aligned() called on C++ type");
    save_free_aligned(name, file, line, ptr);
}
/*! \} */
/*! \endcond */

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
/*! \def sfree
 * \brief
 * Frees memory referenced by \p ptr.
 *
 * \p ptr is allowed to be NULL, in which case nothing is done.
 *
 * \hideinitializer
 */
/*! \def sfree_aligned
 * \brief
 * Frees aligned memory referenced by \p ptr.
 *
 * This must only be called with a pointer obtained through snew_aligned().
 * \p ptr is allowed to be NULL, in which case nothing is done.
 *
 * \hideinitializer
 */

/* C++ implementation */
#define snew(ptr, nelem) gmx_snew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))
#define srenew(ptr, nelem) gmx_srenew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))
#define smalloc(ptr, size) gmx_smalloc_impl(#ptr, __FILE__, __LINE__, (ptr), (size))
#define snew_aligned(ptr, nelem, alignment) \
    gmx_snew_aligned_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem), alignment)
#define sfree(ptr) gmx_sfree_impl(#ptr, __FILE__, __LINE__, (ptr))
#define sfree_aligned(ptr) gmx_sfree_aligned_impl(#ptr, __FILE__, __LINE__, (ptr))

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
constexpr float OVER_ALLOC_FAC = 1.19F;

/*! \brief
 * Turns over allocation for variable size atoms/cg/top arrays on or off,
 * default is off.
 *
 * \todo
 * This is mdrun-specific, so it might be better to put this and
 * over_alloc_dd() much higher up.
 */
void set_over_alloc_dd(bool set);

/*! \brief
 * Returns new allocation count for domain decomposition allocations.
 *
 * Returns n when domain decomposition over allocation is off.
 * Returns OVER_ALLOC_FAC*n + 100 when over allocation in on.
 * This is to avoid frequent reallocation during domain decomposition in mdrun.
 */
int over_alloc_dd(int n);

/** Over allocation for small data types: int, real etc. */
template<typename T>
constexpr T over_alloc_small(T n)
{
    return OVER_ALLOC_FAC * n + 8000;
}

/** Over allocation for large data types: complex structs */
template<typename T>
constexpr T over_alloc_large(T n)
{
    return OVER_ALLOC_FAC * n + 1000;
}

#endif
