/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020,2021, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "smalloc.h"

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>

#include <mutex>

#include <cstring>

#include "gromacs/utility/fatalerror.h"
#ifdef PRINT_ALLOC_KB
#    include "gromacs/utility/basenetwork.h"
#    include "gromacs/utility/gmxmpi.h"
#endif

static bool       g_bOverAllocDD = false;
static std::mutex g_overAllocMutex;

void* save_malloc(const char* name, const char* file, int line, size_t size)
{
    void* p = nullptr;

    if (size == 0)
    {
        p = nullptr;
    }
    else
    {
        if ((p = malloc(size)) == nullptr)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory. Failed to malloc %" PRId64
                      " bytes for %s\n"
                      "(called from file %s, line %d)",
                      static_cast<int64_t>(size),
                      name,
                      file,
                      line);
        }
        (void)memset(p, 0, size);
    }
    return p;
}

void* save_calloc(const char* name, const char* file, int line, size_t nelem, size_t elsize)
{
    void* p = nullptr;

    if ((nelem == 0) || (elsize == 0))
    {
        p = nullptr;
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (nelem * elsize >= PRINT_ALLOC_KB * 1024)
        {
            int rank = gmx_node_rank();
            printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   nelem * elsize / 1048576.0,
                   name,
                   file,
                   line,
                   rank);
        }
#endif
#if GMX_BROKEN_CALLOC
        /* emulate calloc(3) with malloc/memset on machines with
           a broken calloc, e.g. in -lgmalloc on cray xt3. */
        if ((p = malloc((size_t)nelem * (size_t)elsize)) == NULL)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory. Failed to calloc %" PRId64 " elements of size %" PRId64
                      " for %s\n(called from file %s, line %d)",
                      (int64_t)nelem,
                      (int64_t)elsize,
                      name,
                      file,
                      line);
        }
        memset(p, 0, (size_t)(nelem * elsize));
#else
        if ((p = calloc(nelem, elsize)) == nullptr)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory. Failed to calloc %" PRId64 " elements of size %" PRId64
                      " for %s\n(called from file %s, line %d)",
                      static_cast<int64_t>(nelem),
                      static_cast<int64_t>(elsize),
                      name,
                      file,
                      line);
        }
#endif
    }
    return p;
}

void* save_realloc(const char* name, const char* file, int line, void* ptr, size_t nelem, size_t elsize)
{
    void*  p    = nullptr;
    size_t size = nelem * elsize;

    if (size == 0)
    {
        save_free(name, file, line, ptr);
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (size >= PRINT_ALLOC_KB * 1024)
        {
            int rank = gmx_node_rank();
            printf("Reallocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   size / 1048576.0,
                   name,
                   file,
                   line,
                   rank);
        }
#endif
        if (ptr == nullptr)
        {
            p = malloc(size);
        }
        else
        {
            p = realloc(ptr, size);
        }
        if (p == nullptr)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory. Failed to realloc %zu bytes for %s, %s=%p\n"
                      "(called from file %s, line %d)",
                      size,
                      name,
                      name,
                      ptr,
                      file,
                      line);
        }
    }
    return p;
}

void save_free(const char gmx_unused* name, const char gmx_unused* file, int gmx_unused line, void* ptr)
{
    if (ptr != nullptr)
    {
        free(ptr);
    }
}

/* Pointers allocated with this routine should only be freed
 * with save_free_aligned. */
void* save_malloc_aligned(const char* name, const char* file, int line, size_t nelem, size_t elsize, size_t alignment)
{
    void* p = nullptr;

    if (alignment == 0)
    {
        gmx_fatal(errno,
                  __FILE__,
                  __LINE__,
                  "Cannot allocate aligned memory with alignment of zero!\n(called from file %s, "
                  "line %d)",
                  file,
                  line);
    }

    if (nelem == 0 || elsize == 0)
    {
        p = nullptr;
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (nelem * elsize >= PRINT_ALLOC_KB * 1024)
        {
            int rank = gmx_node_rank();
            printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   nelem * elsize / 1048576.0,
                   name,
                   file,
                   line,
                   rank);
        }
#endif

        // Adhere to the implementation requirements. Also avoids false
        // sharing.
        size_t bytes                = nelem * elsize;
        auto   multiplesOfAlignment = (bytes / alignment + 1) * alignment;
        p                           = std::aligned_alloc(alignment, multiplesOfAlignment);

        if (p == nullptr)
        {
            gmx_fatal(errno,
                      __FILE__,
                      __LINE__,
                      "Not enough memory or invalid alignment. Failed to allocate %zu aligned "
                      "elements of size %zu for %s\n(called from file %s, line %d)",
                      nelem,
                      elsize,
                      name,
                      file,
                      line);
        }
    }
    return p;
}

void* save_calloc_aligned(const char* name, const char* file, int line, size_t nelem, size_t elsize, size_t alignment)
{
    void* aligned = save_malloc_aligned(name, file, line, nelem, elsize, alignment);
    if (aligned != nullptr)
    {
        memset(aligned, 0, static_cast<size_t>(nelem * elsize));
    }
    return aligned;
}

/* This routine should only be called with pointers obtained from save_*_aligned family of functions */
void save_free_aligned(const char gmx_unused* name, const char gmx_unused* file, int gmx_unused line, void* ptr)
{
    std::free(ptr);
}

void set_over_alloc_dd(bool set)
{
    std::lock_guard<std::mutex> lock(g_overAllocMutex);
    /* we just make sure that we don't set this at the same time.
       We don't worry too much about reading this rarely-set variable */
    g_bOverAllocDD = set;
}

int over_alloc_dd(int n)
{
    if (g_bOverAllocDD)
    {
        return static_cast<int>(OVER_ALLOC_FAC * n + 100);
    }
    else
    {
        return n;
    }
}
