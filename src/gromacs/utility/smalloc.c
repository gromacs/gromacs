/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif
#ifdef HAVE__ALIGNED_MALLOC
#include <malloc.h>
#endif

#include "thread_mpi/threads.h"

#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/fatalerror.h"
#ifdef PRINT_ALLOC_KB
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxmpi.h"
#endif

static gmx_bool            g_bOverAllocDD     = FALSE;
static tMPI_Thread_mutex_t g_over_alloc_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

#ifdef DEBUG
static void log_action(int bMal, const char *what, const char *file, int line,
                       int nelem, int size, void *ptr)
{
    static int btot = 0;
    char      *NN   = "NULL";
    int        bytes;

    bytes = size*nelem;
    if (!bMal)
    {
        bytes = -bytes;
    }

    tMPI_Thread_mutex_lock(&gmx_logfile_mtx);

    /* This total memory count is not correct, since with realloc
     * it adds the whole size again, not just the increment.
     */
    /* This static variable is protected by the mutex too... */
    btot += bytes;

    bytes /= 1024;
    if (debug && (bytes != 0))
    {
        fprintf(debug, "%s:%d kB (%7d kB) [%s, line %d, nelem %d, size %d]\n",
                what ? what : NN, bytes, btot/1024,
                file ? file : NN, line, nelem, size);
    }
    /* Print to stderr for things larger than 1 MB */
    if (bytes >= 1024 || bytes <= -1024)
    {
        char *fname = NULL;
        if (file)
        {
            fname = strrchr(file, DIR_SEPARATOR);
            if (fname)
            {
                fname++;
            }
            else
            {
                fname = file;
            }
        }
        printf("%s: %.1f MB [%s, line %d, nelem %d, size %d]\n",
               what ? what  : NN, bytes/1024.0,
               file ? fname : NN, line, nelem, size);
    }
    tMPI_Thread_mutex_unlock(&gmx_logfile_mtx);
}
#endif

void *save_malloc(const char *name, const char *file, int line, size_t size)
{
    void *p;

    p = NULL;
    if (size == 0)
    {
        p = NULL;
    }
    else
    {
        if ((p = malloc(size)) == NULL)
        {
            gmx_fatal(errno, __FILE__, __LINE__,
                      "Not enough memory. Failed to malloc %"GMX_PRId64 " bytes for %s\n"
                      "(called from file %s, line %d)",
                      (gmx_int64_t)size, name, file, line);
        }
        (void) memset(p, 0, size);
    }
#ifdef DEBUG
    log_action(1, name, file, line, 1, size, p);
#endif
    return p;
}

void *save_calloc(const char *name, const char *file, int line,
                  size_t nelem, size_t elsize)
{
    void *p;

    p = NULL;
    if ((nelem == 0) || (elsize == 0))
    {
        p = NULL;
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (nelem*elsize >= PRINT_ALLOC_KB*1024)
        {
            int rank = gmx_node_rank();
            printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   nelem*elsize/1048576.0, name, file, line, rank);
        }
#endif
#ifdef GMX_BROKEN_CALLOC
        /* emulate calloc(3) with malloc/memset on machines with
           a broken calloc, e.g. in -lgmalloc on cray xt3. */
        if ((p = malloc((size_t)nelem*(size_t)elsize)) == NULL)
        {
            gmx_fatal(errno, __FILE__, __LINE__,
                      "Not enough memory. Failed to calloc %"GMX_PRId64
                      " elements of size %"GMX_PRId64
                      " for %s\n(called from file %s, line %d)",
                      (gmx_int64_t)nelem, (gmx_int64_t)elsize,
                      name, file, line);
        }
        memset(p, 0, (size_t) (nelem * elsize));
#else
        if ((p = calloc((size_t)nelem, (size_t)elsize)) == NULL)
        {
            gmx_fatal(errno, __FILE__, __LINE__,
                      "Not enough memory. Failed to calloc %"GMX_PRId64
                      " elements of size %"GMX_PRId64
                      " for %s\n(called from file %s, line %d)",
                      (gmx_int64_t)nelem, (gmx_int64_t)elsize, name, file, line);
        }
#endif
    }
#ifdef DEBUG
    log_action(1, name, file, line, nelem, elsize, p);
#endif
    return p;
}

void *save_realloc(const char *name, const char *file, int line, void *ptr,
                   size_t nelem, size_t elsize)
{
    void  *p;
    size_t size = nelem*elsize;

    p = NULL;
    if (size == 0)
    {
        save_free(name, file, line, ptr);
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (size >= PRINT_ALLOC_KB*1024)
        {
            int rank = gmx_node_rank();
            printf("Reallocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   size/1048576.0, name, file, line, rank);
        }
#endif
        if (ptr == NULL)
        {
            p = malloc((size_t)size);
        }
        else
        {
            p = realloc(ptr, (size_t)size);
        }
        if (p == NULL)
        {
            gmx_fatal(errno, __FILE__, __LINE__,
                      "Not enough memory. Failed to realloc %"GMX_PRId64 " bytes for %s, %s=%x\n"
                      "(called from file %s, line %d)",
                      (gmx_int64_t)size, name, name, ptr, file, line);
        }
#ifdef DEBUG
        log_action(1, name, file, line, 1, size, p);
#endif
    }
    return p;
}

void save_free(const char gmx_unused *name, const char gmx_unused *file, int gmx_unused line, void *ptr)
{
#ifdef DEBUG
    log_action(0, name, file, line, 0, 0, ptr);
#endif
    if (ptr != NULL)
    {
        free(ptr);
    }
}

/* If we don't have useful routines for allocating aligned memory,
 * then we have to use the old-style GROMACS approach bitwise-ANDing
 * pointers to ensure alignment. We store the pointer to the originally
 * allocated region in the space before the returned pointer */

/* we create a positive define for the absence of an system-provided memalign */
#if (!defined HAVE_POSIX_MEMALIGN && !defined HAVE_MEMALIGN && \
     !defined HAVE__ALIGNED_MALLOC)
#define GMX_OWN_MEMALIGN
#endif


/* Pointers allocated with this routine should only be freed
 * with save_free_aligned, however this will only matter
 * on systems that lack posix_memalign() and memalign() when
 * freeing memory that needed to be adjusted to achieve
 * the necessary alignment. */
void *save_malloc_aligned(const char *name, const char *file, int line,
                          size_t nelem, size_t elsize, size_t alignment)
{
    void   **aligned  = NULL;
    void    *malloced = NULL;
    gmx_bool allocate_fail;

    if (alignment == 0)
    {
        gmx_fatal(errno, __FILE__, __LINE__,
                  "Cannot allocate aligned memory with alignment of zero!\n(called from file %s, line %d)", file, line);
    }


    if (nelem == 0 || elsize == 0)
    {
        aligned  = NULL;
    }
    else
    {
#ifdef PRINT_ALLOC_KB
        if (nelem*elsize >= PRINT_ALLOC_KB*1024)
        {
            int rank = gmx_node_rank();
            printf("Allocating %.1f MB for %s (called from file %s, line %d on %d)\n",
                   nelem*elsize/1048576.0, name, file, line, rank);
        }
#endif

#ifdef HAVE_POSIX_MEMALIGN
        allocate_fail = (0 != posix_memalign(&malloced, alignment, nelem*elsize));
#elif defined HAVE_MEMALIGN
        allocate_fail = ((malloced = memalign(alignment, nelem*elsize)) == NULL);
#elif defined HAVE__ALIGNED_MALLOC
        allocate_fail = ((malloced = _aligned_malloc(nelem*elsize, alignment))
                         == NULL);
#else
        allocate_fail = ((malloced = malloc(nelem*elsize+alignment+
                                            sizeof(void*))) == NULL);
#endif
        if (allocate_fail)
        {
            gmx_fatal(errno, __FILE__, __LINE__,
                      "Not enough memory. Failed to allocate %u aligned elements of size %u for %s\n(called from file %s, line %d)", nelem, elsize, name, file, line);
        }
        /* we start with the original pointer */
        aligned = (void**)malloced;

#ifdef GMX_OWN_MEMALIGN
        /* Make the aligned pointer, and save the underlying pointer that
         * we're allowed to free(). */

        /* we first make space to store that underlying pointer: */
        aligned = aligned + 1;
        /* then we apply a bit mask */
        aligned = (void *) (((size_t) aligned + alignment - 1) &
                            (~((size_t) (alignment-1))));
        /* and we store the original pointer in the area just before the
           pointer we're going to return */
        aligned[-1] = malloced;
#endif
    }
    return (void*)aligned;
}

void *save_calloc_aligned(const char *name, const char *file, int line,
                          size_t nelem, size_t elsize, size_t alignment)
{
    void *aligned = save_malloc_aligned(name, file, line, nelem, elsize, alignment);
    if (aligned != NULL)
    {
        memset(aligned, 0, (size_t)(nelem * elsize));
    }
    return aligned;
}

/* This routine can NOT be called with any pointer */
void save_free_aligned(const char gmx_unused *name, const char gmx_unused *file, int gmx_unused line, void *ptr)
{
    int   i, j;
    void *free = ptr;

    if (NULL != ptr)
    {
#ifdef GMX_OWN_MEMALIGN
        /* we get the pointer from just before the memaligned pointer */
        free = ((void**)ptr)[-1];
#endif

#ifndef HAVE__ALIGNED_MALLOC
        /* (Now) we're allowed to use a normal free() on this pointer. */
        save_free(name, file, line, free);
#else
        _aligned_free(free);
#endif
    }
}

void set_over_alloc_dd(gmx_bool set)
{
    tMPI_Thread_mutex_lock(&g_over_alloc_mutex);
    /* we just make sure that we don't set this at the same time.
       We don't worry too much about reading this rarely-set variable */
    g_bOverAllocDD = set;
    tMPI_Thread_mutex_unlock(&g_over_alloc_mutex);
}

int over_alloc_dd(int n)
{
    if (g_bOverAllocDD)
    {
        return OVER_ALLOC_FAC*n + 100;
    }
    else
    {
        return n;
    }
}
