/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2010 David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 * Copyright (c) 2012, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Portable implementation of thread-safe sort routines.
 *
 * This module provides a Gromacs version of the qsort() routine defined.
 * It is not highly optimized, but it is thread safe, i.e. multiple threads
 * can simultaneously call gmx_qsort() with different data.
 *
 * The rational is that some implementations of qsort() are not threadsafe.
 * For instance qsort in glibc contains a bug which makes it non-threadsafe:
 * http://sources.redhat.com/bugzilla/show_bug.cgi?id=11655
 * On the other hand, system qsort might be faster than our own.
 */
#ifndef GMX_SORT_H
#define GMX_SORT_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/*! \brief
 * Portable thread-safe sort routine.
 *
 * \param base    Pointer to first element in list to sort
 * \param nmemb   Number of elements in list
 * \param size    Size in bytes of each element
 * \param compar  Comparison function that takes two pointers to elements
 *                being compared as arguments.  The function should return an
 *                integer less than, equal to, or greater than zero if the
 *                first argument is considered to be respectively less than,
 *                equal to, or greater than the second.
 */
void
gmx_qsort(void *           base,
          size_t           nmemb,
          size_t           size,
          int            (*compar)(const void *, const void *));


/*! \def qsort_threadsafe
 * \brief
 * Thread-safe qsort.
 *
 * Expands to gmx_qsort() if Gromacs is built with threading, or system qsort()
 * otherwise.
 */
#ifdef GMX_THREAD_MPI
#define qsort_threadsafe gmx_qsort
#else
#define qsort_threadsafe qsort
#endif

#ifdef __cplusplus
}
#endif

#endif /* GMX_SORT_H */
