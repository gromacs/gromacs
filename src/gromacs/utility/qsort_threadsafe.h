/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2014, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Portable implementation of threadsafe quicksort.
 *
 * This module provides a \Gromacs version of the qsort() routine defined.
 * It is not highly optimized, but it is threadsafe, i.e. multiple threads
 * can simultaneously call gmx_qsort() with different data.
 *
 * The rational is that some implementations of qsort() are not threadsafe.
 * For instance qsort() in glibc contains a bug which makes it not threadsafe:
 * http://sources.redhat.com/bugzilla/show_bug.cgi?id=11655
 * On the other hand, system qsort() might be faster than our own.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_QSORT_THREADSAFE_H
#define GMX_UTILITY_QSORT_THREADSAFE_H

#include <stdlib.h>

/* For GMX_THREAD_MPI */
#include "config.h"

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

/*! \addtogroup module_utility
 * \{
 */

/*! \brief
 * Portable threadsafe sort routine.
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
gmx_qsort(void            *base,
          size_t           nmemb,
          size_t           size,
          int            (*compar)(const void *, const void *));


/*! \def gmx_qsort_threadsafe
 * \brief
 * Threadsafe qsort().
 *
 * Expands to gmx_qsort() if Gromacs is built with threading, or system qsort()
 * otherwise.
 */
#ifdef GMX_THREAD_MPI
#define gmx_qsort_threadsafe gmx_qsort
#else
#define gmx_qsort_threadsafe qsort
#endif

/*! \} */

#ifdef __cplusplus
}
#endif

#endif
