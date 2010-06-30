/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2010
 * David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifndef _GMX_SORT_H_
#define _GMX_SORT_H_

/** @file gmx_sort.h
 *
 *  @brief Portable implementation of thread-safe sort routines.
 *
 *
 *  This module provides a Gromacs version of the qsort() routine defined.
 *  It is not highly optimized, but it is thread safe, i.e. multiple threads
 *  can simultaneously call gmx_qsort with different data.
 */

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif
#if 0
} /* fixes auto-indentation problems */
#endif


/*
 *  @param base    Pointer to first element in list to sort
 *  @param nmemb   Number of elements in list
 *  @param size    Size in bytes of each element
 *  @param compar  Comparison function that takes two pointers to elements
 *                 being compared as arguments. The function should return an
 *                 integer less than, equal to, or greater than zero if the 
 *                 first argument is considered to be respectively less than,
 *                 equal to, or greater than the second.
 */
void
gmx_qsort(void *           base, 
          size_t           nmemb, 
          size_t           size, 
          int            (*compar)(const void *, const void *));


#ifdef GMX_THREADS
/* Some implementations of qsort are not threadsafe.
 * For instance qsort in glibc contains a bug which makes it non-threadsafe:
 * http://sources.redhat.com/bugzilla/show_bug.cgi?id=11655
 */
#define qsort_threadsafe gmx_qsort
#else
/* System qsort might be faster than our own */
#define qsort_threadsafe qsort
#endif


#ifdef __cplusplus
}
#endif


#endif /* _GMX_SORT_H_ */
