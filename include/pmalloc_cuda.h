/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef PMALLOC_CUDA_H
#define PMALLOC_CUDA_H

#ifdef GMX_GPU
#define FUNC_TERM ;
#else
#define FUNC_TERM {}
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! Allocates nbytes of page-locked memory. */
void pmalloc(void **h_ptr, size_t nbytes) FUNC_TERM

/*! Allocates nbytes of page-locked memory with write-combining. */
void pmalloc_wc(void **h_ptr, size_t nbytes) FUNC_TERM

/*! Frees page locked memory allocated with pmalloc. */
void pfree(void *h_ptr) FUNC_TERM

#ifdef __cplusplus
}
#endif
#endif  /* PMALLOC_CUDA_H */


