/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

/* This module provides routines to read/write sparse or full storage
 * matrices from/to files. It is normally used for the Hessian matrix
 * in normal mode analysis.
 */

#ifndef _MTXIO_H_
#define _MTXIO_H_

#include "types/simple.h"
#include "sparsematrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Write a full or sparse matrix to a file.
 *
 * You should provide the filename, dimensions (nrow/ncol), and
 * EITHER a pointer to a full storage matrix or a sparse storage
 * matrix. If both pointers are non-NULL a fatal error will occur.
 */
void
gmx_mtxio_write(const char *             filename,
                int                      nrow,
                int                      ncol,
                real *                   full_matrix,
                gmx_sparsematrix_t *     sparse_matrix);


/* Read a matrix from file.
 *
 * This routine will autodetect the matrix format stored in the file
 * (sparse or full) and set either the full or sparse matrix arguments (ptr to ptr)
 * to a newly allocated matrix structure. Note that the full storage
 * structure is simply nrow*ncol floating-point elements. The sparse                                                  
 * matrix structure should be freed with gmx_sparsematrix_destroy() when you are done.
 *
 * To determine the format you should set full_matrix and sparse_matrix to NULL
 * before calling this routine, and check which one is non-NULL on return.
 */
void
gmx_mtxio_read (const char *            filename,
                int *                   nrow,
                int *                   ncol,
                real **                 full_matrix,
                gmx_sparsematrix_t **   sparse_matrix);

#ifdef __cplusplus
}
#endif

#endif

