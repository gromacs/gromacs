/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

/* This module provides routines to read/write sparse or full storage
 * matrices from/to files. It is normally used for the Hessian matrix
 * in normal mode analysis.
 */

#ifndef GMX_FILEIO_MTXIO_H
#define GMX_FILEIO_MTXIO_H

#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/utility/real.h"

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
 * To determine the format you should set *full_matrix and *sparse_matrix to NULL
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
