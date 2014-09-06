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
#ifndef GMX_LINEARALGEBRA_SPARSEMATRIX_H
#define GMX_LINEARALGEBRA_SPARSEMATRIX_H

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
    gmx_sparsematrix_entry
{
    int      col;
    real     value;
} gmx_sparsematrix_entry_t;

/*! \brief Sparse matrix storage format
 *
 *  This structure specifies a storage format for a sparse matrix.
 *  The memory requirements are only proportional to the number
 *  of nonzero elements, and it provides a reasonably fast way to
 *  perform matrix-vector multiplications.
 *
 *  The data format is very similar to a neighborlist. It is optimized
 *  for fast access, but it is difficult to add entries. If you are
 *  constructing a matrix you should either do it in exactly the order
 *  specified here, or use some other more flexible intermediate structure.
 *
 *  The index array is of size nrow+1. All non-zero matrix elements
 *  on row i are stored in positions index[i] through index[i+1]-1 in
 *  the arrays column and value. The column array contains the column
 *  index for each entry, in ascending order, and the corresponding
 *  position in the value array contains the floating point matrix element.
 *
 *  index[nrow] should be equal to the total number of elements stored.
 *
 *  Thus, to find the value of matrix element [5,4] you should loop
 *  over positions index[5] to index[6]-1 in column until you either find
 *  the value 4, or a higher value (meaning the element was zero).
 *
 *  It is fairly easy to construct the matrix on-the-fly if you can do
 *  it row-by-row.
 *
 *  IMPORTANT:
 *  If compressed_symmetric is set to TRUE, you should only store EITHER the upper OR
 *  lower triangle (and the diagonal), and the other half is assumed to be
 *  symmetric. Otherwise, if compressed_symmetric==FALSE, no symmetry is implied and all
 *  elements should be stored.
 *
 *  The symmetry compression saves us a factor 2 both in storage and
 *  matrix multiplication CPU-time, which can be very useful for huge eigenproblems.
 *
 *  If you are unsure, just set compressed_symmetric to FALSE and list all elements. If
 *  you enable it but still list all elements (both upper and lower triangle) you will be sorry...
 *
 *  Internally, the sparse data is stored as a separate list for each row, where the list
 *  element is a structure with a column and (floating-point) data value. This makes it
 *  possible, although not completely transparent, to update values in random access order.
 *  The drawback is that the structure will allocate nrow memory regions.
 *  The matrix data could be stored in a single contiguous array with indices for each row,
 *  but then we could only insert elements at the end without copying the entire matrix.
 *
 *  After you have
 *
 *  In other words: Not perfect, but it works.
 */
typedef struct
    gmx_sparsematrix
{
    gmx_bool                     compressed_symmetric; /**< Store half elements and assume symmetry. */
    int                          nrow;                 /**< Number of rows in matrix                 */
    int *                        ndata;                /**< Number of entries on each row (list)     */
    int *                        nalloc;               /**< Allocated entry list length for each row */
    gmx_sparsematrix_entry_t **  data;                 /**< data[i] is a list with entries on row i  */
}
gmx_sparsematrix_t;


/*! \brief Allocate a new sparse matrix structure
 *
 *  The number of rows is used to allocate the index array entry. Obviously you
 *  can reallocate these later yourself if necessary - this is a
 *  convenience routine.
 *
 *  By default, the compressed_symmetric flag in the structure will
 *  be FALSE. Set it to TRUE manually if you are only storing either the
 *  upper or lower half of the matrix.
 */
gmx_sparsematrix_t *
gmx_sparsematrix_init            (int                    nrow);


/*! \brief Release all resources used by a sparse matrix structure
 *
 *  All arrays in the structure will be freed, and the structure itself.
 */
void
gmx_sparsematrix_destroy         (gmx_sparsematrix_t *   A);


/*! \brief Print sparse matrix to a stream.
 *
 *  Mainly used for debugging. Be warned that the real sparse matrices used
 *  in Gromacs runs can be HUGE (think 100,000 rows).
 */
void
gmx_sparsematrix_print           (FILE *                 stream,
                                  gmx_sparsematrix_t *   A);

/* Adds value at row,col. If the value did not exist
 * previously it is added, otherwise it is incremented with difference.
 *
 * The column sort order might change, so you need to run fix_sparsematrix
 * once you are done changing the matrix.
 */
real
gmx_sparsematrix_value          (gmx_sparsematrix_t *    A,
                                 int                     row,
                                 int                     col);


/* Adds value at row,col. If the value did not exist
 * previously it is added, otherwise it is incremented with difference.
 *
 * The column sort order might change, so you need to run fix_sparsematrix
 * once you are done changing the matrix.
 */
void
gmx_sparsematrix_increment_value(gmx_sparsematrix_t *    A,
                                 int                     row,
                                 int                     col,
                                 real                    difference);



/*! \brief Sort elements in each column and remove zeros.
 *
 *  Sparse matrix access is faster when the elements are stored in
 *  increasing column order in each row. In some cases previously non-zero
 *  elements will be zero after adding more data, and this routine also removes
 *  those entries to reduce the storage requirements.
 *
 *  It never hurts to run this routine if you have been updating the matrix...
 */
void
gmx_sparsematrix_compress       (gmx_sparsematrix_t *    A);



/*! \brief Sparse matrix vector multiplication
 *
 * Calculate y = A * x for a sparse matrix A.
 */
void
gmx_sparsematrix_vector_multiply(gmx_sparsematrix_t *    A,
                                 real *                  x,
                                 real *                  y);

#ifdef __cplusplus
}
#endif


#endif
