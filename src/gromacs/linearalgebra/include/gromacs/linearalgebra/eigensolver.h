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
#ifndef GMX_LINEARALGEBRA_EIGENSOLVER_H
#define GMX_LINEARALGEBRA_EIGENSOLVER_H

#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/utility/real.h"

/** Calculate eigenvalues/vectors a matrix stored in linear memory (not sparse).
 *
 *  This routine uses lapack to diagonalize a real symmetric matrix efficiently,
 *  and the eigenvalues/vectors will be sorted in ascending order on output.
 *  Gromacs comes with a built-in portable BLAS/LAPACK, but if performance
 *  matters it is advisable to link with an optimized vendor-provided library.
 *  Note that this method uses the LAPACK routine SYEVR, which reduces the
 *  matrix to an upper triangular form and calculates eigenvalues/vectors for
 *  the submatrix, and this assumes symmetry of the input matrix.
 *
 *  \param a            Pointer to matrix data, total size n*n
 *                      The input data in the matrix will be destroyed/changed.
 *  \param n            Side of the matrix to calculate eigenvalues for.
 *  \param index_lower  Index of first eigenvector to determine.
 *  \param index_upper  Last eigenvector determined is index_upper-1.
 *  \param eigenvalues  Array of the eigenvalues on return. The length
 *                      of this array _must_ be n, even if not all
 *                      eigenvectors are calculated, since all eigenvalues
 *                      might be needed as an intermediate step.
 *  \param eigenvec     If this pointer is non-NULL, the eigenvectors
 *                      specified by the indices are returned as rows of
 *                      a matrix, i.e. eigenvector j starts at offset j*n, and
 *                      is of length n.
 */
void eigensolver(real* a, int n, int index_lower, int index_upper, real* eigenvalues, real* eigenvec);


/*! \brief Sparse matrix eigensolver.
 *
 *  This routine is intended for large matrices that might not fit in memory.
 *
 *  It will determine the neig lowest eigenvalues, and if the eigenvectors pointer
 *  is non-NULL also the corresponding eigenvectors.
 *
 *  maxiter=100000 should suffice in most cases!
 */
void sparse_eigensolver(gmx_sparsematrix_t* A, int neig, real* eigenvalues, real* eigenvectors, int maxiter);

#endif
