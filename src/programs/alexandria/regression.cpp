/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "regression.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

extern "C" void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
                        double* b, int* ldb, double* s, double* rcond, int* rank,
                        double* work, int* lwork, int* iwork, int* info );

void multi_regression2(int nrow, double y[], int ncol,
                       double **a, double x[])
{
    /* Query and allocate the optimal workspace */
    int     lwork = -1;
    int     lda   = nrow;
    int     ldb   = nrow;
    int     nrhs  = 1;
    int     rank;
    double  rcond = -1.0;
    double  wkopt;
    double *s;
    snew(s, nrow);
    // Compute length of integer array iwork according to
    // https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
    int  smlsiz = 25;
    int  nlvl   = std::max(0L, std::lround(std::log2(std::min(nrow, ncol)/smlsiz + 1) ) + 1);
    int  liwork = 3*std::min(ncol, nrow)*nlvl + 11*std::min(nrow, ncol);
    int *iwork;
    snew(iwork, liwork);
    int  info;
    dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s, &rcond, &rank, &wkopt, &lwork,
             iwork, &info );
    lwork = (int)wkopt;
    double *work;
    snew(work, lwork);
    /* Solve the equations A*X = B */
    dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, y, &ldb, s, &rcond, &rank, work, &lwork,
             iwork, &info );
    /* Check for convergence */
    if (info > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "The algorithm computing SVD failed to converge and the least squares solution could not be computed. Info = %d.", info);
        GMX_THROW(gmx::InvalidInputError(buf));
    }
    for (int i = 0; i < ncol; i++)
    {
        x[i] = y[i];
    }

    sfree(s);
    sfree(work);
    sfree(iwork);
}
