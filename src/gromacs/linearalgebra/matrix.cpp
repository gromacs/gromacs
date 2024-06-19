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
#include "gmxpre.h"

#include "gromacs/linearalgebra/matrix.h"

#include "config.h"

#include <cstdio>

#include <filesystem>

#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_lapack.h"

double** alloc_matrix(int n, int m)
{
    double** ptr;
    int      i;

    /* There's always time for more pointer arithmetic! */
    /* This is necessary in order to be able to work with LAPACK */
    snew(ptr, n);
    snew(ptr[0], n * m);
    for (i = 1; (i < n); i++)
    {
        ptr[i] = ptr[i - 1] + m;
    }
    return ptr;
}

void free_matrix(double** a)
{
    sfree(a[0]);
    sfree(a);
}

#define DEBUG_MATRIX
void matrix_multiply(FILE* fp, int n, int m, double** x, double** y, double** z)
{
    int i, j, k;

#ifdef DEBUG_MATRIX
    if (fp)
    {
        std::fprintf(fp, "Multiplying %d x %d matrix with a %d x %d matrix\n", n, m, m, n);
    }
    if (fp)
    {
        for (i = 0; (i < n); i++)
        {
            for (j = 0; (j < m); j++)
            {
                std::fprintf(fp, " %7g", x[i][j]);
            }
            std::fprintf(fp, "\n");
        }
    }
#endif
    for (i = 0; (i < m); i++)
    {
        for (j = 0; (j < m); j++)
        {
            z[i][j] = 0;
            for (k = 0; (k < n); k++)
            {
                z[i][j] += x[k][i] * y[j][k];
            }
        }
    }
}

static void dump_matrix(FILE* fp, const char* title, int n, double** a)
{
    double d = 1;
    int    i, j;

    std::fprintf(fp, "%s\n", title);
    for (i = 0; (i < n); i++)
    {
        d = d * a[i][i];
        for (j = 0; (j < n); j++)
        {
            std::fprintf(fp, " %8.2f", a[i][j]);
        }
        std::fprintf(fp, "\n");
    }
    std::fprintf(fp, "Prod a[i][i] = %g\n", d);
}

int matrix_invert(FILE* fp, int n, double** a)
{
    int      i, j, m, lda, *ipiv, lwork, info;
    double **test = nullptr, **id, *work;

#ifdef DEBUG_MATRIX
    if (fp)
    {
        std::fprintf(fp, "Inverting %d square matrix\n", n);
        test = alloc_matrix(n, n);
        for (i = 0; (i < n); i++)
        {
            for (j = 0; (j < n); j++)
            {
                test[i][j] = a[i][j];
            }
        }
        dump_matrix(fp, "before inversion", n, a);
    }
#endif
    snew(ipiv, n);
    lwork = n * n;
    snew(work, lwork);
    m = lda = n;
    info    = 0;
    F77_FUNC(dgetrf, DGETRF)(&n, &m, a[0], &lda, ipiv, &info);
#ifdef DEBUG_MATRIX
    if (fp)
    {
        dump_matrix(fp, "after dgetrf", n, a);
    }
#endif
    if (info != 0)
    {
        return info;
    }
    F77_FUNC(dgetri, DGETRI)(&n, a[0], &lda, ipiv, work, &lwork, &info);
#ifdef DEBUG_MATRIX
    if (fp)
    {
        dump_matrix(fp, "after dgetri", n, a);
    }
#endif
    if (info != 0)
    {
        return info;
    }

#ifdef DEBUG_MATRIX
    if (fp)
    {
        id = alloc_matrix(n, n);
        matrix_multiply(fp, n, n, test, a, id);
        dump_matrix(fp, "And here is the product of A and Ainv", n, id);
        free_matrix(id);
        free_matrix(test);
    }
#endif
    sfree(ipiv);
    sfree(work);

    return 0;
}

double multi_regression(FILE* fp, int nrow, double* y, int ncol, double** xx, double* a0)
{
    int    row, i, j;
    double ax, chi2, **a, **at, **ata, *atx;

    a   = alloc_matrix(nrow, ncol);
    at  = alloc_matrix(ncol, nrow);
    ata = alloc_matrix(ncol, ncol);
    for (i = 0; (i < nrow); i++)
    {
        for (j = 0; (j < ncol); j++)
        {
            at[j][i] = a[i][j] = xx[j][i];
        }
    }
    matrix_multiply(fp, nrow, ncol, a, at, ata);
    if ((row = matrix_invert(fp, ncol, ata)) != 0)
    {
        gmx_fatal(
                FARGS,
                "Matrix inversion failed. Incorrect row = %d.\nThis probably indicates that you do "
                "not have sufficient data points, or that some parameters are linearly dependent.",
                row);
    }
    snew(atx, ncol);

    for (i = 0; (i < ncol); i++)
    {
        atx[i] = 0;
        for (j = 0; (j < nrow); j++)
        {
            atx[i] += at[i][j] * y[j];
        }
    }
    for (i = 0; (i < ncol); i++)
    {
        a0[i] = 0;
        for (j = 0; (j < ncol); j++)
        {
            a0[i] += ata[i][j] * atx[j];
        }
    }
    chi2 = 0;
    for (j = 0; (j < nrow); j++)
    {
        ax = 0;
        for (i = 0; (i < ncol); i++)
        {
            ax += a0[i] * a[j][i];
        }
        chi2 += (y[j] - ax) * (y[j] - ax);
    }

    sfree(atx);
    free_matrix(a);
    free_matrix(at);
    free_matrix(ata);

    return chi2;
}
