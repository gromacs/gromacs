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

#include "gromacs/linearalgebra/eigensolver.h"

#include "config.h"

#include <cstdio>

#include <filesystem>

#include "gromacs/linearalgebra/sparsematrix.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_arpack.h"
#include "gmx_lapack.h"

void eigensolver(real* a, int n, int index_lower, int index_upper, real* eigenvalues, real* eigenvectors)
{
    int*        isuppz;
    int         lwork, liwork;
    int         m, iw0, info;
    real        w0, abstol;
    int*        iwork;
    real*       work;
    real        vl, vu;
    const char* jobz;

    if (index_lower < 0)
    {
        index_lower = 0;
    }

    if (index_upper >= n)
    {
        index_upper = n - 1;
    }

    /* Make jobz point to the character "V" if eigenvectors
     * should be calculated, otherwise "N" (only eigenvalues).
     */
    jobz = (eigenvectors != nullptr) ? "V" : "N";

    /* allocate lapack stuff */
    snew(isuppz, 2 * n);
    vl = vu = 0;

    /* First time we ask the routine how much workspace it needs */
    lwork  = -1;
    liwork = -1;
    abstol = 0;

    /* Convert indices to fortran standard */
    index_lower++;
    index_upper++;

    /* Call LAPACK routine using fortran interface. Note that we use upper storage,
     * but this corresponds to lower storage ("L") in Fortran.
     */
#if GMX_DOUBLE
    F77_FUNC(dsyevr, DSYEVR)
    (jobz,
     "I",
     "L",
     &n,
     a,
     &n,
     &vl,
     &vu,
     &index_lower,
     &index_upper,
     &abstol,
     &m,
     eigenvalues,
     eigenvectors,
     &n,
     isuppz,
     &w0,
     &lwork,
     &iw0,
     &liwork,
     &info);
#else
    F77_FUNC(ssyevr, SSYEVR)
    (jobz,
     "I",
     "L",
     &n,
     a,
     &n,
     &vl,
     &vu,
     &index_lower,
     &index_upper,
     &abstol,
     &m,
     eigenvalues,
     eigenvectors,
     &n,
     isuppz,
     &w0,
     &lwork,
     &iw0,
     &liwork,
     &info);
#endif

    if (info != 0)
    {
        sfree(isuppz);
        gmx_fatal(FARGS, "Internal error in LAPACK diagonalization.");
    }

    lwork  = static_cast<int>(w0);
    liwork = iw0;

    snew(work, lwork);
    snew(iwork, liwork);

    abstol = 0;

#if GMX_DOUBLE
    F77_FUNC(dsyevr, DSYEVR)
    (jobz,
     "I",
     "L",
     &n,
     a,
     &n,
     &vl,
     &vu,
     &index_lower,
     &index_upper,
     &abstol,
     &m,
     eigenvalues,
     eigenvectors,
     &n,
     isuppz,
     work,
     &lwork,
     iwork,
     &liwork,
     &info);
#else
    F77_FUNC(ssyevr, SSYEVR)
    (jobz,
     "I",
     "L",
     &n,
     a,
     &n,
     &vl,
     &vu,
     &index_lower,
     &index_upper,
     &abstol,
     &m,
     eigenvalues,
     eigenvectors,
     &n,
     isuppz,
     work,
     &lwork,
     iwork,
     &liwork,
     &info);
#endif

    sfree(isuppz);
    sfree(work);
    sfree(iwork);

    if (info != 0)
    {
        gmx_fatal(FARGS, "Internal error in LAPACK diagonalization.");
    }
}


#ifdef GMX_MPI_NOT
void sparse_parallel_eigensolver(gmx_sparsematrix_t* A, int neig, real* eigenvalues, real* eigenvectors, int maxiter)
{
    int   iwork[80];
    int   iparam[11];
    int   ipntr[11];
    real* resid;
    real* workd;
    real* workl;
    real* v;
    int   n;
    int   ido, info, lworkl, i, ncv, dovec;
    real  abstol;
    int*  select;
    int   iter;
    int   nnodes, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (eigenvectors != NULL)
    {
        dovec = 1;
    }
    else
    {
        dovec = 0;
    }

    n   = A->nrow;
    ncv = 2 * neig;

    if (ncv > n)
    {
        ncv = n;
    }

    for (i = 0; i < 11; i++)
    {
        iparam[i] = ipntr[i] = 0;
    }

    iparam[0] = 1;       /* Don't use explicit shifts */
    iparam[2] = maxiter; /* Max number of iterations */
    iparam[6] = 1;       /* Standard symmetric eigenproblem */

    lworkl = ncv * (8 + ncv);
    snew(resid, n);
    snew(workd, (3 * n + 4));
    snew(workl, lworkl);
    snew(select, ncv);
    snew(v, n * ncv);

    /* Use machine tolerance - roughly 1e-16 in double precision */
    abstol = 0;

    ido = info = 0;
    fprintf(stderr, "Calculation Ritz values and Lanczos vectors, max %d iterations...\n", maxiter);

    iter = 1;
    do
    {
#    if GMX_DOUBLE
        F77_FUNC(pdsaupd, PDSAUPD)
        (&ido, "I", &n, "SA", &neig, &abstol, resid, &ncv, v, &n, iparam, ipntr, workd, iwork, workl, &lworkl, &info);
#    else
        F77_FUNC(pssaupd, PSSAUPD)
        (&ido, "I", &n, "SA", &neig, &abstol, resid, &ncv, v, &n, iparam, ipntr, workd, iwork, workl, &lworkl, &info);
#    endif
        if (ido == -1 || ido == 1)
        {
            gmx_sparsematrix_vector_multiply(A, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
        }

        fprintf(stderr, "\rIteration %4d: %3d out of %3d Ritz values converged.", iter++, iparam[4], neig);
        fflush(stderr);
    } while (info == 0 && (ido == -1 || ido == 1));

    fprintf(stderr, "\n");
    if (info == 1)
    {
        gmx_fatal(FARGS,
                  "Maximum number of iterations (%d) reached in Arnoldi\n"
                  "diagonalization, but only %d of %d eigenvectors converged.\n",
                  maxiter,
                  iparam[4],
                  neig);
    }
    else if (info != 0)
    {
        gmx_fatal(FARGS, "Unspecified error from Arnoldi diagonalization:%d\n", info);
    }

    info = 0;
    /* Extract eigenvalues and vectors from data */
    fprintf(stderr, "Calculating eigenvalues and eigenvectors...\n");

#    if GMX_DOUBLE
    F77_FUNC(pdseupd, PDSEUPD)
    (&dovec,
     "A",
     select,
     eigenvalues,
     eigenvectors,
     &n,
     NULL,
     "I",
     &n,
     "SA",
     &neig,
     &abstol,
     resid,
     &ncv,
     v,
     &n,
     iparam,
     ipntr,
     workd,
     workl,
     &lworkl,
     &info);
#    else
    F77_FUNC(psseupd, PSSEUPD)
    (&dovec,
     "A",
     select,
     eigenvalues,
     eigenvectors,
     &n,
     NULL,
     "I",
     &n,
     "SA",
     &neig,
     &abstol,
     resid,
     &ncv,
     v,
     &n,
     iparam,
     ipntr,
     workd,
     workl,
     &lworkl,
     &info);
#    endif

    sfree(v);
    sfree(resid);
    sfree(workd);
    sfree(workl);
    sfree(select);
}
#endif


void sparse_eigensolver(gmx_sparsematrix_t* A, int neig, real* eigenvalues, real* eigenvectors, int maxiter)
{
    int   iwork[80];
    int   iparam[11];
    int   ipntr[11];
    real* resid;
    real* workd;
    real* workl;
    real* v;
    int   n;
    int   ido, info, lworkl, i, ncv, dovec;
    real  abstol;
    int*  select;
    int   iter;

#ifdef GMX_MPI_NOT
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    if (n > 1)
    {
        sparse_parallel_eigensolver(A, neig, eigenvalues, eigenvectors, maxiter);
        return;
    }
#endif

    if (eigenvectors != nullptr)
    {
        dovec = 1;
    }
    else
    {
        dovec = 0;
    }

    n   = A->nrow;
    ncv = 2 * neig;

    if (ncv > n)
    {
        ncv = n;
    }

    for (i = 0; i < 11; i++)
    {
        iparam[i] = ipntr[i] = 0;
    }

    iparam[0] = 1;       /* Don't use explicit shifts */
    iparam[2] = maxiter; /* Max number of iterations */
    iparam[6] = 1;       /* Standard symmetric eigenproblem */

    lworkl = ncv * (8 + ncv);
    snew(resid, n);
    snew(workd, (3 * n + 4));
    snew(workl, lworkl);
    snew(select, ncv);
    snew(v, n * ncv);

    /* Use machine tolerance - roughly 1e-16 in double precision */
    abstol = 0;

    ido = info = 0;
    fprintf(stderr, "Calculation Ritz values and Lanczos vectors, max %d iterations...\n", maxiter);

    iter = 1;
    do
    {
#if GMX_DOUBLE
        F77_FUNC(dsaupd, DSAUPD)
        (&ido, "I", &n, "SA", &neig, &abstol, resid, &ncv, v, &n, iparam, ipntr, workd, iwork, workl, &lworkl, &info);
#else
        F77_FUNC(ssaupd, SSAUPD)
        (&ido, "I", &n, "SA", &neig, &abstol, resid, &ncv, v, &n, iparam, ipntr, workd, iwork, workl, &lworkl, &info);
#endif
        if (ido == -1 || ido == 1)
        {
            gmx_sparsematrix_vector_multiply(A, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
        }

        fprintf(stderr, "\rIteration %4d: %3d out of %3d Ritz values converged.", iter++, iparam[4], neig);
        fflush(stderr);
    } while (info == 0 && (ido == -1 || ido == 1));

    fprintf(stderr, "\n");
    if (info == 1)
    {
        gmx_fatal(FARGS,
                  "Maximum number of iterations (%d) reached in Arnoldi\n"
                  "diagonalization, but only %d of %d eigenvectors converged.\n",
                  maxiter,
                  iparam[4],
                  neig);
    }
    else if (info != 0)
    {
        gmx_fatal(FARGS, "Unspecified error from Arnoldi diagonalization:%d\n", info);
    }

    info = 0;
    /* Extract eigenvalues and vectors from data */
    fprintf(stderr, "Calculating eigenvalues and eigenvectors...\n");

#if GMX_DOUBLE
    F77_FUNC(dseupd, DSEUPD)
    (&dovec,
     "A",
     select,
     eigenvalues,
     eigenvectors,
     &n,
     nullptr,
     "I",
     &n,
     "SA",
     &neig,
     &abstol,
     resid,
     &ncv,
     v,
     &n,
     iparam,
     ipntr,
     workd,
     workl,
     &lworkl,
     &info);
#else
    F77_FUNC(sseupd, SSEUPD)
    (&dovec,
     "A",
     select,
     eigenvalues,
     eigenvectors,
     &n,
     nullptr,
     "I",
     &n,
     "SA",
     &neig,
     &abstol,
     resid,
     &ncv,
     v,
     &n,
     iparam,
     ipntr,
     workd,
     workl,
     &lworkl,
     &info);
#endif

    sfree(v);
    sfree(resid);
    sfree(workd);
    sfree(workl);
    sfree(select);
}
