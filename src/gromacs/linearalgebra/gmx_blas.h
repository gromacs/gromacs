/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Header definitions for the standard BLAS library.
 *
 * This is the subset of BLAS routines used for the
 * linear algebra operations in Gromacs.
 * Do NOT use this for other purposes - we only provide this as a
 * simple fallback/reference implementation when no optimized BLAS
 * is present. If you need an implementation for your own code
 * there are several much faster versions out there.
 *
 * All routines are compatible with the BLAS reference implementation,
 * meaning they assume fortran-style matrix row/column organization.
 *
 * There is plenty of documentation for these routines available
 * at http://www.netlib.org/blas , so there is no point in repeating
 * it here.
 */
#ifndef GMX_BLAS_H
#define GMX_BLAS_H

/*! \cond */

#include "config.h"

/* These are not required by this file, but by the internal BLAS
 * implementation.  In principle, they could be included in each file
 * that requires them, but this is simpler.  Since the header is internal
 * to the linearyalgebra/ module, the added complexity may not be worth it. */
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* Double precision versions */
double
    F77_FUNC(dasum, DASUM) (int *n, double *dx, int *incx);

void
    F77_FUNC(daxpy, DAXPY) (int *n, double *da, double *dx, int *incx, double *dy, int *incy);

void
    F77_FUNC(dcopy, DCOPY) (int *n, double *dx, int *incx, double *dy, int *incy);

double
    F77_FUNC(ddot, DDOT) (int *n, double *dx, int *incx, double *dy, int *incy);

void
    F77_FUNC(dgemm, DGEMM) (const char *transa, const char *transb, int *m, int *n, int *k,
                            double *alpha, double *a, int *lda, double *b, int *ldb,
                            double *beta, double *c, int *ldc);

void
    F77_FUNC(dgemv, DGEMV) (const char *trans, int *m, int *n, double *alpha, double *a, int *lda,
                            double *x, int *incx, double *beta, double *y, int *incy);

void
    F77_FUNC(dger, DGER) (int *m, int *n, double *alpha, double *x, int *incx,
                          double *y, int *incy, double *a, int *lda);

double
    F77_FUNC(dnrm2, DNRM2) (int  *n, double *x, int *incx);

void
    F77_FUNC(drot, DROT) (int *n, double *dx, int *incx,
                          double *dy, int *incy, double *c, double *s);

void
    F77_FUNC(dscal, DSCAL) (int *n, double *fact, double *dx, int *incx);

void
    F77_FUNC(dswap, DSWAP) (int *n, double *dx, int *incx, double *dy, int *incy);

void
    F77_FUNC(dsymv, DSYMV) (const char *uplo, int *n, double *alpha, double *a, int *lda,
                            double *x, int *incx, double *beta, double *y, int *incy);

void
    F77_FUNC(dsyr2, DSYR2) (const char *uplo, int *n, double *alpha, double *x, int *incx,
                            double *y, int *incy, double *a, int *lda);

void
    F77_FUNC(dsyr2k, DSYR2K) (const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a,
                              int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

void
    F77_FUNC(dtrmm, DTRMM) (const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
                            double *alpha, double *a, int *lda, double *b, int *ldb);

void
    F77_FUNC(dtrmv, DTRMV) (const char *uplo, const char *trans, const char *diag, int *n,
                            double *a, int *lda, double *x, int *incx);

void
    F77_FUNC(dtrsm, DTRSM) (const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
                            double *alpha, double *a, int *lda, double *b, int *ldb);

int
    F77_FUNC(idamax, IDAMAX) (int *n, double *dx, int *incx);



/* Single precision versions */
float
    F77_FUNC(sasum, SASUM) (int *n, float *dx, int *incx);

void
    F77_FUNC(saxpy, SAXPY) (int *n, float *da, float *dx, int *incx, float *dy, int *incy);

void
    F77_FUNC(scopy, SCOPY) (int *n, float *dx, int *incx, float *dy, int *incy);

float
    F77_FUNC(sdot, SDOT) (int *n, float *dx, int *incx, float *dy, int *incy);

void
    F77_FUNC(sgemm, SGEMM) (const char *transa, const char *transb, int *m, int *n, int *k,
                            float *alpha, float *a, int *lda, float *b, int *ldb,
                            float *beta, float *c, int *ldc);

void
    F77_FUNC(sgemv, SGEMV) (const char *trans, int *m, int *n, float *alpha, float *a, int *lda,
                            float *x, int *incx, float *beta, float *y, int *incy);

void
    F77_FUNC(sger, SGER) (int *m, int *n, float *alpha, float *x, int *incx,
                          float *y, int *incy, float *a, int *lda);

float
    F77_FUNC(snrm2, SNRM2) (int  *n, float *x, int *incx);

void
    F77_FUNC(srot, SROT) (int *n, float *dx, int *incx,
                          float *dy, int *incy, float *c, float *s);

void
    F77_FUNC(sscal, SSCAL) (int *n, float *fact, float *dx, int *incx);

void
    F77_FUNC(sswap, SSWAP) (int *n, float *dx, int *incx, float *dy, int *incy);

void
    F77_FUNC(ssymv, SSYMV) (const char *uplo, int *n, float *alpha, float *a, int *lda,
                            float *x, int *incx, float *beta, float *y, int *incy);

void
    F77_FUNC(ssyr2, SSYR2) (const char *uplo, int *n, float *alpha, float *x, int *incx,
                            float *y, int *incy, float *a, int *lda);

void
    F77_FUNC(ssyr2k, SSYR2K) (const char *uplo, const char *trans, int *n, int *k, float *alpha, float *a,
                              int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

void
    F77_FUNC(strmm, STRMM) (const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
                            float *alpha, float *a, int *lda, float *b, int *ldb);

void
    F77_FUNC(strmv, STRMV) (const char *uplo, const char *trans, const char *diag, int *n,
                            float *a, int *lda, float *x, int *incx);

void
    F77_FUNC(strsm, STRSM) (const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
                            float *alpha, float *a, int *lda, float *b, int *ldb);

int
    F77_FUNC(isamax, ISAMAX) (int *n, float *dx, int *incx);


#ifdef __cplusplus
}
#endif

/*! \endcond */

#endif /* GMX_BLAS_H */
