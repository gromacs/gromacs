/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _GMX_BLAS_H_
#define _GMX_BLAS_H_

#include "types/simple.h"

/* Suppress Cygwin compiler warnings from using newlib version of
 * ctype.h */
#ifdef GMX_CYGWIN
#undef toupper
#endif


/** @file
 *
 *  @brief Header definitions for the standard BLAS library.
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

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#ifndef FortranCInterface_GLOBAL
#define FortranCInterface_GLOBAL(name,NAME) name ## _
#endif


/* Double precision versions */
double
FortranCInterface_GLOBAL(dasum,DASUM)(int *n, double *dx, int *incx);

void
FortranCInterface_GLOBAL(daxpy,DAXPY)(int *n, double *da, double *dx, int *incx, double *dy, int *incy);

void
FortranCInterface_GLOBAL(dcopy,DCOPY)(int *n, double *dx, int *incx, double *dy, int *incy);

double
FortranCInterface_GLOBAL(ddot,DDOT)(int *n, double *dx, int *incx, double *dy, int *incy);

void
FortranCInterface_GLOBAL(dgemm,DGEMM)(const char *transa, const char *transb, int *m, int *n, int *k, 
       double *alpha, double *a, int *lda, double *b, int *ldb, 
       double *beta, double *c, int *ldc);

void
FortranCInterface_GLOBAL(dgemv,DGEMV)(const char *trans, int *m, int *n, double *alpha, double *a, int *lda,
       double *x, int *incx, double *beta, double *y, int *incy);

void
FortranCInterface_GLOBAL(dger,DGER)(int *m, int *n, double *alpha, double *x, int *incx, 
      double *y, int *incy, double *a, int *lda);

double
FortranCInterface_GLOBAL(dnrm2,DNRM2)(int  *n, double *x, int *incx);

void
FortranCInterface_GLOBAL(drot,DROT)(int *n, double *dx, int *incx, 
      double *dy, int *incy, double *c, double *s);

void 
FortranCInterface_GLOBAL(dscal,DSCAL)(int *n, double *fact, double *dx, int *incx);

void
FortranCInterface_GLOBAL(dswap,DSWAP)(int *n, double *dx, int *incx, double *dy, int *incy);

void
FortranCInterface_GLOBAL(dsymv,DSYMV)(const char *uplo, int *n, double *alpha, double *a, int *lda,
       double *x, int *incx, double *beta, double *y, int *incy);

void
FortranCInterface_GLOBAL(dsyr2,DSYR2)(const char *uplo, int *n, double *alpha, double *x, int *incx,
       double *y, int *incy, double *a, int *lda);

void
FortranCInterface_GLOBAL(dsyr2k,DSYR2K)(const char *uplo, const char *trans, int *n, int *k, double *alpha, double *a,
        int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

void 
FortranCInterface_GLOBAL(dtrmm,DTRMM)(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, 
       double *alpha, double *a, int *lda, double *b, int *ldb);

void 
FortranCInterface_GLOBAL(dtrmv,DTRMV)(const char *uplo, const char *trans, const char *diag, int *n, 
       double *a, int *lda, double *x, int *incx);

void
FortranCInterface_GLOBAL(dtrsm,DTRSM)(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
       double *alpha, double *a,int *lda, double *b, int *ldb);

int
FortranCInterface_GLOBAL(idamax,IDAMAX)(int *n, double *dx, int *incx);



/* Single precision versions */
float
FortranCInterface_GLOBAL(sasum,SASUM)(int *n, float *dx, int *incx);

void
FortranCInterface_GLOBAL(saxpy,SAXPY)(int *n, float *da, float *dx, int *incx, float *dy, int *incy);

void
FortranCInterface_GLOBAL(scopy,SCOPY)(int *n, float *dx, int *incx, float *dy, int *incy);

float
FortranCInterface_GLOBAL(sdot,SDOT)(int *n, float *dx, int *incx, float *dy, int *incy);

void
FortranCInterface_GLOBAL(sgemm,SGEMM)(const char *transa, const char *transb, int *m, int *n, int *k, 
       float *alpha, float *a, int *lda, float *b, int *ldb, 
       float *beta, float *c, int *ldc);

void
FortranCInterface_GLOBAL(sgemv,SGEMV)(const char *trans, int *m, int *n, float *alpha, float *a, int *lda,
       float *x, int *incx, float *beta, float *y, int *incy);

void
FortranCInterface_GLOBAL(sger,SGER)(int *m, int *n, float *alpha, float *x, int *incx, 
      float *y, int *incy, float *a, int *lda);

float
FortranCInterface_GLOBAL(snrm2,SNRM2)(int  *n, float *x, int *incx);

void
FortranCInterface_GLOBAL(srot,SROT)(int *n, float *dx, int *incx, 
      float *dy, int *incy, float *c, float *s);

void 
FortranCInterface_GLOBAL(sscal,SSCAL)(int *n, float *fact, float *dx, int *incx);

void
FortranCInterface_GLOBAL(sswap,SSWAP)(int *n, float *dx, int *incx, float *dy, int *incy);

void
FortranCInterface_GLOBAL(ssymv,SSYMV)(const char *uplo, int *n, float *alpha, float *a, int *lda,
       float *x, int *incx, float *beta, float *y, int *incy);

void
FortranCInterface_GLOBAL(ssyr2,SSYR2)(const char *uplo, int *n, float *alpha, float *x, int *incx,
       float *y, int *incy, float *a, int *lda);

void
FortranCInterface_GLOBAL(ssyr2k,SSYR2K)(const char *uplo, const char *trans, int *n, int *k, float *alpha, float *a,
        int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);

void 
FortranCInterface_GLOBAL(strmm,STRMM)(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n, 
       float *alpha, float *a, int *lda, float *b, int *ldb);

void 
FortranCInterface_GLOBAL(strmv,STRMV)(const char *uplo, const char *trans, const char *diag, int *n, 
       float *a, int *lda, float *x, int *incx);

void
FortranCInterface_GLOBAL(strsm,STRSM)(const char *side, const char *uplo, const char *transa, const char *diag, int *m, int *n,
       float *alpha, float *a,int *lda, float *b, int *ldb);

int
FortranCInterface_GLOBAL(isamax,ISAMAX)(int *n, float *dx, int *incx);


#ifdef __cplusplus
}
#endif



#endif /* _BLAS_H_ */
