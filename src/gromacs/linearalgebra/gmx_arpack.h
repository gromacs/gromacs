/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2004 David van der Spoel, Erik Lindahl, University of Groningen.
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
 * Selected routines from ARPACK
 *
 * This file contains a subset of ARPACK functions to perform
 * diagonalization and SVD for sparse matrices in Gromacs.
 *
 * Consult the main ARPACK site for detailed documentation:
 * http://www.caam.rice.edu/software/ARPACK/
 *
 * Below, we just list the options and any specific differences
 * from ARPACK. The code is essentially the same, but the routines
 * have been made thread-safe by using extra workspace arrays.
 */
#ifndef GMX_ARPACK_H
#define GMX_ARPACK_H

#include "config.h"
#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Implicitly Restarted Arnoldi Iteration, double precision.
 *
 *  Reverse communication interface for the Implicitly Restarted Arnoldi
 *  Iteration.  For symmetric problems this reduces to a variant of the
 *  Lanczos method. See the ARPACK site for details.
 *
 *  \param ido     Reverse communication flag. Set to 0 first time.
 *                 Upon return with ido=-1 or ido=1 you should calculate
 *                 Y=A*X and recall the routine. Return with ido=2 means
 *                 Y=B*X should be calculated. ipntr[0] is the pointer in
 *                 workd for X, ipntr[1] is the index for Y.
 *                 Return with ido=99 means it finished.
 *  \param bmat    'I' for standard eigenproblem, 'G' for generalized.
 *  \param n       Order of eigenproblem.
 *  \param which   Which eigenvalues to calculate. 'LA' for largest
 *                 algebraic, 'SA' for smallest algebraic, 'LM' for largest
 *                 magnitude, 'SM' for smallest magnitude, and finally
 *                 'BE' (both ends) to calculate half from each end of
 *                 the spectrum.
 *  \param nev     Number of eigenvalues to calculate. 0<nev<n.
 *  \param tol     Tolerance. Machine precision of it is 0.
 *  \param resid   Optional starting residual vector at input if info=1,
 *                 otherwise a random one is used. Final residual vector on
 *                 return.
 *  \param ncv     Number of columns in matrix v.
 *  \param v       N*NCV matrix. V contain the Lanczos basis vectors.
 *  \param ldv     Leading dimension of v.
 *  \param iparam  Integer array, size 11. Same contents as arpack.
 *  \param ipntr   Integer array, size 11. Points to starting locations
 *                 in the workd/workl arrays. Same contents as arpack.
 *  \param workd   Double precision work array, length 3*n+4.
 *                 Provide the same array for all calls, and don't touch it.
 *                 IMPORTANT: This is 4 units larger than standard ARPACK!
 *  \param iwork   Integer work array, size 80.
 *                 Provide the same array for all calls, and don't touch it.
 *                 IMPORTANT: New argument compared to standard ARPACK!
 *  \param workl   Double precision work array, length lwork.
 *  \param lworkl  Length of the work array workl. Must be at least ncv*(ncv+8)
 *  \param info    Set info to 0 to use random initial residual vector,
 *                 or to 1 if you provide a one. On output, info=0 means
 *                 normal exit, 1 that max number of iterations was reached,
 *                 and 3 that no shifts could be applied. Negative numbers
 *                 correspond to errors in the arguments provided.
 */
void
    F77_FUNC(dsaupd, DSAUPD) (int *     ido,
                              const char *    bmat,
                              int *     n,
                              const char *      which,
                              int *     nev,
                              double *  tol,
                              double *  resid,
                              int *     ncv,
                              double *  v,
                              int *     ldv,
                              int *     iparam,
                              int *     ipntr,
                              double *  workd,
                              int *     iwork,
                              double *  workl,
                              int *     lworkl,
                              int *     info);



/*! \brief Get eigenvalues/vectors after Arnoldi iteration, double prec.
 *
 *  See the ARPACK site for details. You must have finished the interative
 *  part with dsaupd() before calling this function.
 *
 *  \param rvec    1 if you want eigenvectors, 0 if not.
 *  \param howmny  'A' if you want all nvec vectors, 'S' if you
 *                 provide a subset selection in select[].
 *  \param select  Integer array, dimension nev. Indices of the
 *                 eigenvectors to calculate. Fortran code means we
 *                 start counting on 1. This array must be given even in
 *                 howmny is 'A'. (Arpack documentation is wrong on this).
 *  \param d       Double precision array, length nev. Eigenvalues.
 *  \param z       Double precision array, n*nev. Eigenvectors.
 *  \param ldz     Leading dimension of z. Normally n.
 *  \param sigma   Shift if iparam[6] is 3,4, or 5. Ignored otherwise.
 *  \param bmat    Provide the same argument as you did to dsaupd()
 *  \param n       Provide the same argument as you did to dsaupd()
 *  \param which   Provide the same argument as you did to dsaupd()
 *  \param nev     Provide the same argument as you did to dsaupd()
 *  \param tol     Provide the same argument as you did to dsaupd()
 *  \param resid   Provide the same argument as you did to dsaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ncv     Provide the same argument as you did to dsaupd()
 *  \param v       Provide the same argument as you did to dsaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ldv     Provide the same argument as you did to dsaupd()
 *  \param iparam  Provide the same argument as you did to dsaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ipntr   Provide the same argument as you did to dsaupd()
 *                 The array must not be touched between the two function calls!
 *  \param workd   Provide the same argument as you did to dsaupd()
 *                 The array must not be touched between the two function calls!
 *  \param workl   Double precision work array, length lwork.
 *                 The array must not be touched between the two function calls!
 *  \param lworkl  Provide the same argument as you did to dsaupd()
 *  \param info    Provide the same argument as you did to dsaupd()
 */
void
    F77_FUNC(dseupd, DSEUPD) (int *     rvec,
                              const char *    howmny,
                              int *     select,
                              double *  d,
                              double *  z,
                              int *     ldz,
                              double *  sigma,
                              const char *    bmat,
                              int *     n,
                              const char *    which,
                              int *     nev,
                              double *  tol,
                              double *  resid,
                              int *     ncv,
                              double *  v,
                              int *     ldv,
                              int *     iparam,
                              int *     ipntr,
                              double *  workd,
                              double *  workl,
                              int *     lworkl,
                              int *     info);





/*! \brief Implicitly Restarted Arnoldi Iteration, single precision.
 *
 *  Reverse communication interface for the Implicitly Restarted Arnoldi
 *  Iteration.  For symmetric problems this reduces to a variant of the
 *  Lanczos method. See the ARPACK site for details.
 *
 *  \param ido     Reverse communication flag. Set to 0 first time.
 *                 Upon return with ido=-1 or ido=1 you should calculate
 *                 Y=A*X and recall the routine. Return with ido=2 means
 *                 Y=B*X should be calculated. ipntr[0] is the pointer in
 *                 workd for X, ipntr[1] is the index for Y.
 *                 Return with ido=99 means it finished.
 *  \param bmat    'I' for standard eigenproblem, 'G' for generalized.
 *  \param n       Order of eigenproblem.
 *  \param which   Which eigenvalues to calculate. 'LA' for largest
 *                 algebraic, 'SA' for smallest algebraic, 'LM' for largest
 *                 magnitude, 'SM' for smallest magnitude, and finally
 *                 'BE' (both ends) to calculate half from each end of
 *                 the spectrum.
 *  \param nev     Number of eigenvalues to calculate. 0<nev<n.
 *  \param tol     Tolerance. Machine precision of it is 0.
 *  \param resid   Optional starting residual vector at input if info=1,
 *                 otherwise a random one is used. Final residual vector on
 *                 return.
 *  \param ncv     Number of columns in matrix v.
 *  \param v       N*NCV matrix. V contain the Lanczos basis vectors.
 *  \param ldv     Leading dimension of v.
 *  \param iparam  Integer array, size 11. Same contents as arpack.
 *  \param ipntr   Integer array, size 11. Points to starting locations
 *                 in the workd/workl arrays. Same contents as arpack.
 *  \param workd   Single precision work array, length 3*n+4.
 *                 Provide the same array for all calls, and don't touch it.
 *                 IMPORTANT: This is 4 units larger than standard ARPACK!
 *  \param iwork   Integer work array, size 80.
 *                 Provide the same array for all calls, and don't touch it.
 *                 IMPORTANT: New argument compared to standard ARPACK!
 *  \param workl   Single precision work array, length lwork.
 *  \param lworkl  Length of the work array workl. Must be at least ncv*(ncv+8)
 *  \param info    Set info to 0 to use random initial residual vector,
 *                 or to 1 if you provide a one. On output, info=0 means
 *                 normal exit, 1 that max number of iterations was reached,
 *                 and 3 that no shifts could be applied. Negative numbers
 *                 correspond to errors in the arguments provided.
 */
void
    F77_FUNC(ssaupd, SSAUPD) (int *     ido,
                              const char *    bmat,
                              int *     n,
                              const char *    which,
                              int *     nev,
                              float *   tol,
                              float *   resid,
                              int *     ncv,
                              float *   v,
                              int *     ldv,
                              int *     iparam,
                              int *     ipntr,
                              float *   workd,
                              int *     iwork,
                              float *   workl,
                              int *     lworkl,
                              int *     info);





/*! \brief Get eigenvalues/vectors after Arnoldi iteration, single prec.
 *
 *  See the ARPACK site for details. You must have finished the interative
 *  part with ssaupd() before calling this function.
 *
 *  \param rvec    1 if you want eigenvectors, 0 if not.
 *  \param howmny  'A' if you want all nvec vectors, 'S' if you
 *                 provide a subset selection in select[].
 *  \param select  Integer array, dimension nev. Indices of the
 *                 eigenvectors to calculate. Fortran code means we
 *                 start counting on 1. This array must be given even in
 *                 howmny is 'A'. (Arpack documentation is wrong on this).
 *  \param d       Single precision array, length nev. Eigenvalues.
 *  \param z       Single precision array, n*nev. Eigenvectors.
 *  \param ldz     Leading dimension of z. Normally n.
 *  \param sigma   Shift if iparam[6] is 3,4, or 5. Ignored otherwise.
 *  \param bmat    Provide the same argument as you did to ssaupd()
 *  \param n       Provide the same argument as you did to ssaupd()
 *  \param which   Provide the same argument as you did to ssaupd()
 *  \param nev     Provide the same argument as you did to ssaupd()
 *  \param tol     Provide the same argument as you did to ssaupd()
 *  \param resid   Provide the same argument as you did to ssaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ncv     Provide the same argument as you did to ssaupd()
 *  \param v       Provide the same argument as you did to ssaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ldv     Provide the same argument as you did to ssaupd()
 *  \param iparam  Provide the same argument as you did to ssaupd()
 *                 The array must not be touched between the two function calls!
 *  \param ipntr   Provide the same argument as you did to ssaupd()
 *                 The array must not be touched between the two function calls!
 *  \param workd   Provide the same argument as you did to ssaupd()
 *                 The array must not be touched between the two function calls!
 *  \param workl   Single precision work array, length lwork.
 *                 The array must not be touched between the two function calls!
 *  \param lworkl  Provide the same argument as you did to ssaupd()
 *  \param info    Provide the same argument as you did to ssaupd()
 */
void
    F77_FUNC(sseupd, SSEUPD) (int *     rvec,
                              const char *    howmny,
                              int *     select,
                              float *   d,
                              float *   z,
                              int *     ldz,
                              float *   sigma,
                              const char *    bmat,
                              int *     n,
                              const char *    which,
                              int *     nev,
                              float *   tol,
                              float *   resid,
                              int *     ncv,
                              float *   v,
                              int *     ldv,
                              int *     iparam,
                              int *     ipntr,
                              float *   workd,
                              float *   workl,
                              int *     lworkl,
                              int *     info);

#ifdef __cplusplus
}
#endif

#endif
