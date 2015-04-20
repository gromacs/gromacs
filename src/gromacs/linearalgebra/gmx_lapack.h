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
 * Header definitions for the standard LAPACK library.
 *
 * This is the subset of LAPACK routines used for the
 * linear algebra operations in Gromacs. Most of the execution time
 * will be spent in the BLAS routines, which you hopefully have an
 * optimized version of. Gromacs includes reference implementations
 * of both BLAS and LAPACK so it compiles everywhere, but you should
 * really try to find a vendor or otherwise optimized version at least
 * of BLAS for better performance.
 *
 * Do NOT use this code for other purposes - we only provide this as a
 * simple fallback/reference implementation when no optimized BLAS
 * is present. If you need an implementation for your own code
 * there are several much faster versions out there.
 *
 * All routines are compatible with the LAPACK/BLAS reference implementations,
 * meaning they assume fortran-style matrix row/column organization.
 *
 * There is plenty of documentation for these routines available
 * at http://www.netlib.org/lapack , so there is no point in repeating
 * it here.
 */
#ifndef GMX_LAPACK_H
#define GMX_LAPACK_H

/*! \cond */

#include "config.h"

/* These are not required by this file, but by the internal LAPACK
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
/* Double precision */

void
    F77_FUNC(dbdsdc, DBDSDC) (const char *uplo, const char *compq, int *n, double *d, double *e, double *u,
                              int *ldu, double *vt, int *ldvt, double *q, int *iq, double *work,
                              int *iwork, int *info);

void
    F77_FUNC(dgetf2, DGETF2) (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

void
    F77_FUNC(dlamrg, DLAMRG) (int *n1, int *n2, double *a, int *dtrd1, int *dtrd2, int *index);

void
    F77_FUNC(dlarnv, DLARNV) (int *idist, int *iseed, int *n, double *x);

void
    F77_FUNC(dlasd0, DLASD0) (int *n, int *sqre, double *d, double *e, double *u,
                              int *ldu, double *vt, int *ldvt, int *smlsiz, int *iwork,
                              double *work, int *info);

void
    F77_FUNC(dlasda, DLASDA) (int *icompq, int *smlsiz, int *n, int *sqre, double *d, double *e,
                              double *u, int *ldu, double *vt, int *k, double *difl, double *difr,
                              double *z, double *poles, int *givptr, int *givcol, int *ldgcol,
                              int *perm, double *givnum, double *c, double *s,
                              double *work, int *iwork, int *info);

void
    F77_FUNC(dlasq6, DLASQ6) (int *i0, int *n0, double *z, int *pp, double *dmin, double *dmin1,
                              double *dmin2, double *dn, double *dnm1, double *dnm2);

void
    F77_FUNC(dorgl2, DORGL2) (int *m, int *n, int *k, double *a, int *lda,
                              double *tau, double *work, int *info);

void
    F77_FUNC(dbdsqr, DBDSQR) (const char *uplo, int *n, int *ncvt, int *nru, int *ncc, double *d,
                              double *e, double *vt, int *ldvt, double *u, int *ldu,
                              double *c, int *ldc, double *work, int *info);

void
    F77_FUNC(dgetrf, DGETRF) (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

void
    F77_FUNC(dgetri, DGETRI) (int *n, double *a, int *lda, int *ipiv, double *work,
                              int *lwork, int *info);

void
    F77_FUNC(dgetrs, DGETRS) (const char *trans, int *n, int *nrhs,   double *a, int *lda, int *ipiv,
                              double *b, int *ldb, int *info);

void
    F77_FUNC(dtrtri, DTRTRI) (const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

void
    F77_FUNC(dtrti2, DTRTI2) (const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

double
    F77_FUNC(dlange, DLANGE) (const char *norm, int *m, int *n, double *a, int *lda, double *work);

void
    F77_FUNC(dlarrbx, DLARRBX) (int *n, double *d, double *l, double *ld, double *lld, int *ifirst,
                                int *ilast, double *rtol1, double *rtol2, int *offset, double *w,
                                double *wgap, double *werr, double *work, int *iwork, int *info);

void
    F77_FUNC(dlasd1, DLASD1) (int *nl, int *nr, int *sqre, double *d, double *alpha, double *beta,
                              double *u, int *ldu, double *vt, int *ldvt, int *idxq, int *iwork,
                              double *work, int *info);

void
    F77_FUNC(dlasdq, DLASDQ) (const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc,
                              double *d, double *e, double *vt, int *ldvt, double *u, int *ldu,
                              double *c, int *ldc, double *work, int *info);

void
    F77_FUNC(dlasr, DLASR) (const char *side, const char *pivot, const char *direct, int *m, int *n, double *c,
                            double *s, double *a, int *lda);

void
    F77_FUNC(dorglq, DORGLQ) (int *m, int *n, int *k, double *a, int *lda,
                              double *tau, double *work, int *lwork, int *info);

void
    F77_FUNC(dormtr, DORMTR) (const char *side, const char *uplo, const char *trans, int *m, int *n, double *a,
                              int *lda, double *tau, double *c, int *ldc,
                              double *work, int *lwork, int *info);

void
    F77_FUNC(dgebd2, DGEBD2) (int *m, int *n, double *a, int *lda, double *d, double *e,
                              double *tauq, double *taup, double *work, int *info);

void
    F77_FUNC(dlabrd, DLABRD) (int *m, int *n, int *nb, double *a, int *lda, double *d,
                              double *e, double *tauq, double *taup, double *x,
                              int *ldx, double *y, int *ldy);

double
    F77_FUNC(dlanst, DLANST) (const char *norm, int *n, double *d, double *e);

double
    F77_FUNC(dlansy, DLANSY) (const char *norm, const char *uplo, int *n, double *a, int *lda, double *work);

void
    F77_FUNC(dlarrex, DLARREX) (const char *range, int *n, double *vl, double *vu, int *il, int *iu,
                                double *d, double *e, double *tol, int *nsplit,
                                int *isplit, int *m, double *w, int *iblock, int *indexw,
                                double *gersch, double *work, int *iwork, int *info);

void
    F77_FUNC(dlasd2, DLASD2) (int *nl, int *nr, int *sqre, int *k, double *d, double *z,
                              double *alpha, double *beta, double *u, int *ldu, double *vt,
                              int *ldvt, double *dsigma, double *u2, int *ldu2, double *vt2,
                              int *ldvt2, int *idxp, int *idx, int *idxc,
                              int *idxq, int *coltyp, int *info);

void
    F77_FUNC(dlasdt, DLASDT) (int *n, int *lvl, int *nd, int *inode, int *ndiml,
                              int *ndimr, int *msub);

void
    F77_FUNC(dlasrt, DLASRT) (const char *id, int *n, double *d, int *info);

void
    F77_FUNC(dlasrt2, DLASRT2) (const char *id, int *n, double *d, int *key, int *info);

void
    F77_FUNC(ilasrt2, ILASRT2) (const char *id, int *n, int *d, int *key, int *info);

void
    F77_FUNC(dorgqr, DORGQR) (int *m, int *n, int *k, double *a, int *lda, double *tau,
                              double *work, int *lwork, int *info);

void
    F77_FUNC(dstebz, DSTEBZ) (const char *range, const char *order, int *n, double *vl, double *vu,
                              int *il, int *iu, double *abstol, double *d, double *e,
                              int *m, int *nsplit, double *w, int *iblock, int *isplit,
                              double *work, int *iwork, int *info);

void
    F77_FUNC(dsteqr, DSTEQR) (const char *compz, int *n, double *d__, double *e,
                              double *z__,  int *ldz, double *work, int *info);

void
    F77_FUNC(dgebrd, DGEBRD) (int *m, int *n, double *a, int *lda, double *d, double *e,
                              double *tauq, double *taup, double *work, int *lwork, int *info);

void
    F77_FUNC(dlacpy, DLACPY) (const char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);

double
    F77_FUNC(dlapy2, DLAPY2) (double * x, double * y);


void
    F77_FUNC(dlarrfx, DLARRFX) (int *n, double *d, double *l, double *ld, double *lld, int *ifirst,
                                int *ilast, double *w, double *sigma, double *dplus, double *lplus,
                                double *work, int *info);

void
    F77_FUNC(dlasd3, DLASD3) (int *nl, int *nr, int *sqre, int *k, double *d, double *q, int *ldq,
                              double *dsigma, double *u, int *ldu, double *u2, int *ldu2,
                              double *vt, int *ldvt, double *vt2, int *ldvt2, int *idxc,
                              int *ctot, double *z, int *info);

void
    F77_FUNC(dlaset, DLASET) (const char *uplo, int *m, int *n, double *alpha,
                              double *beta, double *a, int *lda);

void
    F77_FUNC(dlassq, DLASSQ) (int *n, double *x, int *incx, double *scale, double *sumsq);

void
    F77_FUNC(dorm2l, DORM2L) (const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda,
                              double *tau, double *c, int *ldc, double *work, int *info);

void
    F77_FUNC(dstegr, DSTEGR) (const char *jobz, const char *range, int *n, double *d, double *e, double *vl,
                              double *vu, int *il, int *iu, double *abstol, int *m, double *w,
                              double *z, int *ldz, int *isuppz, double *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(ssteqr, SSTEQR) (const char *compz, int *n, float *d__, float *e,
                              float *z__,  int *ldz, float *work, int *info);

void
    F77_FUNC(dgelq2, DGELQ2) (int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

void
    F77_FUNC(dlae2, DLAE2) (double *a, double *b, double *c, double *rt1, double *rt2);

void
    F77_FUNC(dlaev2, DLAEV2) (double *a, double *b, double *c, double *rt1, double *rt2,
                              double *cs1, double *cs2);

void
    F77_FUNC(dlar1vx, DLAR1VX) (int *n, int *b1, int *bn, double *sigma, double *d, double *l, double *ld,
                                double *lld, double *eval, double *gersch, double *z, double *ztz, double *mingma,
                                int *r, int *isuppz, double *work);

void
    F77_FUNC(dlarrvx, DLARRVX) (int *n, double *d, double *l, int *isplit, int *m, double *w,
                                int *iblock, int *indexw, double *gersch, double *tol, double *z, int *ldz,
                                int *isuppz, double *work, int *iwork, int *info);

void
    F77_FUNC(dlasd4, DLASD4) (int *n, int *i, double *d, double *z, double *delta,
                              double *rho, double *sigma, double *work, int *info);

void
    F77_FUNC(dlasq1, DLASQ1) (int *n, double *d, double *e, double *work, int *info);


void
    F77_FUNC(dlasv2, DLASV2) (double *f, double *g, double *h, double *ssmin, double *ssmax,
                              double *snr, double *csr, double *snl, double *csl);

void
    F77_FUNC(dorm2r, DORM2R) (const char *side, const char *trans, int *m, int *n, int *k, double *a,
                              int *lda, double *tau, double *c, int *ldc, double *work, int *info);

void
    F77_FUNC(dstein, DSTEIN) (int *n, double *d, double *e, int *m, double *w, int *iblock, int *isplit,
                              double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

void
    F77_FUNC(dgelqf, DGELQF) (int *m, int *n, double *a, int *lda, double *tau,
                              double *work, int *lwork, int *info);

void
    F77_FUNC(dlaebz, DLAEBZ) (int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin,
                              double *abstol, double *reltol, double *pivmin, double *d, double *e,
                              double *e2, int *nval, double *ab, double *c, int *mout, int *nab,
                              double *work, int *iwork, int *info);

void
    F77_FUNC(dlarf, DLARF) (const char *side, int *m, int *n, double *v, int *incv, double *tau,
                            double *c, int *ldc, double *work);

void
    F77_FUNC(dlartg, DLARTG) (double *f, double *g, double *cs, double *sn, double *r);

void
    F77_FUNC(dlasd5, DLASD5) (int *i, double *d, double *z, double *delta,
                              double *rho, double *dsigma, double *work);

void
    F77_FUNC(dlasq2, DLASQ2) (int *n, double *z, int *info);

void
    F77_FUNC(dlasq3, DLASQ3) (int *i0, int *n0, double *z, int *pp, double *dmin,
                              double *sigma, double *desig, double *qmax, int *nfail,
                              int *iter, int *ndiv, int *ieee);

void
    F77_FUNC(dlaswp, DLASWP) (int *n, double *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

void
    F77_FUNC(dormbr, DORMBR) (const char *vect, const char *side, const char *trans, int *m, int *n, int *k,
                              double *a, int *lda, double *tau, double *c, int *ldc, double *work,
                              int *lwork, int *info);

void
    F77_FUNC(dsterf, DSTERF) (int *n, double *d, double *e, int *info);

void
    F77_FUNC(dgeqr2, DGEQR2) (int *m, int *n, double *a, int *lda, double *tau,
                              double *work, int *info);

void
    F77_FUNC(dlaed6, DLAED6) (int *kniter, int *orgati, double *rho, double *d,
                              double *z, double *finit, double *tau, int *info);

void
    F77_FUNC(dlarfb, DLARFB) (const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n,
                              int *k, double *v, int *ldv, double *t, int *ldt, double *c,
                              int *ldc, double *work, int *ldwork);

void
    F77_FUNC(dlaruv, DLARUV) (int *iseed, int *n, double *x);

void
    F77_FUNC(dlasd6, DLASD6) (int *icompq, int *nl, int *nr, int *sqre, double *d, double *vf,
                              double *vl, double *alpha, double *beta, int *idxq, int *perm,
                              int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum,
                              double *poles, double *difl, double *difr, double *z, int *k,
                              double *c, double *s, double *work, int *iwork, int *info);

void
    F77_FUNC(dlatrd, DLATRD) (const char *uplo, int *n, int *nb, double *a, int *lda, double *e,
                              double * tau, double *w, int *ldw);

void
    F77_FUNC(dorml2, DORML2) (const char *side, const char *trans, int *m, int *n, int *k, double *a,
                              int *lda, double *tau, double *c, int *ldc, double *work, int *info);

void
    F77_FUNC(dstevr, DSTEVR) (const char *jobz, const char *range, int *n, double *d, double *e, double *vl,
                              double *vu, int *il, int *iu, double *abstol, int *m, double *w,
                              double *z, int *ldz, int *isuppz, double *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(dsytrd, DSYTRD) (const char *uplo, int *n, double *  a, int *lda, double *d,
                              double *e, double *tau, double *work, int *lwork, int *info);

void
    F77_FUNC(dsyevr, DSYEVR) (const char *jobz, const char *range, const char *uplo, int *n,
                              double *a, int *lda, double *vl, double *vu, int *
                              il, int *iu, double *abstol, int *m, double *w,
                              double *z__, int *ldz, int *isuppz, double *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(dormql, DORMQL) (const char *side, const char *trans, int *m, int *n,
                              int *k, double *a, int *lda, double *tau, double *
                              c, int *ldc, double *work, int *lwork, int *info);

void
    F77_FUNC(dormqr, DORMQR) (const char *side, const char *trans, int *m, int *n, int *k, double *a,
                              int *lda, double *tau, double *c, int *ldc,
                              double *work, int *lwork, int *info);

void
    F77_FUNC(dorgbr, DORGBR) (const char *vect, int *m, int *n, int *k, double *a, int *lda,
                              double *tau, double *work, int *lwork, int *info);

void
    F77_FUNC(dlasq5, DLASQ5) (int *i0, int *n0, double *z, int *pp, double *tau, double *dmin,
                              double *dmin1, double *dmin2, double *dn, double *dnm1,
                              double *dnm2, int *ieee);

void
    F77_FUNC(dlasd8, DLASD8) (int *icompq, int *k, double *d, double *z, double *vf, double *vl,
                              double *difl, double *difr, int *lddifr, double *dsigma,
                              double *work, int *info);

void
    F77_FUNC(dlascl, DLASCL) (const char *type, int *kl, int *ku, double *cfrom, double *cto, int *m,
                              int *n, double *a, int *lda, int *info);

void
    F77_FUNC(dlarft, DLARFT) (const char *direct, const char *storev, int *n, int *k, double *v,
                              int *ldv, double *tau, double *t, int *ldt);

void
    F77_FUNC(dlagts, DLAGTS) (int *job, int *n, double *a, double *b, double *c, double *d,
                              int *in, double *y, double *tol, int *info);

void
    F77_FUNC(dgesdd, DGESDD) (const char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u,
                              int *ldu, double *vt, int *ldvt, double *work, int *lwork,
                              int *iwork, int *info);

void
    F77_FUNC(dsytd2, DSYTD2) (const char *uplo, int *n, double *a, int *lda, double *d,
                              double *e, double *tau, int *info);

void
    F77_FUNC(dormlq, DORMLQ) (const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda,
                              double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

void
    F77_FUNC(dorg2r, DORG2R) (int *m, int *n, int *k, double *a, int *lda, double *tau,
                              double *work, int *info);

void
    F77_FUNC(dlasq4, DLASQ4) (int *i0, int *n0, double *z, int *pp, int *n0in, double *dmin,
                              double *dmin1, double *dmin2, double *dn, double *dn1,
                              double *dn2, double *tau, int *ttype);

void
    F77_FUNC(dlasd7, DLASD7) (int *icompq, int *nl, int *nr, int *sqre, int *k, double *d, double *z,
                              double *zw, double *vf, double *vfw, double *vl, double *vlw,
                              double *alpha, double *beta, double *dsigma, int *idx, int *idxp,
                              int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol,
                              double *givnum, int *ldgnum, double *c, double *s, int *info);

void
    F77_FUNC(dlas2, DLAS2) (double *f, double *g, double *h, double *ssmin, double *ssmax);

void
    F77_FUNC(dlarfg, DLARFG) (int *n, double *alpha, double *x, int *incx, double *tau);

void
    F77_FUNC(dlagtf, DLAGTF) (int *n, double *a, double *lambda, double *b, double *c,
                              double *tol, double *d, int *in, int *info);

void
    F77_FUNC(dgeqrf, DGEQRF) (int *m, int *n, double *a, int *lda, double *tau,
                              double *work, int *lwork, int *info);



/* Single precision */

void
    F77_FUNC(sbdsdc, SBDSDC) (const char *uplo, const char *compq, int *n, float *d, float *e, float *u,
                              int *ldu, float *vt, int *ldvt, float *q, int *iq, float *work,
                              int *iwork, int *info);

void
    F77_FUNC(sgetf2, SGETF2) (int *m, int *n, float *a, int *lda, int *ipiv, int *info);

void
    F77_FUNC(slamrg, SLAMRG) (int *n1, int *n2, float *a, int *dtrd1, int *dtrd2, int *index);

void
    F77_FUNC(slarnv, SLARNV) (int *idist, int *iseed, int *n, float *x);

void
    F77_FUNC(slasd0, SLASD0) (int *n, int *sqre, float *d, float *e, float *u,
                              int *ldu, float *vt, int *ldvt, int *smlsiz, int *iwork,
                              float *work, int *info);

void
    F77_FUNC(slasda, SLASDA) (int *icompq, int *smlsiz, int *n, int *sqre, float *d, float *e,
                              float *u, int *ldu, float *vt, int *k, float *difl, float *difr,
                              float *z, float *poles, int *givptr, int *givcol, int *ldgcol,
                              int *perm, float *givnum, float *c, float *s,
                              float *work, int *iwork, int *info);

void
    F77_FUNC(slasq6, SLASQ6) (int *i0, int *n0, float *z, int *pp, float *dmin, float *dmin1,
                              float *dmin2, float *dn, float *dnm1, float *dnm2);

void
    F77_FUNC(sorgl2, SORGL2) (int *m, int *n, int *k, float *a, int *lda,
                              float *tau, float *work, int *info);

void
    F77_FUNC(sbdsqr, SBDSQR) (const char *uplo, int *n, int *ncvt, int *nru, int *ncc, float *d,
                              float *e, float *vt, int *ldvt, float *u, int *ldu,
                              float *c, int *ldc, float *work, int *info);

void
    F77_FUNC(sgetrf, SGETRF) (int *m, int *n, float *a, int *lda, int *ipiv, int *info);

void
    F77_FUNC(sgetri, SGETRI) (int *n, float *a, int *lda, int *ipiv, float *work,
                              int *lwork, int *info);

void
    F77_FUNC(sgetrs, SGETRS) (const char *trans, int *n, int *nrhs,   float *a, int *lda, int *ipiv,
                              float *b, int *ldb, int *info);

void
    F77_FUNC(strtri, STRTRI) (const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

void
    F77_FUNC(strti2, STRTI2) (const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

float
    F77_FUNC(slange, SLANGE) (const char *norm, int *m, int *n, float *a, int *lda, float *work);

void
    F77_FUNC(slarrbx, SLARRBX) (int *n, float *d, float *l, float *ld, float *lld, int *ifirst,
                                int *ilast, float *rtol1, float *rtol2, int *offset, float *w,
                                float *wgap, float *werr, float *work, int *iwork, int *info);

void
    F77_FUNC(slasd1, SLASD1) (int *nl, int *nr, int *sqre, float *d, float *alpha, float *beta,
                              float *u, int *ldu, float *vt, int *ldvt, int *idxq, int *iwork,
                              float *work, int *info);

void
    F77_FUNC(slasdq, SLASDQ) (const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc,
                              float *d, float *e, float *vt, int *ldvt, float *u, int *ldu,
                              float *c, int *ldc, float *work, int *info);

void
    F77_FUNC(slasr, SLASR) (const char *side, const char *pivot, const char *direct, int *m, int *n, float *c,
                            float *s, float *a, int *lda);

void
    F77_FUNC(sorglq, SORGLQ) (int *m, int *n, int *k, float *a, int *lda,
                              float *tau, float *work, int *lwork, int *info);

void
    F77_FUNC(sormtr, SORMTR) (const char *side, const char *uplo, const char *trans, int *m, int *n, float *a,
                              int *lda, float *tau, float *c, int *ldc,
                              float *work, int *lwork, int *info);

void
    F77_FUNC(sgebd2, SGEBD2) (int *m, int *n, float *a, int *lda, float *d, float *e,
                              float *tauq, float *taup, float *work, int *info);

void
    F77_FUNC(slabrd, SLABRD) (int *m, int *n, int *nb, float *a, int *lda, float *d,
                              float *e, float *tauq, float *taup, float *x,
                              int *ldx, float *y, int *ldy);

float
    F77_FUNC(slanst, SLANST) (const char *norm, int *n, float *d, float *e);

float
    F77_FUNC(slansy, SLANSY) (const char *norm, const char *uplo, int *n, float *a, int *lda, float *work);

void
    F77_FUNC(slarrex, SLARREX) (const char *range, int *n, float *vl, float *vu, int *il, int *iu,
                                float *d, float *e, float *tol, int *nsplit,
                                int *isplit, int *m, float *w, int *iblock, int *indexw,
                                float *gersch, float *work, int *iwork, int *info);

void
    F77_FUNC(slasd2, SLASD2) (int *nl, int *nr, int *sqre, int *k, float *d, float *z,
                              float *alpha, float *beta, float *u, int *ldu, float *vt,
                              int *ldvt, float *dsigma, float *u2, int *ldu2, float *vt2,
                              int *ldvt2, int *idxp, int *idx, int *idxc,
                              int *idxq, int *coltyp, int *info);

void
    F77_FUNC(slasdt, SLASDT) (int *n, int *lvl, int *nd, int *inode, int *ndiml,
                              int *ndimr, int *msub);

void
    F77_FUNC(slasrt, SLASRT) (const char *id, int *n, float *d, int *info);

void
    F77_FUNC(slasrt2, SLASRT2) (const char *id, int *n, float *d, int *key, int *info);

void
    F77_FUNC(sorgqr, SORGQR) (int *m, int *n, int *k, float *a, int *lda, float *tau,
                              float *work, int *lwork, int *info);

void
    F77_FUNC(sstebz, SSTEBZ) (const char *range, const char *order, int *n, float *vl, float *vu,
                              int *il, int *iu, float *abstol, float *d, float *e,
                              int *m, int *nsplit, float *w, int *iblock, int *isplit,
                              float *work, int *iwork, int *info);

void
    F77_FUNC(sgebrd, SGEBRD) (int *m, int *n, float *a, int *lda, float *d, float *e,
                              float *tauq, float *taup, float *work, int *lwork, int *info);

void
    F77_FUNC(slacpy, SLACPY) (const char *uplo, int *m, int *n, float *a, int *lda, float *b, int *ldb);

float
    F77_FUNC(slapy2, SLAPY2) (float * x, float * y);

void
    F77_FUNC(slarrfx, SLARRFX) (int *n, float *d, float *l, float *ld, float *lld, int *ifirst,
                                int *ilast, float *w, float *sigma, float *dplus, float *lplus,
                                float *work, int *info);

void
    F77_FUNC(slasd3, SLASD3) (int *nl, int *nr, int *sqre, int *k, float *d, float *q, int *ldq,
                              float *dsigma, float *u, int *ldu, float *u2, int *ldu2,
                              float *vt, int *ldvt, float *vt2, int *ldvt2, int *idxc,
                              int *ctot, float *z, int *info);

void
    F77_FUNC(slaset, SLASET) (const char *uplo, int *m, int *n, float *alpha,
                              float *beta, float *a, int *lda);

void
    F77_FUNC(slassq, SLASSQ) (int *n, float *x, int *incx, float *scale, float *sumsq);

void
    F77_FUNC(sorm2l, SORM2L) (const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda,
                              float *tau, float *c, int *ldc, float *work, int *info);

void
    F77_FUNC(sstegr, SSTEGR) (const char *jobz, const char *range, int *n, float *d, float *e, float *vl,
                              float *vu, int *il, int *iu, float *abstol, int *m, float *w,
                              float *z, int *ldz, int *isuppz, float *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(sgelq2, SGELQ2) (int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

void
    F77_FUNC(slae2, SLAE2) (float *a, float *b, float *c, float *rt1, float *rt2);

void
    F77_FUNC(slaev2, SLAEV2) (float *a, float *b, float *c, float *rt1, float *rt2,
                              float *cs1, float *cs2);

void
    F77_FUNC(slar1vx, SLAR1VX) (int *n, int *b1, int *bn, float *sigma, float *d, float *l, float *ld,
                                float *lld, float *eval, float *gersch, float *z, float *ztz, float *mingma,
                                int *r, int *isuppz, float *work);

void
    F77_FUNC(slarrvx, SLARRVX) (int *n, float *d, float *l, int *isplit, int *m, float *w,
                                int *iblock, int *indexw, float *gersch, float *tol, float *z, int *ldz,
                                int *isuppz, float *work, int *iwork, int *info);

void
    F77_FUNC(slasd4, SLASD4) (int *n, int *i, float *d, float *z, float *delta,
                              float *rho, float *sigma, float *work, int *info);

void
    F77_FUNC(slasq1, SLASQ1) (int *n, float *d, float *e, float *work, int *info);


void
    F77_FUNC(slasv2, SLASV2) (float *f, float *g, float *h, float *ssmin, float *ssmax,
                              float *snr, float *csr, float *snl, float *csl);

void
    F77_FUNC(sorm2r, SORM2R) (const char *side, const char *trans, int *m, int *n, int *k, float *a,
                              int *lda, float *tau, float *c, int *ldc, float *work, int *info);

void
    F77_FUNC(sstein, SSTEIN) (int *n, float *d, float *e, int *m, float *w, int *iblock, int *isplit,
                              float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

void
    F77_FUNC(sgelqf, SGELQF) (int *m, int *n, float *a, int *lda, float *tau,
                              float *work, int *lwork, int *info);

void
    F77_FUNC(slaebz, SLAEBZ) (int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin,
                              float *abstol, float *reltol, float *pivmin, float *d, float *e,
                              float *e2, int *nval, float *ab, float *c, int *mout, int *nab,
                              float *work, int *iwork, int *info);

void
    F77_FUNC(slarf, SLARF) (const char *side, int *m, int *n, float *v, int *incv, float *tau,
                            float *c, int *ldc, float *work);

void
    F77_FUNC(slartg, SLARTG) (float *f, float *g, float *cs, float *sn, float *r);

void
    F77_FUNC(slasd5, SLASD5) (int *i, float *d, float *z, float *delta,
                              float *rho, float *dsigma, float *work);

void
    F77_FUNC(slasq2, SLASQ2) (int *n, float *z, int *info);

void
    F77_FUNC(slasq3, SLASQ3) (int *i0, int *n0, float *z, int *pp, float *dmin,
                              float *sigma, float *desig, float *qmax, int *nfail,
                              int *iter, int *ndiv, int *ieee);

void
    F77_FUNC(slaswp, SLASWP) (int *n, float *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

void
    F77_FUNC(sormbr, SORMBR) (const char *vect, const char *side, const char *trans, int *m, int *n, int *k,
                              float *a, int *lda, float *tau, float *c, int *ldc, float *work,
                              int *lwork, int *info);

void
    F77_FUNC(ssterf, SSTERF) (int *n, float *d, float *e, int *info);

void
    F77_FUNC(sgeqr2, SGEQR2) (int *m, int *n, float *a, int *lda, float *tau,
                              float *work, int *info);

void
    F77_FUNC(slaed6, SLAED6) (int *kniter, int *orgati, float *rho, float *d,
                              float *z, float *finit, float *tau, int *info);

void
    F77_FUNC(slarfb, SLARFB) (const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n,
                              int *k, float *v, int *ldv, float *t, int *ldt, float *c,
                              int *ldc, float *work, int *ldwork);

void
    F77_FUNC(slaruv, SLARUV) (int *iseed, int *n, float *x);

void
    F77_FUNC(slasd6, SLASD6) (int *icompq, int *nl, int *nr, int *sqre, float *d, float *vf,
                              float *vl, float *alpha, float *beta, int *idxq, int *perm,
                              int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum,
                              float *poles, float *difl, float *difr, float *z, int *k,
                              float *c, float *s, float *work, int *iwork, int *info);

void
    F77_FUNC(slatrd, SLATRD) (const char *uplo, int *n, int *nb, float *a, int *lda, float *e,
                              float * tau, float *w, int *ldw);

void
    F77_FUNC(sorml2, SORML2) (const char *side, const char *trans, int *m, int *n, int *k, float *a,
                              int *lda, float *tau, float *c, int *ldc, float *work, int *info);

void
    F77_FUNC(sstevr, SSTEVR) (const char *jobz, const char *range, int *n, float *d, float *e, float *vl,
                              float *vu, int *il, int *iu, float *abstol, int *m, float *w,
                              float *z, int *ldz, int *isuppz, float *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(ssytrd, SSYTRD) (const char *uplo, int *n, float *  a, int *lda, float *d,
                              float *e, float *tau, float *work, int *lwork, int *info);

void
    F77_FUNC(ssyevr, SSYEVR) (const char *jobz, const char *range, const char *uplo, int *n,
                              float *a, int *lda, float *vl, float *vu, int *
                              il, int *iu, float *abstol, int *m, float *w,
                              float *z__, int *ldz, int *isuppz, float *work,
                              int *lwork, int *iwork, int *liwork, int *info);

void
    F77_FUNC(sormql, SORMQL) (const char *side, const char *trans, int *m, int *n,
                              int *k, float *a, int *lda, float *tau, float *
                              c, int *ldc, float *work, int *lwork, int *info);

void
    F77_FUNC(sormqr, SORMQR) (const char *side, const char *trans, int *m, int *n, int *k, float *a,
                              int *lda, float *tau, float *c, int *ldc,
                              float *work, int *lwork, int *info);

void
    F77_FUNC(sorgbr, SORGBR) (const char *vect, int *m, int *n, int *k, float *a, int *lda,
                              float *tau, float *work, int *lwork, int *info);

void
    F77_FUNC(slasq5, SLASQ5) (int *i0, int *n0, float *z, int *pp, float *tau, float *dmin,
                              float *dmin1, float *dmin2, float *dn, float *dnm1,
                              float *dnm2, int *ieee);

void
    F77_FUNC(slasd8, SLASD8) (int *icompq, int *k, float *d, float *z, float *vf, float *vl,
                              float *difl, float *difr, int *lddifr, float *dsigma,
                              float *work, int *info);

void
    F77_FUNC(slascl, SLASCL) (const char *type, int *kl, int *ku, float *cfrom, float *cto, int *m,
                              int *n, float *a, int *lda, int *info);

void
    F77_FUNC(slarft, SLARFT) (const char *direct, const char *storev, int *n, int *k, float *v,
                              int *ldv, float *tau, float *t, int *ldt);

void
    F77_FUNC(slagts, SLAGTS) (int *job, int *n, float *a, float *b, float *c, float *d,
                              int *in, float *y, float *tol, int *info);

void
    F77_FUNC(sgesdd, SGESDD) (const char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u,
                              int *ldu, float *vt, int *ldvt, float *work, int *lwork,
                              int *iwork, int *info);

void
    F77_FUNC(ssytd2, SSYTD2) (const char *uplo, int *n, float *a, int *lda, float *d,
                              float *e, float *tau, int *info);

void
    F77_FUNC(sormlq, SORMLQ) (const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda,
                              float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

void
    F77_FUNC(sorg2r, SORG2R) (int *m, int *n, int *k, float *a, int *lda, float *tau,
                              float *work, int *info);

void
    F77_FUNC(slasq4, SLASQ4) (int *i0, int *n0, float *z, int *pp, int *n0in, float *dmin,
                              float *dmin1, float *dmin2, float *dn, float *dn1,
                              float *dn2, float *tau, int *ttype);

void
    F77_FUNC(slasd7, SLASD7) (int *icompq, int *nl, int *nr, int *sqre, int *k, float *d, float *z,
                              float *zw, float *vf, float *vfw, float *vl, float *vlw,
                              float *alpha, float *beta, float *dsigma, int *idx, int *idxp,
                              int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol,
                              float *givnum, int *ldgnum, float *c, float *s, int *info);

void
    F77_FUNC(slas2, SLAS2) (float *f, float *g, float *h, float *ssmin, float *ssmax);

void
    F77_FUNC(slarfg, SLARFG) (int *n, float *alpha, float *x, int *incx, float *tau);

void
    F77_FUNC(slagtf, SLAGTF) (int *n, float *a, float *lambda, float *b, float *c,
                              float *tol, float *d, int *in, int *info);

void
    F77_FUNC(sgeqrf, SGEQRF) (int *m, int *n, float *a, int *lda, float *tau,
                              float *work, int *lwork, int *info);


#ifdef __cplusplus
}
#endif

/*! \endcond */

#endif /* GMX_LAPACK_H */
