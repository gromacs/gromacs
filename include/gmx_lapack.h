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

#ifndef _GMX_LAPACK_H_
#define _GMX_LAPACK_H_

#include "visibility.h"
/** @file
 *
 *  @brief Header definitions for the standard LAPACK library.
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



#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#include "types/simple.h"

/* Suppress Cygwin compiler warnings from using newlib version of
 * ctype.h */
#ifdef GMX_CYGWIN
#undef toupper
#endif

#ifndef FortranCInterface_GLOBAL
#define FortranCInterface_GLOBAL(name,NAME) name ## _
#endif



/* Double precision */

void
FortranCInterface_GLOBAL(dbdsdc,DBDSDC)(const char *uplo, const char *compq, int *n, double *d, double *e, double *u, 
	int *ldu, double *vt, int *ldvt, double *q, int *iq, double *work, 
	int *iwork, int *info);

void
FortranCInterface_GLOBAL(dgetf2,DGETF2)(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

void
FortranCInterface_GLOBAL(dlamrg,DLAMRG)(int *n1, int *n2, double *a, int *dtrd1, int *dtrd2, int *index);

void
FortranCInterface_GLOBAL(dlarnv,DLARNV)(int *idist, int *iseed, int *n, double *x);

void 
FortranCInterface_GLOBAL(dlasd0,DLASD0)(int *n, int *sqre, double *d, double *e, double *u, 
	int *ldu, double *vt, int *ldvt, int *smlsiz, int *iwork, 
	double *work, int *info);

void 
FortranCInterface_GLOBAL(dlasda,DLASDA)(int *icompq, int *smlsiz, int *n, int *sqre, double *d, double *e, 
	double *u, int *ldu, double *vt, int *k, double *difl, double *difr, 
	double *z, double *poles, int *givptr, int *givcol, int *ldgcol, 
	int *perm, double *givnum, double *c, double *s, 
	double *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(dlasq6,DLASQ6)(int *i0, int *n0, double *z, int *pp, double *dmin, double *dmin1, 
	double *dmin2, double *dn, double *dnm1, double *dnm2);

void
FortranCInterface_GLOBAL(dorgl2,DORGL2)(int *m,	int *n,	int *k,	double *a, int *lda, 
	double *tau, double *work, int *info);

void
FortranCInterface_GLOBAL(dbdsqr,DBDSQR)(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, double *d, 
	double *e, double *vt, int *ldvt, double *u, int *ldu,
	double *c, int *ldc, double *work, int *info);

void
FortranCInterface_GLOBAL(dgetrf,DGETRF)(int *m,	int *n,	double *a, int *lda, int *ipiv, int *info);

void
FortranCInterface_GLOBAL(dgetri,DGETRI)(int *n,	double *a, int *lda, int *ipiv, double *work, 
	int *lwork, int *info);

void
FortranCInterface_GLOBAL(dgetrs,DGETRS)(const char *trans, int *n, int *nrhs,	double *a, int *lda, int *ipiv,
	double *b, int *ldb, int *info);

void
FortranCInterface_GLOBAL(dtrtri,DTRTRI)(const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

void
FortranCInterface_GLOBAL(dtrti2,DTRTI2)(const char *uplo, const char *diag, int *n, double *a, int *lda, int *info);

double
FortranCInterface_GLOBAL(dlange,DLANGE)(const char *norm, int *m, int *n, double *a, int *lda, double *work);

void
FortranCInterface_GLOBAL(dlarrbx,DLARRBX)(int *n, double *d, double *l, double *ld, double *lld, int *ifirst,
	 int *ilast, double *rtol1, double *rtol2, int *offset, double *w,
	 double *wgap, double *werr, double *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(dlasd1,DLASD1)(int *nl, int *nr, int *sqre, double *d, double *alpha, double *beta, 
	double *u, int *ldu, double *vt, int *ldvt, int *idxq, int *iwork, 
	double *work, int *info);

void
FortranCInterface_GLOBAL(dlasdq,DLASDQ)(const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc,
	double *d, double *e, double *vt, int *ldvt, double *u, int *ldu, 
	double *c, int *ldc, double *work, int *info);

void 
FortranCInterface_GLOBAL(dlasr,DLASR)(const char *side, const char *pivot, const char *direct, int *m, int *n, double *c, 
       double *s, double *a, int *lda);

void 
FortranCInterface_GLOBAL(dorglq,DORGLQ)(int *m, int *n, int *k, double *a, int *lda, 
	double *tau, double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dormtr,DORMTR)(const char *side, const char *uplo, const char *trans, int *m, int *n, double *a, 
	int *lda, double *tau, double *c, int *ldc,
	double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dgebd2,DGEBD2)(int *m, int *n, double *a, int *lda, double *d, double *e,
	double *tauq, double *taup, double *work, int *info);

void 
FortranCInterface_GLOBAL(dlabrd,DLABRD)(int *m, int *n, int *nb, double *a, int *lda, double *d,
	double *e, double *tauq, double *taup, double *x,
	int *ldx, double *y, int *ldy);

double
FortranCInterface_GLOBAL(dlanst,DLANST)(const char *norm, int *n, double *d, double *e);

double
FortranCInterface_GLOBAL(dlansy,DLANSY)(const char *norm, const char *uplo, int *n, double *a, int *lda, double *work);

void
FortranCInterface_GLOBAL(dlarrex,DLARREX)(const char *range, int *n, double *vl, double *vu, int *il, int *iu,
	 double *d, double *e, double *tol, int *nsplit, 
	 int *isplit, int *m, double *w, int *iblock, int *indexw,
	 double *gersch, double *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(dlasd2,DLASD2)(int *nl, int *nr, int *sqre, int *k, double *d, double *z, 
	double *alpha, double *beta, double *u, int *ldu, double *vt, 
	int *ldvt, double *dsigma, double *u2, int *ldu2, double *vt2, 
	int *ldvt2, int *idxp, int *idx, int *idxc, 
	int *idxq, int *coltyp, int *info);

void
FortranCInterface_GLOBAL(dlasdt,DLASDT)(int *n, int *lvl, int *nd, int *inode, int *ndiml, 
	int *ndimr, int *msub);

void 
FortranCInterface_GLOBAL(dlasrt,DLASRT)(const char *id, int *n, double *d, int *info);

void
FortranCInterface_GLOBAL(dlasrt2,DLASRT2)(const char *id, int *n, double *d, int *key, int *info);

void
FortranCInterface_GLOBAL(ilasrt2,ILASRT2)(const char *id, int *n, int *d, int *key, int *info);

void 
FortranCInterface_GLOBAL(dorgqr,DORGQR)(int *m, int *n, int *k, double *a, int *lda, double *tau, 
	double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dstebz,DSTEBZ)(const char *range, const char *order, int *n, double *vl, double *vu, 
	int *il, int *iu, double *abstol, double *d, double *e, 
	int *m, int *nsplit, double *w, int *iblock, int *isplit, 
	double *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(dsteqr,DSTEQR)(const char *compz, int *n, double *d__, double *e, 
        double *z__,  int *ldz, double *work, int *info);

void
FortranCInterface_GLOBAL(dgebrd,DGEBRD)(int *m, int *n, double *a, int *lda, double *d, double *e,
	double *tauq, double *taup, double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dlacpy,DLACPY)(const char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb);

double
FortranCInterface_GLOBAL(dlapy2,DLAPY2)(double * x, double * y);


void
FortranCInterface_GLOBAL(dlarrfx,DLARRFX)(int *n, double *d, double *l, double *ld, double *lld, int *ifirst,
        int *ilast, double *w, double *sigma, double *dplus, double *lplus,
        double *work, int *info);

void 
FortranCInterface_GLOBAL(dlasd3,DLASD3)(int *nl, int *nr, int *sqre, int *k, double *d, double *q, int *ldq, 
	double *dsigma, double *u, int *ldu, double *u2, int *ldu2, 
	double *vt, int *ldvt, double *vt2, int *ldvt2, int *idxc, 
	int *ctot, double *z, int *info);

void
FortranCInterface_GLOBAL(dlaset,DLASET)(const char *uplo, int *m, int *n, double *alpha, 
	double *beta, double *a, int *lda);

void
FortranCInterface_GLOBAL(dlassq,DLASSQ)(int *n, double *x, int *incx, double *scale, double *sumsq);

void
FortranCInterface_GLOBAL(dorm2l,DORM2L)(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, 
	double *tau, double *c, int *ldc, double *work, int *info);

void
FortranCInterface_GLOBAL(dstegr,DSTEGR)(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int *m, double *w, 
	double *z, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork,	int *liwork, int *info);

void
FortranCInterface_GLOBAL(ssteqr,SSTEQR)(const char *compz, int *n, float *d__, float *e, 
        float *z__,  int *ldz, float *work, int *info);

void
FortranCInterface_GLOBAL(dgelq2,DGELQ2)(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);

void
FortranCInterface_GLOBAL(dlae2,DLAE2)(double *a, double *b, double *c, double *rt1, double *rt2);

void
FortranCInterface_GLOBAL(dlaev2,DLAEV2)(double *a, double *b, double *c, double *rt1, double *rt2,
	double *cs1, double *cs2);

void
FortranCInterface_GLOBAL(dlar1vx,DLAR1VX)(int *n, int *b1, int *bn, double *sigma, double *d, double *l, double *ld, 
                          double *lld, double *eval, double *gersch, double *z, double *ztz, double *mingma, 
                          int *r, int *isuppz, double *work);

void
FortranCInterface_GLOBAL(dlarrvx,DLARRVX)(int *n, double *d, double *l, int *isplit, int *m, double *w, 
	 int *iblock, int *indexw, double *gersch, double *tol, double *z, int *ldz, 
	 int *isuppz, double *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(dlasd4,DLASD4)(int *n, int *i, double *d, double *z, double *delta, 
	double *rho, double *sigma, double *work, int *info);

void
FortranCInterface_GLOBAL(dlasq1,DLASQ1)(int *n,	double *d, double *e, double *work, int *info);


void 
FortranCInterface_GLOBAL(dlasv2,DLASV2)(double *f, double *g, double *h, double *ssmin, double *ssmax, 
	double *snr, double *csr, double *snl, double *csl);

void 
FortranCInterface_GLOBAL(dorm2r,DORM2R)(const char *side, const char *trans, int *m, int *n, int *k, double *a, 
	int *lda, double *tau, double *c, int *ldc, double *work, int *info);

void
FortranCInterface_GLOBAL(dstein,DSTEIN)(int *n, double *d, double *e, int *m, double *w, int *iblock, int *isplit, 
	double *z, int *ldz, double *work, int *iwork, int *ifail, int *info);

void
FortranCInterface_GLOBAL(dgelqf,DGELQF)(int *m,	int *n, double *a, int *lda, double *tau,
	double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dlaebz,DLAEBZ)(int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin,
	double *abstol, double *reltol, double *pivmin, double *d, double *e,
	double *e2, int *nval, double *ab, double *c, int *mout, int *nab,
	double *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(dlarf,DLARF)(const char *side, int *m, int *n, double *v, int *incv, double *tau,
       double *c, int *ldc, double *work);

void
FortranCInterface_GLOBAL(dlartg,DLARTG)(double *f, double *g, double *cs, double *sn, double *r);

void 
FortranCInterface_GLOBAL(dlasd5,DLASD5)(int *i, double *d, double *z, double *delta, 
	double *rho, double *dsigma, double *work);

void 
FortranCInterface_GLOBAL(dlasq2,DLASQ2)(int *n, double *z, int *info);

void 
FortranCInterface_GLOBAL(dlasq3,DLASQ3)(int *i0, int *n0, double *z, int *pp, double *dmin, 
	double *sigma, double *desig, double *qmax, int *nfail, 
	int *iter, int *ndiv, int *ieee);

void
FortranCInterface_GLOBAL(dlaswp,DLASWP)(int *n,	double *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

void 
FortranCInterface_GLOBAL(dormbr,DORMBR)(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, 
	double *a, int *lda, double *tau, double *c, int *ldc, double *work,
	int *lwork, int *info);

void
FortranCInterface_GLOBAL(dsterf,DSTERF)(int *n, double *d, double *e, int *info);

void
FortranCInterface_GLOBAL(dgeqr2,DGEQR2)(int *m,	int *n,	double *a, int *lda, double *tau, 
	double *work, int *info);

void 
FortranCInterface_GLOBAL(dlaed6,DLAED6)(int *kniter, int *orgati, double *rho, double *d, 
	double *z, double *finit, double *tau, int *info);

void 
FortranCInterface_GLOBAL(dlarfb,DLARFB)(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, 
	int *k, double *v, int *ldv, double *t, int *ldt, double *c,
	int *ldc, double *work, int *ldwork);

void
FortranCInterface_GLOBAL(dlaruv,DLARUV)(int *iseed, int *n, double *x);

void 
FortranCInterface_GLOBAL(dlasd6,DLASD6)(int *icompq, int *nl, int *nr, int *sqre, double *d, double *vf, 
	double *vl, double *alpha, double *beta, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, 
	double *poles, double *difl, double *difr, double *z, int *k, 
	double *c, double *s, double *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(dlatrd,DLATRD)(const char *uplo, int *n, int *nb, double *a, int *lda, double *e, 
	double * tau, double *w, int *ldw);

void
FortranCInterface_GLOBAL(dorml2,DORML2)(const char *side, const char *trans, int *m, int *n, int *k, double *a,
	int *lda, double *tau, double *c, int *ldc, double *work, int *info);

void
FortranCInterface_GLOBAL(dstevr,DSTEVR)(const char *jobz, const char *range, int *n, double *d, double *e, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int *m, double *w, 
	double *z, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);

void
FortranCInterface_GLOBAL(dsytrd,DSYTRD)(const char *uplo, int *n, double *  a, int *lda, double *d, 
	double *e, double *tau, double *work, int *lwork, int *info);

GMX_LIBGMX_EXPORT
void
FortranCInterface_GLOBAL(dsyevr,DSYEVR)(const char *jobz, const char *range, const char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);

void
FortranCInterface_GLOBAL(dormql,DORMQL)(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c, int *ldc, double *work, int *lwork, int *info);

void 
FortranCInterface_GLOBAL(dormqr,DORMQR)(const char *side, const char *trans, int *m, int *n, int *k, double *a, 
        int *lda, double *tau, double *c, int *ldc, 
        double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dorgbr,DORGBR)(const char *vect, int *m, int *n, int *k, double *a, int *lda,
	double *tau, double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dlasq5,DLASQ5)(int *i0, int *n0, double *z, int *pp, double *tau, double *dmin, 
	double *dmin1, double *dmin2, double *dn, double *dnm1, 
	double *dnm2, int *ieee);

void 
FortranCInterface_GLOBAL(dlasd8,DLASD8)(int *icompq, int *k, double *d, double *z, double *vf, double *vl, 
	double *difl, double *difr, int *lddifr, double *dsigma, 
	double *work, int *info);

void
FortranCInterface_GLOBAL(dlascl,DLASCL)(const char *type, int *kl, int *ku, double *cfrom, double *cto, int *m, 
	int *n, double *a, int *lda, int *info);

void 
FortranCInterface_GLOBAL(dlarft,DLARFT)(const char *direct, const char *storev, int *n, int *k, double *v, 
	int *ldv, double *tau, double *t, int *ldt);

void
FortranCInterface_GLOBAL(dlagts,DLAGTS)(int *job, int *n, double *a, double *b, double *c, double *d, 
	int *in, double *y, double *tol, int *info);

void 
FortranCInterface_GLOBAL(dgesdd,DGESDD)(const char *jobz, int *m, int *n, double *a, int *lda, double *s, double *u, 
	int *ldu, double *vt, int *ldvt, double *work, int *lwork, 
	int *iwork, int *info);

void
FortranCInterface_GLOBAL(dsytd2,DSYTD2)(const char *uplo, int *n, double *a, int *lda, double *d, 
	double *e, double *tau, int *info);

void 
FortranCInterface_GLOBAL(dormlq,DORMLQ)(const char *side, const char *trans, int *m, int *n, int *k, double *a, int *lda, 
	double *tau, double *c, int *ldc, double *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(dorg2r,DORG2R)(int *m, int *n,	int *k,	double *a, int *lda, double *tau,
	double *work, int *info);

void 
FortranCInterface_GLOBAL(dlasq4,DLASQ4)(int *i0, int *n0, double *z, int *pp, int *n0in, double *dmin, 
	double *dmin1, double *dmin2, double *dn, double *dn1, 
	double *dn2, double *tau, int *ttype);

void 
FortranCInterface_GLOBAL(dlasd7,DLASD7)(int *icompq, int *nl, int *nr, int *sqre, int *k, double *d, double *z,
	double *zw, double *vf, double *vfw, double *vl, double *vlw,
	double *alpha, double *beta, double *dsigma, int *idx, int *idxp,
	int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, 
	double *givnum, int *ldgnum, double *c, double *s, int *info);

void
FortranCInterface_GLOBAL(dlas2,DLAS2)(double *f, double *g, double *h, double *ssmin, double *ssmax);

void
FortranCInterface_GLOBAL(dlarfg,DLARFG)(int *n, double *alpha, double *x, int *incx, double *tau);

void
FortranCInterface_GLOBAL(dlagtf,DLAGTF)(int *n, double *a, double *lambda, double *b, double *c, 
	double *tol, double *d, int *in, int *info);

void 
FortranCInterface_GLOBAL(dgeqrf,DGEQRF)(int *m, int *n, double *a, int *lda, double *tau,
	double *work, int *lwork, int *info);



/* Single precision */

void
FortranCInterface_GLOBAL(sbdsdc,SBDSDC)(const char *uplo, const char *compq, int *n, float *d, float *e, float *u, 
	int *ldu, float *vt, int *ldvt, float *q, int *iq, float *work, 
	int *iwork, int *info);

void
FortranCInterface_GLOBAL(sgetf2,SGETF2)(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

void
FortranCInterface_GLOBAL(slamrg,SLAMRG)(int *n1, int *n2, float *a, int *dtrd1, int *dtrd2, int *index);

void
FortranCInterface_GLOBAL(slarnv,SLARNV)(int *idist, int *iseed, int *n, float *x);

void 
FortranCInterface_GLOBAL(slasd0,SLASD0)(int *n, int *sqre, float *d, float *e, float *u, 
	int *ldu, float *vt, int *ldvt, int *smlsiz, int *iwork, 
	float *work, int *info);

void 
FortranCInterface_GLOBAL(slasda,SLASDA)(int *icompq, int *smlsiz, int *n, int *sqre, float *d, float *e, 
	float *u, int *ldu, float *vt, int *k, float *difl, float *difr, 
	float *z, float *poles, int *givptr, int *givcol, int *ldgcol, 
	int *perm, float *givnum, float *c, float *s, 
	float *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(slasq6,SLASQ6)(int *i0, int *n0, float *z, int *pp, float *dmin, float *dmin1, 
	float *dmin2, float *dn, float *dnm1, float *dnm2);

void
FortranCInterface_GLOBAL(sorgl2,SORGL2)(int *m,	int *n,	int *k,	float *a, int *lda, 
	float *tau, float *work, int *info);

void
FortranCInterface_GLOBAL(sbdsqr,SBDSQR)(const char *uplo, int *n, int *ncvt, int *nru, int *ncc, float *d, 
	float *e, float *vt, int *ldvt, float *u, int *ldu,
	float *c, int *ldc, float *work, int *info);

void
FortranCInterface_GLOBAL(sgetrf,SGETRF)(int *m,	int *n,	float *a, int *lda, int *ipiv, int *info);

void
FortranCInterface_GLOBAL(sgetri,SGETRI)(int *n,	float *a, int *lda, int *ipiv, float *work, 
	int *lwork, int *info);

void
FortranCInterface_GLOBAL(sgetrs,SGETRS)(const char *trans, int *n, int *nrhs,	float *a, int *lda, int *ipiv,
	float *b, int *ldb, int *info);

void
FortranCInterface_GLOBAL(strtri,STRTRI)(const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

void
FortranCInterface_GLOBAL(strti2,STRTI2)(const char *uplo, const char *diag, int *n, float *a, int *lda, int *info);

float
FortranCInterface_GLOBAL(slange,SLANGE)(const char *norm, int *m, int *n, float *a, int *lda, float *work);

void
FortranCInterface_GLOBAL(slarrbx,SLARRBX)(int *n, float *d, float *l, float *ld, float *lld, int *ifirst,
	 int *ilast, float *rtol1, float *rtol2, int *offset, float *w,
	 float *wgap, float *werr, float *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(slasd1,SLASD1)(int *nl, int *nr, int *sqre, float *d, float *alpha, float *beta, 
	float *u, int *ldu, float *vt, int *ldvt, int *idxq, int *iwork, 
	float *work, int *info);

void
FortranCInterface_GLOBAL(slasdq,SLASDQ)(const char *uplo, int *sqre, int *n, int *ncvt, int *nru, int *ncc,
	float *d, float *e, float *vt, int *ldvt, float *u, int *ldu, 
	float *c, int *ldc, float *work, int *info);

void 
FortranCInterface_GLOBAL(slasr,SLASR)(const char *side, const char *pivot, const char *direct, int *m, int *n, float *c, 
       float *s, float *a, int *lda);

void 
FortranCInterface_GLOBAL(sorglq,SORGLQ)(int *m, int *n, int *k, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(sormtr,SORMTR)(const char *side, const char *uplo, const char *trans, int *m, int *n, float *a, 
	int *lda, float *tau, float *c, int *ldc,
	float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(sgebd2,SGEBD2)(int *m, int *n, float *a, int *lda, float *d, float *e,
	float *tauq, float *taup, float *work, int *info);

void 
FortranCInterface_GLOBAL(slabrd,SLABRD)(int *m, int *n, int *nb, float *a, int *lda, float *d,
	float *e, float *tauq, float *taup, float *x,
	int *ldx, float *y, int *ldy);

float
FortranCInterface_GLOBAL(slanst,SLANST)(const char *norm, int *n, float *d, float *e);

float
FortranCInterface_GLOBAL(slansy,SLANSY)(const char *norm, const char *uplo, int *n, float *a, int *lda, float *work);

void
FortranCInterface_GLOBAL(slarrex,SLARREX)(const char *range, int *n, float *vl, float *vu, int *il, int *iu,
	 float *d, float *e, float *tol, int *nsplit, 
	 int *isplit, int *m, float *w, int *iblock, int *indexw,
	 float *gersch, float *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(slasd2,SLASD2)(int *nl, int *nr, int *sqre, int *k, float *d, float *z, 
	float *alpha, float *beta, float *u, int *ldu, float *vt, 
	int *ldvt, float *dsigma, float *u2, int *ldu2, float *vt2, 
	int *ldvt2, int *idxp, int *idx, int *idxc, 
	int *idxq, int *coltyp, int *info);

void
FortranCInterface_GLOBAL(slasdt,SLASDT)(int *n, int *lvl, int *nd, int *inode, int *ndiml, 
	int *ndimr, int *msub);

void 
FortranCInterface_GLOBAL(slasrt,SLASRT)(const char *id, int *n, float *d, int *info);

void
FortranCInterface_GLOBAL(slasrt2,SLASRT2)(const char *id, int *n, float *d, int *key, int *info);

void 
FortranCInterface_GLOBAL(sorgqr,SORGQR)(int *m, int *n, int *k, float *a, int *lda, float *tau, 
	float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(sstebz,SSTEBZ)(const char *range, const char *order, int *n, float *vl, float *vu, 
	int *il, int *iu, float *abstol, float *d, float *e, 
	int *m, int *nsplit, float *w, int *iblock, int *isplit, 
	float *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(sgebrd,SGEBRD)(int *m, int *n, float *a, int *lda, float *d, float *e,
	float *tauq, float *taup, float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(slacpy,SLACPY)(const char *uplo, int *m, int *n, float *a, int *lda, float *b, int *ldb);

float
FortranCInterface_GLOBAL(slapy2,SLAPY2)(float * x, float * y);

void
FortranCInterface_GLOBAL(slarrfx,SLARRFX)(int *n, float *d, float *l, float *ld, float *lld, int *ifirst,
        int *ilast, float *w, float *sigma, float *dplus, float *lplus,
        float *work, int *info);

void 
FortranCInterface_GLOBAL(slasd3,SLASD3)(int *nl, int *nr, int *sqre, int *k, float *d, float *q, int *ldq, 
	float *dsigma, float *u, int *ldu, float *u2, int *ldu2, 
	float *vt, int *ldvt, float *vt2, int *ldvt2, int *idxc, 
	int *ctot, float *z, int *info);

void
FortranCInterface_GLOBAL(slaset,SLASET)(const char *uplo, int *m, int *n, float *alpha, 
	float *beta, float *a, int *lda);

void
FortranCInterface_GLOBAL(slassq,SLASSQ)(int *n, float *x, int *incx, float *scale, float *sumsq);

void
FortranCInterface_GLOBAL(sorm2l,SORM2L)(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, 
	float *tau, float *c, int *ldc, float *work, int *info);

void
FortranCInterface_GLOBAL(sstegr,SSTEGR)(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, 
	float *vu, int *il, int *iu, float *abstol, int *m, float *w, 
	float *z, int *ldz, int *isuppz, float *work, 
	int *lwork, int *iwork,	int *liwork, int *info);

void
FortranCInterface_GLOBAL(sgelq2,SGELQ2)(int *m, int *n, float *a, int *lda, float *tau, float *work, int *info);

void
FortranCInterface_GLOBAL(slae2,SLAE2)(float *a, float *b, float *c, float *rt1, float *rt2);

void
FortranCInterface_GLOBAL(slaev2,SLAEV2)(float *a, float *b, float *c, float *rt1, float *rt2,
        float *cs1, float *cs2);

void
FortranCInterface_GLOBAL(slar1vx,SLAR1VX)(int *n, int *b1, int *bn, float *sigma, float *d, float *l, float *ld, 
	float *lld, float *eval, float *gersch, float *z, float *ztz, float *mingma, 
	int *r, int *isuppz, float *work);

void
FortranCInterface_GLOBAL(slarrvx,SLARRVX)(int *n, float *d, float *l, int *isplit, int *m, float *w, 
	 int *iblock, int *indexw, float *gersch, float *tol, float *z, int *ldz, 
	 int *isuppz, float *work, int *iwork, int *info);

void 
FortranCInterface_GLOBAL(slasd4,SLASD4)(int *n, int *i, float *d, float *z, float *delta, 
	float *rho, float *sigma, float *work, int *info);

void
FortranCInterface_GLOBAL(slasq1,SLASQ1)(int *n,	float *d, float *e, float *work, int *info);


void 
FortranCInterface_GLOBAL(slasv2,SLASV2)(float *f, float *g, float *h, float *ssmin, float *ssmax, 
	float *snr, float *csr, float *snl, float *csl);

void 
FortranCInterface_GLOBAL(sorm2r,SORM2R)(const char *side, const char *trans, int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *c, int *ldc, float *work, int *info);

void
FortranCInterface_GLOBAL(sstein,SSTEIN)(int *n, float *d, float *e, int *m, float *w, int *iblock, int *isplit, 
	float *z, int *ldz, float *work, int *iwork, int *ifail, int *info);

void
FortranCInterface_GLOBAL(sgelqf,SGELQF)(int *m,	int *n, float *a, int *lda, float *tau,
	float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(slaebz,SLAEBZ)(int *ijob, int *nitmax, int *n, int *mmax, int *minp, int *nbmin,
	float *abstol, float *reltol, float *pivmin, float *d, float *e,
	float *e2, int *nval, float *ab, float *c, int *mout, int *nab,
	float *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(slarf,SLARF)(const char *side, int *m, int *n, float *v, int *incv, float *tau,
       float *c, int *ldc, float *work);

void
FortranCInterface_GLOBAL(slartg,SLARTG)(float *f, float *g, float *cs, float *sn, float *r);

void 
FortranCInterface_GLOBAL(slasd5,SLASD5)(int *i, float *d, float *z, float *delta, 
	float *rho, float *dsigma, float *work);

void 
FortranCInterface_GLOBAL(slasq2,SLASQ2)(int *n, float *z, int *info);

void 
FortranCInterface_GLOBAL(slasq3,SLASQ3)(int *i0, int *n0, float *z, int *pp, float *dmin, 
	float *sigma, float *desig, float *qmax, int *nfail, 
	int *iter, int *ndiv, int *ieee);

void
FortranCInterface_GLOBAL(slaswp,SLASWP)(int *n,	float *a, int *lda, int *k1, int *k2, int *ipiv, int *incx);

void 
FortranCInterface_GLOBAL(sormbr,SORMBR)(const char *vect, const char *side, const char *trans, int *m, int *n, int *k, 
	float *a, int *lda, float *tau, float *c, int *ldc, float *work,
	int *lwork, int *info);

void
FortranCInterface_GLOBAL(ssterf,SSTERF)(int *n, float *d, float *e, int *info);

void
FortranCInterface_GLOBAL(sgeqr2,SGEQR2)(int *m,	int *n,	float *a, int *lda, float *tau, 
	float *work, int *info);

void 
FortranCInterface_GLOBAL(slaed6,SLAED6)(int *kniter, int *orgati, float *rho, float *d, 
	float *z, float *finit, float *tau, int *info);

void 
FortranCInterface_GLOBAL(slarfb,SLARFB)(const char *side, const char *trans, const char *direct, const char *storev, int *m, int *n, 
	int *k, float *v, int *ldv, float *t, int *ldt, float *c,
	int *ldc, float *work, int *ldwork);

void
FortranCInterface_GLOBAL(slaruv,SLARUV)(int *iseed, int *n, float *x);

void 
FortranCInterface_GLOBAL(slasd6,SLASD6)(int *icompq, int *nl, int *nr, int *sqre, float *d, float *vf, 
	float *vl, float *alpha, float *beta, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, float *givnum, int *ldgnum, 
	float *poles, float *difl, float *difr, float *z, int *k, 
	float *c, float *s, float *work, int *iwork, int *info);

void
FortranCInterface_GLOBAL(slatrd,SLATRD)(const char *uplo, int *n, int *nb, float *a, int *lda, float *e, 
	float * tau, float *w, int *ldw);

void
FortranCInterface_GLOBAL(sorml2,SORML2)(const char *side, const char *trans, int *m, int *n, int *k, float *a,
	int *lda, float *tau, float *c, int *ldc, float *work, int *info);

void
FortranCInterface_GLOBAL(sstevr,SSTEVR)(const char *jobz, const char *range, int *n, float *d, float *e, float *vl, 
	float *vu, int *il, int *iu, float *abstol, int *m, float *w, 
	float *z, int *ldz, int *isuppz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);

void
FortranCInterface_GLOBAL(ssytrd,SSYTRD)(const char *uplo, int *n, float *  a, int *lda, float *d, 
	float *e, float *tau, float *work, int *lwork, int *info);

GMX_LIBGMX_EXPORT
void
FortranCInterface_GLOBAL(ssyevr,SSYEVR)(const char *jobz, const char *range, const char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *
	il, int *iu, float *abstol, int *m, float *w, 
	float *z__, int *ldz, int *isuppz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);

void
FortranCInterface_GLOBAL(sormql,SORMQL)(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *
	c, int *ldc, float *work, int *lwork, int *info);

void 
FortranCInterface_GLOBAL(sormqr,SORMQR)(const char *side, const char *trans, int *m, int *n, int *k, float *a, 
        int *lda, float *tau, float *c, int *ldc, 
        float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(sorgbr,SORGBR)(const char *vect, int *m, int *n, int *k, float *a, int *lda,
	float *tau, float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(slasq5,SLASQ5)(int *i0, int *n0, float *z, int *pp, float *tau, float *dmin, 
	float *dmin1, float *dmin2, float *dn, float *dnm1, 
	float *dnm2, int *ieee);

void 
FortranCInterface_GLOBAL(slasd8,SLASD8)(int *icompq, int *k, float *d, float *z, float *vf, float *vl, 
	float *difl, float *difr, int *lddifr, float *dsigma, 
	float *work, int *info);

void
FortranCInterface_GLOBAL(slascl,SLASCL)(const char *type, int *kl, int *ku, float *cfrom, float *cto, int *m, 
	int *n, float *a, int *lda, int *info);

void 
FortranCInterface_GLOBAL(slarft,SLARFT)(const char *direct, const char *storev, int *n, int *k, float *v, 
	int *ldv, float *tau, float *t, int *ldt);

void
FortranCInterface_GLOBAL(slagts,SLAGTS)(int *job, int *n, float *a, float *b, float *c, float *d, 
	int *in, float *y, float *tol, int *info);

void 
FortranCInterface_GLOBAL(sgesdd,SGESDD)(const char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, 
	int *ldu, float *vt, int *ldvt, float *work, int *lwork, 
	int *iwork, int *info);

void
FortranCInterface_GLOBAL(ssytd2,SSYTD2)(const char *uplo, int *n, float *a, int *lda, float *d, 
	float *e, float *tau, int *info);

void 
FortranCInterface_GLOBAL(sormlq,SORMLQ)(const char *side, const char *trans, int *m, int *n, int *k, float *a, int *lda, 
	float *tau, float *c, int *ldc, float *work, int *lwork, int *info);

void
FortranCInterface_GLOBAL(sorg2r,SORG2R)(int *m, int *n,	int *k,	float *a, int *lda, float *tau,
	float *work, int *info);

void 
FortranCInterface_GLOBAL(slasq4,SLASQ4)(int *i0, int *n0, float *z, int *pp, int *n0in, float *dmin, 
	float *dmin1, float *dmin2, float *dn, float *dn1, 
	float *dn2, float *tau, int *ttype);

void 
FortranCInterface_GLOBAL(slasd7,SLASD7)(int *icompq, int *nl, int *nr, int *sqre, int *k, float *d, float *z,
	float *zw, float *vf, float *vfw, float *vl, float *vlw,
	float *alpha, float *beta, float *dsigma, int *idx, int *idxp,
	int *idxq, int *perm, int *givptr, int *givcol, int *ldgcol, 
	float *givnum, int *ldgnum, float *c, float *s, int *info);

void
FortranCInterface_GLOBAL(slas2,SLAS2)(float *f, float *g, float *h, float *ssmin, float *ssmax);

void
FortranCInterface_GLOBAL(slarfg,SLARFG)(int *n, float *alpha, float *x, int *incx, float *tau);

void
FortranCInterface_GLOBAL(slagtf,SLAGTF)(int *n, float *a, float *lambda, float *b, float *c, 
	float *tol, float *d, int *in, int *info);

void 
FortranCInterface_GLOBAL(sgeqrf,SGEQRF)(int *m, int *n, float *a, int *lda, float *tau,
	float *work, int *lwork, int *info);


#ifdef __cplusplus
}
#endif



#endif /* _LAPACK_H_ */
