#include "../gmx_blas.h"
#include "../gmx_lapack.h"

void 
F77_FUNC(dlasda,DLASDA)(int *icompq, 
	int *smlsiz, 
	int *n, 
	int *sqre, 
	double *d__, 
	double *e, 
	double *u, 
	int *ldu, 
	double *vt, 
	int *k, 
	double *difl, 
	double *difr, 
	double *z__, 
	double *poles, 
	int *givptr, 
	int *givcol, 
	int *ldgcol, 
	int *perm, 
	double *givnum, 
	double *c__, 
	double *s, 
	double *work, 
	int *iwork, 
	int *info)
{
    int givcol_dim1, givcol_offset, perm_dim1, perm_offset, difl_dim1, 
	    difl_offset, difr_dim1, difr_offset, givnum_dim1, givnum_offset, 
	    poles_dim1, poles_offset, u_dim1, u_offset, vt_dim1, vt_offset, 
	    z_dim1, z_offset, i__1, i__2;

    int i__, j, m, i1, ic, lf, nd, ll, nl, vf, nr, vl, im1, ncc, 
	    nlf, nrf, vfi, iwk, vli, lvl, nru, ndb1, nlp1, lvl2, nrp1;
    double beta;
    int idxq, nlvl;
    double alpha;
    int inode, ndiml, ndimr, idxqi, itemp;
    int sqrei;
    int nwork1, nwork2;
    int smlszp;
    int c__0 = 0;
    double zero = 0.0;
    double one = 1.;
    int c__1 = 1;
    --d__;
    --e;
    givnum_dim1 = *ldu;
    givnum_offset = 1 + givnum_dim1;
    givnum -= givnum_offset;
    poles_dim1 = *ldu;
    poles_offset = 1 + poles_dim1;
    poles -= poles_offset;
    z_dim1 = *ldu;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    difr_dim1 = *ldu;
    difr_offset = 1 + difr_dim1;
    difr -= difr_offset;
    difl_dim1 = *ldu;
    difl_offset = 1 + difl_dim1;
    difl -= difl_offset;
    vt_dim1 = *ldu;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --k;
    --givptr;
    perm_dim1 = *ldgcol;
    perm_offset = 1 + perm_dim1;
    perm -= perm_offset;
    givcol_dim1 = *ldgcol;
    givcol_offset = 1 + givcol_dim1;
    givcol -= givcol_offset;
    --c__;
    --s;
    --work;
    --iwork;
    *info = 0;

    m = *n + *sqre;

    if (*n <= *smlsiz) {
	if (*icompq == 0) {
	    F77_FUNC(dlasdq,DLASDQ)("U", sqre, n, &c__0, &c__0, &c__0, &d__[1], &e[1], &vt[
		    vt_offset], ldu, &u[u_offset], ldu, &u[u_offset], ldu, &
		    work[1], info);
	} else {
	    F77_FUNC(dlasdq,DLASDQ)("U", sqre, n, &m, n, &c__0, &d__[1], &e[1], &vt[vt_offset]
		    , ldu, &u[u_offset], ldu, &u[u_offset], ldu, &work[1], 
		    info);
	}
	return;
    }

    inode = 1;
    ndiml = inode + *n;
    ndimr = ndiml + *n;
    idxq = ndimr + *n;
    iwk = idxq + *n;

    ncc = 0;
    nru = 0;

    smlszp = *smlsiz + 1;
    vf = 1;
    vl = vf + m;
    nwork1 = vl + m;
    nwork2 = nwork1 + smlszp * smlszp;

    F77_FUNC(dlasdt,DLASDT)(n, &nlvl, &nd, &iwork[inode], &iwork[ndiml], &iwork[ndimr], 
	    smlsiz);

    ndb1 = (nd + 1) / 2;
    i__1 = nd;
    for (i__ = ndb1; i__ <= i__1; ++i__) {
	i1 = i__ - 1;
	ic = iwork[inode + i1];
	nl = iwork[ndiml + i1];
	nlp1 = nl + 1;
	nr = iwork[ndimr + i1];
	nlf = ic - nl;
	nrf = ic + 1;
	idxqi = idxq + nlf - 2;
	vfi = vf + nlf - 1;
	vli = vl + nlf - 1;
	sqrei = 1;
	if (*icompq == 0) {
	    F77_FUNC(dlaset,DLASET)("A", &nlp1, &nlp1, &zero, &one, &work[nwork1], &smlszp);
	    F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nl, &nlp1, &nru, &ncc, &d__[nlf], &e[nlf], &
		    work[nwork1], &smlszp, &work[nwork2], &nl, &work[nwork2], 
		    &nl, &work[nwork2], info);
	    itemp = nwork1 + nl * smlszp;
	    F77_FUNC(dcopy,DCOPY)(&nlp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    F77_FUNC(dcopy,DCOPY)(&nlp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    F77_FUNC(dlaset,DLASET)("A", &nl, &nl, &zero, &one, &u[nlf + u_dim1], ldu);
	    F77_FUNC(dlaset,DLASET)("A", &nlp1, &nlp1, &zero, &one, &vt[nlf + vt_dim1], 
		    ldu);
	    F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nl, &nlp1, &nl, &ncc, &d__[nlf], &e[nlf], &
		    vt[nlf + vt_dim1], ldu, &u[nlf + u_dim1], ldu, &u[nlf + 
		    u_dim1], ldu, &work[nwork1], info);
	    F77_FUNC(dcopy,DCOPY)(&nlp1, &vt[nlf + vt_dim1], &c__1, &work[vfi], &c__1);
	    F77_FUNC(dcopy,DCOPY)(&nlp1, &vt[nlf + nlp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nl;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
	if (i__ == nd && *sqre == 0) {
	    sqrei = 0;
	} else {
	    sqrei = 1;
	}
	idxqi += nlp1;
	vfi += nlp1;
	vli += nlp1;
	nrp1 = nr + sqrei;
	if (*icompq == 0) {
	    F77_FUNC(dlaset,DLASET)("A", &nrp1, &nrp1, &zero, &one, &work[nwork1], &smlszp);
	    F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nr, &nrp1, &nru, &ncc, &d__[nrf], &e[nrf], &
		    work[nwork1], &smlszp, &work[nwork2], &nr, &work[nwork2], 
		    &nr, &work[nwork2], info);
	    itemp = nwork1 + (nrp1 - 1) * smlszp;
	    F77_FUNC(dcopy,DCOPY)(&nrp1, &work[nwork1], &c__1, &work[vfi], &c__1);
	    F77_FUNC(dcopy,DCOPY)(&nrp1, &work[itemp], &c__1, &work[vli], &c__1);
	} else {
	    F77_FUNC(dlaset,DLASET)("A", &nr, &nr, &zero, &one, &u[nrf + u_dim1], ldu);
	    F77_FUNC(dlaset,DLASET)("A", &nrp1, &nrp1, &zero, &one, &vt[nrf + vt_dim1], 
		    ldu);
	    F77_FUNC(dlasdq,DLASDQ)("U", &sqrei, &nr, &nrp1, &nr, &ncc, &d__[nrf], &e[nrf], &
		    vt[nrf + vt_dim1], ldu, &u[nrf + u_dim1], ldu, &u[nrf + 
		    u_dim1], ldu, &work[nwork1], info);
	    F77_FUNC(dcopy,DCOPY)(&nrp1, &vt[nrf + vt_dim1], &c__1, &work[vfi], &c__1);
	    F77_FUNC(dcopy,DCOPY)(&nrp1, &vt[nrf + nrp1 * vt_dim1], &c__1, &work[vli], &c__1)
		    ;
	}
	if (*info != 0) {
	    return;
	}
	i__2 = nr;
	for (j = 1; j <= i__2; ++j) {
	    iwork[idxqi + j] = j;
	}
    }

    j = (1 << nlvl);

    for (lvl = nlvl; lvl >= 1; --lvl) {
	lvl2 = (lvl << 1) - 1;

	if (lvl == 1) {
	    lf = 1;
	    ll = 1;
	} else {
	    i__1 = lvl - 1;
	    lf = (1 << (lvl-1));
	    ll = (lf << 1) - 1;
	}
	i__1 = ll;
	for (i__ = lf; i__ <= i__1; ++i__) {
	    im1 = i__ - 1;
	    ic = iwork[inode + im1];
	    nl = iwork[ndiml + im1];
	    nr = iwork[ndimr + im1];
	    nlf = ic - nl;
	    nrf = ic + 1;
	    if (i__ == ll) {
		sqrei = *sqre;
	    } else {
		sqrei = 1;
	    }
	    vfi = vf + nlf - 1;
	    vli = vl + nlf - 1;
	    idxqi = idxq + nlf - 1;
	    alpha = d__[ic];
	    beta = e[ic];
	    if (*icompq == 0) {
		F77_FUNC(dlasd6,DLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[
			perm_offset], &givptr[1], &givcol[givcol_offset], 
			ldgcol, &givnum[givnum_offset], ldu, &poles[
			poles_offset], &difl[difl_offset], &difr[difr_offset],
			 &z__[z_offset], &k[1], &c__[1], &s[1], &work[nwork1],
			 &iwork[iwk], info);
	    } else {
		--j;
		F77_FUNC(dlasd6,DLASD6)(icompq, &nl, &nr, &sqrei, &d__[nlf], &work[vfi], &
			work[vli], &alpha, &beta, &iwork[idxqi], &perm[nlf + 
			lvl * perm_dim1], &givptr[j], &givcol[nlf + lvl2 * 
			givcol_dim1], ldgcol, &givnum[nlf + lvl2 * 
			givnum_dim1], ldu, &poles[nlf + lvl2 * poles_dim1], &
			difl[nlf + lvl * difl_dim1], &difr[nlf + lvl2 * 
			difr_dim1], &z__[nlf + lvl * z_dim1], &k[j], &c__[j], 
			&s[j], &work[nwork1], &iwork[iwk], info);
	    }
	    if (*info != 0) {
		return;
	    }
	}
    }

    return;

}


