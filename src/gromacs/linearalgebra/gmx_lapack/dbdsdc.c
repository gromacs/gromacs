#include <ctype.h>
#include <math.h>
#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

#include "gromacs/utility/real.h"

void
F77_FUNC(dbdsdc,DBDSDC)(const char *uplo, 
	const char *compq, 
	int *n,
	double *d__, 
	double *e, 
	double *u, 
	int *ldu,
	double *vt, 
	int *ldvt,
	double *q,
	int *iq,
	double *work, 
	int *iwork, 
	int *info)
{
    int u_dim1, u_offset, vt_dim1, vt_offset, i__1, i__2;
    int i__, j, k;
    double p, r__;
    int z__, ic, ii, kk;
    double cs;
    int is, iu;
    double sn;
    int nm1;
    double eps;
    int ivt, difl, difr, ierr, perm, mlvl, sqre;
    int poles, iuplo, nsize, start;
    int givcol;
    int icompq;
    double orgnrm;
    int givnum, givptr, qstart, smlsiz, wstart, smlszp;
    double zero = 0.0;
    double one = 1.0;
    int c_0 = 0;
    int c_1 = 1;

    --d__;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    vt_dim1 = *ldvt;
    vt_offset = 1 + vt_dim1;
    vt -= vt_offset;
    --q;
    --iq;
    --work;
    --iwork;

    k = iu = z__ = ic = is = ivt = difl = difr = perm = 0;
    poles = givnum = givptr = givcol = 0;
    
    smlsiz = DBDSDC_SMALLSIZE;
    *info = 0;

    iuplo = (*uplo=='U' || *uplo=='u') ? 1 : 2;

    switch(*compq) {
    case 'n':
    case 'N':
      icompq = 0;
      break;
    case 'p':
    case 'P':
      icompq = 1;
      break;
    case 'i':
    case 'I':
      icompq = 2;
      break;
    default:
      return;
    }

    if (*n <= 0) 
	return;
    
    if (*n == 1) {
	if (icompq == 1) {
	  q[1] = (d__[1]>0) ? 1.0 : -1.0;
	  q[smlsiz * *n + 1] = 1.;
	} else if (icompq == 2) {
	  u[u_dim1 + 1] = (d__[1]>0) ? 1.0 : -1.0;
	  vt[vt_dim1 + 1] = 1.;
	}
	d__[1] = fabs(d__[1]);
	return;
    }
    nm1 = *n - 1;
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
	F77_FUNC(dcopy,DCOPY)(n, &d__[1], &c_1, &q[1], &c_1);
	i__1 = *n - 1;
	F77_FUNC(dcopy,DCOPY)(&i__1, &e[1], &c_1, &q[*n + 1], &c_1);
    }
    if (iuplo == 2) {
	qstart = 5;
	wstart = (*n << 1) - 1;
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    F77_FUNC(dlartg,DLARTG)(&d__[i__], &e[i__], &cs, &sn, &r__);
	    d__[i__] = r__;
	    e[i__] = sn * d__[i__ + 1];
	    d__[i__ + 1] = cs * d__[i__ + 1];
	    if (icompq == 1) {
		q[i__ + (*n << 1)] = cs;
		q[i__ + *n * 3] = sn;
	    } else if (icompq == 2) {
		work[i__] = cs;
		work[nm1 + i__] = -sn;
	    }
	}
    }
    if (icompq == 0) {
      F77_FUNC(dlasdq,DLASDQ)("U",&c_0,n,&c_0,&c_0,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
	      &u[u_offset], ldu, &u[u_offset], ldu, &work[wstart], info);
	goto L40;
    }
    if (*n <= smlsiz) {
	if (icompq == 2) {
	    F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	    F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
	    F77_FUNC(dlasdq,DLASDQ)("U",&c_0,n,n,n,&c_0,&d__[1],&e[1],&vt[vt_offset],ldvt,
		    &u[u_offset],ldu,&u[u_offset],ldu,&work[wstart],info);
	} else if (icompq == 1) {
	    iu = 1;
	    ivt = iu + *n;
	    F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &q[iu + (qstart - 1) * *n], n);
	    F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &q[ivt + (qstart - 1) * *n], n);
	    F77_FUNC(dlasdq,DLASDQ)("U", &c_0, n, n, n, &c_0, &d__[1], &e[1], 
		    &q[ivt + (qstart - 1) * *n], n, &q[iu + (qstart - 1) * *n], 
		    n, &q[iu + (qstart - 1) * *n], n, &work[wstart], info);
	}
	goto L40;
    }

    if (icompq == 2) {
	F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &u[u_offset], ldu);
	F77_FUNC(dlaset,DLASET)("A", n, n, &zero, &one, &vt[vt_offset], ldvt);
    }

    orgnrm = F77_FUNC(dlanst,DLANST)("M", n, &d__[1], &e[1]);
    if ( fabs(orgnrm)<GMX_DOUBLE_MIN) {
	return;
    }
    F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &orgnrm, &one, n, &c_1, &d__[1], n, &ierr);
    F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &orgnrm, &one, &nm1, &c_1, &e[1], &nm1, &ierr);

    eps = GMX_DOUBLE_EPS;

    mlvl = (int) (log((double) (*n) / (double) (smlsiz + 1)) / 
	    log(2.)) + 1;
    smlszp = smlsiz + 1;

    if (icompq == 1) {
	iu = 1;
	ivt = smlsiz + 1;
	difl = ivt + smlszp;
	difr = difl + mlvl;
	z__ = difr + (mlvl << 1);
	ic = z__ + mlvl;
	is = ic + 1;
	poles = is + 1;
	givnum = poles + (mlvl << 1);

	k = 1;
	givptr = 2;
	perm = 3;
	givcol = perm + mlvl;
    }

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(d__[i__]) < eps) 
	    d__[i__] = (d__[i__]>0) ? eps : -eps;
    }

    start = 1;
    sqre = 0;

    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (fabs(e[i__]) < eps || i__ == nm1) {
	    if (i__ < nm1) {
		nsize = i__ - start + 1;
	    } else if (fabs(e[i__]) >= eps) {
		nsize = *n - start + 1;
	    } else {
		nsize = i__ - start + 1;
		if (icompq == 2) {
		    u[*n + *n * u_dim1] = (d__[*n]>0) ? 1.0 : -1.0; 
		    vt[*n + *n * vt_dim1] = 1.;
		} else if (icompq == 1) {
		    q[*n + (qstart - 1) * *n] = (d__[*n]>0) ? 1.0 : -1.0; 
		    q[*n + (smlsiz + qstart - 1) * *n] = 1.;
		}
		d__[*n] = fabs(d__[*n]);
	    }
	    if (icompq == 2) {
		F77_FUNC(dlasd0,DLASD0)(&nsize, &sqre, &d__[start], &e[start], 
			&u[start + start * u_dim1], ldu, 
			&vt[start + start * vt_dim1], 
			ldvt, &smlsiz, &iwork[1], &work[wstart], info);
	    } else {
		F77_FUNC(dlasda,DLASDA)(&icompq, &smlsiz, &nsize, &sqre, &d__[start], 
			&e[start], &q[start + (iu + qstart - 2) * *n], n, 
			&q[start + (ivt + qstart - 2) * *n], &iq[start + k * *n],
			&q[start + (difl + qstart - 2) * *n], 
			&q[start + (difr + qstart - 2) * *n], 
			&q[start + (z__ + qstart - 2) * *n], 
			&q[start + (poles + qstart - 2) * *n], 
			&iq[start + givptr * *n], &iq[start + givcol * *n], n, 
			&iq[start + perm * *n], 
			&q[start + (givnum + qstart - 2) * *n], 
			&q[start + (ic + qstart - 2) * *n], 
			&q[start + (is + qstart - 2) * *n], &work[wstart], 
			&iwork[1], info);
		if (*info != 0) {
		    return;
		}
	    }
	    start = i__ + 1;
	}
    }
    F77_FUNC(dlascl,DLASCL)("G", &c_0, &c_0, &one, &orgnrm, n, &c_1, &d__[1], n, &ierr);
L40:
    i__1 = *n;
    for (ii = 2; ii <= i__1; ++ii) {
	i__ = ii - 1;
	kk = i__;
	p = d__[i__];
	i__2 = *n;
	for (j = ii; j <= i__2; ++j) {
	    if (d__[j] > p) {
		kk = j;
		p = d__[j];
	    }
	}
	if (kk != i__) {
	    d__[kk] = d__[i__];
	    d__[i__] = p;
	    if (icompq == 1) {
		iq[i__] = kk;
	    } else if (icompq == 2) {
		F77_FUNC(dswap,DSWAP)(n, &u[i__ * u_dim1 + 1],&c_1,&u[kk*u_dim1+1],&c_1);
		F77_FUNC(dswap,DSWAP)(n, &vt[i__ + vt_dim1], ldvt, &vt[kk + vt_dim1], ldvt);
	    }
	} else if (icompq == 1) {
	    iq[i__] = i__;
	}
    }
    if (icompq == 1) {
	if (iuplo == 1) {
	    iq[*n] = 1;
	} else {
	    iq[*n] = 0;
	}
    }
    if (iuplo == 2 && icompq == 2) {
	F77_FUNC(dlasr,DLASR)("L", "V", "B", n, n, &work[1], &work[*n], &u[u_offset], ldu);
    }

    return;
}
