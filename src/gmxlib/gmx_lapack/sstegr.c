#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"

#include <types/simple.h>

void
F77_FUNC(sstegr,SSTEGR)(char *jobz, 
	char *range, 
	int *n, 
	float *d__, 
	float *e, 
	float *vl, 
	float *vu, 
	int *il, 
	int *iu, 
	float *abstol, 
	int *m, 
	float *w, 
	float *z__, 
	int *ldz, 
	int *isuppz,
	float *work, 
	int *lwork, 
	int *iwork, 
	int *liwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2;
    float d__1, d__2;
    int c__1 = 1;

    int i__, j;
    int jj;
    float eps, tol, tmp, rmin, rmax;
    int itmp;
    float tnrm;
    float scale;
    int iinfo, iindw;
    int lwmin;
    int wantz;
    int iindbl;
    int valeig,alleig,indeig;
    float safmin,minval;
    float bignum;
    int iindwk, indgrs;
    float thresh;
    int iinspl, indwrk, liwmin, nsplit;
    float smlnum;
    int lquery;


    --d__;
    --e;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    lquery = *lwork == -1 || *liwork == -1;
    lwmin = *n * 17;
    liwmin = *n * 10;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (valeig && *n > 0 && *vu <= *vl) {
	*info = -7;
    } else if (indeig && (*il < 1 || *il > *n)) {
	*info = -8;
    } else if (indeig && (*iu < *il || *iu > *n)) {
	*info = -9;
    } else if (*ldz < 1 || (wantz && *ldz < *n)) {
	*info = -14;
    } else if (*lwork < lwmin && ! lquery) {
	*info = -17;
    } else if (*liwork < liwmin && ! lquery) {
	*info = -19;
    }
    if (*info == 0) {
	work[1] = (float) lwmin;
	iwork[1] = liwmin;
    }

    if (*info != 0) {
	i__1 = -(*info);
	return;
    } else if (lquery) {
	return;
    }

    *m = 0;
    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = d__[1];
	} else {
	    if (*vl < d__[1] && *vu >= d__[1]) {
		*m = 1;
		w[1] = d__[1];
	    }
	}
	if (wantz) {
	    z__[z_dim1 + 1] = 1.;
	}
	return;
    }

    minval = GMX_FLOAT_MIN;
    safmin = minval*(1.0+GMX_FLOAT_EPS);
    eps = GMX_FLOAT_EPS;
    smlnum = safmin / eps;
    bignum = 1. / smlnum;
    rmin = sqrt(smlnum);
    d__1 = sqrt(bignum), d__2 = 1. / sqrt(sqrt(safmin));
    rmax = (d__1<d__2) ? d__1 : d__2;
    scale = 1.;
    tnrm = F77_FUNC(slanst,SLANST)("M", n, &d__[1], &e[1]);
    if (tnrm > 0. && tnrm < rmin) {
	scale = rmin / tnrm;
    } else if (tnrm > rmax) {
	scale = rmax / tnrm;
    }
    if ( fabs(scale-1.0)>GMX_FLOAT_EPS) {
	F77_FUNC(sscal,SSCAL)(n, &scale, &d__[1], &c__1);
	i__1 = *n - 1;
	F77_FUNC(sscal,SSCAL)(&i__1, &scale, &e[1], &c__1);
	tnrm *= scale;
    }
    indgrs = 1;
    indwrk = (*n << 1) + 1;

    iinspl = 1;
    iindbl = *n + 1;
    iindw = (*n << 1) + 1;
    iindwk = *n * 3 + 1;

    thresh = eps * tnrm;
    F77_FUNC(slarrex,SLARREX)(range, n, vl, vu, il, iu, &d__[1], &e[1], &thresh, &nsplit, &
	    iwork[iinspl], m, &w[1], &iwork[iindbl], &iwork[iindw], &work[
	    indgrs], &work[indwrk], &iwork[iindwk], &iinfo);
    
    if (iinfo != 0) {
	*info = 1;
	return;
    }

    if (wantz) {
	d__1 = *abstol, d__2 = (float) (*n) * eps;
	tol = (d__1>d__2) ? d__1 : d__2;
	F77_FUNC(slarrvx,SLARRVX)(n, &d__[1], &e[1], &iwork[iinspl], m, &w[1], &iwork[iindbl], &
		iwork[iindw], &work[indgrs], &tol, &z__[z_offset], ldz, &
		isuppz[1], &work[indwrk], &iwork[iindwk], &iinfo);
	if (iinfo != 0) {
	    *info = 2;
	    return;
	}
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	itmp = iwork[iindbl + j - 1];
	w[j] += e[iwork[iinspl + itmp - 1]];
    } 

    if (fabs(scale-1.0)>GMX_FLOAT_EPS) {
	d__1 = 1. / scale;
	F77_FUNC(sscal,SSCAL)(m, &d__1, &w[1], &c__1);
    }
    if (nsplit > 1) {
	i__1 = *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp) {
		    i__ = jj;
		    tmp = w[jj];
		}
	    }
	    if (i__ != 0) {
		w[i__] = w[j];
		w[j] = tmp;
		if (wantz) {
		    F77_FUNC(sswap,SSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 
			    + 1], &c__1);
		    itmp = isuppz[(i__ << 1) - 1];
		    isuppz[(i__ << 1) - 1] = isuppz[(j << 1) - 1];
		    isuppz[(j << 1) - 1] = itmp;
		    itmp = isuppz[i__ * 2];
		    isuppz[i__ * 2] = isuppz[j * 2];
		    isuppz[j * 2] = itmp;
		}
	    }
	}
    }

    work[1] = (float) lwmin;
    iwork[1] = liwmin;
    return;

} 
