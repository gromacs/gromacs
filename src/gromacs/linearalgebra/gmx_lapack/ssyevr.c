#include <math.h>

#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"

void
F77_FUNC(ssyevr,SSYEVR)(const char *jobz, const char *range, const char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *
	il, int *iu, float *abstol, int *m, float *w, 
	float *z__, int *ldz, int *isuppz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    float d__1, d__2;

    /* Local variables */
    int c__1 = 1;
    int i__, j, nb, jj;
    float eps, vll, vuu, tmp1;
    int indd, inde;
    float anrm;
    int imax;
    float rmin, rmax;
    int itmp1, inddd, indee;
    float sigma;
    int iinfo;
    int indwk;
    int lwmin;
    int lower, wantz;
    int alleig, indeig;
    int iscale, ieeeok, indibl, indifl;
    int valeig;
    float safmin,minval;
    float abstll, bignum;
    int indtau;
    int indwkn;
    int liwmin;
    int llwrkn, llwork;
    float smlnum;
    int lwkopt;
    int lquery;
    float posinf,neginf,negzro,newzro;
    float fzero = 0.0;
    float fone = 1.0;

    
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    /* Check for IEEE-compliant FP */
    ieeeok = 1;
    posinf = fone/fzero;
    if(posinf<=1.0)
      ieeeok = 0;
    neginf = -fone/fzero;
    if(neginf>=0.0)
      ieeeok = 0;
    negzro = fone/(neginf+fone);
    if(negzro!=0)
      ieeeok = 0;
    neginf = fone/negzro;
    if(neginf>=0)
      ieeeok = 0;
    newzro = negzro + fzero;
    if(newzro!=fzero)
      ieeeok = 0;
    posinf = fone /newzro;
    if(posinf<=fone)
      ieeeok = 0;
    neginf = neginf*posinf;
    if(neginf>=fzero)
      ieeeok = 0;
    posinf = posinf*posinf;
    if(posinf<=1.0)
      ieeeok = 0;

    lower = (*uplo=='L' || *uplo=='l');
    wantz = (*jobz=='V' || *jobz=='v');
    alleig = (*range=='A' || *range=='a');
    valeig = (*range=='V' || *range=='v');
    indeig = (*range=='I' || *range=='i');

    indibl = 0;
    lquery = *lwork == -1 || *liwork == -1;

    i__1 = 1;
    i__2 = *n * 26;

    if(*n>0) 
      lwmin = *n * 26;
    else
      lwmin = 1;

    if(*n>0) 
      liwmin = *n * 10;
    else
      liwmin = 1;

    *info = 0;
    if (! (wantz || (*jobz=='N' || *jobz=='n'))) {
	*info = -1;
    } else if (! (alleig || valeig || indeig)) {
	*info = -2;
    } else if (! (lower || (*uplo=='U' || *uplo=='u'))) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < ((*n>1) ? *n : 1) ) {
	*info = -6;
    } else {
	if (valeig) {
	    if (*n > 0 && *vu <= *vl) {
		*info = -8;
	    }
	} else if (indeig) {
	  if (*il < 1 || *il > ((*n>1) ? *n : 1)) {
		*info = -9;
	    } else if (*iu < ((*n<*il) ? *n : *il) || *iu > *n) {
		*info = -10;
	    }
	}
    }
    if (*info == 0) {
      if (*ldz < 1 || (wantz && *ldz < *n)) {
	    *info = -15;
	} else if (*lwork < lwmin && ! lquery) {
	    *info = -18;
	} else if (*liwork < liwmin && ! lquery) {
	    *info = -20;
	}
    }

    if (*info == 0) {
      nb = 32;
      /* Computing MAX */
      i__1 = (nb + 1) * *n;
      lwkopt = (i__1>lwmin) ? i__1 : lwmin;
      work[1] = (float) lwkopt;
      iwork[1] = liwmin;
    } else 
      return;

    if (lquery) 
	return;
    
    *m = 0;
    if (*n == 0) {
	work[1] = 1.;
	return;
    }

    if (*n == 1) {
	work[1] = 7.;
	if (alleig || indeig) {
	    *m = 1;
	    w[1] = a[a_dim1 + 1];
	} else {
	    if (*vl < a[a_dim1 + 1] && *vu >= a[a_dim1 + 1]) {
		*m = 1;
		w[1] = a[a_dim1 + 1];
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

    iscale = 0;
    abstll = *abstol;
    vll = *vl;
    vuu = *vu;
    anrm = F77_FUNC(slansy,SLANSY)("M", uplo, n, &a[a_offset], lda, &work[1]);
    if (anrm > 0. && anrm < rmin) {
	iscale = 1;
	sigma = rmin / anrm;
    } else if (anrm > rmax) {
	iscale = 1;
	sigma = rmax / anrm; 
    }
    if (iscale == 1) {
	if (lower) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - j + 1;
		F77_FUNC(sscal,SSCAL)(&i__2, &sigma, &a[j + j * a_dim1], &c__1);
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		F77_FUNC(sscal,SSCAL)(&j, &sigma, &a[j * a_dim1 + 1], &c__1);

	    }
	}
	if (*abstol > 0.) {
	    abstll = *abstol * sigma;
	}
	if (valeig) {
	    vll = *vl * sigma;
	    vuu = *vu * sigma;
	}
    }

    indtau = 1;
    inde = indtau + *n;
    indd = inde + *n;
    indee = indd + *n;
    inddd = indee + *n;
    indifl = inddd + *n;
    indwk = indifl + *n;
    llwork = *lwork - indwk + 1;
    F77_FUNC(ssytrd,SSYTRD)(uplo, n, &a[a_offset], lda, &work[indd], &work[inde], &work[
	    indtau], &work[indwk], &llwork, &iinfo);

    i__1 = *n - 1;
    F77_FUNC(scopy,SCOPY)(&i__1, &work[inde], &c__1, &work[indee], &c__1);
    F77_FUNC(scopy,SCOPY)(n, &work[indd], &c__1, &work[inddd], &c__1);

    F77_FUNC(sstegr,SSTEGR)(jobz, range, n, &work[inddd], &work[indee], vl, vu, il, iu, 
	    abstol, m, &w[1], &z__[z_offset], ldz, &isuppz[1], 
	    &work[indwk], lwork, &iwork[1], liwork, info);
    if (wantz && *info == 0) {
      indwkn = inde;
      llwrkn = *lwork - indwkn + 1;
      F77_FUNC(sormtr,SORMTR)("L", uplo, "N", n, m, &a[a_offset], lda, &work[indtau]
	      , &z__[z_offset], ldz, &work[indwkn], &llwrkn, &iinfo);
    }

    if (*info != 0) 
      return;

    if (iscale == 1) {
	if (*info == 0) {
	    imax = *m;
	} else {
	    imax = *info - 1;
	}
	d__1 = 1. / sigma;
	F77_FUNC(sscal,SSCAL)(&imax, &d__1, &w[1], &c__1);
    }

    if (wantz) {
	i__1 = *m - 1;
	
	for (j = 1; j <= i__1; ++j) {
	    i__ = 0;
	    tmp1 = w[j];
	    i__2 = *m;
	    for (jj = j + 1; jj <= i__2; ++jj) {
		if (w[jj] < tmp1) {
		    i__ = jj;
		    tmp1 = w[jj];
		}
	    }

	    if (i__ != 0) {
		itmp1 = iwork[indibl + i__ - 1];
		w[i__] = w[j];
		iwork[indibl + i__ - 1] = iwork[indibl + j - 1];
		w[j] = tmp1;
		iwork[indibl + j - 1] = itmp1;
		F77_FUNC(sswap,SSWAP)(n, &z__[i__ * z_dim1 + 1], &c__1, &z__[j * z_dim1 + 1],
			 &c__1);
	    }
	}
    }

    work[1] = (float) lwkopt;
    iwork[1] = liwmin;
    return;

}
