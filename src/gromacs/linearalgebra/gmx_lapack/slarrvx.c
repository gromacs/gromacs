#include <math.h>

#include "gromacs/utility/real.h"

#include "../gmx_blas.h"
#include "../gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(slarrvx,SLARRVX)(int *n, 
	float *d__, 
	float *l, 
	int *isplit,
	int *m, 
	float *w,
	int *iblock, 
	int *indexw, 
	float *gersch, 
	float *tol, 
	float *z__, 
	int *ldz, 
	int *isuppz, 
	float *work, 
	int *iwork, 
	int *info)
{
    int z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    float d__1, d__2;
    float c_b5 = 0.;
    int c__1 = 1;
    int c__2 = 2;

    int i__, j, k, p, q;
    int im, in;
    float gap, eps, tmp;
    int zto;
    float ztz;
    int iend, jblk;
    int wend, iter, temp[1], ktot;
    int itmp1, itmp2;
    int indld;
    float sigma;
    int ndone, iinfo, iindr;
    float resid;
    int nomgs;
    int nclus;
    int zfrom, iindc1, iindc2;
    float lambda;
    int ibegin;
    int indgap, indlld;
    float mingma;
    int oldien, oldncl, wbegin;
    float relgap;
    int oldcls;
    int ndepth, inderr, iindwk;
    int newcls, oldfst;
    float minrgp=0.0;
    int indwrk, oldlst;
    float reltol;
    int newfrs, newftt, parity;
    float mgstol, nrminv, rqcorr;
    int newlst, newsiz;


    --d__;
    --l;
    --isplit;
    --w;
    --iblock;
    --indexw;
    --gersch;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --isuppz;
    --work;
    --iwork;

    inderr = *n;
    indld = *n << 1;
    indlld = *n * 3;
    indgap = *n << 2;
    indwrk = *n * 5 + 1;

    iindr = *n;
    iindc1 = *n << 1;
    iindc2 = *n * 3;
    iindwk = (*n << 2) + 1;

    eps = GMX_FLOAT_EPS;

    i__1 = *n << 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iwork[i__] = 0;
    }
    F77_FUNC(slaset,SLASET)("Full", n, m, &c_b5, &c_b5, &z__[z_offset], ldz);
    mgstol = eps * 100.;

    ibegin = 1;
    wbegin = 1;
    i__1 = iblock[*m];
    for (jblk = 1; jblk <= i__1; ++jblk) {
	iend = isplit[jblk];

	wend = wbegin - 1;
L171:
	if (wend < *m) {
	    if (iblock[wend + 1] == jblk) {
		++wend;
		goto L171;
	    }
	}
	if (wend < wbegin) {
	    ibegin = iend + 1;
	    continue;
	}

	if (ibegin == iend) {
	    z__[ibegin + wbegin * z_dim1] = 1.;
	    isuppz[(wbegin << 1) - 1] = ibegin;
	    isuppz[wbegin * 2] = ibegin;
	    ibegin = iend + 1;
	    wbegin = wend + 1;
	    continue;
	}
	oldien = ibegin - 1;
	in = iend - oldien;
	d__1 = .001, d__2 = 1. / (float) in;
	reltol = (d__1<d__2) ? d__1 : d__2;
	im = wend - wbegin + 1;
	F77_FUNC(scopy,SCOPY)(&im, &w[wbegin], &c__1, &work[1], &c__1);
	i__2 = im - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[inderr + i__] = eps * fabs(work[i__]);
	    work[indgap + i__] = work[i__ + 1] - work[i__];
	}
	work[inderr + im] = eps * fabs(work[im]);
	d__2 = fabs(work[im]);
	work[indgap + im] = (d__2>eps) ? d__2 : eps;
	ndone = 0;

	ndepth = 0;
	parity = 1;
	nclus = 1;
	iwork[iindc1 + 1] = 1;
	iwork[iindc1 + 2] = im;

L40:
	if (ndone < im) {
	    oldncl = nclus;
	    nclus = 0;
	    parity = 1 - parity;
	    if (parity == 0) {
		oldcls = iindc1;
		newcls = iindc2;
	    } else {
		oldcls = iindc2;
		newcls = iindc1;
	    }
	    i__2 = oldncl;
	    for (i__ = 1; i__ <= i__2; ++i__) {

		j = oldcls + (i__ << 1);
		oldfst = iwork[j - 1];
		oldlst = iwork[j];
		if (ndepth > 0) {
		    j = wbegin + oldfst - 1;
		    F77_FUNC(scopy,SCOPY)(&in, &z__[ibegin + j * z_dim1], &c__1, &d__[ibegin]
			    , &c__1);
		    i__3 = in - 1;
		    F77_FUNC(scopy,SCOPY)(&i__3, &z__[ibegin + (j + 1) * z_dim1], &c__1, &l[
			    ibegin], &c__1);
		    F77_FUNC(slaset,SLASET)("Full", &in, &c__2, &c_b5, &c_b5, &z__[ibegin + j 
			    * z_dim1], ldz);
		}
		k = ibegin;
		i__3 = in - 1;
		for (j = 1; j <= i__3; ++j) {
		    tmp = d__[k] * l[k];
		    work[indld + j] = tmp;
		    work[indlld + j] = tmp * l[k];
		    ++k;
		}
		if (ndepth > 0) {

		    p = indexw[wbegin - 1 + oldfst];
		    q = indexw[wbegin - 1 + oldlst];
		    d__1 = eps * 4.;
		    i__3 = p - oldfst;
		    F77_FUNC(slarrbx,SLARRBX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 1], &
			    work[indlld + 1], &p, &q, &reltol, &d__1, &i__3, &
			    work[1], &work[indgap + 1], &work[inderr + 1], &
			    work[indwrk + in], &iwork[iindwk], &iinfo);
		}
		newfrs = oldfst;
		i__3 = oldlst;
		for (j = oldfst; j <= i__3; ++j) {
		    if (j == oldlst || work[indgap + j] >= 
			reltol * fabs(work[j])) {
			newlst = j;
		    } else {

			relgap = work[indgap + j] / fabs(work[j]);
			if (j == newfrs) {
			    minrgp = relgap;
			} else {
			    minrgp = (minrgp<relgap) ? minrgp : relgap;
			}
			continue;
		    }
		    newsiz = newlst - newfrs + 1;
		    newftt = wbegin + newfrs - 1;
		    nomgs = newsiz == 1 || newsiz > 1 || minrgp < mgstol;
		    if (newsiz > 1 && nomgs) {

			F77_FUNC(slarrfx,SLARRFX)(&in, &d__[ibegin], &l[ibegin], &work[indld + 
				1], &work[indlld + 1], &newfrs, &newlst, &
				work[1], &sigma, &z__[ibegin + newftt * 
				z_dim1], &z__[ibegin + (newftt + 1) * z_dim1],
				 &work[indwrk], info);
			if (*info == 0) {
			    tmp = eps * fabs(sigma);
			    i__4 = newlst;
			    for (k = newfrs; k <= i__4; ++k) {
				work[k] -= sigma;
				d__1 = work[indgap + k];
				work[indgap + k] = (d__1>tmp) ? d__1 : tmp;
				work[inderr + k] += tmp;
			    }
			    ++nclus;
			    k = newcls + (nclus << 1);
			    iwork[k - 1] = newfrs;
			    iwork[k] = newlst;
			} else {
			    *info = 0;
			    if (minrgp >= mgstol) {
				nomgs = 0;
			    } else {

				work[indwrk] = d__[ibegin];
				i__4 = in - 1;
				for (k = 1; k <= i__4; ++k) {
				    work[indwrk + k] = d__[ibegin + k] + work[
					    indlld + k];
				}
				i__4 = newsiz;
				for (k = 1; k <= i__4; ++k) {
				    iwork[iindwk + k - 1] = 1;
				}
				i__4 = newlst;
				for (k = newfrs; k <= i__4; ++k) {
				    isuppz[2*(oldien + k) - 1] = 1;
				    isuppz[(oldien + k) * 2] = in;
				}
				temp[0] = in;
				F77_FUNC(sstein,SSTEIN)(&in, &work[indwrk], &work[indld + 1], 
					&newsiz, &work[newfrs], &iwork[iindwk]
					, temp, &z__[ibegin + newftt * z_dim1]
					, ldz, &work[indwrk + in], &iwork[
					iindwk + in], &iwork[iindwk + (in*2)], &iinfo);
				if (iinfo != 0) {
				    *info = 2;
				    return;
				}
				ndone += newsiz;
			    }
			}
		    } else {
			ktot = newftt;
			i__4 = newlst;
			for (k = newfrs; k <= i__4; ++k) {
			    iter = 0;
L90:
			    lambda = work[k];

			    F77_FUNC(slar1vx,SLAR1VX)(&in, &c__1, &in, &lambda, &d__[ibegin], &
				    l[ibegin], &work[indld + 1], &work[indlld 
				    + 1], &w[wbegin + k - 1], &gersch[(oldien 
				    << 1) + 1], &z__[ibegin + ktot * z_dim1], 
				    &ztz, &mingma, &iwork[iindr + ktot], &
				    isuppz[(ktot << 1) - 1], &work[indwrk]);
			    tmp = 1. / ztz;
			    nrminv = sqrt(tmp);
			    resid = fabs(mingma) * nrminv;
			    rqcorr = mingma * tmp;
			    if (k == in) {
				gap = work[indgap + k - 1];
			    } else if (k == 1) {
				gap = work[indgap + k];
			    } else {
				d__1 = work[indgap + k - 1], d__2 = work[
					indgap + k];
				gap = (d__1<d__2) ? d__1 : d__2;
			    }
			    ++iter;
			    if (resid > *tol * gap && fabs(rqcorr) > eps * 4. *
				     fabs(lambda)) {
				work[k] = lambda + rqcorr;
				if (iter < 8) {
				    goto L90;
				}
			    }
			    iwork[ktot] = 1;
			    if (newsiz == 1) {
				++ndone;
			    }
			    zfrom = isuppz[(ktot << 1) - 1];
			    zto = isuppz[ktot * 2];
			    i__5 = zto - zfrom + 1;
			    F77_FUNC(sscal,SSCAL)(&i__5, &nrminv, &z__[ibegin + zfrom - 1 + 
				    ktot * z_dim1], &c__1);
			    ++ktot;
			}
			if (newsiz > 1) {
			    itmp1 = isuppz[(newftt << 1) - 1];
			    itmp2 = isuppz[newftt * 2];
			    ktot = oldien + newlst;
			    i__4 = ktot;
			    for (p = newftt + 1; p <= i__4; ++p) {
				i__5 = p - 1;
				for (q = newftt; q <= i__5; ++q) {
				    tmp = -F77_FUNC(sdot,SDOT)(&in, &z__[ibegin + p * 
					    z_dim1], &c__1, &z__[ibegin + q * 
					    z_dim1], &c__1);
				    F77_FUNC(saxpy,SAXPY)(&in, &tmp, &z__[ibegin + q * 
					    z_dim1], &c__1, &z__[ibegin + p * 
					    z_dim1], &c__1);
				}
				tmp = 1. / F77_FUNC(snrm2,SNRM2)(&in, &z__[ibegin + p * 
					z_dim1], &c__1);
				F77_FUNC(sscal,SSCAL)(&in, &tmp, &z__[ibegin + p * z_dim1], &
					c__1);
				i__5 = itmp1, i__6 = isuppz[(p << 1) - 1];
				itmp1 = (i__5<i__6) ? i__5 : i__6;
				i__5 = itmp2, i__6 = isuppz[p * 2];
				itmp2 = (i__5>i__6) ? i__5 : i__6;
			    }
			    i__4 = ktot;
			    for (p = newftt; p <= i__4; ++p) {
				isuppz[(p << 1) - 1] = itmp1;
				isuppz[p * 2] = itmp2;
			    }
			    ndone += newsiz;
			}
		    }
		    newfrs = j + 1;
		}
	    }
	    ++ndepth;
	    goto L40;
	}
	j = wbegin << 1;
	i__2 = wend;
	for (i__ = wbegin; i__ <= i__2; ++i__) {
	    isuppz[j - 1] += oldien;
	    isuppz[j] += oldien;
	    j += 2;

	}
	ibegin = iend + 1;
	wbegin = wend + 1;
    }

    return;

} 
