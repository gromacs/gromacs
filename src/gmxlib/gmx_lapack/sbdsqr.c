#include <ctype.h>
#include <math.h>
#include "gmx_blas.h"
#include "gmx_lapack.h"
#include "lapack_limits.h"


void
F77_FUNC(sbdsqr,SBDSQR)(char *uplo, 
	int *n, 
	int *ncvt, 
	int *nru, 
	int *ncc, 
	float *d, 
	float *e,
	float *vt, 
	int *ldvt, 
	float *u, 
	int *ldu,
	float *c, 
	int *ldc,
	float *work,
	int *info)
{
  char xuplo=toupper(*uplo);
  float f, g, h;
  int i, j, m, i1;
  float r, cs;
  int ll;
  float sn, mu;
  int nm1, nm12, nm13, lll;
  float sll, tol, abse,minval;
  int idir;
  float abss;
  int oldm;
  float cosl;
  int isub, iter;
  float unfl, sinl, cosr, smin, smax, sinr;
  int n1 = *n;
  float oldcs;
  int oldll;
  float shift, sigmn, oldsn;
  int maxit;
  float sminl, sigmx;
  float sminoa, thresh;
  int rotate;
  float sminlo, tolmul,d1,d2;
  const float tolmul1 = pow(LAPACK_EPS_FLOAT,-0.125);
  float minusone = -1.0;

  tolmul = (100.0 < tolmul1) ? 100.0 : tolmul1;
  tolmul = (10.0 > tolmul) ? 10.0 : tolmul;
  tol    = tolmul * LAPACK_EPS_FLOAT;

  *info = 0;

  if (n1<=0) 
    return;
  
  if (n1==1)
    goto L160;

  rotate =(*ncvt > 0 || *nru > 0 || *ncc > 0);
  
  if (! rotate) {
    F77_FUNC(slasq1,SLASQ1)(n,d,e,work,info);
    return;
  }
  
  nm1 = n1 - 1;
  nm12 = nm1 + nm1;
  nm13 = nm12 + nm1;
  idir = 0;
  
  minval = (1.0 + LAPACK_EPS_FLOAT)/LAPACK_MAX_FLOAT;
  unfl = minval / LAPACK_EPS_FLOAT;
  
  if (xuplo=='L') {
    for (i=0;i<(n1-1);i++) {
      F77_FUNC(slartg,SLARTG)(d+i, e+i, &cs, &sn, &r);
      d[i] = r;
      e[i] = sn * d[i + 1];
      d[i + 1] *= cs;
      work[i] = cs;
      work[nm1 + i] = sn;
    }
    
    if (*nru > 0) 
      F77_FUNC(slasr,SLASR)("R","V","F",nru,n,work,work+n1-1,u,ldu);
    if (*ncc > 0) 
      F77_FUNC(slasr,SLASR)("L","V","F",n,ncc,work,work+n1-1,c,ldc);
  }
  
  smax = 0.0;
  for (i=0;i<n1;i++) {
    d1 = fabs(d[i]);
    if(d1>smax)
      smax = d1;
  }
  for (i=0; i<(n1-1);i++) {
    d1 = fabs(e[i]);
    if(d1>smax)
      smax = d1;
  }
  sminl = 0.0;
  if (tol >= 0.0) {
    sminoa = fabs(d[0]);
    if (sminoa!= 0.0) {
      mu = sminoa;
      for (i=1;i<n1 && sminoa!=0.0;i++) {
	mu *= fabs(d[i])/(mu+fabs(e[i-1]));
	if(mu<sminoa)
	  sminoa = mu;
      }
    }
    sminoa /= sqrt((float)n1);
    d1 = tol * sminoa;
    d2 = DBDSQR_MAXITR * n1 * n1 * unfl;
    thresh = (d1>d2) ? d1 : d2;
  } else {
    d1 = fabs(tol) * smax;
    d2 = DBDSQR_MAXITR * n1 * n1 * unfl;
    thresh = (d1>d2) ? d1 : d2;
  }
  maxit = DBDSQR_MAXITR * n1 * n1;
  iter = 0;
  oldll = -1;
  oldm = -1;
  m = n1;

    /* Start of main iteration loop */
 L60:

  if(m<=1)
    goto L160;
  if(iter>maxit) {
    *info = 0;
    for (i=0;i<n1-1;i++) {
      if (e[i] != 0.0) 
	++(*info);
    }
    return;
  }
 
  if (tol < 0.0 && fabs(d[m-1])<=thresh) 
    d[m-1] = 0.0;
  
  smax = fabs(d[m-1]);
  smin = smax;
  for (lll = 1; lll<=(m-1);lll++) {
    ll = m - lll;
    abss = fabs(d[ll-1]);
    abse = fabs(e[ll-1]);
    if (tol < 0.0 && abss <= thresh) 
      d[ll-1] = 0.0;
    if (abse <= thresh) 
      goto L80;
    
    if(abss<smin)
      smin = abss;
    if(abss>smax)
      smax=abss;
    if(abse>smax)
      smax=abse;
  }
  
  ll = 0;
  goto L90;
  
 L80: 
  e[ll-1] = 0.0;
  if (ll == m-1) {
    m--;;
    goto L60;
  }

L90:
    ll++;

    if (ll == m-1) {
	F77_FUNC(slasv2,SLASV2)(d+m-2,e+m-2,d+m-1,&sigmn,&sigmx,&sinr,&cosr,&sinl,&cosl);
	d[m-2] = sigmx;
	e[m-2] = 0.0;
	d[m-1] = sigmn;
	i1 = 1;
	if (*ncvt > 0) 
	  F77_FUNC(srot,SROT)(ncvt,vt+m-2,ldvt,vt+m-1,ldvt,&cosr,&sinr);
	if (*nru > 0) 
	  F77_FUNC(srot,SROT)(nru,u+(m-2)*(*ldu),&i1,u+(m-1)*(*ldu),&i1,&cosl,&sinl);	
	if (*ncc > 0) 
	  F77_FUNC(srot,SROT)(ncc,c+m-2,ldc,c+m-1,ldc,&cosl,&sinl);
	m -= 2;
	goto L60;
    }

    if (ll> oldm || m < oldll) {
      if(fabs(d[ll-1])>=fabs(d[m-1]))
	idir = 1;
      else
	idir = 2;
    }
    if (idir == 1) {
      if((fabs(e[m-2])<=fabs(tol*d[m-1])) ||
	 (tol<0.0 && fabs(e[m-2])<thresh)) {
	e[m - 2] = 0.0;
	goto L60;
      }

      if (tol >= 0.0) {
	mu = fabs(d[ll-1]);
	sminl = mu;
	for (lll = ll; lll<m;lll++) {
	  if(fabs(e[lll-1])<=(tol*mu)) {
	    e[lll-1] = 0.0;
	    goto L60;
	    sminlo = sminl;
	    mu *= fabs(d[lll])/(mu+fabs(e[lll-1]));
	    if(mu<sminl)
	      sminl=mu;
	  }
	}
      }
    } else {
      if((fabs(e[ll-1])<=fabs(tol*d[ll-1])) ||
	 (tol<0.0 && fabs(e[ll-1])<thresh)) {
	e[ll-1] = 0.0;
	goto L60;
      }

      if (tol >= 0.0) {
	mu = fabs(d[m-1]);
	sminl = mu;
	for (lll = m-1; lll>=ll;lll--) {
	  if(fabs(e[lll-1])<=(tol*mu)) {
	    e[lll-1] = 0.0;
	    goto L60;
	    sminlo = sminl;
	    mu *= fabs(d[lll-1])/(mu+fabs(e[lll-1]));
	    if(mu<sminl)
	      sminl=mu;
	  }
	}
      }
    }
    oldll = ll;
    oldm = m;

    d1 = 0.01*tol;
    if(LAPACK_EPS_FLOAT>d1)
      d1 = LAPACK_EPS_FLOAT;

    if(tol>=0.0 && n1*tol*(sminl/smax)<=d1) {
      shift = 0.0;
    } else {
	if (idir == 1) {
	    sll = fabs(d[ll-1]);
	    F77_FUNC(slas2,SLAS2)(d+m-2,e+m-2,d+m-1, &shift, &r);
	} else {
	    sll = fabs(d[m-1]);
	    F77_FUNC(slas2,SLAS2)(d+ll-1,e+ll-1,d+ll, &shift, &r);
	}
	if (sll > 0.0) {
	  d1 = shift / sll;
	  if (d1 * d1 < LAPACK_EPS_FLOAT) 
	    shift = 0.0;
	}
    }
    iter += m - ll;
    if (shift == 0.0) {
	if (idir == 1) {

	    cs = 1.0;
	    oldcs = 1.0;
	    for (i = ll-1; i<(m-1);i++) {
		d1 = d[i] * cs;
		F77_FUNC(slartg,SLARTG)(&d1,e+i, &cs, &sn, &r);
		if ((i+1) > ll) 
		    e[i - 1] = oldsn * r;
		d1 = oldcs * r;
		d2 = d[i + 1] * sn;
		F77_FUNC(slartg,SLARTG)(&d1, &d2, &oldcs, &oldsn,d+i);
		work[i - ll + 1] = cs;
		work[i - ll + 1 + nm1] = sn;
		work[i - ll + 1 + nm12] = oldcs;
		work[i - ll + 1 + nm13] = oldsn;
	    }
	    h = d[m-1] * cs;
	    d[m-1] = h * oldcs;
	    e[m - 2] = h * oldsn;

	    i1 = m-ll+1;
	    if (*ncvt > 0) 
	      F77_FUNC(slasr,SLASR)("L","V","F",&i1,ncvt,work,work+n1-1, vt+ll-1,ldvt);
	    
	    if (*nru > 0) 
	      F77_FUNC(slasr,SLASR)("R","V","F",nru,&i1,work+nm12,work+nm13,u+(ll-1)*(*ldu),ldu);
	    
	    if (*ncc > 0) 
	      F77_FUNC(slasr,SLASR)("L","V","F",&i1,ncc,work+nm12,work+nm13,c+ll-1,ldc);

	    if(fabs(e[m-2])<=thresh)
		e[m - 2] = 0.0;
	} else {

	    cs = 1.0;
	    oldcs = 1.0;
	    for (i = m-1; i>=ll;i--) {
		d1 = d[i] * cs;
		F77_FUNC(slartg,SLARTG)(&d1,e+i-1, &cs, &sn, &r);
		if ((i+1) < m) 
		    e[i] = oldsn * r;
		d1 = oldcs * r;
		d2 = d[i-1] * sn;
		F77_FUNC(slartg,SLARTG)(&d1, &d2, &oldcs, &oldsn,d+i);
		work[i - ll] = cs;
		work[i - ll + nm1] = -sn;
		work[i - ll + nm12] = oldcs;
		work[i - ll + nm13] = -oldsn;
	    }
	    h = d[ll-1] * cs;
	    d[ll-1] = h * oldcs;
	    e[ll-1] = h * oldsn;

	    i1 = m-ll+1;
	    if (*ncvt > 0) 
	      F77_FUNC(slasr,SLASR)("L","V","B",&i1,ncvt,work+nm12,work+nm13, vt+ll-1,ldvt);
	    
	    if (*nru > 0) 
	      F77_FUNC(slasr,SLASR)("R","V","B",nru,&i1,work,work+n1-1,u+(ll-1)*(*ldu),ldu);
	    
	    if (*ncc > 0) 
	      F77_FUNC(slasr,SLASR)("L","V","B",&i1,ncc,work,work+n1-1,c+ll-1,ldc);

	    if(fabs(e[ll-1])<=thresh)
		e[ll-1] = 0.0;
	}
    } else {

      /* nonzero shift */
	if (idir == 1) {
	  d1 = (d[ll-1]>0) ? 1.0 : -1.0;
	  f = ( fabs(d[ll-1])-shift)*(d1+shift/d[ll-1]);
	  g = e[ll-1];
	  for (i = ll-1; i<m-1;++i) {
	    F77_FUNC(slartg,SLARTG)(&f, &g, &cosr, &sinr, &r);
	    if ((i+1) > ll) 
	      e[i - 1] = r;
	    
	    f = cosr * d[i] + sinr * e[i];
	    e[i] = cosr * e[i] - sinr * d[i];
	    g = sinr * d[i + 1];
	    d[i + 1] = cosr * d[i + 1];
	    F77_FUNC(slartg,SLARTG)(&f, &g, &cosl, &sinl, &r);
	    d[i] = r;
	    f = cosl * e[i] + sinl * d[i + 1];
	    d[i + 1] = cosl * d[i + 1] - sinl * e[i];
	    if ((i+1) < m - 1) {
	      g = sinl * e[i + 1];
	      e[i + 1] = cosl * e[i + 1];
	    }
	    work[i - ll + 1] = cosr;
	    work[i - ll + 1 + nm1] = sinr;
	    work[i - ll + 1 + nm12] = cosl;
	    work[i - ll + 1 + nm13] = sinl;
	  }
	  e[m - 2] = f;

	  i1 = m - ll + 1;
	  if (*ncvt > 0)
	    F77_FUNC(slasr,SLASR)("L","V","F",&i1, ncvt,work,work+n1-1,vt+ll-1,ldvt);
	  
	  if (*nru > 0) 
	    F77_FUNC(slasr,SLASR)("R","V","F",nru, &i1,work+nm12,work+nm13,u+(ll-1)*(*ldu),ldu);
	  
	  if (*ncc > 0) 
	    F77_FUNC(slasr,SLASR)("L","V","F",&i1, ncc, work+nm12, work+nm13,c+ll-1,ldc);
	  
	  if(fabs(e[m-2])<=thresh)
		e[m - 2] = 0.0;

	} else {

	  d1 = (d[m-1]>0) ? 1.0 : -1.0;
	  f = ( fabs(d[m-1])-shift)*(d1+shift/d[m-1]);
	  g = e[m-2];
	  for (i = m-1; i>=ll;i--) {
	    F77_FUNC(slartg,SLARTG)(&f, &g, &cosr, &sinr, &r);
	    if ((i+1) < m)
	      e[i] = r;
	    
	    f = cosr * d[i] + sinr * e[i-1];
	    e[i-1] = cosr * e[i-1] - sinr * d[i];
	    g = sinr * d[i - 1];
	    d[i - 1] = cosr * d[i - 1];
	    F77_FUNC(slartg,SLARTG)(&f, &g, &cosl, &sinl, &r);
	    d[i] = r;
	    f = cosl * e[i-1] + sinl * d[i - 1];
	    d[i - 1] = cosl * d[i - 1] - sinl * e[i-1];
	    if (i > ll) {
	      g = sinl * e[i - 2];
	      e[i - 2] = cosl * e[i - 2];
	    }
	    work[i - ll] = cosr;
	    work[i - ll + nm1] = sinr;
	    work[i - ll + nm12] = cosl;
	    work[i - ll + nm13] = sinl;
	  }
	  e[ll-1] = f;


	  if(fabs(e[ll-1])<=thresh)
		e[ll-1] = 0.0;

	  i1 = m - ll + 1;
	  if (*ncvt > 0)
	    F77_FUNC(slasr,SLASR)("L","V","B",&i1, ncvt,work+nm12,work+nm13,vt+ll-1,ldvt);
	  
	  if (*nru > 0) 
	    F77_FUNC(slasr,SLASR)("R","V","B",nru, &i1,work,work+n1-1,u+(ll-1)*(*ldu),ldu);
	  
	  if (*ncc > 0) 
	    F77_FUNC(slasr,SLASR)("L","V","B",&i1, ncc, work, work+n1-1,c+ll-1,ldc);
	  
	}
    }

    goto L60;

L160:
    for (i=0; i<n1;i++) {
	if (d[i] < 0.0) {
	    d[i] = -d[i];

	    if (*ncvt > 0) 
	      F77_FUNC(sscal,SSCAL)(ncvt, &minusone, vt+i, ldvt);
	    
	}
    }
    for (i=0; i<n1-1;i++) {
	isub = 0;
	smin = d[0];
	for (j=1; j<n1-i;j++) {
	    if (d[j] <= smin) {
		isub = j;
		smin = d[j];
	    }
	}
	if (isub != n1 - i - 1) {
	  i1 = 1;
	  d[isub] = d[n1 - i -1];
	  d[n1 - i - 1] = smin;
	  if (*ncvt > 0)
	    F77_FUNC(sswap,SSWAP)(ncvt,vt+isub,ldvt,vt+n1-i-1,ldvt);
	  
	  if (*nru > 0) 
	    F77_FUNC(sswap,SSWAP)(nru,u+isub*(*ldu),&i1,u+(n1-i-1)*(*ldu),&i1);
	  
	  if (*ncc > 0) 
	    F77_FUNC(sswap,SSWAP)(ncc,c+isub,ldc,c+n1-i-1,ldc);
	}
    }
    return;
}

