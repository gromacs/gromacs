/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_levenmar_c = "$Id$";

#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>


#ifdef DOUBLE
typedef double real;
#else
typedef float real;
#endif

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}



real *vector(int nl,int nh)
{
  real *v;
  
  v=(real *)malloc((unsigned) (nh-nl+1)*sizeof(real));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}

int *ivector(int nl, int nh)
{
  int *v;
  
  v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl;
}

double *dvector(int nl, int nh)
{
  double *v;
	
  v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  return v-nl;
}



real **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	real **m;

	m=(real **) malloc((unsigned) (nrh-nrl+1)*sizeof(real*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(real *) malloc((unsigned) (nch-ncl+1)*sizeof(real));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



real **submatrix(real **a, int oldrl, int oldrh, int oldcl, int oldch,
		 int newrl, int newcl)
{
	int i,j;
	real **m;

	m=(real **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(real*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(real *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_dvector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}



void free_matrix(real **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(real **b, int nrl, int nrh, int ncl, int nch)
{
	free((char*) (b+nrl));
}



real **convert_matrix(real *a, int nrl, int nrh, int ncl, int nch)
{
	int i,j,nrow,ncol;
	real **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (real **) malloc((unsigned) (nrow)*sizeof(real*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(real **b, int nrl, int nrh, int ncl, int nch)
{
	free((char*) (b+nrl));
}

#define SWAP(a,b) {real temp=(a);(a)=(b);(b)=temp;}

void gaussj(real **a, int n, real **b, int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol=0,irow=0,j,k,l,ll;
  real big,dum,pivinv;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
	for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}

#undef SWAP


void covsrt(real **covar, int ma, int lista[], int mfit)
{
  int i,j;
  real swap;
  
  for (j=1;j<ma;j++)
    for (i=j+1;i<=ma;i++) covar[i][j]=0.0;
  for (i=1;i<mfit;i++)
    for (j=i+1;j<=mfit;j++) {
      if (lista[j] > lista[i])
	covar[lista[j]][lista[i]]=covar[i][j];
      else
	covar[lista[i]][lista[j]]=covar[i][j];
    }
  swap=covar[1][1];
  for (j=1;j<=ma;j++) {
    covar[1][j]=covar[j][j];
    covar[j][j]=0.0;
  }
  covar[lista[1]][lista[1]]=swap;
  for (j=2;j<=mfit;j++) covar[lista[j]][lista[j]]=covar[1][j];
  for (j=2;j<=ma;j++)
    for (i=1;i<=j-1;i++) covar[i][j]=covar[j][i];
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void covsrt_new(real **covar,int ma, int ia[], int mfit)
     /* Expand in storage the covariance matrix covar, so as to take 
      * into account parameters that are being held fixed. (For the
      * latter, return zero covariances.)
      */
{
  int i,j,k;
  real swap;
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
      for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
      k--;
    }
  }
}
#undef SWAP
	
void mrqcof(real x[], real y[], real sig[], int ndata, real a[], 
	    int ma, int lista[], int mfit, 
	    real **alpha, real beta[], real *chisq,
	    void (*funcs)(real,real *,real *,real *,int)) 
{
  int k,j,i;
  real ymod,wt,sig2i,dy,*dyda;
  
  dyda=vector(1,ma);
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=j;k++) alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  *chisq=0.0;
  for (i=1;i<=ndata;i++) {
    (*funcs)(x[i],a,&ymod,dyda,ma);
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for (j=1;j<=mfit;j++) {
      wt=dyda[lista[j]]*sig2i;
      for (k=1;k<=j;k++)
	alpha[j][k] += wt*dyda[lista[k]];
      beta[j] += dy*wt;
    }
    (*chisq) += dy*dy*sig2i;
  }
  for (j=2;j<=mfit;j++)
    for (k=1;k<=j-1;k++) alpha[k][j]=alpha[j][k];
  free_vector(dyda,1,ma);
}

	
void mrqmin(real x[], real y[], real sig[], int ndata, real a[], 
	    int ma, int lista[], int mfit, 
	    real **covar, real **alpha, real *chisq,
	    void (*funcs)(real,real *,real *,real *,int),
	    real *alamda) 
{
  int k,kk,j,ihit;
  static real *da,*atry,**oneda,*beta,ochisq;
  
  if (*alamda < 0.0) {
    oneda=matrix(1,mfit,1,1);
    atry=vector(1,ma);
    da=vector(1,ma);
    beta=vector(1,ma);
    kk=mfit+1;
    for (j=1;j<=ma;j++) {
      ihit=0;
      for (k=1;k<=mfit;k++)
	if (lista[k] == j) ihit++;
      if (ihit == 0)
	lista[kk++]=j;
      else if (ihit > 1) nrerror("Bad LISTA permutation in MRQMIN-1");
    }
    if (kk != ma+1) nrerror("Bad LISTA permutation in MRQMIN-2");
    *alamda=0.001;
    mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq,funcs);
    ochisq=(*chisq);
  }
  for (j=1;j<=mfit;j++) {
    for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    oneda[j][1]=beta[j];
  }
  gaussj(covar,mfit,oneda,1);
  for (j=1;j<=mfit;j++)
    da[j]=oneda[j][1];
  if (*alamda == 0.0) {
    covsrt(covar,ma,lista,mfit);
    free_vector(beta,1,ma);
    free_vector(da,1,ma);
    free_vector(atry,1,ma);
    free_matrix(oneda,1,mfit,1,1);
    return;
  }
  for (j=1;j<=ma;j++) atry[j]=a[j];
  for (j=1;j<=mfit;j++)
    atry[lista[j]] = a[lista[j]]+da[j];
  mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,chisq,funcs);
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=1;j<=mfit;j++) {
      for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
      beta[j]=da[j];
      a[lista[j]]=atry[lista[j]];
    }
  } else {
    *alamda *= 10.0;
    *chisq=ochisq;
  }
  return;
}


void mrqmin_new(real x[],real y[],real sig[],int ndata,real a[], 
		int ia[],int ma,real **covar,real **alpha,real *chisq, 
		void (*funcs)(real, real [], real *, real [], int), 
		real *alamda)
     /* Levenberg-Marquardt method, attempting to reduce the value Chi^2 
      * of a fit between a set of data points x[1..ndata], y[1..ndata]
      * with individual standard deviations sig[1..ndata], and a nonlinear 
      * function dependent on ma coefficients a[1..ma]. The input array
      * ia[1..ma] indicates by nonzero entries those components of a that
      * should be fitted for, and by zero entries those components that 
      * should be held fixed at their input values. The program returns 
      * current best-fit values for the parameters a[1..ma], and 
      * Chi^2 = chisq. The arrays covar[1..ma][1..ma], alpha[1..ma][1..ma]
      * are used as working space during most iterations. Supply a routine
      * funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, 
      * and its derivatives dyda[1..ma] with respect to the fitting 
      * parameters a at x. On the first call provide an initial guess for 
      * the parameters a, and set alamda < 0 for initialization (which then 
      * sets alamda=.001). If a step succeeds chisq becomes smaller and 
      * alamda de-creases by a factor of 10. If a step fails alamda grows by 
      * a factor of 10. You must call this routine repeatedly until 
      * convergence is achieved. Then, make one final call with alamda=0, 
      * so that covar[1..ma][1..ma] returns the covariance matrix, and alpha 
      * the curvature matrix.
      * (Parameters held fixed will return zero covariances.)
      */
{
  void covsrt(real **covar, int ma, int ia[], int mfit);
  void gaussj(real **a, int n, real **b,int m);
  void mrqcof_new(real x[], real y[], real sig[], int ndata, real a[],
	      int ia[], int ma, real **alpha, real beta[], real *chisq,
	      void (*funcs)(real, real [], real *, real [], int));
  int j,k,l;
  static int mfit;
  static real ochisq,*atry,*beta,*da,**oneda;

  if (*alamda < 0.0) {                    /* Initialization. */
    atry=vector(1,ma);
    beta=vector(1,ma);
    da=vector(1,ma);
    for (mfit=0,j=1;j<=ma;j++)
      if (ia[j]) mfit++;
    oneda=matrix(1,mfit,1,1);
    *alamda=0.001;
    mrqcof_new(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
    ochisq=(*chisq);
    for (j=1;j<=ma;j++)
      atry[j]=a[j];
  }
  for (j=1;j<=mfit;j++) { /* Alter linearized fitting matrix, by augmenting. */
    for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k]; /* diagonal elements. */
    covar[j][j]=alpha[j][j]*(1.0+(*alamda));
    oneda[j][1]=beta[j];
  }
  gaussj(covar,mfit,oneda,1);      /* Matrix solution. */
  for (j=1;j<=mfit;j++) 
    da[j]=oneda[j][1];
  if (*alamda == 0.0) { /* Once converged, evaluate covariance matrix. */
    covsrt_new(covar,ma,ia,mfit);
    free_matrix(oneda,1,mfit,1,1);
    free_vector(da,1,ma);
    free_vector(beta,1,ma);
    free_vector(atry,1,ma);
    return;
  }
  for (j=0,l=1;l<=ma;l++) /* Did the trial succeed? */
    if (ia[l]) atry[l]=a[l]+da[++j];
  mrqcof_new(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
  if (*chisq < ochisq) {
    /* Success, accept the new solution. */
    *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=1;j<=mfit;j++) {
	for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
	beta[j]=da[j];
    }
    for (l=1;l<=ma;l++) a[l]=atry[l];
  } else {               /* Failure, increase alamda and return. */
    *alamda *= 10.0;
    *chisq=ochisq;
  }
}

void mrqcof_new(real x[], real y[], real sig[], int ndata, real a[], 
	    int ia[], int ma, real **alpha, real beta[], real *chisq,
	    void (*funcs)(real, real [], real *, real[], int))
     /* Used by mrqmin to evaluate the linearized fitting matrix alpha, and 
      * vector beta as in (15.5.8), and calculate Chi^2.
      */
{
  int i,j,k,l,m,mfit=0;
  real ymod,wt,sig2i,dy,*dyda;

  dyda=vector(1,ma);
  for (j=1;j<=ma;j++)
    if (ia[j]) mfit++;
  for (j=1;j<=mfit;j++) { /* Initialize (symmetric) alpha), beta. */
    for (k=1;k<=j;k++) alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  *chisq=0.0;
  for (i=1;i<=ndata;i++) { /* Summation loop over all data. */
    (*funcs)(x[i],a,&ymod,dyda,ma);
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	wt=dyda[l]*sig2i;
	for (j++,k=0,m=1;m<=l;m++)
	  if (ia[m]) alpha[j][++k] += wt*dyda[m];
	beta[j] += dy*wt;
      }
    }
    *chisq += dy*dy*sig2i;  /* And find Chi^2. */
  }
  for (j=2;j<=mfit;j++)     /* Fill in the symmetric side. */
    for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
  free_vector(dyda,1,ma);
}
