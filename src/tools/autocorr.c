/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "typedefs.h"
#include "physics.h"
#include "smalloc.h"
#include "xvgr.h"
#include "futil.h"
#include "gstat.h"
#include "names.h"
#include "fatal.h"
#include "vec.h"
#include "string2.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define MODE(x) ((mode & (x)) == (x))

typedef struct {
  unsigned long mode;
  int  nrestart,nout,P,fitfn,nskip;
  bool bFour,bNormalize;
  real tbeginfit,tendfit;
} t_acf;

static bool  bACFinit = FALSE;
static t_acf acf;

typedef real fftreal;

int sffn2effn(char **sffn)
{
  int eFitFn,i;
  
  eFitFn = 0;
  for(i=0; i<effnNR; i++)
    if (sffn[i+1] && strcmp(sffn[0],sffn[i+1])==0)
      eFitFn = i;

  return eFitFn;
}

void four1(fftreal data[],int nn,int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  fftreal tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

#undef SWAP

void realft(fftreal data[],int n,int isign)
{
  int i,i1,i2,i3,i4,n2p3;
  fftreal c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;

  theta=3.141592653589793/(double) n;
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  n2p3=2*n+3;
  for (i=2;i<=n/2;i++) {
    i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n,-1);
  }
}

void twofft(fftreal data1[],fftreal data2[],fftreal fft1[],fftreal fft2[],int n)
{
  int nn3,nn2,jj,j;
  fftreal rep,rem,aip,aim;

  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  four1(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2) {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}

enum { enNorm, enCos, enSin };

void correl(fftreal data1[],fftreal data2[],int n,fftreal ans[])
{
  int     no2,i;
  fftreal dum,*fft;

  snew(fft,2*n+1);
  twofft(data1,data2,fft,ans,n);
  no2=n/2;
  for (i=2;i<=n+2;i+=2) {
    dum      = ans[i-1];
    ans[i-1] = (fft[i-1]*dum+fft[i]*ans[i])/no2;
    ans[i]   = (fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  realft(ans,no2,-1);
  sfree(fft);
}

static void low_do_four_core(int nfour,int nframes,real c1[],fftreal cfour[],
			     int nCos,bool bPadding)
{
  int  i=0;
  fftreal aver,*ans;

  aver = 0.0;
  switch (nCos) {
  case enNorm:  
    for(i=0; (i<nframes); i++) {
      aver+=c1[i];
      cfour[i]=c1[i];
    }
    break;
  case enCos:
    for(i=0; (i<nframes); i++) 
      cfour[i]=cos(c1[i]);
    break;
  case enSin:
    for(i=0; (i<nframes); i++) 
      cfour[i]=sin(c1[i]);
    break;
  default:
    fatal_error(0,"nCos = %d, %s %d",nCos,__FILE__,__LINE__);
  }
  for(   ; (i<nfour); i++)
    cfour[i]= 0.0;
  
  if (bPadding) {
    aver /= nframes;
    /* printf("aver = %g\n",aver); */
    for(i=0; (i<nframes); i++)
      cfour[i] -= aver;
  }
    
  snew(ans,2*nfour);
  correl(cfour-1,cfour-1,nfour,ans-1);
  
  if (bPadding)  
    for (i=0; (i<nfour); i++)
      cfour[i] = ans[i]+sqr(aver);
  else
    for (i=0; (i<nfour); i++)
      cfour[i] = ans[i];
      
  sfree(ans);
}

static void do_ac_core(int nframes,int nout,
		       real corr[],real c1[],int nrestart,
		       unsigned long mode)
{
  int     j,k,j3,jk3,m,n;
  fftreal ccc,c0,cth;
  rvec    xj,xk,rr;

  if (nrestart < 1) {
    printf("WARNING: setting number of restarts to 1\n");
    nrestart = 1;
  }
  if (debug)
    fprintf(debug,
	    "Starting do_ac_core: nframes=%d, nout=%d, nrestart=%d,mode=%lu\n",
	    nframes,nout,nrestart,mode);
  
  for(j=0; (j<nout); j++)
    corr[j]=0;
  
  /* Loop over starting points. */
  for(j=0; (j<nframes); j+=nrestart) {
    j3  = DIM*j;
    
    /* Loop over the correlation length for this starting point */
    for(k=0; (k<nout) && (j+k < nframes); k++) {
      jk3 = DIM*(j+k);
      
      /* Switch over possible ACF types. 
       * It might be more efficient to put the loops inside the switch,
       * but this is more clear, and save development time!
       */      
      if (MODE(eacNormal)) {
	corr[k] += c1[j]*c1[j+k];
      }
      else if (MODE(eacCos)) {
	/* Compute the cos (phi(t)-phi(t+dt)) */
	corr[k] += cos(c1[j]-c1[j+k]);
      }
      else if (MODE(eacP1) || MODE(eacP2) || MODE(eacP3)) {
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	cth=cos_angle(xj,xk);
	
	if (cth-1.0 > 1.0e-15) {
	  printf("j: %d, k: %d, xj:(%g,%g,%g), xk:(%g,%g,%g)\n",
		  j,k,xj[XX],xj[YY],xj[ZZ],xk[XX],xk[YY],xk[ZZ]);
	}
	
	corr[k] += LegendreP(cth,mode);  /* 1.5*cth*cth-0.5; */
      }
      else if (MODE(eacRcross)) {
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	oprod(xj,xk,rr);
	
	corr[k] += iprod(rr,rr);
      }
      else if (MODE(eacVector)) {
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	ccc = iprod(xj,xk);
	
	corr[k] += ccc;
      }
      else
	fatal_error(0,"\nInvalid mode (%d) in do_ac_core",mode);
    }
  }
  /* Correct for the number of points and copy results to the data array */
  for(j=0; (j<nout); j++) {
    n = (nframes-j+(nrestart-1))/nrestart;
    c1[j] = corr[j]/n;
  }
}

void normalize_acf(int nout,real corr[])
{
  int  j;
  real c0;

  if (debug) {
    fprintf(debug,"Before normalization\n");
    for(j=0; (j<nout); j++) 
      fprintf(debug,"%5d  %10f\n",j,corr[j]);
  }

  /* Normalisation makes that c[0] = 1.0 and that other points are scaled
   * accordingly.
   */
  if (corr[0] == 0.0)
    c0 = 1.0;
  else
    c0 = 1.0/corr[0];
  for(j=0; (j<nout); j++)
    corr[j] *= c0;

  if (debug) {
    fprintf(debug,"After normalization\n");
    for(j=0; (j<nout); j++) 
      fprintf(debug,"%5d  %10f\n",j,corr[j]);
  }
}

void average_acf(bool bVerbose,int n,int nitem,real **c1)
{
  real c0;
  int  i,j;
  
  if (bVerbose)
    printf("Averaging correlation functions\n");
  
  for(j=0; (j<n); j++) {
    c0 = 0;
    for(i=0; (i<nitem); i++)
      c0+=c1[i][j];
    c1[0][j] = c0/nitem;
  }
}

void norm_and_scale_vectors(int nframes,real c1[],real scale)
{
  int  j,m;
  real *rij;
  
  for(j=0; (j<nframes); j++) {
    rij = &(c1[j*DIM]);
    unitv(rij,rij);
    for(m=0; (m<DIM); m++)
      rij[m]*=scale;
  }
}

void dump_tmp(char *s,int n,real c[])
{
  FILE *fp;
  int  i;
  
  fp=ffopen(s,"w");
  for(i=0; (i<n); i++)
    fprintf(fp,"%10d  %10g\n",i,c[i]);
  ffclose(fp);
}

real print_and_integrate(FILE *fp,int n,real dt,real c[],real *fit,int nskip)
{
  real c0,sum;
  int  j;
  
  /* Use trapezoidal rule for calculating integral */
  sum = 0.0;
  for(j=0; (j<n); j++) {
    c0 = c[j];
    if (fp && (nskip == 0 || j % nskip == 0))
      fprintf(fp,"%10.3f  %10.5f\n",j*dt,c0);
    if (j > 0)
      sum+=dt*(c0+c[j-1]);
  }
  if (fp) {
    fprintf(fp,"&\n");
    if (fit) {
      for(j=0; (j<n); j++)
	if (nskip == 0 || j % nskip == 0)
	  fprintf(fp,"%10.3f  %10.5f\n",j*dt,fit[j]);
      fprintf(fp,"&\n");
    }
  }
  return sum*0.5;
}

real evaluate_integral(int n,real dx,real y[],real dy[],real aver_start,
		       real *stddev)
{
  double c0,sum,dsum=0,dsum2,sss;
  int    j,ndsum=0;
  
  /* Use trapezoidal rule for calculating integral */
  if (n <= 0)
    fatal_error(0,"Evaluating integral: n = %d (file %s, line %d)",
		n,__FILE__,__LINE__);
  sum  = y[0]+y[n-1];
  if (dy)
    dsum2 = sqr(dy[0]) + sqr(dy[n-1]);
  else
    dsum2 = 0;
  for(j=1; (j<n-1); j++) {
    sum += 2*y[j];
    if (j*dx >= aver_start) {
      sss    = dx*sum*0.5;
      dsum  += sss;
      if (dy) 
	dsum2 += sqr(dy[j]);
      else
	dsum2 += sss*sss;
      ndsum++;
    }
  }
  if (ndsum > 1) {
    dsum2 /= ndsum;
    dsum  /= ndsum;
    *stddev = sqrt(dsum2-dsum*dsum);
  }
  else {
    *stddev = 0.0;
    dsum = sum;
  }
  /* return sum*0.5*dx; */
  return dsum;
}

void do_four_core(unsigned long mode,int nfour,int nf2,int nframes,
		  real c1[],real csum[],real ctmp[])
{
  fftreal *cfour;
  char    buf[32];
  real    fac;
  int     i,j,m,m1;
  
  snew(cfour,nfour);
  
  if (MODE(eacNormal)) {
    /********************************************
     *  N O R M A L
     ********************************************/
    low_do_four_core(nfour,nf2,c1,csum,enNorm,FALSE);
  }
  else if (MODE(eacCos)) {
    /***************************************************
     * C O S I N E
     ***************************************************/
    /* Copy the data to temp array. Since we need it twice
     * we can't overwrite original.
     */
    for(j=0; (j<nf2); j++)
      ctmp[j]=c1[j];
    
    /* Cosine term of AC function */
    low_do_four_core(nfour,nf2,ctmp,cfour,enCos,FALSE);
    for(j=0; (j<nf2); j++)
      c1[j]  = cfour[j];
    
    /* Sine term of AC function */
    low_do_four_core(nfour,nf2,ctmp,cfour,enSin,FALSE);
    for(j=0; (j<nf2); j++) {
      c1[j] += cfour[j];
      csum[j] = c1[j];
    }
  }
  else if (MODE(eacP2)) {
    /***************************************************
     * Legendre polynomials
     ***************************************************/
    /* First normalize the vectors */
    norm_and_scale_vectors(nframes,c1,1.0);
    
    /* For P2 thingies we have to do six FFT based correls 
     * First for XX^2, then for YY^2, then for ZZ^2
     * Then we have to do XY, YZ and XZ (counting these twice)
     * After that we sum them and normalise
     * P2(x) = (3 * cos^2 (x) - 1)/2
     * for unit vectors u and v we compute the cosine as the inner product
     * cos(u,v) = uX vX + uY vY + uZ vZ
     *
     *        oo
     *        /
     * C(t) = |  (3 cos^2(u(t'),u(t'+t)) - 1)/2 dt'
     *        /
     *        0
     *
     * For ACF we need:
     * P2(u(0),u(t)) = [3 * (uX(0) uX(t) + 
     *                       uY(0) uY(t) + 
     *                       uZ(0) uZ(t))^2 - 1]/2 
     *               = [3 * ((uX(0) uX(t))^2 +
     *                       (uY(0) uY(t))^2 +
     *                       (uZ(0) uZ(t))^2 +
     *                 2(uX(0) uY(0) uX(t) uY(t)) +
     *                 2(uX(0) uZ(0) uX(t) uZ(t)) +
     *                 2(uY(0) uZ(0) uY(t) uZ(t))) - 1]/2
     *
     *               = [(3/2) * (<uX^2> + <uY^2> + <uZ^2> +
     *                         2<uXuY> + 2<uXuZ> + 2<uYuZ>) - 0.5]
     *
     */
    
    /* Because of normalization the number of -0.5 to subtract
     * depends on the number of data points!
     */
    for(j=0; (j<nf2); j++) 
      csum[j]  = -0.5*(nf2-j);
    
    /***** DIAGONAL ELEMENTS ************/
    for(m=0; (m<DIM); m++) {
      /* Copy the vector data in a linear array */
      for(j=0; (j<nf2); j++)
	ctmp[j]  = sqr(c1[DIM*j+m]);
      if (debug) {
	sprintf(buf,"c1diag%d.xvg",m);
	dump_tmp(buf,nf2,ctmp);
      }
      
      low_do_four_core(nfour,nf2,ctmp,cfour,enNorm,FALSE);
      
      if (debug) {
	sprintf(buf,"c1dfout%d.xvg",m);
	dump_tmp(buf,nf2,cfour);
      }
      fac = 1.5;
      for(j=0; (j<nf2); j++)
	csum[j] += fac*(cfour[j]);
    }
    /******* OFF-DIAGONAL ELEMENTS **********/
    for(m=0; (m<DIM); m++) {
      /* Copy the vector data in a linear array */
      m1=(m+1) % DIM;
      for(j=0; (j<nf2); j++)
	ctmp[j]=c1[DIM*j+m]*c1[DIM*j+m1];
      
      if (debug) {
	sprintf(buf,"c1off%d.xvg",m);
	dump_tmp(buf,nf2,ctmp);
      }
      low_do_four_core(nfour,nf2,ctmp,cfour,enNorm,FALSE);
      if (debug) { 
	sprintf(buf,"c1ofout%d.xvg",m);
	dump_tmp(buf,nf2,cfour);
      }
      fac = 3.0;
      for(j=0; (j<nf2); j++) {
	csum[j] += fac*cfour[j];
      }
    }
  }
  else if (MODE(eacP1) || MODE(eacVector)) {    
    /***************************************************
     * V E C T O R & P1
     ***************************************************/
    if (MODE(eacP1)) {
      /* First normalize the vectors */
      norm_and_scale_vectors(nframes,c1,1.0);
    }
    
    /* For vector thingies we have to do three FFT based correls 
     * First for XX, then for YY, then for ZZ
     * After that we sum them and normalise
     */
    for(j=0; (j<nf2); j++) {
      csum[j]=0.0;
    }
    for(m=0; (m<DIM); m++) {
      /* Copy the vector data in a linear array */
      for(j=0; (j<nf2); j++)
	ctmp[j]=c1[DIM*j+m];
      low_do_four_core(nfour,nf2,ctmp,cfour,enNorm,FALSE);
      for(j=0; (j<nf2); j++) 
	csum[j] += cfour[j];
    }
  }
  else
    fatal_error(0,"\nUnknown mode in do_autocorr (%d)",mode);
  
  sfree(cfour);
  for(j=0; (j<nf2); j++)
    c1[j] = csum[j]/(real)(nframes-j);
}

static real fit_acf(int ncorr,int fitfn,bool bVerbose,
		    real tbeginfit,real tendfit,real dt,real c1[],real *fit)
{
  real    fitparm[3];
  real    tStart,tail_corr,sum,sumtot=0,ct_estimate,*sig;
  int     i,j,jmax,nf_int;
  bool    bPrint;

  bPrint = bVerbose || bDebugMode();

  if (bPrint) printf("COR:\n");    
  
  if (tendfit <= 0)
    tendfit = ncorr*dt;
  nf_int = min(ncorr,(int)(tendfit/dt));
  sum    = print_and_integrate(debug,nf_int,dt,c1,NULL,1);

  /* Estimate the correlation time for better fitting */
  ct_estimate = 0.5*c1[0];
  for(i=1; (i<ncorr) && (c1[i]>0); i++)
      ct_estimate += c1[i];
  ct_estimate *= dt/c1[0];

  if (bPrint) printf("COR: Correlation time (plain integral from %6.3f to %6.3f ps) = %8.5f ps\n", 
		       0.0,dt*nf_int,sum);
  if (bPrint) printf("COR: Relaxation times are computed as fit to an exponential:\n");
  if (bPrint) printf("COR:   %s\n",longs_ffn[fitfn]);
  if (bPrint) printf("COR: Fit to correlation function from %6.3f ps to %6.3f ps, results in a\n",tbeginfit,min(ncorr*dt,tendfit));
  
  tStart = 0;
  if (bPrint) printf("COR:%11s%11s%11s%11s%11s%11s%11s\n",
		       "Fit from","Integral","Tail Value","Sum (ps)"," a1 (ps)",
		     (nfp_ffn[fitfn]>=2) ? " a2 ()" : "",
		     (nfp_ffn[fitfn]>=3) ? " a3 (ps)" : "");
  if (tbeginfit > 0)
    jmax = 3;
  else
    jmax = 1;
  if (fitfn == effnEXP3) {
    fitparm[0] = 0.002*ncorr*dt;
    fitparm[1] = 0.95;
    fitparm[2] = 0.2*ncorr*dt;
  } else {
    /* Good initial guess, this increases the probability of convergence */
    fitparm[0] = ct_estimate;
    fitparm[1] = 1.0;
    fitparm[2] = 1.0;
  }

  /* Generate more or less appropriate sigma's */
  snew(sig,ncorr);
  for(i=0; i<ncorr; i++)
    sig[i] = sqrt(ct_estimate+dt*i);

  for(j=0; ((j<jmax) && (tStart < tendfit)); j++) {
    /* Use the previous fitparm as starting values for the next fit */
    nf_int = min(ncorr,(int)((tStart+1e-4)/dt));
    sum    = print_and_integrate(debug,nf_int,dt,c1,NULL,1);
    tail_corr = do_lmfit(ncorr,c1,sig,dt,NULL,tStart,tendfit,
			 bDebugMode(),fitfn,fitparm,0);
    sumtot = sum+tail_corr;
    if (fit && ((jmax == 1) || (j == 1)))
      for(i=0; (i<ncorr); i++)
	fit[i] = fit_function(fitfn,fitparm,i*dt);
    if (bPrint) {
      printf("COR:%11.4e%11.4e%11.4e%11.4e",tStart,sum,tail_corr,sumtot);
      for(i=0; (i<nfp_ffn[fitfn]); i++)
	printf(" %11.4e",fitparm[i]);
      printf("\n");
    }
    tStart += tbeginfit;
  }
  sfree(sig);

  return sumtot;
}

void low_do_autocorr(char *fn,char *title,
		     int nframes,int nitem,int nout,real **c1,
		     real dt,unsigned long mode,int nrestart,
		     bool bAver,bool bFour,bool bNormalize,
		     bool bVerbose,real tbeginfit,real tendfit,
		     int eFitFn,int nskip)
{
  FILE    *fp,*gp=NULL;
  int     i,k,nfour;
  fftreal *csum;
  real    *ctmp,*fit;
  real    c0,sum,Ct2av,Ctav;
 
  /* Check flags and parameters */ 
  /*  nout = get_acfnout();*/
  if (nout == -1)
    nout = acf.nout = (nframes+1)/2;
  else if (nout > nframes)
    nout=nframes;
  
  if (MODE(eacCos) && MODE(eacVector))
    fatal_error(0,"Incompatible options bCos && bVector (%s, %d)",
		__FILE__,__LINE__);
  if ((MODE(eacP3) || MODE(eacRcross)) && bFour) {
    fprintf(stderr,"Can't combine mode %lu with FFT, turning off FFT\n",mode);
    bFour = FALSE;
  }
  if (MODE(eacNormal) && MODE(eacVector)) 
    fatal_error(0,"Incompatible mode bits: normal and vector (or Legendre)");
    
  /* Print flags and parameters */
  if (bVerbose) {
    printf("Will calculate %s of %d thingies for %d frames\n",
	   title ? title : "autocorrelation",nitem,nframes);
    printf("bAver = %s, bFour = %s bNormalize= %s\n",
	   bool_names[bAver],bool_names[bFour],bool_names[bNormalize]);
    printf("mode = %lu, dt = %g, nrestart = %d\n",mode,dt,nrestart);
  }
  if (bFour) {  
    c0 = log((double)nframes)/log(2.0);
    k  = c0;
    if (k < c0)
      k++;
    k++;
    nfour = 1<<k;
    if (debug)
      fprintf(debug,"Using FFT to calculate %s, #points for FFT = %d\n",
	      title,nfour);
	
    /* Allocate temp arrays */
    snew(csum,nfour);
    snew(ctmp,nfour);
  } else {
    nfour = 0; /* To keep the compiler happy */
    snew(csum,nframes);
    snew(ctmp,nframes);
  }
  
  /* Loop over items (e.g. molecules or dihedrals) 
   * In this loop the actual correlation functions are computed, but without
   * normalizing them.
   */
  k = max(1,pow(10,(int)(log(nitem)/log(100))));
  for(i=0; i<nitem; i++) {
    if (bVerbose && ((i%k==0 || i==nitem-1)))
      fprintf(stderr,"\rThingie %d",i+1);
    
    if (bFour)
      do_four_core(mode,nfour,nframes,nframes,c1[i],csum,ctmp);
    else 
      do_ac_core(nframes,nout,ctmp,c1[i],nrestart,mode);
  }
  if (bVerbose)
    fprintf(stderr,"\n");
  sfree(ctmp);
  sfree(csum);
  
  if (fn) {
    snew(fit,nout);
    fp=xvgropen(fn,title,"Time (ps)","C(t)");
  } else {
    fit = NULL;
    fp  = NULL;
  }
  if (bAver) {
    if (nitem > 1)
      average_acf(bVerbose,nframes,nitem,c1);
    
    if (bNormalize)
      normalize_acf(nout,c1[0]);
    
    if (eFitFn != effnNONE) {
      fit_acf(nout,eFitFn,fn!=NULL,tbeginfit,tendfit,dt,c1[0],fit);
      sum = print_and_integrate(fp,nout,dt,c1[0],fit,1);
    } else {
      sum = print_and_integrate(fp,nout,dt,c1[0],NULL,1);
      if (bVerbose)
	printf("Correlation time (integral over corrfn): %g (ps)\n",sum);
    }
  } else {
    /* Not averaging. Normalize individual ACFs */
    Ctav = Ct2av = 0;
    if (debug)
      gp = xvgropen("ct-distr.xvg","Correlation times","item","time (ps)");
    for(i=0; i<nitem; i++) {
      if (bNormalize)
	normalize_acf(nout,c1[i]);
      if (eFitFn != effnNONE) {
	fit_acf(nout,eFitFn,fn!=NULL,tbeginfit,tendfit,dt,c1[i],fit);
	sum = print_and_integrate(fp,nout,dt,c1[i],fit,1);
      } else {
	sum = print_and_integrate(fp,nout,dt,c1[i],NULL,1);
	if (debug)
	  fprintf(debug,
		  "CORRelation time (integral over corrfn %d): %g (ps)\n",
		  i,sum);
      }
      Ctav += sum;
      Ct2av += sum*sum;
      if (debug)
	fprintf(gp,"%5d  %.3f\n",i,sum);
    }
    if (debug)
      ffclose(gp);
    if (nitem > 1) {
      Ctav  /= nitem;
      Ct2av /= nitem;
      printf("Average correlation time %.3f Std. Dev. %.3f Error %.3f (ps)\n",
	     Ctav,sqrt((Ct2av - sqr(Ctav))),
	     sqrt((Ct2av - sqr(Ctav))/(nitem-1)));
    }
  }
  if (fp)
    ffclose(fp);
  sfree(fit);
}

static char *Leg[]   = { NULL, "0", "1", "2", "3", NULL };
static char *Nparm[] = { NULL, "1", "2", NULL };

t_pargs *add_acf_pargs(int *npargs,t_pargs *pa)
{
  t_pargs acfpa[] = {
    { "-acflen",     FALSE, etINT,  {&acf.nout},
      "Length of the ACF, default is half the number of frames" },
    { "-normalize",FALSE, etBOOL, {&acf.bNormalize},
      "Normalize ACF" },
    { "-fft",      FALSE, etBOOL, {&acf.bFour},
      "HIDDENUse fast fourier transform for correlation function" },
    { "-nrestart", FALSE, etINT,  {&acf.nrestart},
      "HIDDENNumber of frames between time origins for ACF when no FFT is used" },
    { "-P",        FALSE, etENUM, {Leg},
      "Order of Legendre polynomial for ACF (0 indicates none)" },
    { "-fitfn",    FALSE, etENUM, {s_ffn},
      "Fit function" },
    { "-ncskip",   FALSE, etINT,  {&acf.nskip},
      "Skip N points in the output file of correlation functions" },
    { "-beginfit", FALSE, etREAL, {&acf.tbeginfit},
      "Time where to begin the exponential fit of the correlation function" },
    { "-endfit",   FALSE, etREAL, {&acf.tendfit},
      "Time where to end the exponential fit of the correlation function, -1 is till the end" },
   };
#define NPA asize(acfpa)
  t_pargs *ppa;
  int i;
  
  snew(ppa,*npargs+NPA);
  for(i=0; (i<*npargs); i++)
    ppa[i] = pa[i];
  for(i=0; (i<NPA); i++)
    ppa[*npargs+i] = acfpa[i];
  (*npargs) += NPA;

  acf.mode       = 0;
  acf.nrestart   = 1;
  acf.nout       = -1;
  acf.P          = 0;
  acf.fitfn      = effnEXP1;
  acf.bFour      = TRUE;
  acf.bNormalize = TRUE;
  acf.tbeginfit  = 0.0;
  acf.tendfit    = -1;
  
  bACFinit = TRUE;
    
  return ppa;
}

void do_autocorr(char *fn,char *title,int nframes,int nitem,real **c1,
		 real dt,unsigned long mode,bool bAver)
{
  int i;

  if (!bACFinit) {
    printf("ACF data structures have not been initialised. Call add_acf_pargs\n");
  }

  /* Handle enumerated types */
  sscanf(Leg[0],"%d",&acf.P);
  acf.fitfn = sffn2effn(s_ffn);

  switch (acf.P) {
  case 1:
    mode = mode | eacP1;
    break;
  case 2:
    mode = mode | eacP2;
    break;
  case 3:
    mode = mode | eacP3;
    break;
  default:
    break;
  }
  
  low_do_autocorr(fn,title,nframes,nitem,acf.nout,c1,dt,mode,
		  acf.nrestart,bAver,acf.bFour,acf.bNormalize,
		  bDebugMode(),acf.tbeginfit,acf.tendfit,
		  acf.fitfn,acf.nskip);
}

int get_acfnout(void)
{
  if (!bACFinit)
    fatal_error(0,"ACF data not initialized yet");

  return acf.nout;
}

int get_acffitfn(void)
{
  if (!bACFinit)
    fatal_error(0,"ACF data not initialized yet");

  return sffn2effn(s_ffn);
}
