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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_autocorr_c = "$Id$";

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

typedef struct {
  unsigned long mode;
  int  nrestart,nlag,P,nfitparm;
  bool bFull,bFour,bNormalize;
  real tbeginfit,tendfit;
} t_acf;

static bool  bACFinit = FALSE;
static t_acf acf;

typedef real fftreal;

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
    ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
    ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  realft(ans,no2,-1);
  sfree(fft);
}

static void low_do_four_core(int nfour,int nframes,real c1[],fftreal cfour[],
			     int nCos,bool bPadding)
{
  int  i,no2;
  fftreal aver,dum,*ans;

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
    /* fprintf(stderr,"aver = %g\n",aver); */
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

static void do_ac_core(int nframes,int nlag,
		       real corr[],real c1[],int nrestart,
		       unsigned long mode,bool bFull)
{
  int     j,k,j3,jk3,m,n;
  fftreal ccc,c0,cth;
  rvec    xj,xk,rr;

  if (bFull) {
    if (nrestart != 1) 
      fprintf(stderr,"WARNING: setting number of restarts to 1 for Full ACF\n");
    nrestart = 1;
  }
  else {
    if (nrestart < 1) {
      fprintf(stderr,"WARNING: setting number of restarts to 1\n");
      nrestart = 1;
    }
  }
  if (debug)
    fprintf(debug,"Starting do_ac_core: nframes=%d, nlag=%d, nrestart=%d,mode=%d\n",
	    nframes,nlag,nrestart,mode);
  
  for(j=0; (j<nlag); j++)
    corr[j]=0;
  
  /* Loop over starting points. */
  for(j=0; (j<nframes); j+=nrestart) {
    j3  = DIM*j;
    
    /* Loop over the correlation length for this starting point */
    for(k=0; (k<nlag) && (j+k < nframes); k++) {
      jk3 = DIM*(j+k);
      
      /* Switch over possible ACF types. 
       * It might be more efficient to put the loops inside the switch,
       * but this is more clear, and save development time!
       */      
      switch (mode) {
      case eacNormal:
	corr[k] += c1[j]*c1[j+k];
	break;
    
      case eacCos:
	/* Compute the cos (phi(t)-phi(t+dt)) */
	corr[k] += cos(c1[j]-c1[j+k]);
	break;
    
      case eacP1:
      case eacP2:
      case eacP3:
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	cth=cos_angle(xj,xk);
	
	if (cth-1.0 > 1.0e-15) {
	  fprintf(stderr,"j: %d, k: %d, xj:(%g,%g,%g), xk:(%g,%g,%g)\n",
		  j,k,xj[XX],xj[YY],xj[ZZ],xk[XX],xk[YY],xk[ZZ]);
	}
	
	corr[k] += LegendreP(cth,mode);  /* 1.5*cth*cth-0.5; */
	break;
    
      case eacRcross:
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	oprod(xj,xk,rr);
	
	corr[k] += iprod(rr,rr);
	
	break;
    
      case eacVector:
	for(m=0; (m<DIM); m++) {
	  xj[m] = c1[j3+m];
	  xk[m] = c1[jk3+m];
	}
	ccc = iprod(xj,xk);
	
	corr[k] += ccc;
      
	break;
      default:
	fatal_error(0,"\nInvalid mode (%d) in do_ac_core",mode);
      }
    }
  }
  /* Correct for the number of points and copy results to the data array */
  for(j=0; (j<nlag); j++) {
    n = (nframes-j+(nrestart-1))/nrestart;
    c1[j] = corr[j]/n;
  }
}

void normalize_acf(int nframes,int nf2,int nlag,
		   real corr[],bool bFull,bool bFour,bool bNormalize)
{
  int  j;
  real c0;

  if (debug) {
    fprintf(debug,"Before normalization\n");
    for(j=0; (j<nf2); j++) 
      fprintf(debug,"%5d  %10f\n",j,corr[j]);
  }
  
  /* Normalisation makes that c[0] = 1.0 and that other points are scaled
   * accordingly.
   */
  if (bNormalize) {
    if (corr[0] == 0.0)
      c0 = 1.0;
    else
      c0 = 1.0/corr[0];
    for(j=0; (j<nlag); j++)
      corr[j] *= c0;
  }
  if (debug) {
    fprintf(debug,"After normalization\n");
    for(j=0; (j<nf2); j++) 
      fprintf(debug,"%5d  %10f\n",j,corr[j]);
  }
}

void average_acf(int ncorr,int nitem,real **c1)
{
  real c0;
  int  i,j;
  
  fprintf(stderr,"Averaging correlation functions\n");
  
  for(j=0; (j<ncorr); j++) {
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

real integrate(FILE *fp,int n,real dt,real c[])
{
  real c0,sum;
  int  j;
  
  /* Use trapezoidal rule for calculating integral */
  sum = 0.0;
  for(j=0; (j<n); j++) {
    c0 = c[j];
    if (fp)
      fprintf(fp,"%10.3f  %10.5f\n",j*dt,c0);
    if (j > 0)
      sum+=dt*(c0+c[j-1]);
  }
  return sum*0.5;
}

void do_four_core(unsigned long mode,int nfour,int nf2,int nframes,
		  real c1[],real csum[],real ctmp[])
{
  fftreal *cfour;
  char    buf[32];
  real    fac;
  int     i,j,m,m1;
  
  snew(cfour,nfour);
  
  switch (mode) {
    /********************************************
     *  N O R M A L
     ********************************************/
  case eacNormal:
    /********** F F T ********/
    low_do_four_core(nfour,nf2,c1,cfour,enNorm,FALSE);
    break;
    
    /***************************************************
     * C O S I N E
     ***************************************************/
  case eacCos:
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
    for(j=0; (j<nf2); j++)
      c1[j] += cfour[j];
    
    break;
    
    /***************************************************
     * Legendre polynomials
     ***************************************************/
  case eacP2: 
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
    break;
    
    /***************************************************
     * V E C T O R & P1
     ***************************************************/
  case eacP1:
    /* First normalize the vectors */
    norm_and_scale_vectors(nframes,c1,1.0);
    /* Fall thru, don't break */
    
  case eacVector:
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
    break;
  default:
    fatal_error(0,"\nUnknown mode in do_autocorr (%d)",mode);
  }
  sfree(cfour);
}

void fit_acf(int ncorr,int nfitparm,
	     char *fitfn,char *fittitle,bool bVerbose,
	     real tbeginfit,real tendfit,real dt,real c1[])
{
  real    fitparm[3],fit[3];
  real    tStart,tail_corr,sum,sumtot,*sig;
  int     j,nf_int;
  
  fprintf(stderr,"CORR:\n");    
  
  nf_int = min(ncorr,(int)(tendfit/dt));
  sum    = integrate(debug,ncorr,dt,c1);
  
  fprintf(stderr,"CORR: Correlation time (plain integral from %6.3f to %6.3f ps) = %8.5f ps\n", 
	  0.0,dt*nf_int,sum);
  fprintf(stderr,"CORR: Relaxation times are computed as fit to an exponential:\n");
  if (nfitparm == 1)
    fprintf(stderr,"CORR:    Exp[-t/tau_slope]\n");
  else if (nfitparm == 2)
    fprintf(stderr,"CORR:    A Exp[-t/tau_slope]\n");
  else 
    fatal_error(0,"nparm not set to 1 or 2, %s %d",__FILE__,__LINE__);
  fprintf(stderr,"CORR: Fit to correlation function from %6.3f ps to %6.3f ps, results in a\n",tbeginfit,min(ncorr*dt,tendfit));
    
  tStart = 0;
  fprintf(stderr,"CORR:%12s%12s%12s%12s%12s\n",
	  "Integral to","Value","Tail Value","Sum (ps)","Tau (ps)");
  for(j=0; ((j<5) && (tStart < tendfit)); j++) {
    snew(sig,ncorr);
    fitparm[0]=fitparm[1]=fitparm[2] = 1.0;
    nf_int = min(ncorr,(int)((tStart+1e-4)/dt));
    sum    = integrate(debug,nf_int,dt,c1);
    tail_corr = do_lmfit(ncorr,c1,sig,dt,tStart,tendfit,
			 fitfn,fittitle,bVerbose,nfitparm,
			 NULL,fitparm,NULL);
    sumtot = sum+tail_corr;
    fprintf(stderr,"CORR:%12.5e%12.5e%12.5e%12.5e%12.5e\n",
	    tStart,sum,tail_corr,sumtot,fitparm[0]);
      
    tStart += tbeginfit;
    sfree(sig);
  }
}

void low_do_autocorr(char *fn,char *title,
		     int nframes,int nitem,int nlag,real **c1,
		     real dt,unsigned long mode,int nrestart,
		     bool bFull,bool bAver,bool bFour,bool bNormalize,
		     char *fitfn,char *fittitle,bool bVerbose,
		     real tbeginfit,real tendfit,
		     int nfitparm)
{
  FILE    *fp;
  const   real sqrtsqrt15=sqrt(sqrt(1.5));
  int     i,j,j3,m,m1,k,ncorr,nfour;
  fftreal *csum;
  char    buf[256];
  real    *ctmp,*rij,*sig;
  real    dc,c0,sum,rnorm,fac;
 
  /* Check flags and parameters */ 
  nlag = get_acflag();
  if (bFull || bFour) {
    nlag  = acf.nlag = nframes;
    ncorr = nframes;
  }
  else {
    if (nlag == -1)
      nlag = acf.nlag = (nframes+1)/2;
    else if (nlag > nframes)
      nlag=nframes;
    ncorr = nlag;
  }
  
  if ((mode & eacCos) && (mode & eacVector))
    fatal_error(0,"Incompatible options bCos && bVector (%s, %d)",
		__FILE__,__LINE__);
  if (bFull && bFour) {
    fprintf(stderr,"Turning off FFT! (Can't be done with Full)\n");
    bFour=FALSE;
  }
  if (((mode == eacP3) || (mode == eacRcross)) && bFour) {
    fprintf(stderr,"Cant combine mode %d with FFT, turning off FFT\n",mode);
    bFour = FALSE;
  }
    
  /* Print flags and parameters */
  fprintf(stderr,"Will calculate %s of %d thingies for %d frames\n",
	  title,nitem,ncorr);
  fprintf(stderr,"bFull = %s, bAver = %s, bFour = %s bNormalize= %s\n",
	  bool_names[bFull],bool_names[bAver],
	  bool_names[bFour],bool_names[bNormalize]);
  fprintf(stderr,"mode = %d, dt = %g, nrestart = %d\n",mode,dt,nrestart);
  
  if (bFour) {  
    c0 = log((double)ncorr)/log(2.0);
    k  = c0;
    if (k < c0)
      k++;
    k++;
    nfour = pow(2,k);
    fprintf(stderr,"Using FFT to calculate %s, #points for FFT = %d\n",
	    title,nfour);
	
    /* Allocate temp arrays */
    snew(csum,nfour);
    snew(ctmp,nfour);
  }
  else {
    nfour = 0; /* To keep the compiler happy */
    snew(csum,nframes);
    snew(ctmp,nframes);
  }
  
  /* Loop over items (e.g. molecules or dihedrals) 
   * In this loop the actual correlation functions are computed, but without
   * normalizing them.
   */
  for(i=0; (i<nitem); i++) {
    fprintf(stderr,"\rThingie %d",i);
    
    if (bFour)
      do_four_core(mode,nfour,ncorr,nframes,c1[i],csum,ctmp);
    else 
      do_ac_core(nframes,nlag,ctmp,c1[i],nrestart,mode,bFull);
  }
  fprintf(stderr,"\n");
  sfree(ctmp);
  sfree(csum);
  
  fp=xvgropen(fn,title,"Time (ps)","C(t)");
  if (bAver) {
    average_acf(ncorr,nitem,c1);
    
    normalize_acf(nframes,ncorr,nlag,c1[0],bFull,bFour,bNormalize);

    if (tbeginfit < tendfit)
      fit_acf(ncorr,nfitparm,fitfn,fittitle,bVerbose,
	      tbeginfit,tendfit,dt,c1[0]);
    else {
      sum = integrate(fp,ncorr,dt,c1[0]);
      fprintf(stderr,"Correlation time (integral over corrfn): %g (ps)\n",sum);
    }
  }
  else {
    /* Not averaging. Normalize individual ACFs */
    for(i=0; (i<nitem); i++) 
      normalize_acf(nframes,ncorr,nlag,c1[i],bFull,bFour,bNormalize);
      
    /* Now dump them all */
    for(j=0; (j<ncorr); j++) {
      fprintf(fp,"%10f",j*dt);
      for(i=0; (i<nitem); i++) {
	if (((mode == eacP1) || (mode == eacP2) || (mode == eacP3)) && bFour)
	  fprintf(fp,"  %10.5f",c1[i][j]);
	else
	  fprintf(fp,"  %10.5f",c1[i][j]);
      }
      fprintf(fp,"\n");
    }
  }
  ffclose(fp);
  
}

t_pargs *add_acf_pargs(int *npargs,t_pargs *pa)
{
  t_pargs acfpa[] = {
    { "-fft",      FALSE, etBOOL, &acf.bFour,
      "Use fast fourier transform for correlation function" },
    { "-full",     FALSE, etBOOL,  &acf.bFull,
      "HIDDENCompute full ACF leading to inaccurate tail" },
    { "-normalize",FALSE, etBOOL, &acf.bNormalize,
      "Normalize ACF" },
    { "-nrestart", FALSE, etINT,  &acf.nrestart,
      "Number of frames between time origins for ACF when no FFT is used" },
    { "-acflen",     FALSE, etINT,  &acf.nlag,
      "Length of the ACF when no FFT is used, default is half the number of frames" },
    { "-P",        FALSE, etINT,  &acf.P,
      "Order of Legendre polynomial for ACF (1,2 or 3)" },
    { "-nparm",    FALSE, etINT,  &acf.nfitparm,
      "Number of parameters in exponential fit (1 or 2)" },
    { "-beginfit", FALSE, etREAL, &acf.tbeginfit,
      "Time where to begin the exponential fit of the correlation function" },
    { "-endfit",   FALSE, etREAL, &acf.tendfit,
      "Time where to end the exponential fit of the correlation function" },
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
  acf.nlag       = -1;
  acf.P          = 0;
  acf.nfitparm   = 1;
  acf.bFull      = FALSE;
  acf.bFour      = TRUE;
  acf.bNormalize = TRUE;
  acf.tbeginfit  = 0.0;
  acf.tendfit    = 0.0;
  
  bACFinit = TRUE;
    
  return ppa;
}

void do_autocorr(char *fn,char *title,int nframes,int nitem,real **c1,
		 real dt,unsigned long mode,bool bAver,
		 char *fitfn,char *fittitle)
{
  if (!bACFinit) {
    fprintf(stderr,"ACF data structures have not been initialised. Call add_acf_pargs\n");
  }

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
  
  low_do_autocorr(fn,title,nframes,nitem,acf.nlag,c1,dt,mode,
		  acf.nrestart,acf.bFull,bAver,acf.bFour,acf.bNormalize,
		  fitfn,fittitle,bDebugMode(),acf.tbeginfit,acf.tendfit,
		  acf.nfitparm);
}

int get_acflag(void)
{
  if (!bACFinit)
    fatal_error(0,"ACF data not initialized yet");

  return acf.nlag;
}
