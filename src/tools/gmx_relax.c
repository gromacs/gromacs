/*
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

#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "gstat.h"
#include "tpxio.h"
#include "gmx_ana.h"


typedef struct {
  real re,im;
} complex;

typedef struct {
  int ai,aj;
} t_pair;

typedef struct {
  real    rij_3;
  real    rij_6;
  gmx_bool    bNOE;
  real    tauc,dtauc,S2,dS2;
  complex y2;
  complex Ylm[5];
} t_sij;

void read_shifts(FILE *fp,int nx,real shiftx[],int ny,real shifty[])
{
  int    i;
  double d;
  
  for(i=0; (i<nx); i++) {
    fscanf(fp,"%lf",&d);
    shiftx[i]=d;
    shifty[i]=d;
  }
  /*for(i=0; (i<ny); i++) {
    fscanf(fp,"%lf",&d);
    shifty[i]=d;
    }
    */
}

complex c_sqr(complex c)
{
  complex cs;
  
  cs.re = c.re*c.re-c.im*c.im;
  cs.im = 2.0*c.re*c.im;
  
  return cs;
}

complex c_add(complex c,complex d)
{
  complex cs;
  
  cs.re = c.re+d.re;
  cs.im = c.im+d.im;
  
  return cs;
}

complex calc_ylm(int m,rvec rij,real r2,real r_3,real r_6)
{
  static  gmx_bool bFirst=TRUE;
  static  real y0,y1,y2;
  real    x,y,z,xsq,ysq,rxy,r1,cphi,sphi,cth,sth,fac;
  complex cs;
  
  if (bFirst) {
    y0 = sqrt(5/(4*M_PI));
    y1 = sqrt(45/(24*M_PI));
    y2 = sqrt(45/(96*M_PI));
    bFirst = FALSE;
  }
  r1  = sqrt(r2);
  x   = rij[XX];
  y   = rij[YY];
  z   = rij[ZZ];
  xsq = x*x;
  ysq = y*y;
  rxy = sqrt(xsq+ysq);
  cphi= x/rxy;
  sphi= y/rxy;
  cth = z/r1;
  if (cphi != 0.0)
    sth = x/(r1*cphi);
  else
    sth = y/(r1*sphi);
  
  /* Now calculate the spherical harmonics */
  switch(m) {
  case -2:
  case  2:
    fac    = y2*sth*sth;
    cs.re  = fac*(cphi*cphi-sphi*sphi);
    /* Use the index as a prefactor... check your formulas */
    cs.im  = m*fac*cphi*sphi;
    break;
  case -1:
  case  1:
    fac    = y1*cth*sth;
    cs.re  = -m*fac*cphi;
    cs.im  = -fac*sphi;
    break;
  case 0:
    cs.re = y0*(3*cth*cth-1);
    cs.im = 0;
    break;
  }
  cs.re *= r_3;
  cs.im *= r_3;
  
  return cs;
}

void myfunc(real x,real a[],real *y,real dyda[],int na)
{
  /* Fit to function 
   *
   * y = a1 + (1-a1) exp(-a2 x)
   *
   * where in real life a1 is S^2 and a2 = 1/tauc
   *
   */
   
  real eee,S2,tau1;
  
  S2      = a[1];
  tau1    = 1.0/a[2];
  eee     = exp(-x*tau1);
  *y      = S2 + (1-S2)*eee;
  dyda[1] = 1 - eee;
  dyda[2] = x*tau1*tau1*(1-a[1])*eee;
}

void fit_one(gmx_bool bVerbose,
	     int nframes,real x[],real y[],real dy[],real ftol,
	     real *S2,real *dS2,real *tauc,real *dtauc)
{
  void mrqmin(real x[],real y[],real sig[],int ndata,real a[],
	      int ma,int lista[],int mfit,real **covar,real **alpha,
	      real *chisq,
	      void (*funcs)(real x,real a[],real *y,real dyda[],int na),
	      real *alamda);
	      
  real *a,**covar,**alpha;
  real chisq,ochisq,alamda;
  gmx_bool bCont;
  int  i,j,ma,mfit,*lista;
  
  ma=mfit=2;
  snew(a,ma+1);
  snew(covar,ma+1);
  snew(alpha,ma+1);
  snew(lista,ma+1);
  for(i=0; (i<ma+1); i++) {
    lista[i] = i;
    snew(covar[i],ma+1);
    snew(alpha[i],ma+1);
  }

  a[1]   = 0.99;  /* S^2              */
  a[2]   = 0.1;   /* tauc             */
  alamda = -1;    /* Starting value   */
  chisq  = 1e12;
  j      = 0;      
  do {
    ochisq = chisq;
    mrqmin(x-1,y-1,dy-1,nframes,a,ma,lista,mfit,covar,alpha,
	   &chisq,myfunc,&alamda);
    if (bVerbose)
      fprintf(stderr,"\rFitting %d chisq=%g, alamda=%g, tau=%g, S^2=%g\t\t\n",
	      j,chisq,alamda,1.0/a[2],a[1]);
    j++;
    bCont = (((ochisq - chisq) > ftol*chisq) ||
	     ((ochisq == chisq)));
  } while (bCont && (alamda != 0.0) && (j < 50));
  if (bVerbose)
    fprintf(stderr,"\n");

  /* Now get the covariance matrix out */
  alamda = 0;
  mrqmin(x-1,y-1,dy-1,nframes,a,ma,lista,mfit,covar,alpha,
	 &chisq,myfunc,&alamda);

  *S2    = a[1];
  *dS2   = sqrt(covar[1][1]);
  *tauc  = a[2];
  *dtauc = sqrt(covar[2][2]);
  
  for(i=0; (i<ma+1); i++) {
    sfree(covar[i]);
    sfree(alpha[i]);
  }
  sfree(a);
  sfree(covar);
  sfree(alpha);
  sfree(lista);
}

void calc_tauc(gmx_bool bVerbose,int npair,t_pair pair[],real dt,
	       int nframes,t_sij spec[],real **corr)
{
  FILE *fp;
  char buf[32];
  int  i,j,k,n;
  real S2,S22,tauc,fac;
  real *x,*dy;
  real ftol = 1e-3;
  
  snew(x,nframes);
  snew(dy,nframes);
  for(i=0; (i<nframes); i++) 
    x[i]  = i*dt;
  
  fprintf(stderr,"Fitting correlation function to Lipari&Szabo function\n");
  fac=1.0/((real)nframes);
  for(i=0; (i<npair); i++) {
    if (spec[i].bNOE) {
      /* Use Levenberg-Marquardt method to fit */
      for(j=0; (j<nframes); j++) 
	dy[j] = fac;
      fit_one(bVerbose,
	      nframes,x,corr[i],dy,ftol,
	      &(spec[i].S2),&(spec[i].dS2),
	      &(spec[i].tauc),&(spec[i].dtauc));
      if (bVerbose) {
	sprintf(buf,"test%d.xvg",i);
	fp = ffopen(buf,"w");
	for(j=0; (j<nframes); j++) {
	  fprintf(fp,"%10g  %10g  %10g\n",j*dt,corr[i][j],
		  spec[i].S2 + (1-spec[i].S2)*exp(-j*dt/spec[i].tauc));
	}
	ffclose(fp);
      }
    }
  }
}

void calc_aver(FILE *fp,int nframes,int npair,t_pair pair[],t_sij *spec,
	       real maxdist)
{
  int     i,j,m;
  real    nf_1,fac,md_6;
  complex c1,c2,dc2;
  
  md_6 = pow(maxdist,-6.0);
  fac  = 4*M_PI/5;
  nf_1 = 1.0/nframes;
  for(i=0; (i<npair); i++) {
    c2.re = 0;
    c2.im = 0;
    fprintf(fp,"%5d  %5d",pair[i].ai,pair[i].aj);
    for(m=0; (m<5); m++) {
      c1.re  = spec[i].Ylm[m].re*nf_1;
      c1.im  = spec[i].Ylm[m].im*nf_1;
      dc2    = c_sqr(c1);
      c2     = c_add(dc2,c2);
      
      if (c1.im > 0)
	fprintf(fp,"  %8.3f+i%8.3f",c1.re,c1.im);
      else
	fprintf(fp,"  %8.3f-i%8.3f",c1.re,-c1.im);
    }
    fprintf(fp,"\n");
    spec[i].rij_3 *= nf_1;
    spec[i].rij_6 *= nf_1;
    spec[i].y2.re  = fac*c2.re;
    spec[i].y2.im  = fac*c2.im;
    spec[i].bNOE   = (spec[i].rij_6 > md_6);
  }
}

void plot_spectrum(char *noefn,int npair,t_pair pair[],t_sij *spec,real taum)
{
  FILE    *fp,*out;
  int     i,j,m;
  t_rgb   rlo = { 1,0,0 },rhi = {1,1,1};
  real    Sijmax,Sijmin,pow6,pow3,pp3,pp6,ppy,tauc;
  real    *Sij;
  complex sij;
  
  snew(Sij,npair);
  Sijmax = -1000.0;
  Sijmin =  1000.0;
  fp=xvgropen(noefn,"Cross Relaxation","Pair Index","\\8s\\4\\sij\\N");
  for(i=0; (i<npair); i++) {
    tauc   = spec[i].tauc;
    sij.re = -0.4*((taum-tauc)*spec[i].y2.re + tauc*spec[i].rij_6);
    sij.im = -0.4* (taum-tauc)*spec[i].y2.im;
    Sij[i]=sij.re;
    Sijmax=max(Sijmax,sij.re);
    Sijmin=min(Sijmin,sij.re);
    fprintf(fp,"%5d  %10g\n",i,sij.re);
  }
  ffclose(fp);
  fprintf(stderr,"Sijmin: %g, Sijmax: %g\n",Sijmin,Sijmax);
  out=ffopen("spec.out","w");
  pow6 = -1.0/6.0;
  pow3 = -1.0/3.0;
  fprintf(out,"%5s  %5s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
	  "at i","at j","S^2","Sig S^2","tauc","Sig tauc",
	  "<rij6>","<rij3>","<ylm>","rij3-6","ylm-rij6");
  for(i=0; (i<npair); i++) {
    if (spec[i].bNOE) {
      pp6 = pow(spec[i].rij_6,pow6);
      pp3 = pow(spec[i].rij_3,pow3);
      if (spec[i].y2.re < 0)
	ppy = -pow(-spec[i].y2.re,pow6);
      else
	ppy = pow(spec[i].y2.re,pow6);
      fprintf(out,"%5d  %5d  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n",
	      pair[i].ai,pair[i].aj,
	      spec[i].S2,spec[i].dS2,spec[i].tauc,spec[i].dtauc,
	      pp6,pp3,ppy,pp3-pp6,ppy-pp6);
    }
  }
  ffclose(out);
  
  sfree(Sij);
  
  do_view(noefn,NULL);
}

void spectrum(gmx_bool bVerbose,
	      char *trj,char *shifts,gmx_bool bAbInitio,
	      char *corrfn,char *noefn,
	      int maxframes,gmx_bool bFour,gmx_bool bFit,int nrestart,
	      int npair,t_pair pair[],int nat,real chem_shifts[],
	      real taum,real maxdist,
	      real w_rls[],rvec xp[],t_idef *idef)
{
  FILE   *fp;
  int    i,j,m,ii,jj,natoms,status,nframes;
  rvec   *x,dx;
  matrix box;
  real   t0,t1,t,dt;
  real   r2,r6,r_3,r_6,tauc;
  rvec   **corr;
  real   **Corr;
  t_sij  *spec;
  gmx_rmpbc_t  gpbc=NULL;

  snew(spec,npair);
  
  fprintf(stderr,"There is no kill like overkill! Going to malloc %d bytes\n",
	  npair*maxframes*sizeof(corr[0][0]));
  snew(corr,npair);
  for(i=0; (i<npair); i++)
    snew(corr[i],maxframes);
  nframes = 0;
  natoms  = read_first_x(&status,trj,&t0,&x,box);
  if (natoms > nat)
    gmx_fatal(FARGS,"Not enough atoms in trajectory");
  gpbc = gmx_rmpbc_init(idef,ePBC,natoms,box);

  do {
    if (nframes >= maxframes) {
      fprintf(stderr,"\nThere are more than the %d frames you told me!",
	      maxframes);
      break;
    }
    t1 = t;
    if (bVerbose)
      fprintf(stderr,"\rframe: %d",nframes);
    gmx_rmpbc(gpbc,box,x,x);
    if (bFit)
      do_fit(natoms,w_rls,xp,x);  
    
    for(i=0; (i<npair); i++) {
      ii = pair[i].ai;
      jj = pair[i].aj;
      rvec_sub(x[ii],x[jj],dx);
      copy_rvec(dx,corr[i][nframes]);
      
      r2  = iprod(dx,dx);
      r6  = r2*r2*r2;
      r_3 = gmx_invsqrt(r6);
      r_6 = r_3*r_3;
      spec[i].rij_3 += r_3;
      spec[i].rij_6 += r_6;
      for(m=0; (m<5); m++) {
	spec[i].Ylm[m] = c_add(spec[i].Ylm[m],
			       calc_ylm(m-2,dx,r2,r_3,r_6));
      }
    }
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  close_trj(status);
  if (bVerbose)
    fprintf(stderr,"\n");
 
  gmx_rmpbc_done(gpbc);

  fp=ffopen("ylm.out","w");
  calc_aver(fp,nframes,npair,pair,spec,maxdist);
  ffclose(fp);
 
  /* Select out the pairs that have to be correlated */
  snew(Corr,npair);
  for(i=j=0; (i<npair); i++) {
    if (spec[i].bNOE) {
      Corr[j] = &(corr[i][0][0]);
      j++;
    }
  }
  fprintf(stderr,"There are %d NOEs in your simulation\n",j);
  if (nframes > 1)
    dt = (t1-t0)/(nframes-1);
  else
    dt = 1;
  do_autocorr(corrfn,"Correlation Function for Interproton Vectors",
	      nframes,j,Corr,dt,eacP2,nrestart,FALSE,FALSE,bFour,TRUE);
  
  calc_tauc(bVerbose,npair,pair,dt,nframes/2,spec,(real **)corr);
  
  plot_spectrum(noefn,npair,pair,spec,taum);
}

int gmx_relax(int argc,char *argv[])
{
  const char *desc[] = {
    "g_noe calculates a NOE spectrum"
  };

  int        status;
  t_topology *top;
  int        i,j,k,natoms,nprot,*prot_ind;
  int        ifit;
  char       *gn_fit;
  atom_id    *ind_fit,*all_at;
  real       *w_rls;
  rvec       *xp;
  t_pair     *pair;
  matrix     box;
  int        step,nre;
  real       t,lambda;
  real       *shifts=NULL;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL,     ffREAD },
    { efTPX, "-s", NULL,     ffREAD },
    { efNDX, NULL, NULL,     ffREAD },
    { efDAT, "-d", "shifts", ffREAD },
    { efOUT, "-o","spec",    ffWRITE },
    { efXVG, "-corr", "rij-corr", ffWRITE },
    { efXVG, "-noe", "noesy", ffWRITE }
  };
#define NFILE asize(fnm)
  static real taum      = 0.0, maxdist = 0.6;
  static int  nlevels   = 15;
  static int  nrestart  = 1;
  static int  maxframes = 100;
  static gmx_bool bFFT      = TRUE,bFit = TRUE, bVerbose = TRUE;
  t_pargs pa[] = {
    { "-taum",     FALSE, etREAL, &taum, 
      "Rotational correlation time for your molecule. It is obligatory to pass this option" },
    { "-maxdist",  FALSE, etREAL, &maxdist,
      "Maximum distance to be plotted" },
    { "-nlevels",  FALSE, etINT,  &nlevels,
      "Number of levels for plotting" },
    { "-nframes", FALSE, etINT,  &maxframes,
      "Number of frames in your trajectory. Will stop analysis after this" },
    { "-fft",      FALSE, etBOOL, &bFFT,
      "Use FFT for correlation function" },
    { "-nrestart", FALSE, etINT,  &nrestart,
      "Number of frames between starting point for computation of ACF without FFT" },
    { "-fit",      FALSE, etBOOL, &bFit,
      "Do an optimal superposition on reference structure in [TT].tpx[tt] file" },
    { "-v",        FALSE, etBOOL, &bVerbose,
      "Tell you what I am about to do" }
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  if (taum <= 0)
    gmx_fatal(FARGS,"Please give me a sensible taum!\n");
  if (nlevels > 50) {
    nlevels = 50;
    fprintf(stderr,"Warning: too many levels, setting to %d\n",nlevels);
  }
  		    
  top    = read_top(ftp2fn(efTPX,NFILE,fnm));
  natoms = top->atoms.nr;
  snew(xp,natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,NULL,box,
	   &natoms,xp,NULL,NULL,NULL);

  /* Determine the number of protons, and their index numbers 
   * by checking the mass 
   */
  nprot  = 0;
  snew(prot_ind,natoms);
  for(i=0; (i<natoms); i++)
    if (top->atoms.atom[i].m  < 2) {
      prot_ind[nprot++] = i;
    }
  fprintf(stderr,"There %d protons in your topology\n",nprot);
  snew(pair,(nprot*(nprot-1)/2));
  for(i=k=0; (i<nprot); i++) {
    for(j=i+1; (j<nprot); j++,k++) {
      pair[k].ai = prot_ind[i];
      pair[k].aj = prot_ind[j];
    }
  }
  sfree(prot_ind);
  
  fprintf(stderr,"Select group for root least squares fit\n");
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&ifit,&ind_fit,&gn_fit);
  
  if (ifit < 3) 
    gmx_fatal(FARGS,"Need >= 3 points to fit!\n");

  /* Make an array with weights for fitting */
  snew(w_rls,natoms);
  for(i=0; (i<ifit); i++)
    w_rls[ind_fit[i]]=top->atoms.atom[ind_fit[i]].m;
    
  /* Prepare reference frame */
  snew(all_at,natoms);
  for(j=0; (j<natoms); j++)
    all_at[j]=j;
  rm_pbc(&(top->idef),natoms,box,xp,xp);
  reset_x(ifit,ind_fit,natoms,all_at,xp,w_rls);
  sfree(all_at);
  
  spectrum(bVerbose,
	   ftp2fn(efTRX,NFILE,fnm),ftp2fn(efDAT,NFILE,fnm),
	   ftp2bSet(efDAT,NFILE,fnm),opt2fn("-corr",NFILE,fnm),
	   opt2fn("-noe",NFILE,fnm),
	   maxframes,bFFT,bFit,nrestart,
	   k,pair,natoms,shifts,
	   taum,maxdist,w_rls,xp,&(top->idef));
  
  thanx(stderr);
  
  return 0;
}

