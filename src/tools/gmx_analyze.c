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
#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "readinp.h"
#include "statutil.h"
#include "txtdump.h"
#include "gstat.h"
#include "xvgr.h"

/* must correspond to char *avbar_opt[] declared in main() */
enum { avbarSEL, avbarNONE, avbarSTDDEV, avbarERROR, avbar90, avbarNR };

static real **read_val(char *fn,bool bHaveT,bool bTB,real tb,bool bTE,real te,
		       int nsets_in,int *nset,int *nval,real *dt,real **t,
		       int linelen)
{
  FILE   *fp;
  static int  llmax=0;
  static char *line0=NULL;
  char   *line;
  int    a,narg,n,sin,set,nchar;
  double dbl,tend=0;
  bool   bEndOfSet,bTimeInRange,bFirstLine=TRUE;
  real   **val;
  
  if (linelen > llmax) {
    llmax = linelen;
    srenew(line0,llmax);
  }
  
  val = NULL;
  *t  = NULL;
  *dt = 0;
  fp  = ffopen(fn,"r");
  for(sin=0; sin<nsets_in; sin++) {
    if (nsets_in == 1)
      narg = 0;
    else 
      narg = bHaveT ? 2 : 1;
    n = 0;
    bEndOfSet = FALSE;
    while (!bEndOfSet && fgets(line0,linelen,fp)) {
      line = line0;
      bEndOfSet = (line[0] == '&');
      if (line[0]!='#' && line[0]!='@' && line[0]!='\n' && !bEndOfSet) {
	if (bFirstLine && bHaveT) {
	  /* Check the first line that should contain data */
	  a = sscanf(line,"%lf%lf",&dbl,&dbl);
	  if (a == 0) 
	    gmx_fatal(FARGS,"Expected a number in %s on line:\n%s",fn,line0);
	  else if (a == 1) {
	    fprintf(stderr,"Found only 1 number on line, "
		    "assuming no time is present.\n");
	    bHaveT = FALSE;
	    if (nsets_in > 1)
	      narg = 1;
	  }
	}

	a = 0;
	bTimeInRange = TRUE;
	while ((a<narg || (nsets_in==1 && n==0)) && 
	       sscanf(line,"%lf%n",&dbl,&nchar)==1 && bTimeInRange) {
	  /* Use set=-1 as the time "set" */
	  if (sin) {
	    if (!bHaveT || (a>0))
	      set = sin;
	    else
	      set = -1;
	  } else {
	    if (!bHaveT)
	      set = a;
	    else
	      set = a-1;
	  }
	  if (set==-1 && ((bTB && dbl<tb) || (bTE && dbl>te)))
	    bTimeInRange = FALSE;
	    
	  if (bTimeInRange) {
	    if (n==0) {
	      if (nsets_in == 1)
		narg++;
	      if (set >= 0) {
		*nset = set+1;
		srenew(val,*nset);
		val[set] = NULL;
	      }
	    }
	    if (set == -1) {
	      if (sin == 0) {
		if (n % 100 == 0) 
		  srenew(*t,n+100);
		(*t)[n] = dbl;
	      }
	      /* else we should check the time of the next sets with set 0 */
	    } else {
	      if (n % 100 == 0) 
		srenew(val[set],n+100);
	      val[set][n] = (real)dbl;
	    }
	  }
	  a++;
	  line += nchar;
	}
	if (bTimeInRange) {
	  if (a == 0) {
	    fprintf(stderr,"Ignoring invalid line in %s:\n%s",fn,line0);
	  } else {
	    n++;
	    if (a != narg)
	      fprintf(stderr,"Invalid line in %s:\n%s"
		      "Using zeros for the last %d sets\n",
		      fn,line0,narg-a);
	  }
	}
	if (a > 0)
	  bFirstLine = FALSE;
      }
    }
    if (sin==0) {
      *nval = n;
      if (!bHaveT) {
	snew(*t,n);
	for(a=0; a<n; a++)
	  (*t)[a] = a;
      }
      if (n > 1)
	*dt = (real)((*t)[n-1]-(*t)[0])/(n-1.0);
      else
	*dt = 1;
    } else {
      if (n < *nval) {
	fprintf(stderr,"Set %d is shorter (%d) than the previous set (%d)\n",
		sin+1,n,*nval);
	*nval = n;
	fprintf(stderr,"Will use only the first %d points of every set\n",
		*nval);
      }
    }
  }
  fclose(fp);
  
  return val;
}

static void power_fit(int n,int nset,real **val,real *t)
{
  real *x,*y,quality,a,b;
  int  s,i;

  snew(x,n);
  snew(y,n);
  
  if (t[0]>0) {
    for(i=0; i<n; i++)
      if (t[0]>0)
	x[i] = log(t[i]);
  } else {
    fprintf(stdout,"First time is not larger than 0, using index number as time for power fit\n");
    for(i=0; i<n; i++)
      x[i] = log(i+1);
  }
  
  for(s=0; s<nset; s++) {
    i=0;
    for(i=0; i<n && val[s][i]>=0; i++)
      y[i] = log(val[s][i]);
    if (i < n)
      fprintf(stdout,"Will power fit up to point %d, since it is not larger than 0\n",i);
    quality = lsq_y_ax_b(i,x,y,&a,&b);
    fprintf(stdout,"Power fit set %3d:  error %.3f  a %g  b %g\n",
	    s+1,quality,a,exp(b));
  }
  
  sfree(y); 
  sfree(x);
}

static real cosine_content(int nhp,int n,real *y)
     /* Assumes n equidistant points */
{
  double fac,cosyint,yyint;
  int i;

  if (n < 2)
    return 0;
  
  fac = M_PI*nhp/(n-1);

  cosyint = 0;
  yyint = 0;
  for(i=0; i<n; i++) {
    cosyint += cos(fac*i)*y[i];
    yyint += y[i]*y[i];
  }
    
  return 2*cosyint*cosyint/(n*yyint);
}

static void plot_coscont(char *ccfile,int n,int nset,real **val)
{
  FILE *fp;
  int  s;
  real cc;
  
  fp = xvgropen(ccfile,"Cosine content","set / half periods","cosine content");
  
  for(s=0; s<nset; s++) {
    cc = cosine_content(s+1,n,val[s]);
    fprintf(fp," %d %g\n",s+1,cc);
    fprintf(stdout,"Cosine content of set %d with %.1f periods: %g\n",
	    s+1,0.5*(s+1),cc);
  }
  fprintf(stdout,"\n");
	    
  fclose(fp);
}

void histogram(char *distfile,real binwidth,int n, int nset, real **val)
{
  FILE *fp;
  int  i,s;
  real min,max;
  int  nbin;
  real *histo;

  min=val[0][0];
  max=val[0][0];
  for(s=0; s<nset; s++)
    for(i=0; i<n; i++)
      if (val[s][i] < min)
	min = val[s][i];
      else if (val[s][i] > max)
	max = val[s][i];
  
  if (-min > max)
    max = -min;
  nbin = (int)(max/binwidth)+1;
  fprintf(stderr,"Making distributions with %d bins\n",2*nbin+1);
  snew(histo,2*nbin+1);
  fp = xvgropen(distfile,"Distribution","","");
  for(s=0; s<nset; s++) {
    for(i=0; i<2*nbin+1; i++)
      histo[i] = 0;
    for(i=0; i<n; i++)
      histo[nbin+(int)(floor(val[s][i]/binwidth+0.5))]++;
    for(i=0; i<2*nbin+1; i++)
      fprintf(fp," %g  %g\n",(i-nbin)*binwidth,(real)histo[i]/(n*binwidth));
    if (s<nset-1)
      fprintf(fp,"&\n");
  }
  fclose(fp);
}

static int real_comp(const void *a,const void *b)
{
  real dif = *(real *)a - *(real *)b;

  if (dif < 0)
    return -1;
  else if (dif > 0)
    return 1;
  else
    return 0;
}

static void average(char *avfile,int avbar_opt,
		    int n, int nset,real **val,real *t)
{
  FILE   *fp;
  int    i,s,edge=0;
  double av,var,err;
  real   *tmp=NULL;
  
  fp = ffopen(avfile,"w");
  if ((avbar_opt == avbarERROR) && (nset == 1))
    avbar_opt = avbarNONE;
  if (avbar_opt != avbarNONE) {
    if (avbar_opt == avbar90) {
      snew(tmp,nset);
      fprintf(fp,"@TYPE xydydy\n");
      edge = (int)(nset*0.05+0.5);
      fprintf(stdout,"Errorbars: discarding %d points on both sides: %d%%"
	      " interval\n",edge,(int)(100*(nset-2*edge)/nset+0.5));
    } else
      fprintf(fp,"@TYPE xydy\n");
  }
  
  for(i=0; i<n; i++) {
    av = 0;
    for(s=0; s<nset; s++)
      av += val[s][i];
    av /= nset;
    fprintf(fp," %g %g",t[i],av);
    var = 0;
    if (avbar_opt != avbarNONE) {
      if (avbar_opt == avbar90) {
	for(s=0; s<nset; s++)
	  tmp[s] = val[s][i];
	qsort(tmp,nset,sizeof(tmp[0]),real_comp);
	fprintf(fp," %g %g",tmp[nset-1-edge]-av,av-tmp[edge]);
      } else {
	for(s=0; s<nset; s++)
	  var += sqr(val[s][i]-av);
	if (avbar_opt == avbarSTDDEV)
	  err = sqrt(var/nset);
	else
	  err = sqrt(var/(nset*(nset-1)));
	fprintf(fp," %g",err);
      }
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  if (avbar_opt == avbar90)
    sfree(tmp);
}

static real anal_ee_inf(real *parm,real T)
{
  return sqrt(parm[1]*2*parm[0]/T+parm[3]*2*parm[2]/T);
}

static real anal_ee(real *parm,real T,real t)
{
  real e1,e2;

  if (parm[0])
    e1 = exp(-t/parm[0]);
  else
    e1 = 1;
  if (parm[2])
    e2 = exp(-t/parm[2]);
  else
    e2 = 1;

  return sqrt(parm[1]*2*parm[0]/T*((e1 - 1)*parm[0]/t + 1) +
	      parm[3]*2*parm[2]/T*((e2 - 1)*parm[2]/t + 1));
}

static void estimate_error(char *eefile,int nb_min,int resol,int n,int nset,
			   double *av,double *sig,real **val,real dt,
			   bool bFitAc,bool bSingleExpFit,bool bAllowNegLTCorr)
{
  FILE   *fp;
  int    bs,prev_bs,nbs,nb;
  real   spacing,nbr;
  int    s,i,j;
  double blav,var;
  char   **leg;
  real   *tbs,*ybs,rtmp,*fitsig,fitparm[4];

  fp = xvgropen(eefile,"Error estimates","Block size (time)","Error estimate");
  fprintf(fp,
	  "@ subtitle \"using block averaging, total time %g (%d points)\"\n",
	  (n-1)*dt,n);
  snew(leg,2*nset);
  xvgr_legend(fp,2*nset,leg);
  sfree(leg);

  spacing = pow(2,1.0/resol);
  snew(tbs,n);
  snew(ybs,n);
  snew(fitsig,n);
  for(s=0; s<nset; s++) {
    nbs = 0;
    prev_bs = 0;
    nbr = nb_min;
    while (nbr <= n) {
      bs = n/(int)nbr;
      if (bs != prev_bs) {
	nb = n/bs;
	var = 0;
	for(i=0; i<nb; i++) {
	  blav=0;
	  for (j=0; j<bs; j++)
	    blav += val[s][bs*i+j];
	  var += sqr(av[s] - blav/bs);
	}
	tbs[nbs] = bs*dt;
	ybs[nbs] = sqrt(var/(nb*(nb-1.0))*(n*dt))/sig[s];
	nbs++;
      }
      nbr *= spacing;
      nb = (int)(nbr+0.5);
      prev_bs = bs;
    }

    for(i=0; i<nbs/2; i++) {
      rtmp         = tbs[i];
      tbs[i]       = tbs[nbs-1-i];
      tbs[nbs-1-i] = rtmp;
      rtmp         = ybs[i];
      ybs[i]       = ybs[nbs-1-i];
      ybs[nbs-1-i] = rtmp;
    }
    for(i=0; i<nbs; i++)
      fitsig[i] = sqrt(tbs[i]);

    if (!bSingleExpFit) {
      fitparm[0] = 0.002*n*dt;
      fitparm[1] = 0.95;
      fitparm[2] = 0.2*n*dt;
      do_lmfit(nbs,ybs,fitsig,0,tbs,0,dt*n,bDebugMode(),effnERREST,fitparm,0);
      fitparm[3] = 1-fitparm[1];
    }
    if (bSingleExpFit || fitparm[0]<0 || fitparm[2]<0 || fitparm[1]<0
	|| (fitparm[1]>1 && !bAllowNegLTCorr) || fitparm[2]>(n-1)*dt) {
      if (!bSingleExpFit) {
	if (fitparm[2]>(n-1)*dt)
	  fprintf(stdout,
		  "Warning: tau2 is longer than the length of the data (%g)\n"
		  "         the statistics might be bad\n",
		  (n-1)*dt);
	else
	  fprintf(stdout,"a fitted parameter is negative\n");
	fprintf(stdout,"invalid fit:  e.e. %g  a %g  tau1 %g  tau2 %g\n",
		sig[s]*anal_ee_inf(fitparm,n*dt),
		fitparm[1],fitparm[0],fitparm[2]);
	fprintf(stderr,"Will use a single exponential fit for set %d\n",s+1);
      }
      fitparm[0] = n*dt*0.002;
      fitparm[1] = 1;
      fitparm[2] = 0;
      do_lmfit(nbs,ybs,fitsig,0,tbs,0,dt*n,bDebugMode(),effnERREST,fitparm,6);
      fitparm[3] = 1-fitparm[1];
    }
    fprintf(stdout,"Set %3d:  err.est. %g  a %g  tau1 %g  tau2 %g\n",
	    s+1,sig[s]*anal_ee_inf(fitparm,n*dt),
	    fitparm[1],fitparm[0],fitparm[2]);
    fprintf(fp,"@ legend string %d \"av %f\"\n",2*s,av[s]);
    fprintf(fp,"@ legend string %d \"ee %6g\"\n",
	    2*s+1,sig[s]*anal_ee_inf(fitparm,n*dt));
    for(i=0; i<nbs; i++)
      fprintf(fp,"%g %g %g\n",tbs[i],sig[s]/sqrt(n*dt)*ybs[i],
	      sig[s]/sqrt(n*dt)*fit_function(effnERREST,fitparm,tbs[i]));

    if (bFitAc) {
      real *ac,ac_fit[4];
      
      snew(ac,n);
      for(i=0; i<n; i++) {
	ac[i] = val[s][i] - av[s];
	if (i > 0)
	  fitsig[i] = sqrt(i);
	else
	  fitsig[i] = 1;
      }
      low_do_autocorr(NULL,NULL,n,1,-1,&ac,
		      dt,eacNormal,1,FALSE,TRUE,TRUE,
		      FALSE,0,0,
		      effnEXP3,0);
      
      ac_fit[0] = 0.002*n*dt;
      ac_fit[1] = 0.95;
      ac_fit[2] = 0.2*n*dt;
      do_lmfit(n/nb_min,ac,fitsig,dt,0,0,dt*n/nb_min,
              bDebugMode(),effnEXP3,ac_fit,0);
      ac_fit[3] = 1 - ac_fit[1];

      fprintf(stdout,"Set %3d:  ac erest %g  a %g  tau1 %g  tau2 %g\n",
	    s+1,sig[s]*anal_ee_inf(ac_fit,n*dt),
	    ac_fit[1],ac_fit[0],ac_fit[2]);

      fprintf(fp,"&\n");
      for(i=0; i<nbs; i++)
	fprintf(fp,"%g %g\n",tbs[i],
		sig[s]/sqrt(n*dt)*fit_function(effnERREST,ac_fit,tbs[i]));

      sfree(ac);
    }
    if (s < nset-1)
      fprintf(fp,"&\n");
  }
  sfree(fitsig);
  sfree(ybs);
  sfree(tbs);
  fclose(fp);
}

static void filter(real flen,int n,int nset,real **val,real dt)
{
  int    f,s,i,j;
  double *filt,sum,vf,fluc,fluctot;

  f = (int)(flen/(2*dt));
  snew(filt,f+1);
  filt[0] = 1;
  sum = 1;
  for(i=1; i<=f; i++) {
    filt[i] = cos(M_PI*dt*i/flen);
    sum += 2*filt[i];
  }
  for(i=0; i<=f; i++)
    filt[i] /= sum;
  fprintf(stdout,"Will calculate the fluctuation over %d points\n",n-2*f);
  fprintf(stdout,"  using a filter of length %g of %d points\n",flen,2*f+1);
  fluctot = 0;
  for(s=0; s<nset; s++) {
    fluc = 0;
    for(i=f; i<n-f; i++) {
      vf = filt[0]*val[s][i];
      for(j=1; j<=f; j++)
	vf += filt[j]*(val[s][i-f]+val[s][i+f]);
      fluc += sqr(val[s][i] - vf);
    }
    fluc /= n - 2*f;
    fluctot += fluc;
    fprintf(stdout,"Set %3d filtered fluctuation: %12.6e\n",s+1,sqrt(fluc));
  }
  fprintf(stdout,"Overall filtered fluctuation: %12.6e\n",sqrt(fluctot/nset));
  fprintf(stdout,"\n");

  sfree(filt);
}

static void do_fit(char *fn,int nx,int ny,real *x0,real **val)
{
  FILE *out;
  real *c1=NULL,*sig=NULL,*fitparm;
  real dt=0,tendfit,tbeginfit,integral;
  int  i,efitfn,nparm;
  
  out = ffopen(fn,"w");
  efitfn = get_acffitfn();
  nparm  = nfp_ffn[efitfn];
  fprintf(out,"Will fit to the following function:\n");
  fprintf(out,"%s\n",longs_ffn[efitfn]);
  switch (nx) {
  case 1:
    c1 = val[0];
    snew(sig,ny);
    fprintf(out,"Using only column as y values\n");
    break;
  case 2:
    c1  = val[0];
    sig = val[1];
    fprintf(out,"Using two columns as y and sigma values\n");
    break;
  default:
    gmx_fatal(FARGS,"Dont know what to do with %d columns of data",nx);
  }
  tbeginfit = x0 ? x0[0]    : 0;
  tendfit   = x0 ? x0[ny-1] : (ny-1)*dt;
  
  snew(fitparm,nparm);
  switch(efitfn) {
  case effnEXP1:
    fitparm[0] = 0.5;
    break;
  case effnEXP2:
    fitparm[0] = 0.5;
    fitparm[1] = c1[0];
    break;
  case effnEXP3:
    fitparm[0] = 1.0;
    fitparm[1] = 0.5*c1[0];
    fitparm[2] = 10.0;
    break;
  case effnEXP5:
    fitparm[0] = fitparm[2] = 0.5*c1[0];
    fitparm[1] = 10;
    fitparm[3] = 40;
    fitparm[4] = 0;
    break;
  case effnEXP7:
    fitparm[0] = fitparm[2] = fitparm[4] = 0.33*c1[0];
    fitparm[1] = 1;
    fitparm[3] = 10;
    fitparm[5] = 100;
    fitparm[6] = 0;
    break;
  case effnEXP9:
    fitparm[0] = fitparm[2] = fitparm[4] = fitparm[6] = 0.25*c1[0];
    fitparm[1] = 0.1;
    fitparm[3] = 1;
    fitparm[5] = 10;
    fitparm[7] = 100;
    fitparm[8] = 0;
    break;
  default:
    fprintf(out,"Warning: don't know how to initialize the parameters\n");
    for(i=0; (i<nparm); i++)
      fitparm[i] = 1;
  }
  fprintf(out,"Starting parameters:\n");
  for(i=0; (i<nparm); i++) 
    fprintf(out,"a%-2d = %12.5e\n",i+1,fitparm[i]);
  if (do_lmfit(ny,c1,sig,dt,x0,tbeginfit,tendfit,
	       bDebugMode(),efitfn,fitparm,0)) {
    for(i=0; (i<nparm); i++) 
      fprintf(out,"a%-2d = %12.5e\n",i+1,fitparm[i]);
  }
  else {
    fprintf(out,"No solution was found\n");
  }
  fclose(out);
}

int gmx_analyze(int argc,char *argv[])
{
  static char *desc[] = {
    "g_analyze reads an ascii file and analyzes data sets.",
    "A line in the input file may start with a time",
    "(see option [TT]-time[tt]) and any number of y values may follow.",
    "Multiple sets can also be",
    "read when they are seperated by & (option [TT]-n[tt]),",
    "in this case only one y value is read from each line.",
    "All lines starting with # and @ are skipped.",
    "All analyses can also be done for the derivative of a set",
    "(option [TT]-d[tt]).[PAR]",

    "All options, except for [TT]-av[tt] and [TT]-power[tt] assume that the",
    "points are equidistant in time.[PAR]",

    "g_analyze always shows the average and standard deviation of each",
    "set. For each set it also shows the relative deviation of the third",
    "and forth cumulant from those of a Gaussian distribution with the same",
    "standard deviation.[PAR]",

    "Option [TT]-ac[tt] produces the autocorrelation function(s).[PAR]",
    
    "Option [TT]-cc[tt] plots the resemblance of set i with a cosine of",
    "i/2 periods. The formula is:[BR]"
    "2 (int0-T y(t) cos(pi t/i) dt)^2 / int0-T y(t) y(t) dt[BR]",
    "This is useful for principal components obtained from covariance",
    "analysis, since the principal components of random diffusion are",
    "pure cosines.[PAR]",
    
    "Option [TT]-msd[tt] produces the mean square displacement(s).[PAR]",
    
    "Option [TT]-dist[tt] produces distribution plot(s).[PAR]",
    
    "Option [TT]-av[tt] produces the average over the sets.",
    "Error bars can be added with the option [TT]-errbar[tt].",
    "The errorbars can represent the standard deviation, the error",
    "(assuming the points are independent) or the interval containing",
    "90% of the points, by discarding 5% of the points at the top and",
    "the bottom.[PAR]",
    
    "Option [TT]-ee[tt] produces error estimates using block averaging.",
    "A set is divided in a number of blocks and averages are calculated for",
    "each block. The error for the total average is calculated from",
    "the variance between averages of the m blocks B_i as follows:",
    "error^2 = Sum (B_i - <B>)^2 / (m*(m-1)).",
    "These errors are plotted as a function of the block size.",
    "Also an analytical block average curve is plotted, assuming",
    "that the autocorrelation is a sum of two exponentials.",
    "The analytical curve for the block average BA is:[BR]",
    "BA(t) = sigma sqrt(2/T (  a   (tau1 ((exp(-t/tau1) - 1) tau1/t + 1)) +[BR]",
    "                        (1-a) (tau2 ((exp(-t/tau2) - 1) tau2/t + 1)))),[BR]"
    "where T is the total time.",
    "a, tau1 and tau2 are obtained by fitting BA(t) to the calculated block",
    "average.",
    "When the actual block average is very close to the analytical curve,",
    "the error is sigma*sqrt(2/T (a tau1 + (1-a) tau2)).",
    "The complete derivation is given in",
    "B. Hess, J. Chem. Phys. 116:209-217, 2002.[PAR]",

    "Option [TT]-filter[tt] prints the RMS high-frequency fluctuation",
    "of each set and over all sets with respect to a filtered average.",
    "The filter is proportional to cos(pi t/len) where t goes from -len/2",
    "to len/2. len is supplied with the option [TT]-filter[tt].",
    "This filter reduces oscillations with period len/2 and len by a factor",
    "of 0.79 and 0.33 respectively.[PAR]",
    
    "Option [TT]-power[tt] fits the data to b t^a, which is accomplished",
    "by fitting to a t + b on log-log scale. All points after the first",
    "zero or negative value are ignored."
  };
  static real tb=-1,te=-1,frac=0.5,filtlen=0,binwidth=0.1,aver_start=0;
  static bool bHaveT=TRUE,bDer=FALSE,bSubAv=TRUE,bAverCorr=FALSE,bXYdy=FALSE;
  static bool bEESEF=FALSE,bEENLC=FALSE,bEeFitAc=FALSE,bPower=FALSE;
  static bool bIntegrate=FALSE; 
  static int  linelen=4096,nsets_in=1,d=1,nb_min=4,resol=10;

  /* must correspond to enum avbar* declared at beginning of file */
  static char *avbar_opt[avbarNR+1] = { 
    NULL, "none", "stddev", "error", "90", NULL
  };

  t_pargs pa[] = {
    { "-linelen", FALSE, etINT, {&linelen},
      "HIDDENMaximum input line length" },
    { "-time",    FALSE, etBOOL, {&bHaveT},
      "Expect a time in the input" },
    { "-b",       FALSE, etREAL, {&tb},
      "First time to read from set" },
    { "-e",       FALSE, etREAL, {&te},
      "Last time to read from set" },
    { "-n",       FALSE, etINT, {&nsets_in},
      "Read # sets seperated by &" },
    { "-d",       FALSE, etBOOL, {&bDer},
	"Use the derivative" },
    { "-dp",      FALSE, etINT, {&d}, 
      "HIDDENThe derivative is the difference over # points" },
    { "-bw",      FALSE, etREAL, {&binwidth},
      "Binwidth for the distribution" },
    { "-errbar",  FALSE, etENUM, {avbar_opt},
      "Error bars for -av" },
    { "-integrate",FALSE,etBOOL, {&bIntegrate},
      "Integrate data function(s) numerically using trapezium rule" },
    { "-aver_start",FALSE, etREAL, {&aver_start},
      "Start averaging the integral from here" },
    { "-xydy",    FALSE, etBOOL, {&bXYdy},
      "Interpret second data set as error in the y values for integrating" },
    { "-nbmin",   FALSE, etINT, {&nb_min},
      "HIDDENMinimum number of blocks for block averaging" },
    { "-resol", FALSE, etINT, {&resol},
      "HIDDENResolution for the block averaging, block size increases with"
    " a factor 2^(1/#)" },
    { "-eeexpfit", FALSE, etBOOL, {&bEESEF},
      "HIDDENAlways use a single exponential fit for the error estimate" },
    { "-eenlc", FALSE, etBOOL, {&bEENLC},
      "HIDDENAllow a negative long-time correlation" },
    { "-eefitac", FALSE, etBOOL, {&bEeFitAc},
      "HIDDENAlso plot analytical block average using a autocorrelation fit" },
    { "-filter",  FALSE, etREAL, {&filtlen},
      "Print the high-frequency fluctuation after filtering with a cosine filter of length #" },
    { "-power", FALSE, etBOOL, {&bPower},
      "Fit data to: b t^a" },
    { "-subav", FALSE, etBOOL, {&bSubAv},
      "Subtract the average before autocorrelating" },
    { "-oneacf", FALSE, etBOOL, {&bAverCorr},
      "Calculate one ACF over all sets" }
  };
#define NPA asize(pa)

  FILE     *out;
  int      n,nlast,s,nset,i,j=0;
  real     **val,*t,dt,tot,error;
  double   *av,*sig,cum1,cum2,cum3,cum4,db;
  char     *acfile,*msdfile,*ccfile,*distfile,*avfile,*eefile,*fitfile;
  
  t_filenm fnm[] = { 
    { efXVG, "-f",    "graph",    ffREAD   },
    { efXVG, "-ac",   "autocorr", ffOPTWR  },
    { efXVG, "-msd",  "msd",      ffOPTWR  },
    { efXVG, "-cc",   "coscont",  ffOPTWR  },
    { efXVG, "-dist", "distr",    ffOPTWR  },
    { efXVG, "-av",   "average",  ffOPTWR  },
    { efXVG, "-ee",   "errest",   ffOPTWR  },
    { efLOG, "-g",    "fitlog",   ffOPTWR  }
  }; 
#define NFILE asize(fnm) 

  int     npargs;
  t_pargs *ppa;

  npargs = asize(pa); 
  ppa    = add_acf_pargs(&npargs,pa);
  
  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_BE_NICE ,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL); 

  acfile   = opt2fn_null("-ac",NFILE,fnm);
  msdfile  = opt2fn_null("-msd",NFILE,fnm);
  ccfile   = opt2fn_null("-cc",NFILE,fnm);
  distfile = opt2fn_null("-dist",NFILE,fnm);
  avfile   = opt2fn_null("-av",NFILE,fnm);
  eefile   = opt2fn_null("-ee",NFILE,fnm);
  if (opt2parg_bSet("-fitfn",npargs,ppa)) 
    fitfile  = opt2fn("-g",NFILE,fnm);
  else
    fitfile  = opt2fn_null("-g",NFILE,fnm);
    
  val=read_val(opt2fn("-f",NFILE,fnm),bHaveT,
	       opt2parg_bSet("-b",npargs,ppa),tb,
	       opt2parg_bSet("-e",npargs,ppa),te,
	       nsets_in,&nset,&n,&dt,&t,linelen);
  printf("Read %d sets of %d points, dt = %g\n\n",nset,n,dt);
  
  if (bDer) {
    printf("Calculating the derivative as (f[i+%d]-f[i])/(%d*dt)\n\n",
	    d,d);
    n -= d;
    for(s=0; s<nset; s++)
      for(i=0; (i<n); i++)
	val[s][i] = (val[s][i+d]-val[s][i])/(d*dt);
  }
  if (bIntegrate) {
    real sum,stddev;
    printf("Calculating the integral using the trapezium rule\n");
    
    if (bXYdy) {
      sum = evaluate_integral(n,dt,val[0],val[1],aver_start,&stddev);
      printf("Integral %10.3f +/- %10.5f\n",sum,stddev);
    }
    else {
      for(s=0; s<nset; s++) {
	sum = evaluate_integral(n,dt,val[s],NULL,aver_start,&stddev);
	printf("Integral %d  %10.5f  +/- %10.5f\n",s+1,sum,stddev);
      }
    }
  }

  
  if (fitfile)
    do_fit(fitfile,nset,n,t,val);

  printf("                                      std. dev.    relative deviation of\n");
  printf("                       standard       ---------   cumulants from those of\n");
  printf("set      average       deviation      sqrt(n-1)   a Gaussian distribition\n");
  printf("                                                      cum. 3   cum. 4\n");
  snew(av,nset);
  snew(sig,nset);
  for(s=0; (s<nset); s++) {
    cum1 = 0;
    cum2 = 0;
    cum3 = 0;
    cum4 = 0;
    for(i=0; (i<n); i++)
      cum1 += val[s][i];
    cum1 /= n;
    for(i=0; (i<n); i++) {
      db = val[s][i]-cum1;
      cum2 += db*db;
      cum3 += db*db*db;
      cum4 += db*db*db*db;
    }
    cum2  /= n;
    cum3  /= n;
    cum4  /= n;
    av[s]  = cum1;
    sig[s] = sqrt(cum2);
    if (n > 1)
      error = sqrt(cum2/(n-1));
    else
      error = 0;
    printf("SS%d  %13.6e   %12.6e   %12.6e      %6.3f   %6.3f\n",
	   s+1,av[s],sig[s],error,
	   sig[s] ? cum3/(sig[s]*sig[s]*sig[s]*sqrt(8/M_PI)) : 0,
	   sig[s] ? cum4/(sig[s]*sig[s]*sig[s]*sig[s]*3)-1 : 0); 
  }
  printf("\n");

  if (filtlen)
    filter(filtlen,n,nset,val,dt);
  
  if (msdfile) {
    out=xvgropen(msdfile,"Mean square displacement",
		 "time","MSD (nm\\S2\\N)");
    nlast = (int)(n*frac);
    for(s=0; s<nset; s++) {
      for(j=0; j<=nlast; j++) {
	if (j % 100 == 0)
	  fprintf(stderr,"\r%d",j);
	tot=0;
	for(i=0; i<n-j; i++)
	  tot += sqr(val[s][i]-val[s][i+j]); 
	tot /= (real)(n-j);
	fprintf(out," %g %8g\n",dt*j,tot);
      }
      if (s<nset-1)
	fprintf(out,"&\n");
    }
    fclose(out);
    fprintf(stderr,"\r%d, time=%g\n",j-1,(j-1)*dt);
  }
  if (ccfile)
    plot_coscont(ccfile,n,nset,val);
  
  if (distfile)
    histogram(distfile,binwidth,n,nset,val);
  if (avfile)
    average(avfile,nenum(avbar_opt),n,nset,val,t);
  if (eefile)
    estimate_error(eefile,nb_min,resol,n,nset,av,sig,val,dt,
		   bEeFitAc,bEESEF,bEENLC);
  if (bPower)
    power_fit(n,nset,val,t);
  if (acfile) {
    if (bSubAv) 
      for(s=0; s<nset; s++)
	for(i=0; i<n; i++)
	  val[s][i] -= av[s];
    do_autocorr(acfile,"Autocorrelation",n,nset,val,dt,
		eacNormal,bAverCorr);
  }
  
  view_all(NFILE, fnm);
  
  thanx(stderr);

  return 0;
}
  
