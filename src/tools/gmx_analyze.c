/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "gmx_fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "readinp.h"
#include "statutil.h"
#include "txtdump.h"
#include "gstat.h"
#include "gmx_matrix.h"
#include "gmx_statistics.h"
#include "xvgr.h"
#include "gmx_ana.h"
#include "geminate.h"

/* must correspond to char *avbar_opt[] declared in main() */
enum { avbarSEL, avbarNONE, avbarSTDDEV, avbarERROR, avbar90, avbarNR };

static void power_fit(int n,int nset,real **val,real *t)
{
  real *x,*y,quality,a,b,r;
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
    lsq_y_ax_b(i,x,y,&a,&b,&r,&quality);
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

static void plot_coscont(const char *ccfile,int n,int nset,real **val,
                         const output_env_t oenv)
{
  FILE *fp;
  int  s;
  real cc;
  
  fp = xvgropen(ccfile,"Cosine content","set / half periods","cosine content",
                oenv);
  
  for(s=0; s<nset; s++) {
    cc = cosine_content(s+1,n,val[s]);
    fprintf(fp," %d %g\n",s+1,cc);
    fprintf(stdout,"Cosine content of set %d with %.1f periods: %g\n",
	    s+1,0.5*(s+1),cc);
  }
  fprintf(stdout,"\n");
	    
  ffclose(fp);
}

static void regression_analysis(int n,gmx_bool bXYdy,
                                real *x,int nset,real **val)
{
  real S,chi2,a,b,da,db,r=0;

  if (bXYdy || (nset == 1)) 
  {
      printf("Fitting data to a function f(x) = ax + b\n");
      printf("Minimizing residual chi2 = Sum_i w_i [f(x_i) - y_i]2\n");
      printf("Error estimates will be given if w_i (sigma) values are given\n");
      printf("(use option -xydy).\n\n");
      if (bXYdy) 
          lsq_y_ax_b_error(n,x,val[0],val[1],&a,&b,&da,&db,&r,&S);
      else
          lsq_y_ax_b(n,x,val[0],&a,&b,&r,&S);
      chi2 = sqr((n-2)*S);
      printf("Chi2                    = %g\n",chi2);
      printf("S (Sqrt(Chi2/(n-2))     = %g\n",S);
      printf("Correlation coefficient = %.1f%%\n",100*r);
      printf("\n");
      if (bXYdy) {
          printf("a    = %g +/- %g\n",a,da);
          printf("b    = %g +/- %g\n",b,db);
      }
      else {
          printf("a    = %g\n",a);
          printf("b    = %g\n",b);
      }
  }
  else 
  {
      double chi2,*a,**xx,*y;
      int i,j;
      
      snew(y,n);
      snew(xx,nset-1);
      for(j=0; (j<nset-1); j++)
          snew(xx[j],n);
      for(i=0; (i<n); i++)
      {
          y[i] = val[0][i];
          for(j=1; (j<nset); j++)
              xx[j-1][i] = val[j][i];
      }
      snew(a,nset-1);
      chi2 = multi_regression(NULL,n,y,nset-1,xx,a);
      printf("Fitting %d data points in %d sets\n",n,nset-1);
      printf("chi2 = %g\n",chi2);
      printf("A =");
      for(i=0; (i<nset-1); i++)
      {
          printf("  %g",a[i]);
          sfree(xx[i]);
      }
      printf("\n");
      sfree(xx);
      sfree(y);
      sfree(a);
  }
}

void histogram(const char *distfile,real binwidth,int n, int nset, real **val,
               const output_env_t oenv)
{
  FILE *fp;
  int  i,s;
  double min,max;
  int  nbin;
#if (defined SIZEOF_LONG_LONG_INT) && (SIZEOF_LONG_LONG_INT >= 8)    
  long long int *histo;
#else
  double        *histo;
#endif

  min = val[0][0];
  max = val[0][0];
  for(s=0; s<nset; s++)
    for(i=0; i<n; i++)
      if (val[s][i] < min)
	min = val[s][i];
      else if (val[s][i] > max)
	max = val[s][i];
  
  min = binwidth*floor(min/binwidth);
  max = binwidth*ceil(max/binwidth);
  if (min != 0)
    min -= binwidth;
  max += binwidth;

  nbin = (int)((max - min)/binwidth + 0.5) + 1;
  fprintf(stderr,"Making distributions with %d bins\n",nbin);
  snew(histo,nbin);
  fp = xvgropen(distfile,"Distribution","","",oenv);
  for(s=0; s<nset; s++) {
    for(i=0; i<nbin; i++)
      histo[i] = 0;
    for(i=0; i<n; i++)
      histo[(int)((val[s][i] - min)/binwidth + 0.5)]++;
    for(i=0; i<nbin; i++)
      fprintf(fp," %g  %g\n",min+i*binwidth,(double)histo[i]/(n*binwidth));
    if (s < nset-1)
      fprintf(fp,"&\n");
  }
  ffclose(fp);
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

static void average(const char *avfile,int avbar_opt,
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
  ffclose(fp);
  
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

static void estimate_error(const char *eefile,int nb_min,int resol,int n,
                           int nset, double *av,double *sig,real **val,real dt,
                           gmx_bool bFitAc,gmx_bool bSingleExpFit,gmx_bool bAllowNegLTCorr,
                           const output_env_t oenv)
{
    FILE   *fp;
    int    bs,prev_bs,nbs,nb;
    real   spacing,nbr;
    int    s,i,j;
    double blav,var;
    char   **leg;
    real   *tbs,*ybs,rtmp,dens,*fitsig,twooe,tau1_est,tau_sig;
    real   fitparm[4];
    real   ee,a,tau1,tau2;
    
    if (n < 4)
    {
      fprintf(stdout,"The number of points is smaller than 4, can not make an error estimate\n");
      
      return;
    }
    
    fp = xvgropen(eefile,"Error estimates",
                  "Block size (time)","Error estimate", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp,
                "@ subtitle \"using block averaging, total time %g (%d points)\"\n",
                (n-1)*dt,n);
    }
    snew(leg,2*nset);
    xvgr_legend(fp,2*nset,(const char**)leg,oenv);
    sfree(leg);

    spacing = pow(2,1.0/resol);
    snew(tbs,n);
    snew(ybs,n);
    snew(fitsig,n);
    for(s=0; s<nset; s++)
    {
        nbs = 0;
        prev_bs = 0;
        nbr = nb_min;
        while (nbr <= n)
        {
            bs = n/(int)nbr;
            if (bs != prev_bs)
            {
                nb = n/bs;
                var = 0;
                for(i=0; i<nb; i++)
                {
                    blav=0;
                    for (j=0; j<bs; j++)
                    {
                        blav += val[s][bs*i+j];
                    }
                    var += sqr(av[s] - blav/bs);
                }
                tbs[nbs] = bs*dt;
                if (sig[s] == 0)
                {
                    ybs[nbs] = 0;
                }
                else
                {
                    ybs[nbs] = var/(nb*(nb-1.0))*(n*dt)/(sig[s]*sig[s]);
                }
                nbs++;
            }
            nbr *= spacing;
            nb = (int)(nbr+0.5);
            prev_bs = bs;
        }
        if (sig[s] == 0)
        {
            ee   = 0;
            a    = 1;
            tau1 = 0;
            tau2 = 0;
        }
        else
        {
            for(i=0; i<nbs/2; i++)
            {
                rtmp         = tbs[i];
                tbs[i]       = tbs[nbs-1-i];
                tbs[nbs-1-i] = rtmp;
                rtmp         = ybs[i];
                ybs[i]       = ybs[nbs-1-i];
                ybs[nbs-1-i] = rtmp;
            }
            /* The initial slope of the normalized ybs^2 is 1.
             * For a single exponential autocorrelation: ybs(tau1) = 2/e tau1
             * From this we take our initial guess for tau1.
             */
            twooe = 2/exp(1);
            i = -1;
            do
            {
                i++;
                tau1_est = tbs[i];
            } while (i < nbs - 1 &&
                     (ybs[i] > ybs[i+1] || ybs[i] > twooe*tau1_est));
            
            if (ybs[0] > ybs[1])
            {
                fprintf(stdout,"Data set %d has strange time correlations:\n"
                        "the std. error using single points is larger than that of blocks of 2 points\n"
                        "The error estimate might be inaccurate, check the fit\n",
                        s+1);
                /* Use the total time as tau for the fitting weights */
                tau_sig = (n - 1)*dt;
            }
            else
            {
                tau_sig = tau1_est;
            }
            
            if (debug)
            {
                fprintf(debug,"set %d tau1 estimate %f\n",s+1,tau1_est);
            }
            
            /* Generate more or less appropriate sigma's,
             * also taking the density of points into account.
             */
            for(i=0; i<nbs; i++)
            {
                if (i == 0)
                {
                    dens = tbs[1]/tbs[0] - 1;
                }
                else if (i == nbs-1)
                {
                    dens = tbs[nbs-1]/tbs[nbs-2] - 1;
                }
                else
                {
                    dens = 0.5*(tbs[i+1]/tbs[i-1] - 1);
                }
                fitsig[i] = sqrt((tau_sig + tbs[i])/dens);
            }
            
            if (!bSingleExpFit)
            {
                fitparm[0] = tau1_est;
                fitparm[1] = 0.95;
                /* We set the initial guess for tau2
                 * to halfway between tau1_est and the total time (on log scale).
                 */
                fitparm[2] = sqrt(tau1_est*(n-1)*dt);
                do_lmfit(nbs,ybs,fitsig,0,tbs,0,dt*n,oenv,
                         bDebugMode(),effnERREST,fitparm,0);
                fitparm[3] = 1-fitparm[1];
            }
            if (bSingleExpFit || fitparm[0]<0 || fitparm[2]<0 || fitparm[1]<0
                || (fitparm[1]>1 && !bAllowNegLTCorr) || fitparm[2]>(n-1)*dt)
            {
                if (!bSingleExpFit)
                {
                    if (fitparm[2] > (n-1)*dt)
                    {
                        fprintf(stdout,
                                "Warning: tau2 is longer than the length of the data (%g)\n"
                                "         the statistics might be bad\n",
                                (n-1)*dt);
                    }
                    else
                    {
                        fprintf(stdout,"a fitted parameter is negative\n");
                    }
                    fprintf(stdout,"invalid fit:  e.e. %g  a %g  tau1 %g  tau2 %g\n",
                            sig[s]*anal_ee_inf(fitparm,n*dt),
                            fitparm[1],fitparm[0],fitparm[2]);
                    /* Do a fit with tau2 fixed at the total time.
                     * One could also choose any other large value for tau2.
                     */
                    fitparm[0] = tau1_est;
                    fitparm[1] = 0.95;
                    fitparm[2] = (n-1)*dt;
                    fprintf(stderr,"Will fix tau2 at the total time: %g\n",fitparm[2]);
                    do_lmfit(nbs,ybs,fitsig,0,tbs,0,dt*n,oenv,bDebugMode(),
                             effnERREST,fitparm,4);
                    fitparm[3] = 1-fitparm[1];
                }
                if (bSingleExpFit || fitparm[0]<0 || fitparm[1]<0
                    || (fitparm[1]>1 && !bAllowNegLTCorr))
                {
                    if (!bSingleExpFit) {
                        fprintf(stdout,"a fitted parameter is negative\n");
                        fprintf(stdout,"invalid fit:  e.e. %g  a %g  tau1 %g  tau2 %g\n",
                                sig[s]*anal_ee_inf(fitparm,n*dt),
                                fitparm[1],fitparm[0],fitparm[2]);
                    }
                    /* Do a single exponential fit */
                    fprintf(stderr,"Will use a single exponential fit for set %d\n",s+1);
                    fitparm[0] = tau1_est;
                    fitparm[1] = 1.0;
                    fitparm[2] = 0.0;
                    do_lmfit(nbs,ybs,fitsig,0,tbs,0,dt*n,oenv,bDebugMode(),
                             effnERREST,fitparm,6);
                    fitparm[3] = 1-fitparm[1];
                }
            }
            ee   = sig[s]*anal_ee_inf(fitparm,n*dt);
            a    = fitparm[1];
            tau1 = fitparm[0];
            tau2 = fitparm[2];
        }
        fprintf(stdout,"Set %3d:  err.est. %g  a %g  tau1 %g  tau2 %g\n",
                s+1,ee,a,tau1,tau2);
        fprintf(fp,"@ legend string %d \"av %f\"\n",2*s,av[s]);
        fprintf(fp,"@ legend string %d \"ee %6g\"\n",
                2*s+1,sig[s]*anal_ee_inf(fitparm,n*dt));
        for(i=0; i<nbs; i++)
        {
            fprintf(fp,"%g %g %g\n",tbs[i],sig[s]*sqrt(ybs[i]/(n*dt)),
                    sig[s]*sqrt(fit_function(effnERREST,fitparm,tbs[i])/(n*dt)));
        }
        
        if (bFitAc)
        {
            int fitlen;
            real *ac,acint,ac_fit[4];
            
            snew(ac,n);
            for(i=0; i<n; i++) {
                ac[i] = val[s][i] - av[s];
                if (i > 0)
                    fitsig[i] = sqrt(i);
                else
                    fitsig[i] = 1;
            }
            low_do_autocorr(NULL,oenv,NULL,n,1,-1,&ac,
                            dt,eacNormal,1,FALSE,TRUE,
                            FALSE,0,0,effnNONE,0);
            
            fitlen = n/nb_min;
            
            /* Integrate ACF only up to fitlen/2 to avoid integrating noise */ 
            acint = 0.5*ac[0];
            for(i=1; i<=fitlen/2; i++)
            {
                acint += ac[i];
            }
            acint *= dt;
            
            /* Generate more or less appropriate sigma's */
            for(i=0; i<=fitlen; i++)
            {
                fitsig[i] = sqrt(acint + dt*i);
            }
            
            ac_fit[0] = 0.5*acint;
            ac_fit[1] = 0.95;
            ac_fit[2] = 10*acint;
            do_lmfit(n/nb_min,ac,fitsig,dt,0,0,fitlen*dt,oenv,
                     bDebugMode(),effnEXP3,ac_fit,0);
            ac_fit[3] = 1 - ac_fit[1];
            
            fprintf(stdout,"Set %3d:  ac erest %g  a %g  tau1 %g  tau2 %g\n",
                    s+1,sig[s]*anal_ee_inf(ac_fit,n*dt),
                    ac_fit[1],ac_fit[0],ac_fit[2]);
            
            fprintf(fp,"&\n");
            for(i=0; i<nbs; i++)
            {
                fprintf(fp,"%g %g\n",tbs[i],
                        sig[s]*sqrt(fit_function(effnERREST,ac_fit,tbs[i]))/(n*dt));
            }
            
            sfree(ac);
        }
        if (s < nset-1)
        {
            fprintf(fp,"&\n");
        }
    }
    sfree(fitsig);
    sfree(ybs);
    sfree(tbs);
    ffclose(fp);
}

static void luzar_correl(int nn,real *time,int nset,real **val,real temp,
			 gmx_bool bError,real fit_start,real smooth_tail_start,
                         const output_env_t oenv)
{
  const real tol = 1e-8;
  real *kt;
  real weight,d2;
  int  j;
  
  please_cite(stdout,"Spoel2006b");
  
  /* Compute negative derivative k(t) = -dc(t)/dt */
  if (!bError) {
    snew(kt,nn);
    compute_derivative(nn,time,val[0],kt);
    for(j=0; (j<nn); j++)
      kt[j] = -kt[j];
    if (debug) {
      d2 = 0;
      for(j=0; (j<nn); j++)
	d2 += sqr(kt[j] - val[3][j]);
      fprintf(debug,"RMS difference in derivatives is %g\n",sqrt(d2/nn));
    }
    analyse_corr(nn,time,val[0],val[2],kt,NULL,NULL,NULL,fit_start,
		 temp,smooth_tail_start,oenv);
    sfree(kt);
  }
  else if (nset == 6) {
    analyse_corr(nn,time,val[0],val[2],val[4],
		 val[1],val[3],val[5],fit_start,temp,smooth_tail_start,oenv);
  }
  else {
    printf("Inconsistent input. I need c(t) sigma_c(t) n(t) sigma_n(t) K(t) sigma_K(t)\n");
    printf("Not doing anything. Sorry.\n");
  }
}

static void filter(real flen,int n,int nset,real **val,real dt,
                   const output_env_t oenv)
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

static void do_fit(FILE *out,int n,gmx_bool bYdy,
                   int ny,real *x0,real **val,
                   int npargs,t_pargs *ppa,const output_env_t oenv)
{
  real *c1=NULL,*sig=NULL,*fitparm;
  real tendfit,tbeginfit;
  int  i,efitfn,nparm;
  
  efitfn = get_acffitfn();
  nparm  = nfp_ffn[efitfn];
  fprintf(out,"Will fit to the following function:\n");
  fprintf(out,"%s\n",longs_ffn[efitfn]);
  c1 = val[n];
  if (bYdy) {
    c1  = val[n];
    sig = val[n+1];
    fprintf(out,"Using two columns as y and sigma values\n");
  } else {
    snew(sig,ny);
  }
  if (opt2parg_bSet("-beginfit",npargs,ppa)) {
    tbeginfit = opt2parg_real("-beginfit",npargs,ppa);
  } else {
    tbeginfit = x0[0];
  }
  if (opt2parg_bSet("-endfit",npargs,ppa)) {
    tendfit   = opt2parg_real("-endfit",npargs,ppa);
  } else {
    tendfit   = x0[ny-1];
  }
  
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
  if (do_lmfit(ny,c1,sig,0,x0,tbeginfit,tendfit,
	       oenv,bDebugMode(),efitfn,fitparm,0)) {
    for(i=0; (i<nparm); i++) 
      fprintf(out,"a%-2d = %12.5e\n",i+1,fitparm[i]);
  }
  else {
    fprintf(out,"No solution was found\n");
  }
}

static void do_ballistic(const char *balFile, int nData,
                         real *t, real **val, int nSet,
                         real balTime, int nBalExp,
                         gmx_bool bDerivative,
                         const output_env_t oenv)
{
  double **ctd=NULL, *td=NULL;
  t_gemParams *GP = init_gemParams(0, 0, t, nData, 0, 0, 0, balTime, nBalExp, bDerivative);
  static char *leg[] = {"Ac'(t)"};
  FILE *fp;
  int i, set;
  
  if (GP->ballistic/GP->tDelta >= GP->nExpFit*2+1)
  {
      snew(ctd, nSet);
      snew(td,  nData);

      fp = xvgropen(balFile, "Hydrogen Bond Autocorrelation","Time (ps)","C'(t)", oenv);
      xvgr_legend(fp,asize(leg),(const char**)leg,oenv);
      
      for (set=0; set<nSet; set++)
      {
          snew(ctd[set], nData);
          for (i=0; i<nData; i++) {
              ctd[set][i] = (double)val[set][i];
              if (set==0)
                  td[i] = (double)t[i];
          }
          
          takeAwayBallistic(ctd[set], td, nData, GP->ballistic, GP->nExpFit, GP->bDt);
      }
      
      for (i=0; i<nData; i++)
      {
          fprintf(fp, "  %g",t[i]);
          for (set=0; set<nSet; set++)
          {
              fprintf(fp, "  %g", ctd[set][i]);
          }
          fprintf(fp, "\n");
      }


      for (set=0; set<nSet; set++)
          sfree(ctd[set]);
      sfree(ctd);
      sfree(td);
  }
  else
      printf("Number of data points is less than the number of parameters to fit\n."
             "The system is underdetermined, hence no ballistic term can be found.\n\n");
}

static void do_geminate(const char *gemFile, int nData,
                        real *t, real **val, int nSet,
                        const real D, const real rcut, const real balTime,
                        const int nFitPoints, const real begFit, const real endFit,
                        const output_env_t oenv)
{
    double **ctd=NULL, **ctdGem=NULL, *td=NULL;
    t_gemParams *GP = init_gemParams(rcut, D, t, nData, nFitPoints,
                                     begFit, endFit, balTime, 1, FALSE);
    const char *leg[] = {"Ac\\sgem\\N(t)"};
    FILE *fp;
    int i, set;
    
    snew(ctd,    nSet);
    snew(ctdGem, nSet);
    snew(td,  nData);
    
    fp = xvgropen(gemFile, "Hydrogen Bond Autocorrelation","Time (ps)","C'(t)", oenv);
    xvgr_legend(fp,asize(leg),leg,oenv);
    
    for (set=0; set<nSet; set++)
    {
        snew(ctd[set],    nData);
        snew(ctdGem[set], nData);
        for (i=0; i<nData; i++) {
            ctd[set][i] = (double)val[set][i];
            if (set==0)
                td[i] = (double)t[i];
        }
        fitGemRecomb(ctd[set], td, &(ctd[set]), nData, GP);
    }

    for (i=0; i<nData; i++)
	{
        fprintf(fp, "  %g",t[i]);
        for (set=0; set<nSet; set++)
	    {
            fprintf(fp, "  %g", ctdGem[set][i]);
	    }
        fprintf(fp, "\n");
	}

    for (set=0; set<nSet; set++)
    {
        sfree(ctd[set]);
        sfree(ctdGem[set]);
    }
    sfree(ctd);
    sfree(ctdGem);
    sfree(td);
}

int gmx_analyze(int argc,char *argv[])
{
  static const char *desc[] = {
    "g_analyze reads an ascii file and analyzes data sets.",
    "A line in the input file may start with a time",
    "(see option [TT]-time[tt]) and any number of y values may follow.",
    "Multiple sets can also be",
    "read when they are separated by & (option [TT]-n[tt]),",
    "in this case only one y value is read from each line.",
    "All lines starting with # and @ are skipped.",
    "All analyses can also be done for the derivative of a set",
    "(option [TT]-d[tt]).[PAR]",

    "All options, except for [TT]-av[tt] and [TT]-power[tt] assume that the",
    "points are equidistant in time.[PAR]",

    "g_analyze always shows the average and standard deviation of each",
    "set. For each set it also shows the relative deviation of the third",
    "and fourth cumulant from those of a Gaussian distribution with the same",
    "standard deviation.[PAR]",

    "Option [TT]-ac[tt] produces the autocorrelation function(s).[PAR]",
    
    "Option [TT]-cc[tt] plots the resemblance of set i with a cosine of",
    "i/2 periods. The formula is:[BR]"
    "2 (int0-T y(t) cos(i pi t) dt)^2 / int0-T y(t) y(t) dt[BR]",
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
    "The analytical curve for the block average is:[BR]",
    "f(t) = sigma sqrt(2/T (  a   (tau1 ((exp(-t/tau1) - 1) tau1/t + 1)) +[BR]",
    "                       (1-a) (tau2 ((exp(-t/tau2) - 1) tau2/t + 1)))),[BR]"
    "where T is the total time.",
    "a, tau1 and tau2 are obtained by fitting f^2(t) to error^2.",
    "When the actual block average is very close to the analytical curve,",
    "the error is sigma*sqrt(2/T (a tau1 + (1-a) tau2)).",
    "The complete derivation is given in",
    "B. Hess, J. Chem. Phys. 116:209-217, 2002.[PAR]",

    "Option [TT]-bal[tt] finds and subtracts the ultrafast \"ballistic\"",
    "component from a hydrogen bond autocorrelation function by the fitting",
    "of a sum of exponentials, as described in e.g.",
    "O. Markovitch, J. Chem. Phys. 129:084505, 2008. The fastest term",
    "is the one with the most negative coefficient in the exponential,",
    "or with [TT]-d[tt], the one with most negative time derivative at time 0.",
    "[TT]-nbalexp[tt] sets the number of exponentials to fit.[PAR]",

    "Option [TT]-gem[tt] fits bimolecular rate constants ka and kb",
    "(and optionally kD) to the hydrogen bond autocorrelation function",
    "according to the reversible geminate recombination model. Removal of",
    "the ballistic component first is strongly adviced. The model is presented in",
    "O. Markovitch, J. Chem. Phys. 129:084505, 2008.[PAR]",

    "Option [TT]-filter[tt] prints the RMS high-frequency fluctuation",
    "of each set and over all sets with respect to a filtered average.",
    "The filter is proportional to cos(pi t/len) where t goes from -len/2",
    "to len/2. len is supplied with the option [TT]-filter[tt].",
    "This filter reduces oscillations with period len/2 and len by a factor",
    "of 0.79 and 0.33 respectively.[PAR]",

    "Option [TT]-g[tt] fits the data to the function given with option",
    "[TT]-fitfn[tt].[PAR]",
    
    "Option [TT]-power[tt] fits the data to b t^a, which is accomplished",
    "by fitting to a t + b on log-log scale. All points after the first",
    "zero or negative value are ignored.[PAR]"
    
    "Option [TT]-luzar[tt] performs a Luzar & Chandler kinetics analysis",
    "on output from [TT]g_hbond[tt]. The input file can be taken directly",
    "from [TT]g_hbond -ac[tt], and then the same result should be produced."
  };
  static real tb=-1,te=-1,frac=0.5,filtlen=0,binwidth=0.1,aver_start=0;
  static gmx_bool bHaveT=TRUE,bDer=FALSE,bSubAv=TRUE,bAverCorr=FALSE,bXYdy=FALSE;
  static gmx_bool bEESEF=FALSE,bEENLC=FALSE,bEeFitAc=FALSE,bPower=FALSE;
  static gmx_bool bIntegrate=FALSE,bRegression=FALSE,bLuzar=FALSE,bLuzarError=FALSE; 
  static int  nsets_in=1,d=1,nb_min=4,resol=10, nBalExp=4, nFitPoints=100;
  static real temp=298.15,fit_start=1, fit_end=60, smooth_tail_start=-1, balTime=0.2, diffusion=5e-5,rcut=0.35;
  
  /* must correspond to enum avbar* declared at beginning of file */
  static const char *avbar_opt[avbarNR+1] = { 
    NULL, "none", "stddev", "error", "90", NULL
  };

  t_pargs pa[] = {
    { "-time",    FALSE, etBOOL, {&bHaveT},
      "Expect a time in the input" },
    { "-b",       FALSE, etREAL, {&tb},
      "First time to read from set" },
    { "-e",       FALSE, etREAL, {&te},
      "Last time to read from set" },
    { "-n",       FALSE, etINT, {&nsets_in},
      "Read # sets separated by &" },
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
    { "-regression",FALSE,etBOOL,{&bRegression},
      "Perform a linear regression analysis on the data. If -xydy is set a second set will be interpreted as the error bar in the Y value. Otherwise, if multiple data sets are present a multilinear regression will be performed yielding the constant A that minimize chi^2 = (y - A0 x0 - A1 x1 - ... - AN xN)^2 where now Y is the first data set in the input file and xi the others. Do read the information at the option [TT]-time[tt]." },
    { "-luzar",   FALSE, etBOOL, {&bLuzar},
      "Do a Luzar and Chandler analysis on a correlation function and related as produced by g_hbond. When in addition the -xydy flag is given the second and fourth column will be interpreted as errors in c(t) and n(t)." },
    { "-temp",    FALSE, etREAL, {&temp},
      "Temperature for the Luzar hydrogen bonding kinetics analysis" },
    { "-fitstart", FALSE, etREAL, {&fit_start},
      "Time (ps) from which to start fitting the correlation functions in order to obtain the forward and backward rate constants for HB breaking and formation" }, 
    { "-fitend", FALSE, etREAL, {&fit_end},
      "Time (ps) where to stop fitting the correlation functions in order to obtain the forward and backward rate constants for HB breaking and formation. Only with -gem" }, 
    { "-smooth",FALSE, etREAL, {&smooth_tail_start},
      "If >= 0, the tail of the ACF will be smoothed by fitting it to an exponential function: y = A exp(-x/tau)" },
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
      "Calculate one ACF over all sets" },
    { "-nbalexp", FALSE, etINT, {&nBalExp},
      "HIDDENNumber of exponentials to fit to the ultrafast component" },
    { "-baltime", FALSE, etREAL, {&balTime},
      "HIDDENTime up to which the ballistic component will be fitted" },
/*     { "-gemnp", FALSE, etINT, {&nFitPoints}, */
/*       "HIDDENNumber of data points taken from the ACF to use for fitting to rev. gem. recomb. model."}, */
/*     { "-rcut", FALSE, etREAL, {&rcut}, */
/*       "Cut-off for hydrogen bonds in geminate algorithms" }, */
/*     { "-gemtype", FALSE, etENUM, {gemType}, */
/*       "What type of gminate recombination to use"}, */
/*     { "-D", FALSE, etREAL, {&diffusion}, */
/*       "The self diffusion coefficient which is used for the reversible geminate recombination model."} */
  };
#define NPA asize(pa)

  FILE     *out,*out_fit;
  int      n,nlast,s,nset,i,j=0;
  real     **val,*t,dt,tot,error;
  double   *av,*sig,cum1,cum2,cum3,cum4,db;
  const char     *acfile,*msdfile,*ccfile,*distfile,*avfile,*eefile,*balfile,*gemfile,*fitfile;
  output_env_t oenv;
  
  t_filenm fnm[] = { 
    { efXVG, "-f",    "graph",    ffREAD   },
    { efXVG, "-ac",   "autocorr", ffOPTWR  },
    { efXVG, "-msd",  "msd",      ffOPTWR  },
    { efXVG, "-cc",   "coscont",  ffOPTWR  },
    { efXVG, "-dist", "distr",    ffOPTWR  },
    { efXVG, "-av",   "average",  ffOPTWR  },
    { efXVG, "-ee",   "errest",   ffOPTWR  },
    { efXVG, "-bal",  "ballisitc",ffOPTWR  },
/*     { efXVG, "-gem",  "geminate", ffOPTWR  }, */
    { efLOG, "-g",    "fitlog",   ffOPTWR  }
  }; 
#define NFILE asize(fnm) 

  int     npargs;
  t_pargs *ppa;

  npargs = asize(pa); 
  ppa    = add_acf_pargs(&npargs,pa);
  
  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_VIEW,
		    NFILE,fnm,npargs,ppa,asize(desc),desc,0,NULL,&oenv); 

  acfile   = opt2fn_null("-ac",NFILE,fnm);
  msdfile  = opt2fn_null("-msd",NFILE,fnm);
  ccfile   = opt2fn_null("-cc",NFILE,fnm);
  distfile = opt2fn_null("-dist",NFILE,fnm);
  avfile   = opt2fn_null("-av",NFILE,fnm);
  eefile   = opt2fn_null("-ee",NFILE,fnm);
  balfile  = opt2fn_null("-bal",NFILE,fnm);
/*   gemfile  = opt2fn_null("-gem",NFILE,fnm); */
    /* When doing autocorrelation we don't want a fitlog for fitting
     * the function itself (not the acf) when the user did not ask for it.
     */
    if (opt2parg_bSet("-fitfn",npargs,ppa) && acfile == NULL)
    {
        fitfile  = opt2fn("-g",NFILE,fnm);
    }
    else
    {
        fitfile  = opt2fn_null("-g",NFILE,fnm);
    }
    
    val = read_xvg_time(opt2fn("-f",NFILE,fnm),bHaveT,
                        opt2parg_bSet("-b",npargs,ppa),tb,
                        opt2parg_bSet("-e",npargs,ppa),te,
                        nsets_in,&nset,&n,&dt,&t);
    printf("Read %d sets of %d points, dt = %g\n\n",nset,n,dt);
  
    if (bDer)
    {
        printf("Calculating the derivative as (f[i+%d]-f[i])/(%d*dt)\n\n",
               d,d);
        n -= d;
        for(s=0; s<nset; s++)
        {
            for(i=0; (i<n); i++)
            {
                val[s][i] = (val[s][i+d]-val[s][i])/(d*dt);
            }
        }
    }
    
    if (bIntegrate)
    {
        real sum,stddev;

        printf("Calculating the integral using the trapezium rule\n");
    
        if (bXYdy)
        {
            sum = evaluate_integral(n,t,val[0],val[1],aver_start,&stddev);
            printf("Integral %10.3f +/- %10.5f\n",sum,stddev);
        }
        else
        {
            for(s=0; s<nset; s++)
            {
                sum = evaluate_integral(n,t,val[s],NULL,aver_start,&stddev);
                printf("Integral %d  %10.5f  +/- %10.5f\n",s+1,sum,stddev);
            }
        }
    }

    if (fitfile != NULL)
    {
        out_fit = ffopen(fitfile,"w");
        if (bXYdy && nset >= 2)
        {
            do_fit(out_fit,0,TRUE,n,t,val,npargs,ppa,oenv);
        }
        else
        {
            for(s=0; s<nset; s++)
            {
                do_fit(out_fit,s,FALSE,n,t,val,npargs,ppa,oenv);
            }
        }
        ffclose(out_fit);
    }

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
    filter(filtlen,n,nset,val,dt,oenv);
  
  if (msdfile) {
    out=xvgropen(msdfile,"Mean square displacement",
		 "time","MSD (nm\\S2\\N)",oenv);
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
    ffclose(out);
    fprintf(stderr,"\r%d, time=%g\n",j-1,(j-1)*dt);
  }
  if (ccfile)
    plot_coscont(ccfile,n,nset,val,oenv);
  
  if (distfile)
    histogram(distfile,binwidth,n,nset,val,oenv);
  if (avfile)
    average(avfile,nenum(avbar_opt),n,nset,val,t);
  if (eefile)
    estimate_error(eefile,nb_min,resol,n,nset,av,sig,val,dt,
		   bEeFitAc,bEESEF,bEENLC,oenv);
  if (balfile)
      do_ballistic(balfile,n,t,val,nset,balTime,nBalExp,bDer,oenv);
/*   if (gemfile) */
/*       do_geminate(gemfile,n,t,val,nset,diffusion,rcut,balTime, */
/*                   nFitPoints, fit_start, fit_end, oenv); */
  if (bPower)
    power_fit(n,nset,val,t);

    if (acfile != NULL)
    {
        if (bSubAv)
        {
            for(s=0; s<nset; s++)
            {
                for(i=0; i<n; i++)
                {
                    val[s][i] -= av[s];
                }
            }
        }
        do_autocorr(acfile,oenv,"Autocorrelation",n,nset,val,dt,
                    eacNormal,bAverCorr);
    }

  if (bRegression)
      regression_analysis(n,bXYdy,t,nset,val);

  if (bLuzar) 
    luzar_correl(n,t,nset,val,temp,bXYdy,fit_start,smooth_tail_start,oenv);
    
  view_all(oenv,NFILE, fnm);
  
  thanx(stderr);

  return 0;
}
  
