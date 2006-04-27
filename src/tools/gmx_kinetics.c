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
#include <stdlib.h>
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
#include "xvgr.h"
#include "physics.h"
#include <gsl/gsl_multimin.h>

enum { epAuf, epEuf, epAfu, epEfu, epNR };
enum { eqAif, eqEif, eqAfi, eqEfi, eqAui, eqEui, eqAiu, eqEiu, eqNR };
static char *eep[epNR] = { "Af", "Ef", "Au", "Eu" };
static char *eeq[eqNR] = { "Aif","Eif","Afi","Efi","Aui","Eui","Aiu","Eiu" };

typedef struct {
  int nreplica;   /* Number of replicas in the calculation                   */
  int nframe;     /* Number of time frames                                   */
  int nstate;     /* Number of states the system can be in, e.g. F,I,U       */
  int nparams;    /* Is 2, 4 or 8                                            */
  bool *bMask;    /* Determine whether this replica is part of the d2 comp.  */
  bool bSum;
  int nmask;      /* Number of replicas taken into account                   */
  real dt;        /* Timestep between frames                                 */
  int  j0,j1;     /* Range of frames used in calculating delta               */
  real **temp,**data;
  int  **state;   /* State index running from 0 (F) to nstate-1 (U)          */
  real **beta,**fcalt,**icalt;
  real *time,*sumft,*sumit,*sumfct,*sumict;
  real *params;
  real *d2_replica;
} t_remd_data;

char *itoa(int i)
{
  static char ptr[12];
  
  sprintf(ptr,"%d",i);
  return ptr;
}

char *epnm(int nparams,int index)
{
  static char buf[32],from[8],to[8];
  int    nn,ni,ii;
  
  range_check(index,0,nparams);
  if ((nparams == 2) || (nparams == 4))
    return eep[index];
  else if ((nparams > 4) && (nparams % 4 == 0))
    return eeq[index];
  else 
    gmx_fatal(FARGS,"Don't know how to handle %d parameters",nparams);
    
  return NULL;
}

bool bBack(t_remd_data *d) 
{
  return (d->nparams > 2);
}

real is_folded(t_remd_data *d,int irep,int jframe)
{
  if (d->state[irep][jframe] == 0)
    return 1.0;
  else
    return 0.0;
}

real is_unfolded(t_remd_data *d,int irep,int jframe)
{
  if (d->state[irep][jframe] == d->nstate-1)
    return 1.0;
  else
    return 0.0;
}

real is_intermediate(t_remd_data *d,int irep,int jframe)
{
  if ((d->state[irep][jframe] == 1) && (d->nstate > 2))
    return 1.0;
  else
    return 0.0;
}

void integrate_dfdt(t_remd_data *d)
{
  int    i,j;
  double beta,ddf,ddi,df,db,fac,sumf,sumi;
  
  d->sumfct[0] = 0;
  d->sumict[0] = 0;
  for(i=0; (i<d->nreplica); i++) {
    if (d->bMask[i]) {
      ddf = 0.5*d->dt*is_folded(d,i,0);
      ddi = 0.5*d->dt*is_intermediate(d,i,0);
      d->fcalt[i][0] = ddf;
      d->icalt[i][0] = ddi;
      d->sumfct[0]  += ddf;
      d->sumict[0]  += ddi;
    }
  }
  for(j=1; (j<d->nframe); j++) {
    if (j==d->nframe-1)
      fac = 0.5*d->dt;
    else
      fac = d->dt;
    sumf = sumi = 0;
    for(i=0; (i<d->nreplica); i++) {
      if (d->bMask[i]) {
	beta = d->beta[i][j];
	if (d->nstate <= 2) {
	  df = (d->params[epAuf]*exp(-beta*d->params[epEuf])*
		is_unfolded(d,i,j));
	  if (bBack(d)) {
	    db = (d->params[epAfu]*exp(-beta*d->params[epEfu])*
		  is_folded(d,i,j));
	    ddf = fac*(df-db);
	  }
	  else
	    ddf = fac*df;
	  d->fcalt[i][j] = d->fcalt[i][j-1] + ddf;
	  sumf += ddf;
	}
	else {
	  ddf = fac*((d->params[eqAif]*exp(-beta*d->params[eqEif])*
		      is_intermediate(d,i,j)) -
		     (d->params[eqAfi]*exp(-beta*d->params[eqEfi])*
		      is_folded(d,i,j)));
	  ddi = fac*((d->params[eqAui]*exp(-beta*d->params[eqEui])*
		      is_unfolded(d,i,j)) -
		     (d->params[eqAiu]*exp(-beta*d->params[eqEiu])*
		      is_intermediate(d,i,j)));
	  d->fcalt[i][j] = d->fcalt[i][j-1] + ddf;
	  d->icalt[i][j] = d->icalt[i][j-1] + ddi;
	  sumf += ddf;
	  sumi += ddi;
	}
      }
    }
    d->sumfct[j] = d->sumfct[j-1] + sumf;
    d->sumict[j] = d->sumict[j-1] + sumi;
  }
}

void sum_ft(t_remd_data *d)
{
  int i,j;
  double fac;
    
  for(j=0; (j<d->nframe); j++) {
    d->sumft[j] = 0;
    d->sumit[j] = 0;
    if ((j == 0) || (j==d->nframe-1))
      fac = d->dt*0.5;
    else
      fac = d->dt;
    for(i=0; (i<d->nreplica); i++) {
      if (d->bMask[i]) {
	d->sumft[j] += fac*is_folded(d,i,j);
	d->sumit[j] += fac*is_intermediate(d,i,j);
      }
    }
  }
}

double calc_d2(t_remd_data *d)
{
  int    i,j;
  double dd2,d2=0,dr2,tmp;
  
  integrate_dfdt(d);
  
  if (d->bSum) {
    for(j=d->j0; (j<d->j1); j++) {
      d2  += sqr(d->sumft[j]-d->sumfct[j]);
      if (d->nstate > 2)
	d2 += sqr(d->sumit[j]-d->sumict[j]);
    }
  }
  else {
    for(i=0; (i<d->nreplica); i++) {
      dr2 = 0;
      if (d->bMask[i]) {
	for(j=d->j0; (j<d->j1); j++) {
	  tmp  = sqr(is_folded(d,i,j)-d->fcalt[i][j]);
	  d2  += tmp;
	  dr2 += tmp; 
	  if (d->nstate > 2) {
	    tmp = sqr(is_intermediate(d,i,j)-d->icalt[i][j]);
	    d2  += tmp;
	    dr2 += tmp;
	  }
	}
	d->d2_replica[i] = dr2/(d->j1-d->j0);
      }
    }
  }
  dd2 = (d2/(d->j1-d->j0))/d->nmask;
   
  return dd2;
}

double my_f(const gsl_vector *v,void *params)
{
  t_remd_data *d = (t_remd_data *) params;
  double penalty=0;
  int i;
  
  for(i=0; (i<d->nparams); i++) {
    d->params[i] = gsl_vector_get(v, i);
    if (d->params[i] < 0)
      penalty += 10;
  }
  if (penalty > 0)
    return penalty;
  else
    return calc_d2(d);
}

static void optimize_remd_parameters(FILE *fp,t_remd_data *d,int maxiter,
				     real tol)
{
  real   size,d2;
  size_t iter = 0;
  int    status = 0;
  int    i;

  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  gsl_vector *x,*dx;
  gsl_multimin_function my_func;

  my_func.f      = &my_f;
  my_func.n      = d->nparams;
  my_func.params = (void *) d;

  /* Starting point */
  x = gsl_vector_alloc (my_func.n);
  for(i=0; (i<my_func.n); i++)
    gsl_vector_set (x, i, d->params[i]);
  
  /* Step size, different for each of the parameters */
  dx = gsl_vector_alloc (my_func.n);
  for(i=0; (i<my_func.n); i++)
    gsl_vector_set (dx, i, 0.1*d->params[i]);

  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc (T, my_func.n);

  gsl_multimin_fminimizer_set (s, &my_func, x, dx);
  gsl_vector_free (x);
  gsl_vector_free (dx);

  printf ("%5s","Iter");
  for(i=0; (i<my_func.n); i++) 
    printf(" %8s",epnm(my_func.n,i));
  printf (" %8s %8s\n","NM Size","Chi2");
  
  do  {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);
    
    if (status != 0)
      gmx_fatal(FARGS,"Something went wrong in the iteration in minimizer %s",
		gsl_multimin_fminimizer_name(s));
		
    d2     = gsl_multimin_fminimizer_minimum(s);
    size   = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size,tol);
    
    if (status == GSL_SUCCESS)
      printf ("Minimum found using %s at:\n",
	      gsl_multimin_fminimizer_name(s));
  
    printf ("%5d", iter);
    for(i=0; (i<my_func.n); i++) 
      printf(" %8.2e",gsl_vector_get (s->x,i));
    printf (" %8.2e %8.2e\n",size,d2);
  }
  while ((status == GSL_CONTINUE) && (iter < maxiter));
  
  gsl_multimin_fminimizer_free (s);
}


static void preprocess_remd(FILE *fp,t_remd_data *d,real cutoff,real tref,
			    real ucut,bool bBack,real Euf,real Efu,
			    real Ei,real t0,real t1,bool bSum)
{
  int i,j,ninter;
  real dd,tau_f,tau_u;

  ninter = (ucut > cutoff) ? 1 : 0;
  if (ninter && (ucut <= cutoff))
    gmx_fatal(FARGS,"You have requested an intermediate but the cutoff for intermediates %f is smaller than the normal cutoff(%f)",ucut,cutoff);
    
  if (!bBack) {
    d->nparams = 2;
    d->nstate  = 2;
  }
  else {
    d->nparams = 4*(1+ninter);
    d->nstate  = 2+ninter;
  }
  d->bSum = bSum;
  snew(d->beta, d->nreplica);
  snew(d->state,d->nreplica);
  snew(d->bMask,d->nreplica);
  snew(d->d2_replica,d->nreplica);
  snew(d->sumft,d->nframe);
  snew(d->sumit,d->nframe);
  snew(d->sumfct,d->nframe);
  snew(d->sumict,d->nframe);
  snew(d->params,d->nparams);
  snew(d->fcalt,d->nreplica);
  snew(d->icalt,d->nreplica);

  if (t0 < 0)
    d->j0 = 0;
  else
    for(d->j0=0; (d->j0<d->nframe) && (d->time[d->j0] < t0); d->j0++)
      ;
  if (t1 < 0)
    d->j1=d->nframe;
  else
    for(d->j1=0; (d->j1<d->nframe) && (d->time[d->j1] < t1); d->j1++)
      ;
  if ((d->j1-d->j0) < d->nparams+2)
    gmx_fatal(FARGS,"Start (%f) or end time (%f) for fitting inconsistent. Reduce t0, increase t1 or supply more data",t0,t1);
  fprintf(fp,"Will optimize from %g to %g\n",
	  d->time[d->j0],d->time[d->j1-1]);  
  d->nmask = d->nreplica;
  for(i=0; (i<d->nreplica); i++) {
    snew(d->beta[i],d->nframe);
    snew(d->state[i],d->nframe);
    snew(d->fcalt[i],d->nframe);
    snew(d->icalt[i],d->nframe);
    d->bMask[i] = TRUE;
    for(j=0; (j<d->nframe); j++) {
      d->beta[i][j] = 1.0/(BOLTZ*d->temp[i][j]);
      dd = d->data[i][j];
      if (dd <= cutoff)
	d->state[i][j] = 0;
      else if ((ucut > cutoff) && (dd <= ucut))
	d->state[i][j] = 1;
      else
	d->state[i][j] = d->nstate-1;
    }
  }
  sum_ft(d);
  
  /* Assume forward rate constant is half the total time in this
   * simulation and backward is ten times as long */
  tau_f = d->time[d->nframe-1];
  tau_u = 4*tau_f;
  d->params[epEuf] = Euf;
  d->params[epAuf] = exp(d->params[epEuf]/(BOLTZ*tref))/tau_f;
  if (bBack) {
    d->params[epEfu] = Efu;
    d->params[epAfu] = exp(d->params[epEfu]/(BOLTZ*tref))/tau_u;
    if (ninter >0) {
      d->params[eqEui] = Ei;
      d->params[eqAui] = exp(d->params[eqEui]/(BOLTZ*tref))/tau_u;
      d->params[eqEiu] = Ei;
      d->params[eqAiu] = exp(d->params[eqEiu]/(BOLTZ*tref))/tau_u;
    }
  }
  else {
    d->params[epAfu]  = 0;
    d->params[epEfu]  = 0;
  }
}

real tau(real A,real E,real T)
{
  return exp(E/(BOLTZ*T))/A;
}

real folded_fraction(t_remd_data *d,real tref)
{
  real tauf,taub;
  
  tauf = tau(d->params[epAuf],d->params[epEuf],tref);
  taub = tau(d->params[epAfu],d->params[epEfu],tref);
  
  return (taub/(tauf+taub));
}

void print_tau(FILE *gp,t_remd_data *d,real tref)
{
  real tauf,taub,ddd,fff,DG,DH,TDS,Tm,Tb,Te,Fb,Fe,Fm;
  int  i,np=d->nparams;
  
  ddd = calc_d2(d);
  fprintf(gp,"Final value for Chi2 = %12.5e (%d replicas)\n",ddd,d->nmask);
  tauf = tau(d->params[epAuf],d->params[epEuf],tref);
  fprintf(gp,"%s = %12.5e %s = %12.5e (kJ/mole)\n",
	  epnm(np,epAuf),d->params[epAuf],
	  epnm(np,epEuf),d->params[epEuf]);
  if (bBack(d)) {
    taub = tau(d->params[epAfu],d->params[epEfu],tref);
    fprintf(gp,"%s = %12.5e %s = %12.5e (kJ/mole)\n",
	    epnm(np,epAfu),d->params[epAfu],
	    epnm(np,epEfu),d->params[epEfu]);
    fprintf(gp,"Equilibrium properties at T = %g\n",tref);
    fprintf(gp,"tau_f = %8.3f ns, tau_b = %8.3f ns\n",tauf/1000,taub/1000);
    fff = taub/(tauf+taub);
    DG  = BOLTZ*tref*log(fff/(1-fff));
    DH  = d->params[epEfu]-d->params[epEuf];
    TDS = DH-DG;
    fprintf(gp,"Folded fraction     F = %8.3f\n",fff);
    fprintf(gp,"Unfolding energies: DG = %8.3f  DH = %8.3f TDS = %8.3f\n",
	    DG,DH,TDS);
    Tb = 260;
    Te = 420;
    Fb = folded_fraction(d,Tb);
    Fe = folded_fraction(d,Te);
    while((Te-Tb > 0.001) && (Fm != 0.5)) {
      Tm = 0.5*(Tb+Te);
      Fm = folded_fraction(d,Tm);
      if (Fm > 0.5) {
	Fb = Fm;
	Tb = Tm;
      }
      else if (Fm < 0.5) {
	Te = Tm;
	Fe = Fm;
      }
    }
    if ((Fb-0.5)*(Fe-0.5) <= 0)
      fprintf(gp,"Melting temperature Tm = %8.3f K\n",Tm);
    else 
      fprintf(gp,"No melting temperature detected between 260 and 420K\n");
    if (np > 4) {
      char *ptr;
      fprintf(gp,"Data for intermediates at T = %g\n",tref);
      fprintf(gp,"%8s  %10s  %10s  %10s\n","Name","A","E","tau");
      for(i=0; (i<np/2); i++) {
	tauf = tau(d->params[2*i],d->params[2*i+1],tref);
	ptr = epnm(d->nparams,2*i);
	fprintf(gp,"%8s  %10.3e  %10.3e  %10.3e\n",ptr+1,
		d->params[2*i],d->params[2*i+1],tauf/1000);
      }
    }
  }
  else {
    fprintf(gp,"Equilibrium properties at T = %g\n",tref);
    fprintf(gp,"tau_f = %8.3f\n",tauf);
  }
}

void dump_remd_parameters(FILE *gp,t_remd_data *d,char *fn,char *fn2,char *rfn,
			  char *efn,char *mfn,int skip,real tref)
{
  FILE *fp,*hp;
  int  i,j,np=d->nparams;
  real rhs,tauf,taub,fff,DG;
  real *params;
  char *leg[] = { "Measured", "Fit", "Difference" };
  char *mleg[] = { "Folded fraction","DG (kJ/mole)"};
  char **rleg;
  real fac[] = { 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03 };
#define NFAC asize(fac)
  real d2[NFAC];
  
  integrate_dfdt(d);
  print_tau(gp,d,tref);
  
  if (fn) {
    fp = xvgropen(fn,"Optimized fit to data","Time (ps)","Fraction Folded");
    xvgr_legend(fp,asize(leg),leg);
    for(i=0; (i<d->nframe); i++) {
      if ((skip <= 0) || ((i % skip) == 0))
	fprintf(fp,"%12.5e  %12.5e  %12.5e\n",d->time[i],
		d->sumft[i]/d->nmask,d->sumfct[i]/d->nmask,
		(d->sumft[i]-d->sumfct[i])/d->nmask);
    }
    fclose(fp);
  }
  if (!d->bSum && rfn) {
    snew(rleg,d->nreplica*2);
    for(i=0; (i<d->nreplica); i++) {
      snew(rleg[2*i],32);
      snew(rleg[2*i+1],32);
      sprintf(rleg[2*i],"\\f{4}F(t) %d",i);
      sprintf(rleg[2*i+1],"\\f{12}F \\f{4}(t) %d",i);
    }
    fp = xvgropen(rfn,"Optimized fit to data","Time (ps)","Fraction Folded");
    xvgr_legend(fp,d->nreplica*2,rleg);
    for(j=0; (j<d->nframe); j++) {
      if ((skip <= 0) || ((j % skip) == 0)) {
	fprintf(fp,"%12.5e",d->time[j]);
	for(i=0; (i<d->nreplica); i++) 
	  fprintf(fp,"  %5f  %9.2e",is_folded(d,i,j),d->fcalt[i][j]);
	fprintf(fp,"\n");
      }
    }
    fclose(fp);
  }

  if (fn2 && (d->nstate > 2)) {
    fp = xvgropen(fn2,"Optimized fit to data","Time (ps)",
		  "Fraction Intermediate");
    xvgr_legend(fp,asize(leg),leg);
    for(i=0; (i<d->nframe); i++) {
      if ((skip <= 0) || ((i % skip) == 0))
	fprintf(fp,"%12.5e  %12.5e  %12.5e  %12.5e\n",d->time[i],
		d->sumit[i]/d->nmask,d->sumict[i]/d->nmask,
		(d->sumit[i]-d->sumict[i])/d->nmask);
    }
    fclose(fp);
  }
  if (mfn) {
    if (bBack(d)) {
      fp = xvgropen(mfn,"Melting curve","T (K)","");
      xvgr_legend(fp,asize(mleg),mleg);
      for(i=260; (i<=420); i++) {
	tauf = tau(d->params[epAuf],d->params[epEuf],1.0*i);
	taub = tau(d->params[epAfu],d->params[epEfu],1.0*i);
	fff  = taub/(tauf+taub);
	DG   = BOLTZ*i*log(fff/(1-fff));
	fprintf(fp,"%5d  %8.3f  %8.3f\n",i,fff,DG);
      }
      fclose(fp);
    }
  }
  
  if (efn) {
    snew(params,d->nparams);
    for(i=0; (i<d->nparams); i++)
      params[i] = d->params[i];
    
    hp = xvgropen(efn,"Chi2 as a function of relative parameter",
		  "Fraction","Chi2");
    for(j=0; (j<d->nparams); j++) {
      /* Reset all parameters to optimized values */
      fprintf(hp,"@type xy\n");
      for(i=0; (i<d->nparams); i++) 
	d->params[i] = params[i];
      /* Now modify one of them */
      for(i=0; (i<NFAC); i++) {
	d->params[j] = fac[i]*params[j];
	d2[i] = calc_d2(d);
	fprintf(gp,"%s = %12g  d2 = %12g\n",epnm(np,j),d->params[j],d2[i]);
	fprintf(hp,"%12g  %12g\n",fac[i],d2[i]);
      }
      fprintf(hp,"&\n");
    }
    fclose(hp);
    for(i=0; (i<d->nparams); i++) 
      d->params[i] = params[i];
    sfree(params);
  }
  if (!d->bSum) {
    for(i=0; (i<d->nreplica); i++)
      fprintf(gp,"Chi2[%3d] = %8.2e\n",i,d->d2_replica[i]);
  }
}

int gmx_kinetics(int argc,char *argv[])
{
  static char *desc[] = {
    "g_kinetics reads two xvg files, each one containing data for N replicas.",
    "The first file contains the temperature of each replica at each timestep,"
    "and the second contains real values that can be interpreted as",
    "an indicator for folding. If the value in the file is larger than",
    "the cutoff it is taken to be unfolded and the other way around.[PAR]",
    "From these data an estimate of the forward and backward rate constants",
    "for folding is made at a reference temperature. In addition,"
    "a theoretical melting curve and free energy as a function of temperature"
    "are printed in an xvg file.[PAR]",
    "The user can give a max value to be regarded as intermediate",
    "([TT]-ucut[tt]), which, when given will trigger the use of an intermediate state",
    "in the algorithm to be defined as those structures that have",
    "cutoff < DATA < ucut. Structures with DATA values larger than ucut will",
    "not be regarded as potential folders. In this case 8 parameters are optimized.[PAR]",
    "The average fraction foled is printed in an xvg file together with the fit to it.", 
    "If an intermediate is used a further file will show the build of the intermediate and the fit to that process."
  };
  static int  nreplica = 1;
  static real tref     = 298.15;
  static real cutoff   = 0.2;
  static real ucut     = 0.0;
  static real Euf      = 10;
  static real Efu      = 30;
  static real Ei       = 10;
  static bool bHaveT   = TRUE;
  static real t0       = -1;
  static real t1       = -1;
  static real tb       = 0;
  static real te       = 0;
  static real tol      = 1e-3;
  static int  maxiter  = 100;
  static int  skip     = 0;
  static bool bBack    = TRUE;
  static bool bSplit   = TRUE;
  static bool bSum     = TRUE;
  t_pargs pa[] = {
    { "-time",    FALSE, etBOOL, {&bHaveT},
      "Expect a time in the input" },
    { "-b",       FALSE, etREAL, {&tb},
      "First time to read from set" },
    { "-e",       FALSE, etREAL, {&te},
      "Last time to read from set" },
    { "-bfit",    FALSE, etREAL, {&t0},
      "Time to start the fit from" },
    { "-efit",    FALSE, etREAL, {&t1},
      "Time to end the fit" },
    { "-T",       FALSE, etREAL, {&tref},
      "Reference temperature for computing rate constants" },
    { "-n",       FALSE, etINT, {&nreplica},
      "Read data for # replicas" },
    { "-cut",     FALSE, etREAL, {&cutoff},
      "Cut-off (max) value for regarding a structure as folded" },
    { "-ucut",    FALSE, etREAL, {&ucut},
      "Cut-off (max) value for regarding a structure as intermediate (if not folded)" },
    { "-euf",     FALSE, etREAL, {&Euf},
      "Initial guess for energy of activation for folding (kJ/mole)" },
    { "-efu",     FALSE, etREAL, {&Efu},
      "Initial guess for energy of activation for unfolding (kJ/mole)" },
    { "-ei",      FALSE, etREAL, {&Ei},
      "Initial guess for energy of activation for intermediates (kJ/mole)" },
    { "-maxiter", FALSE, etINT, {&maxiter},
      "Max number of iterations" },
    { "-back",    FALSE, etBOOL, {&bBack},
      "Take the back reaction into account" },
    { "-tol",     FALSE, etREAL,{&tol},
      "Absolute tolerance for convergence of the Nelder and Mead simplex algorithm" },
    { "-skip",    FALSE, etINT, {&skip},
      "Skip points in the output xvg file" },
    { "-split",   FALSE, etBOOL,{&bSplit},
      "Estimate error by splitting the number of replicas in two and refitting" },
    { "-sum",     FALSE, etBOOL,{&bSum},
      "Average folding before computing chi^2" }
  };
#define NPA asize(pa)

  FILE        *fp;
  real        dt_t,dt_d;
  int         nset_t,nset_d,n_t,n_d,i;
  char        *tfile,*dfile;
  t_remd_data remd;  
  
  t_filenm fnm[] = { 
    { efXVG, "-f",    "temp",    ffREAD   },
    { efXVG, "-d",    "data",    ffREAD   },
    { efXVG, "-o",    "ft_all",  ffWRITE  },
    { efXVG, "-o2",   "it_all",  ffOPTWR  },
    { efXVG, "-o3",   "ft_repl", ffOPTWR  },
    { efXVG, "-ee",   "err_est", ffOPTWR  },
    { efLOG, "-g",    "remd",    ffWRITE  },
    { efXVG, "-m",    "melt",    ffWRITE  }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL); 

  if (cutoff < 0)
    gmx_fatal(FARGS,"cutoff should be >= 0 (rather than %f)",cutoff);
		    
  tfile  = opt2fn("-f",NFILE,fnm);
  dfile  = opt2fn("-d",NFILE,fnm);
  
  fp = ffopen(opt2fn("-g",NFILE,fnm),"w");
  
  remd.temp = read_xvg_time(tfile,bHaveT,
			    opt2parg_bSet("-b",NPA,pa),tb,
			    opt2parg_bSet("-e",NPA,pa),te,
			    nreplica,&nset_t,&n_t,&dt_t,&remd.time);
  printf("Read %d sets of %d points in %s, dt = %g\n\n",nset_t,n_t,tfile,dt_t);
  sfree(remd.time);
  
  remd.data = read_xvg_time(dfile,bHaveT,
			    opt2parg_bSet("-b",NPA,pa),tb,
			    opt2parg_bSet("-e",NPA,pa),te,
			    nreplica,&nset_d,&n_d,&dt_d,&remd.time);
  printf("Read %d sets of %d points in %s, dt = %g\n\n",nset_d,n_d,dfile,dt_d);
  
  if ((nset_t != nset_d) || (n_t != n_d) || (dt_t != dt_d))
    gmx_fatal(FARGS,"Files %s and %s are inconsistent. Check log file",
	      tfile,dfile);
  remd.nreplica = nset_d;
  remd.nframe   = n_d;
  remd.dt       = dt_d;
  preprocess_remd(fp,&remd,cutoff,tref,ucut,bBack,Euf,Efu,Ei,t0,t1,bSum);
  
  optimize_remd_parameters(fp,&remd,maxiter,tol);
  
  dump_remd_parameters(fp,&remd,opt2fn("-o",NFILE,fnm),
		       opt2fn_null("-o2",NFILE,fnm),
		       opt2fn_null("-o3",NFILE,fnm),
		       opt2fn_null("-ee",NFILE,fnm),
		       opt2fn("-m",NFILE,fnm),skip,tref);

  if (bSplit) {
    printf("Splitting set of replicas in two halves\n");
    for(i=0; (i<remd.nreplica); i++) 
      remd.bMask[i] = FALSE;
    remd.nmask = 0;
    for(i=0; (i<remd.nreplica); i+=2) {
      remd.bMask[i] = TRUE;
      remd.nmask++;
    }
    sum_ft(&remd);
    optimize_remd_parameters(fp,&remd,maxiter,tol);
    dump_remd_parameters(fp,&remd,"test1.xvg",NULL,NULL,NULL,NULL,skip,tref);
    
    for(i=0; (i<remd.nreplica); i++) 
      remd.bMask[i] = !remd.bMask[i];
    remd.nmask = remd.nreplica - remd.nmask;
    
    sum_ft(&remd);
    optimize_remd_parameters(fp,&remd,maxiter,tol);
    dump_remd_parameters(fp,&remd,"test2.xvg",NULL,NULL,NULL,NULL,skip,tref);
    
    for(i=0; (i<remd.nreplica); i++) 
      remd.bMask[i] = FALSE;
    remd.nmask = 0;
    for(i=0; (i<remd.nreplica/2); i++) {
      remd.bMask[i] = TRUE;
      remd.nmask++;
    }
    sum_ft(&remd);
    optimize_remd_parameters(fp,&remd,maxiter,tol);
    dump_remd_parameters(fp,&remd,"test1.xvg",NULL,NULL,NULL,NULL,skip,tref);
    
    for(i=0; (i<remd.nreplica); i++) 
      remd.bMask[i] = FALSE;
    remd.nmask = 0;
    for(i=remd.nreplica/2; (i<remd.nreplica); i++) {
      remd.bMask[i] = TRUE;
      remd.nmask++;
    }
    sum_ft(&remd);
    optimize_remd_parameters(fp,&remd,maxiter,tol);
    dump_remd_parameters(fp,&remd,"test1.xvg",NULL,NULL,NULL,NULL,skip,tref);
  }
  fclose(fp);
  
  view_all(NFILE, fnm);
  
  thanx(stderr);
  
  return 0;
}
  
