/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */
static char *SRCID_xutils_c = "$Id$";

#include "typedefs.h"
#include "smalloc.h"
#include "strdb.h"
#include "string2.h"
#include "xmdrun.h"
#include "vec.h"
#include "genalg.h"
#include "random.h"

real mol_dipole(int k0,int k1,atom_id ma[],rvec x[],real q[])
{
  int  k,kk,m;
  rvec mu;
  
  clear_rvec(mu);
  for(k=k0; (k<k1); k++) {
    kk = ma[k];
    for(m=0; (m<DIM); m++) 
      mu[m] += q[kk]*x[kk][m];
  }
  return norm(mu);  /* Dipole moment of this molecule in e nm */
}

real calc_mu_aver(t_commrec *cr,t_nsborder *nsb,rvec x[],real q[],rvec mu,
		  t_topology *top,t_mdatoms *md,int gnx,atom_id grpindex[])
{
  int     i,start,end;
  real    mu_ave;
  t_atom  *atom;
  t_block *mols;
  
  start = START(nsb);
  end   = start + HOMENR(nsb);  
  
  atom = top->atoms.atom;
  mols = &(top->blocks[ebMOLS]);
  /*
  clear_rvec(mu);
  for(i=start; (i<end); i++)
    for(m=0; (m<DIM); m++)
      mu[m] += q[i]*x[i][m];
  if (PAR(cr)) {
    gmx_sum(DIM,mu,cr);
  }
  */
  /* I guess we have to parallelise this one! */

  if (gnx > 0) {
    mu_ave = 0.0;
    for(i=0; (i<gnx); i++) {
      int gi = grpindex[i];
      mu_ave += mol_dipole(mols->index[gi],mols->index[gi+1],mols->a,x,q);
    }
    
    return(mu_ave/gnx);
  }
  else
    return 0;
}

/* Lots of global variables! Yummy... */
static real tol     = 0.1;
static real epot    = 0.0;
static real fmax    = 100;
static real npow    = 12.0;
static bool bComb   = TRUE;
static bool bLogEps = FALSE;
static real ratio   = 0.01;

void set_ffvars(real ff_tol,real ff_epot,real ff_npow,bool ff_bComb,
		real ff_fmax,bool ff_bLogEps,real ff_ratio)
{
  tol     = ff_tol;
  epot    = ff_epot;
  npow    = ff_npow;
  bComb   = ff_bComb;
  fmax    = ff_fmax;
  ratio   = ff_ratio;
  bLogEps = ff_bLogEps;
}

real cost(real rmsf,real energy)
{
  return ratio*rmsf+fabs(energy-epot);
}

static char     *esenm[eseNR] = { "SIG", "EPS", "BHAMA", "BHAMB", "BHAMC" };
static int      nparm=0,*param_val=NULL;
static t_range  *range=NULL;
static t_genalg *ga=NULL;

static void init_range(t_range *r,int np,int atype,int ptype,real rmin,real rmax)
{
  if (rmin > rmax)
    fatal_error(0,"rmin (%f) > rmax (%f)",rmin,rmax);
  if (np <= 0)
    fatal_error(0,"np (%d) should be > 0",np);
  if ((rmax > rmin) && (np <= 1))
    fatal_error(0,"If rmax > rmin, np should be > 1");
  if ((ptype < 0) || (ptype >= eseNR))
    fatal_error(0,"ptype (%d) should be < %d",ptype,eseNR);
  r->np    = np;
  r->atype = atype;
  r->ptype = ptype;
  r->rmin  = rmin;
  r->rmax  = rmax;
  r->rval  = rmin;
  r->dr    = r->rmax - r->rmin;
}

static t_range *read_range(char *db,int *nrange)
{
  int     nlines,nr,np,i;
  char    **lines=NULL;
  t_range *range;
  int     atype,ptype;
  double  rmin,rmax;
  
  nlines = get_file(db,&lines);
  snew(range,nlines);
  
  nr=0;
  for(i=0; (i < nlines); i++) {
    strip_comment(lines[i]);
    if (sscanf(lines[i],"%d%d%d%lf%lf",&np,&atype,&ptype,&rmin,&rmax) == 5) {
      if (bLogEps && (ptype == eseEPSILON) && (rmin <= 0))
	fatal_error(0,"When using logarithmic epsilon increments the minimum"
		    "value must be > 0");
      init_range(&range[nr],np,atype,ptype,rmin,rmax);
      nr++;
    }
  }
  fprintf(stderr,"found %d variables to iterate over\n",nr);
  
  *nrange = nr;

  for(nr=0; (nr < nlines); nr++)
    sfree(lines[nr]);
  sfree(lines);
    
  return range;
}

static real value_range(t_range *r,int n)
{
  real logrmin,logrmax;
  
  if ((n < 0) || (n > r->np))
    fatal_error(0,"Value (%d) out of range for value_range (max %d)",n,r->np);

  if (r->np == 1)
    r->rval = r->rmin;
  else {
    if ((r->ptype == eseEPSILON) && bLogEps) {
      logrmin = log(r->rmin);
      logrmax = log(r->rmax);
      r->rval = exp(logrmin + (n*(logrmax-logrmin))/(r->np-1));
    }
    else
      r->rval = r->rmin+(n*(r->dr))/(r->np-1);
  }
  return r->rval;
}

real value_rand(t_range *r,int *seed)
{
  real logrmin,logrmax;
  real mr;
  
  if (r->np == 1)
    r->rval = r->rmin;
  else {
    mr = rando(seed);
    if ((r->ptype == eseEPSILON) && bLogEps) {
      logrmin = log(r->rmin);
      logrmax = log(r->rmax);
      r->rval = exp(logrmin + mr*(logrmax-logrmin));
    }
    else
      r->rval = r->rmin + mr*(r->rmax-r->rmin);
  }
  if (debug)
    fprintf(debug,"type: %s, value: %g\n",esenm[r->ptype],r->rval);
  return r->rval;
}

static void update_ff(t_forcerec *fr,int nparm,t_range range[],int param_val[])
{
  static double *sigma=NULL,*eps=NULL,*c6=NULL,*cn=NULL,*bhama=NULL,*bhamb=NULL,*bhamc=NULL;
  real   val,*nbfp;
  int    i,j,atnr;
  
  atnr = fr->ntype;
  nbfp = fr->nbfp;
  
  if (fr->bBHAM) {
    if (bhama == NULL) {
      snew(bhama,atnr);
      snew(bhamb,atnr);
      snew(bhamc,atnr);
    }
  }
  else {
    if (sigma == NULL) {
      snew(sigma,atnr);
      snew(eps,atnr);
      snew(c6,atnr);
      snew(cn,atnr);
    }
  }
  /* Get current values for everything */
  for(i=0; (i<nparm); i++) {
    if (ga)
      val = range[i].rval;
    else
      val = value_range(&range[i],param_val[i]);
    switch (range[i].ptype) {
    case eseSIGMA:
      sigma[range[i].atype] = val;
      break;
    case eseEPSILON:
      eps[range[i].atype] = val;
      break;
    case eseBHAMA:
      bhama[range[i].atype] = val;
      break;
    case eseBHAMB:
      bhamb[range[i].atype] = val;
      break;
    case eseBHAMC:
      bhamc[range[i].atype] = val;
      break;
    default:
      fatal_error(0,"Unknown ptype");
    }
  }
  if (fr->bBHAM) {
    for(i=0; (i<atnr); i++) {
      for(j=0; (j<=i); j++) {
	BHAMA(nbfp,atnr,i,j) = BHAMA(nbfp,atnr,j,i) = sqrt(bhama[i]*bhama[j]);
	BHAMB(nbfp,atnr,i,j) = BHAMB(nbfp,atnr,j,i) = sqrt(bhamb[i]*bhamb[j]);
	BHAMC(nbfp,atnr,i,j) = BHAMC(nbfp,atnr,j,i) = sqrt(bhamc[i]*bhamc[j]);
      }
    }
  }
  else {  
    /* Now build a new matrix */
    for(i=0; (i<atnr); i++) {
      c6[i] = 4*eps[i]*pow(sigma[i],6.0);
      cn[i] = 4*eps[i]*pow(sigma[i],npow);
      if (debug)
	fprintf(debug,"c6[%d] = %12.5e  c12[%d] = %12.5e\n",i,c6[i],i,cn[i]);
    }
    for(i=0; (i<atnr); i++) {
      for(j=0; (j<=i); j++) {
	C6(nbfp,atnr,i,j)  = C6(nbfp,atnr,j,i)  = sqrt(c6[i]*c6[j]);
	C12(nbfp,atnr,i,j) = C12(nbfp,atnr,j,i) = sqrt(cn[i]*cn[j]);
      }
    }
  }
  
  if (debug) {
    if (!fr->bBHAM) 
      for(i=0; (i<atnr); i++)
	fprintf(debug,"atnr = %2d  sigma = %8.4f  eps = %8.4f\n",i,sigma[i],eps[i]);
    for(i=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++) {
	if (fr->bBHAM)
	  fprintf(debug,"i: %2d  j: %2d  A:  %10.5e  B:  %10.5e  C:  %10.5e\n",i,j,
		  BHAMA(nbfp,atnr,i,j),BHAMB(nbfp,atnr,i,j),BHAMC(nbfp,atnr,i,j));
	else
	  fprintf(debug,"i: %2d  j: %2d  c6:  %10.5e  cn:  %10.5e\n",i,j,
		  C6(nbfp,atnr,i,j),C12(nbfp,atnr,i,j));
	
      }
    }
  }
}

void update_forcefield(int nfile,t_filenm fnm[],t_forcerec *fr)
{
  static int ntry,ntried;
  int    i,j;
  bool   bDone;

  /* First time around we have to read the parameters */  
  if (nparm == 0) {    
    range = read_range(ftp2fn(efDAT,nfile,fnm),&nparm);
    if (nparm == 0) 
      fatal_error(0,"No correct parameter info in %s",ftp2fn(efDAT,nfile,fnm));
    snew(param_val,nparm);

    if (opt2bSet("-ga",nfile,fnm)) {
      /* Genetic algorithm time */
      ga = init_ga(opt2fn("-ga",nfile,fnm),nparm,range);
    }
    else {  
      /* Determine the grid size */
      ntry = 1;
      for(i=0; (i<nparm); i++)
	ntry *= range[i].np;
      ntried = 0;
      
      fprintf(stdlog,"Going to try %d different combinations of %d parameters\n",
	      ntry,nparm);
    }
  }
  else if (ga)
    update_ga(stdlog,range,ga);
  else {
    /* Increment the counter
     * Non-trivial, since this is nparm nested loops in principle 
     */
    for(i=0; (i<nparm); i++) {
      if (param_val[i] < (range[i].np-1)) {
	param_val[i]++;
	for(j=0; (j<i); j++)
	  param_val[j] = 0;
	ntried++;
	break;
      }
    }
    if (i == nparm) {
      fprintf(stdlog,"Finished with %d out of %d iterations\n",ntried+1,ntry);
      if (gmx_parallel)
	gmx_finalize();
      exit(0);
    }
  }

  /* Now do the real updating */
  update_ff(fr,nparm,range,param_val);
}

static void print_range(FILE *fp,real rmsf,real energy)
{
  int  i;
  
  fprintf(fp,"%8.3f  %8.3f  %8.3f",cost(rmsf,energy),rmsf,energy);
  for(i=0; (i<nparm); i++)
    fprintf(fp," %s %10g",esenm[range[i].ptype],range[i].rval);
  fprintf(fp," FF\n");
  fflush(fp);
}

static void print_grid(FILE *fp,real energy,int natoms,rvec f[],rvec fshake[],
		       rvec x[],t_block *mols,real mass[])
{
  static bool bFirst = TRUE;
  static char *desc[] = {
    "------------------------------------------------------------------------",
    "In the output from the forcefield scan we have the potential energy,", 
    "then the root mean square force on the atoms, and finally the parameters",
    "in the order they appear in the input file.",
    "------------------------------------------------------------------------" 
  };
  int  i;
  real msf,rmsf;
  
  if (bFirst) {
    for(i=0; (i<asize(desc)); i++)
      fprintf(fp,"%s\n",desc[i]);
    fflush(fp);
    bFirst = FALSE;
  }
  if ((tol == 0) || (fabs(energy-epot) < tol)) {
    msf=0;
    for(i=0; (i<natoms); i++)
      msf += iprod(f[i],f[i]);
    rmsf = sqrt(msf/natoms);
    if ((fmax == 0) || (rmsf < fmax)) 
      print_range(fp,rmsf,energy);
  }
}

void print_forcefield(FILE *fp,real energy,int natoms,rvec f[],rvec fshake[],
		      rvec x[],t_block *mols,real mass[])
{
  int  i;
  real msf,rmsf;

  if (ga) {
    msf=0;
    for(i=0; (i<natoms); i++)
      msf += iprod(f[i],f[i]);
    rmsf = sqrt(msf/natoms);
    if (print_ga(fp,ga,rmsf,energy,range,tol)) {
      if (gmx_parallel)
	gmx_finalize();
      fprintf(stderr,"\n");
      exit(0);
    }
  }
  else
    print_grid(fp,energy,natoms,f,fshake,x,mols,mass);
}
 
