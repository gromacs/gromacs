/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <math.h>
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "vec.h"
#include "atomprop.h"
#include "xvgr.h"
#include "viewit.h"
#include "gmx_random.h"
#include "gstat.h"
#include "x2top_qgen.h"
#include "x2top_eemprops.h"

typedef struct {
  char    *molname;
  real    dip_ref,qtotal,dip_obs,chieq;
  t_atoms *atoms;
  rvec    *x;
  matrix  box;
} t_mymol; 

typedef struct {
  int     nmol,nparam,eemtype;
  int     *index;
  t_mymol *mymol;
  real    J0_0,Chi0_0,w_0,J0_1,Chi0_1,w_1,fc;
  void    *eem;
  void    *atomprop;
} t_moldip;

static void init_mymol(t_mymol *mymol,char *fn,real dip)
{
  int natoms;
  char title[STRLEN];
  
  /* Read coordinates */
  get_stx_coordnum(fn,&natoms); 
  snew(mymol->atoms,1);
  
  /* make space for all the atoms */
  init_t_atoms(mymol->atoms,natoms,TRUE);
  snew(mymol->x,natoms);              

  read_stx_conf(fn,title,mymol->atoms,mymol->x,NULL,mymol->box);

  mymol->molname = strdup(fn);
  mymol->dip_ref = dip;
  mymol->qtotal  = 0;
}

static void print_mols(FILE *logf,char *xvgfn,int nmol,t_mymol mol[],
		       int eemtype,void *eem)
{
  FILE   *xvgf;
  double d2=0;
  real   a,b,rms,sigma,error,aver;
  int    i,j,nout,*elemcnt;
  char   *resnm,*atomnm;
  t_lsq  lsq;
  
  xvgf  = xvgropen(xvgfn,"Correlation between dipoles",
		   "Experimental","Predicted");
  fprintf(xvgf,"@ s0 linestyle 0\n");
  fprintf(xvgf,"@ s0 symbol 1\n");
  init_lsq(&lsq);
  snew(elemcnt,109);
  for(i=0; (i<nmol); i++) {
    fprintf(logf,"Molecule %s, Dipole %6.2f, should be %6.2f <chi>: %6.2f\n",
	    mol[i].molname,mol[i].dip_obs,mol[i].dip_ref,mol[i].chieq);
    fprintf(xvgf,"%10g  %10g\n",mol[i].dip_ref,mol[i].dip_obs);
    add_lsq(&lsq,mol[i].dip_ref,mol[i].dip_obs);
    d2 += sqr(mol[i].dip_ref-mol[i].dip_obs);
    fprintf(logf,"Res  Atom  q\n");
    for(j=0; (j<mol[i].atoms->nr); j++) {
      resnm  = *(mol[i].atoms->resname[mol[i].atoms->atom[j].resnr]);
      atomnm = *(mol[i].atoms->atomname[j]);
      fprintf(logf,"%-5s%-5s  %8.4f\n",resnm,atomnm,mol[i].atoms->atom[j].q);
      elemcnt[eem_get_elem(eem,eem_get_index(eem,resnm,atomnm,eemtype))]++;
    }
    fprintf(logf,"\n");
  }
  fclose(xvgf);
  fprintf(logf,"Statistics over elements in the test set:\n");
  fprintf(logf,"Element   Number\n");
  for(i=0; (i<109); i++) {
    if (elemcnt[i] > 0)
      fprintf(logf,"%-10d%-8d\n",i,elemcnt[i]);
  }
  
  get_lsq_ab(&lsq,&a,&b);
  rms = sqrt(d2/nmol);
  fprintf(logf,"\nStatistics: fit of %d dipoles Dpred = %.3f Dexp + %3f\n",
	  nmol,a,b); 
  fprintf(logf,"RMSD = %.2f D\n",rms);
  aver = aver_lsq(&lsq);
  sigma = rms/aver;
  nout = 0;
  fprintf(logf,"Overview of outliers (> 1.5 sigma)\n");
  fprintf(logf,"----------------------------------\n");
  fprintf(logf,"%-20s  %12s  %12s  %12s\n",
	  "Name","Predicted","Experimental","Deviation");
  for(i=0; (i<nmol); i++) {
    if ((mol[i].dip_ref > 0) && 
	((mol[i].dip_obs/mol[i].dip_ref-1) > 1.5*sigma)) {
      fprintf(logf,"%-20s  %12.3f  %12.3f  %12.3f\n",
	      mol[i].molname,mol[i].dip_obs,mol[i].dip_ref,
	      fabs(mol[i].dip_obs-mol[i].dip_ref));
      nout ++;
    }
  }  
  if (nout)
    printf("There were %d outliers. See at the very bottom of the log file\n",
	   nout);
  else
    printf("No outliers! Well done.\n");
  do_view(xvgfn,NULL);
  
  done_lsq(&lsq);
}

t_moldip *read_moldip(FILE *logf,char *fn,char *eem_fn,
		      real J0_0,real Chi0_0,real w_0,
		      real J0_1,real Chi0_1,real w_1,
		      real fc,int eemtype)
{
  char     **strings,buf[STRLEN];
  int      i,n,nstrings;
  t_moldip *md;
  double   dip;
  
  nstrings = get_file(fn,&strings);
  snew(md,1);
  snew(md->mymol,nstrings);
  md->nmol = nstrings;
  for(i=0; (i<nstrings); i++) {
    if (sscanf(strings[i],"%s%lf",buf,&dip) != 2) 
      gmx_fatal(FARGS,"Error on line %d of %s",i+1,fn);
    init_mymol(&(md->mymol[i]),buf,dip);
  }
  fprintf(logf,"Read %d sets of molecule coordinates and dipoles\n",nstrings);
  md->eem = read_eemprops(eem_fn,-1);
  
  if ((md->eem == NULL) || ((md->nparam = eem_get_numprops(md->eem,eemtype)) == 0))
    gmx_fatal(FARGS,"Could not read %s, or file does not contain the requested parameters",
	      eem_fn ? eem_fn : "eemprops.dat");
  
  snew(md->index,md->nparam);
  n=0;
  for(i=1; (i<109); i++) {
    if ((md->index[n] = eem_get_elem_index(md->eem,i,eemtype)) != -1)
      n++;
  }
  if (n != md->nparam)
    gmx_fatal(FARGS,"Looking for the different elements, nparam = %d, n = %d",
	      md->nparam,n);
  md->atomprop   = get_atomprop();
  md->eemtype    = eemtype;
  md->J0_0       = J0_0;
  md->Chi0_0     = Chi0_0;
  md->w_0        = w_0;
  md->J0_1       = J0_1;
  md->Chi0_1     = Chi0_1;
  md->w_1        = w_1;
  md->fc         = fc;
  fprintf(logf,"There are %d atom types in the input file %s:\n---\n",
	  md->nparam,eem_fn ? eem_fn : "eemprops.dat");
  write_eemprops(logf,md->eem);
  fprintf(logf,"---\n\n");
	  
  return md;
}

static real mymol_calc_dip(t_mymol *mol)
{
  int i;
  rvec mu,mm;
  
  clear_rvec(mu);
  for(i=0; (i<mol->atoms->nr); i++) {
    svmul(mol->atoms->atom[i].q,mol->x[i],mm);
    rvec_inc(mu,mm);
  }
  return norm(mu)*ENM2DEBYE;
}


#ifdef HAVE_LIBGSL
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>

static double dipole_function(const gsl_vector *v,void *params)
{
  t_moldip *md = (t_moldip *) params;
  int      i,j,k;
  double   chi0,w,J0,rms,qt,qq;
  real     wj;
  
  rms = 0;
  
  /* Set parameters in eem record. There is a penalty if parameters
   * become negative as well.
   */
  k=0;
  for(i=0; (i<md->nparam); i++) {
    J0 = gsl_vector_get(v, k++);
    if (J0 < md->J0_0)
      rms += sqr(J0-md->J0_0);
    if (J0 > md->J0_1)
      rms += sqr(J0-md->J0_1);
    chi0 = gsl_vector_get(v, k++);
    if (chi0 < md->Chi0_0)
      rms += sqr(chi0-md->Chi0_0);
    if (chi0 > md->Chi0_1)
      rms += sqr(chi0-md->Chi0_1);
    if (md->eemtype != eqgSM1) {
      w = gsl_vector_get(v, k++);
      if (w < md->w_0)
	rms += sqr(w-md->w_0);
      if (w > md->w_1)
	rms += sqr(w-md->w_1);
    }
    else
      w = eem_get_w(md->eem,i);
    eem_set_props(md->eem,md->index[i],J0,w,chi0);
  }
  rms = rms*md->fc;
  
  for(i=0; (i<md->nmol); i++) {
    md->mymol[i].chieq = generate_charges_sm(debug,md->mymol[i].molname,
					     md->eem,md->mymol[i].atoms,
					     md->mymol[i].x,1e-4,10000,md->atomprop,
					     md->mymol[i].qtotal,md->eemtype);
    md->mymol[i].dip_obs = mymol_calc_dip(&(md->mymol[i]));
    rms += sqr(md->mymol[i].dip_obs - md->mymol[i].dip_ref);
    for(j=0; (j<md->mymol[i].atoms->nr); j++) {
      qq = md->mymol[i].atoms->atom[j].q;
      if (fabs(qq) >= 1)
	rms += (qq*qq-1);
    }
  }

  return sqrt(rms/md->nmol);
}

static real guess_new_param(real x,real step,real x0,real x1,gmx_rng_t rng,
			    bool bRandom)
{
  real r = gmx_rng_uniform_real(rng);
  
  if (bRandom) 
    x = x0+(x1-x0)*r;
  else
    x = x*(1-step+2*step*r);

  if (x < x0)
    return x0;
  else if (x > x1)
    return x1;
  else
    return x;
}

static real optimize_moldip(FILE *fp,FILE *logf,
			    t_moldip *md,int maxiter,real tol,
			    int reinit,real stepsize,int seed,
			    bool bRandom)
{
  FILE   *out,*xvg;
  real   size,d2,wj;
  int    iter   = 0;
  int    status = 0;
  int    i,k;
  bool   bRand;
  double J00,chi0,w;
  gmx_rng_t rng;
  
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  gsl_vector *x,*dx;
  gsl_multimin_function my_func;

  rng = gmx_rng_init(seed);
  
  my_func.f      = &dipole_function;
  my_func.n      = md->nparam*2;
  if (md->eemtype != eqgSM1)
    my_func.n += md->nparam;
  my_func.params = (void *) md;

  /* Starting point */
  x = gsl_vector_alloc (my_func.n);
  /* Step size, different for each of the parameters */
  dx = gsl_vector_alloc (my_func.n);
  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc (T, my_func.n);

  if (fp)
    fprintf(fp,"%5s %12s %12s\n","Iter","Size","RMS");
  
  do {
    if ((iter == 0) || ((reinit > 0) && ((iter % reinit) == 0))) {
      k=0;
      bRand = bRandom && (iter == 0);
      for(i=0; (i<md->nparam); i++) {
	J00 = lo_get_j00(md->eem,md->index[i],&wj,0);
	J00 = guess_new_param(J00,stepsize,md->J0_0,md->J0_1,rng,bRand);
       	gsl_vector_set (x, k, J00);
	gsl_vector_set (dx, k++, stepsize*J00);
	chi0 = eem_get_chi0(md->eem,md->index[i]);
	chi0 = guess_new_param(chi0,stepsize,md->Chi0_0,md->Chi0_1,rng,bRand);
	gsl_vector_set (x, k, chi0);
	gsl_vector_set (dx, k++, stepsize*chi0);
	
	w = eem_get_w(md->eem,md->index[i]);
	if (md->eemtype != eqgSM1) {
	  w = guess_new_param(w,stepsize,md->w_0,md->w_1,rng,bRand);
	  gsl_vector_set (x, k, w);
	  gsl_vector_set (dx, k++, stepsize*w);
	}
	eem_set_props(md->eem,md->index[i],J00,w,chi0);
      }
      gsl_multimin_fminimizer_set (s, &my_func, x, dx);
      if (0)
	write_eemprops(logf,md->eem);
    }
  
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);
    
    if (status != 0)
      gmx_fatal(FARGS,"Something went wrong in the iteration in minimizer %s",
		gsl_multimin_fminimizer_name(s));
    
    d2     = gsl_multimin_fminimizer_minimum(s);
    size   = gsl_multimin_fminimizer_size(s);
    /*    status = gsl_multimin_test_size(size,tol);
    
    if (status == GSL_SUCCESS)
      if (fp) 
	fprintf(fp,"Minimum found using %s\n",
		gsl_multimin_fminimizer_name(s));
    */    
    if (fp) 
      fprintf(fp,"%5d %12.4e %12.4e\n",iter,size,d2);
  } while (/*(status != GSL_SUCCESS) && */(sqrt(d2) > tol) && (iter < maxiter));
  
  gsl_vector_free (x);
  gsl_vector_free (dx);
  gsl_multimin_fminimizer_free (s);
  

  return d2;
}

static real quality_of_fit(real chi2,int N)
{
  return gsl_sf_gamma_inc_Q((N-2)/2.0,chi2/2.0);
}

#else
static real optimize_moldip(FILE *fp,t_moldip *md,int maxiter,real tol,
			    int reinit,real stepsize)
{
  fprintf(stderr,"This program needs the GNU scientific library to work.\n");
  
  return -1;
}

static real quality_of_fit(real chi2,int N)
{
  fprintf(stderr,"This program needs the GNU scientific library to work.\n");
  
  return -1;
}

#endif

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "tune_dip read a series of molecules and corresponding experimental",
    "dipole moments from a file, and tunes parameters in an algorithm",
    "until the experimental dipole moments are reproduces by the",
    "charge generating algorithm SM in the x2top program.[PAR]",
    "Minima and maxima for the parameters can be set, these are however",
    "not strictly enforced, but rather they are penalized with a harmonic",
    "function, for which the force constant can be set explicitly.[PAR]",
    "At every reinit step parameters are changed by a random amount within",
    "the fraction set by step size, and within the boundaries given",
    "by the minima and maxima."
  };
  
  t_filenm fnm[] = {
    { efDAT, "-f", "moldip",   ffREAD  },
    { efDAT, "-d", "eemprops", ffOPTRD },
    { efDAT, "-o", "molprops", ffWRITE },
    { efLOG, "-g", "charges",  ffWRITE },
    { efXVG, "-x", "dipcorr",  ffWRITE }
  };
#define NFILE asize(fnm)
  static int  maxiter=100,reinit=0,seed=1993;
  static real tol=1e-3;
  static bool bRandom=FALSE;
  static real J0_0=0,Chi0_0=0,w_0=0,step=0.01;
  static real J0_1=30,Chi0_1=30,w_1=1,fc=1.0;
  static char *qgen[] = { NULL, "SM1", "SM2", "SM3", "SM4", NULL };
  t_pargs pa[] = {
    { "-tol",   FALSE, etREAL, {&tol},
      "Tolerance for convergence in optimization" },
    { "-maxiter",FALSE, etINT, {&maxiter},
      "Max number of iterations for optimization" },
    { "-reinit", FALSE, etINT, {&reinit},
      "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
    { "-qgen",   FALSE, etENUM, {qgen},
      "Algorithm used for charge generation" },
    { "-j0",    FALSE, etREAL, {&J0_0},
      "Minimum value that J0 can obtain in fitting" },
    { "-chi0",    FALSE, etREAL, {&Chi0_0},
      "Minimum value that Chi0 can obtain in fitting" },
    { "-w0",    FALSE, etREAL, {&w_0},
      "Minimum value that Radius can obtain in fitting" },
    { "-j1",    FALSE, etREAL, {&J0_1},
      "Maximum value that J0 can obtain in fitting" },
    { "-chi1",    FALSE, etREAL, {&Chi0_1},
      "Maximum value that Chi0 can obtain in fitting" },
    { "-w1",    FALSE, etREAL, {&w_1},
      "Maximum value that Radius can obtain in fitting" },
    { "-fc",    FALSE, etREAL, {&fc},
      "Force constant in the penalty function for going outside the borders given with the above six option." },
    { "-step",  FALSE, etREAL, {&step},
      "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." },
    { "-seed", FALSE, etINT, {&seed},
      "Random number seed for reinit" },
    { "-random", FALSE, etBOOL, {&bRandom},
      "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step only." }
  };
  t_moldip *md;
  int      eemtype;
  FILE     *logf,*out;
      
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  eemtype = eqgSM1;
  if (qgen[0]) {
    eemtype = name2eemtype(qgen[0]);
    if (eemtype == -1)
      eemtype = eqgSM1;
  }
  if (eemtype > eqgSM3)
    gmx_fatal(FARGS,"Only models SM1 and SM2 implemented so far");
    
  logf = fopen(opt2fn("-g",NFILE,fnm),"w");
  md = read_moldip(logf,opt2fn("-f",NFILE,fnm),opt2fn_null("-d",NFILE,fnm),
		   J0_0,Chi0_0,w_0,J0_1,Chi0_1,w_1,fc,eemtype);
  
  (void) optimize_moldip(stdout,logf,md,maxiter,tol,reinit,step,seed,bRandom);
  
  print_mols(logf,opt2fn("-x",NFILE,fnm),md->nmol,md->mymol,
	     eemtype,md->eem);

  fclose(logf);
  
  out = fopen(opt2fn("-o",NFILE,fnm),"w");
  write_eemprops(out,md->eem);
  fclose(out);
  
  thanx(stderr);
  
  return 0;
}
