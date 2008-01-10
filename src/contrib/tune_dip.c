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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
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
#include "grompp.h"
#include "gmx_random.h"
#include "x2top_qgen.h"
#include "x2top_eemprops.h"

typedef struct {
  char    *molname;
  real    dip_ref,qtotal,dip_obs;
  t_atoms *atoms;
  rvec    *x;
  matrix  box;
} t_mymol; 

typedef struct {
  int     nmol,nparam;
  t_mymol *mymol;
  bool    bFitJ0,bFitChi0,bFitRadius;
  real    J0_0,Chi0_0,R_0;
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
  mymol->qtotal = 0;
}

static void print_mol(FILE *fp,t_mymol *mol)
{
  int i;
  
  fprintf(fp,"Molecule %s, Dipole %6.2f, should be %6.2f\n",
	  mol->molname,mol->dip_obs,mol->dip_ref);
  fprintf(fp,"Res  Atom  q\n");
  for(i=0; (i<mol->atoms->nr); i++)
    fprintf(fp,"%-5s%-5s  %8.4f\n",
	    *(mol->atoms->resname[mol->atoms->atom[i].resnr]),
	    *(mol->atoms->atomname[i]),mol->atoms->atom[i].q);
  fprintf(fp,"\n");
}

t_moldip *read_moldip(char *fn,bool bFitJ0,bool bFitChi0,bool bFitRadius,
		      real J0_0,real Chi0_0,real R_0)
{
  char     **strings,buf[STRLEN];
  int      i,nstrings;
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
  printf("Read %d sets of molecular coordinates and dipoles\n",nstrings);
  md->eem = read_eemprops(NULL);
  if (md->eem == NULL)
    gmx_fatal(FARGS,"Could not read eemprops.dat");
  md->nparam     = eem_getnumprops(md->eem);
  md->atomprop   = get_atomprop();
  md->bFitJ0     = bFitJ0;
  md->bFitChi0   = bFitChi0;
  md->bFitRadius = bFitRadius;
  md->J0_0       = J0_0;
  md->Chi0_0     = Chi0_0;
  md->R_0        = R_0;
  
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
  double   chi0,radius,J0,rms,qt,qq;
  real     wj;
  
  rms = 0;
  
  /* Set parameters in eem record. There is a penalty if parameters
   * become negative as well.
   */
  k=0;
  for(i=0; (i<md->nparam); i++) {
    if (md->bFitJ0) {
      J0 = gsl_vector_get(v, k++);
      if (J0 <= md->J0_0)
	rms += sqr(J0-md->J0_0);
    }
    else
      J0 = lo_get_j00(md->eem,i,&wj,0);
    
    if (md->bFitChi0) {
      chi0 = gsl_vector_get(v, k++);
      if (chi0 <= md->Chi0_0)
	rms += sqr(chi0-md->Chi0_0);
    }
    else
      chi0 = eem_get_chi0(md->eem,i);
      
    if (md->bFitRadius) {
      radius = gsl_vector_get(v, k++);
      if (radius <= md->R_0)
	rms += sqr(radius-md->R_0);
    }
    else
      radius = eem_get_radius(md->eem,i);
    eem_set_props(md->eem,i,J0,radius,chi0);
  }
    
  for(i=0; (i<md->nmol); i++) {
    generate_charges_sm(debug,md->eem,md->mymol[i].atoms,
			md->mymol[i].x,1e-4,10000,md->atomprop,
			md->mymol[i].qtotal);
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

static real optimize_moldip(FILE *fp,t_moldip *md,int maxiter,real tol,
			    char *outfn,char *logfn)
{
  FILE   *out;
  real   size,d2,wj;
  int    iter   = 0;
  int    status = 0;
  int    i,k;
  double J00,chi0,radius;
  
  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;

  gsl_vector *x,*dx;
  gsl_multimin_function my_func;

  my_func.f      = &dipole_function;
  my_func.n      = md->nparam*(md->bFitJ0 + md->bFitChi0 + md->bFitRadius);
  my_func.params = (void *) md;

  /* Starting point */
  x = gsl_vector_alloc (my_func.n);
  /* Step size, different for each of the parameters */
  dx = gsl_vector_alloc (my_func.n);
  k=0;
  for(i=0; (i<md->nparam); i++) {
    if (md->bFitJ0) {
      J00 = lo_get_j00(md->eem,i,&wj,0);
      gsl_vector_set (x, k, J00);
      gsl_vector_set (dx, k++, 0.01*J00);
    }
    if (md->bFitChi0) {
      chi0 = eem_get_chi0(md->eem,i);
      gsl_vector_set (x, k, chi0);
      gsl_vector_set (dx, k++, 0.01*chi0);
    }
    if (md->bFitRadius) {
      radius = eem_get_radius(md->eem,i);
      gsl_vector_set (x, k, radius);
      gsl_vector_set (dx, k++, 0.01*radius);
    }
  }
  
  
  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc (T, my_func.n);

  gsl_multimin_fminimizer_set (s, &my_func, x, dx);
  gsl_vector_free (x);
  gsl_vector_free (dx);

  if (fp)
    fprintf(fp,"%5s %12s %12s\n","Iter","Size","RMS");
  
  do  {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);
    
    if (status != 0)
      gmx_fatal(FARGS,"Something went wrong in the iteration in minimizer %s",
		gsl_multimin_fminimizer_name(s));
    
    d2     = gsl_multimin_fminimizer_minimum(s);
    size   = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size,tol);
    
    /*if (status == GSL_SUCCESS)
      if (fp) 
	fprintf(fp,"Minimum found using %s\n",
		gsl_multimin_fminimizer_name(s));
    */
    if (fp) {
      fprintf(fp,"%5d", iter);
      fprintf(fp," %12.4e %12.4e\n",size,d2);
    }
  }
  while (/*(status != GSL_SUCCESS) &&*/ (sqrt(d2) > tol) && (iter < maxiter));
  
  gsl_multimin_fminimizer_free (s);
  
  out = fopen(outfn,"w");
  write_eemprops(out,md->eem);
  fclose(out);

  out = fopen(logfn,"w");
  for(i=0; (i<md->nmol); i++) {
    print_mol(out,&(md->mymol[i]));
  }
  fclose(out);
    
  return d2;
}

static real quality_of_fit(real chi2,int N)
{
  return gsl_sf_gamma_inc_Q((N-2)/2.0,chi2/2.0);
}

#else
static real optimize_moldip(FILE *fp,t_qgen *qgen,int maxiter,real tol)
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
    "charge generating algorithm SM in the x2top program."
  };
  
  t_filenm fnm[] = {
    { efDAT, "-f", "moldip",   ffREAD  },
    { efDAT, "-o", "molprops", ffWRITE },
    { efLOG, "-g", "tune_dip", ffWRITE },
  };
#define NFILE asize(fnm)
  static int  maxiter=100;
  static real tol=1e-3;
  static real J0_0=0,Chi0_0=0,R_0=0;
  static bool bFitJ0=TRUE,bFitChi0=TRUE,bFitRadius=FALSE;
  t_pargs pa[] = {
    { "-tol",   FALSE, etREAL, {&tol},
      "Tolerance for convergence in optimization" },
    { "-maxiter",FALSE, etINT, {&maxiter},
      "Max number of iterations for optimization" },
    { "-fitj0", FALSE, etBOOL, {&bFitJ0},
      "Optimize dipoles by fitting J00" },
    { "-j0",    FALSE, etREAL, {&J0_0},
      "Minimum value that J0 can obtain in fitting" },
    { "-fitChi0", FALSE, etBOOL, {&bFitChi0},
      "Optimize dipoles by fiting Chi0" },
    { "-chi0",    FALSE, etREAL, {&Chi0_0},
      "Minimum value that Chi0 can obtain in fitting" },
    { "-fitradius", FALSE, etBOOL, {&bFitRadius},
      "Optimize dipole by fitting Radius" },
    { "-r0",    FALSE, etREAL, {&R_0},
      "Minimum value that Radius can obtain in fitting" }
  };
  t_moldip *md;
    
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  if (!bFitJ0 && !bFitChi0 && !bFitRadius)
    gmx_fatal(FARGS,"Nothing to fit to!");
    		    
  md = read_moldip(opt2fn("-f",NFILE,fnm),bFitJ0,bFitChi0,bFitRadius,
		   J0_0,Chi0_0,R_0);
  
  (void) optimize_moldip(stdout,md,maxiter,tol,
			 opt2fn("-o",NFILE,fnm),
			 opt2fn("-g",NFILE,fnm));
    
  thanx(stderr);
  
  return 0;
}
