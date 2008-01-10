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
#include "toppush.h"
#include "pdb2top.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "vec.h"
#include "atomprop.h"
#include "grompp.h"
#include "add_par.h"
#include "gmx_random.h"
#include "x2top_nm2type.h"
#include "x2top_qgen.h"
#include "x2top_matrix.h"
#include "x2top_eemprops.h"

typedef struct {
  int  natom,eemtype;
  int  *index; /* In the Yang array */
  real *chi,*chi0,*qq,*wj,qtotal;
  real **Jab;
  rvec *x;
  real chiav;
} t_qgen;

static real coul_slater_slater(real w,real r)
{
  const real a = 11.0/16.0;
  const real b =  3.0/16.0;
  const real c =  1.0/48.0;
  real r_w  = r/w;
  real r_w2 = r_w*r_w;
  real r_w3 = r_w2*r_w;
  
  return (1/r)*(1 - (1+a*r_w+b*r_w2+c*r_w3)*exp(-r_w));
}

static real coul_slater_nucl(real w,real r) 
{
  real r_w  = r/w;
  
  return (1/r)*(1-(1+0.5*r_w)*exp(-r_w));
}

static real coul_nucl_nucl(real w,real r)
{
  return 1/r;
}

static real calc_jab(rvec xi,rvec xj,real wi,real wj)
{
  rvec dx;
  real r,wij;
  real e0=0,e1=0,e2=0;

  rvec_sub(xi,xj,dx);
  r = norm(dx);
  
  e0 = coul_nucl_nucl(0,r);
  if ((wi > 0) && (wj > 0)) {
    wij = 2/(1/wi + 1/wj);  /* Dit kan geoptimaliseerd worden */
    e1 = coul_slater_nucl(wij,r);
    e2 = coul_slater_slater(wij,r);
  }

  return ONE_4PI_EPS0*(e0-2*e1+e2)/ELECTRONVOLT;
}

static real get_chi0(void *atomprop,char *resnm,char *name)
{
  real value;
  
  if (!query_atomprop(atomprop,epropElectroneg,resnm,name,&value))
    if (debug)
      fprintf(debug,"Using default electronegativity value %g for %s:%s\n",
	      value,resnm,name);
  return value;
}

static void solve_q_eem(FILE *fp,t_qgen *qgen,real hardness_factor)
{
  double **a,**b,qtot,chieq;
  int i,j,n,nn;

  n = qgen->natom+1;
  a = alloc_matrix(n,n);
  b = alloc_matrix(n,n);
  for(i=0; (i<n-1); i++) {
    for(j=0; (j<n-1); j++) {
      a[i][j] = qgen->Jab[i][j];
    }
    a[i][i] = hardness_factor*qgen->Jab[i][i];
  }
  for(j=0; (j<n-1); j++)
    a[n-1][j] = 1;
  for(i=0; (i<n-1); i++) 
    a[i][n-1] = -1;
  a[n-1][n-1] = 0;

  for(i=0; (i<n); i++) 
    for(j=0; (j<n); j++) 
      b[i][j] = a[i][j];
  
  mat_inv(fp,n,a);
  qtot = 0;  
  for(i=0; (i<n-1); i++) {
    qgen->qq[i] = 0;
    for(j=0; (j<n-1); j++) {
      qgen->qq[i] += -a[i][j]*qgen->chi0[j];
    }
    qtot += qgen->qq[i];
  }
  chieq = 0;
  for(i=0; (i<n-1); i++) {
    qgen->chi[i] = 0;
    for(j=0; (j<n-1); j++) {
      qgen->chi[i] += b[i][j]*qgen->qq[j];
    }
    chieq += qgen->chi[i];
  }
  qgen->chiav = chieq/qgen->natom;
  if (fabs(qtot - qgen->qtotal) > 1e-3)
    fprintf(stderr,"qtot = %g, it should be %g\n",qtot,qgen->qtotal);
  free_matrix(a,n);
  free_matrix(b,n);
}

static void qgen_calc_Jab(t_qgen *qgen,void *eem)
{
  int    i,j;
  double wi,wj;
  
  for(i=0; (i<qgen->natom); i++) { 
    qgen->Jab[i][i] = lo_get_j00(eem,qgen->index[i],&(qgen->wj[i]),qgen->qq[i]);
  }
  for(i=0; (i<qgen->natom); i++) {
    wi = qgen->wj[i];
    for(j=0; (j<qgen->natom); j++) {
      if (i != j) {
	wj = qgen->wj[j];
	qgen->Jab[i][j] = qgen->Jab[j][i] = 
	  calc_jab(qgen->x[i],qgen->x[j],wi,wj);
      }
    }
  }
}

t_qgen *init_qgen(void *eem,t_atoms *atoms,void *atomprop,rvec *x,int eemtype)
{
  t_qgen *qgen;
  gmx_rng_t rng;
  int i,j;
  
  rng = gmx_rng_init(17);
  snew(qgen,1);
  qgen->natom   = atoms->nr;
  qgen->eemtype = eemtype;
  snew(qgen->chi0,atoms->nr);
  snew(qgen->chi,atoms->nr);
  snew(qgen->Jab,atoms->nr);
  snew(qgen->wj,atoms->nr);
  snew(qgen->index,atoms->nr);
  snew(qgen->qq,atoms->nr);
  qgen->x = x;
  for(i=0; (i<atoms->nr); i++) {
    snew(qgen->Jab[i],atoms->nr);
    qgen->qq[i] = 0.1*(gmx_rng_uniform_real(rng)-0.5);
    qgen->index[i] = eem_getindex(eem,
				  *(atoms->resname[atoms->atom[i].resnr]),
				  *(atoms->atomname[i]),qgen->eemtype);
    if (qgen->index[i] == -1)
      gmx_fatal(FARGS,"Can not find index for %s %s",
		*(atoms->resname[atoms->atom[i].resnr]),
		*(atoms->atomname[i]));
    qgen->chi0[i] = eem_get_chi0(eem,qgen->index[i]);
  }  
  
  return qgen;
}

static void done_qgen(FILE *fp,t_atoms *atoms,t_qgen *qgen) 
{
  int i,j;
  
  if (fp)
    fprintf(fp,"Res Atom   Nr           q         chi        chi0       weight\n");
  for(i=0; (i<atoms->nr); i++) {
    atoms->atom[i].q = qgen->qq[i];
    if (fp)
      fprintf(fp,"%4s%4s%5d  %10g  %10g  %10g  %10g\n",
	      *(atoms->resname[atoms->atom[i].resnr]),
	      *(atoms->atomname[i]),i+1,
	      qgen->qq[i],qgen->chi[i],qgen->chi0[i],qgen->wj[i]);
  }
  if (0 && fp) {
    fprintf(fp,"Jab matrix:\n");
    for(i=0; (i<atoms->nr); i++) {
      for(j=0; (j<atoms->nr); j++) {
	fprintf(fp,"  %6.2f",qgen->Jab[i][j]);
      }
      fprintf(fp,"\n");
    }
  }
  if (fp)
    fprintf(fp,"<chieq> = %10g\n",qgen->chiav);
  sfree(qgen->chi0);
  sfree(qgen->chi);
  sfree(qgen->wj);
  sfree(qgen->qq);
  sfree(qgen->x);
  for(i=0; (i<atoms->nr); i++) 
    sfree(qgen->Jab[i]);
  sfree(qgen->Jab);
}

static void generate_charges_yang(void *eem,t_atoms *atoms,rvec x[],
				  real tol,int maxiter,void *atomprop)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  real   *qq;
  int    i,iter;
  real   rms;
  
  printf("Generating charges using Yang & Sharp algorithm\n");
  qgen = init_qgen(eem,atoms,atomprop,x,eqgYang);
  snew(qq,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    qq[i] = qgen->qq[i];
  iter = 0;
  do {
    qgen_calc_Jab(qgen,eem);
    solve_q_eem(NULL,qgen,2.0);
    rms = 0;
    for(i=0; (i<atoms->nr); i++) {
      rms += sqr(qq[i] - qgen->qq[i]);
      qq[i] = qgen->qq[i];
    }
    rms = sqrt(rms/atoms->nr);
    iter++;
  } while ((rms > tol) && (iter < maxiter));
  if (iter < maxiter)
    printf("Converged to tolerance %g after %d iterations\n",tol,iter);
  else
    printf("Did not converge with %d iterations. RMS = %g\n",maxiter,rms);
    
  done_qgen(stdout,atoms,qgen);
}

void generate_charges_sm(FILE *fp,
			 void *eem,t_atoms *atoms,rvec x[],
			 real tol,int maxiter,void *atomprop,
			 real qtotref)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  real   *qq;
  int    i,iter;
  real   rms,mu;
  
  if (fp)
    fprintf(fp,"Generating charges using Van der Spoel & Van Maaren algorithm\n");
  qgen = init_qgen(eem,atoms,atomprop,x,eqgSM);
  snew(qq,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    qq[i] = qgen->qq[i];
  iter = 0;
  do {
    qgen_calc_Jab(qgen,eem);
    solve_q_eem(NULL,qgen,2.0);
    rms = 0;
    for(i=0; (i<atoms->nr); i++) {
      rms += sqr(qq[i] - qgen->qq[i]);
      qq[i] = qgen->qq[i];
    }
    rms = sqrt(rms/atoms->nr);
    iter++;
  } while ((rms > tol) && (iter < maxiter));
  
  if (fp) {
    if (iter < maxiter)
      fprintf(fp,"Converged to tolerance %g after %d iterations\n",tol,iter);
    else
      fprintf(fp,"Did not converge with %d iterations. RMS = %g\n",maxiter,rms);
  }
  done_qgen(fp,atoms,qgen);
  sfree(qq);
  sfree(qgen);
}

static void generate_charges_bultinck(void *eem,t_atoms *atoms,rvec x[],
				      real tol,int maxiter,void *atomprop)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  int    i;
  real   rms;
  
  printf("Generating charges using Bultinck algorithm\n");
  qgen = init_qgen(eem,atoms,atomprop,x,eqgBultinck);
  
  qgen_calc_Jab(qgen,eem);
  solve_q_eem(NULL,qgen,2.0);
  
  done_qgen(stdout,atoms,qgen);
}

static void generate_charges_linear(t_atoms *atoms,rvec x[],t_params *bonds,
				    real tol,real fac,int maxiter,
				    void *atomprop)
{
  /* Novel electronegativity method */
  real   *chi,*chi0;
  int    i,ai,aj,iter;
  double rmsd,msd,dq,chiav,chisum,chi2sum;
  
  fprintf(stderr,"Generating charges using Linear algorithm\n");
  snew(chi,atoms->nr);
  snew(chi0,atoms->nr);
  for(i=0; (i<atoms->nr); i++) {
    chi0[i] = get_chi0(atomprop,*(atoms->resname[atoms->atom[i].resnr]),
		       *atoms->atomname[i]);
    chi[i]  = chi0[i];
    atoms->atom[i].q = 0;
  }
  /* Perform shake like algorithm to optimize q */ 
  iter = 0;
  do {
    for(i=0; (i<bonds->nr); i++) {
      ai = bonds->param[i].AI;
      aj = bonds->param[i].AJ;
      dq = (chi[aj]-chi[ai])/10;
      atoms->atom[ai].q += dq;
      atoms->atom[aj].q -= dq;
      chi[ai] = chi0[ai] + fac*atoms->atom[ai].q;
      chi[aj] = chi0[aj] + fac*atoms->atom[aj].q;
    }
    chi2sum = chisum = 0;
    for(i=0; (i<atoms->nr); i++) {
      chi2sum += (chi[i]*chi[i]);
      chisum  += chi[i];
    }
    chiav = chisum/atoms->nr;
    msd   = chi2sum/atoms->nr - sqr(chiav);
    rmsd  = sqrt(msd);
    iter++;
    fprintf(stderr,"iter: %5d rms: %g, <chi>: %g, q0: %g, chi[0]: %g\n",
	    iter,rmsd,chiav,atoms->atom[0].q,chi[0]);
  } while ((rmsd > tol) && (iter < maxiter));
  
  sfree(chi);
  sfree(chi0);
}

void assign_charge_alpha(int alg,t_atoms *atoms,rvec x[],
			 t_params *bonds,real tol,real fac,int maxiter,
			 void *atomprop,real qtotref)
{
  int  i;
  void *eem;
  
  eem = read_eemprops(NULL);
  if (debug)
    write_eemprops(debug,eem);
  
  if ((eem == NULL) && (alg > eqgLinear))
    gmx_fatal(FARGS,"Nothing interesting in eemprops.dat");
    
  /* Generate charges */
  switch (alg) {
  case eqgNone:
    for(i=0; (i<atoms->nr); i++) {
      atoms->atom[i].q  = atoms->atom[i].qB = 0;
    }
    break;
  case eqgLinear:
    generate_charges_linear(atoms,x,bonds,tol,fac,maxiter,atomprop);
    break;
  case eqgYang:
    please_cite(stdout,"Yang2006b");
    generate_charges_yang(eem,atoms,x,tol,maxiter,atomprop);
    break;
  case eqgBultinck:
    please_cite(stdout,"Bultinck2002a");
    generate_charges_bultinck(eem,atoms,x,tol,maxiter,atomprop);
    break;
  case eqgSM:
    generate_charges_sm(stdout,eem,atoms,x,tol,maxiter,atomprop,qtotref);
    break;
  default:
    gmx_fatal(FARGS,"Algorithm %d out of range in assign_charge_alpha",alg);
  }
}

