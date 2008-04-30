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
#include "x2top_nm2type.h"
#include "x2top_qgen.h"
#include "x2top_matrix.h"
#include "x2top_eemprops.h"
#include "slater_S_integrals.h"

typedef struct {
  int  natom,eemtype;
  int  *eem_ndx; /* In the atoms, resp. eemprops arrays */
  int  *elem,*row;
  real *chi0,*rhs,*qq,*wj,qtotal,chieq;
  real **Jab;
  rvec *x;
} t_qgen;

static real coul_nucl_nucl(real r)
{
  return 1/r;
}

static real calc_jab(rvec xi,rvec xj,real wi,real wj,
		     int rowi,int rowj,int eemtype)
{
  rvec dx;
  real r;
  real eNN=0,eSS=0,eSN=0,eTot;

  rvec_sub(xi,xj,dx);
  r = norm(dx);
  if (r == 0)
    gmx_fatal(FARGS,"Zero distance between atoms!\n");
  range_check(rowi,1,SLATER_MAX);
  range_check(rowj,1,SLATER_MAX);
  switch (eemtype) {
  case eqgBultinck:
  case eqgSMp:
    eNN = coul_nucl_nucl(r);
    break;
  case eqgSMs:
  case eqgSMps:
    eSS = Slater_SS[rowi-1][rowj-1](r,wi,wj);
    break;
  case eqgRappe:
  case eqgYang:
    eSS = Slater_SS[rowi-1][rowj-1](r,wi,wj);
    break;
  case eqgSMg:
  case eqgSMpg:
  default:
    gmx_fatal(FARGS,"Can not treat algorithm %s yet in calc_jab",
	      get_eemtype_name(eemtype));
  }

  eTot = ONE_4PI_EPS0*(eNN+eSN+eSS)/ELECTRONVOLT;
  
  return eTot;
}

static real calc_j1(rvec xi,rvec xj,real wi,real wj,
		    int rowi,int rowj,int eemtype)
{
  rvec dx;
  real r;
  real eNN=0,eNS=0,eSN=0,eSS=0;

  range_check(rowi,1,SLATER_MAX);
  range_check(rowj,1,SLATER_MAX);
  rvec_sub(xi,xj,dx);
  r = norm(dx);
  if (r == 0)
    gmx_fatal(FARGS,"Zero distance between atoms!\n");
  
  eSN = Slater_NS[rowj-1](r,wj);
  eSS = Slater_SS[rowi-1][rowj-1](r,wi,wj);
  
  return ONE_4PI_EPS0*(eNN+eSN-eSS)/ELECTRONVOLT;
}

static void solve_q_eem(FILE *fp,t_qgen *qgen,real hardness_factor)
{
  double **a,**b,qtot;
  int i,j,n,nn;

  n = qgen->natom+1;
  a = alloc_matrix(n,n);
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

  matrix_invert(fp,n,a);
  qtot = 0;  
  for(i=0; (i<n-1); i++) {
    qgen->qq[i] = 0;
    for(j=0; (j<n-1); j++) {
      qgen->qq[i] += a[i][j]*qgen->rhs[j];
    }
    qtot += qgen->qq[i];
  }
  qgen->chieq = 0;
  for(j=0; (j<n-1); j++) 
    qgen->chieq += a[n-1][j]*qgen->rhs[j];
  
  if (fabs(qtot - qgen->qtotal) > 1e-2)
    fprintf(stderr,"qtot = %g, it should be %g\n",qtot,qgen->qtotal);
  free_matrix(a,n);
}

static void qgen_calc_Jab(t_qgen *qgen,void *eem,int eemtype)
{
  int    i,j;
  double wi,wj;
  
  for(i=0; (i<qgen->natom); i++) { 
    qgen->Jab[i][i] = eem_get_j00(eem,qgen->eem_ndx[i],
				  &(qgen->wj[i]),qgen->qq[i]);
    if (debug)
      fprintf(debug,"CALC_JAB: Atom %d, Weight: %g, J0: %g, chi0: %g\n",
	      i,qgen->wj[i],qgen->Jab[i][i],qgen->chi0[i]);
  }
  for(i=0; (i<qgen->natom); i++) {
    wi           =  qgen->wj[i];
    qgen->rhs[i] = -qgen->chi0[i];
    for(j=0; (j<qgen->natom); j++) {
      if (i != j) {
	wj = qgen->wj[j];
	qgen->Jab[i][j] = calc_jab(qgen->x[i],qgen->x[j],wi,wj,
				   qgen->row[i],qgen->row[j],qgen->eemtype);
	if ((eemtype == eqgSMps) || (eemtype == eqgSMpg))
	  qgen->rhs[i] -= qgen->elem[j]*calc_j1(qgen->x[i],qgen->x[j],
						wi,wj,qgen->row[i],qgen->row[j],
						qgen->eemtype);
      }
    }
  }
}

t_qgen *init_qgen(void *eem,t_atoms *atoms,void *atomprop,rvec *x,int eemtype)
{
  t_qgen *qgen;
  int i,j,atm;
  
  snew(qgen,1);
  qgen->eemtype = eemtype;
  for(i=j=0; (i<atoms->nr); i++) 
    if (atoms->atom[i].ptype == eptAtom) 
      qgen->natom++;
  snew(qgen->chi0,qgen->natom);
  snew(qgen->rhs,qgen->natom);
  snew(qgen->elem,qgen->natom);
  snew(qgen->row,qgen->natom);
  snew(qgen->Jab,qgen->natom);
  snew(qgen->wj,qgen->natom);
  snew(qgen->eem_ndx,qgen->natom);
  snew(qgen->qq,qgen->natom);
  snew(qgen->x,qgen->natom);
  for(i=j=0; (i<atoms->nr); i++) {
    if (atoms->atom[i].ptype == eptAtom) {
      snew(qgen->Jab[j],qgen->natom);
      qgen->eem_ndx[j] = eem_get_index(eem,atoms->atom[i].atomnumber,
				       qgen->eemtype);
      if (qgen->eem_ndx[j] == -1)
	gmx_fatal(FARGS,"No electronegativity data for %s %s and algorithm %s",
		  *(atoms->resname[atoms->atom[i].resnr]),
		  *(atoms->atomname[j]),
		  get_eemtype_name(eemtype));
      qgen->elem[j] = eem_get_elem(eem,qgen->eem_ndx[j]);
      qgen->row[j]  = eem_get_row(eem,qgen->eem_ndx[j]);
      qgen->chi0[j] = eem_get_chi0(eem,qgen->eem_ndx[j]);
      copy_rvec(x[i],qgen->x[j]);
      j++;
    }
  }  
  
  return qgen;
}

static void done_qgen(FILE *fp,t_atoms *atoms,t_qgen *qgen) 
{
  int  i,j,k,l,m;
  rvec mu = {0,0,0};
  
  if (fp)
    fprintf(fp,"Res Atom   Nr Row           q        chi0      weight\n");
  for(i=j=0; (i<atoms->nr); i++) {
    if (atoms->atom[i].ptype == eptAtom) {
      atoms->atom[i].q = qgen->qq[j];
      for(m=0; (m<DIM); m++)
	mu[m] += atoms->atom[i].q * qgen->x[i][m];
      if (fp)
	fprintf(fp,"%4s%4s%5d%4d  %10g  %10g  %10g\n",
		*(atoms->resname[atoms->atom[i].resnr]),
		*(atoms->atomname[i]),i+1,qgen->row[j],
		qgen->qq[j],qgen->chi0[j],qgen->wj[j]);
      j++;
    }
  }
  if (0 && fp) {
    fprintf(fp,"Jab matrix:\n");
    for(i=k=0; (i<atoms->nr); i++) {
      if (atoms->atom[i].ptype == eptAtom) {      
	for(j=l=0; (j<atoms->nr); j++) {
	  if (atoms->atom[i].ptype == eptAtom) {
	    fprintf(fp,"  %6.2f",qgen->Jab[k][l]);
	    l++;
	  }
	}
	k++;
	fprintf(fp,"\n");
      }
    }
  }
  if (fp)
    fprintf(fp,"<chieq> = %10g  <mu> = %10g\n",
	    qgen->chieq,norm(mu)*ENM2DEBYE);
  sfree(qgen->chi0);
  sfree(qgen->wj);
  sfree(qgen->qq);
  sfree(qgen->eem_ndx);
  sfree(qgen->elem);
  sfree(qgen->rhs);
  sfree(qgen->x);
  for(i=j=0; (i<atoms->nr); i++) 
    if (atoms->atom[i].ptype == eptAtom) 
      sfree(qgen->Jab[j++]);
  
  sfree(qgen->Jab);
}

static void generate_charges_yang_rappe(void *eem,t_atoms *atoms,rvec x[],
					real tol,int maxiter,void *atomprop,
					int eqg,real hardness)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  real   *qq;
  int    i,iter;
  real   rms;
  
  if (eqg == eqgYang)
    printf("Generating charges using Yang & Sharp algorithm\n");
  else
    printf("Generating charges using Rappe & Goddard algorithm\n");
  qgen = init_qgen(eem,atoms,atomprop,x,eqg);
  snew(qq,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    qq[i] = qgen->qq[i];
  iter = 0;
  do {
    qgen_calc_Jab(qgen,eem,eqg);
    solve_q_eem(debug,qgen,hardness);
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

real generate_charges_sm(FILE *fp,char *molname,
			 void *eem,t_atoms *atoms,rvec x[],
			 real tol,int maxiter,void *atomprop,
			 real qtotref,int eemtype,
			 real hardness)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  real   *qq,chieq;
  int    i,j,iter;
  real   rms,mu;
  
  if (fp)
    fprintf(fp,"Generating charges using Van der Spoel & Van Maaren algorithm for %s\n",molname);
  qgen = init_qgen(eem,atoms,atomprop,x,eemtype);
  snew(qq,atoms->nr+1);
  for(i=j=0; (i<atoms->nr); i++)
    if (atoms->atom[i].ptype == eptShell) {
      qq[j] = qgen->qq[j];
      j++;
    }
  iter = 0;
  do {
    qgen_calc_Jab(qgen,eem,eemtype);
    solve_q_eem(fp,qgen,hardness);
    rms = 0;
    for(i=j=0; (i<atoms->nr); i++) {
      if (atoms->atom[i].ptype == eptShell) {
	rms += sqr(qq[j] - qgen->qq[j]);
	qq[j] = qgen->qq[j];
	j++;
      }
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
  sfree(qq);
  chieq = qgen->chieq;
  done_qgen(fp,atoms,qgen);
  sfree(qgen);
  
  return chieq;
}

static void generate_charges_bultinck(void *eem,t_atoms *atoms,rvec x[],
				      real tol,int maxiter,void *atomprop,
				      real hardness)
{
  /* Use Rappe and Goddard derivative for now */
  t_qgen *qgen;
  int    i;
  real   rms;
  
  printf("Generating charges using Bultinck algorithm\n");
  qgen = init_qgen(eem,atoms,atomprop,x,eqgBultinck);
  
  qgen_calc_Jab(qgen,eem,eqgBultinck);
  solve_q_eem(NULL,qgen,hardness);
  
  done_qgen(stdout,atoms,qgen);
}

void assign_charges(char *molname,
		    int eemtype,t_atoms *atoms,rvec x[],
		    t_params *bonds,real tol,real fac,int maxiter,
		    void *atomprop,real qtotref,real hardness)
{
  int  i;
  void *eem;
  
  eem = read_eemprops(NULL,-1,atomprop);
  if (debug)
    write_eemprops(debug,eem);
  
  if (eem == NULL)
    gmx_fatal(FARGS,"Nothing interesting in eemprops.dat");
    
  /* Generate charges */
  please_cite(stdout,get_eemtype_reference(eemtype));
  switch (eemtype) {
  case eqgNone:
    for(i=0; (i<atoms->nr); i++) {
      atoms->atom[i].q  = atoms->atom[i].qB = 0;
    }
    break;
  case eqgYang:
  case eqgRappe:
    generate_charges_yang_rappe(eem,atoms,x,tol,maxiter,atomprop,eemtype,
				hardness);
    break;
  case eqgBultinck:
    generate_charges_bultinck(eem,atoms,x,tol,maxiter,atomprop,
			      hardness);
    break;
  case eqgSMp:
  case eqgSMpp:
  case eqgSMs:
  case eqgSMps:
  case eqgSMg:
  case eqgSMpg:
    (void) generate_charges_sm(debug,molname,eem,atoms,x,tol,maxiter,atomprop,
			       qtotref,eemtype,hardness);
    break;
  default:
    gmx_fatal(FARGS,"Algorithm %s out of range in assign_charge_alpha",
	      get_eemtype_name(eemtype));
  }
}

