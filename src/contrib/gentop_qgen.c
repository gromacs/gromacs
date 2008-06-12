/* -*- mode: c; tab-width: 4; indent-tabs-mode: n; c-basic-offset: 4 -*- */
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
//#include "toppush.h"
//#include "pdb2top.h"
//#include "gen_ad.h"
//#include "topexcl.h"
#include "vec.h"
#include "atomprop.h"
#include "grompp.h"
//#include "add_par.h"
#include "gentop_nm2type.h"
#include "gentop_qgen.h"
#include "gentop_matrix.h"
#include "gentop_eemprops.h"
#include "slater_integrals.h"

typedef struct {
	int  natom,eemtype,slater_max;
	int  *eem_ndx; /* In the atoms, resp. eemprops arrays */
	int  *elem,*row;
	bool bWarned;
	real *chi0,*rhs,*qq,*j00,*wj,qtotal,chieq,hfac;
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
	real r,fac=10;
	real eNN=0,eSS=0,eSN=0,eTot;

	rvec_sub(xi,xj,dx);
	r = norm(dx);
	if (r == 0)
		gmx_fatal(FARGS,"Zero distance between atoms!\n");
	if (debug) {
		fprintf(debug,"Modifying input to Slater functions\n");
		r  *= fac;
		wi /= fac;
		wj /= fac;
	}
	switch (eemtype) {
	case eqgBultinck:
	case eqgSMp:
		eNN = coul_nucl_nucl(r);
		break;
	case eqgSMs:
	case eqgSMps:
		eSS = Coulomb_SS(r,rowi,rowj,wi,wj);
		break;
	case eqgRappe:
	case eqgYang:
		eSS = Coulomb_SS(r,rowi,rowj,wi,wj);
		break;
	case eqgSMg:
	case eqgSMpg:
	default:
		gmx_fatal(FARGS,"Can not treat algorithm %s yet in calc_jab",
				  get_eemtype_name(eemtype));
	}

	eTot = ONE_4PI_EPS0*(eNN+eSN+eSS)/ELECTRONVOLT;
	if (debug) {
		eTot *= fac;
	}
	return eTot;
}

static real calc_j1(rvec xi,rvec xj,real wi,real wj,
					int rowi,int rowj,int eemtype)
{
	rvec dx;
	real r;
	real eNN=0,eNS=0,eSN=0,eSS=0;

	rowi--;
	rowj--;
	range_check(rowi,0,SLATER_MAX);
	range_check(rowj,0,SLATER_MAX);
	rvec_sub(xi,xj,dx);
	r = norm(dx);
	if (r == 0)
		gmx_fatal(FARGS,"Zero distance between atoms!\n");
  
	eSN = Nuclear_SS(r,rowj,wj);
	eSS = Coulomb_SS(r,rowi,rowj,wi,wj);
  
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

static void qgen_update_J00(t_qgen *qgen,void *eem,int eemtype)
{
	int    i,j;
	double wi,wj,j0,j00,hj0,qq;
	double zetaH = 1.0698;
  
	for(i=0; (i<qgen->natom); i++) { 
		j0  = qgen->j00[i];
		
		/* Bohr is in nm */
#define BOHR  (0.052917)
		if (eem_get_elem(eem,qgen->eem_ndx[i]) == 1) {
			qq = qgen->qq[i];
			switch (eemtype) {
			case eqgYang:
				j0 = (1+qq/zetaH)*j0;
				break;
			case eqgRappe:
				j0 = (1+qq/zetaH)*j0;
				break;
			case eqgSMs:
			case eqgSMps:
			case eqgSMg:
			case eqgSMpg:
				j0 = (1+qgen->hfac*qq)*j0;
				break;
			default:
				break;
			}
		}
		if (debug && (j0 < 0) && !qgen->bWarned) {
			fprintf(debug,"WARNING: J00 = %g for atom %d. The equations will be instable.\n",j0,i+1);
			qgen->bWarned = TRUE;
		}
		qgen->Jab[i][i] = max(0,j0);
	}
}

static void qgen_debug(FILE *fp,t_qgen *qgen)
{
	int i;
	
	for(i=0; (i<qgen->natom); i++) {
		fprintf(fp,"QGEN: i: %2d chi0: %8g J0: %8g q: %8g\n",
				i+1,qgen->chi0[i],qgen->Jab[i][i],qgen->qq[i]);
	}
}

static real qgen_calc_Sij(t_qgen *qgen,int i,int j)
{
	real dist,dism,Sij = 1.0;
	rvec dx;
	int  l,m,tag;
	
	rvec_sub(qgen->x[i],qgen->x[j],dx);
	dist = norm(dx);
	if ((dist < 0.118) && (qgen->elem[i] != 1) && (qgen->elem[j] != 1)) {
		Sij = Sij*1.64;
	}
	else if ((dist < 0.122) && (qgen->elem[i] != 1) && (qgen->elem[j] != 1)) {
		if ((qgen->elem[i]  !=  8) && (qgen->elem[j]  !=  8)) {
			Sij = Sij*2.23;
		}
		else {
			Sij = Sij*2.23;
		}
	}
	else if (dist < 0.125) {
		tag=0;
		if ((qgen->elem[i] == 6) && (qgen->elem[j] == 8)) {
			tag=i;
		}
		else if ((qgen->elem[i] == 8) && (qgen->elem[j] == 6)) {
			tag=j;
		}
		if (tag != 0) {
			printf("found CO\n");
			for(l=0;(l<qgen->natom); l++) {
				if (qgen->elem[l] == 1) {
					printf("found H\n");
					dism=0.0;
					for(m=0; (m<DIM); m++) 
						dism=dism+sqr(qgen->x[tag][m]-qgen->x[l][m]);
					
					printf("dist: %8.3f\n",sqrt(dism));
					if (sqrt(dism) < 0.105) {
						printf("dist %5d %5d %5d  %5d %8.3f\n",
							   i,l,qgen->elem[tag],qgen->elem[l],sqrt(dism));
						Sij = Sij*1.605;
					}
				}
			}
		}
	}
	else if ((qgen->elem[i] == 6) && (qgen->elem[j] == 8)) 
		Sij = Sij*1.03;
	else if (((qgen->elem[j] == 6) && (qgen->elem[i] == 7) && 
			  (dist < 0.15)) ||
			 ((qgen->elem[i] == 6) && (qgen->elem[j] == 7)  &&
			  (dist < 0.15))) {
		if (qgen->elem[i] == 6) 
			tag=i;
		else
			tag=j;
		for(l=0; (l<qgen->natom); l++) {
			if (qgen->elem[l] == 8) {
				printf("found Oxy\n");
				dism=0.0;
				for(m=0; (m<DIM); m++) 
					dism=dism+sqr(qgen->x[tag][m]-qgen->x[l][m]);
				if (sqrt(dism) < 0.130) {
					printf("found peptide bond\n");
					Sij = Sij*0.66;
				}
				else
					Sij = Sij*1.1;
			}
		}
	}
	return Sij;
}

static void qgen_calc_Jab(t_qgen *qgen,void *eem,int eemtype)
{
	int    i,j;
	double Jab;
	
	for(i=0; (i<qgen->natom); i++) {
		for(j=i+1; (j<qgen->natom); j++) {
			Jab = calc_jab(qgen->x[i],qgen->x[j],qgen->wj[i],qgen->wj[j],
						   qgen->row[i],qgen->row[j],
						   qgen->eemtype);
			if (eemtype == eqgYang)
				Jab = Jab*qgen_calc_Sij(qgen,i,j);
			qgen->Jab[j][i] = qgen->Jab[i][j] = Jab;
		}
	}
}

static void qgen_calc_rhs(t_qgen *qgen,void *eem,int eemtype)
{
	int    i,j;
	double j1;
  
	for(i=0; (i<qgen->natom); i++) 
		qgen->rhs[i] = -qgen->chi0[i];
	if ((eemtype == eqgSMps) || (eemtype == eqgSMpg)) {
		for(i=0; (i<qgen->natom); i++) {
			for(j=i+1; (j<qgen->natom); j++) {
				j1 = calc_j1(qgen->x[i],qgen->x[j],qgen->wj[i],
							 qgen->wj[j],qgen->row[i],qgen->row[j],
							 qgen->eemtype);
				qgen->rhs[i] -= qgen->elem[j]*j1;
				qgen->rhs[j] -= qgen->elem[i]*j1;
			}
		}
	}
}

static t_qgen *qgen_init(void *eem,t_atoms *atoms,void *atomprop,
						 rvec *x,int eemtype,real hfac,
						 int slater_max)
{
	t_qgen *qgen;
	int i,j,atm;
  
	snew(qgen,1);
	qgen->eemtype    = eemtype;
	qgen->hfac       = hfac;
	qgen->slater_max = max(1,min(slater_max,SLATER_MAX));
	for(i=j=0; (i<atoms->nr); i++) 
		if (atoms->atom[i].ptype == eptAtom) 
			qgen->natom++;
	snew(qgen->chi0,qgen->natom);
	snew(qgen->rhs,qgen->natom);
	snew(qgen->elem,qgen->natom);
	snew(qgen->row,qgen->natom);
	snew(qgen->Jab,qgen->natom);
	snew(qgen->wj,qgen->natom);
	snew(qgen->j00,qgen->natom);
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
			qgen->row[j]  = min(qgen->slater_max,
								eem_get_row(eem,qgen->eem_ndx[j]));
			if (qgen->row[j] > SLATER_MAX) {
				if (debug)
					fprintf(debug,"Can not handle higher slaters than %d for atom %s %s\n",
							SLATER_MAX,
							*(atoms->resname[atoms->atom[i].resnr]),
							*(atoms->atomname[j]));
				qgen->row[j] = SLATER_MAX;
			}
			qgen->chi0[j] = eem_get_chi0(eem,qgen->eem_ndx[j]);
			qgen->wj[i]   = eem_get_zeta(eem,qgen->eem_ndx[i]);
			qgen->j00[i]  = eem_get_j00(eem,qgen->eem_ndx[i]);
			copy_rvec(x[i],qgen->x[j]);
			j++;
		}
	}  
  
	return qgen;
}

static void qgen_done(FILE *fp,t_atoms *atoms,t_qgen *qgen,
					  int eQGEN) 
{
	int  i,j,k,l,m;
	rvec mu = { 0, 0, 0 };
  
	if (eQGEN == eQGEN_OK) {
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
		if (fp)
			fprintf(fp,"<chieq> = %10g  <mu> = %10g\n",
					qgen->chieq,norm(mu)*ENM2DEBYE);
	}
	sfree(qgen->chi0);
	sfree(qgen->wj);
	sfree(qgen->j00);
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

void qgen_message(FILE *fp,int eQGEN)
{
	switch (eQGEN) {
	case eQGEN_OK:
		fprintf(fp,"Charge generation finished correctly.\n");
		break;
	case eQGEN_NOTCONVERGED:
		fprintf(fp,"Charge generation did not converge.\n");
		break;
	case eQGEN_ERROR:
	default:
		fprintf(fp,"Unknown status %d in charge generation.\n",eQGEN);
	}
}

int generate_charges_sm(FILE *fp,
						void *eem,t_atoms *atoms,rvec x[],
						real tol,int maxiter,void *atomprop,
						real qtotref,int eemtype,
						real hfac,int slater_max)
{
	t_qgen *qgen;
	real   *qq,chieq;
	int    i,j,iter,eQGEN;
	real   rms,mu;
  
	qgen = qgen_init(eem,atoms,atomprop,x,eemtype,hfac,slater_max);
	snew(qq,atoms->nr+1);
	for(i=j=0; (i<atoms->nr); i++)
		if (atoms->atom[i].ptype == eptShell) {
			qq[j] = qgen->qq[j];
			j++;
		}
	iter = 0;
	qgen_calc_Jab(qgen,eem,eemtype);
	qgen_calc_rhs(qgen,eem,eemtype);
	do {
		qgen_update_J00(qgen,eem,eemtype);
		if (debug)
			qgen_debug(debug,qgen);
		solve_q_eem(debug,qgen,1.0);
		rms = 0;
		for(i=j=0; (i<atoms->nr); i++) {
			if (atoms->atom[i].ptype != eptShell) {
				rms += sqr(qq[j] - qgen->qq[j]);
				qq[j] = qgen->qq[j];
				j++;
			}
		}
		rms = sqrt(rms/atoms->nr);
		iter++;
	} while ((rms > tol) && (iter < maxiter));
	
	if (iter < maxiter)
		eQGEN = eQGEN_OK;
	else
		eQGEN = eQGEN_NOTCONVERGED;
	
	if (fp)	{
		if (eQGEN == eQGEN_OK)
			fprintf(fp,"Converged to tolerance %g after %d iterations\n",
					tol,iter);
		else
			fprintf(fp,"Did not converge within %d iterations. RMS = %g\n",
					maxiter,rms);
	}
    
	qgen_done(fp,atoms,qgen,eQGEN);
	sfree(qgen);
	sfree(qq);
	
	return eQGEN;
}

static int generate_charges_bultinck(FILE *fp,
									 void *eem,t_atoms *atoms,
									 rvec x[],real tol,int maxiter,
									 void *atomprop)
{
	t_qgen *qgen;
	int    i,eQGEN=eQGEN_OK;
	real   rms;
  
	qgen = qgen_init(eem,atoms,atomprop,x,eqgBultinck,0,SLATER_MAX);
  
	qgen_calc_Jab(qgen,eem,eqgBultinck);
	qgen_calc_rhs(qgen,eem,eqgBultinck);
	qgen_update_J00(qgen,eem,eqgBultinck);
	solve_q_eem(debug,qgen,2.0);
  
	qgen_done(fp,atoms,qgen,eQGEN);
	sfree(qgen);
	
	return eQGEN;
}

int generate_charges(FILE *fp,char *molname,
					 int eemtype,t_atoms *atoms,rvec x[],
					 real tol,int maxiter,
					 void *atomprop,real qtotref,real hfac)
{
	int  i,eQGEN;
	void *eem;
	bool bConverged = FALSE;
	
	eem = read_eemprops(NULL,-1,atomprop);
	if (debug)
		write_eemprops(debug,eem);
  
	if (eem == NULL)
		gmx_fatal(FARGS,"Nothing interesting in eemprops.dat");
    
	/* Generate charges */
	please_cite(fp,get_eemtype_reference(eemtype));
	please_cite(stdout,get_eemtype_reference(eemtype));
	
	if (fp) 
		fprintf(fp,"Generating charges for %s using %s algorithm\n",
				molname,get_eemtype_name(eemtype));

	if (eemtype == eqgBultinck)
		eQGEN = generate_charges_bultinck(fp,eem,atoms,x,tol,
										  maxiter,atomprop);
	else
		eQGEN = generate_charges_sm(fp,eem,atoms,x,
									tol,maxiter,atomprop,
									qtotref,eemtype,hfac,SLATER_MAX);

	return eQGEN;
}

