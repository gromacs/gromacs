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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <math.h>
#include "macros.h"
#include "main.h"
#include "smalloc.h"
#include "futil.h"
#include "tgroup.h"
#include "vec.h"
#include "network.h"
#include "smalloc.h"
#include "update.h"
#include "rbin.h"
#include "mtop_util.h"

static void init_grptcstat(int ngtc,t_grp_tcstat tcstat[])
{ 
  int i,j;
  
  for(i=0; (i<ngtc); i++) {
    tcstat[i].T = 0;
    clear_mat(tcstat[i].ekin);
  }
}

static void init_grpstat(FILE *log,
			 gmx_mtop_t *mtop,int ngacc,t_grp_acc gstat[])
{
  gmx_groups_t *groups;
  gmx_mtop_atomloop_all_t aloop;
  int    i,grp;
  t_atom *atom;

  if (ngacc > 0) {
    groups = &mtop->groups;
    aloop = gmx_mtop_atomloop_all_init(mtop);
    while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) {
      grp = ggrpnr(groups,egcACC,i);
      if ((grp < 0) && (grp >= ngacc))
	gmx_incons("Input for acceleration groups wrong");
      gstat[grp].nat++;
      /* This will not work for integrator BD */
      gstat[grp].mA += atom->m;
      gstat[grp].mB += atom->mB;
    }
  }
}

static void init_grpener(FILE *log,int ngener,t_grp_ener *estat)
{
  int i,n2;
  
  n2=ngener*ngener;
#ifdef DEBUG
  fprintf(log,"Creating %d sized group matrix for energies\n",n2);
#endif
  estat->nn=n2;
  for(i=0; (i<egNR); i++) {
    snew(estat->ee[i],n2);
  }
}

void init_t_groups(FILE *log,gmx_mtop_t *mtop,t_grpopts *opts,t_groups *grps)
{
  int i;
#ifdef DEBUG
  fprintf(log,"ngtc: %d, ngacc: %d, ngener: %d\n",opts->ngtc,opts->ngacc,
	  opts->ngener);
#endif
  snew(grps->tcstat,opts->ngtc);
  init_grptcstat(opts->ngtc,grps->tcstat);
   /* Set Berendsen tcoupl lambda's to 1, 
   * so runs without Berendsen coupling are not affected.
   */
  for(i=0; i<opts->ngtc; i++) {
    grps->tcstat[i].lambda = 1.0;
  }
  
  snew(grps->grpstat,opts->ngacc);
  init_grpstat(log,mtop,opts->ngacc,grps->grpstat);
 
  init_grpener(log,opts->ngener,&grps->estat);
}

void dump_estat(FILE *log,t_grp_ener *estat)
{
  int i;
  
  for(i=0; (i<estat->nn); i++) {
    fprintf(log,"%12.5e\n",estat->ee[egLJSR][i]);
  }
}

real rms_ener(t_energy *e,int nsteps)
{
  real erms2;

  erms2=e->e2sum*nsteps-e->esum*e->esum;
  if (erms2 <= 0)
    return 0.0;
  else
    return sqrt(erms2)/nsteps;
}

void accumulate_u(t_commrec *cr,t_grpopts *opts,t_groups *grps)
{
  /* This routine will only be called when it's necessary */
  static t_bin *rb=NULL;
  int          g;

  if (rb == NULL) {
    rb=mk_bin();
  }
  else
    reset_bin(rb);

  for(g=0; (g<opts->ngacc); g++) 
    add_binr(rb,DIM,grps->grpstat[g].u);
    
  sum_bin(rb,cr);
  
  for(g=0; (g<opts->ngacc); g++) 
    extract_binr(rb,DIM*g,DIM,grps->grpstat[g].u);
}       

static void accumulate_ekin(t_commrec *cr,t_grpopts *opts,t_groups *grps)
{
  int g;

  if(PAR(cr))
    for(g=0; (g<opts->ngtc); g++) 
      gmx_sum(DIM*DIM,grps->tcstat[g].ekin[0],cr);
}       

void update_grps(int start,int homenr,t_groups *grps,
		 t_grpopts *opts,rvec v[],t_mdatoms *md,real lambda,
		 bool bNEMD)
{
  int  d,g,n;
  real mv;

  /* calculate mean velocities at whole timestep */ 
  for(g=0; (g<opts->ngtc); g++) {
    grps->tcstat[g].T = 0;
  }

  if (bNEMD) {
    for (g=0; (g<opts->ngacc); g++)
      clear_rvec(grps->grpstat[g].u);
    
    g = 0;
    for(n=start; (n<start+homenr); n++) {
      if (md->cACC)
	g = md->cACC[n];
      for(d=0; (d<DIM);d++) {
	mv = md->massT[n]*v[n][d];
	grps->grpstat[g].u[d] += mv;
      }
    }

    for (g=0; (g < opts->ngacc); g++) {
      for(d=0; (d<DIM);d++) {
	grps->grpstat[g].u[d] /=
	  (1-lambda)*grps->grpstat[g].mA + lambda*grps->grpstat[g].mB;
      }
    }
  }
}

real sum_ekin(bool bFirstStep,
	      t_grpopts *opts,t_groups *grps,
	      tensor ekin,real *dekindlambda)
{
  int          i,j,m,ngtc;
  real         T,ek;
  t_grp_tcstat *tcstat;
  real         nrdf,nd,*ndf;
  
  ngtc = opts->ngtc;
  ndf  = opts->nrdf;
  
  clear_mat(ekin);
  
  T = 0; 
  nrdf = 0;

  for(i=0; (i<ngtc); i++) {
    tcstat = &grps->tcstat[i];
    nd = ndf[i];
    /* Sometimes a group does not have degrees of freedom, e.g.
     * when it consists of shells and virtual sites, then we just
     * set the temperatue to 0 and also neglect the kinetic
     * energy, which should be  zero anyway.
     */
    if (nd > 0) {
      if (bFirstStep) {
	/* This Ekin is only used for reporting the initial temperature
	 * or when doing mdrun -rerun.
	 */
	copy_mat(tcstat->ekinh,tcstat->ekin);
      } else {
	/* Calculate the full step Ekin as the average of the half steps */
	for(j=0; (j<DIM); j++)
	  for(m=0; (m<DIM); m++)
	    tcstat->ekin[j][m] =
	      0.5*(tcstat->ekinh[j][m] + tcstat->ekinh_old[j][m]);
      }
      m_add(tcstat->ekin,ekin,ekin);
      ek = 0;
      for(m=0; (m<DIM); m++)
	ek += tcstat->ekinh[m][m];
      tcstat->Th = calc_temp(ek,nd);
      ek = 0;
      for(m=0; (m<DIM); m++) 
	ek += tcstat->ekin[m][m];
      tcstat->T = calc_temp(ek,nd);
    }
    else {
      tcstat->T  = 0;
      tcstat->Th = 0;
    }
    
    T    += nd*tcstat->T;
    nrdf += nd;
  }
  if (nrdf > 0)
    T/=nrdf;

  if (dekindlambda)
    *dekindlambda = 0.5*(grps->dekindl + grps->dekindl_old);

  return T;
}
