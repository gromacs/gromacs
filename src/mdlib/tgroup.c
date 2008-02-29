/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

static void init_grptcstat(int ngtc,t_grp_tcstat tcstat[])
{ 
  int i,j;
  
  for(i=0; (i<ngtc); i++) {
    tcstat[i].T=0;
    clear_mat(tcstat[i].ekin);
  }
}

static void init_grpstat(FILE *log,
			 t_mdatoms *md,int ngacc,t_grp_acc gstat[])
{
  int i,j,grp;

  /* First reset and count the number of atoms... */
  if (ngacc > 0) {
    for(j=0; (j<md->nr); j++) {
      grp=md->cACC[j];
      if ((grp < 0) && (grp >= ngacc))
	gmx_incons("Input for acceleration groups wrong");
      gstat[grp].nat++;
      gstat[grp].M+=md->massT[j];
    }
    
    for(i=0; (i<ngacc); i++) {
      /* Now allocate memory and fill the groups */
#ifdef DEBUG
      fprintf(log,"gstat[%d] has %d atoms\n",i,gstat[i].nat);
#endif
      snew(gstat[i].aid,gstat[i].nat);
      gstat[i].nat=0;
    }
    for(j=0; (j<md->nr); j++) {
      grp=md->cACC[j];
      gstat[grp].aid[gstat[grp].nat++]=j;
    }
#ifdef DEBUG
    for(i=0; (i<ngacc); i++) 
      fprintf(log,"gstat[%d] has %d atoms\n",i,gstat[i].nat);
#endif
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

void init_groups(FILE *log,t_mdatoms *md,t_grpopts *opts,t_groups *grps)
{
  int i;
#ifdef DEBUG
  fprintf(log,"ngtc: %d, ngacc: %d, ngener: %d\n",opts->ngtc,opts->ngacc,
	  opts->ngener);
#endif
  snew(grps->tcstat,opts->ngtc);
  init_grptcstat(opts->ngtc,grps->tcstat);
  
  snew(grps->grpstat,opts->ngacc);
  init_grpstat(log,md,opts->ngacc,grps->grpstat);
  
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
    add_binr(stdlog,rb,DIM,grps->grpstat[g].u);
    
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
		 t_grpopts *opts,rvec v[],t_mdatoms *md,bool bNEMD)
{
  int  d,g,n;
  real mv;

  /* calculate mean velocities at whole timestep */ 
  for(g=0; (g<opts->ngtc); g++) {
    grps->tcstat[g].T=0;
  }

  if (bNEMD) {
    for (g=0; (g<opts->ngacc); g++)
      clear_rvec(grps->grpstat[g].u);
    
    for(n=start; (n<start+homenr); n++) {
      g=md->cACC[n];
      for(d=0; (d<DIM);d++) {
	mv=md->massT[n]*v[n][d];
	grps->grpstat[g].u[d] += mv;
      }
    }

    for (g=0; (g < opts->ngacc); g++) {
      for(d=0; (d<DIM);d++) {
	grps->grpstat[g].u[d] /= grps->grpstat[g].M;
      }
    }
  }
}

real sum_ekin(bool bFirstStep,
	      t_grpopts *opts,t_groups *grps,tensor ekin,
	      real *dekindlambda)
{
  int          i,j,m,ngtc;
  real         T,ek;
  t_grp_tcstat *tcstat;
  real         nrdf,nd,*ndf;
  
  ngtc=opts->ngtc;
  tcstat=grps->tcstat;
  ndf=opts->nrdf;
  
  clear_mat(ekin);
  
  T=0; 
  nrdf=0;

  for(i=0; (i<ngtc); i++) {
    nd=ndf[i];
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
	copy_mat(tcstat[i].ekinh,tcstat[i].ekin);
      } else {
	/* Calculate the full step Ekin as the average of the half steps */
	for(j=0; (j<DIM); j++)
	  for(m=0; (m<DIM); m++)
	    tcstat[i].ekin[j][m] =
	      0.5*(tcstat[i].ekinh[j][m] + tcstat[i].ekinh_old[j][m]);
      }
      m_add(tcstat[i].ekin,ekin,ekin);
      ek=0;
      for(m=0; (m<DIM); m++)
	ek+=tcstat[i].ekinh[m][m];
      tcstat[i].Th=calc_temp(ek,nd);
      ek = 0;
      for(m=0; (m<DIM); m++) 
	ek+=tcstat[i].ekin[m][m];
      tcstat[i].T=calc_temp(ek,nd);
    }
    else {
      tcstat[i].T=0.0;
      tcstat[i].Th=0.0;
    }
    
    T    += nd*tcstat[i].T;
    nrdf += nd;
  }
  if (nrdf > 0)
    T/=nrdf;

  if (dekindlambda)
    *dekindlambda = 0.5*(grps->dekindl + grps->dekindl_old);
  
  return T;
}

static real sum_v(int n,real v[])
{
  real t;
  int  i;
  
  t=0.0;
  
  for(i=0; (i<n); i++)
    t=t+v[i];
    
  return t;
}

void sum_epot(t_grpopts *opts,t_groups *grps,real epot[])
{
  int i;

  /* Accumulate energies */
  epot[F_COUL_SR]  = sum_v(grps->estat.nn,grps->estat.ee[egCOULSR]);
  epot[F_LJ]       = sum_v(grps->estat.nn,grps->estat.ee[egLJSR]);
  epot[F_LJ14]     = sum_v(grps->estat.nn,grps->estat.ee[egLJ14]);
  epot[F_COUL14]   = sum_v(grps->estat.nn,grps->estat.ee[egCOUL14]);
  epot[F_COUL_LR] += sum_v(grps->estat.nn,grps->estat.ee[egCOULLR]);
  epot[F_LJ_LR]   += sum_v(grps->estat.nn,grps->estat.ee[egLJLR]);
/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
  epot[F_BHAM]     = sum_v(grps->estat.nn,grps->estat.ee[egBHAMSR]);
  epot[F_BHAM_LR]  = sum_v(grps->estat.nn,grps->estat.ee[egBHAMLR]);

  epot[F_EPOT] = 0;
  for(i=0; (i<F_EPOT); i++)
    if (i != F_DISRESVIOL && i != F_ORIRESDEV && i != F_DIHRESVIOL)
      epot[F_EPOT] += epot[i];
}
