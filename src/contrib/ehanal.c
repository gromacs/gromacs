/*
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "gromacs/utility/smalloc.h"
#include "macros.h"
#include "gromacs/utility/fatalerror.h"
#include "random.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "names.h"
#include "ehdata.h"
#include "gromacs/fileio/pdbio.h"

t_histo *init_histo(int np,real minx,real maxx)
{
  t_histo *h;
  
  snew(h,1);
  snew(h->y,np+1);
  snew(h->nh,np+1);
  h->np   = np;
  if (maxx <= minx)
    gmx_fatal(FARGS,"minx (%f) should be less than maxx (%f) in init_histo",minx,maxx);
  h->minx = minx;
  h->maxx = maxx;
  h->dx_1 = np/(maxx-minx);
  h->dx   = (maxx-minx)/np;
  
  return h;
}

void done_histo(t_histo *h)
{
  sfree(h->y);
  sfree(h->nh);
  h->np = 0;
}

void add_histo(t_histo *h,real x,real y)
{
  int n;
  
  n = (x-h->minx)*h->dx_1;
  if ((n < 0) || (n > h->np)) 
    gmx_fatal(FARGS,"Invalid x (%f) in add_histo. Should be in %f - %f",x,h->minx,h->maxx);
  h->y[n] += y;
  h->nh[n]++;
}

void dump_histo(t_histo *h,char *fn,char *title,char *xaxis,char *yaxis,
		int enorm,real norm_fac)
{
  FILE *fp;
  int  i,nn;
  
  for(nn=h->np; (nn > 0); nn--)
    if (h->nh[nn] != 0) 
      break;
  for(i=0; (i<nn); i++)
    if (h->nh[i] > 0)
      break;
  fp = xvgropen(fn,title,xaxis,yaxis);
  for(  ; (i<nn); i++) {
    switch (enorm) {
    case enormNO:
      fprintf(fp,"%12f  %12f  %d\n",h->minx+h->dx*i,h->y[i],h->nh[i]);
      break;
    case enormFAC:
      fprintf(fp,"%12f  %12f  %d\n",h->minx+h->dx*i,h->y[i]*norm_fac,h->nh[i]);
      break;
    case enormNP:
      if (h->nh[i] > 0)
	fprintf(fp,"%12f  %12f  %d\n",
		h->minx+h->dx*i,h->y[i]*norm_fac/h->nh[i],h->nh[i]);
      break;
    default:
      gmx_fatal(FARGS,"Wrong value for enorm (%d)",enorm);
    }
  }
  gmx_ffclose(fp);
}

/*******************************************************************
 *
 * Functions to analyse and monitor scattering
 *
 *******************************************************************/	

void add_scatter_event(t_ana_scat *scatter,rvec pos,gmx_bool bInel,
		       real t,real ekin)
{
  int np = scatter->np;
  
  if (np == scatter->maxp) {
    scatter->maxp += 32;
    srenew(scatter->time,scatter->maxp);
    srenew(scatter->ekin,scatter->maxp);
    srenew(scatter->bInel,scatter->maxp);
    srenew(scatter->pos,scatter->maxp);
  }
  scatter->time[np]  = t;
  scatter->bInel[np] = np;
  scatter->ekin[np]  = ekin;
  copy_rvec(pos,scatter->pos[np]);
  scatter->np++;
}

void reset_ana_scat(t_ana_scat *scatter)
{
  scatter->np = 0;
}

void done_scatter(t_ana_scat *scatter)
{
  scatter->np = 0;
  sfree(scatter->time);
  sfree(scatter->ekin);
  sfree(scatter->bInel);
  sfree(scatter->pos);
}

void analyse_scatter(t_ana_scat *scatter,t_histo *hmfp)
{
  int   i,n;
  rvec  dx;
  
  if (scatter->np > 1) {
    for(i=1; (i<scatter->np); i++) {
      rvec_sub(scatter->pos[i],scatter->pos[i-1],dx);
      add_histo(hmfp,scatter->ekin[i],norm(dx));
    }
  }
}

/*******************************************************************
 *
 * Functions to analyse structural changes
 *
 *******************************************************************/	

t_ana_struct *init_ana_struct(int nstep,int nsave,real timestep,
			      int maxparticle)
{
  t_ana_struct *anal;
  
  snew(anal,1);
  anal->nanal = 1.2*((nstep / nsave)+1);
  anal->index = 0;
  anal->dt    = nsave*timestep;
  snew(anal->t,anal->nanal);
  snew(anal->maxdist,anal->nanal);
  snew(anal->d2_com,anal->nanal);
  snew(anal->d2_origin,anal->nanal);
  snew(anal->nion,anal->nanal);
  anal->nstruct   = 1;
  anal->nparticle = 1;
  anal->maxparticle = maxparticle;
  snew(anal->q,1);
  snew(anal->x,1);
  snew(anal->x[0],maxparticle);
  
  return anal;
}

void done_ana_struct(t_ana_struct *anal)
{
  int i;
  
  sfree(anal->t);
  sfree(anal->maxdist);
  sfree(anal->d2_com);
  sfree(anal->d2_origin);
  sfree(anal->nion);
  sfree(anal->q);
  for(i=0; (i<anal->nstruct); i++)
    sfree(anal->x[i]);
  sfree(anal->x);
}

void reset_ana_struct(t_ana_struct *anal)
{
  int i;
  
  for(i=0; (i<anal->nanal); i++) {
    anal->t[i] = 0;
    anal->maxdist[i] = 0;
    clear_rvec(anal->d2_com[i]);
    clear_rvec(anal->d2_origin[i]);
    anal->nion[i] = 0;
  }
  anal->index = 0;
}

void add_ana_struct(t_ana_struct *total,t_ana_struct *add)
{
  int i,m;
  
  if (total->index == 0)
    total->index = add->index;
  else if (total->index != add->index)
    gmx_fatal(FARGS,"Analysis incompatible (total: %d, add: %d) %s, %d",
		total->index,add->index,__FILE__,__LINE__);
  for(i=0; (i<total->index); i++) {
    if (total->t[i] == 0)
      total->t[i] = add->t[i];
    else if (total->t[i] != add->t[i])
      gmx_fatal(FARGS,"Inconsistent times in analysis (%f-%f) %s, %d",
		  total->t[i],add->t[i],__FILE__,__LINE__);
    if (add->maxdist[i] > total->maxdist[i])
      total->maxdist[i]  = add->maxdist[i];
    for(m=0; (m<DIM); m++) {
      total->d2_com[i][m]    += add->d2_com[i][m]/add->nion[i];
      total->d2_origin[i][m] += add->d2_origin[i][m]/add->nion[i];
    }
    total->nion[i]     += add->nion[i];
  }
}

static void do_add_struct(t_ana_struct *anal,int nparticle,rvec x[])
{
  int i,j;
  
  if (nparticle > anal->nparticle) {
    for(i=0; (i<anal->nstruct); i++) {
      for(j=anal->nparticle; (j<nparticle); j++)
	copy_rvec(x[j],anal->x[i][j]);
    }
  }
  anal->nparticle=nparticle;
  srenew(anal->x,anal->nstruct+1);
  snew(anal->x[anal->nstruct],anal->maxparticle);
  for(j=0; (j<nparticle); j++)
    copy_rvec(x[j],anal->x[anal->nstruct][j]);
  anal->nstruct++;
}

void analyse_structure(t_ana_struct *anal,real t,rvec center,
		       rvec x[],int nparticle,real charge[])
{
  int  i,j,m,nel,n=0;
  rvec dx,com;
  real dx2,dx1;
  
  j = anal->index;
  if (j >= anal->nanal)
    gmx_fatal(FARGS,"Too many points in analyse_structure");
  anal->t[j]       = t;
  anal->maxdist[j] = 0;
  clear_rvec(com);
  nel = 0;
  for(i=0; (i<nparticle); i++) {
    if (charge[i] < 0) {
      rvec_inc(com,x[i]);
      nel++;
    }
  }
  if (nel > 0)
    for(m=0; (m<3); m++)
      com[m] /= nel;
  for(i=0; (i<nparticle); i++) {
    if (charge[i] < 0) {
      rvec_sub(x[i],com,dx);
      for(m=0; (m<DIM); m++) {
	anal->d2_com[j][m]    += sqr(dx[m]);
	anal->d2_origin[j][m] += sqr(x[i][m]);
      }
      dx2 = iprod(x[i],x[i]);
      dx1 = sqrt(dx2);
      if (dx1 > anal->maxdist[j])
	anal->maxdist[j] = dx1;
      n++;
    }
  }
  do_add_struct(anal,nparticle,x);
  anal->nion[j] = n;
  anal->index++;
}

void dump_ana_struct(char *rmax,char *nion,char *gyr_com,char *gyr_origin,
		     t_ana_struct *anal,int nsim)
{
  FILE *fp,*gp,*hp,*kp;
  int  i,j;
  real t,d2;
  char *legend[] = { "Rg", "RgX", "RgY", "RgZ" };
  
  fp = xvgropen(rmax,"rmax","Time (fs)","r (nm)");
  gp = xvgropen(nion,"N ion","Time (fs)","N ions");
  hp = xvgropen(gyr_com,"Radius of gyration wrt. C.O.M.",
		"Time (fs)","Rg (nm)");
  xvgr_legend(hp,asize(legend),legend);
  kp = xvgropen(gyr_origin,"Radius of gyration wrt. Origin",
		"Time (fs)","Rg (nm)");
  xvgr_legend(kp,asize(legend),legend);
  for(i=0; (i<anal->index); i++) {
    t = 1000*anal->t[i];
    fprintf(fp,"%12g  %10.3f\n",t,anal->maxdist[i]);
    fprintf(gp,"%12g  %10.3f\n",t,(1.0*anal->nion[i])/nsim-1);
    d2 = anal->d2_com[i][XX] + anal->d2_com[i][YY] + anal->d2_com[i][ZZ];
    fprintf(hp,"%12g  %10.3f  %10.3f  %10.3f  %10.3f\n",
	    t,sqrt(d2/nsim),
	    sqrt(anal->d2_com[i][XX]/nsim),
	    sqrt(anal->d2_com[i][YY]/nsim),
	    sqrt(anal->d2_com[i][ZZ]/nsim));
    d2 = anal->d2_origin[i][XX] + anal->d2_origin[i][YY] + anal->d2_origin[i][ZZ];
    fprintf(kp,"%12g  %10.3f  %10.3f  %10.3f  %10.3f\n",
	    t,sqrt(d2/nsim),
	    sqrt(anal->d2_origin[i][XX]/nsim),
	    sqrt(anal->d2_origin[i][YY]/nsim),
	    sqrt(anal->d2_origin[i][ZZ]/nsim));
  }
  gmx_ffclose(hp);
  gmx_ffclose(gp);
  gmx_ffclose(fp);
  gmx_ffclose(kp);
}

void dump_as_pdb(char *pdb,t_ana_struct *anal)
{
  FILE *kp;
  int  i,j;
  real t;
  
  kp = gmx_ffopen(pdb,"w");
  for(i=0; (i<anal->nstruct); i++) {
    t = 1000*anal->t[i];
    fprintf(kp,"MODEL  %d  time %g fs\n",i+1,t);
    for(j=0; (j<anal->nparticle); j++) {
      fprintf(kp,get_pdbformat(),"ATOM",i+1,(j < anal->nion[i]) ? "O" : "N",
	      "PLS",' ',1,
	      anal->x[i][j][XX]/100,
	      anal->x[i][j][YY]/100,
	      anal->x[i][j][ZZ]/100);
      fprintf(kp,"\n");
    }
    fprintf(kp,"ENDMDL\n");
  }
  gmx_ffclose(kp);
}

char *enms[eNR] = {
  "Coulomb", "Repulsion", "Potential",
  "EkHole",  "EkElectron", "EkLattice", "Kinetic",
  "Total"
};

void add_ana_ener(t_ana_ener *ae,int nn,real e[])
{
  int i;
 
  /* First time around we are constantly increasing the array size */ 
  if (nn >= ae->nx) {
    if (ae->nx == ae->maxx) {
      ae->maxx += 1024;
      srenew(ae->e,ae->maxx);
    }
    for(i=0; (i<eNR); i++)
      ae->e[ae->nx][i] = e[i];
    ae->nx++;
  }
  else {
    for(i=0; (i<eNR); i++)
      ae->e[nn][i] += e[i];
  }
}

void dump_ana_ener(t_ana_ener *ae,int nsim,real dt,char *edump,
		   t_ana_struct *total)
{
  FILE *fp;
  int  i,j;
  real fac;
  
  fac = 1.0/(nsim*ELECTRONVOLT);
  fp=xvgropen(edump,"Energies","Time (fs)","E (eV)");
  xvgr_legend(fp,eNR,enms);
  fprintf(fp,"@ s%d legend \"Ek/Nelec\"\n",eNR);
  fprintf(fp,"@ type nxy\n");
  for(i=0; (i<ae->nx); i++) {
    fprintf(fp,"%10f",1000.0*dt*i);
    for(j=0; (j<eNR); j++)
      fprintf(fp,"  %8.3f",ae->e[i][j]*fac);
    fprintf(fp,"  %8.3f\n",ae->e[i][eELECTRON]/(ELECTRONVOLT*total->nion[i]));
  }    
  fprintf(fp,"&\n");
  gmx_ffclose(fp);
}

