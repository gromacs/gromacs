/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_tgroup_c = "$Id$";

#include <math.h>
#include "macros.h"
#include "main.h"
#include "smalloc.h"
#include "futil.h"
#include "tgroup.h"
#include "vec.h"
#include "network.h"
#include "smalloc.h"
#include "vveclib.h"
#include "update.h"
#include "rbin.h"
#include "assert.h"

static void init_grptcstat(int ngtc,t_grp_tcstat tcstat[])
{
  int i;
  
  for(i=0; (i<ngtc); i++) {
    tcstat[i].T=0;
    tcstat[i].lambda=1;
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
      assert((grp >= 0) && (grp < ngacc));
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
    fprintf(log,"%12.5e\n",estat->ee[egLJ][i]);
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
  static t_bin *rb=NULL;
  int          g;

  if (opts->ngacc < 1)
    return;
    
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
  int       g;

  for(g=0; (g<opts->ngtc); g++) 
    gmx_sum(DIM*DIM,grps->tcstat[g].ekin[0],cr);
}       

void update_grps(int start,int homenr,t_groups *grps,
		 t_grpopts *opts,rvec v[],t_mdatoms *md)
{
  int  d,g,n;
  real mv;

  /* calculate mean velocities at whole timestep */ 
  for(g=0; (g<opts->ngtc); g++) {
    grps->tcstat[g].T=0;
  }

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

real sum_ekin(t_grpopts *opts,t_groups *grps,tensor ekin,bool bTYZ)
{
  int          i,m,ngtc,nrdf,nd;
  real         T,ndfac,ek;
  rvec         tfac;
  t_grp_tcstat *tcstat;
  int          *ndf;
  
  ngtc=opts->ngtc;
  tcstat=grps->tcstat;
  ndf=opts->nrdf;
  
  clear_mat(ekin);
  
  T=0; 
  nrdf=0;
  if (bTYZ) {
    tfac[XX]=0.0;
    tfac[YY]=tfac[ZZ]=1.0;
    ndfac=2.0/3.0;
  }
  else {
    tfac[XX]=tfac[YY]=tfac[ZZ]=1.0;
    ndfac=1.0;
  }
  
  for(i=0; (i<ngtc); i++) {
    nd=ndf[i]*ndfac;
    /* Sometimes a group does not have degrees of freedom, e.g.
     * when it consists of shells and dummies, then we just
     * set the temperatue to 0 and also neglect the kinetic
     * energy, which should be  zero anyway.
     */
    if (nd > 0) {
      m_add(tcstat[i].ekin,ekin,ekin);
      ek=0;
      for(m=0; (m<DIM); m++) 
	ek+=tfac[m]*tcstat[i].ekin[m][m];
      tcstat[i].T=calc_temp(ek,nd);
    }
    else
      tcstat[i].T=0.0;
      
    T    += nd*tcstat[i].T;
    nrdf += nd;
  }
  if (nrdf > 0)
    T/=nrdf;
  
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
  /* Accumulate energies */
  epot[F_SR]   = sum_v(grps->estat.nn,grps->estat.ee[egCOUL]);
  epot[F_LJ]   = sum_v(grps->estat.nn,grps->estat.ee[egLJ]);
  epot[F_LJ14] = (sum_v(grps->estat.nn,grps->estat.ee[egCOUL14])+
		  sum_v(grps->estat.nn,grps->estat.ee[egLJ14]));
  epot[F_LR]   += sum_v(grps->estat.nn,grps->estat.ee[egLR]);
  epot[F_LJLR] += sum_v(grps->estat.nn,grps->estat.ee[egLJLR]);
/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
  epot[F_BHAM] = sum_v(grps->estat.nn,grps->estat.ee[egBHAM]);
  epot[F_EPOT] = sum_v(F_EPOT,epot);
}
