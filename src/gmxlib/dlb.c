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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_dlb_c = "$Id$";
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "nsb.h"
#include "futil.h"
#include "invblock.h"
#include "network.h"
#include "xvgr.h"
#include "main.h"
#include "nrnb.h"
#include "smalloc.h"

static real low_count_nb(t_nblist *nlist,atom_id cgnr[],
			 real weight,real w2,real nb[])
{
  int    i,ia;
  real   inj,nj,ntot;
  t_nl_i *nli;

  ntot  = 0;
  nli   = nlist->nl_i;
  for(i=0; (i<nlist->nri); i++) {
    ia=nli[i].i_atom;
    inj=nli[i].nj;
    if (nli[i].bWater)
      nj=w2*inj;
    else
      nj=weight*inj;
    ntot+=nj;
    nb[cgnr[ia]] += nj;
  }
  
  return ntot;
}

real calc_bfload(int nbonds,
		 t_iatom forceatoms[],
		 t_functype functype[],
		 t_iparams forceparams[])
                       
{
  int  i,nat,type,ftype,ind;
  real vtot;
  
  vtot=0;
  for (i=0; i<nbonds; ) {
    type=forceatoms[i++];
    ftype=functype[type];
    ind=interaction_function[ftype].nrnb_ind;
    if (ind != -1)
      vtot+=cost_nrnb(ind);
    nat=interaction_function[ftype].nratoms;
    i+=nat;
  }
  return vtot;
}

void count_nb(t_commrec *cr,t_nsborder *nsb,t_block *cgs,int nns,int nlr,
	      t_idef *idef,int ngener)
{
  static FILE    *out=NULL;
  static FILE    *o_ld=NULL;
  static atom_id *cgnr;
  static int     step=0;
  static real    *bfload;
  real           *nb;
  
  real   ntp,ntot,load[MAXNODES+1],load_av,unb,tld,ratio,defcost;
  int    nload[MAXNODES+1];
  int    ncg,nodeid,i,j,k;
  
  if (nsb->nnodes < 2)
    return;
  
  ncg=nsb->cgtotal;
  snew(nb,ncg+1);
   
  nload[0]=0;
  for(i=1; (i<=cr->nnodes); i++)
    nload[i]=nsb->workload[i-1];
    
  /* Initiate the cost per charge group,
   * there are nns particles for this node, 
   */
  defcost=0.1*(nns*cost_nrnb(eNR_NS)+nlr*cost_nrnb(eNR_LR))/
    (nload[cr->nodeid+1]-nload[cr->nodeid]);
  for(i=0; (i<nload[cr->nodeid]); i++)
    nb[i]=0;
  for(; (i<nload[cr->nodeid+1]); i++)
    nb[i]=defcost;
  for(; (i<ncg+1); i++)
    nb[i]=0;
  ntot=(nload[cr->nodeid+1]-nload[cr->nodeid])*defcost;
      
  if (cgnr==NULL) {
    cgnr=make_invblock(cgs,nsb->natoms);
    
    /* Load from bonded forces is constant... */
    snew(bfload,cr->nnodes+1);
    bfload[cr->nodeid]=calc_bfload(idef->il[eilBONDEDS].nr,
				idef->il[eilBONDEDS].iatoms,
				idef->functype,idef->iparams);
    fprintf(stdlog,"My    Bonded Force Load: %.0f\n",bfload[cr->nodeid]);
    bfload[cr->nnodes]=bfload[cr->nodeid];
    gmx_sum(cr->nnodes+1,bfload,cr);
    fprintf(stdlog,"Total Bonded Force Load: %.0f\n",bfload[cr->nnodes]);
  }
  
  for(i=0; (i<ngener); i++) {
    ntot += low_count_nb(&nl_nbfp[i],cgnr,
			 cost_nrnb(eNR_LJC),cost_nrnb(eNR_WATER),nb);
    /*ntot += low_count_nb(&nl_14[i],  cgnr,
			 cost_nrnb(eNR_LJC),0,nb);*/
    ntot += low_count_nb(&nl_c[i],   cgnr,
			 cost_nrnb(eNR_QQ),0,nb);
  }
  
  fprintf(stdlog,"My    Nonbonded Force Load =%.0f\n",ntot);
  nb[ncg]=ntot;
  gmx_sum(ncg+1,nb,cr);
  ntot=nb[ncg];
  fprintf(stdlog,"Total Nonbonded Force Load =%.0f\n",ntot);
  ntot+=bfload[cr->nnodes];
  fprintf(stdlog,"Total           Force Load =%.0f\n",ntot);
  
  load_av=ntot/nsb->nnodes;
  unb=1;
  for(nodeid=0; (nodeid<nsb->nnodes); nodeid++) {
    load[nodeid]=bfload[nodeid];
    for(j=nload[nodeid]; (j<nload[nodeid+1]); j++) 
      load[nodeid]+=nb[j];
#ifdef DEBUG
    fprintf(stdlog,"nodeid: %2d, cg0: %5d, hcgHi: %3d, load: %10g\n",
	    nodeid,nsb->cg0[nodeid],nsb->hcgHi[nodeid],load[nodeid]);
#endif
    ratio=load[nodeid]/load_av;
    unb=max(unb,ratio);
  }
  tld=0;
  nodeid=1;
  nload[0]=0;
  for(i=0; ((i<ncg) && (nodeid < nsb->nnodes)); i++) {
    tld+=nb[i];
    if (tld >= (nodeid*load_av)) {
      nload[nodeid]=i;
      nodeid++;
    }
  }
  if (nodeid == nsb->nnodes)
    nload[nodeid]=ncg;
  else {
    fprintf(stdlog,"nodeid(%d) != nnodes (%d)\n",nodeid,nsb->nnodes);
  }
  if (nsb->nodeid == 0) {
    if (out==NULL) {
      out=xvgropen("load.xvg","Maximum Load","Step","% of average");
      o_ld=ffopen("perfect.ld","w");
    }
    fprintf(out, "%10d  %10g\n",step++,unb*100.0);
    for(i=0; (i<nsb->nnodes); i++)
      fprintf(o_ld,"%3d  ",nload[i+1]-nload[i]);
    fprintf(o_ld,"\n");
  }
  for(i=0; (i<nsb->nnodes); i++)
    nsb->workload[i]=nload[i+1];
  
  /*calc_nsbshift(nsb);*/
  nsb->shift=nsb->nnodes-1;
  nsb->bshift=0;
  
  sfree(nb);
}

