/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
  
  real   ntp,ntot,load[MAXPROC+1],load_av,unb,tld,ratio,defcost;
  int    nload[MAXPROC+1];
  int    ncg,pid,i,j,k;
  
  if (nsb->nprocs < 2)
    return;
  
  ncg=nsb->cgtotal;
  snew(nb,ncg+1);
   
  nload[0]=0;
  for(i=1; (i<=cr->nprocs); i++)
    nload[i]=nsb->workload[i-1];
    
  /* Initiate the cost per charge group,
   * there are nns particles for this proc, 
   */
  defcost=0.1*(nns*cost_nrnb(eNR_NS)+nlr*cost_nrnb(eNR_LR))/
    (nload[cr->pid+1]-nload[cr->pid]);
  for(i=0; (i<nload[cr->pid]); i++)
    nb[i]=0;
  for(; (i<nload[cr->pid+1]); i++)
    nb[i]=defcost;
  for(; (i<ncg+1); i++)
    nb[i]=0;
  ntot=(nload[cr->pid+1]-nload[cr->pid])*defcost;
      
  if (cgnr==NULL) {
    cgnr=make_invblock(cgs,nsb->natoms);
    
    /* Load from bonded forces is constant... */
    snew(bfload,cr->nprocs+1);
    bfload[cr->pid]=calc_bfload(idef->il[eilBONDEDS].nr,
				idef->il[eilBONDEDS].iatoms,
				idef->functype,idef->iparams);
    fprintf(stdlog,"My    Bonded Force Load: %.0f\n",bfload[cr->pid]);
    bfload[cr->nprocs]=bfload[cr->pid];
    gmx_sum(cr->nprocs+1,bfload,cr);
    fprintf(stdlog,"Total Bonded Force Load: %.0f\n",bfload[cr->nprocs]);
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
  ntot+=bfload[cr->nprocs];
  fprintf(stdlog,"Total           Force Load =%.0f\n",ntot);
  
  load_av=ntot/nsb->nprocs;
  unb=1;
  for(pid=0; (pid<nsb->nprocs); pid++) {
    load[pid]=bfload[pid];
    for(j=nload[pid]; (j<nload[pid+1]); j++) 
      load[pid]+=nb[j];
#ifdef DEBUG
    fprintf(stdlog,"pid: %2d, cg0: %5d, hcgHi: %3d, load: %10g\n",
	    pid,nsb->cg0[pid],nsb->hcgHi[pid],load[pid]);
#endif
    ratio=load[pid]/load_av;
    unb=max(unb,ratio);
  }
  tld=0;
  pid=1;
  nload[0]=0;
  for(i=0; ((i<ncg) && (pid < nsb->nprocs)); i++) {
    tld+=nb[i];
    if (tld >= (pid*load_av)) {
      nload[pid]=i;
      pid++;
    }
  }
  if (pid == nsb->nprocs)
    nload[pid]=ncg;
  else {
    fprintf(stdlog,"pid(%d) != nprocs (%d)\n",pid,nsb->nprocs);
  }
  if (nsb->pid == 0) {
    if (out==NULL) {
      out=xvgropen("load.xvg","Maximum Load","Step","% of average");
      o_ld=ffopen("perfect.ld","w");
    }
    fprintf(out, "%10d  %10g\n",step++,unb*100.0);
    for(i=0; (i<nsb->nprocs); i++)
      fprintf(o_ld,"%3d  ",nload[i+1]-nload[i]);
    fprintf(o_ld,"\n");
  }
  for(i=0; (i<nsb->nprocs); i++)
    nsb->workload[i]=nload[i+1];
  
  /*calc_nsbshift(nsb);*/
  nsb->shift=nsb->nprocs-1;
  nsb->bshift=0;
  
  sfree(nb);
}

