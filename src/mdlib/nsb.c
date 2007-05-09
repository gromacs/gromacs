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

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "macros.h"
#include "main.h"
#include "nsb.h"
#include "ns.h"

static int home_cpu(int nnodes,t_nsborder *nsb,int atomid)
{
  int k;
 
  for(k=0; (k<nnodes); k++) {
    if (atomid<nsb->index[k]+nsb->homenr[k])
      return k;
  }
  gmx_fatal(FARGS,"Atomid %d is larger than number of atoms (%d)",
	    atomid+1,nsb->index[nnodes-1]+nsb->homenr[nnodes-1]+1);
	    
  return -1;
}

void calc_nsbshift(FILE *fp,t_nsborder *nsb,t_idef *idef)
{
  int i,j,k;
  int lastcg,targetcg,nshift,naaj;
  int homeid[32];

  nsb->bshift=0;
  for(i=1; (i<nsb->nnodes); i++) {
    targetcg = nsb->workload[i-1];
    for(nshift=i; (nshift > 0) && (nsb->cgload[nshift-1] > targetcg); nshift--)
      ;
    nsb->bshift=max(nsb->bshift,i-nshift);
  }

  nsb->shift=(nsb->nnodes+1)/2;
  for(i=0; (i<nsb->nnodes); i++) {
    lastcg=nsb->cgload[i]-1;
    naaj=calc_naaj(lastcg,nsb->cgtotal);
    targetcg=(lastcg+naaj) % nsb->cgtotal;
    
    /* Search until we find the target charge group */
    for(nshift=0; (nshift < nsb->nnodes) && (targetcg > nsb->cgload[nshift]);
	nshift++)
      ;
    /* Now compute the shift, that is the difference in node index */
    nshift=((nshift-i+nsb->nnodes) % nsb->nnodes);
    
    if (fp)
      fprintf(fp,"CPU=%3d, lastcg=%5d, targetcg=%5d, myshift=%5d\n",
	      i,lastcg,targetcg,nshift);
	    
    /* It's the largest shift that matters */
    nsb->shift=max(nshift,nsb->shift);
  }
  /* Now for the bonded forces */
  for(i=0; (i<F_NRE); i++) {
    if (interaction_function[i].flags & IF_BOND) {
      int nratoms = interaction_function[i].nratoms;
      for(j=0; (j<idef->il[i].nr); j+=nratoms+1) {
	for(k=1; (k<=nratoms); k++)
	  homeid[k-1] = home_cpu(nsb->nnodes,nsb,idef->il[i].iatoms[j+k]);
	for(k=1; (k<nratoms); k++)
	  nsb->shift = max(nsb->shift,abs(homeid[k]-homeid[0]));
      }
    }
  }
  if (fp)
    fprintf(fp,"nsb->shift = %3d, nsb->bshift=%3d\n",
	    nsb->shift,nsb->bshift);
}

void calc_nsb(FILE *fp,t_block *cgs,int nnodes,t_nsborder *nsb,int nstDlb,
	      t_idef *idef)
{
  int  i,cg0;
  
  /* Clean! */
  for(i=0; (i<MAXNODES); i++) 
    nsb->homenr[i]=nsb->index[i]=nsb->cgload[i]=nsb->workload[i]=0;
  
  nsb->nnodes=nnodes;
  nsb->nstDlb=nstDlb;
  nsb->cgtotal=cgs->nr;
  nsb->natoms=cgs->nra;
  for(i=0; (i<nnodes); i++) {
    cg0              = (i > 0) ? cgs->multinr[i-1] : 0;
    nsb->cgload[i]   = cgs->multinr[i];
    nsb->workload[i] = cgs->multinr[i];
    nsb->index[i]    = cgs->index[cg0];
    nsb->homenr[i]   = cgs->index[cgs->multinr[i]]-nsb->index[i];
  }
  calc_nsbshift(fp,nsb,idef);
}

void print_nsb(FILE *fp,char *title,t_nsborder *nsb)
{
  int i;

  fprintf(fp,"%s\n",title);
  fprintf(fp,"nsb->nodeid:     %5d\n",nsb->nodeid);
  fprintf(fp,"nsb->nnodes:  %5d\n",nsb->nnodes);
  fprintf(fp,"nsb->cgtotal: %5d\n",nsb->cgtotal);
  fprintf(fp,"nsb->natoms:  %5d\n",nsb->natoms);
  fprintf(fp,"nsb->shift:   %5d\n",nsb->shift);
  fprintf(fp,"nsb->bshift:  %5d\n",nsb->bshift);
  
  fprintf(fp,"Nodeid   index  homenr  cgload  workload\n");
  for(i=0; (i<nsb->nnodes); i++)
    fprintf(fp,"%6d%8d%8d%8d%10d\n",i,
	    nsb->index[i],nsb->homenr[i],nsb->cgload[i],nsb->workload[i]);
  fprintf(fp,"\n");
}

