/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Good gRace! Old Maple Actually Chews Slate
 */
#include "sysstuff.h"
#include "assert.h"
#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "macros.h"
#include "main.h"
#include "nsb.h"
#include "ns.h"

void calc_nsbshift(t_nsborder *nsb)
{
  int i;
  int firstcg,lastcg,targetcg,nshift,naaj;
  
  nsb->bshift=0;
  for(i=1; (i<nsb->nprocs); i++) {
    targetcg = nsb->workload[i-1];
    for(nshift=i; (nshift > 0) && (nsb->cgload[nshift-1] > targetcg); nshift--)
      ;
    nsb->bshift=max(nsb->bshift,i-nshift);
  }

  nsb->shift=(nsb->nprocs+1)/2;
  for(i=0; (i<nsb->nprocs); i++) {
    lastcg=nsb->cgload[i]-1;
    naaj=calc_naaj(stdlog,lastcg,nsb->cgtotal);
    targetcg=(lastcg+naaj) % nsb->cgtotal;
    
    /* Search until we find the target charge group */
    for(nshift=0; (nshift < nsb->nprocs) && (targetcg > nsb->cgload[nshift]);
	nshift++)
      ;
    /* Now compute the shift, that is the difference in processor index */
    nshift=((nshift-i+nsb->nprocs) % nsb->nprocs);
    
    fprintf(stdlog,"CPU=%3d, lastcg=%5d, targetcg=%5d, myshift=%5d\n",
	    i,lastcg,targetcg,nshift);
	    
    /* It's the largest shift that matters */
    nsb->shift=max(nshift,nsb->shift);
  }
  fprintf(stdlog,"nsb->shift = %3d, nsb->bshift=%3d\n",
	  nsb->shift,nsb->bshift);
}

void calc_nsb(t_block *cgs,int nprocs,t_nsborder *nsb,int nstDlb)
{
  int  i,j,k,cg0,mincg,maxcg,npid;
  bool bDone;
  
  /* Clean! */
  for(i=0; (i<MAXPROC); i++) 
    nsb->homenr[i]=nsb->index[i]=nsb->cgload[i]=nsb->workload[i]=0;
  
  nsb->nprocs=nprocs;
  nsb->nstDlb=nstDlb;
  nsb->cgtotal=cgs->nr;
  nsb->natoms=cgs->nra;
  for(i=0; (i<nprocs); i++) {
    cg0              = (i > 0) ? cgs->multinr[i-1] : 0;
    nsb->cgload[i]   = cgs->multinr[i];
    nsb->workload[i] = cgs->multinr[i];
    nsb->index[i]    = cgs->index[cg0];
    nsb->homenr[i]   = cgs->index[cgs->multinr[i]]-nsb->index[i];
  }
  calc_nsbshift(nsb);
}

void print_nsb(FILE *fp,char *title,t_nsborder *nsb)
{
  int i,tcg,mycg,h;

  fprintf(fp,"%s\n",title);
  fprintf(fp,"nsb->pid:     %5d\n",nsb->pid);
  fprintf(fp,"nsb->nprocs:  %5d\n",nsb->nprocs);
  fprintf(fp,"nsb->cgtotal: %5d\n",nsb->cgtotal);
  fprintf(fp,"nsb->natoms:  %5d\n",nsb->natoms);
  fprintf(fp,"nsb->shift:   %5d\n",nsb->shift);
  fprintf(fp,"nsb->bshift:  %5d\n",nsb->bshift);
  
  fprintf(fp,"pid   index  homenr  cgload  workload\n");
  for(i=0; (i<nsb->nprocs); i++)
    fprintf(fp,"%3d%8d%8d%8d%10d\n",i,
	    nsb->index[i],nsb->homenr[i],nsb->cgload[i],nsb->workload[i]);
  fprintf(fp,"\n");
}

