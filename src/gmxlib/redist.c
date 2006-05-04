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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "tpxio.h"
#include "smalloc.h"
#include "vec.h"
#include "main.h"
#include "mvdata.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "txtdump.h"
#include "splitter.h"
#include "mdatoms.h"
#include "mdrun.h"
#include "statutil.h"
#include "names.h"

static void pr_idef_division(FILE *fp,t_idef *idef,int nnodes)
{
  int i,ftype,nr,nra,m0,m1;
  
  fprintf(stdlog,"Division of bonded forces over processors\n");
  fprintf(stdlog,"%-12s","CPU");
  for(i=0; (i<nnodes); i++) 
    fprintf(stdlog," %5d",i);
  fprintf(stdlog,"\n");
  
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (idef->il[ftype].nr > 0) {
      nr  = idef->il[ftype].nr;
      nra = 1+interaction_function[ftype].nratoms;
      fprintf(stdlog,"%-12s", interaction_function[ftype].name);
      /* Loop over processors */
      for(i=0; (i<nnodes); i++) {
	m0 = (i == 0) ? 0 : idef->il[ftype].multinr[i-1]/nra;
	m1 = idef->il[ftype].multinr[i]/nra;
	fprintf(stdlog," %5d",m1-m0);
      }
      fprintf(stdlog,"\n");
    }
  }
}

static void split_idef_old(t_idef *idef,int nnodes)
{
  int i,ftype,nr,nra,m1;
  
  /* Loop over energy types */
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (idef->il[ftype].nr > 0) {
      nr  = idef->il[ftype].nr;
      nra = 1+interaction_function[ftype].nratoms;
      /* Loop over processors */
      for(i=0; (i<nnodes-1); i++) {
	m1 = ((i+1)*(nr/nra)/nnodes);
	idef->il[ftype].multinr[i] = nra*m1;
      }
      idef->il[ftype].multinr[i]=nr;
    }
  }
}

static void select_my_ilist(FILE *log,t_ilist *il,t_commrec *cr)
{
  t_iatom *ia;
  int     i,start,end,nr;
  
  if (cr->nodeid == 0)
    start=0;
  else
    start=il->multinr[cr->nodeid-1];
  end=il->multinr[cr->nodeid];
  
  nr=end-start;
  if (nr < 0)
    gmx_fatal(FARGS,"Negative number of atoms (%d) on node %d\n"
		"You have probably not used the same value for -np with grompp"
		" and mdrun",
		nr,cr->nodeid);
  snew(ia,nr);

  for(i=0; (i<nr); i++)
    ia[i]=il->iatoms[start+i];

  sfree(il->iatoms);
  il->iatoms=ia;
  
  for(i=0; (i<MAXNODES); i++)
    il->multinr[i]=nr;
  il->nr=nr;
}

static void select_my_idef(FILE *log,t_idef *idef,t_commrec *cr)
{
  int i;
  
  for(i=0; (i<F_NRE); i++)
    select_my_ilist(log,&idef->il[i],cr);
}
	


void split_system_first(FILE *log,t_inputrec *inputrec,t_state *state,
			t_commrec *cr,t_topology *top,t_nsborder *nsb)
{
  int    i,npp;
  real   *capacity;
  double tcap = 0;
  
  /* Time to setup the division of charge groups over processors */
  npp = cr->nnodes-cr->npmenodes;
  snew(capacity,npp);
  for(i=0; (i<npp-1); i++) {
    capacity[i] = 1.0/(double)npp;
    tcap += capacity[i];
  }
  /* Take care that the sum of capacities is 1.0 */
  capacity[npp-1] = 1.0 - tcap;

  /* This computes which entities can be placed on processors */
  split_top(log,npp,top,capacity);
  sfree(capacity);
  calc_nsb(log,&(top->blocks[ebCGS]),cr->nnodes,nsb,0);

  /* This should be fine */
  /*split_idef(&(top->idef),cr->nnodes-cr->npmenodes);*/
  
  select_my_idef(log,&(top->idef),cr);
  
  pr_idef_division(log,&(top->idef),npp);
  
  print_nsb(log,"Workload division",nsb);
}

void split_system_again()
{
}

void create_vsitelist(int nindex, int *list,
		      int *targetn, int **listptr)
{
  int i,j,k,inr;
  int minidx;
  int *newlist;

  /* remove duplicates */
  for(i=0;i<nindex;i++) {
    inr=list[i];
    for(j=i+1;j<nindex;j++) {
      if(list[j]==inr) {
	for(k=j;k<nindex-1;k++)
	  list[k]=list[k+1];
	nindex--;
      }
    }
  }

  *targetn=nindex;
  snew(newlist,nindex);
  
  /* sort into the new array */
  for(i=0;i<nindex;i++) {
    inr=-1;
    for(j=0;j<nindex;j++)
      if(list[j]>0 && (inr==-1 || list[j]<list[inr])) 
	inr=j; /* smallest so far */
    newlist[i]=list[inr];
    list[inr]=-1;
  }
  *listptr=newlist;
}
  

bool setup_parallel_vsites(t_idef *idef,t_commrec *cr,t_nsborder *nsb,
			    t_comm_vsites *vsitecomm)
{
  int i,inr,j,k,ftype;
  int minidx,minhome,ihome;
  int nra,nrd,nconstr;
  bool found=FALSE;
  t_iatom   *ia;
  int *idxprevvsite;
  int *idxnextvsite;
  int *idxprevconstr;
  int *idxnextconstr;
  int  nprevvsite=0,nnextvsite=0;
  int  nprevconstr=0,nnextconstr=0;

#define BUFLEN 100
  
  snew(idxprevvsite,BUFLEN);
  snew(idxnextvsite,BUFLEN);
  snew(idxprevconstr,BUFLEN);
  snew(idxnextconstr,BUFLEN);  

  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (interaction_function[ftype].flags & IF_VSITE) {
      nra    = interaction_function[ftype].nratoms;
      nrd    = idef->il[ftype].nr;
      ia     = idef->il[ftype].iatoms;
      
      for(i=0; (i<nrd); ) {
	
	/* The vsite and constructing atoms */
	if (ftype==F_VSITE2)
	  nconstr=2;
	else if(ftype==F_VSITE4FD)
	  nconstr=4;
	else
	  nconstr=3;
	
	minidx=ia[1];
	for(j=2;j<nconstr+2;j++) 
	  if(ia[j]<minidx)
	    minidx=ia[j];

	minhome=0;
	while(minidx>=(nsb->index[minhome]+nsb->homenr[minhome]))
          minhome++;

	if(minhome==cr->nodeid) {
	  /* This is my vsite interaction - but is the vsite local?
	   * If not, he must be on the next node (error otherwise)
	   * (but we do account for the cyclic ring structure)
	   */
	  if(ia[1]<nsb->index[cr->nodeid] ||
	     ia[1]>=(nsb->index[cr->nodeid]+nsb->homenr[cr->nodeid])) {
	    if((nnextvsite%BUFLEN)==0 && nnextvsite>0)
	      srenew(idxnextvsite,nnextvsite+BUFLEN);
	    idxnextvsite[nnextvsite++]=ia[1];
	    found=TRUE;
	  }
	  for(j=2;j<nconstr+2;j++) {
	    inr=ia[j];
	    ihome=0;
	    while(inr>=(nsb->index[ihome]+nsb->homenr[ihome]))
	      ihome++;
	    if( ihome>(cr->nodeid+1))
	      gmx_fatal(FARGS,"Vsite particle %d and its constructing"
			  " atoms are not on the same or adjacent\n" 
			  " nodes. This is necessary to avoid a lot\n"
			  " of extra communication. The easiest way"
			  " to ensure this is to place vsites\n"
			  " close to the constructing atoms.\n"
			  " Sorry, but you will have to rework your topology!\n",
			  ia[1]);
	    else if(ihome==((cr->nodeid+1)%cr->nnodes)) {
	      if((nnextconstr%BUFLEN)==0 && nnextconstr>0)
		srenew(idxnextconstr,nnextconstr+BUFLEN);
	      idxnextconstr[nnextconstr++]=ia[j];
	      found=TRUE;
	    }
	  }
	} else if(minhome==((cr->nodeid-1+cr->nnodes)%cr->nnodes)) {
	  /* Not our vsite, but we might be involved */
	  if(ia[1]>=nsb->index[cr->nodeid] &&
	     (ia[1]<(nsb->index[cr->nodeid]+nsb->homenr[cr->nodeid]))) {
	    if((nprevvsite%BUFLEN)==0 && nprevvsite>0)
	      srenew(idxprevvsite,nprevvsite+BUFLEN);
	    idxprevvsite[nprevvsite++]=ia[1];
	    found=TRUE;
	  }
	  for(j=2;j<nconstr+2;j++) {
	    inr=ia[j];
	    if(ia[j]>=nsb->index[cr->nodeid] &&
	       (ia[1]<(nsb->index[cr->nodeid]+nsb->homenr[cr->nodeid]))) {
	      if((nprevconstr%BUFLEN)==0 && nprevconstr>0)
		srenew(idxprevconstr,nprevconstr+BUFLEN);
	      idxprevconstr[nprevconstr++]=ia[j];
	      found=TRUE;
	    }
	  }
	}
	/* Increment loop variables */
	i  += nra+1;
	ia += nra+1;
      }
    }
  }

  create_vsitelist(nprevvsite,idxprevvsite,
		   &(vsitecomm->nprevvsite),&(vsitecomm->idxprevvsite));
  create_vsitelist(nnextvsite,idxnextvsite,
		   &(vsitecomm->nnextvsite),&(vsitecomm->idxnextvsite));
  create_vsitelist(nprevconstr,idxprevconstr,
		   &(vsitecomm->nprevconstr),&(vsitecomm->idxprevconstr));
  create_vsitelist(nnextconstr,idxnextconstr,
		   &(vsitecomm->nnextconstr),&(vsitecomm->idxnextconstr));

  sfree(idxprevvsite);
  sfree(idxnextvsite);
  sfree(idxprevconstr);
  sfree(idxnextconstr);

  return found;
#undef BUFLEN
}

