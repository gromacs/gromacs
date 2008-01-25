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

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "macros.h"
#include "main.h"
#include "ns.h"
#include "partdec.h"
#include "splitter.h"

typedef struct gmx_partdec {
  int  neighbor[2];             /* The nodeids of left and right neighb */
  int  *cgindex;                /* The charge group boundaries,         */
				/* size nnodes+1,                       */
                                /* only allocated with particle decomp. */
  int  *index;                  /* The home particle boundaries,        */
				/* size nnodes+1,                       */
                                /* only allocated with particle decomp. */
  int  shift,bshift;		/* Coordinates are shifted left for     */
                                /* 'shift' systolic pulses, and right   */
				/* for 'bshift' pulses. Forces are      */
				/* shifted right for 'shift' pulses     */
				/* and left for 'bshift' pulses         */
				/* This way is not necessary to shift   */
				/* the coordinates over the entire ring */
} gmx_partdec_t;

#ifdef GMX_MPI
static MPI_Request mpi_req_tx=MPI_REQUEST_NULL,mpi_req_rx;
#endif

void gmx_tx(const t_commrec *cr,int dir,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx"); 
#else
  int        nodeid;
  int        tag,flag;
  MPI_Status status;

  nodeid = cr->pd->neighbor[dir];
  
#ifdef DEBUG
  fprintf(stderr,"gmx_tx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
#ifdef MPI_TEST
  /* workaround for crashes encountered with MPI on IRIX 6.5 */
  if (mpi_req_tx != MPI_REQUEST_NULL) {
    MPI_Test(&mpi_req_tx,&flag,&status);
    if (flag==FALSE) {
      fprintf(stdlog,"gmx_tx called before previous send was complete: nodeid=%d, buf=%x, bufsize=%d\n",
	      nodeid,buf,bufsize);
      gmx_tx_wait(nodeid);
    }
  }
#endif
  tag = 0;
  if (MPI_Isend(buf,bufsize,MPI_BYTE,RANK(cr,nodeid),tag,cr->mpi_comm_mygroup,&mpi_req_tx) != 0)
    gmx_comm("MPI_Isend Failed");
#endif
}

void gmx_tx_wait(int dir)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx_wait");
#else
  MPI_Status  status;
  int mpi_result;

  if ((mpi_result=MPI_Wait(&mpi_req_tx,&status)) != 0)
    gmx_fatal(FARGS,"MPI_Wait: result=%d",mpi_result);
#endif
}

void gmx_rx(const t_commrec *cr,int dir,void *buf,int bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_rx");
#else
  int nodeid;
  int tag;

  nodeid = cr->pd->neighbor[dir];
#ifdef DEBUG
  fprintf(stderr,"gmx_rx: nodeid=%d, buf=%x, bufsize=%d\n",
	  nodeid,buf,bufsize);
#endif
  tag = 0;
  if (MPI_Irecv( buf, bufsize, MPI_BYTE, RANK(cr,nodeid), tag, cr->mpi_comm_mygroup, &mpi_req_rx) != 0 )
    gmx_comm("MPI_Recv Failed");
#endif
}

void gmx_rx_wait(int nodeid)
{
#ifndef GMX_MPI
  gmx_call("gmx_rx_wait");
#else
  MPI_Status  status;
  int mpi_result;
  
  if ((mpi_result=MPI_Wait(&mpi_req_rx,&status)) != 0)
    gmx_fatal(FARGS,"MPI_Wait: result=%d",mpi_result);
#endif
}

void gmx_tx_rx_real(const t_commrec *cr,
		    int send_dir,real *send_buf,int send_bufsize,
		    int recv_dir,real *recv_buf,int recv_bufsize)
{
#ifndef GMX_MPI
  gmx_call("gmx_tx_rx_real");
#else
  int send_nodeid,recv_nodeid;
  int tx_tag = 0,rx_tag = 0;
  MPI_Status stat;

  send_nodeid = cr->pd->neighbor[send_dir];
  recv_nodeid = cr->pd->neighbor[recv_dir];

#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif
  
  MPI_Sendrecv(send_buf,send_bufsize,mpi_type,RANK(cr,send_nodeid),tx_tag,
	       recv_buf,recv_bufsize,mpi_type,RANK(cr,recv_nodeid),rx_tag,
	       cr->mpi_comm_mygroup,&stat);
#undef mpi_type
#endif
}
		 
void gmx_wait(int dir_send,int dir_recv)
{
#ifndef GMX_MPI
  gmx_call("gmx_wait");
#else
  gmx_tx_wait(dir_send);
  gmx_rx_wait(dir_recv);
#endif
}

static void set_left_right(t_commrec *cr)
{
  cr->pd->neighbor[GMX_LEFT]  = (cr->nnodes + cr->nodeid - 1) % cr->nnodes;
  cr->pd->neighbor[GMX_RIGHT] = (cr->nodeid + 1) % cr->nnodes;
}

int *pd_cgindex(const t_commrec *cr)
{
  return cr->pd->cgindex;
}

int *pd_index(const t_commrec *cr)
{
  return cr->pd->index;
}

int pd_shift(const t_commrec *cr)
{
  return cr->pd->shift;
}

int pd_bshift(const t_commrec *cr)
{
  return cr->pd->bshift;
}

void pd_cg_range(const t_commrec *cr,int *cg0,int *cg1)
{
  *cg0 = cr->pd->cgindex[cr->nodeid];
  *cg1 = cr->pd->cgindex[cr->nodeid+1];
}

void pd_at_range(const t_commrec *cr,int *at0,int *at1)
{
  *at0 = cr->pd->index[cr->nodeid];
  *at1 = cr->pd->index[cr->nodeid+1];
}

static int home_cpu(int nnodes,gmx_partdec_t *pd,int atomid)
{
  int k;
 
  for(k=0; (k<nnodes); k++) {
    if (atomid<pd->index[k+1])
      return k;
  }
  gmx_fatal(FARGS,"Atomid %d is larger than number of atoms (%d)",
	    atomid+1,pd->index[nnodes]+1);
	    
  return -1;
}

static void calc_nsbshift(FILE *fp,int nnodes,gmx_partdec_t *pd,t_idef *idef)
{
  int i,j,k;
  int lastcg,targetcg,nshift,naaj;
  int homeid[32];
      
  pd->bshift=0;
  for(i=1; (i<nnodes); i++) {
    targetcg = pd->cgindex[i];
    for(nshift=i; (nshift > 0) && (pd->cgindex[nshift] > targetcg); nshift--)
      ;
    pd->bshift=max(pd->bshift,i-nshift);
  }

  pd->shift = (nnodes + 1)/2;
  for(i=0; (i<nnodes); i++) {
    lastcg = pd->cgindex[i+1]-1;
    naaj = calc_naaj(lastcg,pd->cgindex[nnodes]);
    targetcg = (lastcg+naaj) % pd->cgindex[nnodes];
    
    /* Search until we find the target charge group */
    for(nshift=0;
	(nshift < nnodes) && (targetcg > pd->cgindex[nshift+1]);
	nshift++)
      ;
    /* Now compute the shift, that is the difference in node index */
    nshift=((nshift - i + nnodes) % nnodes);
    
    if (fp)
      fprintf(fp,"CPU=%3d, lastcg=%5d, targetcg=%5d, myshift=%5d\n",
	      i,lastcg,targetcg,nshift);
	    
    /* It's the largest shift that matters */
    pd->shift=max(nshift,pd->shift);
  }
  /* Now for the bonded forces */
  for(i=0; (i<F_NRE); i++) {
    if (interaction_function[i].flags & IF_BOND) {
      int nratoms = interaction_function[i].nratoms;
      for(j=0; (j<idef->il[i].nr); j+=nratoms+1) {
	for(k=1; (k<=nratoms); k++)
	  homeid[k-1] = home_cpu(nnodes,pd,idef->il[i].iatoms[j+k]);
	for(k=1; (k<nratoms); k++)
	  pd->shift = max(pd->shift,abs(homeid[k]-homeid[0]));
      }
    }
  }
  if (fp)
    fprintf(fp,"pd->shift = %3d, pd->bshift=%3d\n",
	    pd->shift,pd->bshift);
}

static void init_partdec(FILE *fp,t_commrec *cr,t_block *cgs,int *multinr,
			 t_idef *idef)
{
  int  i;
  gmx_partdec_t *pd;

  snew(pd,1);
  cr->pd = pd;

  set_left_right(cr);
  
  if (cr->nnodes > 1) {
    if (multinr == NULL)
      gmx_fatal(FARGS,"Internal error in init_partdec: multinr = NULL");
    snew(pd->index,cr->nnodes+1);
    snew(pd->cgindex,cr->nnodes+1);
    pd->cgindex[0] = 0;
    pd->index[0]   = 0;
    for(i=0; (i < cr->nnodes); i++) {
      pd->cgindex[i+1] = multinr[i];
      pd->index[i+1]   = cgs->index[multinr[i]];
    }
    calc_nsbshift(fp,cr->nnodes,pd,idef);
    /* This is a hack to do with bugzilla 148 */
    /*pd->shift = cr->nnodes-1;
      pd->bshift = 0;*/
  }
}

static void print_partdec(FILE *fp,char *title,int nnodes,gmx_partdec_t *pd)
{
  int i;

  fprintf(fp,"%s\n",title);
  fprintf(fp,"nnodes:       %5d\n",nnodes);
  fprintf(fp,"pd->shift:    %5d\n",pd->shift);
  fprintf(fp,"pd->bshift:   %5d\n",pd->bshift);
  
  fprintf(fp,"Nodeid   atom0   #atom     cg0       #cg\n");
  for(i=0; (i<nnodes); i++)
    fprintf(fp,"%6d%8d%8d%8d%10d\n",
	    i,
	    pd->index[i],pd->index[i+1]-pd->index[i],
	    pd->cgindex[i],pd->cgindex[i+1]-pd->cgindex[i]);
  fprintf(fp,"\n");
}

static void pr_idef_division(FILE *fp,t_idef *idef,int nnodes,int **multinr)
{
  int i,ftype,nr,nra,m0,m1;

  fprintf(fp,"Division of bonded forces over processors\n");
  fprintf(fp,"%-12s","CPU");
  for(i=0; (i<nnodes); i++) 
    fprintf(fp," %5d",i);
  fprintf(fp,"\n");
  
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if (idef->il[ftype].nr > 0) {
      nr  = idef->il[ftype].nr;
      nra = 1+interaction_function[ftype].nratoms;
      fprintf(fp,"%-12s", interaction_function[ftype].name);
      /* Loop over processors */
      for(i=0; (i<nnodes); i++) {
	m0 = (i == 0) ? 0 : multinr[ftype][i-1]/nra;
	m1 = multinr[ftype][i]/nra;
	fprintf(fp," %5d",m1-m0);
      }
      fprintf(fp,"\n");
    }
  }
}

static void select_my_ilist(FILE *log,t_ilist *il,int *multinr,t_commrec *cr)
{
  t_iatom *ia;
  int     i,start,end,nr;
  
  if (cr->nodeid == 0)
    start=0;
  else
    start=multinr[cr->nodeid-1];
  end=multinr[cr->nodeid];
  
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
  
  il->nr=nr;
}

static void select_my_idef(FILE *log,t_idef *idef,int **multinr,t_commrec *cr)
{
  int i;
  
  for(i=0; (i<F_NRE); i++)
    select_my_ilist(log,&idef->il[i],multinr[i],cr);
}

void split_system(FILE *log,t_inputrec *inputrec,t_state *state,
		  t_commrec *cr,t_topology *top)
{
  int    i,npp,n,nn;
  real   *capacity;
  double tcap = 0,cap;
  int    *multinr_cgs,**multinr_nre;
  char   *cap_env;

  /* Time to setup the division of charge groups over processors */
  npp = cr->nnodes-cr->npmenodes;
  snew(capacity,npp);
  cap_env = getenv("GMX_CAPACITY");
  if (cap_env == NULL) {
    for(i=0; (i<npp-1); i++) {
      capacity[i] = 1.0/(double)npp;
      tcap += capacity[i];
    }
    /* Take care that the sum of capacities is 1.0 */
    capacity[npp-1] = 1.0 - tcap;
  } else {
    tcap = 0;
    nn = 0;
    for(i=0; i<npp; i++) {
      cap = 0;
      sscanf(cap_env+nn,"%lf%n",&cap,&n);
      if (cap == 0)
	gmx_fatal(FARGS,"Incorrect entry or number of entries in GMX_CAPACITY='%s'",cap_env);
      capacity[i] = cap;
      tcap += cap;
      nn += n;
    }
    for(i=0; i<npp; i++)
      capacity[i] /= tcap;
  }

  snew(multinr_cgs,npp);
  snew(multinr_nre,F_NRE);
  for(i=0; i<F_NRE; i++)
    snew(multinr_nre[i],npp);
  
  /* This computes which entities can be placed on processors */
  split_top(log,npp,top,capacity,multinr_cgs,multinr_nre);
  sfree(capacity);
  init_partdec(log,cr,&(top->blocks[ebCGS]),multinr_cgs,&(top->idef));

  /* This should be fine */
  /*split_idef(&(top->idef),cr->nnodes-cr->npmenodes);*/
  
  select_my_idef(log,&(top->idef),multinr_nre,cr);
  
  if (log)
    pr_idef_division(log,&(top->idef),npp,multinr_nre);

  for(i=0; i<F_NRE; i++)
    sfree(multinr_nre[i]);
  sfree(multinr_nre);
  sfree(multinr_cgs);
  
  if (log)
    print_partdec(log,"Workload division",cr->nnodes,cr->pd);
}

static void create_vsitelist(int nindex, int *list,
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
  

bool setup_parallel_vsites(t_idef *idef,t_commrec *cr,
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
  gmx_partdec_t *pd;

#define BUFLEN 100

  pd = cr->pd;

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
	else if(ftype==F_VSITE4FDN || ftype==F_VSITE4FD)
	  nconstr=4;
	else
	  nconstr=3;
	
	minidx=ia[1];
	for(j=2;j<nconstr+2;j++) 
	  if(ia[j]<minidx)
	    minidx=ia[j];

	minhome=0;
	while(minidx>=(pd->index[minhome+1]))
          minhome++;

	if(minhome==cr->nodeid) {
	  /* This is my vsite interaction - but is the vsite local?
	   * If not, he must be on the next node (error otherwise)
	   * (but we do account for the cyclic ring structure)
	   */
	  if(ia[1]<pd->index[cr->nodeid] ||
	     ia[1]>=(pd->index[cr->nodeid+1])) {
	    if((nnextvsite%BUFLEN)==0 && nnextvsite>0)
	      srenew(idxnextvsite,nnextvsite+BUFLEN);
	    idxnextvsite[nnextvsite++]=ia[1];
	    found=TRUE;
	  }
	  for(j=2;j<nconstr+2;j++) {
	    inr=ia[j];
	    ihome=0;
	    while(inr>=(pd->index[ihome+1]))
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
	  if(ia[1]>=pd->index[cr->nodeid] &&
	     (ia[1]<(pd->index[cr->nodeid+1]))) {
	    if((nprevvsite%BUFLEN)==0 && nprevvsite>0)
	      srenew(idxprevvsite,nprevvsite+BUFLEN);
	    idxprevvsite[nprevvsite++]=ia[1];
	    found=TRUE;
	  }
	  for(j=2;j<nconstr+2;j++) {
	    inr=ia[j];
	    if(ia[j]>=pd->index[cr->nodeid] &&
	       (ia[1]<(pd->index[cr->nodeid+1]))) {
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
