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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_splitter_c = "$Id$";

#include <stdio.h>
#include "sysstuff.h"
#include "assert.h"
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "mshift.h"
#include "invblock.h"

typedef struct {
  int     nr;
  t_iatom *ia;
} t_sf;

static void init_sf(int nr,t_sf sf[])
{
  int i;

  for(i=0; (i<nr); i++) {
    sf[i].nr=0;
    sf[i].ia=NULL;
  }
}

static void done_sf(int nr,t_sf sf[])
{
  int i;

  for(i=0; (i<nr); i++) {
    sf[i].nr=0;
    sfree(sf[i].ia);
    sf[i].ia=NULL;
  }
}

static void push_sf(t_sf *sf,int nr,t_iatom ia[])
{
  int i;

  srenew(sf->ia,sf->nr+nr);
  for(i=0; (i<nr); i++)
    sf->ia[sf->nr+i]=ia[i];
  sf->nr+=nr;
}

static int min_pid(int nr,atom_id list[],int hid[])
{
  int i,pid,minpid;

  assert(nr > 0);
  minpid=hid[list[0]];
  for (i=1; (i<nr); i++) 
    if ((pid=hid[list[i]]) < minpid) 
      minpid=pid;

  return minpid;
}


static void split_force2(int nprocs,int hid[],t_idef *idef,t_ilist *ilist)
{
  int     i,type,ftype,pid,nratoms,tnr;
  t_iatom ai;
  t_sf    sf[MAXPROC];
  
  init_sf(MAXPROC,sf);

  /* Walk along all the bonded forces, find the appropriate processor
   * to calc it on, and add it to that processors list.
   */
  for (i=0; i<ilist->nr; i+=(1+nratoms)) {
    type    = ilist->iatoms[i];
    ftype   = idef->functype[type];
    nratoms = interaction_function[ftype].nratoms;

    if (ftype == F_SHAKE) {
      /* SPECIAL CASE: All Atoms must have the same home processor! */
      pid=hid[ilist->iatoms[i+1]];
      if (hid[ilist->iatoms[i+2]] != pid) 
	fatal_error(0,"Shake block crossing processor boundaries\n"
		    "constraint between atoms (%d,%d)",
		    ilist->iatoms[i+1],ilist->iatoms[i+2]);
    }
    else if (ftype == F_SETTLE) {
      /* Only the first particle is stored for settles ... */
      ai=ilist->iatoms[i+1];
      pid=hid[ai];
      if ((pid != hid[ai+1]) ||
	  (pid != hid[ai+2]))
	fatal_error(0,"Settle block crossing processor boundaries\n"
		    "constraint between atoms (%d-%d)",ai,ai+2);
    }
    else 
      pid=min_pid(nratoms,&ilist->iatoms[i+1],hid);

    /* Add it to the list */
    push_sf(&(sf[pid]),nratoms+1,&(ilist->iatoms[i]));
  }
  tnr=0;
  for(pid=0; (pid<MAXPROC); pid++) {
    for (i=0; (i<sf[pid].nr); i++) 
      ilist->iatoms[tnr++]=sf[pid].ia[i];

    ilist->multinr[pid]=(pid==0) ? 0 : ilist->multinr[pid-1];
    ilist->multinr[pid]+=sf[pid].nr;
  }
  assert(tnr==ilist->nr);
  done_sf(MAXPROC,sf);
}

static int *home_index(int nprocs,t_block *cgs)
{
  /* This routine determines the processor id for each particle */
  int *hid;
  int pid,j0,j1,j,k,ak;
  
  snew(hid,cgs->nra);
  /* Initiate to -1 to make it possible to check afterwards,
   * all hid's should be set in the loop below
   */
  for(k=0; (k<cgs->nra); k++)
    hid[k]=-1;
    
  /* loop over processors */
  for(pid=0; (pid<nprocs); pid++) {
    j0 = (pid==0) ? 0 : cgs->multinr[pid-1];
    j1 = cgs->multinr[pid];
    
    /* j0 and j1 are the boundariesin the index array */
    for(j=j0; (j<j1); j++) {
      for(k=cgs->index[j]; (k<cgs->index[j+1]); k++) {
	ak=cgs->a[k];
	hid[ak]=pid;
      }
    }
  }
  /* Now verify that all hid's are not -1 */
  for(k=0; (k<cgs->nra); k++)
    assert(hid[k] != -1);
    
  return hid;
}

static int max_index(int start,t_block *b)
{
  int k,k0,k1;
  int mi=-1;

  if (start < b->nr) {  
    k0=b->index[start];
    k1=b->index[start+1];
    if (k1 > k0) {
      mi=b->a[k0];
      for(k=k0+1; (k<k1); k++)
	mi=max(mi,b->a[k]);
    }
  }
  return mi;
}

static int min_index(int start,t_block *b)
{
  int k,k0,k1;
  int mi=INT_MAX;

  if (start < b->nr) {  
    k0=b->index[start];
    k1=b->index[start+1];
    if (k1 > k0) {
      mi=b->a[k0];
      for(k=k0+1; (k<k1); k++)
	mi=min(mi,b->a[k]);
    }
  }
  return mi;
}

static void split_blocks(bool bVerbose,int nprocs,
			 t_block *cgs,t_block *shakes)
{
  int     maxatom[MAXPROC];
  int     i,ai,sbl;
  int     pid;
  real    load,tload;
  
  bool    bSHK;
  atom_id *shknum;
  
  load  = cgs->nra / (real)nprocs;  
  tload = load;
  
  if ((shakes->nr > 0) && (bVerbose)) {
    fprintf(stderr,
	    "Going to use the WINKELHAAK Algorithm to split over %d cpus\n",
	    nprocs);
  }
  shknum = make_invblock(shakes,cgs->nra+1);
  if (debug)
    for(i=0; (i<cgs->nra); i++)
      fprintf(debug,"i: %5d, shknum: %5d\n",i,shknum[i]);
  
  pid = 0;
  sbl = 0;
  for(i=0; (i<cgs->nr); i++) {
    ai   = cgs->a[cgs->index[i]];
    bSHK = ((i == 0) || 
	    ((shknum[ai] == NO_ATID) || (shknum[ai] != shknum[ai-1])));
    
    if (shknum[ai] != NO_ATID) 
      sbl=max(sbl,shknum[ai]);
    
    if (bSHK && (cgs->a[cgs->index[i+1]] >= tload)) {
      cgs->multinr[pid]    = i;
      shakes->multinr[pid] = sbl;
      tload               += load;
      maxatom[pid]         = ai;
      pid++;
    }
  }
  /* Now the last one... */
  while (pid < nprocs) {
    cgs->multinr[pid]=cgs->nr;
    shakes->multinr[pid]=shakes->nr;
    maxatom[pid]=cgs->nra;
    pid++;
  }
  if (pid != nprocs) {
    fatal_error(0,"pid = %d, nprocs = %d, file %s, line %d",
		pid,nprocs,__FILE__,__LINE__);
  }

  if (bVerbose) {
    for(i=nprocs-1; (i>0); i--)
      maxatom[i]-=maxatom[i-1];
    fprintf(stderr,"Division over processors in atoms:\n");
    for(i=0; (i<nprocs); i++)
      fprintf(stderr,"%6d",maxatom[i]);
    fprintf(stderr,"\n");
  }
  sfree(shknum);
}

static void def_mnr(int nr,int mnr[])
{
  int i;

  for (i=0; (i<MAXPROC); i++) 
    mnr[i]=0;
  mnr[0]=nr;
}

void split_top(bool bVerbose,int nprocs,t_topology *top)
{
  int j;
  int *homeind;
  
  if ((bVerbose) && (nprocs>1))
    fprintf(stderr,"splitting topology...\n");
  
  for(j=0; (j<F_NRE); j++)  
    def_mnr(top->idef.il[j].nr,top->idef.il[j].multinr);
  def_mnr(top->atoms.excl.nr,top->atoms.excl.multinr);

  for (j=0; j<ebNR; j++) 
    def_mnr(top->blocks[j].nr,top->blocks[j].multinr);

  if (nprocs > 1) {
    split_blocks(bVerbose,nprocs,
		 &(top->blocks[ebCGS]),&(top->blocks[ebSBLOCKS]));
    homeind=home_index(nprocs,&(top->blocks[ebCGS]));
    for(j=0; (j<F_NRE); j++)
      split_force2(nprocs,homeind,&top->idef,&top->idef.il[j]);
    sfree(homeind);
  }
}

typedef struct {
  int atom,sid;
} t_sid;

static int sid_comp(const void *a,const void *b)
{
  t_sid *sa,*sb;
  int   dd;
  
  sa=(t_sid *)a;
  sb=(t_sid *)b;
  
  dd=sa->sid-sb->sid;
  if (dd == 0)
    return (sa->atom-sb->atom);
  else
    return dd;
}

typedef enum { egcolWhite, egcolGrey, egcolBlack, egcolNR } egCol;

static int mk_grey(int nnodes,egCol egc[],t_graph *g,int *AtomI,
		   t_sid sid[])
{
  int  j,ng,ai,aj,g0;

  ng=0;
  ai=*AtomI;
  
  g0=g->start;
  /* Loop over all the bonds */
  for(j=0; (j<g->nedge[ai]); j++) {
    aj=g->edge[ai][j]-g0;
    /* If there is a white one, make it gray and set pbc */
    if (egc[aj] == egcolWhite) {
      if (aj < *AtomI)
	*AtomI = aj;
      egc[aj] = egcolGrey;

      /* Check whether this one has been set before... */
      if (sid[aj+g0].sid != -1) 
	fatal_error(0,"sid[%d]=%d, sid[%d]=%d, file %s, line %d",
		    ai,sid[ai+g0].sid,aj,sid[aj+g0].sid,__FILE__,__LINE__);
      else {
	sid[aj+g0].sid  = sid[ai+g0].sid;
	sid[aj+g0].atom = aj+g0;
      }
      ng++;
    }
  }
  return ng;
}

static int first_colour(int fC,egCol Col,t_graph *g,egCol egc[])
/* Return the first node with colour Col starting at fC.
 * return -1 if none found.
 */
{
  int i;
  
  for(i=fC; (i<g->nnodes); i++)
    if ((g->nedge[i] > 0) && (egc[i]==Col))
      return i;
  
  return -1;
}

static int mk_sblocks(bool bVerbose,t_graph *g,t_sid sid[])
{
  int    ng,nnodes;
  int    nW,nG,nB;		/* Number of Grey, Black, White	*/
  int    fW,fG;			/* First of each category	*/
  static egCol *egc=NULL;	/* The colour of each node	*/
  int    g0,nblock;
  
  if (!g->nbound) 
    return 0;

  nblock=0;
  
  nnodes=g->nnodes;
  snew(egc,nnodes);
  
  g0=g->start;
  nW=g->nbound;
  nG=0;
  nB=0;

  fW=0;

  /* We even have a loop invariant:
   * nW+nG+nB == g->nbound
   */
  
  if (bVerbose)
    fprintf(stderr,"Walking down the molecule graph to make shake-blocks\n");

  while (nW > 0) {
    /* Find the first white, this will allways be a larger
     * number than before, because no nodes are made white
     * in the loop
     */
    if ((fW=first_colour(fW,egcolWhite,g,egc)) == -1) 
      fatal_error(0,"No WHITE nodes found while nW=%d\n",nW);
    
    /* Make the first white node grey, and set the block number */
    egc[fW]        = egcolGrey;
    sid[fW+g0].sid = nblock++;
    nG++;
    nW--;

    /* Initial value for the first grey */
    fG=fW;

    if (debug)
      fprintf(debug,"Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n",
	      nW,nG,nB,nW+nG+nB);
	      
    while (nG > 0) {
      if ((fG=first_colour(fG,egcolGrey,g,egc)) == -1)
	fatal_error(0,"No GREY nodes found while nG=%d\n",nG);
      
      /* Make the first grey node black */
      egc[fG]=egcolBlack;
      nB++;
      nG--;

      /* Make all the neighbours of this black node grey
       * and set their block number
       */
      ng=mk_grey(nnodes,egc,g,&fG,sid);
      /* ng is the number of white nodes made grey */
      nG+=ng;
      nW-=ng;
    }
  }
  sfree(egc);
  
  return nblock;
}

void gen_sblocks(bool bVerbose,int natoms,t_idef *idef,t_block *sblock)
{
  t_graph *g;
  int     i,j,k;
  t_sid   *sid;
  int     isid,nsid;
  
  g=mk_graph(idef,natoms,TRUE);
  if (bVerbose && debug)
    p_graph(debug,"Graaf Dracula",g);
  snew(sid,natoms);
  for(i=0; (i<natoms); i++) {
    sid[i].atom =  i;
    sid[i].sid  = -1;
  }
  nsid=mk_sblocks(bVerbose,g,sid);
  
  /* Now sort the shake blocks... */
  qsort(sid,natoms,(size_t)sizeof(sid[0]),sid_comp);
  
  /* Now check how many are NOT -1, i.e. how many have to be shaken */
  for(i=0; (i<natoms); i++)
    if (sid[i].sid > -1)
      break;
  
  /* Fill the sblock struct */    
  sblock->nr  = nsid;
  sblock->nra = natoms-i;
  srenew(sblock->a,sblock->nra);
  srenew(sblock->index,sblock->nr+1);
  
  if (i < natoms) {
    isid = sid[i].sid;
    sblock->index[0]=0;
    if (nsid > 0) {
      k=1;
      for(j=0 ; (i<natoms); i++,j++) {
	sblock->a[j]=sid[i].atom;
	if (sid[i].sid != isid) {
	  sblock->index[k]=j;
	  isid = sid[i].sid;
	  k++;
	}
      }
      sblock->index[k]=j;
      
      if (k != nsid) {
	fatal_error(0,"k=%d, nsid=%d\n",k,nsid);
      }
    }
  }
  sfree(sid);
  done_graph(g);
  sfree(g);
}
