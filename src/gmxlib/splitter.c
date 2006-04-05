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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "sysstuff.h"
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "mshift.h"
#include "invblock.h"
#include "txtdump.h"
#include "math.h"

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

static int min_nodeid(int nr,atom_id list[],int hid[])
{
  int i,nodeid,minnodeid;

  if (nr <= 0)
    gmx_incons("Invalid node number");
  minnodeid=hid[list[0]];
  for (i=1; (i<nr); i++) 
    if ((nodeid=hid[list[i]]) < minnodeid) 
      minnodeid=nodeid;

  return minnodeid;
}

static void split_force2(int nnodes,int hid[],t_idef *idef,t_ilist *ilist)
{
  int     i,j,k,type,ftype,nodeid,nratoms,tnr;
  int     nvsite_constr,tmpid;
  t_iatom ai;
  t_sf    sf[MAXNODES];
  
  init_sf(MAXNODES,sf);

  /* Walk along all the bonded forces, find the appropriate node
   * to calc it on, and add it to that nodes list.
   */
  for (i=0; i<ilist->nr; i+=(1+nratoms)) {
    type    = ilist->iatoms[i];
    ftype   = idef->functype[type];
    nratoms = interaction_function[ftype].nratoms;

    if (ftype == F_CONSTR) {
      /* SPECIAL CASE: All Atoms must have the same home node! */
      nodeid=hid[ilist->iatoms[i+1]];
      if (hid[ilist->iatoms[i+2]] != nodeid) 
	gmx_fatal(FARGS,"Shake block crossing node boundaries\n"
		    "constraint between atoms (%d,%d)",
		    ilist->iatoms[i+1],ilist->iatoms[i+2]);
    }
    else if (ftype == F_SETTLE) {
      /* Only the first particle is stored for settles ... */
      ai=ilist->iatoms[i+1];
      nodeid=hid[ai];
      if ((nodeid != hid[ai+1]) ||
	  (nodeid != hid[ai+2]))
	gmx_fatal(FARGS,"Settle block crossing node boundaries\n"
		    "constraint between atoms (%d-%d)",ai,ai+2);
    }
    else if(interaction_function[ftype].flags & IF_VSITE) {
      /* Virtual sites should be constructed on the home node */
  
      if (ftype==F_VSITE2)
	nvsite_constr=2;
      else if(ftype==F_VSITE4FD)
	nvsite_constr=4;
      else
	nvsite_constr=3;
      
      tmpid=hid[ilist->iatoms[i+1]];
      
      for(k=2;k<nvsite_constr+2;k++) {
	if(hid[ilist->iatoms[i+k]]<(tmpid-1) ||
	   hid[ilist->iatoms[i+k]]>(tmpid+1))
	  gmx_fatal(FARGS,"Virtual site %d and its constructing"
		      " atoms are not on the same or adjacent\n" 
		      " nodes. This is necessary to avoid a lot\n"
		      " of extra communication. The easiest way"
		      " to ensure this is to place virtual sites\n"
		      " close to the constructing atoms.\n"
		      " Sorry, but you will have to rework your topology!\n",
		      ilist->iatoms[i+1]);
      }
      nodeid=min_nodeid(nratoms,&ilist->iatoms[i+1],hid);
    } else
      nodeid=min_nodeid(nratoms,&ilist->iatoms[i+1],hid);

    /* Add it to the list */
    push_sf(&(sf[nodeid]),nratoms+1,&(ilist->iatoms[i]));
  }
  tnr=0;
  for(nodeid=0; (nodeid<MAXNODES); nodeid++) {
    for (i=0; (i<sf[nodeid].nr); i++) 
      ilist->iatoms[tnr++]=sf[nodeid].ia[i];

    ilist->multinr[nodeid]=(nodeid==0) ? 0 : ilist->multinr[nodeid-1];
    ilist->multinr[nodeid]+=sf[nodeid].nr;
  }
  if (tnr != ilist->nr)
    gmx_incons("Splitting forces over processors");
  done_sf(MAXNODES,sf);
}

static int *home_index(int nnodes,t_block *cgs)
{
  /* This routine determines the node id for each particle */
  int *hid;
  int nodeid,j0,j1,j,k,ak;
  
  snew(hid,cgs->nra);
  /* Initiate to -1 to make it possible to check afterwards,
   * all hid's should be set in the loop below
   */
  for(k=0; (k<cgs->nra); k++)
    hid[k]=-1;
    
  /* loop over nodes */
  for(nodeid=0; (nodeid<nnodes); nodeid++) {
    j0 = (nodeid==0) ? 0 : cgs->multinr[nodeid-1];
    j1 = cgs->multinr[nodeid];
    
    /* j0 and j1 are the boundariesin the index array */
    for(j=j0; (j<j1); j++) {
      for(k=cgs->index[j]; (k<cgs->index[j+1]); k++) {
	ak=cgs->a[k];
	hid[ak]=nodeid;
      }
    }
  }
  /* Now verify that all hid's are not -1 */
  for(k=0; (k<cgs->nra); k++)
    if (hid[k] == -1)
      gmx_fatal(FARGS,"hid[%d] = -1, cgs->nr = %d, cgs->nra = %d",
		  k,cgs->nr,cgs->nra);
  
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

typedef struct {
  int atom,ic,is;
} t_border;

void set_bor(t_border *b,int atom,int ic,int is)
{
  if (debug)
    fprintf(debug,"border @ atom %5d [ ic = %5d,  is = %5d ]\n",atom,ic,is);
  b->atom = atom;
  b->ic   = ic;
  b->is   = is;
}

static bool is_bor(atom_id ai[],int i)
{
  return ((ai[i] != ai[i-1]) || ((ai[i] == NO_ATID) && (ai[i-1] == NO_ATID)));
}

static t_border *mk_border(FILE *fp,int natom,atom_id *invcgs,
			   atom_id *invshk,int *nb)
{
  t_border *bor;
  atom_id  *sbor,*cbor;
  int      i,j,is,ic,ns,nc,nbor;

  if (debug) {
    for(i=0; (i<natom); i++) {
      fprintf(debug,"atom: %6d  cgindex: %6d  shkindex: %6d\n",
	      i,invcgs[i],invshk[i]);
    }
  }
    
  snew(sbor,natom+1);
  snew(cbor,natom+1);
  ns = nc = 1;
  for(i=1; (i<natom); i++) {
    if (is_bor(invcgs,i)) 
      cbor[nc++] = i;
    if (is_bor(invshk,i)) 
      sbor[ns++] = i;
  }
  sbor[ns] = 0;
  cbor[nc] = 0;
  fprintf(fp,"There are %d charge group borders and %d shake borders\n",
	  nc,ns);
  snew(bor,max(nc,ns));
  ic = is = nbor = 0;
  while ((ic < nc) || (is < ns)) {
    if (sbor[is] == cbor[ic]) {
      set_bor(&(bor[nbor]),cbor[ic],ic,is);
      nbor++;
      if (ic < nc) ic++;
      if (is < ns) is++;
    }
    else if (cbor[ic] > sbor[is]) {
      if (is == ns) {
	set_bor(&(bor[nbor]),cbor[ic],ic,is);
	nbor++;
	if (ic < nc) ic++;
      }
      else if (is < ns) 
	is++;
    }
    else if (ic < nc)
      ic++;
    else
      is++;/*gmx_fatal(FARGS,"Can't happen is=%d, ic=%d (%s, %d)",
	     is,ic,__FILE__,__LINE__);*/
  }
  fprintf(fp,"There are %d total borders\n",nbor);

  if (debug) {
    fprintf(debug,"There are %d actual bor entries\n",nbor);
    for(i=0; (i<nbor); i++) 
      fprintf(debug,"bor[%5d] = atom: %d  ic: %d  is: %d\n",i,
	      bor[i].atom,bor[i].ic,bor[i].is);
  }
  
  *nb = nbor;
  
  return bor;
}

static void split_blocks(FILE *fp,int nnodes,
			 t_block *cgs,t_block *sblock,real capacity[])
{
  int      maxatom[MAXNODES];
  int      i,ii,ai,b0,b1;
  int      nodeid,last_shk,nbor;
  t_border *border;
  double   tload,tcap;
  
  bool    bSHK;
  atom_id *shknum,*cgsnum;
  
  if (debug) {
    pr_block(debug,0,"cgs",cgs,TRUE);
    pr_block(debug,0,"sblock",sblock,TRUE);
  }

  shknum = make_invblock(sblock,cgs->nra+1);
  cgsnum = make_invblock(cgs,cgs->nra+1);
  border = mk_border(fp,cgs->nra,cgsnum,shknum,&nbor);

  tload  = capacity[0]*cgs->nra;
  tcap   = 1.0;
  nodeid = 0;
  /* Start at bor is 1, to force the first block on the first processor */
  for(i=0; (i<nbor) && (tload < cgs->nra); i++) {
    if(i<(nbor-1)) 
      b1=border[i+1].atom;
    else
      b1=cgs->nra;

    b0 = border[i].atom;
    
    if ((fabs(b0-tload)<fabs(b1-tload))) {
      /* New nodeid time */
      cgs->multinr[nodeid]    = border[i].ic;
      /* Store the atom number here, has to be processed later */
      sblock->multinr[nodeid] = border[i].atom;
      maxatom[nodeid]         = b0;
      tcap -= capacity[nodeid];
      nodeid++;
      
      /* Recompute target load */
      tload = b0 +
	(cgs->nra-b0)*capacity[nodeid]/tcap;

      if (debug)
	fprintf(debug,"tload: %g tcap: %g  nodeid: %d\n",tload,tcap,nodeid);
    } 
  }
  /* Now the last one... */
  while (nodeid < nnodes) {
    cgs->multinr[nodeid]    = cgs->nr;
    /* Store atom number, see above */
    sblock->multinr[nodeid] = cgs->nra;
    maxatom[nodeid]         = cgs->nra;
    nodeid++;
  }
  if (nodeid != nnodes) {
    gmx_fatal(FARGS,"nodeid = %d, nnodes = %d, file %s, line %d",
		nodeid,nnodes,__FILE__,__LINE__);
  }

  for(i=nnodes-1; (i>0); i--)
    maxatom[i]-=maxatom[i-1];
  fprintf(fp,"Division over nodes in atoms:\n");
  for(i=0; (i<nnodes); i++)
    fprintf(fp," %7d",maxatom[i]);
  fprintf(fp,"\n");
  
  sfree(shknum);
  sfree(cgsnum);
  sfree(border);
}

static void def_mnr(int nr,int mnr[])
{
  int i;

  for (i=0; (i<MAXNODES); i++) 
    mnr[i]=0;
  mnr[0]=nr;
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
	gmx_fatal(FARGS,"sid[%d]=%d, sid[%d]=%d, file %s, line %d",
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

static int mk_sblocks(FILE *fp,t_graph *g,t_sid sid[])
{
  int    ng,nnodes;
  int    nW,nG,nB;		/* Number of Grey, Black, White	*/
  int    fW,fG;			/* First of each category	*/
  egCol  *egc=NULL;	        /* The colour of each node	*/
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
  
  fprintf(fp,"Walking down the molecule graph to make shake-blocks\n");

  while (nW > 0) {
    /* Find the first white, this will allways be a larger
     * number than before, because no nodes are made white
     * in the loop
     */
    if ((fW=first_colour(fW,egcolWhite,g,egc)) == -1) 
      gmx_fatal(FARGS,"No WHITE nodes found while nW=%d\n",nW);
    
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
	gmx_fatal(FARGS,"No GREY nodes found while nG=%d\n",nG);
      
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

  if (debug)
    fprintf(debug,"Found %d shake blocks\n",nblock);
    
  return nblock;
}

void gen_sblocks(FILE *fp,int natoms,t_idef *idef,t_block *sblock,
		 bool bSettle)
{
  t_graph *g;
  int     i,j,k;
  t_sid   *sid;
  int     isid,nsid;
  
  g=mk_graph(idef,natoms,TRUE,bSettle);
  if (debug)
    p_graph(debug,"Graaf Dracula",g);
  snew(sid,natoms);
  for(i=0; (i<natoms); i++) {
    sid[i].atom =  i;
    sid[i].sid  = -1;
  }
  nsid=mk_sblocks(fp,g,sid);

  if (!nsid)
    return;
    
  /* Now sort the shake blocks... */
  qsort(sid,natoms,(size_t)sizeof(sid[0]),sid_comp);
  
  if (debug) {
    fprintf(debug,"Sorted shake block\n");
    for(i=0; (i<natoms); i++) 
      fprintf(debug,"sid[%5d] = atom:%5d sid:%5d\n",i,sid[i].atom,sid[i].sid);
  }
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
	gmx_fatal(FARGS,"k=%d, nsid=%d\n",k,nsid);
      }
    }
  }
  sfree(sid);
  /* Due to unknown reason this free generates a problem sometimes */
  /* done_graph(g);
     sfree(g); */
  if (debug)
    fprintf(debug,"Done gen_sblocks\n");
}

void reset_multinr(t_topology *top)
{
  int j;
  
  for(j=0; (j<F_NRE); j++)  
    def_mnr(top->idef.il[j].nr,top->idef.il[j].multinr);
  
  for (j=0; j<ebNR; j++) 
    def_mnr(top->blocks[j].nr,top->blocks[j].multinr);
}

void split_top(FILE *fp,int nnodes,t_topology *top,real *capacity)
{
  int     i,j,k,mj,atom,maxatom;
  t_block sblock;
  int     *homeind;
  atom_id *sblinv;
  int ftype,nvsite_constr,nra,nrd;
  t_iatom   *ia;
  int minhome,ihome,minidx;
  
  if (nnodes <= 1)
    return;
    
  fprintf(fp,"splitting topology...\n");
  
  reset_multinr(top);
  
  /* Make a special shake block that includes settles */
  init_block(&sblock);
  gen_sblocks(fp,top->atoms.nr,&top->idef,&sblock,TRUE);
  
  split_blocks(fp,nnodes,&(top->blocks[ebCGS]),&sblock,capacity);
  
  /* Now transform atom numbers to real inverted shake blocks */
  sblinv = make_invblock(&(top->blocks[ebSBLOCKS]),top->atoms.nr+1);
  for(j=0; (j<MAXNODES); j++) {
    atom = sblock.multinr[j];
    mj   = NO_ATID;
    for(k=(j == 0) ? 0 : sblock.multinr[j-1]; (k<atom); k++)
      if (sblinv[k] != NO_ATID)
	mj = max(mj,(int)sblinv[k]);
    if (mj == NO_ATID) 
      mj = (j == 0) ? -1 : top->blocks[ebSBLOCKS].multinr[j-1]-1;
    
    top->blocks[ebSBLOCKS].multinr[j] = mj+1;
  }
  sfree(sblinv);
    
  homeind=home_index(nnodes,&(top->blocks[ebCGS]));
  
  for(j=0; (j<F_NRE); j++)
    split_force2(nnodes,homeind,&top->idef,&top->idef.il[j]);
  sfree(homeind);
  done_block(&sblock);
}

