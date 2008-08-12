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
#include <string.h>

#include "sysstuff.h"
#include "macros.h"
#include "smalloc.h"
#include "typedefs.h"
#include "mshift.h"
#include "invblock.h"
#include "txtdump.h"
#include "math.h"
#include "assert.h"
#include "splitter.h"

typedef struct {
  int     nr;
  t_iatom *ia;
} t_sf;

static t_sf *init_sf(int nr)
{
  t_sf *sf;
  int  i;

  snew(sf,nr);
  for(i=0; (i<nr); i++) {
    sf[i].nr=0;
    sf[i].ia=NULL;
  }

  return sf;
}

static void done_sf(int nr,t_sf *sf)
{
  int i;

  for(i=0; (i<nr); i++) {
    sf[i].nr=0;
    sfree(sf[i].ia);
    sf[i].ia=NULL;
  }
  sfree(sf);
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

static void split_force2(int nnodes,int hid[],int ftype,t_ilist *ilist,
			 int *multinr)
{
  int     i,j,k,type,nodeid,nratoms,tnr;
  int     nvsite_constr,tmpid;
  t_iatom ai;
  t_sf    *sf;
  
  sf = init_sf(nnodes);

  /* Walk along all the bonded forces, find the appropriate node
   * to calc it on, and add it to that nodes list.
   */
  for (i=0; i<ilist->nr; i+=(1+nratoms)) {
    type    = ilist->iatoms[i];
    nratoms = interaction_function[ftype].nratoms;

    if (ftype == F_CONSTR) {
      /* SPECIAL CASE: All Atoms must have the same home node! */
      nodeid=hid[ilist->iatoms[i+1]];
      if (hid[ilist->iatoms[i+2]] != nodeid) 
	gmx_fatal(FARGS,"Shake block crossing node boundaries\n"
		  "constraint between atoms (%d,%d)",
		  ilist->iatoms[i+1]+1,ilist->iatoms[i+2]+1);
    }
    else if (ftype == F_SETTLE) {
      /* Only the first particle is stored for settles ... */
      ai=ilist->iatoms[i+1];
      nodeid=hid[ai];
      if ((nodeid != hid[ai+1]) ||
	  (nodeid != hid[ai+2]))
	gmx_fatal(FARGS,"Settle block crossing node boundaries\n"
		  "constraint between atoms (%d-%d)",ai+1,ai+2+1);
    }
    else if(interaction_function[ftype].flags & IF_VSITE) {
      /* Virtual sites should be constructed on the home node */
  
      if (ftype==F_VSITE2)
	nvsite_constr=2;
      else if(ftype==F_VSITE4FD || ftype==F_VSITE4FDN)
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
  for(nodeid=0; (nodeid<nnodes); nodeid++) {
    for (i=0; (i<sf[nodeid].nr); i++) 
      ilist->iatoms[tnr++]=sf[nodeid].ia[i];

    multinr[nodeid]=(nodeid==0) ? 0 : multinr[nodeid-1];
    multinr[nodeid]+=sf[nodeid].nr;
  }
  if (tnr != ilist->nr)
    gmx_incons("Splitting forces over processors");
  done_sf(nnodes,sf);
}

static int *home_index(int nnodes,t_block *cgs,int *multinr)
{
  /* This routine determines the node id for each particle */
  int *hid;
  int nodeid,j0,j1,j,k;
  
  snew(hid,cgs->index[cgs->nr]);
  /* Initiate to -1 to make it possible to check afterwards,
   * all hid's should be set in the loop below
   */
  for(k=0; (k<cgs->index[cgs->nr]); k++)
    hid[k]=-1;
    
  /* loop over nodes */
  for(nodeid=0; (nodeid<nnodes); nodeid++) {
    j0 = (nodeid==0) ? 0 : multinr[nodeid-1];
    j1 = multinr[nodeid];
    
    /* j0 and j1 are the boundariesin the index array */
    for(j=j0; (j<j1); j++) {
      for(k=cgs->index[j]; (k<cgs->index[j+1]); k++) {
	hid[k]=nodeid;
      }
    }
  }
  /* Now verify that all hid's are not -1 */
  for(k=0; (k<cgs->index[cgs->nr]); k++)
    if (hid[k] == -1)
      gmx_fatal(FARGS,"hid[%d] = -1, cgs->nr = %d, natoms = %d",
		  k,cgs->nr,cgs->index[cgs->nr]);
  
  return hid;
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
  if (fp)
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
  if (fp)
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
			 t_block *cgs,t_blocka *sblock,real capacity[],
			 int *multinr_cgs,int *multinr_shk)
{
  int      natoms,*maxatom;
  int      i,ii,ai,b0,b1;
  int      nodeid,last_shk,nbor;
  t_border *border;
  double   tload,tcap;
  
  bool    bSHK;
  atom_id *shknum,*cgsnum;
  
  natoms = cgs->index[cgs->nr];

  if (debug) {
    pr_block(debug,0,"cgs",cgs,TRUE);
    pr_blocka(debug,0,"sblock",sblock,TRUE);
  }

  shknum = make_invblocka(sblock,natoms+1);
  cgsnum = make_invblock(cgs,natoms+1);
  border = mk_border(fp,natoms,cgsnum,shknum,&nbor);

  snew(maxatom,nnodes);
  tload  = capacity[0]*natoms;
  tcap   = 1.0;
  nodeid = 0;
  /* Start at bor is 1, to force the first block on the first processor */
  for(i=0; (i<nbor) && (tload < natoms); i++) {
    if(i<(nbor-1)) 
      b1=border[i+1].atom;
    else
      b1=natoms;

    b0 = border[i].atom;
    
    if ((fabs(b0-tload)<fabs(b1-tload))) {
      /* New nodeid time */
      multinr_cgs[nodeid] = border[i].ic;
      /* Store the atom number here, has to be processed later */
      multinr_shk[nodeid] = border[i].atom;
      maxatom[nodeid]     = b0;
      tcap -= capacity[nodeid];
      nodeid++;
      
      /* Recompute target load */
      tload = b0 +
	(natoms-b0)*capacity[nodeid]/tcap;

      if (debug)
	fprintf(debug,"tload: %g tcap: %g  nodeid: %d\n",tload,tcap,nodeid);
    } 
  }
  /* Now the last one... */
  while (nodeid < nnodes) {
    multinr_cgs[nodeid] = cgs->nr;
    /* Store atom number, see above */
    multinr_shk[nodeid] = natoms;
    maxatom[nodeid]     = natoms;
    nodeid++;
  }
  if (nodeid != nnodes) {
    gmx_fatal(FARGS,"nodeid = %d, nnodes = %d, file %s, line %d",
		nodeid,nnodes,__FILE__,__LINE__);
  }

  for(i=nnodes-1; (i>0); i--)
    maxatom[i]-=maxatom[i-1];
  if (fp) {
    fprintf(fp,"Division over nodes in atoms:\n");
    for(i=0; (i<nnodes); i++)
      fprintf(fp," %7d",maxatom[i]);
    fprintf(fp,"\n");
  }

  sfree(maxatom);
  sfree(shknum);
  sfree(cgsnum);
  sfree(border);
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
		   int maxsid,t_sid sid[])
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
      range_check(aj+g0,0,maxsid);
      range_check(ai+g0,0,maxsid);
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

static int mk_sblocks(FILE *fp,t_graph *g,int maxsid,t_sid sid[])
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
  
  if (fp)
    fprintf(fp,"Walking down the molecule graph to make constraint-blocks\n");

  while (nW > 0) {
    /* Find the first white, this will allways be a larger
     * number than before, because no nodes are made white
     * in the loop
     */
    if ((fW=first_colour(fW,egcolWhite,g,egc)) == -1) 
      gmx_fatal(FARGS,"No WHITE nodes found while nW=%d\n",nW);
    
    /* Make the first white node grey, and set the block number */
    egc[fW]        = egcolGrey;
    range_check(fW+g0,0,maxsid);
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
      ng=mk_grey(nnodes,egc,g,&fG,maxsid,sid);
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


typedef struct {
  int first,last,sid;
} t_merge_sid;

static int ms_comp(const void *a, const void *b)
{
  t_merge_sid *ma = (t_merge_sid *)a;
  t_merge_sid *mb = (t_merge_sid *)b;
  int d;
  
  d = ma->first-mb->first;
  if (d == 0)
    return ma->last-mb->last;
  else
    return d;
}

static int merge_sid(int i0,int at_start,int at_end,int nsid,t_sid sid[],
		     t_blocka *sblock)
{
  int  i,j,k,n,isid,ndel;
  t_merge_sid *ms;
  int  nChanged;

  /* We try to remdy the following problem:
   * Atom: 1  2  3  4  5 6 7 8 9 10
   * Sid:  0 -1  1  0 -1 1 1 1 1 1
   */
   
  /* Determine first and last atom in each shake ID */
  snew(ms,nsid);

  for(k=0; (k<nsid); k++) {
    ms[k].first = at_end+1;
    ms[k].last  = -1;
    ms[k].sid   = k;
  }
  for(i=at_start; (i<at_end); i++) {
    isid = sid[i].sid;
    range_check(isid,-1,nsid);
    if (nsid >= 0) {
      ms[isid].first = min(ms[isid].first,sid[i].atom);
      ms[isid].last  = max(ms[isid].last,sid[i].atom);
    }
  }
  qsort(ms,nsid,sizeof(ms[0]),ms_comp);
  
  /* Now merge the overlapping ones */
  ndel = 0;
  for(k=0; (k<nsid); ) {
    for(j=k+1; (j<nsid); ) {
      if (ms[j].first <= ms[k].last) {
	ms[k].last  = max(ms[k].last,ms[j].last);
	ms[k].first = min(ms[k].first,ms[j].first);
	ms[j].sid = -1;
	ndel++;
	j++;
      }
      else {
	k = j;
	j = k+1;
      }
    }
    if (j == nsid)
      k++;
  }
  for(k=0; (k<nsid); k++) {
    while ((k < nsid-1) && (ms[k].sid == -1)) {
      for(j=k+1; (j<nsid); j++) {
	memcpy(&(ms[j-1]),&(ms[j]),sizeof(ms[0]));
      }
      nsid--;
    }
  }

  for(k=at_start; (k<at_end); k++) {
    sid[k].atom = k;
    sid[k].sid = -1;
  }
  sblock->nr  = nsid;
  sblock->nalloc_index = sblock->nr+1;
  snew(sblock->index,sblock->nalloc_index);
  sblock->nra = at_end - at_start;
  sblock->nalloc_a = sblock->nra;
  snew(sblock->a,sblock->nalloc_a);
  sblock->index[0] = 0;
  for(k=n=0; (k<nsid); k++) {
    sblock->index[k+1] = sblock->index[k] + ms[k].last - ms[k].first+1;
    for(j=ms[k].first; (j<=ms[k].last); j++) {
      range_check(n,0,sblock->nra);
      sblock->a[n++] = j;
      range_check(j,0,at_end);
      if (sid[j].sid == -1)
	sid[j].sid = k;
      else 
	fprintf(stderr,"Double sids (%d, %d) for atom %d\n",sid[j].sid,k,j);
    }
  }
  assert(k == nsid);
  /* Removed 2007-09-04
  sblock->index[k+1] = natoms;
  for(k=0; (k<natoms); k++)
    if (sid[k].sid == -1)
      sblock->a[n++] = k;
  assert(n == natoms);
  */
  sblock->nra = n;
  assert(sblock->index[k] == sblock->nra);
  sfree(ms);
  
  return nsid;
}

void gen_sblocks(FILE *fp,int at_start,int at_end,
		 t_idef *idef,t_blocka *sblock,
		 bool bSettle)
{
  t_graph *g;
  int     i,i0,j,k,istart,n;
  t_sid   *sid;
  int     isid,nsid;
  
  g=mk_graph(NULL,idef,at_start,at_end,TRUE,bSettle);
  if (debug)
    p_graph(debug,"Graaf Dracula",g);
  snew(sid,at_end);
  for(i=at_start; (i<at_end); i++) {
    sid[i].atom =  i;
    sid[i].sid  = -1;
  }
  nsid=mk_sblocks(fp,g,at_end,sid);

  if (!nsid)
    return;
    
  /* Now sort the shake blocks... */
  qsort(sid+at_start,at_end-at_start,(size_t)sizeof(sid[0]),sid_comp);
  
  if (debug) {
    fprintf(debug,"Sorted shake block\n");
    for(i=at_start; (i<at_end); i++) 
      fprintf(debug,"sid[%5d] = atom:%5d sid:%5d\n",i,sid[i].atom,sid[i].sid);
  }
  /* Now check how many are NOT -1, i.e. how many have to be shaken */
  for(i0=at_start; (i0<at_end); i0++)
    if (sid[i0].sid > -1)
      break;
  
  /* Now we have the sids that have to be shaken. We'll check the min and
   * max atom numbers and this determines the shake block. DvdS 2007-07-19.
   * For the purpose of making boundaries all atoms in between need to be
   * part of the shake block too. There may be cases where blocks overlap
   * and they will have to be merged.
   */
  nsid = merge_sid(i0,at_start,at_end,nsid,sid,sblock);
  /* Now sort the shake blocks again... */
  /*qsort(sid,natoms,(size_t)sizeof(sid[0]),sid_comp);*/
  
  /* Fill the sblock struct */    
  /*  sblock->nr  = nsid;
  sblock->nra = natoms;
  srenew(sblock->a,sblock->nra);
  srenew(sblock->index,sblock->nr+1);
  
  i    = i0;
  isid = sid[i].sid;
  n    = k = 0;
  sblock->index[n++]=k;
  while (i < natoms) {
    istart = sid[i].atom;
    while ((i<natoms-1) && (sid[i+1].sid == isid))
    i++;*/
    /* After while: we found a new block, or are thru with the atoms */
  /*    for(j=istart; (j<=sid[i].atom); j++,k++)
      sblock->a[k]=j;
    sblock->index[n] = k;
    if (i < natoms-1)
      n++;
    if (n > nsid)
      gmx_fatal(FARGS,"Death Horror: nsid = %d, n= %d",nsid,n);
    i++;
    isid = sid[i].sid;
  }
  */
  sfree(sid);
  /* Due to unknown reason this free generates a problem sometimes */
  done_graph(g);
  sfree(g); 
  if (debug)
    fprintf(debug,"Done gen_sblocks\n");
}

static t_blocka block2blocka(t_block *block)
{
  t_blocka blocka;
  int i;

  blocka.nr = block->nr;
  blocka.nalloc_index = blocka.nr + 1;
  snew(blocka.index,blocka.nalloc_index);
  for(i=0; i<=block->nr; i++)
    blocka.index[i] = block->index[i];
  blocka.nra = block->index[block->nr];
  blocka.nalloc_a = blocka.nra;
  snew(blocka.a,blocka.nalloc_a);
  for(i=0; i<blocka.nra; i++)
    blocka.a[i] = i;

  return blocka;
}

void split_top(FILE *fp,int nnodes,gmx_localtop_t *top,t_block *mols,
	       real *capacity,int *multinr_cgs,int **multinr_nre)
{
  int     natoms,i,j,k,mj,atom,maxatom,sstart,send,bstart,nodeid;
  int     *multinr_shk;
  t_blocka sblock,shakeblock;
  int     *homeind;
  atom_id *sblinv;
  int ftype,nvsite_constr,nra,nrd;
  t_iatom   *ia;
  int minhome,ihome,minidx;
  
  if (nnodes <= 1)
    return;
  
  natoms = mols->index[mols->nr];

  if (fp)
    fprintf(fp,"splitting topology...\n");
  
#define MOL_BORDER
#ifndef MOL_BORDER
  /* Make a special shake block that includes settles */
  init_block(&sblock);
  gen_sblocks(fp,0,natoms,&top->idef,&sblock,TRUE);
#else
  sblock = block2blocka(mols);
#endif  

  snew(multinr_shk,nnodes);
  split_blocks(fp,nnodes,&top->cgs,&sblock,capacity,multinr_cgs,multinr_shk);
  
  /* Now transform atom numbers to real inverted shake blocks */
  init_blocka(&shakeblock);
  gen_sblocks(fp,0,natoms,&top->idef,&shakeblock,TRUE);
  sblinv = make_invblocka(&shakeblock,natoms+1);
  done_blocka(&shakeblock);
  for(j=0; (j<nnodes); j++) {
    atom = multinr_shk[j];
    mj   = NO_ATID;
    for(k=(j == 0) ? 0 : multinr_shk[j-1]; (k<atom); k++)
      if (sblinv[k] != NO_ATID)
	mj = max(mj,(int)sblinv[k]);
    if (mj == NO_ATID) 
      mj = (j == 0) ? -1 : multinr_shk[j-1]-1;
    
    multinr_shk[j] = mj+1;
  }
  homeind=home_index(nnodes,&top->cgs,multinr_cgs);
  
  for(j=0; (j<F_NRE); j++)
    split_force2(nnodes,homeind,j,&top->idef.il[j],multinr_nre[j]);

  sfree(sblinv);
  sfree(multinr_shk);
  sfree(homeind);
  done_blocka(&sblock);
}

