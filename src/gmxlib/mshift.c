/*
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "vec.h"
#include "futil.h"
#include "copyrite.h"
#include "mshift.h"
#include "main.h"
#include "pbc.h"

/************************************************************
 *
 *      S H I F T   U T I L I T I E S
 *
 ************************************************************/
 

/************************************************************
 *
 *      G R A P H   G E N E R A T I O N    C O D E
 *
 ************************************************************/

static void add_gbond(t_graph *g,atom_id a0,atom_id a1)
{
  int     i;
  atom_id inda0,inda1;
  gmx_bool    bFound;

  inda0 = a0 - g->start;
  inda1 = a1 - g->start;
  bFound = FALSE;
  /* Search for a direct edge between a0 and a1.
   * All egdes are bidirectional, so we only need to search one way.
   */
  for(i=0; (i<g->nedge[inda0] && !bFound); i++) {
    bFound = (g->edge[inda0][i] == a1);
  }

  if (!bFound) {
    g->edge[inda0][g->nedge[inda0]++] = a1;
    g->edge[inda1][g->nedge[inda1]++] = a0;
  }
}

static void mk_igraph(t_graph *g,int ftype,t_ilist *il,
		      int at_start,int at_end,
		      int *part)
{
  t_iatom *ia;
  int     i,j,np;
  int     end;

  end=il->nr;
  ia=il->iatoms;

  i = 0;
  while (i < end) {
    np = interaction_function[ftype].nratoms;
    
    if (ia[1] >= at_start && ia[1] < at_end) {
      if (ia[np] >= at_end)
	gmx_fatal(FARGS,
		  "Molecule in topology has atom numbers below and "
		  "above natoms (%d).\n"
		  "You are probably trying to use a trajectory which does "
		  "not match the first %d atoms of the run input file.\n"
		  "You can make a matching run input file with tpbconv.",
		  at_end,at_end);
      if (ftype == F_SETTLE) {
	/* Bond all the atoms in the settle */
	add_gbond(g,ia[1],ia[1]+1);
	add_gbond(g,ia[1],ia[1]+2);
      } else if (part == NULL) {
	/* Simply add this bond */
	for(j=1; j<np; j++) {
	  add_gbond(g,ia[j],ia[j+1]);
	}
      } else {
	/* Add this bond when it connects two unlinked parts of the graph */
	for(j=1; j<np; j++) {
	  if (part[ia[j]] != part[ia[j+1]]) {
	    add_gbond(g,ia[j],ia[j+1]);
	  }
	}
      }
    }
    ia+=np+1;
    i+=np+1;
  }
}

static void g_error(int line,const char *file)
{
  gmx_fatal(FARGS,"Tring to print non existant graph (file %s, line %d)",
	      file,line);
}
#define GCHECK(g) if (g == NULL) g_error(__LINE__,__FILE__)

void p_graph(FILE *log,const char *title,t_graph *g)
{
  int  i,j;
  const char *cc[egcolNR] = { "W", "G", "B" };
  
  GCHECK(g);
  fprintf(log,"graph:  %s\n",title);
  fprintf(log,"nnodes: %d\n",g->nnodes);
  fprintf(log,"nbound: %d\n",g->nbound);
  fprintf(log,"start:  %d\n",g->start);
  fprintf(log,"end:    %d\n",g->end);
  fprintf(log," atom shiftx shifty shiftz C nedg    e1    e2 etc.\n");
  for(i=0; (i<g->nnodes); i++)
    if (g->nedge[i] > 0) {
      fprintf(log,"%5d%7d%7d%7d %1s%5d",g->start+i+1,
	      g->ishift[i][XX],g->ishift[i][YY],
	      g->ishift[i][ZZ],
	      (g->negc > 0) ? cc[g->egc[i]] : " ",
	      g->nedge[i]);
      for(j=0; (j<g->nedge[i]); j++)
	fprintf(log," %5u",g->edge[i][j]+1);
      fprintf(log,"\n");
    }
  fflush(log);
}

static void calc_1se(t_graph *g,int ftype,t_ilist *il,
		     int nbond[],int at_start,int at_end)
{
  int     k,nratoms,end,j;
  t_iatom *ia,iaa;

  end=il->nr;

  ia=il->iatoms;
  for(j=0; (j<end); j+=nratoms+1,ia+=nratoms+1) {
    nratoms = interaction_function[ftype].nratoms;
    
    if (ftype == F_SETTLE) {
      iaa          = ia[1];
      if (iaa >= at_start && iaa < at_end) {
	nbond[iaa]   += 2;
	nbond[iaa+1] += 1;
	nbond[iaa+2] += 1;
	g->start      = min(g->start,iaa);
	g->end        = max(g->end,iaa+2);
      }
    } else {
      for(k=1; (k<=nratoms); k++) {
	iaa=ia[k];
	if (iaa >= at_start && iaa < at_end) {
	  g->start=min(g->start,iaa);
	  g->end  =max(g->end,  iaa);
	  /* When making the graph we (might) link all atoms in an interaction
	   * sequentially. Therefore the end atoms add 1 to the count,
	   * the middle atoms 2.
	   */
	  if (k == 1 || k == nratoms) {
	    nbond[iaa] += 1;
	  } else {
	    nbond[iaa] += 2;
	  }
	}
      }
    }
  }
}

static int calc_start_end(FILE *fplog,t_graph *g,t_ilist il[],
			  int at_start,int at_end,
			  int nbond[])
{
  int   i,nnb,nbtot;
  
  g->start=at_end;
  g->end=0;

  /* First add all the real bonds: they should determine the molecular 
   * graph.
   */
  for(i=0; (i<F_NRE); i++)
    if (interaction_function[i].flags & IF_CHEMBOND)
      calc_1se(g,i,&il[i],nbond,at_start,at_end);
  /* Then add all the other interactions in fixed lists, but first
   * check to see what's there already.
   */
  for(i=0; (i<F_NRE); i++) {
    if (!(interaction_function[i].flags & IF_CHEMBOND)) {
      calc_1se(g,i,&il[i],nbond,at_start,at_end);
    }
  }
  
  nnb   = 0;
  nbtot = 0;
  for(i=g->start; (i<=g->end); i++) {
    nbtot += nbond[i];
    nnb    = max(nnb,nbond[i]);
  }
  if (fplog) {
    fprintf(fplog,"Max number of connections per atom is %d\n",nnb);
    fprintf(fplog,"Total number of connections is %d\n",nbtot);
  }
  return nbtot;
}



static void compact_graph(FILE *fplog,t_graph *g)
{
  int i,j,n,max_nedge;
  atom_id *e;

  max_nedge = 0;
  n = 0;
  for(i=0; i<g->nnodes; i++) {
    for(j=0; j<g->nedge[i]; j++) {
      g->edge[0][n++] = g->edge[i][j];
    }
    max_nedge = max(max_nedge,g->nedge[i]);
  }
  srenew(g->edge[0],n);
  /* set pointers after srenew because edge[0] might move */
  for(i=1; i<g->nnodes; i++) {
    g->edge[i] = g->edge[i-1] + g->nedge[i-1];
  }

  if (fplog) {
    fprintf(fplog,"Max number of graph edges per atom is %d\n",
	    max_nedge);
    fprintf(fplog,"Total number of graph edges is %d\n",n);
  }
}

static gmx_bool determine_graph_parts(t_graph *g,int *part)
{
  int  i,e;
  int  nchanged;
  atom_id at_i,*at_i2;
  gmx_bool bMultiPart;

  /* Initialize the part array with all entries different */
  for(at_i=g->start; at_i<g->end; at_i++) {
    part[at_i] = at_i;
  }

  /* Loop over the graph until the part array is fixed */
  do {
    bMultiPart = FALSE;
    nchanged = 0;
    for(i=0; (i<g->nnodes); i++) {
      at_i  = g->start + i;
      at_i2 = g->edge[i];
      for(e=0; e<g->nedge[i]; e++) {
	/* Set part for both nodes to the minimum */
	if (part[at_i2[e]] > part[at_i]) {
	  part[at_i2[e]] = part[at_i];
	  nchanged++;
	} else if (part[at_i2[e]] < part[at_i]) {
	  part[at_i] = part[at_i2[e]];
	  nchanged++;
	}
      }
      if (part[at_i] != part[g->start]) {
	bMultiPart = TRUE;
      }
    }
    if (debug) {
      fprintf(debug,"graph part[] nchanged=%d, bMultiPart=%d\n",
	      nchanged,bMultiPart);
    }
  } while (nchanged > 0);

  return bMultiPart;
}

void mk_graph_ilist(FILE *fplog,
		    t_ilist *ilist,int at_start,int at_end,
		    gmx_bool bShakeOnly,gmx_bool bSettle,
		    t_graph *g)
{
  int     *nbond;
  int     i,nbtot;
  gmx_bool    bMultiPart;

  snew(nbond,at_end);
  nbtot = calc_start_end(fplog,g,ilist,at_start,at_end,nbond);
  
  if (g->start >= g->end) {
    g->nnodes = 0;
    g->nbound = 0;
  }
  else {
    g->nnodes = g->end - g->start + 1;
    snew(g->ishift,g->nnodes);
    snew(g->nedge,g->nnodes);
    snew(g->edge,g->nnodes);
    /* Allocate a single array and set pointers into it */
    snew(g->edge[0],nbtot);
    for(i=1; (i<g->nnodes); i++) {
      g->edge[i] = g->edge[i-1] + nbond[g->start+i-1];
    }

    if (!bShakeOnly) {
      /* First add all the real bonds: they should determine the molecular 
       * graph.
       */
      for(i=0; (i<F_NRE); i++)
	if (interaction_function[i].flags & IF_CHEMBOND)
	  mk_igraph(g,i,&(ilist[i]),at_start,at_end,NULL);

      /* Determine of which separated parts the IF_CHEMBOND graph consists.
       * Store the parts in the nbond array.
       */
      bMultiPart = determine_graph_parts(g,nbond);

      if (bMultiPart) {
	/* Then add all the other interactions in fixed lists,
	 * but only when they connect parts of the graph
	 * that are not connected through IF_CHEMBOND interactions.
	 */	 
	for(i=0; (i<F_NRE); i++) {
	  if (!(interaction_function[i].flags & IF_CHEMBOND)) {
	    mk_igraph(g,i,&(ilist[i]),at_start,at_end,nbond);
	  }
	}
      }
      
      /* Removed all the unused space from the edge array */
      compact_graph(fplog,g);
    }
    else {
      /* This is a special thing used in splitter.c to generate shake-blocks */
      mk_igraph(g,F_CONSTR,&(ilist[F_CONSTR]),at_start,at_end,NULL);
      if (bSettle)
	mk_igraph(g,F_SETTLE,&(ilist[F_SETTLE]),at_start,at_end,NULL);
    }
    g->nbound = 0;
    for(i=0; (i<g->nnodes); i++)
      if (g->nedge[i] > 0)
        g->nbound++;
  }

  g->negc = 0;
  g->egc = NULL;

  sfree(nbond);

  if (gmx_debug_at)
    p_graph(debug,"graph",g);
}

t_graph *mk_graph(FILE *fplog,
		  t_idef *idef,int at_start,int at_end,
		  gmx_bool bShakeOnly,gmx_bool bSettle)
{
  t_graph *g;

  snew(g,1);

  mk_graph_ilist(fplog,idef->il,at_start,at_end,bShakeOnly,bSettle,g);

  return g;
}

void done_graph(t_graph *g)
{
  int i;
  
  GCHECK(g);
  if (g->nnodes > 0) {
    sfree(g->ishift);
    sfree(g->nedge);
    sfree(g->edge[0]);
    sfree(g->edge);
    sfree(g->egc);
  }
}

/************************************************************
 *
 *      S H I F T   C A L C U L A T I O N   C O D E
 *
 ************************************************************/

static void mk_1shift_tric(int npbcdim,matrix box,rvec hbox,
			   rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for triclinic box... */
  int  m,d;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  mj[ZZ] = 0;
  for(m=npbcdim-1; (m>=0); m--) {
    /* If dx < hbox, then xj will be reduced by box, so that
     * xi - xj will be bigger
     */
    if (dx[m] < -hbox[m]) {
      mj[m]=mi[m]-1;
      for(d=m-1; d>=0; d--)
	dx[d]+=box[m][d];
    } else if (dx[m] >= hbox[m]) {
      mj[m]=mi[m]+1;
      for(d=m-1; d>=0; d--)
	dx[d]-=box[m][d];
    } else
      mj[m]=mi[m];
  }
}

static void mk_1shift(int npbcdim,rvec hbox,rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for rectangular box... */
  int  m;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  mj[ZZ] = 0;
  for(m=0; (m<npbcdim); m++) {
    /* If dx < hbox, then xj will be reduced by box, so that
     * xi - xj will be bigger
     */
    if (dx[m] < -hbox[m])
      mj[m]=mi[m]-1;
    else if (dx[m] >= hbox[m])
      mj[m]=mi[m]+1;
    else
      mj[m]=mi[m];
  }
}

static void mk_1shift_screw(matrix box,rvec hbox,
			    rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for rectangular box... */
  int  signi,m;
  rvec dx;

  if ((mi[XX] > 0 &&  mi[XX] % 2 == 1) ||
      (mi[XX] < 0 && -mi[XX] % 2 == 1)) {
    signi = -1;
  } else {
    signi =  1;
  }

  rvec_sub(xi,xj,dx);

  if (dx[XX] < -hbox[XX])
    mj[XX] = mi[XX] - 1;
  else if (dx[XX] >= hbox[XX])
    mj[XX] = mi[XX] + 1;
  else
    mj[XX] = mi[XX];
  if (mj[XX] != mi[XX]) {
    /* Rotate */
    dx[YY] = xi[YY] - (box[YY][YY] + box[ZZ][YY] - xj[YY]);
    dx[ZZ] = xi[ZZ] - (box[ZZ][ZZ]               - xj[ZZ]);
  }
  for(m=1; (m<DIM); m++) {
    /* The signs are taken such that we can first shift x and rotate
     * and then shift y and z.
     */
    if (dx[m] < -hbox[m])
      mj[m] = mi[m] - signi;
    else if (dx[m] >= hbox[m])
      mj[m] = mi[m] + signi;
    else
      mj[m] = mi[m];
  }
}

static int mk_grey(FILE *log,int nnodes,egCol egc[],t_graph *g,int *AtomI,
		   int npbcdim,matrix box,rvec x[],int *nerror)
{
  int      m,j,ng,ai,aj,g0;
  rvec     dx,hbox;
  gmx_bool     bTriclinic;
  ivec     is_aj;
  t_pbc    pbc;
   
  for(m=0; (m<DIM); m++)
    hbox[m]=box[m][m]*0.5;
  bTriclinic = TRICLINIC(box);
  
  ng=0;
  ai=*AtomI;
  
  g0=g->start;
  /* Loop over all the bonds */
  for(j=0; (j<g->nedge[ai]); j++) {
    aj=g->edge[ai][j]-g0;
    /* If there is a white one, make it grey and set pbc */
    if (g->bScrewPBC)
      mk_1shift_screw(box,hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    else if (bTriclinic)
      mk_1shift_tric(npbcdim,box,hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    else
      mk_1shift(npbcdim,hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    
    if (egc[aj] == egcolWhite) {
      if (aj < *AtomI)
	*AtomI = aj;
      egc[aj] = egcolGrey;
      
      copy_ivec(is_aj,g->ishift[aj]);

      ng++;
    }
    else if ((is_aj[XX] != g->ishift[aj][XX]) ||
	     (is_aj[YY] != g->ishift[aj][YY]) ||
	     (is_aj[ZZ] != g->ishift[aj][ZZ])) {
      if (gmx_debug_at) {
	set_pbc(&pbc,-1,box);
	pbc_dx(&pbc,x[g0+ai],x[g0+aj],dx);
	fprintf(debug,"mk_grey: shifts for atom %d due to atom %d\n"
		"are (%d,%d,%d), should be (%d,%d,%d)\n"
		"dx = (%g,%g,%g)\n",
		aj+g0+1,ai+g0+1,is_aj[XX],is_aj[YY],is_aj[ZZ],
		g->ishift[aj][XX],g->ishift[aj][YY],g->ishift[aj][ZZ],
		dx[XX],dx[YY],dx[ZZ]);
      }
      (*nerror)++;
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

void mk_mshift(FILE *log,t_graph *g,int ePBC,matrix box,rvec x[])
{
  static int nerror_tot = 0;
  int    npbcdim;
  int    ng,nnodes,i;
  int    nW,nG,nB;		/* Number of Grey, Black, White	*/
  int    fW,fG;			/* First of each category	*/
  int    nerror=0;

  g->bScrewPBC = (ePBC == epbcSCREW);

  if (ePBC == epbcXY)
    npbcdim = 2;
  else
    npbcdim = 3;

  GCHECK(g);
  /* This puts everything in the central box, that is does not move it 
   * at all. If we return without doing this for a system without bonds
   * (i.e. only settles) all water molecules are moved to the opposite octant
   */
  for(i=0; (i<g->nnodes); i++) {
      g->ishift[i][XX]=g->ishift[i][YY]=g->ishift[i][ZZ]=0;
  }
    
  if (!g->nbound)
    return;

  nnodes=g->nnodes;
  if (nnodes > g->negc) {
    g->negc = nnodes;
    srenew(g->egc,g->negc);
  }
  memset(g->egc,0,(size_t)(nnodes*sizeof(g->egc[0])));

  nW=g->nbound;
  nG=0;
  nB=0;

  fW=0;

  /* We even have a loop invariant:
   * nW+nG+nB == g->nbound
   */
#ifdef DEBUG2
  fprintf(log,"Starting W loop\n");
#endif
  while (nW > 0) {
    /* Find the first white, this will allways be a larger
     * number than before, because no nodes are made white
     * in the loop
     */
    if ((fW=first_colour(fW,egcolWhite,g,g->egc)) == -1) 
      gmx_fatal(FARGS,"No WHITE nodes found while nW=%d\n",nW);
    
    /* Make the first white node grey */
    g->egc[fW]=egcolGrey;
    nG++;
    nW--;

    /* Initial value for the first grey */
    fG=fW;
#ifdef DEBUG2
    fprintf(log,"Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n",
	    nW,nG,nB,nW+nG+nB);
#endif
    while (nG > 0) {
      if ((fG=first_colour(fG,egcolGrey,g,g->egc)) == -1)
	gmx_fatal(FARGS,"No GREY nodes found while nG=%d\n",nG);
      
      /* Make the first grey node black */
      g->egc[fG]=egcolBlack;
      nB++;
      nG--;

      /* Make all the neighbours of this black node grey
       * and set their periodicity 
       */
      ng=mk_grey(log,nnodes,g->egc,g,&fG,npbcdim,box,x,&nerror);
      /* ng is the number of white nodes made grey */
      nG+=ng;
      nW-=ng;
    }
  }
  if (nerror > 0) {
    nerror_tot++;
    if (nerror_tot <= 100) {
      fprintf(stderr,"There were %d inconsistent shifts. Check your topology\n",
	      nerror);
      if (log) {
	fprintf(log,"There were %d inconsistent shifts. Check your topology\n",
		nerror);
      }
    }
    if (nerror_tot == 100) {
      fprintf(stderr,"Will stop reporting inconsistent shifts\n");
      if (log) {
	fprintf(log,"Will stop reporting inconsistent shifts\n");
      }
    }
  }
}

/************************************************************
 *
 *      A C T U A L   S H I F T   C O D E
 *
 ************************************************************/

void shift_x(t_graph *g,matrix box,rvec x[],rvec x_s[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  GCHECK(g);
  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  if (g->bScrewPBC) {
    for(i=0,j=g0; (i<gn); i++,j++) { 
      tx=is[i][XX];
      ty=is[i][YY];
      tz=is[i][ZZ];
      
      if ((tx > 0 && tx % 2 == 1) ||
	  (tx < 0 && -tx %2 == 1)) {
	x_s[j][XX] = x[j][XX] + tx*box[XX][XX];
	x_s[j][YY] = box[YY][YY] + box[ZZ][YY] - x[j][YY];
	x_s[j][ZZ] = box[ZZ][ZZ]               - x[j][ZZ];
      } else {
	x_s[j][XX] = x[j][XX];
      }
      x_s[j][YY] = x[j][YY] + ty*box[YY][YY] + tz*box[ZZ][YY];
      x_s[j][ZZ] = x[j][ZZ] + tz*box[ZZ][ZZ];
    }
  } else if (TRICLINIC(box)) {
     for(i=0,j=g0; (i<gn); i++,j++) { 
	 tx=is[i][XX];
	 ty=is[i][YY];
	 tz=is[i][ZZ];
	 
	 x_s[j][XX]=x[j][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
	 x_s[j][YY]=x[j][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
	 x_s[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
     }
  } else {
     for(i=0,j=g0; (i<gn); i++,j++) { 
	 tx=is[i][XX];
	 ty=is[i][YY];
	 tz=is[i][ZZ];
	 
	 x_s[j][XX]=x[j][XX]+tx*box[XX][XX];
	 x_s[j][YY]=x[j][YY]+ty*box[YY][YY];
	 x_s[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
     }
  }       
     
}

void shift_self(t_graph *g,matrix box,rvec x[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  if (g->bScrewPBC)
    gmx_incons("screw pbc not implemented for shift_self");

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;

#ifdef DEBUG
  fprintf(stderr,"Shifting atoms %d to %d\n",g0,g0+gn);
#endif
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) { 
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
	  x[j][YY]=x[j][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
	  x[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) { 
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]+tx*box[XX][XX];
	  x[j][YY]=x[j][YY]+ty*box[YY][YY];
	  x[j][ZZ]=x[j][ZZ]+tz*box[ZZ][ZZ];
      }
  }       
  
}

void unshift_x(t_graph *g,matrix box,rvec x[],rvec x_s[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  if (g->bScrewPBC)
    gmx_incons("screw pbc not implemented for unshift_x");

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x_s[j][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
	  x[j][YY]=x_s[j][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
	  x[j][ZZ]=x_s[j][ZZ]-tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x_s[j][XX]-tx*box[XX][XX];
	  x[j][YY]=x_s[j][YY]-ty*box[YY][YY];
	  x[j][ZZ]=x_s[j][ZZ]-tz*box[ZZ][ZZ];
      }
  }
}

void unshift_self(t_graph *g,matrix box,rvec x[])
{
  ivec *is;
  int g0,gn;
  int i,j,tx,ty,tz;

  if (g->bScrewPBC)
    gmx_incons("screw pbc not implemented for unshift_self");

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  if(TRICLINIC(box)) {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
	  x[j][YY]=x[j][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
	  x[j][ZZ]=x[j][ZZ]-tz*box[ZZ][ZZ];
      }
  } else {
      for(i=0,j=g0; (i<gn); i++,j++) {
	  tx=is[i][XX];
	  ty=is[i][YY];
	  tz=is[i][ZZ];
	  
	  x[j][XX]=x[j][XX]-tx*box[XX][XX];
	  x[j][YY]=x[j][YY]-ty*box[YY][YY];
	  x[j][ZZ]=x[j][ZZ]-tz*box[ZZ][ZZ];
      }
  }
}
#undef GCHECK
