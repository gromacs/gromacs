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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "smalloc.h"
#include "fatal.h"
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
 
static void add_gbond(t_graph *g,t_iatom ia[],int np)
{
  int     j,k,l;
  atom_id aa,inda;

  for(j=0; (j<np); j++) {
    aa=ia[j];
    inda=aa-g->start;
    for(k=0; (k<np); k++)
      if (j != k) {
	for(l=0; (l<g->nedge[inda]); l++)
	  if (g->edge[inda][l] == ia[k])
	    break;
	if (l == g->nedge[inda]) {
	  if (g->nedge[inda] == g->maxedge)
	    gmx_fatal(FARGS,"More than %d graph edges per atom (atom %d)\n",
			g->maxedge,aa+1);
	  g->edge[inda][g->nedge[inda]++]=ia[k];
	}
      }
  }
}

static void mk_igraph(t_graph *g,t_functype ftype[],t_ilist *il,
		      int natoms,bool bAll)
{
  t_iatom *ia,waterh[3],*iap;
  t_iatom tp;
  int     i,j,np,nbonded;
  int     end;

  end=il->nr;
  ia=il->iatoms;

  i = 0;
  while (i < end) {
    tp=ftype[ia[0]];
    np=interaction_function[tp].nratoms;
    
    if ((ia[1] < natoms) && (interaction_function[tp].flags & IF_GRAPH)) {
      if (ia[np] >= natoms)
	gmx_fatal(FARGS,"Molecule in topology has atom numbers below and "
		    "above natoms (%d).\n"
		    "You are probably trying to use a trajectory which does "
		    "not match the first %d atoms of the run input file.\n"
		    "You can make a matching run input file with tpbconv.",
		    natoms,natoms);
      if (tp == F_SETTLE) {
	/* Bond all the atoms in the settle */
	nbonded = 3;
	waterh[0] = ia[1];
	waterh[1] = ia[1]+1;
	waterh[2] = ia[1]+2;
	iap = waterh;
      } else {
	if (interaction_function[tp].flags & IF_VSITE)
	  /* Bond a virtual site only to the first constructing atom */
	  nbonded = 2;
	else
	  nbonded = np;
	iap = &(ia[1]);
      }
      if (bAll)
	add_gbond(g,iap,nbonded);
      else {
	/* Check whether all atoms are bonded now! */
	for(j=0; j<nbonded; j++)
	  if (g->nedge[iap[j]-g->start] == 0)
	    break;
	if (j < nbonded)
	  add_gbond(g,iap,nbonded);
      }
    }
    ia+=np+1;
    i+=np+1;
  }
}

static void g_error(int line,char *file)
{
  gmx_fatal(FARGS,"Tring to print non existant graph (file %s, line %d)",
	      file,line);
}
#define GCHECK(g) if (g == NULL) g_error(__LINE__,__FILE__)

void p_graph(FILE *log,char *title,t_graph *g)
{
  int i,j;
  
  GCHECK(g);
  fprintf(log,"graph:  %s\n",title);
  fprintf(log,"maxedge:%d\n",g->maxedge);
  fprintf(log,"nnodes: %d\n",g->nnodes);
  fprintf(log,"nbound: %d\n",g->nbound);
  fprintf(log,"start:  %d\n",g->start);
  fprintf(log,"end:    %d\n",g->end);
  fprintf(log," atom shiftx shifty shiftz nedg   e1   e2 etc.\n");
  for(i=0; (i<g->nnodes); i++)
    if (g->nedge[i] > 0) {
	fprintf(log,"%5d%7d%7d%7d%5d",g->start+i+1,g->ishift[i][XX],g->ishift[i][YY],
		g->ishift[i][ZZ],g->nedge[i]);
      for(j=0; (j<g->nedge[i]); j++)
	fprintf(log,"%5u",g->edge[i][j]+1);
      fprintf(log,"\n");
    }
  fflush(log);
}

static void calc_1se(t_graph *g,t_ilist *il,t_functype ftype[],
		     short nbond[],int natoms)
{
  int     k,nratoms,nbonded,tp,end,j;
  t_iatom *ia,iaa;

  end=il->nr;

  ia=il->iatoms;
  for(j=0; (j<end); j+=nratoms+1,ia+=nratoms+1) {
    tp      = ftype[ia[0]];
    nratoms = interaction_function[tp].nratoms;
    
    if (tp == F_SETTLE) {
      iaa          = ia[1];
      if (iaa<natoms) {
	nbond[iaa]   = 2;
	nbond[iaa+1] = 2;
	nbond[iaa+2] = 2;
	g->start     = min(g->start,iaa);
	g->end       = max(g->end,iaa+2);
      }
    } else {
      if (interaction_function[tp].flags & IF_VSITE)
	/* Bond a virtual site only to the first constructing atom */
	nbonded = 2;
      else
	nbonded = nratoms;
      for(k=0; (k<nbonded); k++) {
	iaa=ia[k+1];
	if (iaa<natoms) {
	  g->start=min(g->start,iaa);
	  g->end  =max(g->end,  iaa);
	  if (interaction_function[tp].flags & IF_GRAPH)
	    nbond[iaa]++;
	}
      }
    }
  }
}

static void calc_start_end(t_graph *g,t_idef *idef,int natoms)
{
  short *nbond;
  int   i,nnb;
  
  g->start=natoms;
  g->end=0;

  snew(nbond,natoms);
  for(i=0; (i<F_NRE); i++) {
    if (interaction_function[i].flags & IF_GRAPH)
      calc_1se(g,&idef->il[i],idef->functype,nbond,natoms);
  }
  
  nnb=0;
  for(i=g->start; (i<=g->end); i++)
    nnb=max(nnb,nbond[i]);
  if (stdlog)
    fprintf(stdlog,"Max number of graph edges per atom is %d\n",nnb);
  
  sfree(nbond);
  
  g->maxedge=nnb+6;
}

t_graph *mk_graph(t_idef *idef,int natoms,bool bShakeOnly,bool bSettle)
{
  t_graph *g;
  int     i;
  
  snew(g,1);

  calc_start_end(g,idef,natoms);
  
  if (g->start >= g->end) {
    g->nnodes=0;
  }
  else {
    g->nnodes=g->end-g->start+1;
    snew(g->ishift,g->nnodes);
    snew(g->nedge,g->nnodes);
  
    /* To prevent malloc problems with many small arrays, using realloc,
     * we allocate some more memory, and divide it ourselves.
     * We calculate pointers... (Yuck Yuck)
     */
    if (debug)
      fprintf(debug,"MSHIFT: nnodes=%d, maxedge=%d\n",g->nnodes,g->maxedge);
    snew(g->edge,g->nnodes);
    snew(g->edge[0],g->maxedge*g->nnodes);

    for(i=1; (i<g->nnodes); i++)
      g->edge[i]=g->edge[i-1]+g->maxedge;

    if (!bShakeOnly) {
      /* First add all the real bonds: they should determine the molecular 
       * graph.
       */
      for(i=0; (i<F_NRE); i++)
	if (interaction_function[i].flags & IF_GRAPH)
	  mk_igraph(g,idef->functype,&(idef->il[i]),natoms,TRUE);
      /* Then add all the other interactions in fixed lists, but first
       * check to see what's there already.
       */
      for(i=0; (i<F_NRE); i++) {
	if (!(interaction_function[i].flags & IF_GRAPH)) {
	  mk_igraph(g,idef->functype,&(idef->il[i]),natoms,FALSE);
	}
      }
    }
    else {
      /* This is a special thing used in grompp to generate shake-blocks */
      mk_igraph(g,idef->functype,&(idef->il[F_SHAKE]),natoms,TRUE);
      if (bSettle)
	mk_igraph(g,idef->functype,&(idef->il[F_SETTLE]),natoms,TRUE);
    }
    g->nbound=0;
    for(i=0; (i<g->nnodes); i++)
      if (g->nedge[i] > 0)
        g->nbound++;
  }
#ifdef DEBUG
  p_graph(stdlog,"graph",g);
#endif

  g->negc = 0;
  g->egc = NULL;
  
  return g;
}

void done_graph(t_graph *g)
{
  GCHECK(g);
  if (g->nnodes > 0) {
    sfree(g->ishift);
    sfree(g->nedge);
    /* This is malloced in a NASTY way, see above */
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

static void mk_1shift_tric(matrix box,rvec hbox,rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for triclinic box... */
  int  m,d;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  for(m=DIM-1; (m>=0); m--) {
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

static void mk_1shift(rvec hbox,rvec xi,rvec xj,int *mi,int *mj)
{
  /* Calculate periodicity for rectangular box... */
  int  m;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  for(m=0; (m<DIM); m++) {
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

static int mk_grey(FILE *log,int nnodes,egCol egc[],t_graph *g,int *AtomI,
		   matrix box,rvec x[],int *nerror)
{
  int      m,j,ng,ai,aj,g0;
  rvec     hbox;
  bool     bTriclinic;
  ivec     is_aj;
  
  for(m=0; (m<DIM); m++)
    hbox[m]=box[m][m]*0.5;
  bTriclinic = TRICLINIC(box);
  
  ng=0;
  ai=*AtomI;
  
  g0=g->start;
  /* Loop over all the bonds */
  for(j=0; (j<g->nedge[ai]); j++) {
    aj=g->edge[ai][j]-g0;
    /* If there is a white one, make it gray and set pbc */
    if (bTriclinic)
      mk_1shift_tric(box,hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    else
      mk_1shift(hbox,x[g0+ai],x[g0+aj],g->ishift[ai],is_aj);
    
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

void mk_mshift(FILE *log,t_graph *g,matrix box,rvec x[])
{
  int    ng,nnodes,i;
  int    nW,nG,nB;		/* Number of Grey, Black, White	*/
  int    fW,fG;			/* First of each category	*/
  int    nerror=0;

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
      ng=mk_grey(log,nnodes,g->egc,g,&fG,box,x,&nerror);
      /* ng is the number of white nodes made grey */
      nG+=ng;
      nW-=ng;
    }
  }
  if (nerror > 0)
    fprintf(log,"There were %d inconsistent shifts. Check your topology\n",
	    nerror);
}

/************************************************************
 *
 *      A C T U A L   S H I F T   C O D E
 *
 ************************************************************/
 
static void shift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai)
{
  int tx,ty,tz;
  
  tx=(g->ishift[ai-g->start])[XX];
  ty=(g->ishift[ai-g->start])[YY];
  tz=(g->ishift[ai-g->start])[ZZ];

  x_s[ai][XX]=x[ai][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
  x_s[ai][YY]=x[ai][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
  x_s[ai][ZZ]=x[ai][ZZ]+tz*box[ZZ][ZZ];
}
 
static void unshift_atom(t_graph *g,matrix box,rvec x[],rvec x_s[],atom_id ai)
{
  int tx,ty,tz;
  
  tx=(g->ishift[ai-g->start])[XX];
  ty=(g->ishift[ai-g->start])[YY];
  tz=(g->ishift[ai-g->start])[ZZ];

  x_s[ai][XX]=x[ai][XX]-tx*box[XX][XX]-ty*box[YY][XX]-tz*box[ZZ][XX];
  x_s[ai][YY]=x[ai][YY]-ty*box[YY][YY]-tz*box[ZZ][YY];
  x_s[ai][ZZ]=x[ai][ZZ]-tz*box[ZZ][ZZ];
}

void shift_x(t_graph *g,matrix box,rvec x[],rvec x_s[])
{
  ivec *is;
  int      g0,gn;
  int      i,j,tx,ty,tz;

  GCHECK(g);
  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  if(TRICLINIC(box)) {
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

#ifdef DEBUGMSHIFT
void main(int argc,char *argv[])
{
  FILE         *out;
  t_args       targ;
  t_topology   top;
  t_statheader sh;
  rvec         *x;
  ivec         *mshift;
  matrix       box;

  t_graph      *g;
  int          i,idum,pid;
  real         rdum;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,&targ,PCA_NEED_INOUT,NULL);
  if (argc > 1)
    pid=atoi(argv[1]);
  else
    pid=0;
  
  read_status_header(targ.infile,&sh);
  snew(x,sh.natoms);
  snew(mshift,sh.natoms);

  fprintf(stderr,"Reading Status %s\n",targ.infile);
  read_status(targ.infile,&idum,&rdum,&rdum,NULL,
	      box,NULL,NULL,&idum,x,NULL,NULL,&idum,NULL,&top);

  fprintf(stderr,"Making Graph Structure...\n");
  g=mk_graph(&(top.idef),top.atoms.nr,FALSE,FALSE);

  out=ffopen(targ.outfile,"w");

  fprintf(stderr,"Making Shift...\n");
  mk_mshift(out,g,box,x,mshift);

  p_graph(out,"In Den Haag daar woont een graaf...",g);
  fclose(out);
}
#endif

