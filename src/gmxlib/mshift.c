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
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_mshift_c = "$Id$";

#include <string.h>
#include "assert.h"
#include "bondf.h"
#include "smalloc.h"
#include "fatal.h"
#include "macros.h"
#include "vec.h"
#include "statusio.h"
#include "statutil.h"
#include "futil.h"
#include "copyrite.h"
#include "mshift.h"
#include "main.h"

/************************************************************
 *
 *      S H I F T   U T I L I T I E S
 *
 ************************************************************/
 
static ivec iv_shift[N_IVEC];

static void mk_ivshift(void)
{
  int i,j,k,n;
  
  n=0;
  for(i=-D_BOX; (i<=D_BOX); i++)
    for(j=-D_BOX; (j<=D_BOX); j++)
      for(k=-D_BOX; (k<=D_BOX); k++,n++) {
	iv_shift[n][XX]=i;
	iv_shift[n][YY]=j;
	iv_shift[n][ZZ]=k;
      }
  assert(n==N_IVEC);
}

/************************************************************
 *
 *      G R A P H   G E N E R A T I O N    C O D E
 *
 ************************************************************/
 
typedef enum { egcolWhite, egcolGrey, egcolBlack, egcolNR } egCol;

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
	  if (g->nedge[inda] == g->maxbond)
	    fatal_error(0,"More than %d bonds per atom (atom %d)\n",
			g->maxbond,aa+1);
	  g->edge[inda][g->nedge[inda]++]=ia[k];
	}
      }
  }
}

static void mk_igraph(t_graph *g,t_functype ftype[],t_ilist *il,bool bAll)
{
  t_iatom *ia;
  t_iatom settle[3];
  t_iatom tp;
  int     i,j,np;
  int     end;

  end=il->nr;
  ia=il->iatoms;
  for(i=0; (i < end); ) {
    tp=ftype[ia[0]];
    np=interaction_function[tp].nratoms;

        
    switch (tp) {
    case F_SETTLE:
      if (!bAll) {
	settle[0]=ia[1];
	settle[1]=settle[0]+1;
	settle[2]=settle[0]+2;
	add_gbond(g,settle,3);
      }
      break;
      /*case F_BONDS:
	case F_MORSE:
	case F_SHAKE:
	add_gbond(g,&(ia[1]),np);
	break;*/
    default:
      if (interaction_function[tp].flags & (IF_BOND | IF_SHAKE)) {
	if (bAll) {
	  add_gbond(g,&(ia[1]),np);
	}
	else {
	  /* Check whether all atoms are bonded now! */
	  for(j=0; (j<np); j++)
	    if (g->nedge[ia[1+j]-g->start] == 0)
	      break;
	  if (j < np)
	    add_gbond(g,&(ia[1]),np);
	}
      }
      break;
    }
    ia+=np+1;
    i+=np+1;
  }
}

void p_graph(FILE *log,char *title,t_graph *g)
{
  int i,j;

  fprintf(log,"graph:  %s\n",title);
  fprintf(log,"maxbond:%d\n",g->maxbond);
  fprintf(log,"nnodes: %d\n",g->nnodes);
  fprintf(log,"nbound: %d\n",g->nbound);
  fprintf(log,"start:  %d\n",g->start);
  fprintf(log,"end:    %d\n",g->end);
  fprintf(log," atom shift nedg   e1   e2 etc.\n");
  for(i=0; (i<g->nnodes); i++)
    if (g->nedge[i] > 0) {
      fprintf(log,"%5d%6d%5d",g->start+i+1,g->ishift[i],g->nedge[i]);
      for(j=0; (j<g->nedge[i]); j++)
	fprintf(log,"%5d",g->edge[i][j]+1);
      fprintf(log,"\n");
    }
  fflush(log);
}

static void calc_1se(t_graph *g,t_ilist *il,t_functype ftype[],
		     short nbond[])
{
  int     k,nratoms,tp,end,j;
  t_iatom *ia,iaa;

  end=il->nr;

  ia=il->iatoms;
  for(j=0; (j<end); j+=nratoms+1,ia+=nratoms+1) {
    tp      = ftype[ia[0]];
    nratoms = interaction_function[tp].nratoms;
    
    if (tp == F_SETTLE) {
      iaa          = ia[1];
      nbond[iaa]   = 2;
      nbond[iaa+1] = 2;
      nbond[iaa+2] = 2;
      g->start     = min(g->start,iaa);
      g->end       = max(g->end,iaa+2);
    }
    else {
      for(k=0; (k<nratoms); k++) {
	iaa=ia[k+1];
	g->start=min(g->start,iaa);
	g->end  =max(g->end,  iaa);
	
	switch (tp) {
	case F_BONDS:
	case F_MORSE:
	case F_SHAKE:
	  nbond[iaa]++;
	  break;
	default:
	  break;
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
    if (interaction_function[i].flags & (IF_BOND | IF_SHAKE))
      calc_1se(g,&idef->il[i],idef->functype,nbond);
  }
  
  nnb=0;
  for(i=g->start; (i<=g->end); i++)
    nnb=max(nnb,nbond[i]);
  if (stdlog)
    fprintf(stdlog,"Max number of bonds per atom is %d\n",nnb);
  
  sfree(nbond);
  
  g->maxbond=nnb+6;
}

t_graph *mk_graph(t_idef *idef,int natoms,bool bShakeOnly)
{
  t_graph *g;
  int     i;
  mk_ivshift();
  
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
    snew(g->edge,g->nnodes);
    snew(g->edge[0],g->maxbond*g->nnodes);

    for(i=1; (i<g->nnodes); i++)
      g->edge[i]=g->edge[i-1]+g->maxbond;

    if (!bShakeOnly) {
      for(i=0; (i<F_NRE); i++)
	if (interaction_function[i].flags & (IF_CONNECT))
	  mk_igraph(g,idef->functype,&(idef->il[i]),TRUE);
      for(i=0; (i<F_NRE); i++)
	if (interaction_function[i].flags & (IF_BOND | ~IF_CONNECT))
	  mk_igraph(g,idef->functype,&(idef->il[i]),FALSE);
      
      /*mk_igraph(g,idef->functype,&(idef->il[F_MORSE]));
	mk_igraph(g,idef->functype,&(idef->il[F_SETTLE]));
	mk_igraph(g,idef->functype,&(idef->il[F_LJ14]));*/
    }
    else
      mk_igraph(g,idef->functype,&(idef->il[F_SHAKE]),TRUE);

    g->nbound=0;
    for(i=0; (i<g->nnodes); i++)
      if (g->nedge[i] > 0)
        g->nbound++;
  }
#ifdef DEBUG
  p_graph(stdlog,"graph",g);
#endif
  return g;
}

void done_graph(t_graph *g)
{
  if (g->nnodes > 0) {
    sfree(g->ishift);
    sfree(g->nedge);
    /* This is malloced in a NASTY way, see above */
    sfree(g->edge[0]);
    sfree(g->edge);
  }
}

/************************************************************
 *
 *      S H I F T   C A L C U L A T I O N   C O D E
 *
 ************************************************************/

static int mk_1shift(rvec hbox,rvec xi,t_ishift mi,rvec xj)
{
  /* Calculate periodicity for rectangular box... */
  int  m;
  ivec mj;
  rvec dx;
  
  rvec_sub(xi,xj,dx);

  for(m=0; (m<DIM); m++) {
    /* If dx < hbox, then xj will be reduced by box, so that
     * xi - xj will be bigger
     */
    if (dx[m] < -hbox[m])
      mj[m]=iv_shift[mi][m]-1;
    else if (dx[m] >= hbox[m])
      mj[m]=iv_shift[mi][m]+1;
    else
      mj[m]=iv_shift[mi][m];
  }
  return IVEC2IS(mj);
}

static int mk_grey(FILE *log,int nnodes,egCol egc[],t_graph *g,int *AtomI,
		   matrix box,rvec x[],int *nerror)
{
  int      m,j,ng,ai,aj,g0,is_aj;
  rvec     hbox;
  
  for(m=0; (m<DIM); m++)
    hbox[m]=box[m][m]*0.5;

  ng=0;
  ai=*AtomI;
  
  g0=g->start;
  /* Loop over all the bonds */
  for(j=0; (j<g->nedge[ai]); j++) {
    aj=g->edge[ai][j]-g0;
    /* If there is a white one, make it gray and set pbc */
    is_aj=mk_1shift(hbox,x[g0+ai],g->ishift[ai],x[g0+aj]);
    if ((is_aj >= N_IVEC) || (is_aj < 0))
      fatal_error(0,"is_aj out of range (%d)",is_aj);
    
    if (egc[aj] == egcolWhite) {
      if (aj < *AtomI)
	*AtomI = aj;
      egc[aj] = egcolGrey;
      
      g->ishift[aj]=is_aj;

      ng++;
    }
    else if (is_aj != g->ishift[aj]) {
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
  static egCol *egc=NULL;	/* The colour of each node	*/

  if (!g->nbound)
    return;

  for(i=0; (i<g->nnodes); i++) {
    g->ishift[i]=CENTRAL;
  }
    
  nnodes=g->nnodes;
  if (egc == NULL)
    snew(egc,nnodes);
  else
    memset(egc,0,nnodes*sizeof(egc[0]));
  
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
    if ((fW=first_colour(fW,egcolWhite,g,egc)) == -1) 
      fatal_error(0,"No WHITE nodes found while nW=%d\n",nW);
    
    /* Make the first white node grey */
    egc[fW]=egcolGrey;
    nG++;
    nW--;

    /* Initial value for the first grey */
    fG=fW;
#ifdef DEBUG2
    fprintf(log,"Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n",
	    nW,nG,nB,nW+nG+nB);
#endif
    while (nG > 0) {
      if ((fG=first_colour(fG,egcolGrey,g,egc)) == -1)
	fatal_error(0,"No GREY nodes found while nG=%d\n",nG);
      
      /* Make the first grey node black */
      egc[fG]=egcolBlack;
      nB++;
      nG--;

      /* Make all the neighbours of this black node grey
       * and set their periodicity 
       */
      ng=mk_grey(log,nnodes,egc,g,&fG,box,x,&nerror);
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
 
void shift_atom(t_graph *g,rvec sv[],rvec x[],rvec x_s[],atom_id ai)
{
  int t;
  
  t=g->ishift[ai-g->start];

  x_s[ai][XX]=x[ai][XX]+sv[t][XX];
  x_s[ai][YY]=x[ai][YY]+sv[t][YY];
  x_s[ai][ZZ]=x[ai][ZZ]+sv[t][ZZ];
}
 
void unshift_atom(t_graph *g,rvec sv[],rvec x[],rvec x_s[],atom_id ai)
{
  int t;
  
  t=g->ishift[ai-g->start];
  
  x_s[ai][XX]=x[ai][XX]-sv[t][XX];
  x_s[ai][YY]=x[ai][YY]-sv[t][YY];
  x_s[ai][ZZ]=x[ai][ZZ]-sv[t][ZZ];
}

void shift_x(t_graph *g,rvec sv[],rvec x[],rvec x_s[])
{
  t_ishift *is;
  int      g0,gn;
  int      i,j,isd;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  for(i=0,j=g0; (i<gn); i++,j++) { 
    isd=is[i];
#ifdef DEBUG
    fprintf(stdlog,"x[%d]=(%e,%e,%e), sv[is[%d]]=(%e,%e,%e)\n",
	    j,x[j][XX],x[j][YY],x[j][ZZ],
	    i,sv[isd][XX],sv[isd][YY],sv[isd][ZZ]);
#endif
    rvec_add(x[j],sv[isd],x_s[j]);
  }
}

void shift_self(t_graph *g,rvec sv[],rvec x[])
{
  t_ishift *is;
  int      g0,gn;
  int      i,j,isd;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;

#ifdef DEBUG
  fprintf(stderr,"Shifting atoms %d to %d\n",g0,g0+gn);
#endif
  for(i=0,j=g0; (i<gn); i++,j++) { 
    isd=is[i];
#ifdef DEBUG
    fprintf(stdlog,"x[%d]=(%e,%e,%e), sv[is[%d]]=(%e,%e,%e)\n",
	    j,x[j][XX],x[j][YY],x[j][ZZ],
	    i,sv[isd][XX],sv[isd][YY],sv[isd][ZZ]);
#endif
    rvec_inc(x[j],sv[isd]);
  }
}

void unshift_x(t_graph *g,rvec sv[],rvec x[],rvec x_s[])
{
  t_ishift *is;
  int      g0,gn;
  int      i,j,isd;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  for(i=0,j=g0; (i<gn); i++,j++) {
    isd=is[i];
#ifdef DEBUG
    fprintf(stdlog,"x_s[%d]=(%e,%e,%e), sv[is[%d]]=(%e,%e,%e)\n",
	    j,x_s[j][XX],x_s[j][YY],x_s[j][ZZ],
	    i,sv[isd][XX],sv[isd][YY],sv[isd][ZZ]);
#endif
    rvec_sub(x_s[j],sv[isd],x[j]);
  }
}

void unshift_self(t_graph *g,rvec sv[],rvec x[])
{
  t_ishift *is;
  int      g0,gn;
  int      i,j,isd;

  g0=g->start;
  gn=g->nnodes;
  is=g->ishift;
  
  for(i=0,j=g0; (i<gn); i++,j++) {
    isd=is[i];
#ifdef DEBUG
    fprintf(stdlog,"x_s[%d]=(%e,%e,%e), sv[is[%d]]=(%e,%e,%e)\n",
	    j,x_s[j][XX],x_s[j][YY],x_s[j][ZZ],
	    i,sv[isd][XX],sv[isd][YY],sv[isd][ZZ]);
#endif
    rvec_dec(x[j],sv[isd]);
  }
}

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
  g=mk_graph(&(top.idef),top.atoms.nr,pid);

  out=ffopen(targ.outfile,"w");

  fprintf(stderr,"Making Shift...\n");
  mk_mshift(out,g,box,x,mshift);

  p_graph(out,"In Den Haag daar woont een graaf...",g);
  fclose(out);
}
#endif

