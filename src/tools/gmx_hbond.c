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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>

#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "futil.h"
#include "tpxio.h"
#include "physics.h"
#include "macros.h"
#include "fatal.h"
#include "index.h"
#include "smalloc.h"
#include "assert.h"
#include "vec.h"
#include "xvgr.h"
#include "gstat.h"
#include "matio.h"
#include "string2.h"
#include "pbc.h"

#define max_hx 7
typedef int t_hx[max_hx];
#define NRHXTYPES max_hx
char *hxtypenames[NRHXTYPES]=
{"n-n","n-n+1","n-n+2","n-n+3","n-n+4","n-n+5","n-n>6"};

enum { gr0, gr1, grI, grNR };
#define OGRP (gr1-grp)
  
static char *grpnames[grNR] = {"0","1","I" };

#define HB_NO 0
#define HB_YES 1<<0
#define HB_INS 1<<1
#define HB_YESINS HB_YES|HB_INS
#define HB_NR (1<<2)
#define MAXHYDRO 4

typedef struct {
  int nr;
  int maxnr;
  atom_id *atoms;
} t_ncell;

typedef struct {
  t_ncell d[grNR];
  t_ncell a[grNR];
} t_gridcell;

typedef int     t_icell[grNR];
typedef atom_id h_id[MAXHYDRO];
 
typedef struct {
  bool     bExisted[MAXHYDRO]; /* has this hbond existed ever? */
  bool     bInsert;  /* has this hbond been invaded? */
  /* Bitmask array which tells whether a hbond is present
   * at a given time. Either of these may be NULL 
   */
  unsigned int **exist; 
} t_hbond;

typedef struct {
  int     nra,max_nra;
  atom_id *acc;             /* Atom numbers of the acceptors     */
  int     *aptr;            /* Map atom number to acceptor index */
} t_acceptors;

typedef struct {
  int      nrd,max_nrd;
  int      *d;                 /* Atom numbers of the donors         */
  int      *nhydro;            /* Number of hydrogens for each donor */
  h_id     *hydro;             /* The atom numbers of the hydrogens  */
  int      *dptr;              /* Map atom number to donor index     */
} t_donors;

typedef struct {
  bool        bHBmap,bDAnr;
  int         wordlen;
  /* The following arrays are nframes long */
  int         nframes,max_frames,maxhydro;
  int         *nhb;
  real        *time;
  t_icell     *danr;
  t_hx        *nhx;
  /* These structures are initialized from the topology at start up */
  t_donors    *(d[grNR]);
  t_acceptors *(a[grNR]);
  /* This holds a matrix with all possible hydrogen bonds */
  int         nd,na,nrhb;
  t_hbond     **hbmap;
} t_hbdata;

static t_hbdata *mk_hbdata(bool bHBmap,bool bDAnr,bool bMerge)
{
  t_hbdata *hb;
  
  snew(hb,1);
  hb->wordlen = 8*sizeof(unsigned int);
  hb->bHBmap  = bHBmap;
  hb->bDAnr   = bDAnr;
  if (bMerge)
    hb->maxhydro = 1;
  else
    hb->maxhydro = MAXHYDRO;
    

  return hb;
}

static void mk_hbmap(t_hbdata *hb,bool bTwo,bool bInsert)
{
  int i,j;
  
  if (bInsert) {
    hb->nd = hb->d[grI]->nrd;
    hb->na = hb->a[grI]->nra;
  }
  else {
    hb->nd = hb->d[gr0]->nrd;
    hb->na = hb->a[gr0]->nra;
    if (bTwo) {
      hb->nd += hb->d[gr1]->nrd;
      hb->na += hb->a[gr1]->nra;
    }
  }
  fprintf(stderr,"Going to allocate %d kb of memory,  and that's only the beginning\n",(int)((hb->nd*hb->na*sizeof(hb->hbmap[0][0]))/1024.0));
  
  snew(hb->hbmap,hb->nd);
  for(i=0; (i<hb->nd); i++) {
    snew(hb->hbmap[i],hb->na);
    for(j=0; (j<hb->na); j++) 
      snew(hb->hbmap[i][j].exist,hb->maxhydro);
  }
}

static void add_frames(t_hbdata *hb,int nframes)
{
  int i,j,k,l;
  
  if (nframes >= hb->max_frames) {
    hb->max_frames += 12000;
    srenew(hb->time,hb->max_frames);
    srenew(hb->nhb,hb->max_frames);
    srenew(hb->nhx,hb->max_frames);
    if (hb->bDAnr)
      srenew(hb->danr,hb->max_frames);
    if (hb->bHBmap) {
      for(i=0; (i<hb->nd); i++) 
	for(j=0; (j<hb->na); j++) 
	  for(k=0; (k<hb->maxhydro); k++)
	    if (hb->hbmap[i][j].exist[k]) {
	      srenew(hb->hbmap[i][j].exist[k],hb->max_frames/hb->wordlen);
	      for(l=hb->nframes/hb->wordlen; (l<hb->max_frames/hb->wordlen); l++)
		hb->hbmap[i][j].exist[k][l] = 0;
	    }
    }
  }
  hb->nframes=nframes;
}

#define OFFSET(frame) (frame / 32)
#define MASK(frame)   (1 << (frame % 32))

static void set_hb(unsigned int hbexist[],unsigned int frame,bool bValue)
{
  if (bValue)
    hbexist[OFFSET(frame)] |= MASK(frame);
  else
    hbexist[OFFSET(frame)] &= ~MASK(frame);
}

static bool is_hb(unsigned int hbexist[],int frame)
{
  return ((hbexist[OFFSET(frame)] & MASK(frame)) != 0) ? 1 : 0;
}

static void add_hbond(t_hbdata *hb,int d,int a,int h,int grpd,int grpa,
		      int frame,bool bInsert,bool bMerge)
{ 
  int k,id,ia;
  
  if ((id = hb->d[grpd]->dptr[d]) == NOTSET)
    fatal_error(0,"No donor atom %d (%s %d)",d+1,__FILE__,__LINE__);
  if ((ia = hb->a[grpa]->aptr[a]) == NOTSET)
    fatal_error(0,"No acceptor atom %d (%s %d)",a+1,__FILE__,__LINE__);
  
  /* Loop over hydrogens */
  for(k=0; (k<hb->d[grpd]->nhydro[id]); k++) 
    if (hb->d[grpd]->hydro[id][k] == h)
      break;
  if (k == hb->d[grpd]->nhydro[id])
    fatal_error(0,"Donor %d does not have hydrogen %d  %s %d",d+1,h+1,
		__FILE__,__LINE__);;
    
  if (bMerge)
    k = 0;
  if (hb->bHBmap) {
    if (!hb->hbmap[id][ia].exist[k])
      snew(hb->hbmap[id][ia].exist[k],hb->max_frames/hb->wordlen);
    if (frame >= 0)
      set_hb(hb->hbmap[id][ia].exist[k],frame,TRUE);
  }
  hb->hbmap[id][ia].bInsert = bInsert || hb->hbmap[id][ia].bInsert;

  if (frame >= 0) {
    hb->nhb[frame]++;
    if (!hb->hbmap[id][ia].bExisted[k]) {
      hb->hbmap[id][ia].bExisted[k] = TRUE;
      hb->nrhb++;
    }
  }
}

static bool in_list(atom_id selection,int isize,atom_id *index)
{
  int i;
  bool bFound;
  
  bFound=FALSE;
  for(i=0; (i<isize) && !bFound; i++)
    if(selection == index[i])
      bFound=TRUE;
  
  return bFound;
}

static char *mkatomname(t_atoms *atoms,int i)
{
  static char buf[32];
  int rnr;
  
  rnr = atoms->atom[i].resnr;
  sprintf(buf,"%4s%d%-4s",*atoms->resname[rnr],rnr+1,*atoms->atomname[i]);
  
  return buf;
}

static void add_acc(t_acceptors *a,int ia)
{
  if (a->nra >= a->max_nra) {
    a->max_nra += 16;
    srenew(a->acc,a->max_nra);
  }
  a->acc[a->nra++] = ia;
}

static t_acceptors *search_acceptors(t_topology *top,int isize, 
				     atom_id *index,bool bNitAcc,
				     bool bContact,bool bDoIt)
{
  t_acceptors *a;
  int i;
  
  snew(a,1);
  for (i=0; (i<top->atoms.nr); i++) {
    if (bDoIt) {
      if ((bContact ||
	   (((*top->atoms.atomname[i])[0] == 'O') || 
	    (bNitAcc && ((*top->atoms.atomname[i])[0] == 'N')))) &&
	  in_list(i,isize,index)) {
	add_acc(a,i);
      }
    }
  }
  snew(a->aptr,top->atoms.nr);
  for(i=0; (i<top->atoms.nr); i++)
    a->aptr[i] = NOTSET;
  for(i=0; (i<a->nra); i++)
    a->aptr[a->acc[i]] = i;
    
  return a;
}

static int acceptor_index(t_acceptors *a,atom_id i)
{
  return a->aptr[i];
}

static int donor_index(t_donors *d,atom_id i)
{
  return d->dptr[i];
}

static void add_h2d(int id,int ih,t_donors *ddd)
{
  int i;
  
  for(i=0; (i<ddd->nhydro[id]); i++) 
    if (ddd->hydro[id][i] == ih) {
      fprintf(stderr,"Hm. This isn't first time I find this donor (%d,%d)\n",
	      ddd->d[id],ih);
      break;
    }
  if (i == ddd->nhydro[id]) {
    if (ddd->nhydro[id] >= MAXHYDRO)
      fatal_error(0,"Donor %d has more than %d hydrogens! %s %d",
		  ddd->d[id],MAXHYDRO,__FILE__,__LINE__);
    ddd->hydro[id][i] = ih;
    ddd->nhydro[id]++;
  }
}
  
static void add_dh(t_donors *ddd,int id,int ih)
{
  int i;
  
  for(i=0; (i<ddd->nrd); i++) 
    if (ddd->d[i] == id) {
      add_h2d(i,ih,ddd);
      break;
    }
  if (i == ddd->nrd) {
    if (ddd->nrd >= ddd->max_nrd) {
      ddd->max_nrd += 128;
      srenew(ddd->d,ddd->max_nrd);
      srenew(ddd->nhydro,ddd->max_nrd);
      srenew(ddd->hydro,ddd->max_nrd);
    }
    ddd->d[ddd->nrd] = id;
    ddd->nhydro[ddd->nrd] = 0;
    add_h2d(ddd->nrd,ih,ddd);
    ddd->nrd++;
  }
}

static t_donors *search_donors(t_topology *top, int isize, atom_id *index,
			       bool bContact,bool bDoIt)
{
  int        i,j,nra;
  t_functype func_type;
  t_ilist    *interaction;
  atom_id    nr1,nr2;
  bool       stop;
  t_donors   *ddd;
  
  snew(ddd,1);
  
  if (bContact) {
    if (bDoIt)
      for(i=0; (i<isize); i++)
	add_dh(ddd,index[i],-1);
  }
  else {
    for(func_type=0; (func_type < F_NRE); func_type++) {
      interaction=&(top->idef.il[func_type]);
      for(i=0; i < interaction->nr; 
	  i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1) {
	/* next function */
	if (func_type != top->idef.functype[interaction->iatoms[i]]) {
	  fatal_error(0,"Error in %s,%d func_type %s",__FILE__,__LINE__,
		      interaction_function[func_type].longname);
	}
	
	/* check out this functype */
	if (func_type == F_SETTLE) {
	  nr1=interaction->iatoms[i+1];
	  
	  if (in_list(nr1,  isize,index)) {
	    if (in_list(nr1+1,isize,index))
	      add_dh(ddd,nr1,nr1+1);
	    if (in_list(nr1+2,isize,index))
	      add_dh(ddd,nr1,nr1+2);
	  }
	} 
	else if (IS_CHEMBOND(func_type)) {
	  for (j=0; j<2; j++) {
	    nr1=interaction->iatoms[i+1+j];
	    nr2=interaction->iatoms[i+2-j];
	    if ((*top->atoms.atomname[nr1][0] == 'H') && 
		((*top->atoms.atomname[nr2][0] == 'O') ||
		 (*top->atoms.atomname[nr2][0] == 'N')) &&
		in_list(nr1,isize,index) && in_list(nr2,isize,index))
	      add_dh(ddd,nr2,nr1);
	  }
	}
      }
    }
    for(func_type=0; func_type < F_NRE; func_type++) {
      interaction=&top->idef.il[func_type];
      for(i=0; i < interaction->nr; 
	  i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1) {
	/* next function */
	assert(func_type == top->idef.functype[interaction->iatoms[i]]);
	
	if ( interaction_function[func_type].flags & IF_DUMMY ) {
	  nr1=interaction->iatoms[i+1];
	  if ( *top->atoms.atomname[nr1][0]  == 'H') {
	    nr2=nr1-1;
	    stop=FALSE;
	    while (!stop && ( *top->atoms.atomname[nr2][0] == 'H'))
	      if (nr2)
		nr2--;
	      else
		stop=TRUE;
	    if ( !stop && ( ( *top->atoms.atomname[nr2][0] == 'O') ||
			    ( *top->atoms.atomname[nr2][0] == 'N') ) &&
		 in_list(nr1,isize,index) && in_list(nr2,isize,index) )
	      add_dh(ddd,nr2,nr1);
	  }
	}
      }
    }
  }
  snew(ddd->dptr,top->atoms.nr);
  for(i=0; (i<top->atoms.nr); i++)
    ddd->dptr[i] = NOTSET;
  for(i=0; (i<ddd->nrd); i++)
    ddd->dptr[ddd->d[i]] = i;
  
  return ddd;    
}

static t_gridcell ***init_grid(bool bBox,rvec box[],real rcut,ivec ngrid)
{
  t_gridcell ***grid;
  int i,y,z;
  
  if (bBox)
    for(i=0; i<DIM; i++)
      ngrid[i]=(box[i][i]/(1.2*rcut));
  
  if ( !bBox || (ngrid[XX]<3) || (ngrid[YY]<3) || (ngrid[ZZ]<3) )
    for(i=0; i<DIM; i++)
      ngrid[i]=1;
  printf("\nWill do grid-seach on %dx%dx%d grid, rcut=%g\n",
	 ngrid[XX],ngrid[YY],ngrid[ZZ],rcut);
  snew(grid,ngrid[ZZ]);
  for (z=0; z<ngrid[ZZ]; z++) {
    snew((grid)[z],ngrid[YY]);
    for (y=0; y<ngrid[YY]; y++)
      snew((grid)[z][y],ngrid[XX]);
  }
  return grid;
}

static void build_grid(t_hbdata *hb,rvec x[], rvec xshell,
		       bool bBox, matrix box, rvec hbox,
		       real rcut, real rshell,
		       ivec ngrid, t_gridcell ***grid)
{
  int     i,m,gr,xi,yi,zi,nr;
  atom_id *ad;
  ivec    grididx;
  rvec    invdelta,dshell;
  t_ncell *newgrid;
  bool    bDoRshell,bInShell,bAcc;
  real    rshell2=0;
  int     ggg;
  
  bDoRshell = (rshell > 0);
  rshell2   = sqr(rshell);
  bInShell  = TRUE;
  
  for(m=0; m<DIM; m++) {
    hbox[m]=box[m][m]*0.5;
    if (bBox) {
      invdelta[m]=ngrid[m]/box[m][m];
      if (1/invdelta[m] < rcut)
	fatal_error(0,"box shrank too much to keep using this grid\n");
    } else
      invdelta[m]=0;
  }
  grididx[XX]=0;
  grididx[YY]=0;
  grididx[ZZ]=0;
  for(gr=0; (gr<grNR); gr++) {
    /* resetting atom counts */
    for (zi=0; zi<ngrid[ZZ]; zi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (xi=0; xi<ngrid[XX]; xi++) {
	  grid[zi][yi][xi].d[gr].nr=0;
	  grid[zi][yi][xi].a[gr].nr=0;
	}
    /*if (debug)
      fprintf(debug,"Putting %d %s atoms into grid\n",nr[gr],grpnames[gr]);*/
    /* put atoms in grid cells */
    for(bAcc=FALSE; (bAcc<=TRUE); bAcc++) {
      if (bAcc) {
	nr = hb->a[gr]->nra;
	ad = hb->a[gr]->acc;
      }
      else {
	nr = hb->d[gr]->nrd;
	ad = hb->d[gr]->d;
      }
      for(i=0; (i<nr); i++) {
	/* check if we are inside the shell */
	/* if bDoRshell=FALSE then bInShell=TRUE always */
	if ( bDoRshell ) {
	  bInShell=TRUE;
	  rvec_sub(x[ad[i]],xshell,dshell);
	  if (bBox) 
	    for(m=DIM-1; m>=0 && bInShell; m--) {
	      if ( dshell[m] < -hbox[m] )
		rvec_inc(dshell,box[m]);
	      else if ( dshell[m] >= hbox[m] ) 
		dshell[m] -= 2*hbox[m];
	      /* if we're outside the cube, we're outside the sphere also! */
	      if ( (dshell[m]>rshell) || (-dshell[m]>rshell) )
		bInShell=FALSE;
	    }
	  /* if we're inside the cube, check if we're inside the sphere */
	  if (bInShell)
	    bInShell = norm2(dshell) < rshell2;
	}
	if ( bInShell ) {
	  if (bBox) 
	    for(m=DIM-1; m>=0; m--) {
	      /* put atom in the box */
	      while( x[ad[i]][m] < 0 ) 
		rvec_inc(x[ad[i]],box[m]);
	      while( x[ad[i]][m] >= box[m][m] ) 
		rvec_dec(x[ad[i]],box[m]);
	      /* determine grid index of atom */
	      ggg        = x[ad[i]][m]*invdelta[m];
	      grididx[m] = (ggg+ngrid[m]) % ngrid[m];
	    }
	  /* add atom to grid cell */
	  if (bAcc)
	    newgrid=&(grid[grididx[ZZ]][grididx[YY]][grididx[XX]].a[gr]);
	  else
	    newgrid=&(grid[grididx[ZZ]][grididx[YY]][grididx[XX]].d[gr]);
	  if (newgrid->nr >= newgrid->maxnr) {
	    newgrid->maxnr+=10;
	    srenew(newgrid->atoms, newgrid->maxnr);
	  }
	  newgrid->atoms[newgrid->nr]=ad[i];
	  newgrid->nr++;
	}
      }
    }
  } 
}

static void count_da_grid(ivec ngrid, t_gridcell ***grid, t_icell danr)
{
  int gr,xi,yi,zi;
  
  for(gr=0; gr<grNR; gr++) {
    danr[gr]=0;
    for (zi=0; zi<ngrid[ZZ]; zi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (xi=0; xi<ngrid[XX]; xi++) {
	  danr[gr] += grid[zi][yi][xi].d[gr].nr;
	}
  }
}

/* The grid loop.
 * Without a box, the grid is 1x1x1, so all loops are 1 long.
 * With a rectangular box (bTric==FALSE) all loops are 3 long.
 * With a triclinic box all loops are 3 long, except when a cell is
 * located next to one of the box edges which is not parallel to the
 * x/y-plane, in that case all cells in a line or layer are searched.
 * This could be implemented slightly more efficient, but the code
 * would get much more complicated.
 */
#define B(n,x,bTric,bEdge) ((n==1) ? x : bTric&&(bEdge) ? 0   : (x-1))
#define E(n,x,bTric,bEdge) ((n==1) ? x : bTric&&(bEdge) ? n-1 : (x+1))
#define GRIDMOD(j,n) (j+n)%(n)
#define LOOPGRIDINNER(x,y,z,xx,yy,zz,xo,yo,zo,n,bTric)\
     for(zz=B(n[ZZ],zo,bTric,FALSE); zz<=E(n[ZZ],zo,bTric,FALSE); zz++) {\
       z=GRIDMOD(zz,n[ZZ]);\
       for(yy=B(n[YY],yo,bTric,z==0||z==n[ZZ]-1); \
	  yy<=E(n[YY],yo,bTric,z==0||z==n[ZZ]-1); yy++) {\
	 y=GRIDMOD(yy,n[YY]);\
	 for(xx=B(n[XX],xo,bTric,y==0||y==n[YY]-1||z==0||z==n[ZZ]-1); \
	     xx<=E(n[XX],xo,bTric,y==0||y==n[YY]-1||z==0||z==n[ZZ]-1); xx++) {\
	   x=GRIDMOD(xx,n[XX]);
#define ENDLOOPGRIDINNER\
	 }\
       }\
     }\

static void dump_grid(FILE *fp, ivec ngrid, t_gridcell ***grid)
{
  int gr,x,y,z,sum[grNR];
  
  fprintf(fp,"grid %dx%dx%d\n",ngrid[XX],ngrid[YY],ngrid[ZZ]);
  for (gr=0; gr<grNR; gr++) {
    sum[gr]=0;
    fprintf(fp,"GROUP %d (%s)\n",gr,grpnames[gr]);
    for (z=0; z<ngrid[ZZ]; z+=2) {
      fprintf(fp,"Z=%d,%d\n",z,z+1);
      for (y=0; y<ngrid[YY]; y++) {
	for (x=0; x<ngrid[XX]; x++) {
	  fprintf(fp,"%3d",grid[x][y][z].d[gr].nr);
	  sum[gr]+=grid[z][y][x].d[gr].nr;
	  fprintf(fp,"%3d",grid[x][y][z].a[gr].nr);
	  sum[gr]+=grid[z][y][x].a[gr].nr;
	  
	}
	fprintf(fp," | ");
	if ( (z+1) < ngrid[ZZ] )
	  for (x=0; x<ngrid[XX]; x++) {
	    fprintf(fp,"%3d",grid[z+1][y][x].d[gr].nr);
	    sum[gr]+=grid[z+1][y][x].d[gr].nr;
	    fprintf(fp,"%3d",grid[z+1][y][x].a[gr].nr);
	    sum[gr]+=grid[z+1][y][x].a[gr].nr;
	  }
	fprintf(fp,"\n");
      }
    }
  }
  fprintf(fp,"TOTALS:");
  for (gr=0; gr<grNR; gr++)
    fprintf(fp," %d=%d",gr,sum[gr]);
  fprintf(fp,"\n");
}

/* New GMX record! 5 * in a row. Congratulations! 
 * Sorry, only four left.
 */
static void free_grid(ivec ngrid, t_gridcell ****grid)
{
  int y,z;
  t_gridcell ***g = *grid;
  
  for (z=0; z<ngrid[ZZ]; z++) {
    for (y=0; y<ngrid[YY]; y++) {
      sfree(g[z][y]);
    }
    sfree(g[z]);
  }
  sfree(g);
  g=NULL;
}

static void pbc_correct(rvec dx,matrix box,rvec hbox)
{
  int m;
  
  for(m=DIM-1; m>=0; m--) {
    if ( dx[m] < -hbox[m] )
      rvec_inc(dx,box[m]);
    else if ( dx[m] >= hbox[m] )
      rvec_dec(dx,box[m]);
  }
}

static bool low_is_hbond(int d,int a,int h,
			 rvec r_da,real rcut2, real ccut, 
			 rvec x[], bool bBox, matrix box,rvec hbox,
			 real *d_ha, real *ang,real d2)
{
  rvec r_dh,r_ha;
  real ca;
  
  if (d2 == 0) {
    rvec_sub(x[h],x[a],r_ha);
    if (bBox) 
      pbc_correct(r_ha,box,hbox);    
    d2 = iprod(r_ha,r_ha);
  }
  
  if ( d2 <= rcut2 ) {
    rvec_sub(x[d],x[h],r_dh);
    if (bBox)
      pbc_correct(r_dh,box,hbox);
    
    ca = cos_angle(r_dh,r_da);
    /* if angle is smaller, cos is larger */
    if (ca >= ccut) {
      *d_ha = sqrt(d2);
      *ang = acos(ca);
      return TRUE;
    }
  }
  return FALSE;
}

static bool is_hbond(t_hbdata *hb,int grpd,int grpa,int d,int a,
		     real rcut, real ccut, 
		     rvec x[], bool bBox, matrix box,rvec hbox,
		     real *d_ha, real *ang,bool bDA,int *hhh,
		     bool bContact)
{
  int  h,hh,id,ja;
  rvec r_da;
  real rc2,d2;
  
  if (d == a)
    return FALSE;

  id = donor_index(hb->d[grpd],d);
  ja = acceptor_index(hb->a[grpa],a);
  
  rvec_sub(x[d],x[a],r_da);
  if (bBox) 
    pbc_correct(r_da,box,hbox);    
  d2 = iprod(r_da,r_da);
  rc2 = rcut*rcut;
  
  if (d2 < rc2) {
    if (bContact) {
      *hhh = -1;
      return TRUE;
    }
    if (!bDA)
      d2 = 0;
    
    for(h=0; (h<hb->d[grpd]->nhydro[id]); h++) {
      hh = hb->d[grpd]->hydro[id][h];
      if (low_is_hbond(d,a,hh,r_da,rc2,ccut,x,bBox,box,hbox,d_ha,ang,d2)) {
	*hhh = hh;
	return TRUE;
      }
    }
  }
  return FALSE;
}

static void merge_hb(t_hbdata *hb,bool bTwo)
{
  int  i,inew,j,ii,jj,m,id,ia,ia0,id0,grp,ogrp,itest;
  int  *gNrD,*gNrA;
  
  inew = hb->nrhb;
    
  /* Check whether donors are also acceptors */
  fprintf(stderr,"Merging hbonds with Acceptor and Donor swapped\n");
  /* In the two arrays below we store the group number+1 (i.e. either
   * gr0 or gr1). This is necessary in case we have two groups.
   */
  snew(gNrD,hb->na);
  snew(gNrA,hb->nd);
  ia0 = id0 = 0;
  for(grp=gr0; (grp <= bTwo?gr1:gr0); grp++) {
    for(i=0; (i<hb->d[grp]->nrd); i++) {
      id = hb->d[grp]->d[i];
      ia = hb->a[grp]->aptr[id];
      gNrA[ia0++] = (ia != NOTSET) ? 1+grp : 0;
    }
    for(i=0; (i<hb->a[grp]->nra); i++) {
      ia = hb->a[grp]->acc[i];
      id = hb->d[grp]->dptr[ia];
      gNrD[id0++] = (id != NOTSET) ? 1+grp : 0;
    }
  }
  /* Consistency check */
  if ((ia0 != hb->nd) || (id0 != hb->na))
    fatal_error(0,"Unexepected inconsistency: ia0=%d, na=%d, id0=%d, nd=%d (%s,%d)",ia0,hb->na,id0,hb->nd,__FILE__,__LINE__);
    
  itest = 0;
  for(grp=gr0; (grp <= bTwo?gr1:gr0); grp++) {
    ogrp = bTwo ? OGRP : grp;
    for(i=0; (i<hb->nd); i++) {
      if (gNrA[i] == (1+grp)) {
	id = hb->d[grp]->d[i];
	ii = hb->a[grp]->aptr[id];
	for(j=0; (j<hb->na); j++) {
	  if (gNrD[j] == (1+ogrp)) {
	    ia = hb->a[ogrp]->acc[j];
	    jj = hb->d[ogrp]->dptr[ia];
	    if (id != ia) {
	      if (hb->hbmap[i][j].exist[0] && hb->hbmap[jj][ii].exist[0]) {
		for(m=0; (m<hb->max_frames/hb->wordlen); m++)
		  hb->hbmap[i][j].exist[0][m] |= hb->hbmap[jj][ii].exist[0][m];
		sfree(hb->hbmap[jj][ii].exist[0]);
		hb->hbmap[jj][ii].exist[0] = NULL;
		inew--;
	      }
	      itest++;
	    }
	  }
	}
      }
    }
  }
  sfree(gNrD);
  sfree(gNrA);
  fprintf(stderr,"Reduced number of hbonds from %d to %d\n",hb->nrhb,inew);
  fprintf(stderr,"tested %d pairs\n",itest);
  hb->nrhb = inew;
}

static void do_hbac(char *fn,t_hbdata *hb,real aver_nhb,bool bDump,bool bMerge)
{
  FILE *fp;
  int  i,j,k,m,ihb;
  bool bNorm=FALSE;
  double nhb = 0;
  real *rhbex;
  real *ct,tail,tail2,dtail,ct0;
  const real tol = 1e-3;
  int   nframes = hb->nframes;
  unsigned int **exist;
  int   nh,nhbonds,nhydro;
  
  /* build hbexist matrix in reals for autocorr */
  /* Allocate memory for computing ACF (rhbex) and aggregating the ACF (ct) */
  snew(rhbex,nframes);
  snew(ct,nframes);
  snew(exist,hb->maxhydro);
  
  /* Loop over hbonds */
  if (bDump) {
    fp = ffopen("debug-ac.xvg","w");
    for(j=0; (j<nframes); j++) {
      fprintf(fp,"%10g",hb->time[j]);
      for(i=0; (i<hb->nd); i++) {
	for(k=0; (k<hb->na); k++) {
	  if (bMerge) {
	    if (hb->hbmap[i][k].exist[0]) 
	      ihb = is_hb(hb->hbmap[i][k].exist[0],j);
	    else
	      ihb = 0;
	  }
	  else {
	    ihb = 0;
	    for(m=0; (m<hb->maxhydro) && !ihb ; m++)
	      ihb = ihb || ((hb->hbmap[i][k].exist[m]) && 
			    is_hb(hb->hbmap[i][k].exist[m],j));
	  }
	  fprintf(fp,"  %3d",ihb);
	}
	fprintf(fp,"\n");
      }
      ffclose(fp);
    }
  }
  /* Total number of hbonds analyzed here */
  nhbonds = 0;
  for(i=0; (i<hb->nd); i++) {
    for(k=0; (k<hb->na); k++) {
      if (bMerge) {
	if (hb->hbmap[i][k].exist[0]) {
	  exist[0] = hb->hbmap[i][k].exist[0];
	  nhydro = 1;
	}
	else
	  nhydro = 0;
      }
      else {
	nhydro = 0;
	for(m=0; (m<hb->maxhydro); m++)
	  if (hb->hbmap[i][k].exist[m]) {
	    exist[nhydro++] = hb->hbmap[i][k].exist[m];
	  }
      }
      for(nh=0; (nh<nhydro); nh++) {
	if ((((nhbonds+1) % 10) == 0) || (nhbonds+1 == hb->nrhb))
	  fprintf(stderr,"\rACF %d/%d",nhbonds+1,hb->nrhb);
	nhbonds++;
	for(j=0; (j<nframes); j++) {
	  ihb      = is_hb(exist[nh],j);
	  rhbex[j] = ihb-aver_nhb;
	  nhb     += ihb; 
	}
    
	/* The autocorrelation function is normalized after summation only */
	low_do_autocorr(NULL,NULL,
			nframes,1,-1,&rhbex,hb->time[1]-hb->time[0],eacNormal,1,
			FALSE,TRUE,bNorm,FALSE,0,-1,0,1);
	for(j=0; (j<nframes/2); j++) {
	  ct[j] += rhbex[j];
	}
      }
    }
  }
  fprintf(stderr,"\n");
  sfree(rhbex);
  
  /* Normalize */
  ct0 = ct[0];
  for(j=0; (j<nframes/2); j++)
    ct[j] /= ct0;
  
  /* Determine tail value for statistics */
  tail  = 0;
  tail2 = 0;
  for(j=nframes/4; (j<nframes/2); j++) {
    tail += ct[j];
    tail2 += ct[j]*ct[j];
  }
  tail  /= (nframes/2 - nframes/4);
  tail2 /= (nframes/2 - nframes/4);
  dtail  = sqrt(tail2-tail*tail);
  
  /* Check whether the ACF is long enough */
  if (dtail > tol) {
    printf("\nWARNING: Correlation function is probably not long enough\n"
	   "because the standard deviation in the tail of C(t) > %g\n"
	   "Total number of hbonds found: %g\n"
	   "Tail value (average C(t) over second half of acf): %g +/- %g\n",
	   tol,nhb,tail,dtail);
  }
  fp = xvgropen(fn, "Hydrogen Bond Autocorrelation","Time (ps)","C(t)");
  for(j=0; (j<nframes/2); j++)
    fprintf(fp,"%10g  %10g  %10g\n",
	    hb->time[j]-hb->time[0],(ct[j]-tail)/(1-tail),ct[j]);
  
  fclose(fp);
  do_view(fn,NULL);
  sfree(ct);
}

static void init_hbframe(t_hbdata *hb,int nframes,real t)
{
  int i,j,m;
  
  hb->time[nframes] = t;
  hb->nhb[nframes]  = 0;
  for (i=0; (i<max_hx); i++)
    hb->nhx[nframes][i]=0;
  if (hb->bHBmap)
    for (i=0; (i<hb->nd); i++)
      for (j=0; (j<hb->na); j++)
	for (m=0; (m<hb->maxhydro); m++)
	  if (hb->hbmap[i][j].exist[m])
	    set_hb(hb->hbmap[i][j].exist[m],nframes,HB_NO);
}


int gmx_hbond(int argc,char *argv[])
{
  static char *desc[] = {
    "g_hbond computes and analyzes hydrogen bonds. Hydrogen bonds are",
    "determined based on cutoffs for the angle Donor - Hydrogen - Acceptor",
    "(zero is extended) and the distance Hydrogen - Acceptor.",
    "OH and NH groups are regarded as donors, O is an acceptor always,",
    "N is an acceptor by default, but this can be switched using",
    "[TT]-nitacc[tt]. Dummy hydrogen atoms are assumed to be connected",
    "to the first preceding non-hydrogen atom.[PAR]",
    
    "You need to specify two groups for analysis, which must be either",
    "identical or non-overlapping. All hydrogen bonds between the two",
    "groups are analyzed.[PAR]",
    
    "If you set -shell, you will be asked for an additional index group",
    "which should contain exactly one atom. In this case, only hydrogen",
    "bonds between atoms within the shell distance from the one atom are",
    "considered.[PAR]"
    
    "It is also possible to analyse specific hydrogen bonds with",
    "[TT]-sel[tt]. This index file must contain a group of atom triplets",
    "Donor Hydrogen Acceptor, in the following way:[PAR]",
    
    "[TT]",
    "[ selected ][BR]",
    "     20    21    24[BR]",
    "     25    26    29[BR]",
    "      1     3     6[BR]",
    "[tt][BR]",
    "Note that the triplets need not be on separate lines.",
    "Each atom triplet specifies a hydrogen bond to be analyzed,",
    "note also that no check is made for the types of atoms.[PAR]",
    
    "[TT]-ins[tt] turns on computing solvent insertion into hydrogen bonds.",
    "In this case an additional group must be selected, specifying the",
    "solvent molecules.[PAR]",
    
    "[BB]Output:[bb][BR]",
    "[TT]-num[tt]:  number of hydrogen bonds as a function of time.[BR]",
    "[TT]-ac[tt]:   average over all autocorrelations of the existence",
    "functions (either 0 or 1) of all hydrogen bonds.[BR]",
    "[TT]-dist[tt]: distance distribution of all hydrogen bonds.[BR]",
    "[TT]-ang[tt]:  angle distribution of all hydrogen bonds.[BR]",
    "[TT]-hx[tt]:   the number of n-n+i hydrogen bonds as a function of time",
    "where n and n+i stand for residue numbers and i ranges from 0 to 6.",
    "This includes the n-n+3, n-n+4 and n-n+5 hydrogen bonds associated",
    "with helices in proteins.[BR]",
    "[TT]-hbn[tt]:  all selected groups, donors, hydrogens and acceptors",
    "for selected groups, all hydrogen bonded atoms from all groups and",
    "all solvent atoms involved in insertion.[BR]",
    "[TT]-hbm[tt]:  existence matrix for all hydrogen bonds over all",
    "frames, this also contains information on solvent insertion",
    "into hydrogen bonds. Ordering is identical to that in [TT]-hbn[tt]",
    "index file.[BR]",
    "[TT]-dan[tt]: write out the number of donors and acceptors analyzed for",
    "each timeframe. This is especially usefull when using [TT]-shell[tt]."
  };
  
  static real acut=30, abin=1, rcut=0.35, rbin=0.005, rshell=-1,maxnhb=0;
  static bool bNitAcc=TRUE,bInsert=FALSE,bDA=TRUE,bDump=FALSE,bMerge=TRUE;
  static bool bContact=FALSE;
  /* options */
  t_pargs pa [] = {
    { "-ins",  FALSE,  etBOOL, {&bInsert},
      "Analyze solvent insertion" },
    { "-a",    FALSE,  etREAL, {&acut},
      "Cutoff angle (degrees, Donor - Hydrogen - Acceptor)" },
    { "-r",    FALSE,  etREAL, {&rcut},
      "Cutoff radius (nm, X - Acceptor, see next option)" },
    { "-da",   FALSE,  etBOOL, {&bDA},
      "Use distance Donor-Acceptor (if TRUE) or Hydrogen-Acceptor (FALSE)" },
    { "-abin", FALSE,  etREAL, {&abin},
      "Binwidth angle distribution (degrees)" },
    { "-rbin", FALSE,  etREAL, {&rbin},
      "Binwidth distance distribution (nm)" },
    { "-nitacc",FALSE, etBOOL, {&bNitAcc},
      "Regard nitrogen atoms as acceptors" },
    { "-contact",FALSE,etBOOL, {&bContact},
      "Do not look for hydrogen bonds, but merely for contacts within the cut-off distance" },
    { "-shell", FALSE, etREAL, {&rshell},
      "when > 0, only calculate hydrogen bonds within # nm shell around "
      "one particle" },
    { "-dump",  FALSE, etBOOL, {&bDump},
      "Dump all hydrogen bond ACFs (maximum 1000) in a single xvg file for debugging" },
    { "-max_hb",FALSE, etREAL, {&maxnhb},
      "Theoretical maximum number of hydrogen bonds used for normalizing HB autocorrelation function. Can be useful in case the program estimates it wrongly" },
    { "-merge", FALSE, etBOOL, {&bMerge},
      "H-bonds between the same donor and accepter, but with different hydrogen are treated as a single H-bond. Mainly important for the ACF." }
  };

  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD  },
    { efTPX, NULL,   NULL,     ffREAD  },
    { efNDX, NULL,   NULL,     ffOPTRD },
    { efLOG, "-g",   "hbond",  ffOPTWR },
    { efNDX, "-sel", "select", ffOPTRD },
    { efXVG, "-num", "hbnum",  ffWRITE },
    { efXVG, "-ac",  "hbac",   ffOPTWR },
    { efXVG, "-dist","hbdist", ffOPTWR },
    { efXVG, "-ang", "hbang",  ffOPTWR },
    { efXVG, "-hx",  "hbhelix",ffOPTWR },
    { efNDX, "-hbn", "hbond",  ffOPTWR },
    { efXPM, "-hbm", "hbmap",  ffOPTWR },
    { efXVG, "-dan",  "danum",  ffOPTWR }
  };
#define NFILE asize(fnm)
  
  char  hbmap [HB_NR]={ ' ',    'o',      '-',       '*' };
  char *hbdesc[HB_NR]={ "None", "Present", "Inserted", "Present & Inserted" };
  t_rgb hbrgb [HB_NR]={ {1,1,1},{1,0,0},   {0,0,1},    {1,0,1} };
  
  int     status;
  t_topology top;
  t_inputrec ir;
  int     natoms,nframes=0,shatom;
  int     *isize;
  char    **grpnames;
  atom_id **index;
  rvec    *x,hbox;
  matrix  box;
  real    t,ccut,dist,ang;
  double  max_nhb,aver_nhb;
  int     h,i,j,k,l,start,end,id,ja;
  int     xi,yi,zi,ai;
  int     xj,yj,zj,aj,xjj,yjj,zjj;
  int     xk,yk,zk,ak,xkk,ykk,zkk;
  bool    bSelected,bStop,bTwo,new,was,bBox,bTric;
  int     *adist,*rdist;
  int        grp,nabin,nrbin,bin,resdist;
  t_hbdata   *hb;
  FILE       *fp,*fpins=NULL,*fplog;
  t_gridcell ***grid;
  t_ncell    *icell,*jcell,*kcell;
  ivec       ngrid;
    
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,NFILE,fnm,asize(pa),
		    pa,asize(desc),desc,0,NULL);

  /* process input */
  bSelected = opt2bSet("-sel",NFILE,fnm);
  ccut = cos(acut*DEG2RAD);
  
  if (bContact) {
    if (bSelected)
      fatal_error(0,"Can not analyze selected contacts: turn off -sel");
    if (bInsert)
      fatal_error(0,"Can not analyze inserted contacts: turn off -ins");
  }
  
  /* Initiate main data structure! */
  hb = mk_hbdata(opt2bSet("-hbm",NFILE,fnm) || opt2bSet("-ac",NFILE,fnm),
		 opt2bSet("-da",NFILE,fnm),bMerge);
  
  /* get topology */
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&i,&t,&t,
	   &ir,box,&natoms,NULL,NULL,NULL,&top);
  
  snew(grpnames,grNR);
  snew(index,grNR);
  snew(isize,grNR);
  if (bSelected) {
    /* analyze selected hydrogen bonds */
    fprintf(stderr,"Select group with selected atoms:\n");
    get_index(&(top.atoms),opt2fn("-sel",NFILE,fnm),
	      1,isize,index,grpnames);
    if (isize[0] % 3)
      fatal_error(0,"Number of atoms in group '%s' not a multiple of 3\n"
		  "and therefore cannot contain triplets of "
		  "Donor-Hydrogen-Acceptor",grpnames[0]);
    bTwo=FALSE;
    
    snew(hb->d[gr0],1);
    snew(hb->a[gr0],1);
    for(i=0; (i<isize[0]); i+=3) {
      int dd = index[0][i];
      int hh = index[0][i+1];
      int aa = index[0][i+2];
      add_dh (hb->d[gr0],dd,hh);
      add_acc(hb->a[gr0],aa);
      /* Should this be here ? */
      add_hbond(hb,dd,hh,aa,gr0,gr0,-1,FALSE,bMerge);
    }
    fprintf(stderr,"Analyzing %d selected hydrogen bonds from '%s'\n",
	    hb->nrhb,grpnames[0]);
  } 
  else {
    /* analyze all hydrogen bonds: get group(s) */
    fprintf(stderr,"Specify 2 groups to analyze:\n");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      2,isize,index,grpnames);
    
    /* check if we have two identical or two non-overlapping groups */
    bTwo = isize[0] != isize[1];
    for (i=0; (i<isize[0]) && !bTwo; i++)
      bTwo = index[0][i] != index[1][i];
    if (bTwo) {
      fprintf(stderr,"Checking for overlap...\n");
      for (i=0; i<isize[0]; i++)
	for (j=0; j<isize[1]; j++)
	  if (index[0][i] == index[1][j]) 
	    fatal_error(0,"Partial overlap between groups '%s' and '%s'",
			grpnames[0],grpnames[1]);
    }
    if (bContact && !bTwo)
      fatal_error(0,"I need two different groups for a contact analysis");
    if (bTwo)
      fprintf(stderr,"Calculating %s "
	      "between two groups of %d and %d atoms\n",
	      bContact ? "contacts" : "hydrogen bonds",isize[0],isize[1]);
    else
      fprintf(stderr,"Calculating hydrogen bonds in one group of %d atoms\n",
	      isize[0]);
  }
  if (bInsert) {
    fprintf(stderr,"Specify group for insertion analysis:\n");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      1,&(isize[grI]),&(index[grI]),&(grpnames[grI]));
    fprintf(stderr,"Checking for overlap...\n");
    for (i=0; i<isize[grI]; i++)
      for (grp=0; grp<(bTwo?2:1); grp++)
	for (j=0; j<isize[grp]; j++)
	  if (index[grI][i] == index[grp][j]) 
	    fatal_error(0,"Partial overlap between groups '%s' and '%s'",
			grpnames[grp],grpnames[grI]);
    fpins=ffopen("insert.dat","w");
    fprintf(fpins,"%4s: %15s -> %15s (%7s) - %15s (%7s)\n",
	    "time","insert","donor","distang","acceptor","distang");
  }
  
  /* search donors and acceptors in groups */
  for (i=0; (i<grNR); i++)
    if ( ((i==gr0) && !bSelected ) ||
	 ((i==gr1) && bTwo ) ||
	 ((i==grI) && bInsert ) ) {
      if (bContact) {
	hb->a[i] = search_acceptors(&top,isize[i],index[i],bNitAcc,
				    TRUE,(i==gr0));
	hb->d[i] = search_donors   (&top,isize[i],index[i],
				    TRUE,(i==gr1));
      }
      else {
	hb->a[i] = search_acceptors(&top,isize[i],index[i],bNitAcc,FALSE,TRUE);
	hb->d[i] = search_donors   (&top,isize[i],index[i],FALSE,TRUE);
      }
      fprintf(stderr,"Found %d donors and %d acceptors in group '%s'\n",
	      hb->d[i]->nrd,hb->a[i]->nra,grpnames[i]);
    }
  /*if (bSelected)
    snew(donors[gr0D], dons[gr0D].nrd);*/
  if (!bTwo) {
    hb->d[gr1]   = hb->d[gr0];
    hb->a[gr1]   = hb->a[gr0];
    /*donors[gr1D] = donors[gr0D];*/
  }
  if (!bInsert) {
    /* Initiate structures with zero entries */
    snew(hb->d[grI],1);
    snew(hb->a[grI],1);
  }
  /* Generate hbond data structure */
  mk_hbmap(hb,bTwo,bInsert);
  
  /* check input */
  bStop=FALSE;
  for (grp=gr0; (grp <= (bTwo ? gr1 : gr0)); grp++)
    if (hb->d[grp]->nrd + hb->a[grp]->nra == 0) {
      fprintf(stderr,"No Donors or Acceptors found in group '%s'\n",
	      grpnames[grp]);
      bStop=TRUE;
    }
  if (!bStop) {
    if (hb->d[gr0]->nrd + hb->d[gr1]->nrd == 0) {
      fprintf(stderr,"No Donors found\n");
      bStop=TRUE;
    }
    if (hb->a[gr0]->nra + hb->a[gr1]->nra == 0) {
      fprintf(stderr,"No Acceptors found\n");
      bStop=TRUE;
    }
  }
  if (bStop)
    fatal_error(0,"Nothing to be done");
  if ( bInsert && ((hb->d[grI]->nrd == 0) || (hb->a[grI]->nra == 0)) )
    fatal_error(0,"No %s%s%s found in insertion group '%s'\n",
		(hb->d[grI]->nrd == 0) ? "donors" : "",
		((hb->d[grI]->nrd == 0) && 
		 (hb->a[grI]->nra == 0)) ? " and " : "",
		(hb->a[grI]->nra == 0) ? "acceptors" : "",
		grpnames[grI]);

  shatom=0;
  if (rshell > 0) {
    int shisz;
    atom_id *shidx;
    char *shgrpnm;
    /* get index group with atom for shell */
    do {
      fprintf(stderr,"Select atom for shell (1 atom):\n");
      get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
		1,&shisz,&shidx,&shgrpnm);
      if (shisz != 1)
	fprintf(stderr,"group contains %d atoms, should be 1 (one)\n",shisz);
    } while(shisz != 1);
    shatom = shidx[0];
    fprintf(stderr,"Will calculate hydrogen bonds within a shell "
	    "of %g nm around atom %i\n",rshell,shatom+1);
  }

  /* Analyze trajectory */
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( natoms > top.atoms.nr )
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);
		
  bBox  = ir.ePBC!=epbcNONE;
  grid  = init_grid(bBox, box, rcut, ngrid);
  nabin = acut/abin;
  nrbin = rcut/rbin;
  snew(adist,nabin);
  snew(rdist,nrbin);
  
  do {
    bTric = bBox && TRICLINIC(box);
    build_grid(hb,x,x[shatom], bBox,box,hbox, rcut, rshell, ngrid,grid);
    if (debug)
      dump_grid(debug, ngrid, grid);
    
    add_frames(hb,nframes);
    init_hbframe(hb,nframes,t);

    if (hb->bDAnr)
      count_da_grid(ngrid, grid, hb->danr[nframes]);
    /* loop over all gridcells (xi,yi,zi)      */
    /* Removed confusing macro, DvdS 27/12/98  */
    for(xi=0; (xi<ngrid[XX]); xi++)
      for(yi=0; (yi<ngrid[YY]); yi++)
	for(zi=0; (zi<ngrid[ZZ]); zi++) {
	
	  /* loop over donor groups gr0 (always) and gr1 (if necessary) */
	  for (grp=gr0; (grp <= (bTwo?gr1:gr0)); grp++) {
	    icell=&(grid[zi][yi][xi].d[grp]);
	    
	    /* loop over all hydrogen atoms from group (grp) 
	     * in this gridcell (icell) 
	     */
	    for (ai=0; (ai<icell->nr); ai++) {
	      i  = icell->atoms[ai];
	      
	      /* loop over all adjacent gridcells (xj,yj,zj) */
	      /* This is a macro!!! */
	      LOOPGRIDINNER(xj,yj,zj,xjj,yjj,zjj,xi,yi,zi,ngrid,bTric) {
		jcell=&(grid[zj][yj][xj].a[OGRP]);
		/* loop over acceptor atoms from other group (OGRP) 
		 * in this adjacent gridcell (jcell) 
		 */
		for (aj=0; (aj<jcell->nr); aj++) {
		  j = jcell->atoms[aj];
		  
		  if ((bSelected && (j == i)) || (!bSelected)) {
		    /* check if this once was a h-bond */
		    if (is_hbond(hb,grp,OGRP,i,j,rcut,ccut,x,
				 bBox,box,hbox,&dist,&ang,bDA,&h,bContact)) {
		      /* add to index if not already there */
		      /* Add a hbond */
		      add_hbond(hb,i,j,h,grp,OGRP,nframes,FALSE,bMerge);

		      /* make angle and distance distributions */
		      ang*=RAD2DEG;
		      adist[(int)( ang/abin)]++;
		      rdist[(int)(dist/rbin)]++;
		      
		      if (!bTwo) {
			int id = donor_index(hb->d[grp],i);
			int ia = acceptor_index(hb->a[OGRP],j);
			resdist=abs(top.atoms.atom[id].resnr-
				    top.atoms.atom[ia].resnr);
			if (resdist >= max_hx) 
			  resdist = max_hx-1;
			hb->nhx[nframes][resdist]++;
		      }
		    }
		    if (bInsert && bSelected) {
		      /* this has been a h-bond, or we are analyzing 
			 selected bonds: check for inserted */
		      bool ins_d, ins_a;
		      real ins_d_dist, ins_d_ang, ins_a_dist, ins_a_ang;
		      int  ins_d_k=0,ins_a_k=0;
		      
		      ins_d=ins_a=FALSE;
		      ins_d_dist=ins_d_ang=ins_a_dist=ins_a_ang=1e6;
		      
		      /* loop over gridcells adjacent to i (xk,yk,zk) */
		      LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xi,yi,zi,ngrid,bTric){
			kcell=&(grid[zk][yk][xk].a[grI]);
			/* loop over acceptor atoms from ins group 
			   in this adjacent gridcell (kcell) */
			for (ak=0; (ak<kcell->nr); ak++) {
			  k=kcell->atoms[ak];
			  if (is_hbond(hb,grp,grI,i,k,rcut,ccut,x,
				       bBox,box,hbox,&dist,&ang,bDA,&h,
				       bContact)) {
			    if (dist < ins_d_dist) {
			      ins_d=TRUE;
			      ins_d_dist=dist;
			      ins_d_ang =ang ;
			      ins_d_k   =k   ;
			    }
			  }
			}
		      }
		      ENDLOOPGRIDINNER;
			/* loop over gridcells adjacent to j (xk,yk,zk) */
		      LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xj,yj,zj,ngrid,bTric){
			kcell=&grid[zk][yk][xk].d[grI];
			/* loop over hydrogen atoms from ins group 
			   in this adjacent gridcell (kcell) */
			for (ak=0; ak<kcell->nr; ak++) {
			  k=kcell->atoms[ak];
			  if (is_hbond(hb,grI,OGRP,k,j,rcut,ccut,x,
				       bBox,box,hbox,&dist,&ang,bDA,&h,
				       bContact) ) {
			    if (dist<ins_a_dist) {
			      ins_a=TRUE;
			      ins_a_dist=dist;
			      ins_a_ang =ang ;
			      ins_a_k   =k   ;
			    }
			  }
			}
		      }
		      ENDLOOPGRIDINNER;
			
		      {
			if (ins_d && ins_a && 
			    is_hbond(hb,grI,grI,ins_d_k,ins_a_k,rcut,ccut,x,
				     bBox,box,hbox,&dist,&ang,bDA,&h,
				     bContact)) {
			  /* add to hbond index if not already there */
			  add_hbond(hb,ins_d_k,ins_a_k,h,grI,OGRP,
				    nframes,TRUE,bMerge);
			  
			  /* print insertion info to file */
			  /*fprintf(fpins,
			    "%4g: %4u:%3.3s%4d%3.3s -> "
			    "%4u:%3.3s%4d%3.3s (%4.2f,%2.0f) - "
			    "%4u:%3.3s%4d%3.3s (%4.2f,%2.0f)\n",t,
			    a[grIA][ins_d_k]+1,
			    *top.atoms.resname[top.atoms.atom[a[grIA][ins_d_k]].resnr],
			    top.atoms.atom[a[grIA][ins_d_k]].resnr+1,
			    *top.atoms.atomname[a[grIA][ins_d_k]],
			    a[grp+grD][i]+1,
			    *top.atoms.resname[top.atoms.atom[a[grp+grD][i]].resnr],
			    top.atoms.atom[a[grp+grD][i]].resnr+1,
			    *top.atoms.atomname[a[grp+grD][i]],
			    ins_d_dist,ins_d_ang*RAD2DEG,
			    a[OGRP+grA][j]+1,
			    *top.atoms.resname[top.atoms.atom[a[OGRP+grA][j]].resnr],
			    top.atoms.atom[a[OGRP+grA][j]].resnr+1,
			    *top.atoms.atomname[a[OGRP+grA][j]],
			    ins_a_dist,ins_a_ang*RAD2DEG);*/
			}
		      }
		    }
		  }
		} /* for aj  */
	      }
	      ENDLOOPGRIDINNER;
	    } /* for ai  */
	  } /* for grp */
	} /* for xi,yi,zi */
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  
  free_grid(ngrid,&grid);
  
  close_trj(status);
  if (bInsert)
    fclose(fpins);
  
  if (hb->nrhb==0)
    fprintf(stderr,"No hydrogen bonds found!!\n");
  else {
    fprintf(stderr,"Found %d different hydrogen bonds in trajectory\n",
	    hb->nrhb);
    
    if (bMerge)
      merge_hb(hb,bTwo);

    aver_nhb = 0;    
    fp = xvgropen(opt2fn("-num",NFILE,fnm),bContact ? "Contacts" :
		  "Hydrogen Bonds","Time","Number");
    for(i=0; (i<nframes); i++) {
      fprintf(fp,"%10g %10d\n",hb->time[i],hb->nhb[i]);
      aver_nhb += hb->nhb[i];
    }
    fclose(fp);
    aver_nhb /= nframes;
    
    if (opt2bSet("-dist",NFILE,fnm)) {
      long sum;
      
      sum=0;
      for(i=0; i<nrbin; i++)
	sum+=rdist[i];
      
      fp = xvgropen(opt2fn("-dist",NFILE,fnm),
		    "Hydrogen Bond Distribution",
		    "Hydrogen - Acceptor Distance (nm)","");
      for(i=0; i<nrbin; i++)
	fprintf(fp,"%10g %10g\n",(i+0.5)*rbin,rdist[i]/(rbin*(real)sum));
      fclose(fp);
    }
  
    if (opt2bSet("-ang",NFILE,fnm)) {
      long sum;
      
      sum=0;
      for(i=0; i<nabin; i++)
	sum+=adist[i];
      
      fp = xvgropen(opt2fn("-ang",NFILE,fnm),
		    "Hydrogen Bond Distribution",
		    "Donor - Hydrogen - Acceptor Angle (\\SO\\N)","");
      for(i=0; i<nabin; i++)
	fprintf(fp,"%10g %10g\n",(i+0.5)*abin,adist[i]/(abin*(real)sum));
      fclose(fp);
    }
    
    if (opt2bSet("-hx",NFILE,fnm)) {
      fp = xvgropen(opt2fn("-hx",NFILE,fnm),
		    "Hydrogen Bonds","Time(ps)","Count");
      xvgr_legend(fp,NRHXTYPES,hxtypenames);
      for(i=0; i<nframes; i++) {
	fprintf(fp,"%10g",hb->time[i]);
	for (j=0; j<max_hx; j++)
	  fprintf(fp," %6d",hb->nhx[i][j]);
	fprintf(fp,"\n");
      }
      fclose(fp);
    }
    
    /* Compute maximum possible number of different hbonds */
    if (maxnhb > 0)
      max_nhb = maxnhb;
    else {
      if (bTwo) 
	max_nhb = 0.5*(hb->d[gr0]->nrd*hb->a[gr1]->nra +
		       hb->d[gr1]->nrd*hb->a[gr0]->nra);
      else 
	max_nhb = 0.5*(hb->d[gr0]->nrd*hb->a[gr0]->nra);
    }
    printf("Average number of hbonds per timeframe %.3f out of %g possible\n",
	   aver_nhb,max_nhb);
    if (opt2bSet("-ac",NFILE,fnm))
      do_hbac(opt2fn("-ac",NFILE,fnm),hb,aver_nhb/max_nhb,bDump,bMerge);
    
    if (opt2bSet("-hbm",NFILE,fnm)) {
      t_matrix mat;
      int id,ia,hh,x,y;
      
      mat.nx=nframes;
      mat.ny=hb->nrhb;
      snew(mat.matrix,mat.nx);
      for(x=0; (x<mat.nx); x++){
	snew(mat.matrix[x],mat.ny);
	y=0;
	for(id=0; (id<hb->nd); id++) 
	  for(ia=0; (ia<hb->na); ia++) 
	    for(hh=0; (hh<hb->maxhydro); hh++)
	      if (hb->hbmap[id][ia].exist[hh])
		for(j=0; (j<mat.ny); j++)
		  mat.matrix[x][y++]=is_hb(hb->hbmap[id][ia].exist[hh],x);
	
	mat.axis_x=hb->time;
	snew(mat.axis_y,mat.ny);
	for(j=0; j<mat.ny; j++)
	  mat.axis_y[j]=j;
	sprintf(mat.title,"Hydrogen Bond Existence Map");
	sprintf(mat.legend,"Hydrogen Bonds");
	sprintf(mat.label_x,"time (ps)");
	sprintf(mat.label_y,"Hydrogen Bond Index");
	mat.bDiscrete=TRUE;
	if (bInsert)
	  mat.nmap=HB_NR;
	else
	mat.nmap=2;
	snew(mat.map,mat.nmap);
	for(i=0; i<mat.nmap; i++) {
	  mat.map[i].code.c1=hbmap[i];
	  mat.map[i].desc=hbdesc[i];
	  mat.map[i].rgb=hbrgb[i];
	}
	fp = opt2FILE("-hbm",NFILE,fnm,"w");
	write_xpm_m(fp, mat);
	fclose(fp);
	for(x=0; x<mat.nx; x++)
	  sfree(mat.matrix[x]);
	sfree(mat.axis_y);
	sfree(mat.matrix);
	sfree(mat.map);
      }
    }
    
    if (opt2bSet("-hbn",NFILE,fnm)) {
      fp = opt2FILE("-hbn",NFILE,fnm,"w");
      if (opt2bSet("-g",NFILE,fnm)) {
	fplog = ffopen(opt2fn("-g",NFILE,fnm),"w");
	fprintf(fplog,"# %10s  %12s  %12s\n","Donor","Hydrogen","Acceptor");
      }
      else
	fplog = NULL;
      for (grp=gr0; grp<=(bTwo?gr1:gr0); grp++) {
	fprintf(fp,"[ %s ]",grpnames[grp]);
	for (i=0; i<isize[grp]; i++) {
	  fprintf(fp,(i%15)?" ":"\n");
	  fprintf(fp,"%4u",index[grp][i]+1);
	}
	fprintf(fp,"\n");
	fprintf(fp,"[ donors_hydrogens_%s ]",grpnames[grp]);
	for (i=0; (i<hb->d[grp]->nrd); i++) {
	  for(j=0; (j<hb->d[grp]->nhydro[i]); j++)
	    fprintf(fp,"%4u %4u",hb->d[grp]->d[i]+1,
		    hb->d[grp]->hydro[i][j]+1);
	  fprintf(fp,"\n");
	}
	fprintf(fp,"[ acceptors_%s ]",grpnames[grp]);
	for (i=0; (i<hb->a[grp]->nra); i++) {
	  fprintf(fp,(i%15)?" ":"\n");
	  fprintf(fp,"%4u",hb->a[grp]->acc[i]+1);
	}
	fprintf(fp,"\n");
      }
      if (bTwo)
	fprintf(fp,"[ hbonds_%s-%s ]\n",grpnames[0],grpnames[1]);
      else
	fprintf(fp,"[ hbonds_%s ]\n",grpnames[0]);
      for(i=0; (i<hb->nrhb); i++) {
	int  ddd,hhh,aaa;
	char ds[32],hs[32],as[32];
	
	/*	ddd = hb->hb[i].d;
	hhh = hb->hb[i].h;
	aaa = hb->hb[i].a;
	sprintf(ds,"%s",mkatomname(&top.atoms,ddd));
	    sprintf(hs,"%s",mkatomname(&top.atoms,hhh));
	    sprintf(as,"%s",mkatomname(&top.atoms,aaa));
	    fprintf(fp,"%6u %6u %6u\n",ddd+1,hhh+1,aaa+1);
	    if (fplog) 
	      fprintf(fplog,"%12s  %12s  %12s\n",ds,hs,as);
	*/
      }
      if (bInsert) {
	if (bTwo)
	  fprintf(fp,"[ insert_%s->%s-%s ]",
		  grpnames[2],grpnames[0],grpnames[1]);
	else
	  fprintf(fp,"[ insert_%s->%s ]",grpnames[2],grpnames[0]);
	j=0;
	/*	for(i=0; (i<hb->nrhb); i++) {
	  if (hb->hb[i].bInsert) {
	    fprintf(fp,(j%15)?" ":"\n");
	    fprintf(fp,"%4d",i+1);
	    j++;
	  }
	  }*/
	fprintf(fp,"\n");
      }
      fclose(fp);
      if (fplog)
	fclose(fplog);
    }
    
    if (hb->bDAnr) {
      int  i,j,nleg;
      char **legnames;
      char buf[STRLEN];
      
#define USE_THIS_GROUP(j) ( (j == gr0) || (bTwo && (j == gr1)) || (bInsert && (j == grI)) )
      
      fp = xvgropen(opt2fn("-da",NFILE,fnm),
		    "Donors and Acceptors","Time(ps)","Count");
      nleg = (bTwo?2:1)*2 + (bInsert?2:0);
      snew(legnames,nleg);
      i=0;
      for(j=0; j<grNR; j++)
	if ( USE_THIS_GROUP(j) ) {
	  sprintf(buf,"Donors %s",grpnames[j]);
	  legnames[i++] = strdup(buf);
	  sprintf(buf,"Acceptors %s",grpnames[j]);
	  legnames[i++] = strdup(buf);
	}
      assert(i==nleg);
      xvgr_legend(fp,nleg,legnames);
      for(i=0; i<nframes; i++) {
	fprintf(fp,"%10g",hb->time[i]);
	for (j=0; (j<grNR); j++)
	  if ( USE_THIS_GROUP(j) )
	    fprintf(fp," %6d",hb->danr[i][j]);
	fprintf(fp,"\n");
      }
      fclose(fp);
    }
  }
  
  thanx(stderr);
  
  return 0;
}
  
