/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gromacs Runs One Microsecond At Cannonball Speeds
 */

#include <math.h>
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
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

enum { gr0D, gr0H, gr0A, gr1D, gr1H, gr1A, grID, grIH, grIA, grNR };
#define grD gr0D
#define grH gr0H
#define grA gr0A
#define gr0 gr0D
#define gr1 gr1D
#define grI grID
#define grINC (gr1-gr0)

typedef struct {
  int nr;
  int maxnr;
  atom_id *atoms;
} t_ncell;

typedef t_ncell t_gridcell[grNR];
typedef int     t_icell[grNR];

typedef struct {
  int a;
  int nr;
} t_acceptor;

typedef struct {
  int nrhb;
  int maxnr;
  t_acceptor *hb;
} t_donor;

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

static void add_acc(int i, int *max_nr,int *nr,atom_id **a)
{
  (*nr)++;
  if ( *nr > *max_nr ) {
    (*max_nr)+=10;
    srenew(*a,*max_nr);
  }
  (*a)[*nr-1]=i;
}

static void search_acceptors(t_topology *top, int isize, atom_id *index,
		      int *nr_a, atom_id **a, bool bNitAcc)
{
  int i,max_nr_a;
  
  max_nr_a=*nr_a;
  for (i=0; i<top->atoms.nr; i++)
    if ( ( *top->atoms.atomname[i][0] == 'O' || 
	   ( bNitAcc && ( *top->atoms.atomname[i][0] == 'N' ) ) ) &&
	 in_list(i,isize,index) )
      add_acc(i,&max_nr_a,nr_a,a);
  srenew(*a,*nr_a);
}

static void add_dh(int id, int ih, int *max_nr,int *nr,atom_id **d,atom_id **h,
		   bool *bDonor)
{
  if (!bDonor[ih]) {
    (*nr)++;
    if ( *nr > *max_nr ) {
      (*max_nr)+=10;
      srenew(*d,*max_nr);
      srenew(*h,*max_nr);
    }
    (*d)[*nr-1]=id;
    (*h)[*nr-1]=ih;

    bDonor[ih] = TRUE;
  }
}

static void search_donors(t_topology *top, int isize, atom_id *index,
			  int *nr_d, atom_id **d, atom_id **h)
{
  int        i,j,max_nr_d;
  t_functype func_type;
  t_ilist    *interaction;
  atom_id    nr1,nr2;
  bool       *bDonor,stop;
  
  snew(bDonor,top->atoms.nr);

  max_nr_d=*nr_d;
  for(func_type=0; func_type < F_NRE; func_type++) {
    interaction=&top->idef.il[func_type];
    for(i=0; i < interaction->nr; i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1 /* next function */) {
      assert(func_type == top->idef.functype[interaction->iatoms[i]]);
      
      /* check out this functype */
      if (func_type == F_SETTLE) {
	nr1=interaction->iatoms[i+1];
	
	if (in_list(nr1,  isize,index)) {
	  if (in_list(nr1+1,isize,index))
	    add_dh(nr1,nr1+1,&max_nr_d,nr_d,d,h,bDonor);
	  if (in_list(nr1+2,isize,index))
	    add_dh(nr1,nr1+2,&max_nr_d,nr_d,d,h,bDonor);
	}
      } else if (IS_CHEMBOND(func_type)) {
	for (j=0; j<2; j++) {
	  nr1=interaction->iatoms[i+1+j];
	  nr2=interaction->iatoms[i+2-j];
	  if ( ( *top->atoms.atomname[nr1][0] == 'H' ) && 
	       ( ( *top->atoms.atomname[nr2][0] == 'O' ) ||
		 ( *top->atoms.atomname[nr2][0] == 'N' )) &&
	       in_list(nr1,isize,index) && in_list(nr2,isize,index))
	    add_dh(nr2,nr1,&max_nr_d,nr_d,d,h,bDonor);
	}
      }
    }
  }
  for(func_type=0; func_type < F_NRE; func_type++) {
    interaction=&top->idef.il[func_type];
    for(i=0; i < interaction->nr; i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1 /* next function */) {
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
	    add_dh(nr2,nr1,&max_nr_d,nr_d,d,h,bDonor);
	}
      }
    }
  }
  sfree(bDonor);
  
  srenew(*d,*nr_d);
  srenew(*h,*nr_d);
}

static void init_grid(bool bBox, matrix box, real rcut, 
		      ivec ngrid, t_gridcell ****grid)
{
  int i,x,y,z;
  
  if (bBox)
    for(i=0; i<DIM; i++)
      ngrid[i]=box[i][i]/(1.2*rcut);
  
  if ( !bBox || (ngrid[XX]<3) || (ngrid[YY]<3) || (ngrid[ZZ]<3) )
    for(i=0; i<DIM; i++)
      ngrid[i]=1;
  if (debug) 
    fprintf(debug,"Will do grid-seach on %dx%dx%d grid\n",
	    ngrid[XX],ngrid[YY],ngrid[ZZ]);
  snew(*grid,ngrid[ZZ]);
  for (z=0; z<ngrid[ZZ]; z++) {
    snew((*grid)[z],ngrid[YY]);
    for (y=0; y<ngrid[YY]; y++)
      snew((*grid)[z][y],ngrid[XX]);
  }
}

char *grpnames[grNR] = {"0D","0H","0A","1D","1H","1A","iD","iH","iA"};

static void build_grid(int *nr, atom_id **a, rvec x[], rvec xshell,
		       bool bBox, matrix box, rvec hbox,
		       real rcut, real rshell,
		       ivec ngrid, t_gridcell ***grid)
{
  int     i,m,gr,xi,yi,zi;
  ivec    grididx;
  rvec    invdelta,dshell;
  t_ncell *newgrid;
  bool    bDoRshell,bInShell;
  real    rshell2=0;
  
  bDoRshell = rshell > 0;
  if (bDoRshell)
    rshell2 = sqr(rshell);
  bInShell = TRUE;
  
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
  for(gr=0; gr<grNR; gr++) {
    /* resetting atom counts */
    for (zi=0; zi<ngrid[ZZ]; zi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (xi=0; xi<ngrid[XX]; xi++)
	  grid[zi][yi][xi][gr].nr=0;
    if (debug)
      fprintf(debug,"Putting %d %s atoms into grid\n",nr[gr],grpnames[gr]);
    /* put atoms in grid cells */
    for(i=0; i<nr[gr]; i++) {
      /* check if we are inside the shell */
      /* if bDoRshell=FALSE then bInShell=TRUE always */
      if ( bDoRshell ) {
	bInShell=TRUE;
	rvec_sub(x[a[gr][i]],xshell,dshell);
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
	    while( x[a[gr][i]][m] < 0 ) 
	      rvec_inc(x[a[gr][i]],box[m]);
	    while( x[a[gr][i]][m] >= box[m][m] ) 
	      rvec_dec(x[a[gr][i]],box[m]);
	    /* get grid index of atom */
	    grididx[m]=x[a[gr][i]][m]*invdelta[m];
	  }
	/* add atom to grid cell */
	newgrid=&(grid[grididx[ZZ]][grididx[YY]][grididx[XX]][gr]);
	if (newgrid->nr >= newgrid->maxnr) {
	  newgrid->maxnr+=10;
	  srenew(newgrid->atoms, newgrid->maxnr);
	}
	newgrid->atoms[newgrid->nr]=i;
	newgrid->nr++;
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
	for (xi=0; xi<ngrid[XX]; xi++)
	  danr[gr] += grid[zi][yi][xi][gr].nr;
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
	  fprintf(fp,"%3d",grid[x][y][z][gr].nr);
	  sum[gr]+=grid[z][y][x][gr].nr;
	}
	fprintf(fp," | ");
	if ( (z+1) < ngrid[ZZ] )
	  for (x=0; x<ngrid[XX]; x++) {
	    fprintf(fp,"%3d",grid[z+1][y][x][gr].nr);
	    sum[gr]+=grid[z+1][y][x][gr].nr;
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

static pbc_correct(rvec dx,matrix box,rvec hbox)
{
  int m;
  
  for(m=DIM-1; m>=0; m--) {
    if ( dx[m] < -hbox[m] )
      rvec_inc(dx,box[m]);
    else if ( dx[m] >= hbox[m] )
      rvec_dec(dx,box[m]);
  }
}

static bool is_hbond(atom_id d, atom_id h, atom_id a, 
		     real rcut, real ccut, 
		     rvec x[], bool bBox, matrix box,rvec hbox,
		     real *d_ha, real *ang,bool bDA)
{
  rvec dist,r_dh,r_ha;
  real d2,ca;
  int m;

  if (d == a)
    return FALSE;
  if (bDA)
    rvec_sub(x[a],x[d],dist);
  else 
    rvec_sub(x[a],x[h],dist);

  if (bBox) 
    pbc_correct(dist,box,hbox);    

  d2 = iprod(dist,dist);
  if ( d2 <= rcut*rcut ) {
    rvec_sub(x[d],x[h],r_dh);
    rvec_sub(x[h],x[a],r_ha);
    if (bBox) 
      pbc_correct(r_dh,box,hbox);    
    if (bBox) 
      pbc_correct(r_ha,box,hbox);    

    ca = cos_angle(r_dh,r_ha);
    /* if angle is smaller, cos is larger */
    if (ca >= ccut) {
      *d_ha = sqrt(d2);
      *ang = acos(ca);
      return TRUE;
    }
  }
  return FALSE;
}

static void sort_dha(int nr_d, atom_id *d, atom_id *h, int nr_a, atom_id *a)
{
  int i,j;
  atom_id temp;
  
  for(i=0; i<nr_d; i++)
    for(j=i+1; j<nr_d; j++)
      if ( (d[i] > d[j]) || 
	   ( (d[i] == d[j]) && (h[i] > h[j]) ) ) {
	temp=d[i];
	d[i]=d[j];
	d[j]=temp;
	temp=h[i];
	h[i]=h[j];
	h[j]=temp;
      }
  for(i=0; i<nr_a; i++)
    for(j=i+1; j<nr_a; j++)
      if (a[i] > a[j]) {
	temp=a[i];
	a[i]=a[j];
	a[j]=temp;
      }
}

static void sort_hb(int *nr_a, t_donor **donors)
{
  int gr,i,j,k;
  t_acceptor ta;
  
  fprintf(stderr,"Sorting hydrogen bonds\n");
  for (gr=0; gr<grNR; gr+=grINC)
    for (i=0; i<nr_a[gr+grD]; i++)
      for (j=0; j<donors[gr+grD][i].nrhb; j++)
	for (k=j+1; k<donors[gr+grD][i].nrhb; k++)
	  if (donors[gr+grD][i].hb[j].a > donors[gr+grD][i].hb[k].a) {
	    ta                      = donors[gr+grD][i].hb[k];
	    donors[gr+grD][i].hb[k] = donors[gr+grD][i].hb[j];
	    donors[gr+grD][i].hb[j] = ta;
	  }
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

static void do_hbac(char *fn,unsigned int **hbexist,int nrhb,int nframes,
		    real time[])
{
  FILE *fp;
  int  i,j,j0;
  real *rhbex;
  real *ct,tail;
  int  *nct;
  
  /* build hbexist matrix in reals for autocorr */
  fprintf(stderr,"Will allocate %.2f Mb of memory\n",
	  nrhb*nframes*sizeof(real)/(1024.0*1024.0));
  snew(rhbex,nframes);
  snew(ct,nframes);
  snew(nct,nframes);
  for(i=0; i<nrhb; i++) {
    j0=0;
    /* while ((j0 < nframes) && (!is_hb(hbexist[i],j0)))
       j0++; */
    for(j=j0; (j<nframes); j++)
      rhbex[j-j0]=is_hb(hbexist[i],j);
      
    low_do_autocorr(NULL,NULL,
		    nframes-j0,1,-1,&rhbex,time[1]-time[0],eacNormal,1,
		    TRUE,TRUE,TRUE,FALSE,0,-1,0,1);
    for(j=0; (j<(nframes-j0)/2); j++) {
      ct[j] += rhbex[j];
      nct[j]++;
    }
  }
  sfree(rhbex);

  /* Normalize */
  for(j=0; ((j<nframes/2) && (nct[j] > 0)); j++)
    ct[j] /= nct[j];
  nframes = j;
  
  /* Determine final value for normalizing */
  tail = 0;
  for(j=nframes/4; (j<nframes/2); j++)
    tail += ct[j];
  tail /= (nframes/2 - nframes/4);
  
  printf("Tail value is %g\n",tail);
  fp = xvgropen(fn, "Hydrogen Bond Autocorrelation","Time (ps)","C(t)");
  for(j=0; (j<nframes); j++)
    fprintf(fp,"%10g  %10g  %10d\n",time[j]-time[0],ct[j],nct[j]);
  /*(ct[j]-tail)/(nrhb-tail));*/
  fclose(fp);
  do_view(fn,NULL);
  sfree(ct);
  sfree(nct);
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
    "[TT]-da[tt]: write out the number of donors and acceptors analyzed for",
    "each timeframe. This is especially usefull when using [TT]-shell[tt]."
  };
  
  static real acut=30, abin=1, rcut=0.35, rbin=0.005, rshell=-1;
  static bool bNitAcc=TRUE,bInsert=FALSE,bDA=TRUE;
  /* options */
  t_pargs pa [] = {
    { "-ins",  FALSE, etBOOL, {&bInsert},
      "Analyze solvent insertion" },
    { "-a",    FALSE, etREAL, {&acut},
      "Cutoff angle (degrees, Donor - Hydrogen - Acceptor)" },
    { "-r",    FALSE, etREAL, {&rcut},
      "Cutoff radius (nm, X - Acceptor, see next option)" },
    { "-dh",   FALSE, etBOOL, {&bDA},
      "Use distance Donor-Acceptor (if TRUE) or Hydrogen-Acceptor (FALSE)" },
    { "-abin", FALSE, etREAL, {&abin},
      "Binwidth angle distribution (degrees)" },
    { "-rbin", FALSE, etREAL, {&rbin},
      "Binwidth distance distribution (nm)" },
    { "-nitacc",FALSE,etBOOL, {&bNitAcc},
      "Regard nitrogen atoms as acceptors" },
    { "-shell", FALSE,etREAL, {&rshell},
      "when > 0, only calculate hydrogen bonds within # nm shell around "
      "one particle" }
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
    { efXVG, "-da",  "danum",  ffOPTWR }
  };
#define NFILE asize(fnm)
  
#define FRINC 8192
#define HBINC 100
  
#define NRGRPS 3
#define OGRP (gr1-grp)
#define INSGRP 2
  
#define HB_NO 0
#define HB_YES 1<<0
#define HB_INS 1<<1
#define HB_YESINS HB_YES|HB_INS
#define HB_NR (1<<2)
  char  hbmap [HB_NR]={ ' ',    'o',      '-',       '*' };
  char *hbdesc[HB_NR]={ "None", "Present", "Inserted", "Present & Inserted" };
  t_rgb hbrgb [HB_NR]={ {1,1,1},{1,0,0},   {0,0,1},    {1,0,1} };
  
  int     status;
  t_topology top;
  t_inputrec ir;
  int     natoms,nframes,max_nframes,shatom;
  int     *isize;
  char    **grpnames;
  atom_id **index;
  rvec    *x,hbox;
  matrix  box;
  real    t,ccut,dist,ang;
  real    *time;
  int     i,j,k,l,start,end;
  int     xi,yi,zi,ai;
  int     xj,yj,zj,aj,xjj,yjj,zjj;
  int     xk,yk,zk,ak,xkk,ykk,zkk;
  bool    bSelected,bStop,bTwo,bHBMap,new,was,bBox,bTric,bDAnr;
  bool    *insert=NULL;
  int     nr_a[grNR];
  atom_id *a[grNR];
  int     grp;
  int     *nhb,*adist,*rdist;
  t_hx    *nhx;
  int     max_nrhb,nrhb,nabin,nrbin,bin,resdist,idx;
  unsigned int **hbexist;
  FILE       *fp,*fpins=NULL,*fplog;
  t_gridcell ***grid;
  t_ncell    *icell,*jcell,*kcell;
  t_icell    *danr;
  ivec       ngrid;
  t_donor    *donors[grNR];
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_BE_NICE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  /* process input */
  bHBMap = opt2bSet("-hbm",NFILE,fnm) || opt2bSet("-ac",NFILE,fnm);
  bDAnr  = opt2bSet("-da",NFILE,fnm);
  bSelected = opt2bSet("-sel",NFILE,fnm);
  ccut   = cos(acut*DEG2RAD);
  
  /* get topology */
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&i,&t,&t,
	   &ir,box,&natoms,NULL,NULL,NULL,&top);
  
  /* initialize h-bond atom groups: */
  for (i=gr0; i<grNR; i+=grINC) {
    nr_a[i+grD] = nr_a[i+grH] = nr_a[i+grA] = 0;
    a[i+grD]    = a[i+grH]    = a[i+grA]    = NULL;
  }
  
  snew(grpnames,NRGRPS);
  snew(index,NRGRPS);
  snew(isize,NRGRPS);
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
    for(grp=gr0; grp<gr0+grINC; grp++) {
      nr_a[grp]=isize[0]/3;
      snew(a[grp],nr_a[grp]);
      for(i=0; i<nr_a[grp]; i++)
	a[grp][i]=index[0][i*3+grp];
    }
    fprintf(stderr,"Analyzing %d selected hydrogen bonds from '%s'\n",
	    nr_a[gr0D],grpnames[0]);
  } else {
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
    if (bTwo)
      fprintf(stderr,"Calculating hydrogen bonds "
	      "between two groups of %d and %d atoms\n",
	      isize[0],isize[1]);
    else
      fprintf(stderr,"Calculating hydrogen bonds in one group of %d atoms\n",
	      isize[0]);
  }
  if (bInsert) {
    fprintf(stderr,"Specify group for insertion analysis:\n");
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      1,&(isize[INSGRP]),&(index[INSGRP]),&(grpnames[INSGRP]));
    fprintf(stderr,"Checking for overlap...\n");
    for (i=0; i<isize[INSGRP]; i++)
      for (grp=0; grp<(bTwo?2:1); grp++)
	for (j=0; j<isize[grp]; j++)
	  if (index[INSGRP][i] == index[grp][j]) 
	    fatal_error(0,"Partial overlap between groups '%s' and '%s'",
			grpnames[grp],grpnames[INSGRP]);
    fpins=ffopen("insert.dat","w");
    fprintf(fpins,"%4s: %15s -> %15s (%7s) - %15s (%7s)\n",
	    "time","insert","donor","distang","acceptor","distang");
  }
  
  /* search donors and acceptors in groups */
  for (i=gr0; i<grNR; i+=grINC)
    if ( ((i==gr0) && !bSelected ) ||
	 ((i==gr1) && bTwo ) ||
	 ((i==grI) && bInsert ) ) {
      search_acceptors(&top, isize[i/grINC], index[i/grINC],
		       &nr_a[i+grA], &a[i+grA], bNitAcc);
      search_donors   (&top, isize[i/grINC], index[i/grINC],
		       &nr_a[i+grH], &a[i+grD], &a[i+grH]);
      nr_a[i+grD] = nr_a[i+grH];
      fprintf(stderr,"Found %d donors and %d acceptors in group '%s'\n",
	      nr_a[i+grD], nr_a[i+grA], grpnames[i/grINC]);
      snew(donors[i+grD], nr_a[i+grD] );
      sort_dha(nr_a[i+grD], a[i+grD], a[i+grH], nr_a[i+grA], a[i+grA]);
    }
  if (bSelected)
    snew(donors[gr0D], nr_a[gr0D]);
  if (!bTwo) {
    for(i=0; i<grINC; i++) {
      nr_a[gr1+i]   = nr_a[gr0+i];
      a[gr1+i]      = a[gr0+i];
    }
    donors[gr1D] = donors[gr0D];
  }
  
  /* check input */
  bStop=FALSE;
  for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC)
    if (nr_a[grp+grD]+nr_a[grp+grA]==0) {
      fprintf(stderr,"No Donors or Acceptors found in group '%s'\n",
	      grpnames[grp/grINC]);
      bStop=TRUE;
    }
  if (!bStop) {
    if (nr_a[gr0D]+nr_a[gr1D]==0) {
      fprintf(stderr,"No Donors found\n");
      bStop=TRUE;
    }
    if (nr_a[gr0A]+nr_a[gr1A]==0) {
      fprintf(stderr,"No Acceptors found\n");
      bStop=TRUE;
    }
  }
  if (bStop)
    fatal_error(0,"Nothing to be done");
  if ( bInsert && ((nr_a[grID]==0) || (nr_a[grIA]==0)) )
    fatal_error(0,"No %s%s%s found in insertion group '%s'\n",
		(nr_a[grID]==0)?"donors":"",
		((nr_a[grID]==0)&&(nr_a[grIA]==0))?" and ":"",
		(nr_a[grIA]==0)?"acceptors":"",
		grpnames[grI/grINC]);
  
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

  /* analyze trajectory */
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( natoms > top.atoms.nr )
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);
  bBox = ir.ePBC!=epbcNONE;
  init_grid(bBox, box, rcut, ngrid, &grid);
  max_nframes = nframes = 0;
  max_nrhb    = nrhb    = 0;
  hbexist     = NULL;
  nhb         = NULL;
  nhx         = NULL;
  time        = NULL;
  danr        = NULL;
  nabin       = acut/abin;
  nrbin       = rcut/rbin;
  snew(adist,nabin);
  snew(rdist,nrbin);
  if (bInsert)
    snew(insert,natoms);
  do {
    bTric = bBox && TRICLINIC(box);
    build_grid(nr_a,a, x,x[shatom], bBox,box,hbox, rcut, rshell, ngrid,grid);
    if (debug)
      dump_grid(debug, ngrid, grid);
    
    if (nframes >= max_nframes) {
      max_nframes += FRINC;
      srenew(nhb,max_nframes);
      srenew(nhx,max_nframes);
      srenew(time,max_nframes);
      if (bHBMap)
	for (i=0; i<max_nrhb; i++)
	  srenew(hbexist[i],max_nframes/32);
      if (bDAnr)
	srenew(danr,max_nframes);
    }
    time[nframes] = t;
    nhb[nframes]  = 0;
    for (i=0; i<max_hx; i++)
      nhx[nframes][i]=0;
    if (bHBMap)
      for (i=0; i<max_nrhb; i++)
	set_hb(hbexist[i],nframes,HB_NO);
    if (bDAnr)
      count_da_grid(ngrid, grid, danr[nframes]);
    /* loop over all gridcells (xi,yi,zi)      */
    /* Removed confusing macro, DvdS 27/12/98  */
    for(xi=0; (xi<ngrid[XX]); xi++)
      for(yi=0; (yi<ngrid[YY]); yi++)
	for(zi=0; (zi<ngrid[ZZ]); zi++) {
	  /* loop over groups gr0 (always) and gr1 (if necessary) */
	  for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC) {
	    icell=&grid[zi][yi][xi][grp+grH];
	    /* loop over all hydrogen atoms from group (grp) 
	     * in this gridcell (icell) 
	     */
	    for (ai=0; ai<icell->nr; ai++) {
	      i=icell->atoms[ai];
	      /* loop over all adjacent gridcells (xj,yj,zj) */
	      /* This is a macro!!! */
	      LOOPGRIDINNER(xj,yj,zj,xjj,yjj,zjj,xi,yi,zi,ngrid,bTric) {
		jcell=&grid[zj][yj][xj][OGRP+grA];
		/* loop over acceptor atoms from other group (OGRP) 
		 * in this adjacent gridcell (jcell) 
		 */
		for (aj=0; aj<jcell->nr; aj++) {
		  j=jcell->atoms[aj];
		  if ( (bSelected && (j==i)) || (!bSelected) ) {
		    /* check if this once was a h-bond */
		    idx=NOTSET;
		    for (k=0; (k<donors[grp][i].nrhb) && (idx==NOTSET); k++)
		      if (j == donors[grp][i].hb[k].a)
			idx=k;
		    if ( is_hbond(a[ grp+grD][i],a[ grp+grH][i],a[OGRP+grA][j],
				  rcut,ccut,x,bBox,box,hbox,&dist,&ang,bDA) ) {
		      /* add to index if not already there */
		      if (idx==NOTSET) {
			if (donors[grp][i].nrhb>=donors[grp][i].maxnr) {
			  donors[grp][i].maxnr+=10;
			  srenew(donors[grp][i].hb,donors[grp][i].maxnr);
			}
			/* Add a donor atom in a hbond */
			donors[grp][i].hb[donors[grp][i].nrhb].a=j;
			donors[grp][i].hb[donors[grp][i].nrhb].nr=nrhb;
			idx = donors[grp][i].nrhb;
			donors[grp][i].nrhb++;
			if (bHBMap && (nrhb >= max_nrhb) ) {
			  max_nrhb+=HBINC;
			  srenew(hbexist,max_nrhb);
			  for (l=max_nrhb-HBINC; l<max_nrhb; l++)
			    snew(hbexist[l],max_nframes/32);
			}
			nrhb++;
		      }
		      if (bHBMap)
			/* update matrix */
			set_hb(hbexist[donors[grp][i].hb[idx].nr],
			       nframes,HB_YES);
		      
		      /* count number of hbonds per frame */
		      nhb[nframes]++;
		      
		      /* make angle and distance distributions */
		      ang*=RAD2DEG;
		      adist[(int)( ang/abin)]++;
		      rdist[(int)(dist/rbin)]++;
		      
		      if (!bTwo) {
			resdist=abs(top.atoms.atom[a[ grp+grD][i]].resnr-
				    top.atoms.atom[a[OGRP+grA][j]].resnr);
			if (resdist >= max_hx) 
			  resdist = max_hx-1;
			nhx[nframes][resdist]++;
		      }
		    } 
		    if (bInsert && ( (idx!=NOTSET) || bSelected ) ) {
		      /* this has been a h-bond, or we are analyzing 
			 selected bonds: check for inserted */
		      bool ins_d, ins_a;
		      real ins_d_dist, ins_d_ang, ins_a_dist, ins_a_ang;
		      int  ins_d_k=0,ins_a_k=0;
		      
		      ins_d=ins_a=FALSE;
		      ins_d_dist=ins_d_ang=ins_a_dist=ins_a_ang=1e6;
		      
		      /* loop over gridcells adjacent to i (xk,yk,zk) */
		      LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xi,yi,zi,ngrid,bTric){
			kcell=&grid[zk][yk][xk][grIA];
			/* loop over acceptor atoms from ins group 
			   in this adjacent gridcell (kcell) */
			for (ak=0; ak<kcell->nr; ak++) {
			  k=kcell->atoms[ak];
			  if (is_hbond(a[grp+grD][i],
				       a[grp+grH][i],
				       a[   grIA][k],
				       rcut,ccut,x,bBox,box,hbox,
				       &dist,&ang,bDA))
			    if (dist<ins_d_dist) {
			      ins_d=TRUE;
			      ins_d_dist=dist;
			      ins_d_ang =ang ;
			      ins_d_k   =k   ;
			    }
			}
		      }
		      ENDLOOPGRIDINNER;
		      /* loop over gridcells adjacent to j (xk,yk,zk) */
		      LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xj,yj,zj,ngrid,bTric){
			kcell=&grid[zk][yk][xk][grIH];
			/* loop over hydrogen atoms from ins group 
			   in this adjacent gridcell (kcell) */
			for (ak=0; ak<kcell->nr; ak++) {
			  k=kcell->atoms[ak];
			  if (is_hbond(a[    grID][k],
				       a[    grIH][k],
				       a[OGRP+grA][j],
				       rcut,ccut,x,bBox,box,hbox,
				       &dist,&ang,bDA))
			    if (dist<ins_a_dist) {
			      ins_a=TRUE;
			      ins_a_dist=dist;
			      ins_a_ang =ang ;
			      ins_a_k   =k   ;
			    }
			}
		      }
		      ENDLOOPGRIDINNER;
		      
		      if (ins_d && ins_a && 
			  (a[grIA][ins_d_k] == a[grID][ins_a_k])) {
			/* set insertion index */
			insert[a[grID][ins_a_k]]=TRUE;
			insert[a[grIH][ins_a_k]]=TRUE;
			insert[a[grIA][ins_d_k]]=TRUE;
			
			/* add to hbond index if not already there */
			if (idx==NOTSET) {
			  if (donors[grp][i].nrhb>=donors[grp][i].maxnr) {
			    donors[grp][i].maxnr+=10;
			    srenew(donors[grp][i].hb,donors[grp][i].maxnr);
			  }
			  donors[grp][i].hb[donors[grp][i].nrhb].a=j;
			  donors[grp][i].hb[donors[grp][i].nrhb].nr=nrhb;
			  idx=donors[grp][i].nrhb;
			  donors[grp][i].nrhb++;
			  if (bHBMap && (nrhb>=max_nrhb) ) {
			    max_nrhb+=HBINC;
			    srenew(hbexist,max_nrhb);
			    for (l=max_nrhb-HBINC; l<max_nrhb; l++)
			      snew(hbexist[l],max_nframes/32);
			  }
			  nrhb++;
			}
			
			/* mark insertion in hbond index */
			if (bHBMap)
			  set_hb(hbexist[donors[grp][i].hb[idx].nr],
				 nframes,HB_INS);
			
			/* print insertion info to file */
			fprintf(fpins,
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
				ins_a_dist,ins_a_ang*RAD2DEG);
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
  
  if (nrhb==0)
    fprintf(stderr,"No hydrogen bonds found!!\n");
  else {
    fprintf(stderr,"Found %d different hydrogen bonds in trajectory\n",nrhb);
    
    /* Dump everything to output */
    sort_hb(nr_a,donors);
    
    fp = xvgropen(opt2fn("-num",NFILE,fnm),"Hydrogen Bonds","Time","Number");
    for(i=0; i<nframes; i++)
      fprintf(fp,"%10g %10d\n",time[i],nhb[i]);
    fclose(fp);
    
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
	fprintf(fp,"%10g",time[i]);
	for (j=0; j<max_hx; j++)
	  fprintf(fp," %6d",nhx[i][j]);
	fprintf(fp,"\n");
      }
      fclose(fp);
    }
    
    if (opt2bSet("-ac",NFILE,fnm))
      do_hbac(opt2fn("-ac",NFILE,fnm),hbexist,nrhb,nframes,time);
    
    if (opt2bSet("-hbm",NFILE,fnm)) {
      t_matrix mat;
      int x,y;
      
      mat.nx=nframes;
      mat.ny=nrhb;
      snew(mat.matrix,mat.nx);
      for(x=0; x<mat.nx; x++){
	snew(mat.matrix[x],mat.ny);
	y=0;
	for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC)
	  for(i=0; i<nr_a[grp+grD]; i++)
	    for(j=0; j<donors[grp][i].nrhb; j++)
	      mat.matrix[x][y++]=hbexist[donors[grp][i].hb[j].nr][x];
      }
      mat.axis_x=time;
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
    for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC) {
      fprintf(fp,"[ %s ]",grpnames[grp/grINC]);
      for (i=0; i<isize[grp/grINC]; i++) {
	fprintf(fp,(i%15)?" ":"\n");
	fprintf(fp,"%4u",index[grp/grINC][i]+1);
      }
      fprintf(fp,"\n");
      fprintf(fp,"[ donors_hydrogens_%s ]",grpnames[grp/grINC]);
      for (i=0; i<nr_a[grp+grD]; i++) {
	fprintf(fp,(i%6)?"   ":"\n");
	fprintf(fp,"%4u %4u",a[grp+grD][i]+1,a[grp+grH][i]+1);
      }
      fprintf(fp,"\n");
      fprintf(fp,"[ acceptors_%s ]",grpnames[grp/grINC]);
      for (i=0; i<nr_a[grp+grA]; i++) {
	fprintf(fp,(i%15)?" ":"\n");
	fprintf(fp,"%4u",a[grp+grA][i]+1);
      }
      fprintf(fp,"\n");
    }
    if (bTwo)
      fprintf(fp,"[ hbonds_%s-%s ]\n",grpnames[0],grpnames[1]);
    else
      fprintf(fp,"[ hbonds_%s ]\n",grpnames[0]);
    for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC)
      for(i=0; i<nr_a[grp+grD]; i++)
	for(j=0; j<donors[grp][i].nrhb; j++) {
	  int  ddd,hhh,aaa;
	  char ds[32],hs[32],as[32];
	  
	  ddd = a[ grp+grD][i];
	  hhh = a[ grp+grH][i];
	  aaa = a[OGRP+grA][donors[grp][i].hb[j].a];
	  sprintf(ds,"%s",mkatomname(&top.atoms,ddd));
	  sprintf(hs,"%s",mkatomname(&top.atoms,hhh));
	  sprintf(as,"%s",mkatomname(&top.atoms,aaa));
	  fprintf(fp,"%6u %6u %6u\n",ddd+1,hhh+1,aaa+1);
	  if (fplog) 
	    fprintf(fplog,"%12s  %12s  %12s\n",ds,hs,as);
	}
    if (bInsert) {
      if (bTwo)
	fprintf(fp,"[ insert_%s->%s-%s ]",
		grpnames[2],grpnames[0],grpnames[1]);
      else
	fprintf(fp,"[ insert_%s->%s ]",grpnames[2],grpnames[0]);
      j=0;
      for(i=0; i<natoms; i++) 
	if (insert[i]) {
	  fprintf(fp,(j%15)?" ":"\n");
	  fprintf(fp,"%4d",i+1);
	  j++;
	}
      fprintf(fp,"\n");
    }
    fclose(fp);
    if (fplog)
      fclose(fplog);
  }
  
  if (bDAnr) {
    int  i,j,nleg;
    char **legnames;
    char buf[STRLEN];
    char *danames[grINC] = {"Donors", NULL, "Acceptors"};
    
#define USE_THIS_GROUP(j) ( ( j % grINC != grH ) && ( bTwo || ( j / grINC != 1 ) ) && ( bInsert || ( j / grINC != 2 ) ) )
    
    fp = xvgropen(opt2fn("-da",NFILE,fnm),
		  "Donors and Acceptors","Time(ps)","Count");
    nleg = (bTwo?2:1)*2 + (bInsert?2:0);
    snew(legnames,nleg);
    i=0;
    for(j=0; j<grNR; j++)
      if ( USE_THIS_GROUP(j) ) {
	sprintf(buf,"%s %s",danames[j % grINC],grpnames[j / grINC]);
	legnames[i++] = strdup(buf);
      }
    assert(i==nleg);
    xvgr_legend(fp,nleg,legnames);
    for(i=0; i<nframes; i++) {
      fprintf(fp,"%10g",time[i]);
      for (j=0; j<grNR; j++)
	if ( USE_THIS_GROUP(j) )
	  fprintf(fp," %6d",danr[i][j]);
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  thanx(stderr);
  
  return 0;
}
