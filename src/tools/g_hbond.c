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
 * GRowing Old MAkes el Chrono Sweat
 */

static char *SRCID_g_hbond_c = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) g_hbond.cc 1.29 9/30/97"
#endif /* HAVE_IDENT */

#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "futil.h"
#include "tpxio.h"
#include "physics.h"
#include "macros.h"
#include "fatal.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "assert.h"
#include "vec.h"
#include "xvgr.h"
#include "gstat.h"
#include "matio.h"
#include <math.h>

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
#define grINC gr1-gr0

typedef struct {
  int nr;
  int maxnr;
  atom_id *atoms;
} t_gridcell;

typedef struct {
  int a;
  int nr;
} t_acceptor;

typedef struct {
  int nrhb;
  int maxnr;
  t_acceptor *hb;
} t_donor;

bool in_list(atom_id selection,int isize,atom_id *index)
{
  int i;
  bool bFound;
  
  bFound=FALSE;
  for(i=0;(i<isize) && !bFound;i++)
    if(selection == index[i])
      bFound=TRUE;
  
  return bFound;
}

void add_acc(int i, int *max_nr,int *nr,atom_id **a)
{
  (*nr)++;
  if ((*nr) > (*max_nr)) {
    (*max_nr)+=10;
    srenew((*a),(*max_nr));
  }
  (*a)[(*nr)-1]=i;
}

void search_acceptors(t_topology *top, int isize, atom_id *index,
		      int *nr_a, atom_id **a, bool bNitAcc)
{
  int i,max_nr_a;
  
  max_nr_a=*nr_a;
  for (i=0;(i<top->atoms.nr);i++)
    if (*top->atoms.atomname[i][0] == 'O' || 
	(bNitAcc && (*top->atoms.atomname[i][0] == 'N')))
      if (in_list(i,isize,index))
	add_acc(i,&max_nr_a,nr_a,a);
  srenew((*a),(*nr_a));
}

void add_dh(int id, int ih, int *max_nr,int *nr,atom_id **d,atom_id **h)
{
  (*nr)++;
  if ((*nr) > (*max_nr)) {
    (*max_nr)+=10;
    srenew((*d),(*max_nr));
    srenew((*h),(*max_nr));
  }
  (*d)[(*nr)-1]=id;
  (*h)[(*nr)-1]=ih;
}

void search_donors(t_topology *top, int isize, atom_id *index,
		   int *nr_d, atom_id **d, atom_id **h, bool bDumConn)
{
  int i,j,max_nr_d;
  t_functype func_type;
  t_ilist *interaction;
  atom_id nr1,nr2;
  bool stop;
  
  max_nr_d=*nr_d;
  for(func_type=0; func_type < F_NRE; func_type++) {
    interaction=&top->idef.il[func_type];
    for(i=0; i<interaction->nr; i+=interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1 /* next function */) {
      assert(func_type == top->idef.functype[interaction->iatoms[i]]);
      
      /* check out this functype */
      if (func_type == F_SETTLE) {
	nr1=interaction->iatoms[i+1];
	
	if (in_list(nr1,  isize,index) && 
	    in_list(nr1+1,isize,index) && 
	    in_list(nr1+2,isize,index) ) {
	  add_dh(nr1,nr1+1,&max_nr_d,nr_d,d,h);
	  add_dh(nr1,nr1+2,&max_nr_d,nr_d,d,h);
	}
      } else if ( interaction_function[func_type].flags & IF_CONNECT ) {
	for (j=0; j<2; j++) {
	  nr1=interaction->iatoms[i+1+j];
	  nr2=interaction->iatoms[i+2-j];
	  if ( ((*top->atoms.atomname[nr1][0]) == 'H') && 
	       ( (*top->atoms.atomname[nr2][0] == 'O') ||
		 (*top->atoms.atomname[nr2][0] == 'N')) &&
	       in_list(nr1,isize,index) && in_list(nr2,isize,index))
	    add_dh(nr2,nr1,&max_nr_d,nr_d,d,h);
	}
      } else if ( bDumConn &&
		  (interaction_function[func_type].flags & IF_DUMMY) ) {
	nr1=interaction->iatoms[i+1];
	if ((*top->atoms.atomname[nr1][0]) == 'H') {
	  nr2=nr1-1;
	  stop=FALSE;
	  while (!stop && ((*top->atoms.atomname[nr2][0]) == 'H'))
	    if (nr2)
	      nr2--;
	    else
	      stop=TRUE;
	  if ( !stop && ( (*top->atoms.atomname[nr2][0] == 'O') ||
			  (*top->atoms.atomname[nr2][0] == 'N') ) &&
	       in_list(nr1,isize,index) && in_list(nr2,isize,index) )
	    add_dh(nr2,nr1,&max_nr_d,nr_d,d,h);
	}
      }
    }
  }
  srenew((*d),(*nr_d));
  srenew((*h),(*nr_d));
}

void init_grid(bool bBox, matrix box, real rcut, 
	       ivec ngrid, t_gridcell *****grid)
{
  int i,x,y,z;
  
  if (bBox)
    for(i=0; i<DIM; i++)
      ngrid[i]=box[i][i]/(1.2*rcut);
  
  if ( !bBox || (ngrid[XX]<3) || (ngrid[YY]<3) || (ngrid[ZZ]<3) )
    for(i=0; i<DIM; i++)
      ngrid[i]=1;
  if (debug) fprintf(stderr,"\nWill do grid-seach on %dx%dx%d grid\n",
		     ngrid[XX],ngrid[YY],ngrid[ZZ]);
  snew(*grid,ngrid[XX]);
  for (x=0; x<ngrid[XX]; x++) {
    snew((*grid)[x],ngrid[YY]);
    for (y=0; y<ngrid[YY]; y++) {
      snew((*grid)[x][y],ngrid[ZZ]);
      for (z=0; z<ngrid[ZZ]; z++)
	snew((*grid)[x][y][z],grNR);
    }
  }
}

void build_grid(int *nr, atom_id **a, rvec x[], 
		bool bBox, matrix box, rvec hbox, real rcut,
		ivec ngrid, t_gridcell ****grid)
{
  int i,m,gr,xi,yi,zi;
  ivec grididx;
  rvec invdelta;
  t_gridcell *newgrid;
  
  for(m=0; m<DIM; m++) {
    hbox[m]=box[m][m]*0.5;
    if (bBox) {
      invdelta[m]=ngrid[m]/box[m][m];
      if (1/invdelta[m] < rcut)
	fatal_error(0,"box shrank too much to keep using this grid\n");
    } else
      invdelta[m]=0;
  }
  for(gr=0; gr<grNR; gr++) {
    /* resetting atom counts */
    for (xi=0; xi<ngrid[XX]; xi++)
      for (yi=0; yi<ngrid[YY]; yi++)
	for (zi=0; zi<ngrid[ZZ]; zi++)
	  grid[xi][yi][zi][gr].nr=0;
    /* put atoms in grid cells */
    for(i=0; i<nr[gr]; i++) {
      for(m=0; m<DIM; m++) {
	/* put atom in the box */
	if (bBox) {
	  while(x[a[gr][i]][m]<0) 
	    x[a[gr][i]][m]+=box[m][m];
	  while(x[a[gr][i]][m]>=box[m][m]) 
	    x[a[gr][i]][m]-=box[m][m];
	}
	/* get grid index of atom */
	grididx[m]=x[a[gr][i]][m]*invdelta[m];
      }
      /* add atom to grid cell */
      newgrid=&(grid[grididx[XX]][grididx[YY]][grididx[ZZ]][gr]);
      if (newgrid->nr >= newgrid->maxnr) {
	newgrid->maxnr+=10;
	srenew(newgrid->atoms, newgrid->maxnr);
      }
      newgrid->atoms[newgrid->nr]=i;
      newgrid->nr++;
    }
  }
}

#define LOOPGRIDOUTER(x,y,z,n)\
  for(x=0; x<n[XX]; x++)\
    for(y=0; y<n[YY]; y++)\
      for(z=0; z<n[ZZ]; z++)
    
#define B(n,x)((n!=1)?(x-1):x)
#define E(n,x)((n!=1)?(x+1):x)
#define GRIDMOD(j,n) (j+n)%(n)
#define LOOPGRIDINNER(x,y,z,xx,yy,zz,xo,yo,zo,n)\
     for(xx=B(n[XX],xo); xx<=E(n[XX],xo); xx++) {\
       x=GRIDMOD(xx,n[XX]);\
       for(yy=B(n[YY],yo); yy<=E(n[YY],yo); yy++) {\
	 y=GRIDMOD(yy,n[YY]);\
	 for(zz=B(n[ZZ],zo); zz<=E(n[ZZ],zo); zz++) {\
	   z=GRIDMOD(zz,n[ZZ]);
#define ENDLOOPGRIDINNER\
	 }\
       }\
     }\

void dump_grid(FILE *fp, ivec ngrid, t_gridcell ****grid)
{
  int gr,x,y,z,sum[grNR];
  
  fprintf(fp,"grid %dx%dx%d\n",ngrid[XX],ngrid[YY],ngrid[ZZ]);
  for (gr=0; gr<grNR; gr++) {
    sum[gr]=0;
    fprintf(fp,"GROUP %d\n",gr);
    for (z=0; z<ngrid[ZZ]; z+=2) {
      fprintf(fp,"Z=%d,%d\n",z,z+1);
      for (y=0; y<ngrid[YY]; y++) {
	for (x=0; x<ngrid[XX]; x++) {
	  fprintf(fp,"%3d",grid[x][y][z][gr].nr);
	  sum[gr]+=grid[x][y][z][gr].nr;
	}
	fprintf(fp," | ");
	if ( (z+1) < ngrid[ZZ] )
	  for (x=0; x<ngrid[XX]; x++) {
	    fprintf(fp,"%3d",grid[x][y][z+1][gr].nr);
	    sum[gr]+=grid[x][y][z+1][gr].nr;
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

void free_grid(ivec ngrid, t_gridcell *****grid)
{
  int x,y,z;
  
  for (x=0; x<ngrid[XX]; x++) {
    for (y=0; y<ngrid[YY]; y++) {
      for (z=0; z<ngrid[ZZ]; z++)
	sfree((*grid)[x][y][z]);
      sfree((*grid)[x][y]);
    }
    sfree((*grid)[x]);
  }
  sfree(*grid);
}

bool is_hbond(atom_id d, atom_id h, atom_id a, 
	      real rcut, real ccut, 
	      rvec x[], bool bBox, rvec hbox, real *d_da, real *ang)
{
  rvec r_da,r_dh;
  real d_da2,ca;
  int m;
  
  for(m=0; m<DIM; m++) {
    r_da[m]=x[d][m]-x[a][m];
    if (bBox)
      if (r_da[m]<-hbox[m]) 
	r_da[m]+=2*hbox[m];
      else if (r_da[m]>=hbox[m]) 
	r_da[m]-=2*hbox[m];
    if ( (r_da[m]>rcut) || (-r_da[m]>rcut) )
      return FALSE;
  }
  d_da2 = iprod(r_da,r_da);
  if ( d_da2 <= rcut*rcut ) {
    rvec_sub(x[d],x[h],r_dh);
    for(m=0; m<DIM; m++) {
      if (bBox)
	if (r_dh[m]<-hbox[m]) 
	  r_dh[m]+=2*hbox[m];
	else if (r_dh[m]>=hbox[m]) 
	  r_dh[m]-=2*hbox[m];
    }
    ca = cos_angle(r_da,r_dh);
    /* if angle is smaller, cos is larger */
    if ( ca>=ccut ) {
      *d_da = sqrt(d_da2);
      *ang = acos(ca);
      return TRUE;
    }
  }
  return FALSE;
}

void sort_dha(int nr_d, atom_id *d, atom_id *h, int nr_a, atom_id *a)
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

void sort_hb(int nrhb, int *nr_a, t_donor **donors)
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_hbond computes and analyzes hydrogen bonds. Hydrogen bonds are",
    "determined based on cutoffs for the angle Hydrogen - Donor - Acceptor",
    "(zero is optimum) and the distance Donor - Acceptor.[PAR]",
    "You need to specify two groups for analysis, which must be either",
    "identical or non-overlapping. All hydrogen bonds between the two",
    "groups are analyzed.[PAR]",
    "It is also possible to analyse specific hydrogen bonds with",
    "[TT]-sel[tt]. This index file must contain a group of atom triplets",
    "Donor Hydrogen Acceptor, in the following way:[PAR]",
    "[TT]",
    "[ selected ][BR]",
    "     20    21    24[BR]",
    "     25    26    29[BR]",
    "      0     3     6[BR]",
    "[tt][BR]",
    "(Note that the triplets need not be on separate lines.)",
    "Each atom triplet specifies a hydrogen bond to be analyzed.",
    "Note that no check is made for the types of atoms.[PAR]",
    "[TT]-ins[tt] turns on computing solvent insertion into hydrogen bonds.",
    "In this case an additional group must be selected, specifying the",
    "solvent molecules.[PAR]",
    "[TT]-dumconn[tt] makes g_hbond assume a covalent bond exists between",
    "any dummy atom and the first preceding (in sequence) heavy atom.",
    "This is used in searching Donor-Hydrogen pairs.[PAR]",
    "The following output files are/can be generated:[BR]",
    "[TT]-ac[tt]:   average over all autocorrelations of the existence",
    "functions (either 0 or 1) of all hydrogen bonds.[BR]",
    "[TT]-num[tt]:  number of hydrogen bonds as a function of time.[BR]",
    "[TT]-dist[tt]: distance distribution of all hydrogen bonds.[BR]",
    "[TT]-ang[tt]:  angle distribution of all hydrogen bonds.[BR]",
    "[BB]Output:[bb][BR]",
    "[TT]-hx[tt]:   the number of n-n+i hydrogen bonds as a function of time",
    "where n and n+i stand for residue numbers and i ranges from 0 to 6.",
    "This includes the n-n+3, n-n+4 and n-n+5 hydrogen bonds associated",
    "with helices in proteins.[BR]",
    "[TT]-hbn[tt]:  all selected groups, donors, hydrogens and acceptors",
    "for selected groups, all hydrogen bonded atoms from all groups and",
    "all solvent atoms involved in insertion.[BR]",
    "[TT]-hbm[tt]:  existence matrix for all hydrogen bonds over all",
    "frames, this also contains information on solvent insertion",
    "into hydrogen bonds.[BR]",
  };
  
  static real acut=60, abin=1, rcut=0.35, rbin=0.005;
  static bool bNitAcc=TRUE,bDumConn=TRUE,bInsert=FALSE;
  /* options */
  t_pargs pa [] = {
    { "-ins",  FALSE, etBOOL, &bInsert,
      "analyze solvent insertion" },
    { "-a",    FALSE, etREAL, &acut,
      "cutoff angle (degrees, Hydrogen - Donor - Acceptor)" },
    { "-abin", FALSE, etREAL, &abin,
      "binwidth angle distribution (degrees)" },
    { "-r",    FALSE, etREAL, &rcut,
      "cutoff radius (nm, Donor - Acceptor)" },
    { "-rbin", FALSE, etREAL, &rbin,
      "binwidth distance distribution (nm)" },
    { "-nitacc",FALSE,etBOOL, &bNitAcc,
      "regard nitrogen atoms as acceptors" },
    { "-dumconn",FALSE,etBOOL, &bDumConn,
      "dummy atoms bond to first preceding heavy atom" }
  };

  t_filenm fnm[] = {
    { efTRX, "-f",   NULL,     ffREAD  },
    { efTPX, NULL,   NULL,     ffREAD  },
    { efNDX, NULL,   NULL,     ffOPTRD },
    { efNDX, "-sel", "select", ffOPTRD },
    { efXVG, "-num", "hbnum",  ffWRITE },
    { efXVG, "-ac",  "hbac",   ffOPTWR },
    { efXVG, "-dist","hbdist", ffOPTWR },
    { efXVG, "-ang", "hbang",  ffOPTWR },
    { efXVG, "-hx",  "hbhelix",ffOPTWR },
    { efNDX, "-hbn", "hbond",  ffOPTWR },
    { efXPM, "-hbm", "hbmap",  ffOPTWR }
  };
#define NFILE asize(fnm)
  
#define FRINC 100
#define HBINC 100
  
#define ngrps 3
#define ogrp (gr1-grp)
#define insgrp grI/grINC
  
#define hbNO 0
#define hbYES 1<<0
#define hbINS 1<<1
#define hbYESINS hbYES|hbINS
#define hbNR (1<<2)
  char  hbmap [hbNR]={ ' ',    'o',      '-',       '*' };
  char *hbdesc[hbNR]={ "None", "Present","Inserted","Present & Inserted" };
  t_rgb hbrgb [hbNR]={ {1,1,1},{1,0,0},  {0,0,1},   {1,0,1} };
  
  int     status;
  t_topology top;
  t_inputrec ir;
  int     natoms,nframes,max_nframes;
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
  bool    bSelected,bStop,bTwo,new,was,bBox;
  bool    *insert;
  int     nr_a[grNR];
  atom_id *a[grNR];
  int     grp;
  int     *nhb,*adist,*rdist;
  t_hx    *nhx;
  int     max_nrhb,nrhb,nabin,nrbin,bin,resdist,idx;
  unsigned char **hbexist;
  FILE    *fp,*fpins;
  t_gridcell ****grid,*icell,*jcell,*kcell;
  ivec    ngrid;
  t_donor *donors[grNR];
  real    **rhbex;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,TRUE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  /* process input */
  bSelected = opt2bSet("-sel",NFILE,fnm);
  ccut = cos(acut*DEG2RAD);
  
  /* get topology */
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&i,&t,&t,
	   &ir,box,&natoms,NULL,NULL,NULL,&top);
  
  if (bSelected) {
    /* analyze selected hydrogen bonds */
    fprintf(stderr,"Select group with selected atoms:\n");
    snew(grpnames,ngrps);
    snew(index,ngrps);
    snew(isize,ngrps);
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
    snew(grpnames,ngrps);
    snew(index,ngrps);
    snew(isize,ngrps);
    get_index(&(top.atoms),ftp2fn_null(efNDX,NFILE,fnm),
	      2,isize,index,grpnames);
    
    /* check if we have two identical or two non-overlapping groups */
    bTwo=!(isize[0] == isize[1]);
    for (i=0; (i<isize[0]) && !bTwo; i++)
      bTwo=(index[0][i] != index[1][i]);
    if (bTwo) {
      fprintf(stderr,"Checking for overlap...\n");
      for (i=0; (i<isize[0]); i++)
	for (j=0; (j<isize[1]); j++)
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
	      1,&(isize[insgrp]),&(index[insgrp]),&(grpnames[insgrp]));
    fprintf(stderr,"Checking for overlap...\n");
    for (i=0; (i<isize[insgrp]); i++)
      for (grp=0; grp<(bTwo?2:1); grp++)
	for (j=0; (j<isize[grp]); j++)
	  if (index[insgrp][i] == index[grp][j]) 
	    fatal_error(0,"Partial overlap between groups '%s' and '%s'",
			grpnames[grp],grpnames[insgrp]);
    fpins=ffopen("insert.dat","w");
    fprintf(fpins,"%4s: %15s -> %15s (%7s) - %15s (%7s)\n",
	    "time","insert","donor","distang","acceptor","distang");
  }
  
  /* search donors and acceptors in groups */
  for (i=gr0; i<grNR; i+=grINC) 
    if ( ((i==gr0) && !bSelected ) ||
	 ((i==gr1) && bTwo ) ||
	 ((i==grI) && bInsert ) ) {
      nr_a[i+grD]=nr_a[i+grA]=0;
      a[i+grD]=a[i+grH]=a[i+grA]=NULL;
      search_acceptors(&top,isize[i/grINC],index[i/grINC],
		       &nr_a[i+grA],&a[i+grA],bNitAcc);
      search_donors(&top,isize[i/grINC],index[i/grINC],
		    &nr_a[i+grD],&a[i+grD],&a[i+grH],bDumConn);
      fprintf(stderr,"Found %d donors and %d acceptors in group '%s'\n",
	      nr_a[i+grD],nr_a[i+grA],grpnames[i/grINC]);
      snew(donors[i+grD],nr_a[i+grD]);
      sort_dha(nr_a[i+grD],a[i+grD],a[i+grH],nr_a[i+grA],a[i+grA]);
    }
  if (bSelected)
    snew(donors[gr0D],nr_a[gr0D]);
  if (!bTwo)
    for(i=0; i<grINC; i++) {
      nr_a[gr1+i]=nr_a[gr0+i];
      a[gr1+i]=a[gr0+i];
      donors[gr1+i]=donors[gr0+i];
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

  /* analyze trajectory */
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( natoms > top.atoms.nr )
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);
  bBox=ir.eBox!=ebtNONE;
  init_grid(bBox, box, rcut, ngrid, &grid);
  max_nframes=nframes=0;
  max_nrhb=nrhb=0;
  hbexist=NULL;
  nhb=NULL;
  nhx=NULL;
  time=NULL;
  nabin=acut/abin;
  nrbin=rcut/rbin;
  snew(adist,nabin);
  snew(rdist,nrbin);
  if (bInsert)
    snew(insert,natoms);
  do {
    build_grid(nr_a, a, x, bBox, box, hbox, rcut, ngrid, grid);
    if (nframes>=max_nframes) {
      max_nframes+=FRINC;
      srenew(nhb,max_nframes);
      srenew(nhx,max_nframes);
      srenew(time,max_nframes);
      for (i=0; i<max_nrhb; i++)
	srenew(hbexist[i],max_nframes);
    }
    time[nframes]=t;
    nhb[nframes]=0;
    for (i=0; i<max_hx; i++)
      nhx[nframes][i]=0;
    for (i=0; i<max_nrhb; i++)
      hbexist[i][nframes]=hbNO;
    /* loop over all gridcells (xi,yi,zi) */
    LOOPGRIDOUTER(xi,yi,zi,ngrid) {
      /* loop over groups gr0 (always) and gr1 (if necessary) */
      for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC) {
	icell=&grid[xi][yi][zi][grp+grD];
	/* loop over al donor atoms from group (grp) 
	   in this gridcell (icell) */
	for (ai=0; ai<icell->nr; ai++) {
	  i=icell->atoms[ai];
	  /* loop over all adjacent gridcells (xj,yj,zj) */
	  LOOPGRIDINNER(xj,yj,zj,xjj,yjj,zjj,xi,yi,zi,ngrid) {
	    jcell=&grid[xj][yj][zj][ogrp+grA];
	    /* loop over acceptor atoms from other group (ogrp) 
	       in this adjacent gridcell (jcell) */
	    for (aj=0; aj<jcell->nr; aj++) {
	      j=jcell->atoms[aj];
	      if ( (bSelected && (j==i)) || (!bSelected) ) {
		/* check if this once was a h-bond */
		idx=NOTSET;
		for (k=0; (k<donors[grp][i].nrhb) && (idx==NOTSET); k++)
		  if (j==donors[grp][i].hb[k].a)
		    idx=k;
		if ( is_hbond(a[ grp+grD][i],
			      a[ grp+grH][i],
			      a[ogrp+grA][j],
			      rcut,ccut,x,bBox,hbox,&dist,&ang) ) {
		  /* add to index if not already there */
		  if (idx==NOTSET) {
		    if (donors[grp][i].nrhb>=donors[grp][i].maxnr) {
		      donors[grp][i].maxnr+=10;
		      srenew(donors[grp][i].hb,donors[grp][i].maxnr);
		    }
		    donors[grp][i].hb[donors[grp][i].nrhb].a=j;
		    donors[grp][i].hb[donors[grp][i].nrhb].nr=nrhb;
		    idx=donors[grp][i].nrhb;
		    donors[grp][i].nrhb++;
		    if (nrhb>=max_nrhb) {
		      max_nrhb+=HBINC;
		      srenew(hbexist,max_nrhb);
		      for (l=max_nrhb-HBINC; l<max_nrhb; l++)
			snew(hbexist[l],max_nframes);
		    }
		    nrhb++;
		  }
		  /* update matrix */
		  hbexist[donors[grp][i].hb[idx].nr][nframes] |= hbYES;
		  
		  /* count number of hbonds per frame */
		  nhb[nframes]++;
		  
		  /* make angle and distance distributions */
		  ang*=RAD2DEG;
		  adist[(int)( ang/abin)]++;
		  rdist[(int)(dist/rbin)]++;
		  
		  if (!bTwo) {
		    resdist=abs(top.atoms.atom[a[ grp+grD][i]].resnr-
				top.atoms.atom[a[ogrp+grA][j]].resnr);
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
		  int  ins_d_k,ins_a_k;
		  
		  ins_d=ins_a=FALSE;
		  ins_d_dist=ins_d_ang=ins_a_dist=ins_a_ang=1e6;
		  
		  /* loop over gridcells adjacent to i (xk,yk,zk) */
		  LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xi,yi,zi,ngrid) {
		    kcell=&grid[xk][yk][zk][grIA];
		    /* loop over acceptor atoms from ins group 
		       in this adjacent gridcell (kcell) */
		    for (ak=0; ak<kcell->nr; ak++) {
		      k=kcell->atoms[ak];
		      if (is_hbond(a[grp+grD][i],
				   a[grp+grH][i],
				   a[   grIA][k],
				   rcut,ccut,x,bBox,hbox,&dist,&ang))
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
		  LOOPGRIDINNER(xk,yk,zk,xkk,ykk,zkk,xj,yj,zj,ngrid) {
		    kcell=&grid[xk][yk][zk][grID];
		    /* loop over donor atoms from ins group 
		       in this adjacent gridcell (kcell) */
		    for (ak=0; ak<kcell->nr; ak++) {
		      k=kcell->atoms[ak];
		      if (is_hbond(a[    grID][k],
				   a[    grIH][k],
				   a[ogrp+grA][j],
				   rcut,ccut,x,bBox,hbox,&dist,&ang))
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
		      if (nrhb>=max_nrhb) {
			max_nrhb+=HBINC;
			srenew(hbexist,max_nrhb);
			for (l=max_nrhb-HBINC; l<max_nrhb; l++)
			  snew(hbexist[l],max_nframes);
		      }
		      nrhb++;
		    }
		    
		    /* mark insertion in hbond index */
		    hbexist[donors[grp][i].hb[idx].nr][nframes] |= hbINS;
		    
		    /* print insertion info to file */
		    fprintf(fpins,
			    "%4g: %4d:%3.3s%4d%3.3s -> "
			    "%4d:%3.3s%4d%3.3s (%4.2f,%2.0f) - "
			    "%4d:%3.3s%4d%3.3s (%4.2f,%2.0f)\n",t,
			    a[grIA][ins_d_k]+1,
			    *top.atoms.resname[top.atoms.atom[a[grIA][ins_d_k]].resnr],
			    top.atoms.atom[a[grIA][ins_d_k]].resnr+1,
			    *top.atoms.atomname[a[grIA][ins_d_k]],
			    a[grp+grD][i]+1,
			    *top.atoms.resname[top.atoms.atom[a[grp+grD][i]].resnr],
			    top.atoms.atom[a[grp+grD][i]].resnr+1,
			    *top.atoms.atomname[a[grp+grD][i]],
			    ins_d_dist,ins_d_ang*RAD2DEG,
			    a[ogrp+grA][j]+1,
			    *top.atoms.resname[top.atoms.atom[a[ogrp+grA][j]].resnr],
			    top.atoms.atom[a[ogrp+grA][j]].resnr+1,
			    *top.atoms.atomname[a[ogrp+grA][j]],
			    ins_a_dist,ins_a_ang*RAD2DEG);
		  }
		}
	      }
	    } /* for aj  */
	  }
	  ENDLOOPGRIDINNER;
	} /* for ai  */
      } /* for grp */
    } /* LOOPGRIDOUTER */
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
    sort_hb(nrhb,nr_a,donors);
    
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
		    "Hydrogen Bond Distribution","Distance (nm)","");
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
		    "Hydrogen Bond Distribution","Angle (\\SO\\N)","");
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
    
    
    if (opt2bSet("-ac",NFILE,fnm)) {
      /* build hbexist matrix in reals for autocorr */
      snew(rhbex,nrhb);
      for(i=0; i<nrhb; i++) {
	snew(rhbex[i],nframes);
	for(j=0; j<nframes; j++)
	  rhbex[i][j]=(hbexist[i][j] & hbYES);
      }
      low_do_autocorr(opt2fn("-ac",NFILE,fnm), "Hydrogen Bond Autocorrelation",
		      nframes,nrhb,-1,rhbex,time[1]-time[0],eacNormal,1,
		      TRUE,TRUE,TRUE,NULL,NULL,FALSE,0,0,0);
      for(i=0; i<nrhb; i++)
	sfree(rhbex[i]);
      sfree(rhbex);
    }
    
    if (opt2bSet("-hbm",NFILE,fnm) || bInsert) {
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
	mat.nmap=hbNR;
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
      sfree(mat.matrix);
      sfree(mat.map);
    }
  }
  
  if (opt2bSet("-hbn",NFILE,fnm)) {
    fp = opt2FILE("-hbn",NFILE,fnm,"w");
    for (grp=gr0; grp<=(bTwo?gr1:gr0); grp+=grINC) {
      fprintf(fp,"[ %s ]",grpnames[grp/grINC]);
      for (i=0; i<isize[grp/grINC]; i++) {
	fprintf(fp,(i%15)?" ":"\n");
	fprintf(fp,"%4d",index[grp/grINC][i]+1);
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
	for(j=0; j<donors[grp][i].nrhb; j++)
	  fprintf(fp,"%6d %6d %6d\n",
		  a[ grp+grD][i]+1,
		  a[ grp+grH][i]+1,
		  a[ogrp+grA][donors[grp][i].hb[j].a]+1);
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
	  fprintf(fp,"%4u",i+1);
	  j++;
	}
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  thanx(stdout);
  
  return 0;
}
