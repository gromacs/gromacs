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
static char *SRCID_g_msd_cc = "$Id$";

 
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "maths.h"
#include "futil.h"
#include "rdgroup.h"
#include "copyrite.h"
#include "typedefs.h"
#include "xvgr.h"
#include "corrutil.h"
#include "gstat.h"

/* NORMAL = total diffusion coefficient (default). X,Y,Z is diffusion 
   coefficient in X,Y,Z direction. LATERAL is diffusion coefficient in
   plane perpendicular to axis
*/
#define NORMAL 0     
#define X 1
#define Y 2
#define Z 3
#define LATERAL 4

static int type = NORMAL;
static int axis = 2;            // default is perpendicular to z-axis

class c_msd: public c_corr {
public:
  c_msd(int nrgrp,int nrframes);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
  void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		 matrix box);
};

class c_msd_m: public c_msd {
public:
  real *mass;
  c_msd_m(t_topology *top,int nrgrp,int nrframes);
  real calc_one(int ix,int nx0,rvec xc[],real *tm);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
};

class c_msd_mol: public c_msd_m {
private:
  t_block *mols;
  int     nnx;
  t_lsq   **lsq;
public:
  c_msd_mol(t_topology *top,int nrgrp,int nrframes);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
  void printdist(char *fn,char *diff);
  void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		 matrix box);
};

c_msd::c_msd(int nrgrp,int nrframes): c_corr(nrgrp,nrframes)
{
}

real c_msd::calc1(int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i,ix,m;
  real g,r,r2;
  
  g=0.0;
  
  for(i=0; (i<nx); i++) {
    ix=index[i];
    r2=0.0;
    switch (type) {
    case NORMAL:
      for(m=0; (m<DIM); m++) {
	r=x0[nx0][ix][m]-xc[ix][m];
	r2+=r*r;
      }
      break;
    case X:
    case Y:
    case Z:
      r=x0[nx0][ix][type-1]-xc[ix][type-1];
      break;
    case LATERAL:
      for(m=0; (m<DIM); m++) {
	if (m != axis) {
	  r=x0[nx0][ix][m]-xc[ix][m];
	  r2+=r*r;
	}
      }
      break;
    default:
      fprintf(stderr,"Error: did not expect option value %d", type);
      exit(1);
    }
    g+=r2;
  }
  g/=nx;
  
  return g;
}

void c_msd_mol::prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
			  matrix box)
{
  int  i,j,k,m,ind;
  rvec hbox;

  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m]=0.5*box[m][m];
  for(i=0; (i<gnx); i++) {
    ind=index[i];
    for(j=mols->index[ind]; (j<mols->index[ind+1]); j++) {
      k=mols->a[j];
      for(m=0; (m<DIM); m++) {
	while(xcur[k][m]-xprev[k][m] <= hbox[m])
	  xcur[k][m] += box[m][m];
	while(xcur[k][m]-xprev[k][m] >  hbox[m])
	  xcur[k][m] -= box[m][m];
      }      
    }
  }
}

real c_msd_m::calc_one(int ix,int nx0,rvec xc[],real *tm)
{
  real r2,r,mm;
  int  m;
  
  mm=mass[ix];
  if (mm < 1)
    return 0;
  (*tm)+=mass[ix];
  r2=0.0;
    switch (type)
    {
    case NORMAL:
      for(m=0; (m<DIM); m++) {
	r=x0[nx0][ix][m]-xc[ix][m];
	r2+=mass[ix]*r*r;
      }
      break;
    case X:
    case Y:
    case Z:
      r=x0[nx0][ix][type-1]-xc[ix][type-1];
      r2 = mass[ix]*r*r;
      break;
    case LATERAL:
      for(m=0; (m<DIM); m++) {
	if (m != axis) {
	  r=x0[nx0][ix][m]-xc[ix][m];
	  r2+=mass[ix]*r*r;
	}
      }
      break;
    default:
      fprintf(stderr,"Error: options got screwd. ");
      fprintf(stderr,"Did not expect value %d\n",type);
      exit(1);
    } /* end switch */
  return r2;
}

real c_msd_m::calc1(int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i;
  real g,tm;
  
  g=tm=0.0;
  for(i=0; (i<nx); i++) 
    g+=calc_one(index[i],nx0,xc,&tm);
  
  g/=tm;
  
  return g;
}

c_msd_m::c_msd_m(t_topology *top,int nrgrp,int nrframes): c_msd(nrgrp,nrframes)
{
  int i;
  t_atoms *atoms;
  
  atoms=&(top->atoms);
  mass=(real *)calloc(atoms->nr,sizeof(mass[0]));
  for(i=0; (i<atoms->nr); i++) {
    mass[i]=atoms->atom[i].m;
  }
}

void c_msd::prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		      matrix box)
{
  int  i,m,ind;
  rvec hbox;

  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m]=0.5*box[m][m];
  for(i=0; (i<gnx); i++) {
    ind=index[i];
    for(m=0; (m<DIM); m++) {
      while(xcur[ind][m]-xprev[ind][m] <= hbox[m])
	xcur[ind][m] += box[m][m];
      while(xcur[ind][m]-xprev[ind][m] >  hbox[m])
	xcur[ind][m] -= box[m][m];
    }      
  }
}

c_msd_mol::c_msd_mol(t_topology *top,int nrgrp,int nrframes): c_msd_m(top,nrgrp,nrframes)
{
  mols=&(top->blocks[ebMOLS]);
  lsq=NULL;
}

real c_msd_mol::calc1(int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i,ii,j;
  real g,mm,gtot,tt;

  if (lsq==NULL) {
    nnx=nx;
    lsq=(t_lsq **)calloc(nrestart,sizeof(lsq[0]));
    for(i=0; (i<nrestart); i++)
      lsq[i]=(t_lsq *)calloc(nx,sizeof(lsq[0][0]));
  }
  tt=time[in_data(nx0)];
  gtot=0;
  for(i=0; (i<nx); i++) {
    ii=index[i];
    g=mm=0;
    for(j=mols->index[ii]; (j<mols->index[ii+1]); j++) 
      g+=calc_one(mols->a[j],nx0,xc,&mm);
    
    g/=mm;
    gtot+=g;
    add_lsq(&(lsq[nx0][i]),tt,g);
  }
  return gtot/nx;
}

void c_msd_mol::printdist(char *fn,char *difn)
{
#define NDIST 100
  FILE  *out;
  int   *ndist;
  real  *diff,mind=0,maxd=0,dav,a,b;
  t_lsq lsq1;
  int   i,j;
  
  out=xvgropen(difn,"Diffusion Coefficients / Molecule","Molecule","D");
  
  dav=0;
  diff=(real *)calloc(nnx,sizeof(real));  
  for(i=0; (i<nnx); i++) {
    init_lsq(&lsq1);
    for(j=0; (j<nrestart); j++) {
      lsq1.sx+=lsq[j][i].sx;
      lsq1.sy+=lsq[j][i].sy;
      lsq1.xx+=lsq[j][i].xx;
      lsq1.yx+=lsq[j][i].yx;
      lsq1.np+=lsq[j][i].np;
    }
    get_lsq_ab(&lsq1,&a,&b);
    diff[i]=a*FACTOR/dim_factor;
    fprintf(out,"%10d  %10g\n",i,diff[i]);
    if (dav == 0) {
      mind=maxd=dav=diff[i];
    }
    else {
      mind=min(mind,diff[i]);
      maxd=max(maxd,diff[i]);
      dav+=diff[i];
    }
  }
  dav/=nnx;
  fclose(out);
  xvgr_file(difn,"-graphtype bar");
  
  ndist=(int *)calloc(NDIST+1,sizeof(*ndist));
  for(i=0; i<NDIST+1; i++)
    ndist[i]=0;
  for(i=0; (i<nnx); i++) {
    int index=(int)(0.5+NDIST*(diff[i]-mind)/(maxd-mind));
    if ((index >= 0) && (index <= NDIST))
      ndist[index]++;
  }
  out=xvgropen(fn,"Distribution of Diffusion Constants","D ()","Number");
  for(i=0; (i<=NDIST); i++)
    fprintf(out,"%10g  %10d\n",
	    mind+(i*(maxd-mind))/NDIST,ndist[i]);
  fclose(out);
  xvgr_file(fn,NULL);
  
}

void do_corr(int NFILE, t_filenm fnm[],c_msd *msd,int nrgrp)
{
  t_first_x    *fx;
  t_next_x     *nx;
  int          *gnx;
  atom_id      **index;
  char         **grpname;
  char         title[256];
  
  fx=read_first_x;
  nx=read_next_x;

  gnx     = (int *)calloc(nrgrp,sizeof(int));
  index   = (atom_id **)calloc(nrgrp,sizeof(atom_id *));
  grpname = (char **)calloc(nrgrp,sizeof(char *));
    
  rd_index(ftp2fn(efNDX,NFILE,fnm),nrgrp,gnx,index,grpname);
  sprintf(title,"Mean  Square Displacement");

  msd->loop(opt2fn("-f",NFILE,fnm),gnx,index,fx,nx);
  
  msd->print(opt2fn("-o",NFILE,fnm),title,
	     "MSD (nm\\S2\\N)",TRUE);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_msd computes the mean square displacement of atoms from",
    "their initial positions. This provides an easy way to compute",
    "the diffusion constant using the Einstein relation.[PAR]",
    "If the -m option is specified, the index is interpreted as containing",
    " molecule numbers and the diffusion coefficient per molecule is ",
    "computed."
  };
  static char *normtype=NULL;
  static char *axtitle =NULL;
  static int  ngroup   = 1;
  static int  nrframes = 10;
  static int  bMW      = 0;
  t_pargs pa[] = {
    { "-type",    FALSE, etSTR, &normtype,
      "If x, y or z is specified, the diffusion coefficient in the x, y, or "
      "z directions is computed." },
    { "-lateral", FALSE, etSTR, &axtitle, 
      "Calculate the lateral diffusion in a plane perpendicular to X, Y, Z" },
    { "-ngroup",  FALSE, etINT, &ngroup,
      "Number of groups to calculate MSD for" },
    { "-nframes", FALSE, etINT, &nrframes,
      "Number of frames in your trajectory" },
    { "-mw",      FALSE, etBOOL,&bMW,
      "Use mass weighted MSD" }
  };
  static char *bugs[] = {
    "The diffusion constant given in the title of the graph for lateral"
    " diffusion has to be multiplied by 6/4"
  };

  t_filenm fnm[] = { 
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffREAD },
    { efTPX, NULL, NULL,  ffREAD },
    { efXVG, NULL, NULL,  ffWRITE },
    { efXVG, "-m", "mol", ffOPTWR },
    { efXVG, "-d", "diff",ffOPTWR }
  };
#define NFILE asize(fnm)

  t_topology *top;
  char        cc;

  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (ngroup < 1)
    fatal_error(0,"Must have at least 1 group (now %d)",ngroup);
    
  if (opt2parg_bSet("-type",asize(pa),pa)) {
    type=(toupper(normtype[0]) - 'X')+1; /* See defines above */
    dim_factor=2.0;
    if ((type < X) || (type > Z)) {
      fprintf(stderr,"Invalid type (%s) given!\n",normtype);
      exit(1);
    }
  }
  else if (opt2parg_bSet("-lateral",asize(pa),pa)) {
    type=LATERAL;
    dim_factor = 4.0;
    axis = (toupper(axtitle[0]) - 'X'); /* See defines above */
    if ((axis < XX) || (type > ZZ)) {
      fprintf(stderr,"Invalid lateral (%s) given!\n",axtitle);
      exit(1);
    }
  }
  
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  
  if (opt2bSet("-m",NFILE,fnm)) {
    c_msd_mol cm(top,ngroup,nrframes);
    
    do_corr(NFILE,fnm,&cm,ngroup);
    cm.printdist(opt2fn("-m",NFILE,fnm),opt2fn("-d",NFILE,fnm));
  }
  else {
    /* MSD class */
    if (bMW) {
      top=read_top(ftp2fn(efTPX,NFILE,fnm));
      do_corr(NFILE,fnm,new c_msd_m(top,ngroup,nrframes),ngroup);
    }
    else {
      do_corr(NFILE,fnm,new c_msd(ngroup,nrframes),ngroup);
    }
  }
  thanx(stdout);
  
  return 0;
}
