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
#include "gstat.h"
#include "tpxio.h"

#define FACTOR  1000.0	/* Convert nm^2/ps to 10e-5 cm^2/s */

extern float dim_factor;

class c_corr {
public:
  int  natoms,nrestart;
  real t0,delta_t;
  int  nframes,maxframes,nlast,ngrp;
  int  *n_offs;
  int  **ndata;
  real **data,*time,*data_x,*data_y,*data_z,*data_xy;
  rvec **x0;
public:
  /* Constructor, destructor */
  c_corr(int nrgrp);
  virtual ~c_corr();
  
  /* Normal methods */
  void normalise();
  void init_restart();
  void print(char *fn,char *title,char *yaxis,bool bXvgr);
  void subtitle(FILE *out);
  void calc(int nr,int nx,atom_id index[],rvec xc[]);
  real thistime() { return time[nframes]; };
  void loop(char *fn,int gnx[],atom_id *index[],
	    t_first_x *fx,t_next_x *nx);
  int  in_data(int nx00) { return nframes-n_offs[nx00]; };
  
  /* Virtual methods, may/must be overridden by children */
  virtual real calc1(int nx,atom_id index[],int nx0,rvec xc[]) = 0;
  virtual void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
			 matrix box) = 0;
};

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

float dim_factor = 6.0;         // 6 for 3D, 4 for 2D, 2 for 1D

c_corr::c_corr(int nrgrp)
{
  int i;

  maxframes = 0;
  ngrp      = nrgrp;
  nframes   = 0;
  nlast     = 0;

  ndata=(int **)calloc(nrgrp,sizeof(ndata[0]));
  data=(real **)calloc(nrgrp,sizeof(data[0]));
  for(i=0; (i<nrgrp); i++) {
    ndata[i]=NULL;
    data[i]=NULL;
  }
  time=NULL;
}

c_corr::~c_corr()
{
  int i;
  
  free(n_offs);
  for(i=0; (i<nrestart); i++)
    free(x0[i]);
  free(x0);
}

void c_corr::init_restart(void)
{
  char *intro[] = {
    "\n",
    "Mean Square Displacement calculations and Correlation functions",
    "can be calculated more accurately, when using multiple starting",
    "points (see also Gromacs Manual). You can select the number of",
    "starting points, and the interval (in picoseconds) between starting",
    "points. More starting points implies more CPU time."
    };
  int    i;
  double dt;
  
  for(i=0; (i<asize(intro)); i++)
    printf("%s\n",intro[i]);
  do {
    printf("\nNumber of starting points ( >= 1) ? "); fflush(stdout);
    scanf("%d",&nrestart);
  } while (nrestart <= 0);
  if (nrestart > 1)
    do {
      printf("\nTime interval (> 0) ? "); fflush(stdout);
      scanf("%lf",&dt);
      delta_t=dt;
    } while (delta_t <= 0.0);
  else
    delta_t = 0.0;
  
  printf("\nThe number of starting points you requested takes %d"
	 " bytes memory\n\n",natoms*nrestart*sizeof(rvec));

  x0=(rvec **)calloc(nrestart,sizeof(rvec *));
  for(i=0; (i<nrestart); i++)
    x0[i]=(rvec *)calloc(natoms,sizeof(rvec));
  n_offs=(int *)calloc(nrestart,sizeof(int));
}

void lsq_y_ax_b2(int n, real x[], real y[], real *a, real *b,
		 real *sa,real *sb)
{
  int    i;
  double yx,xx,sx,sy,sa2;

  yx=xx=sx=sy=0.0;
  for (i=0; (i < n); i++) {
    yx+=y[i]*x[i];
    xx+=x[i]*x[i];
    sx+=x[i];
    sy+=y[i];
  }
  *a=(n*yx-sy*sx)/(n*xx-sx*sx);
  *b=(sy-(*a)*sx)/n;
  sa2=0.0;
  for(i=0; (i<n); i++) {
    real tmp=(y[i]-(*a)*x[i]-(*b));
    sa2+=tmp*tmp;
  }
  *sa=sqrt((sa2/(n*(n-2)))*(1.0/(xx/n-(sx/n)*(sx/n))));
  *sb=(*sa)*sqrt(xx/n);
}


void c_corr::subtitle(FILE *out)
{
  int i;
  real aa,bb,da,db,D;
  
  for(i=0; (i<ngrp); i++) {
    lsq_y_ax_b2(nframes,time,data[i],&aa,&bb,&da,&db);
    D=aa*FACTOR/dim_factor;
    fprintf(out,"# D[%d] = %.4f (1e-5 cm^2 s^-1)\n",i,D);
  }
}

void c_corr::print(char *fn,char *title,char *yaxis,bool bXvgr)
{
  FILE *out;
  int  i,j;
  real aa,bb,da,db;
  
  /*  real err_fac; */
  
  /*  err_fac=sqrt(2.0/3.0*sqrt(M_PI)); */
  for(j=0; (j<ngrp); j++)
    for(i=0; (i<nframes); i++)
      data[j][i]/=ndata[j][i];
  out=xvgropen(fn,title,"Time (ps)",yaxis);
  if (bXvgr)
    subtitle(out);
    
  lsq_y_ax_b2(nframes,time,data[0],&aa,&bb,&da,&db);
  for(i=0; (i<nframes); i++) {
    fprintf(out,"%10g",time[i]);
    for(j=0; (j<ngrp); j++)
      fprintf(out,"  %10g",data[j][i]);
    fprintf(out,"\n");
  }
  fclose(out);
}

// called from c_corr::loop, to do the main calculations
void c_corr::calc(int nr,int nx,atom_id index[],rvec xc[])
{
  int  nx0;
  real g;
  
  /* Check for new starting point */
  if (nlast < nrestart) {
    if ((thistime() >= (t0+nlast*delta_t)) && (nr==0)) {
      printf("New starting point\n");
      memcpy(x0[nlast],xc,natoms*sizeof(xc[0]));
      n_offs[nlast]=nframes;
      nlast++;
    }
  }

  // nx0 appears to be the number of new starting points,
  // so for all starting points, call calc1. 
  for(nx0=0; (nx0<nlast); nx0++) {
    g=calc1(nx,index,nx0,xc);
#ifdef DEBUG2
    printf("g[%d]=%g\n",nx0,g);
#endif
    data [nr][in_data(nx0)]+=g;
    ndata[nr][in_data(nx0)]++;
  }
}

// this is the mean loop for the correlation type functions 
// fx and nx are file pointers to things like read_first_x and
// read_next_x
void c_corr::loop(char *fn,int gnx[],atom_id *index[],
		  t_first_x *fx,t_next_x *nx)
{
  rvec         *x[2];
  real         t;
  int          i,j,status,cur=0,maxframes=0;
#define        prev (1-cur)
  matrix       box;
  
  natoms=fx(&status,fn,&t0,&(x[prev]),box);
#ifdef DEBUG
  fprintf(stderr,"Read %d atoms for first frame\n",natoms);
#endif


  x[cur]=(rvec *)calloc(natoms,sizeof(x[cur][0]));

  init_restart();
  memcpy(x[cur],x[prev],natoms*sizeof(x[prev][0]));
  t=t0;
  do {
    if (nframes >= maxframes-1) {
      if (maxframes==0) {
	for(i=0; (i<ngrp); i++) {
	  ndata[i]=(int *)calloc(0,sizeof(**ndata));
	  data[i]=(real *)calloc(0,sizeof(**data));
	}
	time=(real *)calloc(0,sizeof(*time)); 
      }
      maxframes+=10;
      for(i=0; (i<ngrp); i++) {
	ndata[i]=(int *)realloc(ndata[i],maxframes*sizeof(**ndata));
	data[i]=(real *)realloc(data[i],maxframes*sizeof(**data));
	for(j=maxframes-10; j<maxframes; j++) {
	  ndata[i][j]=0;
	  data[i][j]=0;
	}
      }
      time=(real *)realloc(time,maxframes*sizeof(*time)); 
    }

    time[nframes]=t;
    
    // loop over all groups in index file
    for(i=0; (i<ngrp); i++) {
      // nice for putting things in boxes and such
      prep_data(gnx[i],index[i],x[cur],x[prev],box);
      // calculate something usefull, like mean square displacements
      calc(i,gnx[i],index[i],x[cur]);
    }
    cur=prev;
    
    nframes++;
  } while (nx(status,&t,natoms,x[cur],box));
  fprintf(stderr,"\n");
  
  close_trj(status);
}

class c_msd: public c_corr {
public:
  c_msd(int nrgrp);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
  void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		 matrix box);
};

class c_msd_m: public c_msd {
public:
  real *mass;
  c_msd_m(t_topology *top,int nrgrp);
  real calc_one(int ix,int nx0,rvec xc[],real *tm);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
};

class c_msd_mol: public c_msd_m {
private:
  t_block *mols;
  int     nnx;
  t_lsq   **lsq;
public:
  c_msd_mol(t_topology *top,int nrgrp);
  real calc1(int nx,atom_id index[],int nx0,rvec xc[]);
  void printdist(char *fn,char *diff);
  void prep_data(int gnx,atom_id index[],rvec xcur[],rvec xprev[],
		 matrix box);
};

c_msd::c_msd(int nrgrp): c_corr(nrgrp)
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

c_msd_m::c_msd_m(t_topology *top,int nrgrp): c_msd(nrgrp)
{
  t_atoms *atoms;
  int i;
  
  atoms=&top->atoms;
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

c_msd_mol::c_msd_mol(t_topology *top,int nrgrp): c_msd_m(top,nrgrp)
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

void do_corr(int NFILE, t_filenm fnm[],c_msd *msd,int nrgrp,
	     t_topology *top)
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
    
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),nrgrp,gnx,index,grpname);
  sprintf(title,"Mean Square Displacement");

  msd->loop(ftp2fn(efTRX,NFILE,fnm),gnx,index,fx,nx);
  
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
  static char *normtype[]= { "no","x","y","z",NULL };
  static char *axtitle[] = { "no","x","y","z",NULL };
  static int  ngroup     = 1;
  static int  bMW        = 1;
  t_pargs pa[] = {
    { "-type",    FALSE, etENUM, &normtype,
      "Compute diffusion coefficient in one direction" },
    { "-lateral", FALSE, etENUM, &axtitle, 
      "Calculate the lateral diffusion in a plane perpendicular to" },
    { "-ngroup",  FALSE, etINT, &ngroup,
      "Number of groups to calculate MSD for" },
    { "-mw",      FALSE, etBOOL,&bMW,
      "Mass weighted MSD" }
  };
  static char *bugs[] = {
    "The diffusion constant given in the title of the graph for lateral"
    " diffusion has to be multiplied by 6/4"
  };

  t_filenm fnm[] = { 
    { efTRX, NULL, NULL,  ffREAD },
    { efTPS, NULL, NULL,  ffREAD }, 
    { efNDX, NULL, NULL,  ffOPTRD },
    { efXVG, NULL, "msd",  ffWRITE },
    { efXVG, "-m", "mol", ffOPTWR },
    { efXVG, "-d", "diff",ffOPTWR }
  };
#define NFILE asize(fnm)

  t_topology  top;
  matrix      box;
  char        cc,title[256];
  rvec        *xdum;
  int         bTop;

  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (ngroup < 1)
    fatal_error(0,"Must have at least 1 group (now %d)",ngroup);
    
  if (normtype[0][0]!='n') {
    type=normtype[0][0] - 'x'+1; /* See defines above */
    dim_factor=2.0;
  }
  if ((type==NORMAL) && (axtitle[0][0]!='n')) {
    type=LATERAL;
    dim_factor = 4.0;
    axis = (axtitle[0][0] - 'x'); /* See defines above */
  }

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xdum,NULL,box,bMW); 
  
  if (opt2bSet("-m",NFILE,fnm)) {
    if (!bTop)
      fatal_error(0,"need a topology for molecules");
    c_msd_mol cm(&top,ngroup);
    
    do_corr(NFILE,fnm,&cm,ngroup,&top);
    cm.printdist(opt2fn("-m",NFILE,fnm),opt2fn("-d",NFILE,fnm));
  }
  else {
    /* MSD class */
    if (bMW)
      do_corr(NFILE,fnm,new c_msd_m(&top,ngroup),ngroup,&top);
    else
      do_corr(NFILE,fnm,new c_msd(ngroup),ngroup,&top);
  }
  thanx(stdout);
  
  return 0;
}
