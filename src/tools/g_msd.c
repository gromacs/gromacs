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
/* NORMAL = total diffusion coefficient (default). X,Y,Z is diffusion 
   coefficient in X,Y,Z direction. LATERAL is diffusion coefficient in
   plane perpendicular to axis
*/
enum { NOT_USED, NORMAL, X, Y, Z, LATERAL };

typedef struct {
  real    t0,delta_t,dim_factor;
  real    **data,*time,*data_x,*data_y,*data_z,*data_xy,*mass;
  rvec    **x0;
  t_block *mols;
  t_lsq   **lsq;
  int     type,axis,natoms,nrestart,nnx,nframes,nlast,ngrp;
  int     *n_offs;
  int     **ndata;
} t_corr;

typedef real t_calc_func(t_corr *,int,atom_id[],int,rvec[]);
typedef void t_prep_data_func(t_corr *this,int gnx,atom_id index[],
			      rvec xcur[],rvec xprev[],matrix box);
			      
static real thistime(t_corr *this) 
{
  return this->time[this->nframes]; 
}

static bool in_data(t_corr *this,int nx00) 
{ 
  return this->nframes-this->n_offs[nx00]; 
}

t_corr *init_corr(int nrgrp,int type,int axis,real dim_factor,
		  bool bMass,bool bMol,t_topology *top)
{
  t_corr  *this;
  t_atoms *atoms;
  int     i;

  snew(this,1);
  this->type      = type;
  this->ngrp      = nrgrp;
  this->nframes   = 0;
  this->nlast     = 0;

  snew(this->ndata,nrgrp);
  snew(this->data,nrgrp);
  for(i=0; (i<nrgrp); i++) {
    this->ndata[i] = NULL;
    this->data[i]  = NULL;
  }
  this->time = NULL;
  this->lsq  = NULL;
  if (bMass) {
    atoms=&top->atoms;
    snew(this->mass,atoms->nr);
    for(i=0; (i<atoms->nr); i++) {
      this->mass[i]=atoms->atom[i].m;
    }
  }
  else if (bMol) {
    this->mols = &(top->blocks[ebMOLS]);
  }
  
  return this;
}

static void done_corr(t_corr *this)
{
  int i;
  
  sfree(this->n_offs);
  for(i=0; (i<this->nrestart); i++)
    sfree(this->x0[i]);
  sfree(this->x0);
}

static void init_restart(t_corr *this)
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
    scanf("%d",&this->nrestart);
  } while (this->nrestart <= 0);
  if (this->nrestart > 1)
    do {
      printf("\nTime interval (> 0) ? "); fflush(stdout);
      scanf("%lf",&dt);
      this->delta_t=dt;
    } while (this->delta_t <= 0.0);
  else
    this->delta_t = 0.0;
  
  printf("\nThe number of starting points you requested takes %d"
	 " bytes memory\n\n",this->natoms*this->nrestart*sizeof(rvec));

  snew(this->x0,this->nrestart);
  for(i=0; (i<this->nrestart); i++)
    snew(this->x0[i],this->natoms);
  snew(this->n_offs,this->nrestart);
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

static void subtitle(FILE *out,t_corr *this)
{
  int i;
  real aa,bb,da,db,D;
  
  for(i=0; (i<this->ngrp); i++) {
    lsq_y_ax_b2(this->nframes,this->time,this->data[i],&aa,&bb,&da,&db);
    D=aa*FACTOR/this->dim_factor;
    fprintf(out,"# D[%d] = %.4f (1e-5 cm^2 s^-1)\n",i,D);
  }
}

static void corr_print(t_corr *this,char *fn,char *title,char *yaxis,bool bXvgr)
{
  FILE *out;
  int  i,j;
  real aa,bb,da,db;
  
  /*  real err_fac; */
  
  /*  err_fac=sqrt(2.0/3.0*sqrt(M_PI)); */
  for(j=0; (j<this->ngrp); j++)
    for(i=0; (i<this->nframes); i++)
      this->data[j][i]/=this->ndata[j][i];
  out=xvgropen(fn,title,"Time (ps)",yaxis);
  if (bXvgr)
    subtitle(out,this);
    
  lsq_y_ax_b2(this->nframes,this->time,this->data[0],&aa,&bb,&da,&db);
  for(i=0; (i<this->nframes); i++) {
    fprintf(out,"%10g",this->time[i]-this->time[0]);
    for(j=0; (j<this->ngrp); j++)
      fprintf(out,"  %10g",this->data[j][i]);
    fprintf(out,"\n");
  }
  fclose(out);
}

/* called from corr_loop, to do the main calculations */
static void calc_corr(t_corr *this,int nr,int nx,atom_id index[],rvec xc[],
		      t_calc_func *calc1)
{
  int  nx0;
  real g;
  
  /* Check for new starting point */
  if (this->nlast < this->nrestart) {
    if ((thistime(this) >= (this->t0+this->nlast*this->delta_t)) && (nr==0)) {
      printf("New starting point\n");
      memcpy(this->x0[this->nlast],xc,this->natoms*sizeof(xc[0]));
      this->n_offs[this->nlast]=this->nframes;
      this->nlast++;
    }
  }

  /* nx0 appears to be the number of new starting points,
   * so for all starting points, call calc1. 
   */
  for(nx0=0; (nx0<this->nlast); nx0++) {
    g = calc1(this,nx,index,nx0,xc);
#ifdef DEBUG2
    printf("g[%d]=%g\n",nx0,g);
#endif
    this->data [nr][in_data(this,nx0)]+=g;
    this->ndata[nr][in_data(this,nx0)]++;
  }
}

static real calc1_norm(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i,ix,m;
  real g,r,r2;
  
  g=0.0;
  
  for(i=0; (i<nx); i++) {
    ix=index[i];
    r2=0.0;
    switch (this->type) {
    case NORMAL:
      for(m=0; (m<DIM); m++) {
	r   = this->x0[nx0][ix][m]-xc[ix][m];
	r2 += r*r;
      }
      break;
    case X:
    case Y:
    case Z:
      r = this->x0[nx0][ix][this->type-1]-xc[ix][this->type-1];
      break;
    case LATERAL:
      for(m=0; (m<DIM); m++) {
	if (m != this->axis) {
	  r   = this->x0[nx0][ix][m]-xc[ix][m];
	  r2 += r*r;
	}
      }
      break;
    default:
      fprintf(stderr,"Error: did not expect option value %d",this->type);
      exit(1);
    }
    g+=r2;
  }
  g/=nx;
  
  return g;
}

/* this is the mean loop for the correlation type functions 
 * fx and nx are file pointers to things like read_first_x and
 * read_next_x
 */
void corr_loop(t_corr *this,char *fn,int gnx[],atom_id *index[],
	       t_calc_func *calc1,t_prep_data_func *prep1,
	       t_first_x *fx,t_next_x *nx)
{
  rvec         *x[2];
  real         t;
  int          i,j,status,cur=0,maxframes=0;
#define        prev (1-cur)
  matrix       box;
  
  this->natoms=fx(&status,fn,&this->t0,&(x[prev]),box);
#ifdef DEBUG
  fprintf(stderr,"Read %d atoms for first frame\n",this->natoms);
#endif

  snew(x[cur],this->natoms);

  init_restart(this);
  memcpy(x[cur],x[prev],this->natoms*sizeof(x[prev][0]));
  t=this->t0;
  do {
    if (this->nframes >= maxframes-1) {
      if (maxframes==0) {
	for(i=0; (i<this->ngrp); i++) {
	  this->ndata[i] = NULL;
	  this->data[i]  = NULL;
	}
	this->time = NULL;
      }
      maxframes+=10;
      for(i=0; (i<this->ngrp); i++) {
	srenew(this->ndata[i],maxframes);
	srenew(this->data[i],maxframes);
	for(j=maxframes-10; j<maxframes; j++) {
	  this->ndata[i][j]=0;
	  this->data[i][j]=0;
	}
      }
      srenew(this->time,maxframes);
    }

    this->time[this->nframes]=t;
    
    /* loop over all groups in index file */
    for(i=0; (i<this->ngrp); i++) {
      /* nice for putting things in boxes and such */
      prep1(this,gnx[i],index[i],x[cur],x[prev],box);
      /* calculate something useful, like mean square displacements */
      calc_corr(this,i,gnx[i],index[i],x[cur],calc1);
    }
    cur=prev;
    
    this->nframes++;
  } while (nx(status,&t,this->natoms,x[cur],box));
  fprintf(stderr,"\n");
  
  close_trj(status);
}

static void prep_data_mol(t_corr *this,int gnx,atom_id index[],
			  rvec xcur[],rvec xprev[],matrix box)
{
  int  i,j,k,m,ind;
  rvec hbox;

  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m]=0.5*box[m][m];
  for(i=0; (i<gnx); i++) {
    ind=index[i];
    for(j=this->mols->index[ind]; (j<this->mols->index[ind+1]); j++) {
      k=this->mols->a[j];
      for(m=0; (m<DIM); m++) {
	while(xcur[k][m]-xprev[k][m] <= hbox[m])
	  xcur[k][m] += box[m][m];
	while(xcur[k][m]-xprev[k][m] >  hbox[m])
	  xcur[k][m] -= box[m][m];
      }      
    }
  }
}

static real calc_one_mw(t_corr *this,int ix,int nx0,rvec xc[],real *tm)
{
  real r2,r,mm;
  int  m;
  
  mm=this->mass[ix];
  if (mm < 1)
    return 0;
  (*tm)+=this->mass[ix];
  r2=0.0;
  switch (this->type) {
  case NORMAL:
    for(m=0; (m<DIM); m++) {
      r   = this->x0[nx0][ix][m]-xc[ix][m];
      r2 += this->mass[ix]*r*r;
    }
    break;
  case X:
  case Y:
  case Z:
    r  = this->x0[nx0][ix][this->type-1]-xc[ix][this->type-1];
    r2 = this->mass[ix]*r*r;
      break;
  case LATERAL:
    for(m=0; (m<DIM); m++) {
      if (m != this->axis) {
	r   = this->x0[nx0][ix][m]-xc[ix][m];
	r2 += this->mass[ix]*r*r;
      }
    }
    break;
  default:
    fprintf(stderr,"Error: options got screwed. ");
    fprintf(stderr,"Did not expect value %d\n",this->type);
    exit(1);
  } /* end switch */
  return r2;
}

static real calc1_mw(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i;
  real g,tm;
  
  g=tm=0.0;
  for(i=0; (i<nx); i++) 
    g+=calc_one_mw(this,index[i],nx0,xc,&tm);
  
  g/=tm;
  
  return g;
}

static void prep_data_norm(t_corr *this,int gnx,atom_id index[],
			   rvec xcur[],rvec xprev[],matrix box)
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

static real calc1_mol(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[])
{
  int  i,ii,j;
  real g,mm,gtot,tt;

  if (this->lsq == NULL) {
    this->nnx = nx;
    snew(this->lsq,this->nrestart);
    for(i=0; (i<this->nrestart); i++)
      snew(this->lsq[i],nx);
  }
  tt=this->time[in_data(this,nx0)];
  gtot=0;
  for(i=0; (i<nx); i++) {
    ii=index[i];
    g=mm=0;
    for(j=this->mols->index[ii]; (j<this->mols->index[ii+1]); j++) 
      g += calc_one_mw(this,this->mols->a[j],nx0,xc,&mm);
    
    g/=mm;
    gtot+=g;
    add_lsq(&(this->lsq[nx0][i]),tt,g);
  }
  return gtot/nx;
}

void printdist(t_corr *this,char *fn,char *difn)
{
#define NDIST 100
  FILE  *out;
  int   *ndist;
  real  *diff,mind=0,maxd=0,dav,a,b;
  t_lsq lsq1;
  int   i,j;
  
  out=xvgropen(difn,"Diffusion Coefficients / Molecule","Molecule","D");
  
  dav=0;
  snew(diff,this->nnx);
  for(i=0; (i<this->nnx); i++) {
    init_lsq(&lsq1);
    for(j=0; (j<this->nrestart); j++) {
      lsq1.sx+=this->lsq[j][i].sx;
      lsq1.sy+=this->lsq[j][i].sy;
      lsq1.xx+=this->lsq[j][i].xx;
      lsq1.yx+=this->lsq[j][i].yx;
      lsq1.np+=this->lsq[j][i].np;
    }
    get_lsq_ab(&lsq1,&a,&b);
    diff[i]=a*FACTOR/this->dim_factor;
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
  dav/=this->nnx;
  fclose(out);
  xvgr_file(difn,"-graphtype bar");
  
  ndist=(int *)calloc(NDIST+1,sizeof(*ndist));
  for(i=0; i<NDIST+1; i++)
    ndist[i]=0;
  for(i=0; (i<this->nnx); i++) {
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

void do_corr(int NFILE, t_filenm fnm[],int nrgrp,
	     t_topology *top,bool bMol,bool bMW,
	     int type,real dim_factor,int axis)
{
  t_corr       *msd;
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

  msd = init_corr(nrgrp,type,axis,dim_factor,bMW,bMol,top);
  
  corr_loop(msd,ftp2fn(efTRX,NFILE,fnm),gnx,index,
	    bMW ? calc1_mw : (bMol ? calc1_mol : calc1_norm), 
	    bMol ? prep_data_mol : prep_data_norm,
	    fx,nx);
  
  corr_print(msd,opt2fn("-o",NFILE,fnm),title,
	     "MSD (nm\\S2\\N)",TRUE);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  if (bMol) 
    printdist(msd,opt2fn("-m",NFILE,fnm),opt2fn("-d",NFILE,fnm));
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_msd computes the mean square displacement of atoms from",
    "their initial positions. This provides an easy way to compute",
    "the diffusion constant using the Einstein relation."
  };
  static char *normtype[]= { NULL,"no","x","y","z",NULL };
  static char *axtitle[] = { NULL,"no","x","y","z",NULL };
  static int  ngroup     = 1;
  static bool bMW        = TRUE;
  t_pargs pa[] = {
    { "-type",    FALSE, etENUM, normtype,
      "Compute diffusion coefficient in one direction" },
    { "-lateral", FALSE, etENUM, axtitle, 
      "Calculate the lateral diffusion in a plane perpendicular to" },
    { "-ngroup",  FALSE, etINT,  &ngroup,
      "Number of groups to calculate MSD for" },
    { "-mw",      FALSE, etBOOL, &bMW,
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
  bool        bMol,bTop;
  int         axis,type;
  real        dim_factor;

  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  if (ngroup < 1)
    fatal_error(0,"Must have at least 1 group (now %d)",ngroup);

  bMol = opt2bSet("-m",NFILE,fnm);
          
  if (normtype[0][0]!='n') {
    type = normtype[0][0] - 'x'+1; /* See defines above */
    dim_factor = 2.0;
  }
  else {
    type       = NORMAL;
    dim_factor = 6.0;
  }
  if ((type==NORMAL) && (axtitle[0][0]!='n')) {
    type=LATERAL;
    dim_factor = 4.0;
    axis = (axtitle[0][0] - 'x'); /* See defines above */
  }

  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xdum,NULL,box,bMW); 
  if (bMol && !bTop)
    fatal_error(0,"Could not read a topology form %s. Try a tpr file instead.",
		ftp2fn(efTPS,NFILE,fnm));
    
  do_corr(NFILE,fnm,ngroup,&top,bMol,bMW,type,dim_factor,axis);

  thanx(stdout);
  
  return 0;
}
