/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_msd_c = "$Id$";

 
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
#include "pbc.h"

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

t_corr *init_corr(int nrgrp,int type,real dim_factor,
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
  this->dim_factor = dim_factor;
  
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
  if (bMol)
    this->mols = &(top->blocks[ebMOLS]);
  
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

static void init_restart(t_corr *this,int nrestart,real dt)
{
  int i;
  
  this->nrestart = max(1,nrestart);
  this->delta_t  = (nrestart > 1) ? dt : 0;
  
  printf("\nThe number of starting points you requested takes %lu"
	 " bytes memory\n\n",
	 (unsigned long) this->natoms*this->nrestart*sizeof(rvec));

  snew(this->x0,this->nrestart);
  for(i=0; (i<this->nrestart); i++)
    snew(this->x0[i],this->natoms);
  snew(this->n_offs,this->nrestart);
}

static void corr_print(t_corr *this,char *fn,char *title,char *yaxis,
		       real beginfit,real endfit,
		       real *DD,real *SigmaD,char *grpname[])
{
  FILE *out;
  int  i,j;
  
  out=xvgropen(fn,title,"Time (ps)",yaxis);
  if (DD) {
    fprintf(out,"# Diffusion constants fitted from time %g to %g (ps)\n",
	    beginfit,endfit);
    for(i=0; i<this->ngrp; i++) 
      fprintf(out,"# D[%10s] = %.3f (+/- %.3f) (1e-5 cm^2/s)\n",
	      grpname[i],DD[i],SigmaD[i]);
  }
  for(i=0; i<this->nframes; i++) {
    fprintf(out,"%10g",this->time[i]);
    for(j=0; j<this->ngrp; j++)
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
    if ((thistime(this) >= (this->nlast*this->delta_t)) && (nr==0)) {
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
      fatal_error(0,"Error: did not expect option value %d",this->type);
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
	       int nrestart,real dt)
{
  rvec         *x[2];
  real         t;
  int          i,j,status,cur=0,maxframes=0;
#define        prev (1-cur)
  matrix       box;
  
  this->natoms=read_first_x(&status,fn,&this->t0,&(x[prev]),box);
#ifdef DEBUG
  fprintf(stderr,"Read %d atoms for first frame\n",this->natoms);
#endif

  snew(x[cur],this->natoms);

  init_restart(this,nrestart,dt);
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

    this->time[this->nframes] = t - this->t0;
    
    /* loop over all groups in index file */
    for(i=0; (i<this->ngrp); i++) {
      /* nice for putting things in boxes and such */
      prep1(this,gnx[i],index[i],x[cur],x[prev],box);
      /* calculate something useful, like mean square displacements */
      calc_corr(this,i,gnx[i],index[i],x[cur],calc1);
    }
    cur=prev;
    
    this->nframes++;
  } while (read_next_x(status,&t,this->natoms,x[cur],box));
  fprintf(stderr,"\n");
  
  close_trj(status);
}

static void prep_data_mol(t_corr *this,int gnx,atom_id index[],
			  rvec xcur[],rvec xprev[],matrix box)
{
  int  i,j,k,m,d,ind;
  rvec hbox;

  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m]=0.5*box[m][m];
  for(i=0; i<gnx; i++) {
    ind=index[i];
    for(j=this->mols->index[ind]; j<this->mols->index[ind+1]; j++) {
      k=this->mols->a[j];
      for(m=DIM-1; m>=0; m--) {
	while (xcur[k][m]-xprev[k][m] <= -hbox[m])
	  for(d=0; d<=m; d++)
	    xcur[k][d] += box[m][d];
	while (xcur[k][m]-xprev[k][m] >  hbox[m])
	  for(d=0; d<=m; d++)
	    xcur[k][d] -= box[m][d];
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
    fatal_error(0,"Options got screwed. Did not expect value %d\n",this->type);
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

void printmol(t_corr *this,char *fn)
{
#define NDIST 100
  FILE  *out;
  t_lsq lsq1;
  int   i,j;
  real  a,b;
  
  out=xvgropen(fn,"Diffusion Coefficients / Molecule","Molecule","D");
  
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
    fprintf(out,"%10d  %10g\n",i,a*FACTOR/this->dim_factor);
  }
  fclose(out);
  do_view(fn,"-graphtype bar");
}

void do_corr(int NFILE, t_filenm fnm[],int nrgrp,
	     t_topology *top,bool bMol,bool bMW,
	     int type,real dim_factor,int axis,
	     int nrestart,real dt,real beginfit,real endfit)
{
  t_corr       *msd;
  int          *gnx;
  atom_id      **index;
  char         **grpname;
  int          i,i0,i1,j,N;
  real         *DD,*SigmaD,a,a2,b;
  
  gnx     = (int *)calloc(nrgrp,sizeof(int));
  index   = (atom_id **)calloc(nrgrp,sizeof(atom_id *));
  grpname = (char **)calloc(nrgrp,sizeof(char *));
    
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),nrgrp,gnx,index,grpname);

  msd = init_corr(nrgrp,type,dim_factor,bMW,bMol,top);
  
  corr_loop(msd,ftp2fn(efTRX,NFILE,fnm),gnx,index,
	    bMol ? calc1_mol : (bMW ? calc1_mw : calc1_norm),
	    bMol ? prep_data_mol : prep_data_norm,nrestart,dt);
  
  /* Correct for the number of points */
  for(j=0; (j<msd->ngrp); j++)
    for(i=0; (i<msd->nframes); i++)
      msd->data[j][i] /= msd->ndata[j][i];

  if (bMol) 
    printmol(msd,opt2fn("-mol",NFILE,fnm));

  DD     = NULL;
  SigmaD = NULL;
  for(i0=0; i0<msd->nframes && msd->time[i0]<beginfit; i0++) 
    ;
  if (endfit == -1) {
      i1 = msd->nframes;
      endfit = msd->time[i1-1];
  } else
    for(i1=i0; i1<msd->nframes && msd->time[i1]<=endfit; i1++)
		      ;
  N = i1-i0;
  if (N <= 2) {
    fprintf(stdout,"Not enough points for fitting (%d).\n"
	    "Can not determine the diffusion constant.\n",N);
  } else {
    snew(DD,msd->ngrp);
    snew(SigmaD,msd->ngrp);
    for(j=0; j<msd->ngrp; j++) {
      if (N >= 4) {
	lsq_y_ax_b(N/2,&(msd->time[i0]),&(msd->data[j][i0]),&a,&b);
	lsq_y_ax_b(N/2,&(msd->time[i0+N/2]),&(msd->data[j][i0+N/2]),&a2,&b);
	SigmaD[j] = fabs(a-a2);
      } else
	SigmaD[j] = 0;
      lsq_y_ax_b(N,&(msd->time[i0]),&(msd->data[j][i0]),&(DD[j]),&b);
      DD[j]     *= FACTOR/msd->dim_factor;
      SigmaD[j] *= FACTOR/msd->dim_factor;
      fprintf(stdout,"D[%10s] %.3f (+/- %.3f) 1e-5 cm^2/s\n",
		grpname[j],DD[j],SigmaD[j]);
    }
  }
  /* Print and show mean square displacement */
  corr_print(msd,opt2fn("-o",NFILE,fnm),
	     "Mean Square Displacement",
	     "MSD (nm\\S2\\N)",
	     beginfit,endfit,DD,SigmaD,grpname);
  do_view(opt2fn("-o",NFILE,fnm),NULL);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_msd computes the mean square displacement (MSD) of atoms from",
    "their initial positions. This provides an easy way to compute",
    "the diffusion constant using the Einstein relation.",
    "The diffusion constant is calculated by least squares fitting a",
    "straight line through the MSD from [TT]-beginfit[tt] to",
    "[TT]-endfit[tt]. An error estimate given, which is the difference",
    "of the diffusion coefficients obtained from fits over the two halfs",
    "of the fit interval.[PAR]",
    "Option [TT]-mol[tt] plots the MSD for molecules, this implies",
    "[TT]-mw[tt].[PAR]",
    "Mean Square Displacement calculations and Correlation functions",
    "can be calculated more accurately, when using multiple starting",
    "points (see also Gromacs Manual). You can select the number of",
    "starting points, and the interval (in picoseconds) between starting",
    "points. More starting points implies more CPU time."
  };
  static char *normtype[]= { NULL,"no","x","y","z",NULL };
  static char *axtitle[] = { NULL,"no","x","y","z",NULL };
  static int  ngroup     = 1;
  static int  nrestart   = 1;
  static real dt         = 0; 
  static real beginfit   = 0; 
  static real endfit     = -1; 
  static bool bMW        = TRUE;
  t_pargs pa[] = {
    { "-type",    FALSE, etENUM, {normtype},
      "Compute diffusion coefficient in one direction" },
    { "-lateral", FALSE, etENUM, {axtitle}, 
      "Calculate the lateral diffusion in a plane perpendicular to" },
    { "-ngroup",  FALSE, etINT,  {&ngroup},
      "Number of groups to calculate MSD for" },
    { "-mw",      FALSE, etBOOL, {&bMW},
      "Mass weighted MSD" },
    { "-nrestart",FALSE, etINT,  {&nrestart},
      "Number of restarting points in trajectory" },
    { "-trestart",FALSE, etREAL, {&dt},
      "Time between restarting points in trajectory" },
    { "-beginfit",FALSE, etREAL, {&beginfit},
      "Start time for fitting the MSD" },
    { "-endfit",FALSE, etREAL, {&endfit},
      "End time for fitting the MSD (-1 is to the end)" }
  };

  t_filenm fnm[] = { 
    { efTRX, NULL, NULL,  ffREAD },
    { efTPS, NULL, NULL,  ffREAD }, 
    { efNDX, NULL, NULL,  ffOPTRD },
    { efXVG, NULL, "msd", ffWRITE },
    { efXVG, "-mol", "diff_mol",ffOPTWR }
  };
#define NFILE asize(fnm)

  t_topology  top;
  matrix      box;
  char        title[256];
  rvec        *xdum;
  bool        bTop,bMol;
  int         axis,type;
  real        dim_factor;

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  if (ngroup < 1)
    fatal_error(0,"Must have at least 1 group (now %d)",ngroup);

  bMol = opt2bSet("-mol",NFILE,fnm);
  if (bMol) {
    bMW  = TRUE;
    fprintf(stderr,"Calculating diffusion coefficients for molecules.\n");
  }

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
  else
    axis = 0;
  fprintf(stdout,"type = %d, axis = %d, dim_factor = %g\n",
	  type,axis,dim_factor);
  bTop=read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xdum,NULL,box,bMW); 
  if (bMol && !bTop)
    fatal_error(0,"Could not read a topology from %s. Try a tpr file instead.",
		ftp2fn(efTPS,NFILE,fnm));
    
  do_corr(NFILE,fnm,ngroup,
	  &top,bMol,bMW,type,dim_factor,axis,
	  nrestart,dt,beginfit,endfit);

  thanx(stderr);
  
  return 0;
}
