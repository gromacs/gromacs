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

#include <string.h>
#include <ctype.h>
#include <math.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "statutil.h"
#include "maths.h"
#include "futil.h"
#include "index.h"
#include "copyrite.h"
#include "typedefs.h"
#include "xvgr.h"
#include "gstat.h"
#include "tpxio.h"
#include "pbc.h"
#include "vec.h"

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
		  bool bMass,bool bMol,real dt,t_topology *top)
{
  t_corr  *this;
  t_atoms *atoms;
  int     i;

  snew(this,1);
  this->type      = type;
  this->axis      = axis;
  this->ngrp      = nrgrp;
  this->nrestart  = 0;
  this->delta_t   = dt;
  this->x0        = NULL;
  this->n_offs    = NULL;
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

static void init_restart(t_corr *this)
{
  int i;
  
  this->nrestart++;
  
  srenew(this->x0,this->nrestart);
  snew(this->x0[this->nrestart-1],this->natoms);
  srenew(this->n_offs,this->nrestart);
}

static void corr_print(t_corr *this,char *fn,char *title,char *yaxis,
		       real beginfit,real endfit,
		       real *DD,real *SigmaD,char *grpname[])
{
  FILE *out;
  int  i,j;
  
  out=xvgropen(fn,title,xvgr_tlabel(),yaxis);
  if (DD) {
    fprintf(out,"# Diffusion constants fitted from time %g to %g (%s)\n",
	    beginfit,endfit,time_unit());
    for(i=0; i<this->ngrp; i++) 
      fprintf(out,"# D[%10s] = %.3f (+/- %.3f) (1e-5 cm^2/s)\n",
	      grpname[i],DD[i],SigmaD[i]);
  }
  for(i=0; i<this->nframes; i++) {
    fprintf(out,"%10g",convert_time(this->time[i]));
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
      r = this->x0[nx0][ix][this->type-X]-xc[ix][this->type-X];
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
      gmx_fatal(FARGS,"Error: did not expect option value %d",this->type);
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
	       t_calc_func *calc1,t_prep_data_func *prep1,real dt)
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

  /* init_restart(this,nrestart,dt); */
  memcpy(x[cur],x[prev],this->natoms*sizeof(x[prev][0]));
  t=this->t0;
  do {
    if (bRmod(t,this->t0,dt))
      init_restart(this);
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
  fprintf(stderr,"\nUsed %d restart points spaced %g %s over %g %s\n\n", 
	  this->nrestart, 
	  convert_time(dt), time_unit(),
	  convert_time(this->time[this->nframes-1]), time_unit() );
  
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
    r  = this->x0[nx0][ix][this->type-X]-xc[ix][this->type-X];
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
    gmx_fatal(FARGS,"Options got screwed. Did not expect value %d\n",this->type);
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
    for(m=DIM-1; m>=0; m--) {
      while(xcur[ind][m]-xprev[ind][m] <= -hbox[m])
	rvec_inc(xcur[ind],box[m]);
      while(xcur[ind][m]-xprev[ind][m] >  hbox[m])
	rvec_dec(xcur[ind],box[m]);
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
  real  a,b,*D,Dav,D2av,VarD;
  
  out=xvgropen(fn,"Diffusion Coefficients / Molecule","Molecule","D");
  
  snew(D,this->nnx);
  Dav = D2av = 0;
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
    D[i]  = a*FACTOR/this->dim_factor;
    Dav  += D[i];
    D2av += sqr(D[i]);
    fprintf(out,"%10d  %10g\n",i,D[i]);
  }
  fclose(out);
  do_view(fn,"-graphtype bar");
  
  /* Compute variance, stddev and error */
  Dav  /= this->nnx;
  D2av /= this->nnx;
  VarD  = D2av - sqr(Dav);
  printf("<D> = %.3f Std. Dev. = %.3f Error = %.3f\n",
	 Dav,sqrt(VarD),sqrt(VarD/this->nnx));
  
  sfree(D);
}

void do_corr(char *trx_file, char *ndx_file, char *msd_file, char *mol_file,
	     int nrgrp, t_topology *top,bool bMW,
	     int type,real dim_factor,int axis,
	     real dt,real beginfit,real endfit)
{
  t_corr       *msd;
  int          *gnx;
  atom_id      **index;
  char         **grpname;
  int          i,i0,i1,j,N;
  real         *DD,*SigmaD,a,a2,b;
  
  snew(gnx,nrgrp);
  snew(index,nrgrp);
  snew(grpname,nrgrp);
    
  if (mol_file && !ndx_file) {
    gnx[0] = top->blocks[ebMOLS].nr;
    snew(index[0],gnx[0]);
    grpname[0] = "Molecules";
    for(i=0; (i<gnx[0]); i++)
      index[0][i] = i;
  }
  else
    get_index(&top->atoms,ndx_file,nrgrp,gnx,index,grpname);

  msd = init_corr(nrgrp,type,axis,dim_factor,bMW,(mol_file!=NULL),dt,top);
  
  corr_loop(msd,trx_file,gnx,index,
	    (mol_file!=NULL) ? calc1_mol : (bMW ? calc1_mw : calc1_norm),
	    (mol_file!=NULL) ? prep_data_mol : prep_data_norm,dt);
  
  /* Correct for the number of points */
  for(j=0; (j<msd->ngrp); j++)
    for(i=0; (i<msd->nframes); i++)
      msd->data[j][i] /= msd->ndata[j][i];

  if (mol_file) 
    printmol(msd,mol_file);

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
  /* Print mean square displacement */
  corr_print(msd,msd_file,
	     "Mean Square Displacement",
	     "MSD (nm\\S2\\N)",
	     beginfit,endfit,DD,SigmaD,grpname);
}

int gmx_msd(int argc,char *argv[])
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
    "[TT]-mw[tt], i.e. for each inidividual molecule an diffusion constant",
    "is computed. When using an index file, it should contain molecule",
    "numbers instead of atom numbers.",
    "Using this option one also gets an accurate error estimate",
    "based on the statistics between individual molecules. Since one usually",
    "is interested in self-diffusion at infinite dilution this is probably",
    "the most useful number.[PAR]",
  };
  static char *normtype[]= { NULL,"no","x","y","z",NULL };
  static char *axtitle[] = { NULL,"no","x","y","z",NULL };
  static int  ngroup     = 1;
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
    { "-trestart",FALSE, etTIME, {&dt},
      "Time between restarting points in trajectory (%t)" },
    { "-beginfit",FALSE, etTIME, {&beginfit},
      "Start time for fitting the MSD (%t)" },
    { "-endfit",FALSE, etTIME, {&endfit},
      "End time for fitting the MSD (%t), -1 is till end" }
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
  char        *trx_file, *tps_file, *ndx_file, *msd_file, *mol_file;
  rvec        *xdum;
  bool        bTop;
  int         axis,type;
  real        dim_factor;

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  trx_file = ftp2fn_null(efTRX,NFILE,fnm);
  tps_file = ftp2fn_null(efTPS,NFILE,fnm);
  ndx_file = ftp2fn_null(efNDX,NFILE,fnm);
  msd_file = ftp2fn_null(efXVG,NFILE,fnm);
  mol_file = opt2fn_null("-mol",NFILE,fnm);
  
  if (ngroup < 1)
    gmx_fatal(FARGS,"Must have at least 1 group (now %d)",ngroup);

  if (mol_file) {
    bMW  = TRUE;
    fprintf(stderr,"Calculating diffusion coefficients for molecules.\n");
  }

  if (normtype[0][0]!='n') {
    type = normtype[0][0] - 'x' + X; /* See defines above */
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
  bTop=read_tps_conf(tps_file,title,&top,&xdum,NULL,box,bMW); 
  if (mol_file && !bTop)
    gmx_fatal(FARGS,"Could not read a topology from %s. Try a tpr file instead.",
		tps_file);
    
  do_corr(trx_file,ndx_file,msd_file,mol_file,ngroup,
	  &top,bMW,type,dim_factor,axis,dt,beginfit,endfit);
  
  view_all(NFILE, fnm);
  
  thanx(stderr);
  
  return 0;
}
