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
  rvec    *com;
  t_lsq   **lsq;
  int     type,axis,ncoords,nrestart,nmol,nframes,nlast,ngrp;
  int     *n_offs;
  int     **ndata;
} t_corr;

typedef real t_calc_func(t_corr *,int,atom_id[],int,rvec[],rvec);
			      
static real thistime(t_corr *this) 
{
  return this->time[this->nframes]; 
}

static bool in_data(t_corr *this,int nx00) 
{ 
  return this->nframes-this->n_offs[nx00]; 
}

t_corr *init_corr(int nrgrp,int type,int axis,real dim_factor,
		  bool bMol,bool bMass,real dt,t_topology *top)
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
  if (bMol) {
    this->nmol = top->blocks[ebMOLS].nr;
    snew(this->mass,this->nmol);
    for(i=0; i<this->nmol; i++)
      this->mass[i] = 1;
  }
  else {
    this->nmol  = 0;
    if (bMass) {
      atoms = &top->atoms;
      snew(this->mass,atoms->nr);
      for(i=0; (i<atoms->nr); i++) {
	this->mass[i] = atoms->atom[i].m;
      }
    }
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

static void corr_print(t_corr *this,char *fn,char *title,char *yaxis,
		       real msdtime,real beginfit,real endfit,
		       real *DD,real *SigmaD,char *grpname[])
{
  FILE *out;
  int  i,j;
  
  out=xvgropen(fn,title,xvgr_tlabel(),yaxis);
  if (DD) {
    fprintf(out,"# MSD gathered over %g %s with %d restarts\n",
	    msdtime,time_unit(),this->nrestart);
    fprintf(out,"# Diffusion constants fitted from time %g to %g %s\n",
	    beginfit,endfit,time_unit());
    for(i=0; i<this->ngrp; i++) 
      fprintf(out,"# D[%10s] = %.4f (+/- %.4f) (1e-5 cm^2/s)\n",
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
		      bool bRmCOMM,rvec com,t_calc_func *calc1)
{
  int  nx0;
  real g;
  rvec dcom;
  
  /* Check for new starting point */
  if (this->nlast < this->nrestart) {
    if ((thistime(this) >= (this->nlast*this->delta_t)) && (nr==0)) {
      memcpy(this->x0[this->nlast],xc,this->ncoords*sizeof(xc[0]));
      this->n_offs[this->nlast]=this->nframes;
      copy_rvec(com,this->com[this->nlast]);
      this->nlast++;
    }
  }

  /* nx0 appears to be the number of new starting points,
   * so for all starting points, call calc1. 
   */
  for(nx0=0; (nx0<this->nlast); nx0++) {
    if (bRmCOMM) {
      rvec_sub(this->com[nx0],com,dcom);
    } else {
      clear_rvec(dcom);
    }
    g = calc1(this,nx,index,nx0,xc,dcom);
#ifdef DEBUG2
    printf("g[%d]=%g\n",nx0,g);
#endif
    this->data [nr][in_data(this,nx0)]+=g;
    this->ndata[nr][in_data(this,nx0)]++;
  }
}

static real calc1_norm(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[],
		       rvec dcom)
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
	r   = this->x0[nx0][ix][m] - xc[ix][m] - dcom[m];
	r2 += r*r;
      }
      break;
    case X:
    case Y:
    case Z:
      r = this->x0[nx0][ix][this->type-X] - xc[ix][this->type-X]
	- dcom[this->type-X];
      r2 += r*r;
      break;
    case LATERAL:
      for(m=0; (m<DIM); m++) {
	if (m != this->axis) {
	  r   = this->x0[nx0][ix][m] - xc[ix][m] - dcom[m];
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

static void calc_mol_com(t_block *mols,t_atoms *atoms,rvec *x,rvec *xa)
{
  int  mol,i,j,d;
  rvec xm;
  real m,mtot;

  for(mol=0; mol<mols->nr; mol++) {
    clear_rvec(xm);
    mtot = 0;
    for(i=mols->index[mol]; i<mols->index[mol+1]; i++) {
      j = mols->a[i];
      m = atoms->atom[j].m;
      for(d=0; d<DIM; d++)
	xm[d] += m*x[j][d];
      mtot += m;
    }
    svmul(1/mtot,xm,xa[mol]);
  }
}

static real calc_one_mw(t_corr *this,int ix,int nx0,rvec xc[],real *tm,
			rvec dcom)
{
  real r2,r,mm;
  int  m;
  
  mm=this->mass[ix];
  if (mm == 0)
    return 0;
  (*tm) += mm;
  r2     = 0.0;
  switch (this->type) {
  case NORMAL:
    for(m=0; (m<DIM); m++) {
      r   = this->x0[nx0][ix][m] - xc[ix][m] - dcom[m];
      r2 += mm*r*r;
    }
    break;
  case X:
  case Y:
  case Z:
    r  = this->x0[nx0][ix][this->type-X] - xc[ix][this->type-X]
      - dcom[this->type-X];
    r2 = mm*r*r;
      break;
  case LATERAL:
    for(m=0; (m<DIM); m++) {
      if (m != this->axis) {
	r   = this->x0[nx0][ix][m] - xc[ix][m] - dcom[m];
	r2 += mm*r*r;
      }
    }
    break;
  default:
    gmx_fatal(FARGS,"Options got screwed. Did not expect value %d\n",this->type);
  } /* end switch */
  return r2;
}

static real calc1_mw(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[],
		     rvec dcom)
{
  int  i;
  real g,tm;
  
  g=tm=0.0;
  for(i=0; (i<nx); i++) 
    g += calc_one_mw(this,index[i],nx0,xc,&tm,dcom);
  
  g/=tm;
  
  return g;
}

static void remove_pbc(int natoms,t_atoms *atoms,
		       rvec xcur[],rvec xprev[],matrix box,
		       rvec com)
{
  int  i,m;
  rvec hbox;
  real mass;
  dvec sx;
  double tmass;

  /* Remove periodicity */
  for(m=0; (m<DIM); m++)
    hbox[m] = 0.5*box[m][m];
  clear_dvec(sx);
  tmass = 0;
  for(i=0; (i<natoms); i++) {
      for(m=DIM-1; m>=0; m--) {
	while(xcur[i][m]-xprev[i][m] <= -hbox[m])
	  rvec_inc(xcur[i],box[m]);
	while(xcur[i][m]-xprev[i][m] >  hbox[m])
	  rvec_dec(xcur[i],box[m]);
    }
    mass = atoms->atom[i].m;
    for(m=0; m<DIM; m++)
      sx[m] += mass*xcur[i][m];
    tmass += mass;
  }
  for(m=0; m<DIM; m++)
    com[m] = sx[m]/tmass;
}

static void prep_data(int gnx,atom_id index[],
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

static real calc1_mol(t_corr *this,int nx,atom_id index[],int nx0,rvec xc[],
		      rvec dcom)
{
  int  i;
  real g,mm,gtot,tt;

  tt = this->time[in_data(this,nx0)];
  gtot = 0;
  for(i=0; (i<nx); i++) {
    mm = 0;
    g = calc_one_mw(this,index[i],nx0,xc,&mm,dcom);
    /* We don't need to normalize as the mass was set to 1 */
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
  real  a,b,r,*D,Dav,D2av,VarD;
  
  out=xvgropen(fn,"Diffusion Coefficients / Molecule","Molecule","D");
  
  snew(D,this->nmol);
  Dav = D2av = 0;
  for(i=0; (i<this->nmol); i++) {
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
  Dav  /= this->nmol;
  D2av /= this->nmol;
  VarD  = D2av - sqr(Dav);
  printf("<D> = %.4f Std. Dev. = %.4f Error = %.4f\n",
	 Dav,sqrt(VarD),sqrt(VarD/this->nmol));
  
  sfree(D);
}

/* this is the mean loop for the correlation type functions 
 * fx and nx are file pointers to things like read_first_x and
 * read_next_x
 */
void corr_loop(t_corr *this,char *fn,t_topology *top,
	       int nmol,int gnx[],atom_id *index[],
	       t_calc_func *calc1,bool bRmCOMM,real dt)
{
  rvec         *x[2],*xa[2],com;
  real         t;
  int          natoms,i,j,status,cur=0,maxframes=0;
#define        prev (1-cur)
  matrix       box;
  bool         bFirst;

  natoms = read_first_x(&status,fn,&this->t0,&(x[cur]),box);
#ifdef DEBUG
  fprintf(stderr,"Read %d atoms for first frame\n",natoms);
#endif
  if (bRmCOMM && natoms < top->atoms.nr)
    fprintf(stderr,"WARNING: The trajectory only contains part of the system (%d of %d atoms) and therefore the COM motion of only this part of the system will be removed\n",natoms,top->atoms.nr);

  snew(x[prev],natoms);

  if (nmol > 0) {
    this->ncoords = nmol;
    snew(xa[0],nmol);
    snew(xa[1],nmol);
  } else {
    this->ncoords = natoms;
    xa[0] = x[0];
    xa[1] = x[1];
  }

  bFirst = TRUE;
  t=this->t0;
  do {
    if (bRmod(t,this->t0,dt)) {
      this->nrestart++;
  
      srenew(this->x0,this->nrestart);
      snew(this->x0[this->nrestart-1],this->ncoords);
      srenew(this->com,this->nrestart);
      srenew(this->n_offs,this->nrestart);
      srenew(this->lsq,this->nrestart);
      snew(this->lsq[this->nrestart-1],this->nmol);
      if (debug)
	fprintf(debug,"Extended data structures because of new restart %d\n",
		this->nrestart);
    }
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

    if (bFirst) {
      memcpy(xa[prev],xa[cur],this->ncoords*sizeof(xa[prev][0]));
      bFirst = FALSE;
    }

    if (bRmCOMM)
      remove_pbc(natoms,&top->atoms,xa[cur],xa[prev],box,com);

    if (nmol) {
      /*rm_pbc(&top->idef,natoms,box,x[cur],x[cur]);*/
      calc_mol_com(&top->blocks[ebMOLS],&top->atoms,x[cur],xa[cur]);
    }
    
    /* loop over all groups in index file */
    for(i=0; (i<this->ngrp); i++) {
      /* nice for putting things in boxes and such */
      if (!bRmCOMM)
	prep_data(gnx[i],index[i],xa[cur],xa[prev],box);
      /* calculate something useful, like mean square displacements */
      calc_corr(this,i,gnx[i],index[i],xa[cur],bRmCOMM,com,calc1);
    }
    cur=prev;
    
    this->nframes++;
  } while (read_next_x(status,&t,natoms,x[cur],box));
  fprintf(stderr,"\nUsed %d restart points spaced %g %s over %g %s\n\n", 
	  this->nrestart, 
	  convert_time(dt), time_unit(),
	  convert_time(this->time[this->nframes-1]), time_unit() );
  
  close_trj(status);
}

void do_corr(char *trx_file, char *ndx_file, char *msd_file, char *mol_file,
	     int nrgrp, t_topology *top,bool bMW,bool bRmCOMM,
	     int type,real dim_factor,int axis,
	     real dt,real beginfit,real endfit)
{
  t_corr       *msd;
  int          *gnx;
  atom_id      **index;
  char         **grpname;
  int          i,i0,i1,j,N;
  real         *DD,*SigmaD,a,a2,b,r;
  
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

  msd = init_corr(nrgrp,type,axis,dim_factor,(mol_file!=NULL),bMW,dt,top);
  
  corr_loop(msd,trx_file,top,mol_file ? gnx[0] : 0,gnx,index,
	    (mol_file!=NULL) ? calc1_mol : (bMW ? calc1_mw : calc1_norm),
	    bRmCOMM,dt);
  
  /* Correct for the number of points */
  for(j=0; (j<msd->ngrp); j++)
    for(i=0; (i<msd->nframes); i++)
      msd->data[j][i] /= msd->ndata[j][i];

  if (mol_file) 
    printmol(msd,mol_file);

  DD     = NULL;
  SigmaD = NULL;

  if (beginfit == -1) {
    i0 = (int)(0.1*(msd->nframes - 1) + 0.5);
    beginfit = msd->time[i0];
  } else
  for(i0=0; i0<msd->nframes && msd->time[i0]<beginfit; i0++) 
    ;

  if (endfit == -1) {
    i1 = (int)(0.9*(msd->nframes - 1) + 0.5) + 1;
    endfit = msd->time[i1-1];
  } else
    for(i1=i0; i1<msd->nframes && msd->time[i1]<=endfit; i1++)
		      ;
  fprintf(stdout,"Fitting from %g to %g %s\n\n",beginfit,endfit,time_unit());

  N = i1-i0;
  if (N <= 2) {
    fprintf(stdout,"Not enough points for fitting (%d).\n"
	    "Can not determine the diffusion constant.\n",N);
  } else {
    snew(DD,msd->ngrp);
    snew(SigmaD,msd->ngrp);
    for(j=0; j<msd->ngrp; j++) {
      if (N >= 4) {
	lsq_y_ax_b(N/2,&(msd->time[i0]),&(msd->data[j][i0]),&a,&b,&r);
	lsq_y_ax_b(N/2,&(msd->time[i0+N/2]),&(msd->data[j][i0+N/2]),&a2,&b,&r);
	SigmaD[j] = fabs(a-a2);
      } else
	SigmaD[j] = 0;
      lsq_y_ax_b(N,&(msd->time[i0]),&(msd->data[j][i0]),&(DD[j]),&b,&r);
      DD[j]     *= FACTOR/msd->dim_factor;
      SigmaD[j] *= FACTOR/msd->dim_factor;
      fprintf(stdout,"D[%10s] %.4f (+/- %.4f) 1e-5 cm^2/s\n",
	      grpname[j],DD[j],SigmaD[j]);
    }
  }
  /* Print mean square displacement */
  corr_print(msd,msd_file,
	     "Mean Square Displacement",
	     "MSD (nm\\S2\\N)",
	     msd->time[msd->nframes-1],beginfit,endfit,DD,SigmaD,grpname);
}

int gmx_msd(int argc,char *argv[])
{
  static char *desc[] = {
    "g_msd computes the mean square displacement (MSD) of atoms from",
    "their initial positions. This provides an easy way to compute",
    "the diffusion constant using the Einstein relation.",
    "The time between additional starting points for the MSD calculation",
    "is set with [TT]-trestart[tt].",
    "The diffusion constant is calculated by least squares fitting a",
    "straight line through the MSD from [TT]-beginfit[tt] to",
    "[TT]-endfit[tt]. An error estimate given, which is the difference",
    "of the diffusion coefficients obtained from fits over the two halfs",
    "of the fit interval.[PAR]",
    "With option [TT]-rmcomm[tt] center of mass motion can be removed.",
    "For trajectories produced with GROMACS this is usually not necessary",
    "as mdrun usually already removes the center of mass motion.",
    "When you use this option be sure that the whole system is stored",
    "in the trajectory file.[PAR]",
    "Option [TT]-mol[tt] plots the MSD for molecules, this implies",
    "[TT]-mw[tt], i.e. for each inidividual molecule an diffusion constant",
    "is computed for its center of mass. When using an index file,",
    "it should contain molecule numbers instead of atom numbers.",
    "Using this option one also gets an accurate error estimate",
    "based on the statistics between individual molecules. Since one usually",
    "is interested in self-diffusion at infinite dilution this is probably",
    "the most useful number.[PAR]",
  };
  static char *normtype[]= { NULL,"no","x","y","z",NULL };
  static char *axtitle[] = { NULL,"no","x","y","z",NULL };
  static int  ngroup     = 1;
  static real dt         = 10; 
  static real beginfit   = -1; 
  static real endfit     = -1; 
  static bool bMW        = TRUE;
  static bool bRmCOMM    = FALSE;
  t_pargs pa[] = {
    { "-type",    FALSE, etENUM, {normtype},
      "Compute diffusion coefficient in one direction" },
    { "-lateral", FALSE, etENUM, {axtitle}, 
      "Calculate the lateral diffusion in a plane perpendicular to" },
    { "-ngroup",  FALSE, etINT,  {&ngroup},
      "Number of groups to calculate MSD for" },
    { "-mw",      FALSE, etBOOL, {&bMW},
      "Mass weighted MSD" },
    { "-rmcomm",      FALSE, etBOOL, {&bRmCOMM},
      "Remove center of mass motion" },
    { "-trestart",FALSE, etTIME, {&dt},
      "Time between restarting points in trajectory (%t)" },
    { "-beginfit",FALSE, etTIME, {&beginfit},
      "Start time for fitting the MSD (%t), -1 is 10%" },
    { "-endfit",FALSE, etTIME, {&endfit},
      "End time for fitting the MSD (%t), -1 is 90%" }
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
  bTop = read_tps_conf(tps_file,title,&top,&xdum,NULL,box,bMW || bRmCOMM); 
  if (mol_file && !bTop)
    gmx_fatal(FARGS,"Could not read a topology from %s. Try a tpr file instead.",
		tps_file);
    
  do_corr(trx_file,ndx_file,msd_file,mol_file,ngroup,
	  &top,bMW,bRmCOMM,type,dim_factor,axis,dt,beginfit,endfit);
  
  view_all(NFILE, fnm);
  
  thanx(stderr);
  
  return 0;
}
