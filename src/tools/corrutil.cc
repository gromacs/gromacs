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
static char *SRCID_corrutil_cc = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) corrutil.cc 1.8 31 Oct 1994"
#endif /* HAVE_IDENT */
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "statutil.h"
#include "futil.h"
#include "rdgroup.h"
#include "copyrite.h"
#include "typedefs.h"
#include "statusio.h"
#include "corrutil.h"
#include "xvgr.h"

float dim_factor = 6.0;         // 6 for 3D, 4 for 2D, 2 for 1D

c_corr::c_corr(int nrgrp,int nrframes)
{
  int i;
  
  maxframes = nrframes;
  ngrp      = nrgrp;
  nframes   = 0;
  nlast     = 0;
  
  ndata=(int **)calloc(nrgrp,sizeof(ndata[0]));
  data=(real **)calloc(nrgrp,sizeof(data[0]));
  for(i=0; (i<nrgrp); i++) {
    ndata[i]=(int *)calloc(maxframes,sizeof(ndata[0][0]));
    data[i]=(real *)calloc(maxframes,sizeof(data[0][0]));
  }
  time=(real *)calloc(maxframes,sizeof(time[0]));
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

void c_corr::normalise()
{
  int i,j;
  
  for(j=0; (j<ngrp); j++)
    for(i=0; (i<nframes); i++)
      data[j][i]/=ndata[j][i];
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
  real aa,bb,da,db,err_fac;
  
  fprintf(stderr,"new print\n");
  err_fac=sqrt(2.0/3.0*sqrt(M_PI));
  normalise();
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
  int          i,status,cur=0;
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
    time[nframes]=t;
    
    // loop over all groups in index file
    for(i=0; (i<ngrp); i++) {
      // nice for putting things in boxes and such
      prep_data(gnx[i],index[i],x[cur],x[prev],box);
      // calculate something usefull, like mean square displacements
      calc(i,gnx[i],index[i],x[cur]);
    }
    cur=prev;
    
    if (++nframes >= maxframes) {
      fprintf(stderr,"\nToo many frames in trajectory (> %d)\n",maxframes);
      break;
    }

  } while (nx(status,&t,natoms,x[cur],box));
  fprintf(stderr,"\n");
  
  close_trj(status);
}


