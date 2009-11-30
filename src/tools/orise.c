/*
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


#include "typedefs.h"
#include "maths.h"
#include "string2.h"
#include "hxprops.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "smalloc.h"
#include "readev.h"
#include "macros.h"
#include "confio.h"
#include "copyrite.h"
#include "statutil.h"
#include "mcprop.h"
#include "orise.h"

void optimize(FILE *fp,int nx,t_propfunc *f,t_pinp *p)
{
  real *fact,sf;
  int  i;
  
  snew(fact,nx);
  sf=0.0;
  for(i=0; (i<nx); i++) {
    fact[i]=1.0/(i+2.0);
    sf+=fact[i]*fact[i];
  }
  for(i=0; (i<nx); i++) 
    fact[i]/=sqrt(sf);

  do_mc(fp,nx,fact,p->step,p->v0,p->tol,p->nsteps,f);
  
  fprintf(fp,"MC Done\n");
  fprintf(fp,"Components are:\n");
  sf=0.0;
  for(i=0; (i<nx); i++) {
    fprintf(fp,"EV %5d:  Fac: %8.3f\n",i,fact[i]);
    sf+=fact[i]*fact[i];
  }
  fprintf(fp,"Norm of vector: %8.3f\n",sqrt(sf));
}

static rvec    *xav;
static rvec    **EV;
static real    **evprj;
static int     nca,natoms,nframes;
static atom_id *ca_index;

real risefunc(int nx,real x[])
{
  static rvec *xxx=NULL;
  real   cx,cy,cz;
  double rrr,rav2,rav;
  int    i,j,m,n,ai;
  
  if (xxx == NULL) 
    snew(xxx,natoms);
  
  rav2=0,rav=0;
  
  for(j=0; (j<nframes); j++) {
    /* Make a structure, we only have to do Z */
    for(i=0; (i<nca); i++) {
      ai=ca_index[i];
      /*cx=xav[ai][XX];
      cy=xav[ai][YY];*/
      cz=xav[ai][ZZ];
      for(n=0; (n<nx); n++) {
	/*cx+=EV[n][ai][XX]*x[n]*evprj[n][j];
	cy+=EV[n][ai][YY]*x[n]*evprj[n][j];*/
	cz+=EV[n][ai][ZZ]*x[n]*evprj[n][j];
      }
      /*xxx[ai][XX]=cx;
      xxx[ai][YY]=cy;*/
      xxx[ai][ZZ]=cz;
    }
    rrr   = rise(nca,ca_index,xxx);
    rav  += rrr;
    rav2 += rrr*rrr;
  }
  rav/=nframes;
  rrr=sqrt(rav2/nframes-rav*rav);
  
  return -rrr;
}

real radfunc(int nx,real x[])
{
  static rvec *xxx=NULL;
  real   cx,cy,cz;
  double rrr,rav2,rav;
  int    i,j,m,n,ai;
  
  if (xxx == NULL) 
    snew(xxx,natoms);
  
  rav2=0,rav=0;
  
  for(j=0; (j<nframes); j++) {
    /* Make a structure, we only have to do X & Y */
    for(i=0; (i<nca); i++) {
      ai=ca_index[i];
      cx=xav[ai][XX];
      cy=xav[ai][YY];
      cz=xav[ai][ZZ];
      for(n=0; (n<nx); n++) {
	cx+=EV[n][ai][XX]*x[n]*evprj[n][j];
	cy+=EV[n][ai][YY]*x[n]*evprj[n][j];
	cz+=EV[n][ai][ZZ]*x[n]*evprj[n][j];
      }
      xxx[ai][XX]=cx;
      xxx[ai][YY]=cy;
      xxx[ai][ZZ]=cz;
    }
    rrr   = radius(NULL,nca,ca_index,xxx);
    rav  += rrr;
    rav2 += rrr*rrr;
  }
  rav/=nframes;
  rrr=sqrt(rav2/nframes-rav*rav);
  
  return -rrr;
}

void init_optim(int nx,rvec *xxav,rvec **EEV,
		real **eevprj,int nnatoms,
		int nnca,atom_id *cca_index,
		t_pinp *p)
{
  xav      = xxav;
  EV       = EEV;
  evprj    = eevprj;
  natoms   = nnatoms;
  nframes  = p->nframes;
  nca      = nnca;
  ca_index = cca_index;
}

void optim_rise(int nx,rvec *xxav,rvec **EEV,
		real **eevprj,int nnatoms,
		int nnca,atom_id *cca_index,
		t_pinp *p)
{
  FILE *fp;
  
  fp       = ffopen("rise.log","w");
  init_optim(nx,xxav,EEV,eevprj,nnatoms,nnca,cca_index,p);
  optimize(fp,nx,risefunc,p);
  ffclose(fp);
}

void optim_radius(int nx,rvec *xxav,rvec **EEV,
		  real **eevprj,int nnatoms,
		  int nnca,atom_id *cca_index,
		  t_pinp *p)
{
  FILE *fp;
  
  fp       = ffopen("radius.log","w");
  
  init_optim(nx,xxav,EEV,eevprj,nnatoms,nnca,cca_index,p);
  optimize(fp,nx,radfunc,p);
  ffclose(fp);
}

