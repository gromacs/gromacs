/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * S  C  A  M  O  R  G
 */
static char *SRCID_sfac_c = "$Id$";
#include <stdlib.h>
#include <math.h>
#include "assert.h"
#include "typedefs.h"
#include "complex.h"
#include "smalloc.h"
#include "network.h"
#include "txtdump.h"
#include "random.h"
#include "sfac.h"
#include "vec.h"
#include "futil.h"
#include "matio.h"

void print_hkl(char *hklfn,int nh,int nk,int nl,int l0,real Intensity[],
	       int nlevels)
{
  FILE *out;
  real *t_x,*t_y,**matrix;
  real lo,hi,val;
  t_rgb rlo = { 1.0, 1.0, 0.0 };
  t_rgb rhi = { 0.0, 0.0, 1.0 };
  int  n_x,n_y;
  int  i,j,k,hkl;

  fprintf(stderr,"nh=%d, nk=%d, nl=%d\n",nh,nk,nl);  
  n_x = (2*nh+1);
  n_y = (2*nk+1);
  
  assert ((l0>= -nl) && (l0<=nl));
  out=ffopen(hklfn,"w");
  snew(t_x,n_x);
  snew(t_y,n_y);
  for(i=-nh; (i<=nh); i++)
    t_x[i+nh] = i;
  for(i=-nk; (i<=nk); i++)
    t_y[i+nk] = i;
  lo = 1000;
  hi = -1000;
  snew(matrix,n_x);
  hkl = 0;
  for(i=-nh; (i<=nh); i++) {
    snew(matrix[i+nh],n_y);
    for(j=-nk; (j<=nk); j++) {
      if ((i == 0) && (j == 0) && (l0 == 0))
	val = 0.0;
      else
	val = Intensity[hkl+nl+l0];
      lo  = min(lo,val);
      hi  = max(hi,val);
      matrix[i+nh][j+nk] =  val;
      hkl += (2*nl+1);
    }
  }
  write_xpm(out,"HKL l=0","electrons","h","k",n_x,n_y,t_x,t_y,
	    matrix,lo,hi,rlo,rhi,&nlevels);
  for(i=0; (i<=n_x); i++)
    sfree(matrix[i]);
  sfree(matrix);
  sfree(t_x);
  sfree(t_y);
  fclose(out);
}

void sum_Fhkl(t_commrec *cr,int nhkl,t_complex Fhkl[])
{
  if (PAR(cr))
    gmx_sum(nhkl*2,(real *)Fhkl,cr);
}

void calc_Fhkl(FILE *log,int nh,int nk,int nl,
	       t_complex Fhkl[],int start,int homenr,rvec x[],
	       real fatom[])
{
  /* Calculate the structure factor F(hkl) from the coordinates 
   * x[start] .. x[start+homenr-1]
   * Fhkl is a linear array of length (2*nh+1)*(2*nk+1)*(2*nl+1)
   */
  int        h,k,l,n,n0,nhkl;
  real       rh,rk,rl,inproduct,fatomn,F000,F000_1;
  const real twopi=2.0*M_PI;
  t_complex  shkl;
  
  /* Counter into structure factor array */
  nhkl       = 0;  
  
  for(h=-nh; (h<=nh); h++) {
    rh = h;
    for(k=-nk; (k<=nk); k++) {
      rk = k;
      for(l=-nl; (l<=nl); l++) {
	rl = l;
	shkl = cnul;
	for(n=0; (n<homenr); n++) {
	  n0         = n+start;
	  fatomn     = fatom[n];
	  inproduct  = twopi*(rh*x[n0][XX] + rk*x[n0][YY] + rl*x[n0][ZZ]);
	  shkl.re   += fatomn*cos(inproduct);
	  shkl.im   += fatomn*sin(inproduct);
	}
	if ((h==0) && (k==0) && (l==0))
	  F000 = shkl.re;
	Fhkl[nhkl++] = shkl;
      }
    }
  }
  /* Normalize structrue factor */
  F000_1 = 1.0/F000;
  /*for(n=0; (n<nhkl); n++) {
    Fhkl[n].re *= F000_1;
    Fhkl[n].im *= F000_1;
  }*/
  
}

void calc_I0(int nhkl,t_complex Fhkl[],real I0[])
{
  int  i;
  
  for(i=0; (i<nhkl); i++) {
    I0[i] = (sqr(Fhkl[i].re) + sqr(Fhkl[i].im));
  }
}

real calc_xray(FILE *log,int nh,int nk,int nl,
	       t_complex Fhkl[],real I0[],real kI,
	       int start,int homenr,rvec x[],rvec f_xray[],
	       real fatom[])
{
  int  h,k,l,n,n0,nhkl;
  real rh,rk,rl,kkk,VVV,DeltaI;
  real fhkl_re,fhkl_im,fhkl2,FFF;
  real inproduct,cosrS,sinrS;
  const real twopi=2.0*M_PI;
  
  /* Now compute the potential and forces */
  nhkl = 0;
  VVV  = 0;
  for(h=-nh; (h<=nh); h++) {
    rh = h;
    for(k=-nk; (k<=nk); k++) {
      rk = k;
      for(l=-nl; (l<=nl); l++) {
	rl = l;

	fhkl_re  = Fhkl[nhkl].re;
	fhkl_im  = Fhkl[nhkl].im;
	fhkl2    = sqr(fhkl_re) + sqr(fhkl_im);
        DeltaI   = fhkl2-I0[nhkl];
        VVV     += 0.5*kI*DeltaI*DeltaI;
	kkk      = kI*DeltaI*2.0*twopi;
	
	for(n=0; (n<homenr); n++) {
	  n0        = n+start;
	  inproduct = twopi*(rh*x[n0][XX] + rk*x[n0][YY] + rl*x[n0][ZZ]);
	  cosrS     = cos(inproduct);
	  sinrS     = sin(inproduct);
	  FFF       = kkk*fatom[n]*(fhkl_re*sinrS - fhkl_im*cosrS);
	  
	  f_xray[n][XX] += rh*FFF;
	  f_xray[n][YY] += rk*FFF;
	  f_xray[n][ZZ] += rl*FFF;
	}
	nhkl++;
      }
    }
  }
  return VVV;
}

void rand_coords(int start,int homenr,rvec x[])
{
  int n,m;
  int seed=1993;
  
  for(n=start; (n<start+homenr); n++) {
    for(m=0; (m<DIM); m++)
      x[n][m] = 1.0*rando(&seed);
  }
}

real *calc_nelectrons(t_atoms *atoms)
{
  real *fatom,nel,nel_1;
  int  i;
  char c;
  
  snew(fatom,atoms->nr);
  
  nel = 0;
  for(i=0; (i<atoms->nr); i++) {
    c=toupper((*(atoms->atomname[i]))[0]);
    switch (c) {
    case 'C':
      fatom[i] = 6;
      break;
    case 'N':
      fatom[i] = 7;
      break;
    case 'O':
      fatom[i] = 8;
      break;
    case 'P':
      fatom[i] = 15;
      break;
    case 'S':
      fatom[i] = 16;
      break;
    case 'D':
      fatom[i] = 2;
      break;
    case 'H':
      fatom[i] = 1;
      break;
    default:
      fprintf(stderr,"Atom: %s, first char %c\n",*(atoms->atomname[0]),c);
      fatom[i] = 1;
    }
    nel+=fatom[i];
  }
  nel_1 = 1.0/nel;
  for(i=0; (i<atoms->nr); i++)
    fatom[i] *= nel_1;
    
  return fatom;
}

real lo_calc_sfac(FILE *log,int start,int homenr,rvec x[],
		  int nh,int nk,int nl,
		  real kI,real I0[],real *n_elec,
		  rvec f[],t_commrec *cr)
{
  static     bool bFirst = TRUE,bDump=FALSE;
  static     t_complex *Fhkl=NULL;
  static     rvec *f_xray=NULL;
  char       *ptr;
  rvec       f_xraysum;
  int        n,nhkl;
  real       VVV,mass;
  
  if (Fhkl == NULL) {
    fprintf(log,"Initiating data (%s,%d)\n",__FILE__,__LINE__);
    nhkl  = (2*nh+1)*(2*nk+1)*(2*nl+1);
    
    fprintf(log,"Going to malloc %d bytes for F(hkl)\n",
	    nhkl*sizeof(Fhkl[0]));
    snew(Fhkl,nhkl);
    snew(f_xray,homenr);
    bDump = (getenv("DUMP") != NULL);
  }
  
  /* First calculate the structure factor F(hkl) */
  calc_Fhkl(log,nh,nk,nl,Fhkl,start,homenr,x,n_elec);
  
  /* Now do global summing of F(hkl) if we are in parallel */
  sum_Fhkl(cr,nhkl,Fhkl);

  if (bFirst) {
    if (getenv("RANDOM") != NULL)
      rand_coords(start,homenr,x);
    bFirst = FALSE;
  }
  
  /* Reset imaginary forces */
  for(n=0; (n<homenr); n++) {
    clear_rvec(f_xray[n]);
  }

  /* Compute potential and forces */
  VVV = calc_xray(log,nh,nk,nl,Fhkl,I0,kI,
		  start,homenr,x,f_xray,n_elec);
  
  /* Add the forces to f array */
  clear_rvec(f_xraysum);
  for(n=0; (n<homenr); n++) {
    rvec_inc(f[start+n],f_xray[n]);
    rvec_inc(f_xraysum,f_xray[n]);
  }
  
  /* Dump the forces */
  if (bDump) {
    pr_rvecs(log,0,"f_tot",f,homenr);
    pr_rvecs(log,0,"sum   ",&f_xraysum,1);
  }
  return VVV;
}
	
real calc_sfac(FILE *log,int start,int homenr,rvec x[],
	       rvec f[],t_commrec *cr)
{
  static real *I0=NULL,*n_el=NULL,kI;
  static int  nh,nk,nl;
  static bool bXray=TRUE;
  char        *ptr;
  int         i;
  
  if (!bXray)
    return 0.0;
    
  if (I0 == NULL) {
    if ((ptr = getenv("SFAC")) == NULL) {
      bXray = FALSE;
      fprintf(stderr,"No SFAC environment variable: no refinement\n");
      return 0.0;
    }
    read_diffraction(ptr,&nh,&nk,&nl,&I0);
    if ((ptr = getenv("NELEC")) == NULL) {
      bXray = FALSE;
      fprintf(stderr,"No NELEC environment variable: no refinement\n");
      return 0.0;
    }
    read_nelectron(ptr,&i,&n_el);
    assert(i == homenr);

    fprintf(stderr,"The Intensity grid is %dx%dx%d\n",nh,nk,nl);
        
    if ((ptr = getenv("KI")) != NULL)
      kI = atof(ptr);
    else
      kI = 1.0;
    fprintf(log,"kI = %g\n",kI);
  }
  
  return lo_calc_sfac(log,start,homenr,x,nh,nk,nl,kI,I0,n_el,f,cr);
}

void write_diffraction(char *fn,int nh,int nk,int nl,real I[])
{
  FILE *out;
  int  hkl,nhkl;
  
  out=ffopen(fn,"w");
  
  fprintf(out,"%8d%8d%8d\n",nh,nk,nl);
  nhkl = (2*nh+1)*(2*nk+1)*(2*nl+1);
  for(hkl = 0; (hkl<nhkl); hkl++)
    fprintf(out,"%15e\n",I[hkl]);
    
  ffclose(out);
}

void read_diffraction(char *fn,int *nh,int *nk,int *nl,real **I)
{
  FILE   *in;
  int    hkl,nhkl;
  double d;
  
  in=ffopen(fn,"r");
  
  fscanf(in,"%8d%8d%8d\n",nh,nk,nl);
  nhkl = (2* *nh + 1)*(2* *nk + 1)*(2* *nl + 1);
  snew(*I,nhkl);
  for(hkl=0; (hkl<nhkl); hkl++) {
    fscanf(in,"%lf\n",&d);
    (*I)[hkl] = d;
  }
  
  ffclose(in);
}

void write_nelectron(char *fn,int natom,real n_el[])
{
  FILE *out;
  int  i;
  
  out=ffopen(fn,"w");
  
  fprintf(out,"%d\n",natom);
  for(i = 0; (i<natom); i++)
    fprintf(out,"%g\n",n_el[i]);
  
  ffclose(out);
}

void read_nelectron(char *fn,int *natom,real **n_el)
{
  FILE   *in;
  int    i;
  double d;
  
  in=ffopen(fn,"r");
  
  fscanf(in,"%8d\n",natom);
  snew(*n_el,*natom);
  for(i=0; (i<*natom); i++) {
    fscanf(in,"%lf\n",&d);
    (*n_el)[i] = d;
  }
  
  ffclose(in);
}

