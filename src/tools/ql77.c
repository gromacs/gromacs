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

#include <math.h>
#include "typedefs.h"
#include "vec.h"
#include "smalloc.h"

void ql77 (int n,real *x,real *d)
     /*c
       c     matrix diagonalization routine
       c
       c     code derived from eispack
       c     update : fmp 16.2.86 (conversion to fortran 77)
       c     update : fmp 09.3.89 (get rid of double precision)
       c 
       c     update : bh  16.5.99 (now works, real C code iso f2c)     
       c
       c
       c     input
       c
       c     n           order of the matrix, not larger than nmax
       c     x(n*n)      the real symmetric matrix to be diagonalized as a
       c                 two-indexed square array of which only the lower left
       c                 triangle is used and need be supplied
       c
       c     output
       c
       c     d(n)        computed eigenvalues in descending order 
       c     x(n*n)      eigenvectors (in columns) in the same order as the
       c                 eigenvalues
       c
       */
{
  int i,j,k,l,ni;
  real *e,h,g,f,b,s,p,r,c,absp;
  real totwork,work;

#define pr_pr(a,b,c) fprintf(stderr,"\rreduction: %g%%  accumulation:  %g%%  accumulation: %g%%",a,b,c)

  const real eps=7.e-14,tol=1.e-30;

  snew(e,n);

  /*
   *  householder reduction to tridiagonal form                         
   */
  
  totwork = 0;
  for(ni=1; ni<n; ni++)
    totwork += pow(n-ni,3);

  work=0;
  pr_pr(0.,0.,0.);
  for(ni=1; (ni < n); ni++) {
    i=n-ni;
    l=i-1;
    h=0.0;                                                          
    g=x[i+n*(i-1)];
    if (l > 0) {
      for(k=0; (k<l); k++) 
	h=h+sqr(x[i+n*k]);
      s=h+g*g;
      if (s < tol)
	h=0.0;
      else if (h > 0.0) { 
	l=l+1;
	f=g;
	g=sqrt(s);
	g=-g;
	h=s-f*g;                                                           
	x[i+n*(i-1)]=f-g;
	f=0.0;
	for(j=0; (j<l); j++) {
	  x[j+n*i]=x[i+n*j]/h;
	  s=0.0;
	  for(k=0; (k<=j); k++)
	    s=s+x[j+n*k]*x[i+n*k];
	  for(k=j+1; (k<l); k++)
	    s=s+x[k+n*j]*x[i+n*k];
	  e[j]=s/h;
	  f=f+s*x[j+n*i];
	}
	f=f/(h+h);
	for(j=0; (j<l); j++) 
	  e[j]=e[j]-f*x[i+n*j];
	for(j=0; (j<l); j++) {
	  f=x[i+n*j];
	  s=e[j];
	  for(k=0; (k<=j); k++)
	    x[j+n*k]=x[j+n*k]-f*e[k]-x[i+n*k]*s;
	}
      }
    }
    d[i]=h;
    e[i-1]=g;

    work += pow(n-ni,3);
    if (ni % 5 == 0)
      pr_pr(floor(100*work/totwork+0.5),0.,0.);
  }

  /*
   *  accumulation of transformation matrix and intermediate d vector
   */
  
  d[0]=x[0];
  x[0]=1.0;

  work=0;
  pr_pr(100.,0.,0.);
  for(i=1; (i<n); i++) {
    if (d[i] > 0.0) {
      for(j=0; (j<i); j++) {
	s=0.0;
	for(k=0; (k<i); k++) 
	  s=s+x[i+n*k]*x[k+n*j];
	for(k=0; (k<i); k++)
	  x[k+n*j]=x[k+n*j]-s*x[k+n*i];
      }
    }
    d[i]=x[i+n*i];
    x[i+n*i]=1.0;
    for(j=0; (j<i); j++) {
      x[i+n*j]=0.0;
      x[j+n*i]=0.0;
    }
    work += pow(i,3);
    if (i % 5 == 0)
      pr_pr(100.,floor(100*work/totwork+0.5),0.);
  }

  /*
   *  ql iterates
   */

  b=0.0;
  f=0.0;
  e[n-1]=0.0;
  totwork += pow(n,3);
  work=0;
  pr_pr(100.,100.,0.);
  for(l=0; (l<n); l++) {
    h=eps*(fabs(d[l])+fabs(e[l]));
    if (h > b) 
      b=h;                                                   
    for(j=l; (j<n); j++) {
      if(fabs(e[j]) <= b) 
	break;
    }
    if (j != l) { 
      do {
	g=d[l];
	p=(d[l+1]-g)*0.5/e[l];
	r=sqrt(p*p+1.0);
	if(p < 0.0)
	  p=p-r;
	else
	  p=p+r;
	d[l]=e[l]/p;
	h=g-d[l];                                                     
	for(i=l+1; (i<n); i++)
	  d[i]=d[i]-h;                                                       
	f=f+h;                                                             
	p=d[j];
	c=1.0;
	s=0.0;
	for(ni=l; (ni<j); ni++) {
	  i=l+j-1-ni;
	  g=c*e[i];
	  h=c*p;
	  if(fabs(p) >= fabs(e[i])) {
	    c=e[i]/p;
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*p*r;
	    s=c/r;
	    c=1.0/r;
	  } else {
	    c=p/e[i];
	    r=sqrt(c*c+1.0);
	    e[i+1]=s*e[i]*r;
	    s=1.0/r;
	    c=c/r;
	  }
	  p=c*d[i]-s*g;
	  d[i+1]=h+s*(c*g+s*d[i]);
	  for(k=0; (k<n); k++) {
	    h=x[k+n*(i+1)];
	    x[k+n*(i+1)]=x[k+n*i]*s+h*c;
	    x[k+n*i]=x[k+n*i]*c-h*s;
	  }
	}
	e[l]=s*p;
	d[l]=c*p;
      } while (fabs(e[l]) > b); 
    }
    d[l]=d[l]+f;

    work += pow(n-l,3);
    if (l % 5 == 0)
      pr_pr(100.,100.,floor(100*work/totwork+0.5));
  }
  fprintf(stderr,"\n");

  /*
   *  put eigenvalues and eigenvectors in 
   *  desired ascending order
   */
  
  for(i=0; (i<n-1); i++) {
    k    = i;
    p    = d[i];
    absp = fabs(d[i]);
    for(j=i+1; (j<n); j++) {
      if(fabs(d[j]) < absp) {
	k    = j;
	p    = d[j];
	absp = fabs(d[j]);
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for(j=0; (j<n); j++) {
	p        = x[j+n*i];
	x[j+n*i] = x[j+n*k];
	x[j+n*k] = p;
      }
    }
  }

  sfree(e);

  /*c
    c     fmp
    c     last but not least i have to confess that there is an original    
    c     g. binsch remark on this routine: 'ql is purported to be one of
    c     the fastest and most compact routines of its kind presently known
    c     to mankind.'
    c
    
    SIC! DvdS
  
    */

}
