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
 * Grunge ROck MAChoS
 */
static char *SRCID_ql77_c = "$Id$";

#include <math.h>
#include "typedefs.h"
#include "vec.h"

void ql77 (int n,real **x,real *d,real *e,int nmax)
     /*c
       c     matrix diagonalization routine
       c
       c     code derived from eispack
       c     update : fmp 16.2.86 (conversion to fortran 77)
       c     update : fmp 09.3.89 (get rid of double precision)
       c
       c
       c     input
       c
       c     n           order of the matrix, not larger than nmax
       c     x(n,n)      the real symmetric matrix to be diagonalized as a
       c                 two-indexed square array of which only the lower left
       c                 triangle is used and need be supplied
       c     nmax        maximum value of n as declared in the calling routine
       c                 dimension statement
       c     e           working space of length nmax
       c
       c     output
       c
       c     d(n)        computed eigenvalues in descending order 
       c     x(n,n)      eigenvectors (in columns) in the same order as the
       c                 eigenvalues
       c
       */
{
  int i,j,k,l,j1,ni;
  real h,g,f,b,s,p,r,c,absp;
  
  const real eps=7.e-14,tol=1.e-70;
  /*c
    c     handle special case                                               
    c*/
  
  if (n == 1) {
    d[1]=x[1][1];
    x[1][1]=1.0;
    return;
  }
  /*c
    c     householder reduction to tridiagonal form                         
    c*/
    
  for(ni=2; (ni <= n); ni++) {
    i=n+2-ni;
    l=i-2;
    h=0.0;                                                          
    g=x[i][i-1];
    if (l <= 0) 
      goto L140;
    for(k=1; (k<=l); k++) 
      h=h+sqr(x[i][k]);
    s=h+g*g;
    if (s < tol) {
      h=0.0;
      goto L140;
    }
    if (h <= 0.0) 
      goto L140;
    
    l=l+1;
    f=g;
    g=sqrt(s);
    /*if (f.le.0.0)goto 70*/
    g=-g;
    h=s-f*g;                                                           
    x[i][i-1]=f-g;
    f=0.0;
    for(j=1; (j<=l); l++) {
      x[j][i]=x[i][j]/h;
      s=0.0;
      for(k=1; (k<=j); k++)
	s=s+x[j][k]*x[i][k];
      j1=j+1;
      if (j1 <= l) 
	for(k=j1; (k<=l); k++)
	  s=s+x[k][j]*x[i][k];
      e[j]=s/h;
      f=f+s*x[j][i];
    }
    f=f/(h+h);
    for(j=1; (j<=l); j++) 
      e[j]=e[j]-f*x[i][j];
    for(j=1; (j<=l); j++) {
      f=x[i][j];
      s=e[j];
      for(k=1; (k<=j); k++)
	x[j][k]=x[j][k]-f*e[k]-x[i][k]*s;
    }
  L140:
    d[i]=h;
    e[i-1]=g;
  }
  /*c
    c     accumulation of transformation matrix and intermediate d vector
    c*/
  
  d[1]=x[1][1];
  x[1][1]=1.0;
  for(i=2; (i<=n); i++) {
    l=i-1;
    if (d[i] > 0.0) {
      for(j=1; (j<=l); j++) {
	s=0.0;
	for(k=1; (k<=l); k++) 
	  s=s+x[i][k]*x[k][j];
	for(k=1; (k<=l); k++)
	  x[k][j]=x[k][j]-s*x[k][i];
      }
      d[i]=x[i][i];
    }
    x[i][i]=1.0;
    for(j=1; (j<=l); j++) {
      x[i][j]=0.0;
      x[j][i]=0.0;
    }
  }
  /*
    c
    c     ql iterates
    c
    */
  b=0.0;
  f=0.0;
  e[n]=0.0;
  for(l=1; (l<=n); l++) {
    h=eps*(fabs(d[l])+fabs(e[l]));
    if (h > b) 
      b=h;                                                   
    for(j=l; (j<=n); j++) {
      if(fabs(e[j]) <= b) 
	break;
    }
    if(j == l) 
      continue;
  L260:
    g=d[l];
    p=(d[l+1]-g)*0.5/e[l];
    r=sqrt(p*p+1.0);
    if(p < 0.0)
      p=p-r;
    else
      p=p+r;
    d[l]=e[l]/p;
    h=g-d[l];                                                     
    k=l+1;
    for(i=k; (i<=n); i++)
      d[i]=d[i]-h;                                                       
    f=f+h;                                                             
    p=d[j];
    c=1.0;
    s=0.0;
    j1=j-1;
    for(ni=l; (ni<=j1); ni++) {
      i=l+j1-ni;
      g=c*e[i];
      h=c*p;
      if(fabs(p) >= fabs(e[i])) {
	c=e[i]/p;
	r=sqrt(c*c+1.0);
	e[i+1]=s*p*r;
	s=c/r;
	c=1.0/r;
      }
      else {
	c=p/e[i];
	r=sqrt(c*c+1.0);
	e[i+1]=s*e[i]*r;
	s=1.0/r;
	c=c/r;
      }
      p=c*d[i]-s*g;
      d[i+1]=h+s*(c*g+s*d[i]);
      for(k=1; (k<=n); k++) {
	h=x[k][i+1];
	x[k][i+1]=x[k][i]*s+h*c;
	x[k][i]=x[k][i]*c-h*s;
      }
      e[l]=s*p;
      d[l]=c*p;
      if(fabs(e[l]) > b) 
	goto L260;
      d[l]=d[l]+f;
    }
  }
  /*
    c
    c**** put eigenvalues and eigenvectors in 
    c**** desired ascending order
    c
    */
  ni = n-1;
  
  for(i=1; (i<=ni); i++) {
    k    = i;
    p    = d[i];
    absp = fabs(d[i]);
    j1   = i+1;
    
    for(j=j1; (j<=n); j++) {
      if(fabs(d[j]) < absp) {
	k    = j;
	p    = d[j];
	absp = fabs(d[j]);
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for(j=1; (j<=n); j++) {
	p      = x[j][i];
	x[j][i] = x[j][k];
	x[j][k] = p;
      }
    }
  }
  
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
