/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Gyas ROwers Mature At Cryogenic Speed
 */
static char *SRCID_nrjac_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gstat.h"

/* These routines were taken from NUMERICAL RECIPES */

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
  a[k][l]=h+s*(g-h*tau);
	
static void nrerror(char *error_text)
{
  fprintf(stderr,"Numerical Recipes run_time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to the system...\n");
  exit(1);
}

static void free_vector(double *v,int nl,int nh)
{
  free((char*)(v+nl));
}

static double *vector(int nl,int nh)
{
  double *v;
  v=(double*)malloc((unsigned)(nh-nl+1)*sizeof(double));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}

void jacobi(double a[7][7],int n,double d[7],double v[7][7],int *nrot)
{
  int j,i;
  int iq_1,ip_1;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=vector(1,n);
  z=vector(1,n);
  for (ip_1=1;ip_1<=n;ip_1++) {
    for (iq_1=1;iq_1<=n;iq_1++) v[ip_1][iq_1]=0.0;
    v[ip_1][ip_1]=1.0;
  }
  for (ip_1=1;ip_1<=n;ip_1++) {
    b[ip_1]=d[ip_1]=a[ip_1][ip_1];
    z[ip_1]=0.0;
  }
  *nrot=0;
  for (i=1; (i<=50); i++) {
    sm=0.0;
    for (ip_1=1;ip_1<=n-1;ip_1++) {
      for (iq_1=ip_1+1;iq_1<=n;iq_1++)
        sm += fabs(a[ip_1][iq_1]);
    }
    if (sm == 0.0) {
      free_vector(z,1,n);
      free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip_1=1;ip_1<=n-1;ip_1++) {
      for (iq_1=ip_1+1;iq_1<=n;iq_1++) {
        g=100.0*fabs(a[ip_1][iq_1]);
        if (i > 4 && fabs(d[ip_1])+g == fabs(d[ip_1])
            && fabs(d[iq_1])+g == fabs(d[iq_1]))
          a[ip_1][iq_1]=0.0;
        else if (fabs(a[ip_1][iq_1]) > tresh) {
          h=d[iq_1]-d[ip_1];
          if (fabs(h)+g == fabs(h))
            t=(a[ip_1][iq_1])/h;
          else {
            theta=0.5*h/(a[ip_1][iq_1]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip_1][iq_1];
          z[ip_1] -= h;
          z[iq_1] += h;
          d[ip_1] -= h;
          d[iq_1] += h;
          a[ip_1][iq_1]=0.0;
          for (j=1;j<=ip_1-1;j++) {
            ROTATE(a,j,ip_1,j,iq_1)
	  }
          for (j=ip_1+1;j<=iq_1-1;j++) {
            ROTATE(a,ip_1,j,j,iq_1)
            }
          for (j=iq_1+1;j<=n;j++) {
            ROTATE(a,ip_1,j,iq_1,j)
            }
          for (j=1;j<=n;j++) {
            ROTATE(v,j,ip_1,j,iq_1)
            }
          ++(*nrot);
        }
      }
    }
    for (ip_1=1; (ip_1<=n); ip_1++) {
      b[ip_1] +=  z[ip_1];
      d[ip_1]  =  b[ip_1];
      z[ip_1]  =  0.0;
    }
  }
  nrerror("Too many iterations in routine JACOBI");
}

