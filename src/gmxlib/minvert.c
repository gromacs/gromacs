/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Gyas ROwers Mature At Cryogenic Speed
 */

/* This file is completely threadsafe - keep it that way! */

#include "minvert.h"
#include "nr.h"
#include "smalloc.h"

void mat_mult(int n,real **A,real **B,real **C)
{
  int i,j,k;

  for(i=1; (i<=n); i++)
    for(j=1; (j<=n); j++) {
      C[i][j]=0.0;
      for(k=1; (k<=n); k++)
	C[i][j]+=A[i][k]*B[k][j];
    }
}

real **mk_mat(int n)
{
  real **a;
  int  i;

  snew(a,n+1);
  for(i=0; (i<=n); i++)
    snew(a[i],n+1);

  return a;
}

real **mk_mat2(int nrow,int ncol)
{
  real **a;
  int  i;

  snew(a,nrow+1);
  for(i=0; (i<=nrow); i++)
    snew(a[i],ncol+1);

  return a;
}

void cp_mat(int n,real **src,real **dest)
{
  int i,j;
  
  for(i=1; (i<=n); i++)
    for(j=1; (j<=n); j++)
      dest[i][j]=src[i][j];
}

void print_mat(FILE *fp,char *title,int n,real **a,int *indx)
{
  int i,j,ii;
  
  fprintf(fp,"%s\n",title);
  for(i=1; (i<=n); i++) {
    if (indx)
      ii=indx[i];
    else
      ii=i;
    for(j=1; (j<=n); j++)
      fprintf(fp," %8.3f",a[ii][j]);
    fprintf(fp,"\n");
  }
}

void invert_mat(int n,real **A,real **Ainv)
{
  real *col=NULL;
  int  *indx=NULL;
  
  real d;
  int  i,j;
  
  snew(col,n+1);
  snew(indx,n+1);
  
  ludcmp(A,n,indx,&d);
  
  /* A is an LU-decomposed matrix, now back substitute to get Ainv */
  for(j=1; (j<=n); j++) {
    for(i=1; (i<=n); i++)
      col[i]=0;
    col[j]=1.0;
    lubksb(A,n,indx,col);
    for(i=1; (i<=n); i++)
      Ainv[i][j]=col[i];
  }
  sfree(col);
  sfree(indx);
}

