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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_minvert_c = "$Id$";

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
  static real *col=NULL;
  static int  *indx=NULL;
  static int  N=0;
  
  real d;
  int  i,j;
  
  if (n > N) {
    srenew(col,n+1);
    srenew(indx,n+1);
    N=n;
  }
  ludcmp(A,N,indx,&d);
  
  /* A is an LU-decomposed matrix, now back substitute to get Ainv */
  for(j=1; (j<=n); j++) {
    for(i=1; (i<=n); i++)
      col[i]=0;
    col[j]=1.0;
    lubksb(A,n,indx,col);
    for(i=1; (i<=n); i++)
      Ainv[i][j]=col[i];
  }
}

