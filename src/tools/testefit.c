/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * S  C  A  M  O  R  G
 */
#include <sysstuff.h>
#include <string.h>
#include <math.h>
#include "smalloc.h"
#include "vec.h"
#include "statutil.h"
#include "futil.h"
#include "rdgroup.h"
#include "copyrite.h"
#include "typedefs.h"
#include "statusio.h"
#include "gstat.h"

void main()
{
  real *x,*y,*Dy;
  real a,b,sa,sb,A,tau,dA,dtau;
  int  i,n;
  FILE *out;
  
  printf("Give the number of points : ");
  scanf("%d",&n);
  printf("Give A                    : ");
  scanf("%f",&A);
  printf("Give tau                  : ");
  scanf("%f",&tau);
  printf("\n");
  
  printf("\tA = \t%10.5f\ttau = \t%10.5f\n",A,tau);
  /*  n=1000; */
  snew(x,n);
  snew(y,n);
  snew(Dy,n);
  
  out=ffopen("explist.dat","w");
  for(i=0;(i<n);i++){
    x[i]=0.01*i;
    y[i]=A*exp(-x[i]/tau);
    /*    Dy[i]=1.0/sqrt(n-i); */
    Dy[i]=1.0;
    fprintf(out,"%f %f\n",x[i],y[i]);
  }
  fclose(out);
  
  printf("\nhallo world\n");
  expfit(n,x,y,Dy,&a,&sa,&b,&sb);
  A=exp(a);
  dA=exp(a)*sa;
  tau=-1.0/b;
  dtau=sb/sqr(b);

  printf("Fitted to y=exp(a+bx):\n");
  printf("a = %10.5f\t b = %10.5f",a,b);
  printf("\n");
  printf("Fitted to y=Aexp(-x/tau):\n");
  printf("A  = %10.5f\t tau  = %10.5f\n",A,tau);
  printf("dA = %10.5f\t dtau = %10.5f\n",dA,dtau);
}




