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
 * GROtesk MACabre and Sinister
 */
#include "readev.h"
#include "futil.h"
#include "smalloc.h"
	
rvec **read_ev(char *fn,int natoms)
{
  FILE   *in;
  rvec   **ev;
  double xx,yy,zz;
  int    i,j,k,ei,ai;

  snew(ev,DIM*natoms);
  for(i=0; (i<DIM*natoms); i++)
    snew(ev[i],natoms);
  in=ffopen(fn,"r");
  for(i=0; (i<DIM*natoms); i++) {
    for(j=k=0; (j<natoms); j++) {
      fscanf(in,"%d%d%lf%lf%lf",&ei,&ai,&xx,&yy,&zz);
      ev[i][j][XX]=xx;
      ev[i][j][YY]=yy;
      ev[i][j][ZZ]=zz;
    }
  }
  fclose(in);
  
  return ev;
}

real **read_proj(int nev,int nframes,char *base)
{
  FILE   *in;
  real   **evprj;
  char   buf[256];
  int    i,j,d;
  double x;
  
  snew(evprj,nev);
  for(i=0; (i<nev); i++) {
    fprintf(stderr,"\rEV %d",i);
    snew(evprj[i],nframes);
    sprintf(buf,"%s%d",base,i+1);
    in=ffopen(buf,"r");
    for(j=0; (j<nframes); j++) {
      fscanf(in,"%d%lf",&d,&x);
      evprj[i][j]=x;
    }
    fclose(in);
  }
  fprintf(stderr,"\rSuccesfully read eigenvector projections\n");
  
  return evprj;
}
