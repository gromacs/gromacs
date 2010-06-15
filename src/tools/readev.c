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
  ffclose(in);
  
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
    ffclose(in);
  }
  fprintf(stderr,"\rSuccessfully read eigenvector projections\n");
  
  return evprj;
}
