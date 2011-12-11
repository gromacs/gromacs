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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "gstat.h"
#include "wgms.h"

/*static int     n=0;*/
#define FPL 10

void write_gms(FILE *fp,int natoms,rvec x[],matrix box)
{
  int i,j,n;

  n=0;
  for(i=0;(i<natoms);i++)
    for(j=0;(j<3);j++) {
      fprintf(fp,"%8.3f",x[i][j]);
      n++;
      if (n==FPL) {
	fprintf(fp,"\n");
	n=0;
      }
    }
  if (n != 0) 
    fprintf(fp,"\n");
  if (box != NULL)
    fprintf(fp,"%8.3f%8.3f%8.3f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
}

void write_gms_ndx(FILE *fp,int isize,atom_id index[],rvec x[],matrix box)
{
  int i,j,n;

  n=0;
  for(i=0;(i<isize);i++)
    for(j=0;(j<3);j++) {
      fprintf(fp,"%8.3f",x[index[i]][j]);
      n++;
      if (n==FPL) {
	fprintf(fp,"\n");
	n=0;
      }
    }
  if (n != 0) 
    fprintf(fp,"\n");
  if (box != NULL)
    fprintf(fp,"%8.3f%8.3f%8.3f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
}

