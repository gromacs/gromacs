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
 * Gromacs Runs On Most of All Computer Systems
 */
#include <stdio.h>
#include "gstat.h"

static int     n=0;
#define FPL 10

void write_gms(FILE *fp,int natoms,rvec x[],matrix box)
{
  int i,j;

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
  int i,j;

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

