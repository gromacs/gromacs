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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"

#define MDIM 18

real **read_shifts(char *fn)
{
  FILE   *fp;
  real   **ptr;
  int    i,j;
  double d;
  
  fp=libopen(fn,"r");
  snew(ptr,MDIM);
  for(i=0; (i<MDIM); i++) {
    snew(ptr[i],MDIM);
    for(j=0; (j<MDIM); j++) {
      fscanf(fp,"%lf",&d);
      ptr[i][j] = d;
    }
  }
  fclose(fp);
  
  return ptr;
}

void calc_shifts(real phi,real psi,real *ca,real *cb,real *ha,real *co)
{
  static bool bFirst=TRUE;
  static real **ca_shifts,**cb_shifts,**ha_shifts,**co_shifts;
  real   phi_ind,psi_ind;
  int    iphi,ipsi;
  
  if (bFirst) {
    ca_shifts=read_shifts("ca-shifts.dat");
    cb_shifts=read_shifts("cb-shifts.dat");
    ha_shifts=read_shifts("ha-shifts.dat");
    co_shifts=read_shifts("co-shifts.dat");
    bFirst = FALSE;
  }
  
  phi    += M_PI;
  psi    += M_PI;
  phi_ind = (phi*MDIM)/(2*M_PI);
  psi_ind = (psi*MDIM)/(2*M_PI);
  iphi    = phi_ind;
  ipsi    = psi_ind;
  
}
