/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_virtual_c = "$Id$";

#include "cdist.h"

int set_virtual (int *ATMS,int N,real margin,t_dist *d,int natoms)
{
  /* Routine to add distances to virtual particle. 
     The virtual particle is placed 10A outside
     the plane common to all atoms, straight above
     the center of mass, which is calculated assuming 
     unit mass for all atoms.

     Adam Kirrander 990211                        */
  
  int ndist=0,i,j;
  real CONST=0.0,Ki,tmp,len,lb,ub;
  
  /* Calculate the constant part */
  for (i=1 ; i<N ; i++ ) {
    for (j=0 ; j<i ; j++) {
      tmp = d_len(d,natoms,ATMS[i],ATMS[j]);
      CONST += tmp*tmp;
    }
  }
  CONST = CONST/N;
  
  /* Calculate and set distances */
  for (i=0 ; i<N ; i++ ) {
    Ki = 0.0;
    for (j=0 ; j<N ; j++) {                      /*Calc. variable part*/
      if (i == j) continue;
      tmp = d_len(d,natoms,ATMS[i],ATMS[j]); 
      Ki += tmp*tmp;
    }
    len = sqrt(64.0+((Ki-CONST)/N));              /*Pythagoras*/
    lb  = (1.0-margin)*len;
    ub  = (1.0+margin)*len; 
    set_dist(d,natoms,ATMS[i],ATMS[N],lb,ub,len); 
    ndist++;
  }
  
  /* If number of virtual dist. should correspond to nr. atoms */
  if (ndist != N) fprintf(stderr,"Check routine set_virtual!\n");
  
  return ndist;
} 

