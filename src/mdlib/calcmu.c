/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
static char *SRCID_calcmu_c = "$Id$";
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "network.h"
#include "vec.h"
#include "fatal.h"
#include "physics.h"
#include "nsb.h"
#include "main.h"

void calc_mu_and_q(t_nsborder *nsb,rvec x[],real q[],rvec mu,real *qsum)
{
  int i,start,end,m;
  /* temporary double prec. to maintain precision */
  double tmpmu[3];
  double tmpq;
  
  start = START(nsb);
  end   = start + HOMENR(nsb);  

  tmpq=0;
  for(m=0;m<DIM;m++)
    tmpmu[m]=0;

  for(i=start; (i<end); i++) {
    for(m=0; (m<DIM); m++) {
      tmpmu[m] += q[i]*x[i][m];
    }
    tmpq+=q[i];
  }
  for(m=0; (m<DIM); m++)
    mu[m] = tmpmu[m] * ENM2DEBYE;

  *qsum=tmpq;
}

bool read_mu(FILE *fp,rvec mu,real *vol)
{
  /* For backward compatibility */
  real mmm[4];
  
  if (fread(mmm,(size_t)(4*sizeof(mmm)),1,fp) != 1)
    return FALSE;
    
  copy_rvec(mmm,mu);
  *vol = mmm[3];
  
  return TRUE;
}

