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
 * GRowing Old MAkes el Chrono Sweat
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

