/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_calcmu_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "network.h"
#include "vec.h"
#include "fatal.h"
#include "nsb.h"

void calc_mu(t_commrec *cr,t_nsborder *nsb,rvec x[],real q[],rvec mu)
{
  int i,start,end,m;
  
  start = START(nsb);
  end   = start + HOMENR(nsb);  
  
  clear_rvec(mu);
  for(i=start; (i<end); i++)
    for(m=0; (m<DIM); m++)
      mu[m] += q[i]*x[i][m];
  if (PAR(cr)) {
    gmx_sum(DIM,mu,cr);
  }
}

void write_mu(FILE *fp,rvec mu,matrix box)
{
  real mmm[4];
  
  copy_rvec(mu,mmm);
  mmm[3] = det(box);
  if (fwrite(mmm,4*sizeof(mmm[0]),1,fp) != 1)
    fatal_error(0,"Writing mu!");
}

bool read_mu(FILE *fp,rvec mu,real *vol)
{
  real mmm[4];
  
  if (fread(mmm,4*sizeof(mmm),1,fp) != 1)
    return FALSE;
    
  copy_rvec(mmm,mu);
  *vol = mmm[3];
  
  return TRUE;
}
