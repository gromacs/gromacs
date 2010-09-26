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
 * GROwing Monsters And Cloning Shrimps
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "network.h"
#include "vec.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "main.h"
#include "calcmu.h"

void calc_mu(int start,int homenr,rvec x[],real q[],real qB[],
	     int nChargePerturbed,
	     dvec mu,dvec mu_B)
{
  int i,end,m;
  
  end   = start + homenr;  
  
  clear_dvec(mu);
  for(i=start; (i<end); i++)
    for(m=0; (m<DIM); m++)
      mu[m] += q[i]*x[i][m];
  
  for(m=0; (m<DIM); m++)
    mu[m] *= ENM2DEBYE;
  
  if (nChargePerturbed) {
    clear_dvec(mu_B);
    for(i=start; (i<end); i++)
      for(m=0; (m<DIM); m++)
	mu_B[m] += qB[i]*x[i][m];
    
    for(m=0; (m<DIM); m++)
      mu_B[m] *= ENM2DEBYE;
  } else {
    copy_dvec(mu,mu_B);
  }
}

gmx_bool read_mu(FILE *fp,rvec mu,real *vol)
{
  /* For backward compatibility */
  real mmm[4];
  
  if (fread(mmm,(size_t)(4*sizeof(real)),1,fp) != 1)
    return FALSE;
    
  copy_rvec(mmm,mu);
  *vol = mmm[3];
  
  return TRUE;
}

