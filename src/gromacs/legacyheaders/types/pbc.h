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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifndef _pbc_h
#define _pbc_h


#include "simple.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Maximum number of combinations of single triclinic box vectors
 * required to shift atoms that are within a brick of the size of
 * the diagonal of the box to within the maximum cut-off distance.
 */
#define MAX_NTRICVEC 12

typedef struct {
  int    ndim_ePBC;
  int    ePBCDX;
  int    dim;
  matrix box;
  rvec   fbox_diag;
  rvec   hbox_diag;
  rvec   mhbox_diag;
  real   max_cutoff2;
  gmx_bool   bLimitDistance;
  real   limit_distance2;
  int    ntric_vec;
  ivec   tric_shift[MAX_NTRICVEC];
  rvec   tric_vec[MAX_NTRICVEC];
} t_pbc;

#ifdef __cplusplus
}
#endif

#endif
