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

#ifndef _mdatom_h
#define _mdatom_h

#include "simple.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  real          tmassA,tmassB,tmass;
  int           nr;
  int           nalloc;
  int           nenergrp;
  gmx_bool          bVCMgrps;
  int           nPerturbed;
  int           nMassPerturbed;
  int           nChargePerturbed;
  gmx_bool          bOrires;
  real          *massA,*massB,*massT,*invmass;
  real          *chargeA,*chargeB;
  gmx_bool          *bPerturbed;
  int           *typeA,*typeB;
  unsigned short        *ptype;
  unsigned short        *cTC,*cENER,*cACC,*cFREEZE,*cVCM;
  unsigned short        *cU1,*cU2,*cORF;
  /* for QMMM, atomnumber contains atomic number of the atoms */
  gmx_bool          *bQM;
  /* The range of home atoms */
  int           start;
  int           homenr;
  /* The lambda value used to create the contents of the struct */
  real          lambda;
} t_mdatoms;

#ifdef __cplusplus
}
#endif


#endif
