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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _chargegroup_h
#define _chargegroup_h

#include "sysstuff.h"
#include "typedefs.h"

#ifdef __cplusplus
extern "C" { 
#endif

  void calc_chargegroup_radii(const gmx_mtop_t *mtop,rvec *x,
				     real *rvdw1,real *rvdw2,
				     real *rcoul1,real *rcoul2);
  /* This routine calculates the two largest charge group radii in x,
   * separately for VdW and Coulomb interactions.
   */

  void calc_cgcm(FILE *log,int cg0,int cg1,t_block *cgs,
			rvec pos[],rvec cg_cm[]);
  /* Routine to compute centers of geometry of charge groups. No periodicity
   * is used.
   */
  
  void put_charge_groups_in_box (FILE *log,int cg0,int cg1,
					int ePBC,matrix box,t_block *cgs,
					rvec pos[],
					rvec cg_cm[]);
  /* This routine puts charge groups in the periodic box, keeping them
   * together.
   */
  
#ifdef __cplusplus
}
#endif

#endif	/* _chargegroup_h */
