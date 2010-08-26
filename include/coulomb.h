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

#ifndef _coulomb_h
#define _coulomb_h

#include <stdio.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Ewald related stuff */

void 
init_ewald_tab(ewald_tab_t *et, const t_commrec *cr, const t_inputrec *ir,
                   FILE *fp);
/* initialize the ewald table (as found in the t_forcerec) */

real 
calc_ewaldcoeff(real rc,real dtol);
/* Determines the Ewald parameter, both for Ewald and PME */


real
do_ewald(FILE *log,       gmx_bool bVerbose,
         t_inputrec *ir,
         rvec x[],        rvec f[],
         real chargeA[],  real chargeB[],
         rvec box,
         t_commrec *cr,  int natoms,
         matrix lrvir,   real ewaldcoeff,
         real lambda,    real *dvdlambda,
         ewald_tab_t et);
/* Do an Ewald calculation for the long range electrostatics. */
 
real
ewald_LRcorrection(FILE *fp,
			       int start,int end,
			       t_commrec *cr,t_forcerec *fr,
			       real *chargeA,real *chargeB,
			       t_blocka *excl,rvec x[],
			       matrix box,rvec mu_tot[],
			       int ewald_geometry,real epsilon_surface,
			       real lambda,real *dvdlambda,
			       real *vdip,real *vcharge);
/* Calculate the Long range correction to ewald, due to 
 * 1-4 interactions, surface dipole term and charge terms
 */

/* Routines to set global constants for speeding up the calculation
 * of potentials and forces.
 */
void 
set_shift_consts(FILE *log,real r1,real rc,rvec box,
			     t_forcerec *fr);

real
shift_LRcorrection(FILE *fp,int start,int natoms,
			       t_commrec *cr,t_forcerec *fr,
			       real charge[],t_blocka *excl,rvec x[],
			       gmx_bool bOld,matrix box,matrix lrvir);
/* Calculate the self energy and forces
 * when using long range electrostatics methods.
 * Part of this is a constant, it is computed only once and stored in
 * a local variable. The remainder is computed every step.
 * PBC is taken into account. (Erik L.) 
 */

#ifdef __cplusplus
}
#endif
 
#endif


