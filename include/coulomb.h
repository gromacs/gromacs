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
#include "types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Ewald related stuff */

void 
init_ewald_tab(ewald_tab_t *et, const t_commrec *cr, const t_inputrec *ir,
                   FILE *fp);
/* initialize the ewald table (as found in the t_forcerec) */

GMX_LIBGMX_EXPORT real 
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
 
GMX_LIBGMX_EXPORT real
ewald_LRcorrection(FILE *fp,
		   int start,int end,
		   t_commrec *cr,int thread,t_forcerec *fr,
		   real *chargeA,real *chargeB,
		   gmx_bool calc_excl_corr,
		   t_blocka *excl,rvec x[],
		   matrix box,rvec mu_tot[],
		   int ewald_geometry,real epsilon_surface,
		   rvec *f,tensor vir,
		   real lambda,real *dvdlambda);
/* Calculate the Long range correction to the Ewald sum,
 * due to excluded pairs and/or surface dipole terms.
 */

GMX_LIBGMX_EXPORT real
ewald_charge_correction(t_commrec *cr,t_forcerec *fr,real lambda,matrix box,
			real *dvdlambda,tensor vir);
/* Calculate the Long range correction to the Ewald sum,
 * due to a net system charge.
 * Should only be called on one thread.
 */

/* Routines to set global constants for speeding up the calculation
 * of potentials and forces.
 */
GMX_LIBGMX_EXPORT void 
set_shift_consts(FILE *log,real r1,real rc,rvec box,
			     t_forcerec *fr);

#ifdef __cplusplus
}
#endif
 
#endif


