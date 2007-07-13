/*
 * $Id$
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

#ifndef _pme_h
#define _pme_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "gmxcomplex.h"
#include "fftgrid.h"
#include "gmx_wallcycle.h"

typedef real *splinevec[DIM];

enum { GMX_SUM_QGRID_FORWARD, GMX_SUM_QGRID_BACKWARD };

extern int gmx_pme_init(gmx_pme_t *pmedata,t_commrec *cr,
			t_inputrec *ir,int homenr,
			bool bFreeEnergy);
			
extern int gmx_pme_destroy(FILE *log,gmx_pme_t *pmedata);
/* Initialize and destroy the pme data structures resepectively.
 * Return value 0 indicates all well, non zero is an error code.
 */

extern int gmx_pme_do(gmx_pme_t pme,
		      int start,       int homenr,
		      rvec x[],        rvec f[],
		      real chargeA[],  real chargeB[],
		      matrix box,      t_commrec *cr,
		      t_nrnb *nrnb,
		      matrix lrvir,    real ewaldcoeff,
		      real *energy,    real lambda,    
		      real *dvdlambda, bool bGatherOnly);
/* Do a PME calculation for the long range electrostatics. 
 * If bGatherOnly is set, the energy from the last computation will be used, 
 * and the forces will be interpolated at the new positions. No new solving 
 * is done then.
 * Return value 0 indicates all well, non zero is an error code.
 */

extern int gmx_pmeonly(gmx_pme_t pme,
                       t_commrec *cr,     t_nrnb *mynrnb,
		       gmx_wallcycle_t wcycle,
		       real ewaldcoeff,   bool bGatherOnly);
/* Called on the nodes that do PME exclusively (as slaves) 
 */

extern void gmx_sum_qgrid(gmx_pme_t pme,t_commrec *cr,t_fftgrid *grid,
			  int direction);

/* The following two routines are for PME/PP node splitting: */
extern void gmx_pme_send_x_q(t_commrec *cr,
			     matrix box, rvec *x, 
			     real *chargeA, real *chargeB, 
			     bool bFreeEnergy, real lambda,
			     bool bLastStep);
/* Send the particle coordinates and/or charges from the PP to the PME nodes
 */

extern void gmx_pme_receive_f(t_commrec *cr,
			      rvec f[], matrix vir, 
			      real *energy, real *dvdlambda,
			      float *pme_cycles);
/* PP nodes receive the long range forces from the PME nodes */

#endif
