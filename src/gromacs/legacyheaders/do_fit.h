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

#ifndef _do_fit_h
#define _do_fit_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

real calc_similar_ind(gmx_bool bRho,int nind,atom_id *index,real mass[],
			     rvec x[],rvec xp[]);
/* Returns RMSD or Rho (depending on bRho) over all atoms in index */

real rmsdev_ind(int nind,atom_id index[],real mass[],
		       rvec x[],rvec xp[]);
/* Returns the RMS Deviation betweem x and xp over all atoms in index */

real rmsdev(int natoms,real mass[],rvec x[],rvec xp[]);
/* Returns the RMS Deviation betweem x and xp over all atoms */

real rhodev_ind(int nind,atom_id index[],real mass[],rvec x[],rvec xp[]);
/* Returns size-independent Rho similarity parameter over all atoms in index
 * Maiorov & Crippen, PROTEINS 22, 273 (1995).
 */
 
real rhodev(int natoms,real mass[],rvec x[],rvec xp[]);
/* Returns size-independent Rho similarity parameter over all atoms
 * Maiorov & Crippen, PROTEINS 22, 273 (1995).
 */

void calc_fit_R(int ndim,int natoms,real *w_rls,rvec *xp,rvec *x,
		       matrix R);
/* Calculates the rotation matrix R for which
 * sum_i w_rls_i (xp_i - R x_i).(xp_i - R x_i)
 * is minimal. ndim=3 gives full fit, ndim=2 gives xy fit.
 * This matrix is also used do_fit.
 * x_rotated[i] = sum R[i][j]*x[j]
 */

void do_fit_ndim(int ndim,int natoms,real *w_rls,rvec *xp,rvec *x);
/* Do a least squares fit of x to xp. Atoms which have zero mass
 * (w_rls[i]) are not taken into account in fitting.
 * This makes is possible to fit eg. on Calpha atoms and orient
 * all atoms. The routine only fits the rotational part,
 * therefore both xp and x should be centered round the origin.
 */

void do_fit(int natoms,real *w_rls,rvec *xp,rvec *x);
/* Calls do_fit with ndim=3, thus fitting in 3D */

void reset_x_ndim(int ndim,int ncm,const atom_id *ind_cm,
			 int nreset,const atom_id *ind_reset,
			 rvec x[],const real mass[]);
/* Put the center of mass of atoms in the origin for dimensions 0 to ndim.
 * The center of mass is computed from the index ind_cm.
 * When ind_cm!=NULL the COM is determined using ind_cm.
 * When ind_cm==NULL the COM is determined for atoms 0 to ncm.
 * When ind_reset!=NULL the coordinates indexed by ind_reset are reset.
 * When ind_reset==NULL the coordinates up to nreset are reset.
 */

void reset_x(int ncm,const atom_id *ind_cm,
		    int nreset,const atom_id *ind_reset,
		    rvec x[],const real mass[]);
/* Calls reset_x with ndim=3, thus resetting all dimesions */

#ifdef __cplusplus
}
#endif

#endif	/* _do_fit_h */
