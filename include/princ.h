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

#ifndef _princ_h
#define _princ_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void rotate_atoms(int gnx,atom_id index[],rvec x[],matrix trans);
/* Rotate all atoms in index using matrix trans */

void principal_comp(int n,atom_id index[],t_atom atom[],rvec x[],
			   matrix trans,rvec d);
/* Calculate the principal components of atoms in index. Atoms are
 * mass weighted. It is assumed that the center of mass is in the origin!
 */

void orient_princ(t_atoms *atoms, int isize, atom_id *index,
			 int natoms, rvec x[], rvec *v, rvec d);
/* rotates molecule to align principal axes with coordinate axes */

real calc_xcm(rvec x[],int gnx,atom_id *index,t_atom *atom,rvec xcm,
		     gmx_bool bQ);
/* Calculate the center of mass of the atoms in index. if bQ then the atoms
 * will be charge weighted rather than mass weighted.
 * Returns the total mass/charge.
 */

real sub_xcm(rvec x[],int gnx,atom_id *index,t_atom atom[],rvec xcm,
		    gmx_bool bQ);
/* Calc. the center of mass and subtract it from all coordinates.
 * Returns the original center of mass in xcm
 * Returns the total mass
 */

void add_xcm(rvec x[],int gnx,atom_id *index,rvec xcm);
/* Increment all atoms in index with xcm */

#ifdef __cplusplus
}
#endif

#endif
