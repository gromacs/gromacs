/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * GRowing Old MAkes el Chrono Sweat
 */
#ifndef _princ_h
#define _princ_h

#include "typedefs.h"

extern void pr_jacobi(real **a,int n,real d[],real **v,int *nrot);

extern void rotate_atoms(int gnx,atom_id index[],rvec x[],matrix trans);
/* Rotate all atoms in index using matrix trans */

extern void principal_comp(int n,atom_id index[],t_atom atom[],rvec x[],
			   matrix trans,rvec d);
/* Calculate the principal components of atoms in index. Atoms are
 * mass weighted. It is assumed that the center of mass is in the origin!
 */
			   
extern real calc_xcm(rvec x[],int gnx,atom_id index[],t_atom atom[],rvec xcm,
		     bool bQ);
/* Calculate the center of mass of the atoms in index. if bQ then the atoms
 * will be charge weighted rather than mass weighted.
 * Returns the total mass/charge.
 */
 
extern real sub_xcm(rvec x[],int gnx,atom_id index[],t_atom atom[],rvec xcm,
		    bool bQ);
/* Calc. the center of mass and subtract it from all coordinates.
 * Returns the original center of mass in xcm
 * Returns the total mass
 */
 
extern void add_xcm(rvec x[],int gnx,atom_id index[],rvec xcm);
/* Increment all atoms in index with xcm */

#endif
