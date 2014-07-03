/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_GMXANA_PRINC_H
#define GMX_GMXANA_PRINC_H

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void rotate_atoms(int gnx, atom_id index[], rvec x[], matrix trans);
/* Rotate all atoms in index using matrix trans */

void principal_comp(int n, atom_id index[], t_atom atom[], rvec x[],
                    matrix trans, rvec d);
/* Calculate the principal components of atoms in index. Atoms are
 * mass weighted. It is assumed that the center of mass is in the origin!
 */

void orient_princ(t_atoms *atoms, int isize, atom_id *index,
                  int natoms, rvec x[], rvec *v, rvec d);
/* rotates molecule to align principal axes with coordinate axes */

real calc_xcm(rvec x[], int gnx, atom_id *index, t_atom *atom, rvec xcm,
              gmx_bool bQ);
/* Calculate the center of mass of the atoms in index. if bQ then the atoms
 * will be charge weighted rather than mass weighted.
 * Returns the total mass/charge.
 */

real sub_xcm(rvec x[], int gnx, atom_id *index, t_atom atom[], rvec xcm,
             gmx_bool bQ);
/* Calc. the center of mass and subtract it from all coordinates.
 * Returns the original center of mass in xcm
 * Returns the total mass
 */

void add_xcm(rvec x[], int gnx, atom_id *index, rvec xcm);
/* Increment all atoms in index with xcm */

#ifdef __cplusplus
}
#endif

#endif
