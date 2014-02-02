/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#ifndef GMX_GMXPREPROCESS_ADDCONF_H
#define GMX_GMXPREPROCESS_ADDCONF_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern
void add_conf(t_atoms *atoms, rvec **x, rvec **v, real **r, gmx_bool bSrenew,
              int ePBC, matrix box, gmx_bool bInsert,
              t_atoms *atoms_solvt, rvec *x_solvt, rvec *v_solvt, real *r_solvt,
              gmx_bool bVerbose, real rshell, int max_sol, const output_env_t oenv);
/* Add two conformations together, without generating overlap.
 * When not inserting, don't check overlap in the middle of the box.
 * If rshell > 0, keep all the residues around the protein (0..natoms_prot-1)
 * that are within rshell distance.
 * If max_sol > 0, add max max_sol solvent molecules.
 */

#ifdef __cplusplus
}
#endif

#endif
