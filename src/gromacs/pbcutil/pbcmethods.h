/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#ifndef GMX_PBCUTIL_PBCMETHODS_H
#define GMX_PBCUTIL_PBCMETHODS_H

struct t_topology;
struct t_block;
struct t_atom;

#include "gmxpre.h"

#include "gromacs/math/vec.h"

enum
{
    euSel,
    euRect,
    euTric,
    euCompact,
    euNR
};

void calc_pbc_cluster(int ecenter, int nrefat, t_topology* top, int ePBC, rvec x[], const int index[], matrix box);


void put_molecule_com_in_box(int      unitcell_enum,
                             int      ecenter,
                             t_block* mols,
                             int      natoms,
                             t_atom   atom[],
                             int      ePBC,
                             matrix   box,
                             rvec     x[]);

void put_residue_com_in_box(int unitcell_enum, int ecenter, int natoms, t_atom atom[], int ePBC, matrix box, rvec x[]);

void center_x(int ecenter, rvec x[], matrix box, int n, int nc, const int ci[]);

#endif
