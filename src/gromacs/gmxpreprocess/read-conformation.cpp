/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include "read-conformation.h"

#include "gromacs/fileio/confio.h"
#include "atomprop.h"
#include "types/simple.h"
#include "types/atoms.h"
#include "smalloc.h"

void mk_vdw(t_atoms *a, real rvdw[], gmx_atomprop_t aps,
            real r_distance, real r_scale)
{
    int i;

    /* initialise van der waals arrays of configuration */
    fprintf(stderr, "Initialising van der waals distances...\n");
    for (i = 0; (i < a->nr); i++)
    {
        if (!gmx_atomprop_query(aps, epropVDW,
                                *(a->resinfo[a->atom[i].resind].name),
                                *(a->atomname[i]), &(rvdw[i])))
        {
            rvdw[i] = r_distance;
        }
        else
        {
            rvdw[i] *= r_scale;
        }
    }
}

char *read_conformation(const char *confin, t_atoms *atoms, rvec **x, rvec **v,
                        real **r, int *ePBC, matrix box, gmx_atomprop_t aps,
                        real r_distance, real r_scale)
{
    char *title;
    int   natoms;

    snew(title, STRLEN);
    get_stx_coordnum(confin, &natoms);

    /* allocate memory for atom coordinates of configuration */
    snew(*x, natoms);
    if (v)
    {
        snew(*v, natoms);
    }
    snew(*r, natoms);
    init_t_atoms(atoms, natoms, FALSE);

    /* read residue number, residue names, atomnames, coordinates etc. */
    fprintf(stderr, "Reading solute configuration%s\n", v ? " and velocities" : "");
    read_stx_conf(confin, title, atoms, *x, v ? *v : NULL, ePBC, box);
    fprintf(stderr, "%s\nContaining %d atoms in %d residues\n",
            title, atoms->nr, atoms->nres);

    mk_vdw(atoms, *r, aps, r_distance, r_scale);

    return title;
}
