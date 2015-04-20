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
#include "gmxpre.h"

#include "read-conformation.h"

#include "gromacs/fileio/confio.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

real *makeExclusionDistances(const t_atoms *a, gmx_atomprop_t aps,
                             real defaultDistance, real scaleFactor)
{
    int   i;
    real *exclusionDistances;

    snew(exclusionDistances, a->nr);
    /* initialise arrays with distances usually based on van der Waals
       radii */
    for (i = 0; (i < a->nr); i++)
    {
        if (!gmx_atomprop_query(aps, epropVDW,
                                *(a->resinfo[a->atom[i].resind].name),
                                *(a->atomname[i]), &(exclusionDistances[i])))
        {
            exclusionDistances[i] = defaultDistance;
        }
        else
        {
            exclusionDistances[i] *= scaleFactor;
        }
    }
    return exclusionDistances;
}

char *readConformation(const char *confin, t_atoms *atoms, rvec **x, rvec **v,
                       int *ePBC, matrix box, const char *statusTitle)
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
    init_t_atoms(atoms, natoms, FALSE);

    /* read residue number, residue names, atomnames, coordinates etc. */
    fprintf(stderr, "Reading %s configuration%s\n", statusTitle, v ? " and velocities" : "");
    read_stx_conf(confin, title, atoms, *x, v ? *v : NULL, ePBC, box);
    fprintf(stderr, "%s\nContaining %d atoms in %d residues\n",
            title, atoms->nr, atoms->nres);

    return title;
}
