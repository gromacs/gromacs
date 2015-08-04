/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015, by the GROMACS development team, led by
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

#include "gpp_bond_atomtype.h"

#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int              nr;       /* The number of atomtypes		*/
    char          ***atomname; /* Names of the atomtypes		*/
} gpp_bond_atomtype;

int get_bond_atomtype_type(char *str, t_bond_atomtype at)
{
    gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;

    int                i;

    for (i = 0; (i < ga->nr); i++)
    {
        /* Atom types are always case sensitive */
        if (strcmp(str, *(ga->atomname[i])) == 0)
        {
            return i;
        }
    }

    return NOTSET;
}

char *get_bond_atomtype_name(int nt, t_bond_atomtype at)
{
    gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;

    if ((nt < 0) || (nt >= ga->nr))
    {
        return NULL;
    }

    return *(ga->atomname[nt]);
}

t_bond_atomtype init_bond_atomtype(void)
{
    gpp_bond_atomtype *ga;

    snew(ga, 1);

    return (t_bond_atomtype ) ga;
}

void add_bond_atomtype(t_bond_atomtype at, t_symtab *tab,
                       char *name)
{
    gpp_bond_atomtype *ga = (gpp_bond_atomtype *) at;

    ga->nr++;
    srenew(ga->atomname, ga->nr);
    ga->atomname[ga->nr-1] = put_symtab(tab, name);
}

void done_bond_atomtype(t_bond_atomtype *at)
{
    gpp_bond_atomtype *ga = (gpp_bond_atomtype *) *at;

    sfree(ga->atomname);
    ga->nr = 0;
    sfree(ga);

    *at = NULL;
}
