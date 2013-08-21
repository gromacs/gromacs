/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "topdirs.h"

#include <stdarg.h>
#include <stdio.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* Must correspond to the directive enum in grompp-impl.h */
static const char *directive_names[d_maxdir+1] = {
    "defaults",
    "atomtypes",
    "bondtypes",
    "constrainttypes",
    "pairtypes",
    "angletypes",
    "dihedraltypes",
    "nonbond_params",
    "implicit_genborn_params",
    "implicit_surface_params",
    "cmaptypes",
    /* All the directives above can not appear after moleculetype */
    "moleculetype",
    "atoms",
    "virtual_sites2",
    "virtual_sites3",
    "virtual_sites4",
    "virtual_sitesn",
    "bonds",
    "exclusions",
    "pairs",
    "pairs_nb",
    "angles",
    "dihedrals",
    "constraints",
    "settles",
    "polarization",
    "water_polarization",
    "thole_polarization",
    "system",
    "molecules",
    "position_restraints",
    "angle_restraints",
    "angle_restraints_z",
    "distance_restraints",
    "orientation_restraints",
    "dihedral_restraints",
    "cmap",
    "intermolecular_interactions",
    "invalid"
};

int ifunc_index(directive d, int type)
{
    switch (d)
    {
        case d_bondtypes:
        case d_bonds:
            switch (type)
            {
                case 1:
                    return F_BONDS;
                case 2:
                    return F_G96BONDS;
                case 3:
                    return F_MORSE;
                case 4:
                    return F_CUBICBONDS;
                case 5:
                    return F_CONNBONDS;
                case 6:
                    return F_HARMONIC;
                case 7:
                    return F_FENEBONDS;
                case 8:
                    return F_TABBONDS;
                case 9:
                    return F_TABBONDSNC;
                case 10:
                    return F_RESTRBONDS;
                default:
                    gmx_fatal(FARGS, "Invalid bond type %d", type);
                    break;
            }
            break;
        case d_angles:
        case d_angletypes:
            switch (type)
            {
                case 1:
                    return F_ANGLES;
                case 2:
                    return F_G96ANGLES;
                case 3:
                    return F_CROSS_BOND_BONDS;
                case 4:
                    return F_CROSS_BOND_ANGLES;
                case 5:
                    return F_UREY_BRADLEY;
                case 6:
                    return F_QUARTIC_ANGLES;
                case 8:
                    return F_TABANGLES;
                case 9:
                    return F_LINEAR_ANGLES;
                case 10:
                    return F_RESTRANGLES;
                default:
                    gmx_fatal(FARGS, "Invalid angle type %d", type);
                    break;
            }
            break;
        case d_pairs:
        case d_pairtypes:
            if (type == 1 || (d == d_pairtypes && type == 2))
            {
                return F_LJ14;
            }
            else if (type == 2)
            {
                return F_LJC14_Q;
            }
            else
            {
                gmx_fatal(FARGS, "Invalid pairs type %d", type);
            }
            break;
        case d_pairs_nb:
            return F_LJC_PAIRS_NB;
        case d_dihedrals:
        case d_dihedraltypes:
            switch (type)
            {
                case 1:
                    return F_PDIHS;
                case 2:
                    return F_IDIHS;
                case 3:
                    return F_RBDIHS;
                case 4:
                    return F_PIDIHS;
                case 5:
                    return F_FOURDIHS;
                case 8:
                    return F_TABDIHS;
                case 9:
                    return F_PDIHS; /* proper dihedrals where we allow multiple terms over single bond */
                case 10:
                    return F_RESTRDIHS;
                case 11:
                    return F_CBTDIHS;
                default:
                    gmx_fatal(FARGS, "Invalid dihedral type %d", type);
            }
            break;
        case d_cmaptypes:
        case d_cmap:
            return F_CMAP;

        case d_nonbond_params:
            if (type == 1)
            {
                return F_LJ;
            }
            else
            {
                return F_BHAM;
            }
        case d_vsites2:
            return F_VSITE2;
        case d_vsites3:
            switch (type)
            {
                case 1:
                    return F_VSITE3;
                case 2:
                    return F_VSITE3FD;
                case 3:
                    return F_VSITE3FAD;
                case 4:
                    return F_VSITE3OUT;
                default:
                    gmx_fatal(FARGS, "Invalid vsites3 type %d", type);
            }
            break;
        case d_vsites4:
            switch (type)
            {
                case 1:
                    return F_VSITE4FD;
                case 2:
                    return F_VSITE4FDN;
                default:
                    gmx_fatal(FARGS, "Invalid vsites4 type %d", type);
            }
            break;
        case d_vsitesn:
            return F_VSITEN;
        case d_constraints:
        case d_constrainttypes:
            switch (type)
            {
                case 1:
                    return F_CONSTR;
                case 2:
                    return F_CONSTRNC;
                default:
                    gmx_fatal(FARGS, "Invalid constraints type %d", type);
            }
            break;
        case d_settles:
            return F_SETTLE;
        case d_position_restraints:
            switch (type)
            {
                case 1:
                    return F_POSRES;
                case 2:
                    return F_FBPOSRES;
                    break;
                default:
                    gmx_fatal(FARGS, "Invalid position restraint type %d", type);
            }
            break;
        case d_polarization:
            switch (type)
            {
                case 1:
                    return F_POLARIZATION;
                case 2:
                    return F_ANHARM_POL;
                default:
                    gmx_fatal(FARGS, "Invalid polarization type %d", type);
            }
            break;
        case d_thole_polarization:
            return F_THOLE_POL;
        case d_water_polarization:
            return F_WATER_POL;
        case d_angle_restraints:
            return F_ANGRES;
        case d_angle_restraints_z:
            return F_ANGRESZ;
        case d_distance_restraints:
            return F_DISRES;
        case d_orientation_restraints:
            return F_ORIRES;
        case d_dihedral_restraints:
            return F_DIHRES;
        default:
            gmx_fatal(FARGS, "invalid directive %s in ifunc_index (%s:%s)",
                      dir2str(d), __FILE__, __LINE__);
    }
    return -1;
}

const char *dir2str (directive d)
{
    if (d < d_maxdir)
    {
        return directive_names[d];
    }
    else
    {
        return directive_names[d_maxdir];
    }
}

directive str2dir (char *dstr)
{
    int  d;
    char buf[STRLEN], *ptr;

    /* Hack to be able to read old topologies */
    if (gmx_strncasecmp_min(dstr, "dummies", 7) == 0)
    {
        sprintf(buf, "virtual_sites%s", dstr+7);
        ptr = buf;
    }
    else
    {
        ptr = dstr;
    }

    for (d = 0; (d < d_maxdir); d++)
    {
        if (gmx_strcasecmp_min(ptr, dir2str((directive)d)) == 0)
        {
            return (directive)d;
        }
    }

    return d_invalid;
}

static directive **necessary = NULL;

static void set_nec(directive **n, ...)
/* Must always have at least one extra argument */
{
    va_list   ap;
    int       ind = 0;
    directive d;

    va_start(ap, n);
    do
    {
        d = (directive)va_arg(ap, int);
        srenew(*n, ++ind);
        (*n)[ind-1] = d;
    }
    while (d != d_none);
    va_end(ap);
}

void DS_Init(DirStack **DS)
{
    if (necessary == NULL)
    {
        int i;

        snew(necessary, d_maxdir);
        set_nec(&(necessary[d_defaults]), d_none);
        set_nec(&(necessary[d_atomtypes]), d_defaults, d_none);
        set_nec(&(necessary[d_bondtypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_constrainttypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_pairtypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_angletypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_dihedraltypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_nonbond_params]), d_atomtypes, d_none);
        set_nec(&(necessary[d_implicit_genborn_params]), d_atomtypes, d_none);
        set_nec(&(necessary[d_implicit_surface_params]), d_atomtypes, d_none);
        set_nec(&(necessary[d_cmaptypes]), d_atomtypes, d_none);
        set_nec(&(necessary[d_moleculetype]), d_atomtypes, d_none);
        set_nec(&(necessary[d_atoms]), d_moleculetype, d_none);
        set_nec(&(necessary[d_vsites2]), d_atoms, d_none);
        set_nec(&(necessary[d_vsites3]), d_atoms, d_none);
        set_nec(&(necessary[d_vsites4]), d_atoms, d_none);
        set_nec(&(necessary[d_vsitesn]), d_atoms, d_none);
        set_nec(&(necessary[d_bonds]), d_atoms, d_none);
        set_nec(&(necessary[d_exclusions]), d_bonds, d_constraints, d_settles, d_none);
        set_nec(&(necessary[d_pairs]), d_atoms, d_none);
        set_nec(&(necessary[d_pairs_nb]), d_atoms, d_none);
        set_nec(&(necessary[d_angles]), d_atoms, d_none);
        set_nec(&(necessary[d_polarization]), d_atoms, d_none);
        set_nec(&(necessary[d_water_polarization]), d_atoms, d_none);
        set_nec(&(necessary[d_thole_polarization]), d_atoms, d_none);
        set_nec(&(necessary[d_dihedrals]), d_atoms, d_none);
        set_nec(&(necessary[d_constraints]), d_atoms, d_none);
        set_nec(&(necessary[d_settles]), d_atoms, d_none);
        set_nec(&(necessary[d_system]), d_moleculetype, d_none);
        set_nec(&(necessary[d_molecules]), d_system, d_none);
        set_nec(&(necessary[d_position_restraints]), d_atoms, d_none);
        set_nec(&(necessary[d_angle_restraints]), d_atoms, d_none);
        set_nec(&(necessary[d_angle_restraints_z]), d_atoms, d_none);
        set_nec(&(necessary[d_distance_restraints]), d_atoms, d_none);
        set_nec(&(necessary[d_orientation_restraints]), d_atoms, d_none);
        set_nec(&(necessary[d_dihedral_restraints]), d_atoms, d_none);
        set_nec(&(necessary[d_cmap]), d_atoms, d_none);
        set_nec(&(necessary[d_intermolecular_interactions]), d_molecules, d_none);

        for (i = 0; (i < d_maxdir); i++)
        {
            if (debug)
            {
                fprintf(debug, "%20s:  ", dir2str((directive)i));
            }
            if (necessary[i])
            {
                directive d;
                int       j = 0;
                do
                {
                    d = necessary[i][j++];
                    if (debug)
                    {
                        fprintf(debug, "%20s  ", dir2str(d));
                    }
                }
                while (d != d_none);
            }
            if (debug)
            {
                fprintf(debug, "\n");
            }
        }
    }
    *DS = NULL;

}

void DS_Done (DirStack **DS)
{
    DirStack *D;

    while (*DS != NULL)
    {
        D   = *DS;
        *DS = (*DS)->prev;
        sfree (D);
    }
}

void DS_Push (DirStack **DS, directive d)
{
    DirStack *D;

    snew(D, 1);
    D->d    = d;
    D->prev = *DS;
    *DS     = D;
}

int DS_Search(DirStack *DS, directive d)
{
    DirStack *D;

    D = DS;
    while ((D != NULL) && (D->d != d))
    {
        D = D->prev;
    }

    return (D != NULL);
}

int DS_Check_Order(DirStack *DS, directive d)
{
    directive d0;
    int       i = 0;

    /* Check if parameter definitions appear after a moleculetype directive */
    if (d < d_moleculetype && DS_Search(DS, d_moleculetype))
    {
        return FALSE;
    }

    /* Check if all the necessary directives have appeared before directive d */
    if (necessary[d][0] == d_none)
    {
        return TRUE;
    }
    else
    {
        do
        {
            d0 = necessary[d][i++];
            if (DS_Search(DS, d0))
            {
                return TRUE;
            }
        }
        while (d0 != d_none);
    }
    return FALSE;
}
