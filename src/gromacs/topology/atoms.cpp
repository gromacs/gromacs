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
#include "gmxpre.h"

#include "atoms.h"

#include <cstring>

#include "gromacs/topology/symtab.h"
#include "gromacs/utility/smalloc.h"

void init_atom(t_atoms *at)
{
    at->nr        = 0;
    at->nres      = 0;
    at->atom      = NULL;
    at->resinfo   = NULL;
    at->atomname  = NULL;
    at->atomtype  = NULL;
    at->atomtypeB = NULL;
    at->pdbinfo   = NULL;
}

void init_atomtypes(t_atomtypes *at)
{
    at->nr         = 0;
    at->radius     = NULL;
    at->vol        = NULL;
    at->atomnumber = NULL;
    at->gb_radius  = NULL;
    at->S_hct      = NULL;
}

void done_atom(t_atoms *at)
{
    at->nr       = 0;
    at->nres     = 0;
    sfree(at->atom);
    sfree(at->resinfo);
    sfree(at->atomname);
    sfree(at->atomtype);
    sfree(at->atomtypeB);
    if (at->pdbinfo)
    {
        sfree(at->pdbinfo);
    }
}

void done_atomtypes(t_atomtypes *atype)
{
    atype->nr = 0;
    sfree(atype->radius);
    sfree(atype->vol);
    sfree(atype->surftens);
    sfree(atype->atomnumber);
    sfree(atype->gb_radius);
    sfree(atype->S_hct);
}

void add_t_atoms(t_atoms *atoms, int natom_extra, int nres_extra)
{
    int i;

    if (natom_extra > 0)
    {
        srenew(atoms->atomname, atoms->nr+natom_extra);
        srenew(atoms->atom, atoms->nr+natom_extra);
        if (NULL != atoms->pdbinfo)
        {
            srenew(atoms->pdbinfo, atoms->nr+natom_extra);
        }
        if (NULL != atoms->atomtype)
        {
            srenew(atoms->atomtype, atoms->nr+natom_extra);
        }
        if (NULL != atoms->atomtypeB)
        {
            srenew(atoms->atomtypeB, atoms->nr+natom_extra);
        }
        for (i = atoms->nr; (i < atoms->nr+natom_extra); i++)
        {
            atoms->atomname[i] = NULL;
            memset(&atoms->atom[i], 0, sizeof(atoms->atom[i]));
            if (NULL != atoms->pdbinfo)
            {
                std::memset(&atoms->pdbinfo[i], 0, sizeof(atoms->pdbinfo[i]));
            }
            if (NULL != atoms->atomtype)
            {
                atoms->atomtype[i] = NULL;
            }
            if (NULL != atoms->atomtypeB)
            {
                atoms->atomtypeB[i] = NULL;
            }
        }
        atoms->nr += natom_extra;
    }
    if (nres_extra > 0)
    {
        srenew(atoms->resinfo, atoms->nres+nres_extra);
        for (i = atoms->nres; (i < atoms->nres+nres_extra); i++)
        {
            std::memset(&atoms->resinfo[i], 0, sizeof(atoms->resinfo[i]));
        }
        atoms->nres += nres_extra;
    }
}

void init_t_atoms(t_atoms *atoms, int natoms, gmx_bool bPdbinfo)
{
    atoms->nr   = natoms;
    atoms->nres = 0;
    snew(atoms->atomname, natoms);
    atoms->atomtype  = NULL;
    atoms->atomtypeB = NULL;
    snew(atoms->resinfo, natoms);
    snew(atoms->atom, natoms);
    if (bPdbinfo)
    {
        snew(atoms->pdbinfo, natoms);
    }
    else
    {
        atoms->pdbinfo = NULL;
    }
}

t_atoms *copy_t_atoms(t_atoms *src)
{
    t_atoms *dst;
    int      i;

    snew(dst, 1);
    init_t_atoms(dst, src->nr, (NULL != src->pdbinfo));
    dst->nr = src->nr;
    if (NULL != src->atomname)
    {
        snew(dst->atomname, src->nr);
    }
    if (NULL != src->atomtype)
    {
        snew(dst->atomtype, src->nr);
    }
    if (NULL != src->atomtypeB)
    {
        snew(dst->atomtypeB, src->nr);
    }
    for (i = 0; (i < src->nr); i++)
    {
        dst->atom[i] = src->atom[i];
        if (NULL != src->pdbinfo)
        {
            dst->pdbinfo[i] = src->pdbinfo[i];
        }
        if (NULL != src->atomname)
        {
            dst->atomname[i]  = src->atomname[i];
        }
        if (NULL != src->atomtype)
        {
            dst->atomtype[i] = src->atomtype[i];
        }
        if (NULL != src->atomtypeB)
        {
            dst->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    dst->nres = src->nres;
    for (i = 0; (i < src->nres); i++)
    {
        dst->resinfo[i] = src->resinfo[i];
    }
    return dst;
}

void t_atoms_set_resinfo(t_atoms *atoms, int atom_ind, t_symtab *symtab,
                         const char *resname, int resnr, unsigned char ic,
                         int chainnum, char chainid)
{
    t_resinfo *ri;

    ri           = &atoms->resinfo[atoms->atom[atom_ind].resind];
    ri->name     = put_symtab(symtab, resname);
    ri->rtp      = NULL;
    ri->nr       = resnr;
    ri->ic       = ic;
    ri->chainnum = chainnum;
    ri->chainid  = chainid;
}

void free_t_atoms(t_atoms *atoms, gmx_bool bFreeNames)
{
    int i;

    if (bFreeNames && atoms->atomname != NULL)
    {
        for (i = 0; i < atoms->nr; i++)
        {
            if (atoms->atomname[i] != NULL)
            {
                sfree(*atoms->atomname[i]);
                *atoms->atomname[i] = NULL;
            }
        }
    }
    if (bFreeNames && atoms->resinfo != NULL)
    {
        for (i = 0; i < atoms->nres; i++)
        {
            if (atoms->resinfo[i].name != NULL)
            {
                sfree(*atoms->resinfo[i].name);
                *atoms->resinfo[i].name = NULL;
            }
        }
    }
    if (bFreeNames && atoms->atomtype != NULL)
    {
        for (i = 0; i < atoms->nr; i++)
        {
            if (atoms->atomtype[i] != NULL)
            {
                sfree(*atoms->atomtype[i]);
                *atoms->atomtype[i] = NULL;
            }
        }
    }
    if (bFreeNames && atoms->atomtypeB != NULL)
    {
        for (i = 0; i < atoms->nr; i++)
        {
            if (atoms->atomtypeB[i] != NULL)
            {
                sfree(*atoms->atomtypeB[i]);
                *atoms->atomtypeB[i] = NULL;
            }
        }
    }
    sfree(atoms->atomname);
    sfree(atoms->atomtype);
    sfree(atoms->atomtypeB);
    sfree(atoms->resinfo);
    sfree(atoms->atom);
    sfree(atoms->pdbinfo);
    atoms->nr        = 0;
    atoms->nres      = 0;
    atoms->atomname  = NULL;
    atoms->atomtype  = NULL;
    atoms->atomtypeB = NULL;
    atoms->resinfo   = NULL;
    atoms->atom      = NULL;
    atoms->pdbinfo   = NULL;
}
