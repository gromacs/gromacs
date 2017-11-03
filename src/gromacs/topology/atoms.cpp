/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <cstdio>
#include <cstring>

#include <algorithm>

#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

const char *ptype_str[eptNR+1] = {
    "Atom", "Nucleus", "Shell", "Bond", "VSite", nullptr
};

void init_atom(t_atoms *at)
{
    at->nr          = 0;
    at->nres        = 0;
    at->atom        = nullptr;
    at->resinfo     = nullptr;
    at->atomname    = nullptr;
    at->atomtype    = nullptr;
    at->atomtypeB   = nullptr;
    at->pdbinfo     = nullptr;
    at->haveMass    = FALSE;
    at->haveCharge  = FALSE;
    at->haveType    = FALSE;
    at->haveBState  = FALSE;
    at->havePdbInfo = FALSE;
}

void init_atomtypes(t_atomtypes *at)
{
    at->nr         = 0;
    at->radius     = nullptr;
    at->vol        = nullptr;
    at->atomnumber = nullptr;
    at->gb_radius  = nullptr;
    at->S_hct      = nullptr;
}

void done_atom(t_atoms *at)
{
    sfree(at->atom);
    sfree(at->resinfo);
    sfree(at->atomname);
    sfree(at->atomtype);
    sfree(at->atomtypeB);
    sfree(at->pdbinfo);
    init_atom(at);
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
        if (nullptr != atoms->pdbinfo)
        {
            srenew(atoms->pdbinfo, atoms->nr+natom_extra);
        }
        if (nullptr != atoms->atomtype)
        {
            srenew(atoms->atomtype, atoms->nr+natom_extra);
        }
        if (nullptr != atoms->atomtypeB)
        {
            srenew(atoms->atomtypeB, atoms->nr+natom_extra);
        }
        for (i = atoms->nr; (i < atoms->nr+natom_extra); i++)
        {
            atoms->atomname[i] = nullptr;
            memset(&atoms->atom[i], 0, sizeof(atoms->atom[i]));
            if (nullptr != atoms->pdbinfo)
            {
                std::memset(&atoms->pdbinfo[i], 0, sizeof(atoms->pdbinfo[i]));
            }
            if (nullptr != atoms->atomtype)
            {
                atoms->atomtype[i] = nullptr;
            }
            if (nullptr != atoms->atomtypeB)
            {
                atoms->atomtypeB[i] = nullptr;
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
    atoms->atomtype  = nullptr;
    atoms->atomtypeB = nullptr;
    snew(atoms->resinfo, natoms);
    snew(atoms->atom, natoms);
    atoms->haveMass    = FALSE;
    atoms->haveCharge  = FALSE;
    atoms->haveType    = FALSE;
    atoms->haveBState  = FALSE;
    atoms->havePdbInfo = bPdbinfo;
    if (atoms->havePdbInfo)
    {
        snew(atoms->pdbinfo, natoms);
    }
    else
    {
        atoms->pdbinfo = nullptr;
    }
}

void gmx_pdbinfo_init_default(t_pdbinfo *pdbinfo)
{
    pdbinfo->type         = epdbATOM;
    pdbinfo->atomnr       = 0;
    pdbinfo->altloc       = ' ';
    pdbinfo->atomnm[0]    = '\0';
    pdbinfo->occup        = 1.0;
    pdbinfo->bfac         = 0.0;
    pdbinfo->bAnisotropic = FALSE;
    std::fill(pdbinfo->uij, pdbinfo->uij+6, 0.0);
}

t_atoms *copy_t_atoms(const t_atoms *src)
{
    t_atoms *dst;
    int      i;

    snew(dst, 1);
    init_t_atoms(dst, src->nr, (nullptr != src->pdbinfo));
    dst->nr = src->nr;
    if (nullptr != src->atomname)
    {
        snew(dst->atomname, src->nr);
    }
    if (nullptr != src->atomtype)
    {
        snew(dst->atomtype, src->nr);
    }
    if (nullptr != src->atomtypeB)
    {
        snew(dst->atomtypeB, src->nr);
    }
    for (i = 0; (i < src->nr); i++)
    {
        dst->atom[i] = src->atom[i];
        if (nullptr != src->pdbinfo)
        {
            dst->pdbinfo[i] = src->pdbinfo[i];
        }
        if (nullptr != src->atomname)
        {
            dst->atomname[i]  = src->atomname[i];
        }
        if (nullptr != src->atomtype)
        {
            dst->atomtype[i] = src->atomtype[i];
        }
        if (nullptr != src->atomtypeB)
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
    ri->rtp      = nullptr;
    ri->nr       = resnr;
    ri->ic       = ic;
    ri->chainnum = chainnum;
    ri->chainid  = chainid;
}

static void pr_atom(FILE *fp, int indent, const char *title, const t_atom *atom, int n)
{
    int i;

    if (available(fp, atom, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%6d]={type=%3hu, typeB=%3hu, ptype=%8s, m=%12.5e, "
                    "q=%12.5e, mB=%12.5e, qB=%12.5e, resind=%5d, atomnumber=%3d}\n",
                    title, i, atom[i].type, atom[i].typeB, ptype_str[atom[i].ptype],
                    atom[i].m, atom[i].q, atom[i].mB, atom[i].qB,
                    atom[i].resind, atom[i].atomnumber);
        }
    }
}

static void pr_strings2(FILE *fp, int indent, const char *title,
                        char ***nm, char ***nmB, int n, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, nm, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\",nameB=\"%s\"}\n",
                    title, bShowNumbers ? i : -1, *(nm[i]), *(nmB[i]));
        }
    }
}

static void pr_resinfo(FILE *fp, int indent, const char *title, const t_resinfo *resinfo, int n,
                       gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, resinfo, indent, title))
    {
        indent = pr_title_n(fp, indent, title, n);
        for (i = 0; i < n; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\", nr=%d, ic='%c'}\n",
                    title, bShowNumbers ? i : -1,
                    *(resinfo[i].name), resinfo[i].nr,
                    (resinfo[i].ic == '\0') ? ' ' : resinfo[i].ic);
        }
    }
}

void pr_atoms(FILE *fp, int indent, const char *title, const t_atoms *atoms,
              gmx_bool bShownumbers)
{
    if (available(fp, atoms, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_atom(fp, indent, "atom", atoms->atom, atoms->nr);
        pr_strings(fp, indent, "atom", atoms->atomname, atoms->nr, bShownumbers);
        pr_strings2(fp, indent, "type", atoms->atomtype, atoms->atomtypeB, atoms->nr, bShownumbers);
        pr_resinfo(fp, indent, "residue", atoms->resinfo, atoms->nres, bShownumbers);
    }
}


void pr_atomtypes(FILE *fp, int indent, const char *title, const t_atomtypes *atomtypes,
                  gmx_bool bShowNumbers)
{
    int i;
    if (available(fp, atomtypes, indent, title))
    {
        indent = pr_title(fp, indent, title);
        for (i = 0; i < atomtypes->nr; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp,
                    "atomtype[%3d]={radius=%12.5e, volume=%12.5e, gb_radius=%12.5e, surftens=%12.5e, atomnumber=%4d, S_hct=%12.5e)}\n",
                    bShowNumbers ? i : -1, atomtypes->radius[i], atomtypes->vol[i],
                    atomtypes->gb_radius[i],
                    atomtypes->surftens[i], atomtypes->atomnumber[i], atomtypes->S_hct[i]);
        }
    }
}

static void cmp_atom(FILE *fp, int index, const t_atom *a1, const t_atom *a2, real ftol, real abstol)
{
    if (a2)
    {
        cmp_us(fp, "atom.type", index, a1->type, a2->type);
        cmp_us(fp, "atom.ptype", index, a1->ptype, a2->ptype);
        cmp_int(fp, "atom.resind", index, a1->resind, a2->resind);
        cmp_int(fp, "atom.atomnumber", index, a1->atomnumber, a2->atomnumber);
        cmp_real(fp, "atom.m", index, a1->m, a2->m, ftol, abstol);
        cmp_real(fp, "atom.q", index, a1->q, a2->q, ftol, abstol);
        cmp_us(fp, "atom.typeB", index, a1->typeB, a2->typeB);
        cmp_real(fp, "atom.mB", index, a1->mB, a2->mB, ftol, abstol);
        cmp_real(fp, "atom.qB", index, a1->qB, a2->qB, ftol, abstol);
    }
    else
    {
        cmp_us(fp, "atom.type", index, a1->type, a1->typeB);
        cmp_real(fp, "atom.m", index, a1->m, a1->mB, ftol, abstol);
        cmp_real(fp, "atom.q", index, a1->q, a1->qB, ftol, abstol);
    }
}

void cmp_atoms(FILE *fp, const t_atoms *a1, const t_atoms *a2, real ftol, real abstol)
{
    int i;

    fprintf(fp, "comparing atoms\n");

    if (a2)
    {
        cmp_int(fp, "atoms->nr", -1, a1->nr, a2->nr);
        for (i = 0; i < std::min(a1->nr, a2->nr); i++)
        {
            cmp_atom(fp, i, &(a1->atom[i]), &(a2->atom[i]), ftol, abstol);
        }
    }
    else
    {
        for (i = 0; (i < a1->nr); i++)
        {
            cmp_atom(fp, i, &(a1->atom[i]), nullptr, ftol, abstol);
        }
    }
}

void atomsSetMassesBasedOnNames(t_atoms *atoms, gmx_bool printMissingMasses)
{
    if (atoms->haveMass)
    {
        /* We could decide to anyhow assign then or generate a fatal error,
         * but it's probably most useful to keep the masses we have.
         */
        return;
    }

    int            maxWarn  = (printMissingMasses ? 10 : 0);
    int            numWarn  = 0;

    gmx_atomprop_t aps      = gmx_atomprop_init();

    gmx_bool       haveMass = TRUE;
    for (int i = 0; i < atoms->nr; i++)
    {
        if (!gmx_atomprop_query(aps, epropMass,
                                *atoms->resinfo[atoms->atom[i].resind].name,
                                *atoms->atomname[i],
                                &atoms->atom[i].m))
        {
            haveMass = FALSE;

            if (numWarn < maxWarn)
            {
                fprintf(stderr, "Can not find mass in database for atom %s in residue %d %s\n",
                        *atoms->atomname[i],
                        atoms->resinfo[atoms->atom[i].resind].nr,
                        *atoms->resinfo[atoms->atom[i].resind].name);
                numWarn++;
            }
            else
            {
                break;
            }
        }
    }
    atoms->haveMass = haveMass;

    gmx_atomprop_destroy(aps);
}
