/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "gromacs/compat/make_unique.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

const char *ptype_str[eptNR+1] = {
    "Atom", "Nucleus", "Shell", "Bond", "VSite", nullptr
};

void init_atomtypes(t_atomtypes *at)
{
    at->nr         = 0;
    at->atomnumber = nullptr;
}

void done_atomtypes(t_atomtypes *atype)
{
    atype->nr = 0;
    sfree(atype->atomnumber);
}

static void printAtomData(FILE *fp,
                          int indent,
                          const char *title,
                          gmx::ArrayRef<const AtomInfo> atoms)
{
    if (available(fp, &atoms, indent, title))
    {
        indent = pr_title_n(fp, indent, title, atoms.size());
        int i = 0;
        for (const auto &a : atoms)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%6d]={type=%3hu, typeB=%3hu, ptype=%8s, m=%12.5e, "
                    "q=%12.5e, mB=%12.5e, qB=%12.5e, resind=%5d, atomnumber=%3d}\n",
                    title, i, a.type_, a.typeB_,
                    ptype_str[a.ptype_],
                    a.m_, a.q_, a.mB_, a.qB_,
                    a.resind_, a.atomnumber_);
            i++;
        }
    }
}

static void printAtomInfoTypeStrings(FILE *fp,
                                     int indent,
                                     const char *title,
                                     gmx::ArrayRef<const AtomInfo> atoms,
                                     gmx_bool bShowNumbers)
{
    if (available(fp, &atoms, indent, title))
    {
        indent = pr_title_n(fp, indent, title, atoms.size());
        int i = 0;
        for (const auto &a : atoms)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\"}\n",
                    title, bShowNumbers ? i : -1, *(a.atomname));
            i++;
        }
    }
}

static void printAtomNames(FILE *fp,
                           int indent,
                           const char *title,
                           gmx::ArrayRef<const AtomInfo> atoms,
                           gmx_bool bShowNumbers)
{
    if (available(fp, &atoms, indent, title))
    {
        indent = pr_title_n(fp, indent, title, atoms.size());
        int i = 0; 
        for (const auto &a : atoms)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\",nameB=\"%s\"}\n",
                    title, bShowNumbers ? i : -1, *(a.atomtype), *(a.atomtypeB));
            i++;
        }
    }
}

void printResidues(FILE *fp,
                   int indent,
                   const char *title,
                   gmx::ArrayRef<const Residue> resinfo,
                   gmx_bool bShowNumbers)
{
    if (available(fp, &resinfo, indent, title))
    {
        indent = pr_title_n(fp, indent, title, resinfo.size());
        int i = 0;
        for (const auto &r : resinfo)
        {
            pr_indent(fp, indent);
            fprintf(fp, "%s[%d]={name=\"%s\", nr=%d, ic='%c'}\n",
                    title, bShowNumbers ? i : -1,
                    *(r.name_), r.nr_,
                    (r.ic_ == '\0') ? ' ' : r.ic_);
            i++;
        }
    }
}

void printAtoms(FILE *fp,
                int indent,
                const char *title,
                gmx::ArrayRef<const AtomInfo> atoms,
                gmx_bool bShownumbers)
{
    if (available(fp, &atoms, indent, title))
    {
        indent = pr_title(fp, indent, title);
        printAtomData(fp, indent, "atom", atoms);
        printAtomNames(fp, indent, "atom", atoms, bShownumbers);
        printAtomInfoTypeStrings(fp, indent, "type", atoms, bShownumbers);
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
                    "atomtype[%3d]={atomnumber=%4d}\n",
                    bShowNumbers ? i : -1, atomtypes->atomnumber[i]);
        }
    }
}

static void compareAtomData(FILE *fp,
                            int index,
                            const AtomInfo &a1,
                            const AtomInfo &a2,
                            real ftol,
                            real abstol)
{
    cmp_us(fp, "atom.type", index, a1.type_, a2.type_);
    cmp_us(fp, "atom.ptype", index, a1.ptype_, a2.ptype_);
    cmp_int(fp, "atom.resind", index, a1.resind_, a2.resind_);
    cmp_int(fp, "atom.atomnumber", index, a1.atomnumber_, a2.atomnumber_);
    cmp_real(fp, "atom.m", index, a1.m_, a2.m_, ftol, abstol);
    cmp_real(fp, "atom.q", index, a1.q_, a2.q_, ftol, abstol);
    cmp_us(fp, "atom.typeB", index, a1.typeB_, a2.typeB_);
    cmp_real(fp, "atom.mB", index, a1.mB_, a2.mB_, ftol, abstol);
    cmp_real(fp, "atom.qB", index, a1.qB_, a2.qB_, ftol, abstol);
}

void compareAtomFEPData(FILE *fp,
                        gmx::ArrayRef<const AtomInfo> atoms,
                        real ftol,
                        real abstol)
{
    fprintf(fp, "comparing atoms\n");
    int i = 0;
    for (const auto &a : atoms)
    {
        cmp_us(fp, "atom.type", i, a.type_, a.typeB_);
        cmp_real(fp, "atom.m", i, a.m_, a.mB_, ftol, abstol);
        cmp_real(fp, "atom.q", i, a.q_, a.qB_, ftol, abstol);
        i++;
    }
}


void compareAtomInfo(FILE *fp,
                     gmx::ArrayRef<const AtomInfo> a1,
                     gmx::ArrayRef<const AtomInfo> a2,
                     real ftol,
                     real abstol)
{
    fprintf(fp, "comparing atoms\n");

    cmp_int(fp, "atoms.size()", -1, a1.size(), a2.size());
    for (int i = 0; i < std::min(a1.size(), a2.size()); i++)
    {
        compareAtomData(fp, i, a1[i], a2[i], ftol, abstol);
    }
}

void atomsSetMassesBasedOnNames(gmx::ArrayRef<AtomInfo> atoms,
                                gmx::ArrayRef<const Residue> resinfo,
                                bool printMissingMasses)
{
    if (allAtomsHaveMass(atoms))
    {
        /* We could decide to anyhow assign then or generate a fatal error,
         * but it's probably most useful to keep the masses we have.
         */
        return;
    }

    int               maxWarn  = (printMissingMasses ? 10 : 0);
    int               numWarn  = 0;

    AtomPropertiesPtr aps      = gmx::compat::make_unique<AtomProperties>();

    for (int i = 0; i < atoms.size(); i++)
    {
        if (!setAtomProperty(aps.get(), epropMass,
                             *resinfo[atoms[i].resind_].name_,
                             *atoms[i].atomname,
                             &atoms[i].m_))
        {
            atoms[i].haveMass_ = false;

            if (numWarn < maxWarn)
            {
                fprintf(stderr, "Can not find mass in database for atom %s in residue %d %s\n",
                        *atoms[i].atomname,
                        resinfo[atoms[i].resind_].nr_,
                        *resinfo[atoms[i].resind_].name_);
                numWarn++;
            }
            else
            {
                break;
            }
        }
    }
}
