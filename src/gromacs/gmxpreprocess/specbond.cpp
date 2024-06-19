/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "specbond.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

struct SpecialBond
{
    std::string firstResidue, secondResidue;
    std::string firstAtomName, secondAtomName;
    std::string newFirstResidue, newSecondResidue;
    real        length;
};

static bool yesno()
{
    char c;

    do
    {
        c = toupper(fgetc(stdin));
    } while ((c != 'Y') && (c != 'N'));

    return (c == 'Y');
}

std::vector<SpecialBond> generateSpecialBonds()
{
    const char* sbfile = "specbond.dat";

    std::vector<SpecialBond> specialBonds;
    char                     r1buf[32], r2buf[32], a1buf[32], a2buf[32], nr1buf[32], nr2buf[32];
    double                   length;
    int                      nb1, nb2;
    char**                   lines;

    int nlines = get_lines(sbfile, &lines);
    for (int i = 0; (i < nlines); i++)
    {
        if (sscanf(lines[i], "%s%s%d%s%s%d%lf%s%s", r1buf, a1buf, &nb1, r2buf, a2buf, &nb2, &length, nr1buf, nr2buf)
            != 9)
        {
            fprintf(stderr, "Invalid line '%s' in %s\n", lines[i], sbfile);
        }
        else
        {
            SpecialBond newBond;

            newBond.firstResidue     = r1buf;
            newBond.secondResidue    = r2buf;
            newBond.newFirstResidue  = nr1buf;
            newBond.newSecondResidue = nr2buf;
            newBond.firstAtomName    = a1buf;
            newBond.secondAtomName   = a2buf;
            newBond.length           = length;
            specialBonds.push_back(newBond);
        }
        sfree(lines[i]);
    }
    if (nlines > 0)
    {
        sfree(lines);
    }
    fprintf(stderr, "%zu out of %d lines of %s converted successfully\n", specialBonds.size(), nlines, sbfile);

    return specialBonds;
}

static bool is_special(gmx::ArrayRef<const SpecialBond> sb, const char* res, const char* atom)
{
    return std::any_of(sb.begin(), sb.end(), [res, atom](const auto& bond) {
        return (((strncmp(bond.firstResidue.c_str(), res, 3) == 0)
                 && (gmx::equalCaseInsensitive(bond.firstAtomName, atom)))
                || ((strncmp(bond.secondResidue.c_str(), res, 3) == 0)
                    && (gmx::equalCaseInsensitive(bond.secondAtomName, atom))));
    });
}

static bool is_bond(gmx::ArrayRef<const SpecialBond> sb, t_atoms* pdba, int a1, int a2, real d, int* index_sb, bool* bSwap)
{
    const char* at1  = *pdba->atomname[a1];
    const char* at2  = *pdba->atomname[a2];
    const char* res1 = *pdba->resinfo[pdba->atom[a1].resind].name;
    const char* res2 = *pdba->resinfo[pdba->atom[a2].resind].name;

    int i = 0;
    for (const auto& bond : sb)
    {
        *index_sb = i;
        if (((strncmp(bond.firstResidue.c_str(), res1, 3) == 0)
             && (gmx::equalCaseInsensitive(bond.firstAtomName, at1))
             && (strncmp(bond.secondResidue.c_str(), res2, 3) == 0)
             && (gmx::equalCaseInsensitive(bond.secondAtomName, at2))))
        {
            *bSwap = FALSE;
            if ((0.9 * bond.length < d) && (1.1 * bond.length > d))
            {
                return TRUE;
            }
        }
        if (((strncmp(bond.firstResidue.c_str(), res2, 3) == 0)
             && (gmx::equalCaseInsensitive(bond.firstAtomName, at2))
             && (strncmp(bond.secondResidue.c_str(), res1, 3) == 0)
             && (gmx::equalCaseInsensitive(bond.secondAtomName, at1))))
        {
            *bSwap = TRUE;
            if ((0.9 * bond.length < d) && (1.1 * bond.length > d))
            {
                return TRUE;
            }
        }
        i++;
    }
    return FALSE;
}

static void rename_1res(t_atoms* pdba, t_symtab* symtab, int resind, const char* newres, bool bVerbose)
{
    if (bVerbose)
    {
        printf("Using rtp entry %s for %s %d\n",
               newres,
               *pdba->resinfo[resind].name,
               pdba->resinfo[resind].nr);
    }
    pdba->resinfo[resind].rtp = put_symtab(symtab, newres);
}

std::vector<DisulfideBond>
makeDisulfideBonds(t_atoms* pdba, t_symtab* symtab, rvec x[], bool bInteractive, bool bVerbose)
{
    std::vector<SpecialBond>   specialBonds = generateSpecialBonds();
    std::vector<DisulfideBond> bonds;
    bool                       bSwap;
    int                        index_sb;
    char                       buf[10];


    if (!specialBonds.empty())
    {
        std::vector<int> specialBondResIdxs;
        std::vector<int> specialBondAtomIdxs;

        for (int i = 0; (i < pdba->nr); i++)
        {
            /* Check if this atom is special and if it is not a double atom
             * in the input that still needs to be removed.
             */
            int prevAtom = -1;
            if (!specialBondAtomIdxs.empty())
            {
                prevAtom = specialBondAtomIdxs.back();
            }
            if (is_special(specialBonds, *pdba->resinfo[pdba->atom[i].resind].name, *pdba->atomname[i])
                && !(!specialBondAtomIdxs.empty() && pdba->atom[prevAtom].resind == pdba->atom[i].resind
                     && gmx_strcasecmp(*pdba->atomname[prevAtom], *pdba->atomname[i]) == 0))
            {
                specialBondResIdxs.push_back(pdba->atom[i].resind);
                specialBondAtomIdxs.push_back(i);
            }
        }
        /* distance matrix d[nspec][nspec] */
        int                            nspec = specialBondAtomIdxs.size();
        std::vector<std::vector<real>> d(nspec);
        for (int i = 0; (i < nspec); i++)
        {
            d[i].resize(nspec);
            int ai = specialBondAtomIdxs[i];
            for (int j = 0; (j < nspec); j++)
            {
                int aj  = specialBondAtomIdxs[j];
                d[i][j] = std::sqrt(distance2(x[ai], x[aj]));
            }
        }
        if (nspec > 1)
        {
#define MAXCOL 7
            fprintf(stderr, "Special Atom Distance matrix:\n");
            for (int b = 0; (b < nspec); b += MAXCOL)
            {
                /* print resname/number column headings */
                fprintf(stderr, "%8s%8s", "", "");
                int e = std::min(b + MAXCOL, nspec - 1);
                for (int i = b; (i < e); i++)
                {
                    sprintf(buf,
                            "%s%d",
                            *pdba->resinfo[pdba->atom[specialBondAtomIdxs[i]].resind].name,
                            pdba->resinfo[specialBondResIdxs[i]].nr);
                    fprintf(stderr, "%8s", buf);
                }
                fprintf(stderr, "\n");
                /* print atomname/number column headings */
                fprintf(stderr, "%8s%8s", "", "");
                e = std::min(b + MAXCOL, nspec - 1);
                for (int i = b; (i < e); i++)
                {
                    std::string buf = gmx::formatString(
                            "%s%d", *pdba->atomname[specialBondAtomIdxs[i]], specialBondAtomIdxs[i] + 1);
                    fprintf(stderr, "%8s", buf.c_str());
                }
                fprintf(stderr, "\n");
                /* print matrix */
                e = std::min(b + MAXCOL, nspec);
                for (int i = b + 1; (i < nspec); i++)
                {
                    std::string buf = gmx::formatString(
                            "%s%d",
                            *pdba->resinfo[pdba->atom[specialBondAtomIdxs[i]].resind].name,
                            pdba->resinfo[specialBondResIdxs[i]].nr);
                    fprintf(stderr, "%8s", buf.c_str());
                    buf = gmx::formatString(
                            "%s%d", *pdba->atomname[specialBondAtomIdxs[i]], specialBondAtomIdxs[i] + 1);
                    fprintf(stderr, "%8s", buf.c_str());
                    int e2 = std::min(i, e);
                    for (int j = b; (j < e2); j++)
                    {
                        fprintf(stderr, " %7.3f", d[i][j]);
                    }
                    fprintf(stderr, "\n");
                }
            }
        }

        for (int i = 0; (i < nspec); i++)
        {
            int ai = specialBondAtomIdxs[i];
            for (int j = i + 1; (j < nspec); j++)
            {
                int aj = specialBondAtomIdxs[j];
                /* Ensure creation of at most nspec special bonds to avoid overflowing bonds[] */
                if (bonds.size() < specialBondAtomIdxs.size()
                    && is_bond(specialBonds, pdba, ai, aj, d[i][j], &index_sb, &bSwap))
                {
                    fprintf(stderr,
                            "%s %s-%d %s-%d and %s-%d %s-%d%s",
                            bInteractive ? "Link" : "Linking",
                            *pdba->resinfo[pdba->atom[ai].resind].name,
                            pdba->resinfo[specialBondResIdxs[i]].nr,
                            *pdba->atomname[ai],
                            ai + 1,
                            *pdba->resinfo[pdba->atom[aj].resind].name,
                            pdba->resinfo[specialBondResIdxs[j]].nr,
                            *pdba->atomname[aj],
                            aj + 1,
                            bInteractive ? " (y/n) ?" : "...\n");
                    bool bDoit = bInteractive ? yesno() : true;

                    if (bDoit)
                    {
                        DisulfideBond newBond;
                        /* Store the residue numbers in the bonds array */
                        newBond.firstResidue  = specialBondResIdxs[i];
                        newBond.secondResidue = specialBondResIdxs[j];
                        newBond.firstAtom     = *pdba->atomname[ai];
                        newBond.secondAtom    = *pdba->atomname[aj];
                        bonds.push_back(newBond);
                        /* rename residues */
                        if (bSwap)
                        {
                            rename_1res(pdba,
                                        symtab,
                                        specialBondResIdxs[i],
                                        specialBonds[index_sb].newSecondResidue.c_str(),
                                        bVerbose);
                            rename_1res(pdba,
                                        symtab,
                                        specialBondResIdxs[j],
                                        specialBonds[index_sb].newFirstResidue.c_str(),
                                        bVerbose);
                        }
                        else
                        {
                            rename_1res(pdba,
                                        symtab,
                                        specialBondResIdxs[i],
                                        specialBonds[index_sb].newFirstResidue.c_str(),
                                        bVerbose);
                            rename_1res(pdba,
                                        symtab,
                                        specialBondResIdxs[j],
                                        specialBonds[index_sb].newSecondResidue.c_str(),
                                        bVerbose);
                        }
                    }
                }
            }
        }
    }

    return bonds;
}
