/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "hackblock.h"

#include <cstring>

#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/* these MUST correspond to the enum in hackblock.h */
const char* btsNames[ebtsNR] = { "bonds", "angles", "dihedrals", "impropers", "exclusions", "cmap" };
const int   btsNiatoms[ebtsNR] = { 2, 3, 4, 4, 2, 5 };

MoleculePatchType MoleculePatch::type() const
{
    if (oname.empty() && !nname.empty())
    {
        return MoleculePatchType::Add;
    }
    else if (!oname.empty() && nname.empty())
    {
        return MoleculePatchType::Delete;
    }
    else if (!oname.empty() && !nname.empty())
    {
        return MoleculePatchType::Replace;
    }
    else
    {
        GMX_THROW(gmx::InvalidInputError("Unknown type of atom modification"));
    }
}

void clearModificationBlock(MoleculePatchDatabase* globalPatches)
{
    globalPatches->name.clear();
    globalPatches->hack.clear();
    for (int i = 0; i < ebtsNR; i++)
    {
        globalPatches->rb[i].b.clear();
    }
}

#define safe_strdup(str) (((str) != NULL) ? gmx_strdup(str) : NULL)

static bool contains_char(const BondedInteraction& s, char c)
{

    bool bRet = false;
    for (int i = 0; i < MAXATOMLIST; i++)
    {
        if (!s.a[i].empty() && s.a[i][0] == c)
        {
            bRet = true;
        }
    }

    return bRet;
}

static int rbonded_find_atoms_in_list(const BondedInteraction&               b,
                                      gmx::ArrayRef<const BondedInteraction> blist,
                                      int                                    natoms)
{
    int foundPos = -1;

    for (auto it = blist.begin(); (it != blist.end()) && (foundPos < 0); it++)
    {
        bool atomsMatch = true;
        for (int k = 0; k < natoms && atomsMatch; k++)
        {
            atomsMatch = atomsMatch && (b.a[k] == it->a[k]);
        }
        /* Try reverse if forward match did not work */
        if (!atomsMatch)
        {
            atomsMatch = true;
            for (int k = 0; k < natoms && atomsMatch; k++)
            {
                atomsMatch = atomsMatch && (b.a[k] == it->a[natoms - 1 - k]);
            }
        }
        if (atomsMatch)
        {
            foundPos = std::distance(blist.begin(), it);
            /* If all the atoms AND all the parameters match, it is likely that
             * the user made a copy-and-paste mistake (since it would be much cheaper
             * to just bump the force constant 2x if you really want it twice).
             * Since we only have the unparsed string here we can only detect
             * EXACT matches (including identical whitespace).
             */
            if (b.s != it->s)
            {
                gmx_warning("Duplicate line found in or between hackblock and rtp entries");
            }
        }
    }
    return foundPos;
}

bool mergeBondedInteractionList(gmx::ArrayRef<const BondedInteractionList> s,
                                gmx::ArrayRef<BondedInteractionList>       d,
                                bool                                       bMin,
                                bool                                       bPlus)
{
    bool bBondsRemoved = false;
    for (int i = 0; i < ebtsNR; i++)
    {
        if (!s[i].b.empty())
        {
            /* Record how many bonds we have in the destination when we start.
             *
             * If an entry is present in the hackblock (destination), we will
             * not add the one from the main rtp, since the point is for hackblocks
             * to overwrite it. However, if there is no hackblock entry we do
             * allow multiple main rtp entries since some forcefield insist on that.
             *
             * We accomplish this by checking the position we find an entry in,
             * rather than merely checking whether it exists at all.
             * If that index is larger than the original (hackblock) destination
             * size, it was added from the main rtp, and then we will allow more
             * such entries. In contrast, if the entry found has a lower index
             * it is a hackblock entry meant to override the main rtp, and then
             * we don't add the main rtp one.
             */
            int nbHackblockStart = d[i].b.size();

            for (const auto& b : s[i].b)
            {
                /* Check if this bonded string already exists before adding.
                 * We are merging from the main RTP to the hackblocks, so this
                 * will mean the hackblocks overwrite the man RTP, as intended.
                 */
                int index = rbonded_find_atoms_in_list(b, d[i].b, btsNiatoms[i]);
                /* - If we did not find this interaction at all, the index will be -1,
                 *   and then we should definitely add it to the merged hackblock and rtp.
                 *
                 * Alternatively, if it was found, index will be >=0.
                 * - In case this index is lower than the original number of entries,
                 *   it is already present as a *hackblock* entry, and those should
                 *   always override whatever we have listed in the RTP. Thus, we
                 *   should just keep that one and not add anything from the RTP.
                 * - Finally, if it was found, but with an index higher than
                 *   the original number of entries, it comes from the RTP rather
                 *   than hackblock, and then we must have added it ourselves
                 *   in a previous iteration. In that case it is a matter of
                 *   several entries for the same sequence of atoms, and we allow
                 *   that in the RTP. In this case we should simply copy all of
                 *   them, including this one.
                 */
                if (index < 0 || index >= nbHackblockStart)
                {
                    if (!(bMin && contains_char(b, '-')) && !(bPlus && contains_char(b, '+')))
                    {
                        d[i].b.push_back(b);
                    }
                    else if (i == ebtsBONDS)
                    {
                        bBondsRemoved = true;
                    }
                }
                else
                {
                    /* This is the common case where a hackblock entry simply
                     * overrides the RTP, so we cannot warn here.
                     */
                }
            }
        }
    }
    return bBondsRemoved;
}

void copyPreprocessResidues(const PreprocessResidue& s, PreprocessResidue* d, t_symtab* symtab)
{
    *d = s;
    d->atom.clear();
    for (const auto& a : s.atom)
    {
        d->atom.push_back(a);
    }
    d->atomname.clear();
    for (const auto& a : s.atomname)
    {
        d->atomname.push_back(put_symtab(symtab, *a));
    }
    d->cgnr.clear();
    for (const auto& c : s.cgnr)
    {
        d->cgnr.push_back(c);
    }
    for (int i = 0; i < ebtsNR; i++)
    {
        d->rb[i].type = s.rb[i].type;
        d->rb[i].b.clear();
    }
    mergeBondedInteractionList(s.rb, d->rb, FALSE, FALSE);
}

void mergeAtomModifications(const MoleculePatchDatabase& s, MoleculePatchDatabase* d)
{
    for (const auto& h : s.hack)
    {
        d->hack.push_back(h);
    }
}

void mergeAtomAndBondModifications(const MoleculePatchDatabase& s, MoleculePatchDatabase* d)
{
    mergeAtomModifications(s, d);
    mergeBondedInteractionList(s.rb, d->rb, FALSE, FALSE);
}

void copyModificationBlocks(const MoleculePatchDatabase& s, MoleculePatchDatabase* d)
{
    *d      = s;
    d->name = s.name;
    d->hack.clear();
    for (int i = 0; i < ebtsNR; i++)
    {
        d->rb[i].b.clear();
    }
    mergeAtomAndBondModifications(s, d);
}

#undef safe_strdup
