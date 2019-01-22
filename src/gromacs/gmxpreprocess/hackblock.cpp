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
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* these MUST correspond to the enum in hackblock.h */
const char *btsNames[ebtsNR] = { "bonds", "angles", "dihedrals", "impropers", "exclusions", "cmap" };
const int   btsNiatoms[ebtsNR] = { 2,       3,        4,           4,           2,             5 };

void free_t_restp(gmx::ArrayRef<t_restp> rtp)
{
    for (auto it = rtp.begin(); it != rtp.end(); it++)
    {
        for (int j = 0; j < it->natom(); j++)
        {
            sfree(*it->atomname[j]);
            sfree(it->atomname[j]);
        }
        it->atom.clear();
        it->atomname.clear();
        it->cgnr.clear();
    }
}

void clear_t_hackblock(t_hackblock *hb)
{
    hb->name.clear();
    hb->hack.clear();
    for (int i = 0; i < ebtsNR; i++)
    {
        hb->rb[i].b.clear();
    }
}

void clear_t_hack(t_hack *hack)
{
    hack->oname.clear();
    hack->nname.clear();
    hack->atom.clear();
    hack->cgnr  = NOTSET;
    hack->tp    = 0;
    hack->nctl  = 0;
    for (int i = 0; i < 4; i++)
    {
        hack->a[i].clear();
    }
    for (int i = 0; i < DIM; i++)
    {
        hack->newx[i] = NOTSET;
    }
}

#define safe_strdup(str) (((str) != NULL) ? gmx_strdup(str) : NULL)

static void copy_t_rbonded(const t_rbonded &s, t_rbonded *d)
{
    for (int i = 0; i < MAXATOMLIST; i++)
    {
        d->a[i] = s.a[i];
    }
    d->s     = s.s;
    d->match = s.match;
}


static bool contains_char(const t_rbonded &s, char c)
{
    bool     bRet;

    bRet = FALSE;
    for (int i = 0; i < MAXATOMLIST; i++)
    {
        if (!s.a[i].empty() && s.a[i][0] == c)
        {
            bRet = TRUE;
        }
    }

    return bRet;
}

static int
rbonded_find_atoms_in_list(const t_rbonded               &b,
                           gmx::ArrayRef<const t_rbonded> blist,
                           int                            natoms)
{
    int      foundPos = -1;
    bool     atomsMatch;

    for (int i = 0; i < blist.size() && foundPos < 0; i++)
    {
        atomsMatch = TRUE;
        for (int k = 0; k < natoms && atomsMatch; k++)
        {
            atomsMatch = atomsMatch && (b.a[k] ==  blist[i].a[k]);
        }
        /* Try reverse if forward match did not work */
        if (!atomsMatch)
        {
            atomsMatch = TRUE;
            for (int k = 0; k < natoms && atomsMatch; k++)
            {
                atomsMatch = atomsMatch && (b.a[k] == blist[i].a[natoms-1-k]);
            }
        }
        if (atomsMatch)
        {
            foundPos = i;
            /* If all the atoms AND all the parameters match, it is likely that
             * the user made a copy-and-paste mistake (since it would be much cheaper
             * to just bump the force constant 2x if you really want it twice).
             * Since we only have the unparsed string here we can only detect
             * EXACT matches (including identical whitespace).
             */
            if (b.s == blist[i].s)
            {
                gmx_warning("Duplicate line found in or between hackblock and rtp entries");
            }
        }
    }
    return foundPos;
}

bool merge_t_bondeds(const std::array<t_rbondeds, ebtsNR> &s,
                     std::array<t_rbondeds, ebtsNR> *d,
                     bool bMin, bool bPlus)
{
    bool     bBondsRemoved;
    int      nbHackblockStart;
    int      index;

    bBondsRemoved = FALSE;
    for (int i = 0; i < ebtsNR; i++)
    {
        if (s[i].nb() > 0)
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
            nbHackblockStart = d->at(i).nb();

            /* make space */
            for (int j = 0; j < s[i].nb(); j++)
            {
                /* Check if this bonded string already exists before adding.
                 * We are merging from the main RTP to the hackblocks, so this
                 * will mean the hackblocks overwrite the man RTP, as intended.
                 */
                index = rbonded_find_atoms_in_list(s[i].b[j], d->at(i).b, btsNiatoms[i]);
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
                    if (!(bMin && contains_char(s[i].b[j], '-'))
                        && !(bPlus && contains_char(s[i].b[j], '+')))
                    {
                        d->at(i).b.push_back(t_rbonded());
                        copy_t_rbonded(s[i].b[j], &d->at(i).b.back());
                    }
                    else if (i == ebtsBONDS)
                    {
                        bBondsRemoved = TRUE;
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

void copy_t_restp(const t_restp *s, t_restp *d)
{
    *d = *s;
    d->atom.clear();
    for (const auto &a : s->atom)
    {
        d->atom.push_back(a);
    }
    d->atomname.clear();
    for (const auto &anm : s->atomname)
    {
        char **tmp = nullptr;
        snew(tmp, 1);
        *tmp = safe_strdup(*anm);
        d->atomname.push_back(tmp);
    }
    d->cgnr.clear();
    for (const auto &c : s->cgnr)
    {
        d->cgnr.push_back(c);
    }
    for (int i = 0; i < ebtsNR; i++)
    {
        d->rb[i].type = s->rb[i].type;
    }
    merge_t_bondeds(s->rb, &d->rb, FALSE, FALSE);
}

void copy_t_hack(const t_hack &s, t_hack *d)
{
    *d       = s;
    d->oname = s.oname;
    d->nname = s.nname;
    d->atom.clear();
    for (const auto &a : s.atom)
    {
        d->atom.push_back(a);
    }
    for (int i = 0; i < 4; i++)
    {
        d->a[i] = s.a[i];
    }
    copy_rvec(s.newx, d->newx);
}

void merge_hacks_lo(gmx::ArrayRef<const t_hack> s, std::vector<t_hack> *d)
{
    for (const auto &h : s)
    {
        d->push_back(t_hack());
        copy_t_hack(h, &d->back());
    }
}

void merge_hacks(const t_hackblock &s, t_hackblock *d)
{
    merge_hacks_lo(s.hack, &d->hack);
}

void merge_t_hackblock(const t_hackblock &s, t_hackblock *d)
{
    merge_hacks(s, d);
    merge_t_bondeds(s.rb, &d->rb, FALSE, FALSE);
}

void copy_t_hackblock(const t_hackblock &s, t_hackblock *d)
{
    *d       = s;
    d->name  = s.name;
    d->hack.clear();
    for (int i = 0; i < ebtsNR; i++)
    {
        d->rb[i].b.clear();
    }
    merge_t_hackblock(s, d);
}

#undef safe_strdup

void dump_hb(FILE *out, int nres, t_hackblock hb[])
{
#define SS(s) (s) ? (s) : "-"
#define SA(s) (s) ? "+" : ""
    fprintf(out, "t_hackblock\n");
    for (int i = 0; i < nres; i++)
    {
        fprintf(out, "%3d %4s %2d\n",
                i, SS(hb[i].name.c_str()), hb[i].nhack());
        if (hb[i].nhack())
        {
            for (int j = 0; j < hb[i].nhack(); j++)
            {
                fprintf(out, "%d: %d %4s %4s %1s %2d %d %4s %4s %4s %4s\n",
                        j, hb[i].hack[j].nr(),
                        SS(hb[i].hack[j].oname.c_str()), SS(hb[i].hack[j].nname.c_str()),
                        SA(hb[i].hack[j].atom.data()), hb[i].hack[j].tp, hb[i].hack[j].cgnr,
                        SS(hb[i].hack[j].ai()), SS(hb[i].hack[j].aj()),
                        SS(hb[i].hack[j].ak()), SS(hb[i].hack[j].al()) );
            }
        }
        for (int j = 0; j < ebtsNR; j++)
        {
            if (hb[i].rb[j].nb())
            {
                fprintf(out, " %c %d:", btsNames[j][0], hb[i].rb[j].nb());
                int k;
                for (k = 0; k < hb[i].rb[j].nb(); k++)
                {
                    fprintf(out, " [");
                    for (int l = 0; l < btsNiatoms[j]; l++)
                    {
                        fprintf(out, " %s", hb[i].rb[j].b[k].a[l].c_str());
                    }
                    fprintf(out, " %s]", SS(hb[i].rb[j].b[k].s.c_str()));
                }
                fprintf(out, " Entry matched: %s\n", yesno_names[hb[i].rb[j].b[k].match]);
            }
        }
        fprintf(out, "\n");
    }
#undef SS
#undef SA
}
