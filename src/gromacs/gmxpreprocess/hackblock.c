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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "hackblock.h"

#include <string.h>

#include "gromacs/legacyheaders/names.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/* these MUST correspond to the enum in hackblock.h */
const char *btsNames[ebtsNR] = { "bonds", "angles", "dihedrals", "impropers", "exclusions", "cmap" };
const int   btsNiatoms[ebtsNR] = { 2,       3,        4,           4,           2,             5 };

static void free_t_bonded(t_rbonded *rb)
{
    int i;

    for (i = 0; i < MAXATOMLIST; i++)
    {
        sfree(rb->a[i]);
    }
    sfree(rb->s);
}

static void free_t_bondeds(t_rbondeds *rbs)
{
    int i;

    for (i = 0; i < rbs->nb; i++)
    {
        free_t_bonded(&rbs->b[i]);
    }
    sfree(rbs->b);
    rbs->b  = NULL;
    rbs->nb = 0;
}

void free_t_restp(int nrtp, t_restp **rtp)
{
    int i, j;

    for (i = 0; i < nrtp; i++)
    {
        sfree((*rtp)[i].resname);
        sfree((*rtp)[i].atom);
        for (j = 0; j < (*rtp)[i].natom; j++)
        {
            sfree(*(*rtp)[i].atomname[j]);
            sfree((*rtp)[i].atomname[j]);
        }
        sfree((*rtp)[i].atomname);
        sfree((*rtp)[i].cgnr);
        for (j = 0; j < ebtsNR; j++)
        {
            free_t_bondeds(&(*rtp)[i].rb[j]);
        }
    }
    sfree(*rtp);
}

void free_t_hack(int nh, t_hack **h)
{
    int i, j;

    for (i = 0; i < nh; i++)
    {
        sfree((*h)[i].oname);
        sfree((*h)[i].nname);
        sfree((*h)[i].atom);
        for (j = 0; j < 4; j++)
        {
            sfree((*h)[i].a[j]);
        }
    }
    sfree(*h);
    *h = NULL;
}

void free_t_hackblock(int nhb, t_hackblock **hb)
{
    int i, j;

    for (i = 0; i < nhb; i++)
    {
        sfree((*hb)[i].name);
        free_t_hack((*hb)[i].nhack, &(*hb)[i].hack);
        for (j = 0; j < ebtsNR; j++)
        {
            free_t_bondeds(&(*hb)[i].rb[j]);
        }
    }
    sfree(*hb);
}

void clear_t_hackblock(t_hackblock *hb)
{
    int i;

    hb->name    = NULL;
    hb->nhack   = 0;
    hb->maxhack = 0;
    hb->hack    = NULL;
    for (i = 0; i < ebtsNR; i++)
    {
        hb->rb[i].nb = 0;
        hb->rb[i].b  = NULL;
    }
}

void clear_t_hack(t_hack *hack)
{
    int i;

    hack->nr    = 0;
    hack->oname = NULL;
    hack->nname = NULL;
    hack->atom  = NULL;
    hack->cgnr  = NOTSET;
    hack->tp    = 0;
    hack->nctl  = 0;
    for (i = 0; i < 4; i++)
    {
        hack->a[i]  = NULL;
    }
    for (i = 0; i < DIM; i++)
    {
        hack->newx[i] = NOTSET;
    }
}

#define safe_strdup(str) ((str != NULL) ? gmx_strdup(str) : NULL)

static void copy_t_rbonded(t_rbonded *s, t_rbonded *d)
{
    int i;

    for (i = 0; i < MAXATOMLIST; i++)
    {
        d->a[i] = safe_strdup(s->a[i]);
    }
    d->s     = safe_strdup(s->s);
    d->match = s->match;
}

static gmx_bool contains_char(t_rbonded *s, char c)
{
    int      i;
    gmx_bool bRet;

    bRet = FALSE;
    for (i = 0; i < MAXATOMLIST; i++)
    {
        if (s->a[i] && s->a[i][0] == c)
        {
            bRet = TRUE;
        }
    }

    return bRet;
}

int
rbonded_find_atoms_in_list(t_rbonded *b, t_rbonded blist[], int nlist, int natoms)
{
    int      i, k;
    int      foundPos = -1;
    gmx_bool atomsMatch;

    for (i = 0; i < nlist && foundPos < 0; i++)
    {
        atomsMatch = TRUE;
        for (k = 0; k < natoms && atomsMatch; k++)
        {
            atomsMatch = atomsMatch && !strcmp(b->a[k], blist[i].a[k]);
        }
        /* Try reverse if forward match did not work */
        if (!atomsMatch)
        {
            atomsMatch = TRUE;
            for (k = 0; k < natoms && atomsMatch; k++)
            {
                atomsMatch = atomsMatch && !strcmp(b->a[k], blist[i].a[natoms-1-k]);
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
            if (!strcmp(b->s, blist[i].s))
            {
                gmx_warning("Duplicate line found in or between hackblock and rtp entries");
            }
        }
    }
    return foundPos;
}

gmx_bool merge_t_bondeds(t_rbondeds s[], t_rbondeds d[], gmx_bool bMin, gmx_bool bPlus)
{
    int      i, j;
    gmx_bool bBondsRemoved;
    int      nbHackblockStart;
    int      index;

    bBondsRemoved = FALSE;
    for (i = 0; i < ebtsNR; i++)
    {
        if (s[i].nb > 0)
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
            nbHackblockStart = d[i].nb;

            /* make space */
            srenew(d[i].b, d[i].nb + s[i].nb);
            for (j = 0; j < s[i].nb; j++)
            {
                /* Check if this bonded string already exists before adding.
                 * We are merging from the main RTP to the hackblocks, so this
                 * will mean the hackblocks overwrite the man RTP, as intended.
                 */
                index = rbonded_find_atoms_in_list(&s[i].b[j], d[i].b, d[i].nb, btsNiatoms[i]);
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
                    if (!(bMin && contains_char(&s[i].b[j], '-'))
                        && !(bPlus && contains_char(&s[i].b[j], '+')))
                    {
                        copy_t_rbonded(&s[i].b[j], &d[i].b[ d[i].nb ]);
                        d[i].nb++;
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

void copy_t_restp(t_restp *s, t_restp *d)
{
    int i;

    *d         = *s;
    d->resname = safe_strdup(s->resname);
    snew(d->atom, s->natom);
    for (i = 0; i < s->natom; i++)
    {
        d->atom[i] = s->atom[i];
    }
    snew(d->atomname, s->natom);
    for (i = 0; i < s->natom; i++)
    {
        snew(d->atomname[i], 1);
        *d->atomname[i] = safe_strdup(*s->atomname[i]);
    }
    snew(d->cgnr, s->natom);
    for (i = 0; i < s->natom; i++)
    {
        d->cgnr[i] = s->cgnr[i];
    }
    for (i = 0; i < ebtsNR; i++)
    {
        d->rb[i].type = s->rb[i].type;
        d->rb[i].nb   = 0;
        d->rb[i].b    = NULL;
    }
    merge_t_bondeds(s->rb, d->rb, FALSE, FALSE);
}

void copy_t_hack(t_hack *s, t_hack *d)
{
    int i;

    *d       = *s;
    d->oname = safe_strdup(s->oname);
    d->nname = safe_strdup(s->nname);
    if (s->atom)
    {
        snew(d->atom, 1);
        *(d->atom) = *(s->atom);
    }
    else
    {
        d->atom = NULL;
    }
    for (i = 0; i < 4; i++)
    {
        d->a[i] = safe_strdup(s->a[i]);
    }
    copy_rvec(s->newx, d->newx);
}

void merge_hacks_lo(int ns, t_hack *s, int *nd, t_hack **d)
{
    int i;

    if (ns)
    {
        srenew(*d, *nd + ns);
        for (i = 0; i < ns; i++)
        {
            copy_t_hack(&s[i], &(*d)[*nd + i]);
        }
        (*nd) += ns;
    }
}

void merge_hacks(t_hackblock *s, t_hackblock *d)
{
    merge_hacks_lo(s->nhack, s->hack, &d->nhack, &d->hack);
}

void merge_t_hackblock(t_hackblock *s, t_hackblock *d)
{
    merge_hacks(s, d);
    merge_t_bondeds(s->rb, d->rb, FALSE, FALSE);
}

void copy_t_hackblock(t_hackblock *s, t_hackblock *d)
{
    int i;

    *d       = *s;
    d->name  = safe_strdup(s->name);
    d->nhack = 0;
    d->hack  = NULL;
    for (i = 0; i < ebtsNR; i++)
    {
        d->rb[i].nb = 0;
        d->rb[i].b  = NULL;
    }
    merge_t_hackblock(s, d);
}

#undef safe_strdup

void dump_hb(FILE *out, int nres, t_hackblock hb[])
{
    int i, j, k, l;

#define SS(s) (s) ? (s) : "-"
#define SA(s) (s) ? "+" : ""
    fprintf(out, "t_hackblock\n");
    for (i = 0; i < nres; i++)
    {
        fprintf(out, "%3d %4s %2d %2d\n",
                i, SS(hb[i].name), hb[i].nhack, hb[i].maxhack);
        if (hb[i].nhack)
        {
            for (j = 0; j < hb[i].nhack; j++)
            {
                fprintf(out, "%d: %d %4s %4s %1s %2d %d %4s %4s %4s %4s\n",
                        j, hb[i].hack[j].nr,
                        SS(hb[i].hack[j].oname), SS(hb[i].hack[j].nname),
                        SA(hb[i].hack[j].atom), hb[i].hack[j].tp, hb[i].hack[j].cgnr,
                        SS(hb[i].hack[j].AI), SS(hb[i].hack[j].AJ),
                        SS(hb[i].hack[j].AK), SS(hb[i].hack[j].AL) );
            }
        }
        for (j = 0; j < ebtsNR; j++)
        {
            if (hb[i].rb[j].nb)
            {
                fprintf(out, " %c %d:", btsNames[j][0], hb[i].rb[j].nb);
                for (k = 0; k < hb[i].rb[j].nb; k++)
                {
                    fprintf(out, " [");
                    for (l = 0; l < btsNiatoms[j]; l++)
                    {
                        fprintf(out, " %s", hb[i].rb[j].b[k].a[l]);
                    }
                    fprintf(out, " %s]", SS(hb[i].rb[j].b[k].s));
                }
                fprintf(out, " Entry matched: %s\n", yesno_names[hb[i].rb[j].b[k].match]);
            }
        }
        fprintf(out, "\n");
    }
#undef SS
#undef SA
}

void init_t_protonate(t_protonate *protonate)
{
    protonate->bInit = FALSE;
}
