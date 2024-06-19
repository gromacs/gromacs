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

#include "splitter.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <vector>

#include "gromacs/pbcutil/mshift.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"

typedef struct
{
    int atom, sid;
} t_sid;

static bool sid_comp(const t_sid& sa, const t_sid& sb)
{
    if (sa.sid == sb.sid)
    {
        return sa.atom < sb.atom;
    }
    else
    {
        return sa.sid < sb.sid;
    }
}

static int mk_grey(gmx::ArrayRef<egCol> edgeColor, const t_graph* g, int* AtomI, int maxsid, gmx::ArrayRef<t_sid> sid)
{
    int ng, ai, g0;

    ng = 0;
    ai = *AtomI;

    g0 = g->edgeAtomBegin;
    /* Loop over all the bonds */
    for (int aj : g->edges[ai])
    {
        aj -= g0;
        /* If there is a white one, make it gray and set pbc */
        if (edgeColor[aj] == egcolWhite)
        {
            if (aj < *AtomI)
            {
                *AtomI = aj;
            }
            edgeColor[aj] = egcolGrey;

            /* Check whether this one has been set before... */
            range_check(aj + g0, 0, maxsid);
            range_check(ai + g0, 0, maxsid);
            if (sid[aj + g0].sid != -1)
            {
                gmx_fatal(FARGS,
                          "sid[%d]=%d, sid[%d]=%d, file %s, line %d",
                          ai,
                          sid[ai + g0].sid,
                          aj,
                          sid[aj + g0].sid,
                          __FILE__,
                          __LINE__);
            }
            else
            {
                sid[aj + g0].sid  = sid[ai + g0].sid;
                sid[aj + g0].atom = aj + g0;
            }
            ng++;
        }
    }
    return ng;
}

static int first_colour(const int fC, const egCol Col, const t_graph* g, gmx::ArrayRef<const egCol> edgeColor)
/* Return the first node with colour Col starting at fC.
 * return -1 if none found.
 */
{
    int i;

    for (i = fC; i < int(g->edges.size()); i++)
    {
        if (!g->edges[i].empty() && edgeColor[i] == Col)
        {
            return i;
        }
    }

    return -1;
}

static int mk_sblocks(FILE* fp, t_graph* g, int maxsid, gmx::ArrayRef<t_sid> sid)
{
    int ng;
    int nW, nG, nB; /* Number of Grey, Black, White	*/
    int fW, fG;     /* First of each category	*/
    int g0, nblock;

    if (!g->numConnectedAtoms)
    {
        return 0;
    }

    nblock = 0;

    std::vector<egCol> edgeColor(g->edges.size(), egcolWhite);

    g0 = g->edgeAtomBegin;
    nW = g->numConnectedAtoms;
    nG = 0;
    nB = 0;

    fW = 0;

    /* We even have a loop invariant:
     * nW+nG+nB == g->nbound
     */

    if (fp)
    {
        fprintf(fp, "Walking down the molecule graph to make constraint-blocks\n");
    }

    while (nW > 0)
    {
        /* Find the first white, this will always be a larger
         * number than before, because no nodes are made white
         * in the loop
         */
        if ((fW = first_colour(fW, egcolWhite, g, edgeColor)) == -1)
        {
            gmx_fatal(FARGS, "No WHITE nodes found while nW=%d\n", nW);
        }

        /* Make the first white node grey, and set the block number */
        edgeColor[fW] = egcolGrey;
        range_check(fW + g0, 0, maxsid);
        sid[fW + g0].sid = nblock++;
        nG++;
        nW--;

        /* Initial value for the first grey */
        fG = fW;

        if (debug)
        {
            fprintf(debug, "Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n", nW, nG, nB, nW + nG + nB);
        }

        while (nG > 0)
        {
            if ((fG = first_colour(fG, egcolGrey, g, edgeColor)) == -1)
            {
                gmx_fatal(FARGS, "No GREY nodes found while nG=%d\n", nG);
            }

            /* Make the first grey node black */
            edgeColor[fG] = egcolBlack;
            nB++;
            nG--;

            /* Make all the neighbours of this black node grey
             * and set their block number
             */
            ng = mk_grey(edgeColor, g, &fG, maxsid, sid);
            /* ng is the number of white nodes made grey */
            nG += ng;
            nW -= ng;
        }
    }

    if (debug)
    {
        fprintf(debug, "Found %d shake blocks\n", nblock);
    }

    return nblock;
}


typedef struct
{
    int first, last, sid;
} t_merge_sid;

static int ms_comp(const void* a, const void* b)
{
    const t_merge_sid* ma = reinterpret_cast<const t_merge_sid*>(a);
    const t_merge_sid* mb = reinterpret_cast<const t_merge_sid*>(b);
    int                d;

    d = ma->first - mb->first;
    if (d == 0)
    {
        return ma->last - mb->last;
    }
    else
    {
        return d;
    }
}

static gmx::ListOfLists<int> merge_sid(int at_start, int at_end, int nsid, gmx::ArrayRef<t_sid> sid)
{
    int i, j, isid;

    /* We try to remdy the following problem:
     * Atom: 1  2  3  4  5 6 7 8 9 10
     * Sid:  0 -1  1  0 -1 1 1 1 1 1
     */

    /* Determine first and last atom in each shake ID */
    std::vector<t_merge_sid> ms;
    ms.reserve(nsid);

    for (int k = 0; k < nsid; k++)
    {
        ms.push_back({ at_end + 1, -1, k });
    }
    for (i = at_start; (i < at_end); i++)
    {
        isid = sid[i].sid;
        range_check(isid, -1, nsid);
        if (isid >= 0)
        {
            ms[isid].first = std::min(ms[isid].first, sid[i].atom);
            ms[isid].last  = std::max(ms[isid].last, sid[i].atom);
        }
    }
    qsort(ms.data(), gmx::ssize(ms), sizeof(ms[0]), ms_comp);

    /* Now merge the overlapping ones */
    for (int k = 0; k < nsid;)
    {
        for (j = k + 1; (j < nsid);)
        {
            if (ms[j].first <= ms[k].last)
            {
                ms[k].last  = std::max(ms[k].last, ms[j].last);
                ms[k].first = std::min(ms[k].first, ms[j].first);
                ms[j].sid   = -1;
                j++;
            }
            else
            {
                k = j;
                j = k + 1;
            }
        }
        if (j == nsid)
        {
            k++;
        }
    }
    for (int k = 0; k < nsid; k++)
    {
        while ((k < nsid - 1) && (ms[k].sid == -1))
        {
            for (j = k + 1; (j < nsid); j++)
            {
                std::memcpy(&(ms[j - 1]), &(ms[j]), sizeof(ms[0]));
            }
            nsid--;
        }
    }

    for (int k = at_start; k < at_end; k++)
    {
        sid[k].atom = k;
        sid[k].sid  = -1;
    }
    gmx::ListOfLists<int> sblock;
    std::vector<int>      listForSid;
    for (int k = 0; k < nsid; k++)
    {
        listForSid.clear();
        for (j = ms[k].first; (j <= ms[k].last); j++)
        {
            listForSid.push_back(j);
            GMX_RELEASE_ASSERT(sid[j].sid == -1, "Can not have double sids for an atom");
            sid[j].sid = k;
        }
        sblock.pushBack(listForSid);
    }

    return sblock;
}

gmx::ListOfLists<int> gen_sblocks(FILE* fp, int at_end, const InteractionDefinitions& idef, const bool useSettles)
{
    t_graph* g;
    int      i, i0;
    int      nsid;

    g = mk_graph(nullptr, idef, at_end, TRUE, useSettles);
    if (debug)
    {
        p_graph(debug, "Graaf Dracula", g);
    }
    std::vector<t_sid> sid;
    sid.reserve(at_end);
    for (i = 0; i < at_end; i++)
    {
        sid.push_back({ i, -1 });
    }
    nsid = mk_sblocks(fp, g, at_end, sid);

    if (!nsid)
    {
        return {};
    }

    /* Now sort the shake blocks... */
    std::sort(sid.data(), sid.data() + at_end, sid_comp);

    if (debug)
    {
        fprintf(debug, "Sorted shake block\n");
        for (i = 0; i < at_end; i++)
        {
            fprintf(debug, "sid[%5d] = atom:%5d sid:%5d\n", i, sid[i].atom, sid[i].sid);
        }
    }
    /* Now check how many are NOT -1, i.e. how many have to be shaken */
    for (i0 = 0; i0 < at_end; i0++)
    {
        if (sid[i0].sid > -1)
        {
            break;
        }
    }

    /* Now we have the sids that have to be shaken. We'll check the min and
     * max atom numbers and this determines the shake block. DvdS 2007-07-19.
     * For the purpose of making boundaries all atoms in between need to be
     * part of the shake block too. There may be cases where blocks overlap
     * and they will have to be merged.
     */
    gmx::ListOfLists<int> sblock = merge_sid(0, at_end, nsid, sid);
    /* Due to unknown reason this free generates a problem sometimes */
    done_graph(g);
    if (debug)
    {
        fprintf(debug, "Done gen_sblocks\n");
    }

    return sblock;
}
