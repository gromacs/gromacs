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

#include "gromacs/legacyheaders/splitter.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

typedef struct {
    int atom, sid;
} t_sid;

static int sid_comp(const void *a, const void *b)
{
    t_sid *sa, *sb;
    int    dd;

    sa = (t_sid *)a;
    sb = (t_sid *)b;

    dd = sa->sid-sb->sid;
    if (dd == 0)
    {
        return (sa->atom-sb->atom);
    }
    else
    {
        return dd;
    }
}

static int mk_grey(egCol egc[], t_graph *g, int *AtomI,
                   int maxsid, t_sid sid[])
{
    int  j, ng, ai, aj, g0;

    ng = 0;
    ai = *AtomI;

    g0 = g->at_start;
    /* Loop over all the bonds */
    for (j = 0; (j < g->nedge[ai]); j++)
    {
        aj = g->edge[ai][j]-g0;
        /* If there is a white one, make it gray and set pbc */
        if (egc[aj] == egcolWhite)
        {
            if (aj < *AtomI)
            {
                *AtomI = aj;
            }
            egc[aj] = egcolGrey;

            /* Check whether this one has been set before... */
            range_check(aj+g0, 0, maxsid);
            range_check(ai+g0, 0, maxsid);
            if (sid[aj+g0].sid != -1)
            {
                gmx_fatal(FARGS, "sid[%d]=%d, sid[%d]=%d, file %s, line %d",
                          ai, sid[ai+g0].sid, aj, sid[aj+g0].sid, __FILE__, __LINE__);
            }
            else
            {
                sid[aj+g0].sid  = sid[ai+g0].sid;
                sid[aj+g0].atom = aj+g0;
            }
            ng++;
        }
    }
    return ng;
}

static int first_colour(int fC, egCol Col, t_graph *g, egCol egc[])
/* Return the first node with colour Col starting at fC.
 * return -1 if none found.
 */
{
    int i;

    for (i = fC; (i < g->nnodes); i++)
    {
        if ((g->nedge[i] > 0) && (egc[i] == Col))
        {
            return i;
        }
    }

    return -1;
}

static int mk_sblocks(FILE *fp, t_graph *g, int maxsid, t_sid sid[])
{
    int     ng, nnodes;
    int     nW, nG, nB; /* Number of Grey, Black, White	*/
    int     fW, fG;     /* First of each category	*/
    egCol  *egc = NULL; /* The colour of each node	*/
    int     g0, nblock;

    if (!g->nbound)
    {
        return 0;
    }

    nblock = 0;

    nnodes = g->nnodes;
    snew(egc, nnodes);

    g0 = g->at_start;
    nW = g->nbound;
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
        /* Find the first white, this will allways be a larger
         * number than before, because no nodes are made white
         * in the loop
         */
        if ((fW = first_colour(fW, egcolWhite, g, egc)) == -1)
        {
            gmx_fatal(FARGS, "No WHITE nodes found while nW=%d\n", nW);
        }

        /* Make the first white node grey, and set the block number */
        egc[fW]        = egcolGrey;
        range_check(fW+g0, 0, maxsid);
        sid[fW+g0].sid = nblock++;
        nG++;
        nW--;

        /* Initial value for the first grey */
        fG = fW;

        if (debug)
        {
            fprintf(debug, "Starting G loop (nW=%d, nG=%d, nB=%d, total %d)\n",
                    nW, nG, nB, nW+nG+nB);
        }

        while (nG > 0)
        {
            if ((fG = first_colour(fG, egcolGrey, g, egc)) == -1)
            {
                gmx_fatal(FARGS, "No GREY nodes found while nG=%d\n", nG);
            }

            /* Make the first grey node black */
            egc[fG] = egcolBlack;
            nB++;
            nG--;

            /* Make all the neighbours of this black node grey
             * and set their block number
             */
            ng = mk_grey(egc, g, &fG, maxsid, sid);
            /* ng is the number of white nodes made grey */
            nG += ng;
            nW -= ng;
        }
    }
    sfree(egc);

    if (debug)
    {
        fprintf(debug, "Found %d shake blocks\n", nblock);
    }

    return nblock;
}


typedef struct {
    int first, last, sid;
} t_merge_sid;

static int ms_comp(const void *a, const void *b)
{
    t_merge_sid *ma = (t_merge_sid *)a;
    t_merge_sid *mb = (t_merge_sid *)b;
    int          d;

    d = ma->first-mb->first;
    if (d == 0)
    {
        return ma->last-mb->last;
    }
    else
    {
        return d;
    }
}

static int merge_sid(int at_start, int at_end, int nsid, t_sid sid[],
                     t_blocka *sblock)
{
    int          i, j, k, n, isid, ndel;
    t_merge_sid *ms;
    int          nChanged;

    /* We try to remdy the following problem:
     * Atom: 1  2  3  4  5 6 7 8 9 10
     * Sid:  0 -1  1  0 -1 1 1 1 1 1
     */

    /* Determine first and last atom in each shake ID */
    snew(ms, nsid);

    for (k = 0; (k < nsid); k++)
    {
        ms[k].first = at_end+1;
        ms[k].last  = -1;
        ms[k].sid   = k;
    }
    for (i = at_start; (i < at_end); i++)
    {
        isid = sid[i].sid;
        range_check(isid, -1, nsid);
        if (isid >= 0)
        {
            ms[isid].first = min(ms[isid].first, sid[i].atom);
            ms[isid].last  = max(ms[isid].last, sid[i].atom);
        }
    }
    qsort(ms, nsid, sizeof(ms[0]), ms_comp);

    /* Now merge the overlapping ones */
    ndel = 0;
    for (k = 0; (k < nsid); )
    {
        for (j = k+1; (j < nsid); )
        {
            if (ms[j].first <= ms[k].last)
            {
                ms[k].last  = max(ms[k].last, ms[j].last);
                ms[k].first = min(ms[k].first, ms[j].first);
                ms[j].sid   = -1;
                ndel++;
                j++;
            }
            else
            {
                k = j;
                j = k+1;
            }
        }
        if (j == nsid)
        {
            k++;
        }
    }
    for (k = 0; (k < nsid); k++)
    {
        while ((k < nsid-1) && (ms[k].sid == -1))
        {
            for (j = k+1; (j < nsid); j++)
            {
                memcpy(&(ms[j-1]), &(ms[j]), sizeof(ms[0]));
            }
            nsid--;
        }
    }

    for (k = at_start; (k < at_end); k++)
    {
        sid[k].atom = k;
        sid[k].sid  = -1;
    }
    sblock->nr           = nsid;
    sblock->nalloc_index = sblock->nr+1;
    snew(sblock->index, sblock->nalloc_index);
    sblock->nra      = at_end - at_start;
    sblock->nalloc_a = sblock->nra;
    snew(sblock->a, sblock->nalloc_a);
    sblock->index[0] = 0;
    for (k = n = 0; (k < nsid); k++)
    {
        sblock->index[k+1] = sblock->index[k] + ms[k].last - ms[k].first+1;
        for (j = ms[k].first; (j <= ms[k].last); j++)
        {
            range_check(n, 0, sblock->nra);
            sblock->a[n++] = j;
            range_check(j, 0, at_end);
            if (sid[j].sid == -1)
            {
                sid[j].sid = k;
            }
            else
            {
                fprintf(stderr, "Double sids (%d, %d) for atom %d\n", sid[j].sid, k, j);
            }
        }
    }
    assert(k == nsid);
    /* Removed 2007-09-04
       sblock->index[k+1] = natoms;
       for(k=0; (k<natoms); k++)
       if (sid[k].sid == -1)
        sblock->a[n++] = k;
       assert(n == natoms);
     */
    sblock->nra = n;
    assert(sblock->index[k] == sblock->nra);
    sfree(ms);

    return nsid;
}

void gen_sblocks(FILE *fp, int at_start, int at_end,
                 t_idef *idef, t_blocka *sblock,
                 gmx_bool bSettle)
{
    t_graph *g;
    int      i, i0, j, k, istart, n;
    t_sid   *sid;
    int      isid, nsid;

    g = mk_graph(NULL, idef, at_start, at_end, TRUE, bSettle);
    if (debug)
    {
        p_graph(debug, "Graaf Dracula", g);
    }
    snew(sid, at_end);
    for (i = at_start; (i < at_end); i++)
    {
        sid[i].atom =  i;
        sid[i].sid  = -1;
    }
    nsid = mk_sblocks(fp, g, at_end, sid);

    if (!nsid)
    {
        return;
    }

    /* Now sort the shake blocks... */
    qsort(sid+at_start, at_end-at_start, (size_t)sizeof(sid[0]), sid_comp);

    if (debug)
    {
        fprintf(debug, "Sorted shake block\n");
        for (i = at_start; (i < at_end); i++)
        {
            fprintf(debug, "sid[%5d] = atom:%5d sid:%5d\n", i, sid[i].atom, sid[i].sid);
        }
    }
    /* Now check how many are NOT -1, i.e. how many have to be shaken */
    for (i0 = at_start; (i0 < at_end); i0++)
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
    nsid = merge_sid(at_start, at_end, nsid, sid, sblock);
    /* Now sort the shake blocks again... */
    /*qsort(sid,natoms,(size_t)sizeof(sid[0]),sid_comp);*/

    /* Fill the sblock struct */
    /*  sblock->nr  = nsid;
       sblock->nra = natoms;
       srenew(sblock->a,sblock->nra);
       srenew(sblock->index,sblock->nr+1);

       i    = i0;
       isid = sid[i].sid;
       n    = k = 0;
       sblock->index[n++]=k;
       while (i < natoms) {
       istart = sid[i].atom;
       while ((i<natoms-1) && (sid[i+1].sid == isid))
       i++;*/
    /* After while: we found a new block, or are thru with the atoms */
    /*    for(j=istart; (j<=sid[i].atom); j++,k++)
        sblock->a[k]=j;
       sblock->index[n] = k;
       if (i < natoms-1)
        n++;
       if (n > nsid)
        gmx_fatal(FARGS,"Death Horror: nsid = %d, n= %d",nsid,n);
       i++;
       isid = sid[i].sid;
       }
     */
    sfree(sid);
    /* Due to unknown reason this free generates a problem sometimes */
    done_graph(g);
    sfree(g);
    if (debug)
    {
        fprintf(debug, "Done gen_sblocks\n");
    }
}
