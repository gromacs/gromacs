/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "block.h"

#include <cstdio>

#include <algorithm>

#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

void init_block(t_block *block)
{
    block->nr           = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0]     = 0;
}

void init_blocka(t_blocka *block)
{
    block->nr           = 0;
    block->nra          = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0]     = 0;
    block->nalloc_a     = 0;
    block->a            = nullptr;
}

t_blocka *new_blocka(void)
{
    t_blocka *block;

    snew(block, 1);
    snew(block->index, 1);

    return block;
}

void done_block(t_block *block)
{
    block->nr    = 0;
    sfree(block->index);
    block->nalloc_index = 0;
}

void done_blocka(t_blocka *block)
{
    block->nr    = 0;
    block->nra   = 0;
    sfree(block->index);
    sfree(block->a);
    block->index        = nullptr;
    block->a            = nullptr;
    block->nalloc_index = 0;
    block->nalloc_a     = 0;
}

void stupid_fill_block(t_block *grp, int natom, gmx_bool bOneIndexGroup)
{
    if (bOneIndexGroup)
    {
        grp->nalloc_index = 2;
        snew(grp->index, grp->nalloc_index);
        grp->index[0] = 0;
        grp->index[1] = natom;
        grp->nr       = 1;
    }
    else
    {
        grp->nalloc_index = natom+1;
        snew(grp->index, grp->nalloc_index);
        snew(grp->index, natom+1);
        for (int i = 0; i <= natom; ++i)
        {
            grp->index[i] = i;
        }
        grp->nr = natom;
    }
}

void stupid_fill_blocka(t_blocka *grp, int natom)
{
    grp->nalloc_a = natom;
    snew(grp->a, grp->nalloc_a);
    for (int i = 0; i < natom; ++i)
    {
        grp->a[i] = i;
    }
    grp->nra = natom;

    grp->nalloc_index = natom + 1;
    snew(grp->index, grp->nalloc_index);
    for (int i = 0; i <= natom; ++i)
    {
        grp->index[i] = i;
    }
    grp->nr = natom;
}

void copy_blocka(const t_blocka *src, t_blocka *dest)
{
    dest->nr           = src->nr;
    /* Workaround for inconsistent handling of nalloc_index in
     * other parts of the code. Often nalloc_index and nalloc_a
     * are not set.
     */
    dest->nalloc_index = std::max(src->nalloc_index, dest->nr + 1);
    snew(dest->index, dest->nalloc_index);
    for (int i = 0; i < dest->nr+1; ++i)
    {
        dest->index[i] = src->index[i];
    }
    dest->nra      = src->nra;
    /* See above. */
    dest->nalloc_a = std::max(src->nalloc_a, dest->nra);
    snew(dest->a, dest->nalloc_a);
    for (int i = 0; i < dest->nra; ++i)
    {
        dest->a[i] = src->a[i];
    }
}

static int pr_block_title(FILE *fp, int indent, const char *title, const t_block *block)
{
    if (available(fp, block, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "nr=%d\n", block->nr);
    }
    return indent;
}

static int pr_blocka_title(FILE *fp, int indent, const char *title, const t_blocka *block)
{
    if (available(fp, block, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "nr=%d\n", block->nr);
        pr_indent(fp, indent);
        fprintf(fp, "nra=%d\n", block->nra);
    }
    return indent;
}

static void low_pr_blocka(FILE *fp, int indent, const char *title, const t_blocka *block, gmx_bool bShowNumbers)
{
    int i;

    if (available(fp, block, indent, title))
    {
        indent = pr_blocka_title(fp, indent, title, block);
        for (i = 0; i <= block->nr; i++)
        {
            pr_indent(fp, indent+INDENT);
            fprintf(fp, "%s->index[%d]=%d\n",
                    title, bShowNumbers ? i : -1, block->index[i]);
        }
        for (i = 0; i < block->nra; i++)
        {
            pr_indent(fp, indent+INDENT);
            fprintf(fp, "%s->a[%d]=%d\n",
                    title, bShowNumbers ? i : -1, block->a[i]);
        }
    }
}

void pr_block(FILE *fp, int indent, const char *title, const t_block *block, gmx_bool bShowNumbers)
{
    int i, start;

    if (available(fp, block, indent, title))
    {
        indent = pr_block_title(fp, indent, title, block);
        start  = 0;
        if (block->index[start] != 0)
        {
            fprintf(fp, "block->index[%d] should be 0\n", start);
        }
        else
        {
            for (i = 0; i < block->nr; i++)
            {
                int end  = block->index[i+1];
                pr_indent(fp, indent);
                if (end <= start)
                {
                    fprintf(fp, "%s[%d]={}\n", title, i);
                }
                else
                {
                    fprintf(fp, "%s[%d]={%d..%d}\n",
                            title, bShowNumbers ? i : -1,
                            bShowNumbers ? start : -1, bShowNumbers ? end-1 : -1);
                }
                start = end;
            }
        }
    }
}

void pr_blocka(FILE *fp, int indent, const char *title, const t_blocka *block, gmx_bool bShowNumbers)
{
    int i, j, ok, size, start, end;

    if (available(fp, block, indent, title))
    {
        indent = pr_blocka_title(fp, indent, title, block);
        start  = 0;
        end    = start;
        if ((ok = (block->index[start] == 0)) == 0)
        {
            fprintf(fp, "block->index[%d] should be 0\n", start);
        }
        else
        {
            for (i = 0; i < block->nr; i++)
            {
                end  = block->index[i+1];
                size = pr_indent(fp, indent);
                if (end <= start)
                {
                    size += fprintf(fp, "%s[%d]={", title, i);
                }
                else
                {
                    size += fprintf(fp, "%s[%d][%d..%d]={",
                                    title, bShowNumbers ? i : -1,
                                    bShowNumbers ? start : -1, bShowNumbers ? end-1 : -1);
                }
                for (j = start; j < end; j++)
                {
                    if (j > start)
                    {
                        size += fprintf(fp, ", ");
                    }
                    if ((size) > (USE_WIDTH))
                    {
                        fprintf(fp, "\n");
                        size = pr_indent(fp, indent+INDENT);
                    }
                    size += fprintf(fp, "%d", block->a[j]);
                }
                fprintf(fp, "}\n");
                start = end;
            }
        }
        if ((end != block->nra) || (!ok))
        {
            pr_indent(fp, indent);
            fprintf(fp, "tables inconsistent, dumping complete tables:\n");
            low_pr_blocka(fp, indent, title, block, bShowNumbers);
        }
    }
}
