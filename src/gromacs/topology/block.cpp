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

#include "gromacs/topology/block.h"

#include <cstdio>

#include <algorithm>
#include <vector>

#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"


void init_block(t_block* block)
{
    block->nr           = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0] = 0;
}

void init_blocka(t_blocka* block)
{
    block->nr           = 0;
    block->nra          = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0] = 0;
    block->nalloc_a = 0;
    block->a        = nullptr;
}

void done_block(t_block* block)
{
    block->nr = 0;
    sfree(block->index);
    block->index        = nullptr;
    block->nalloc_index = 0;
}

void done_blocka(t_blocka* block)
{
    block->nr  = 0;
    block->nra = 0;
    sfree(block->index);
    sfree(block->a);
    block->index        = nullptr;
    block->a            = nullptr;
    block->nalloc_index = 0;
    block->nalloc_a     = 0;
}

void stupid_fill_block(t_block* grp, int natom, gmx_bool bOneIndexGroup)
{
    if (bOneIndexGroup)
    {
        grp->nalloc_index = 2;
        srenew(grp->index, grp->nalloc_index);
        grp->index[0] = 0;
        grp->index[1] = natom;
        grp->nr       = 1;
    }
    else
    {
        grp->nalloc_index = natom + 1;
        srenew(grp->index, grp->nalloc_index);
        for (int i = 0; i <= natom; ++i)
        {
            grp->index[i] = i;
        }
        grp->nr = natom;
    }
}

static int pr_block_title(FILE* fp, int indent, const char* title, const t_block* block)
{
    if (available(fp, block, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "nr=%d\n", block->nr);
    }
    return indent;
}

static int pr_blocka_title(FILE* fp, int indent, const char* title, const int numBlocks)
{
    indent = pr_title(fp, indent, title);
    pr_indent(fp, indent);
    fprintf(fp, "nr=%d\n", numBlocks);

    return indent;
}

static int pr_listoflists_title(FILE* fp, int indent, const char* title, const gmx::ListOfLists<int>* lists)
{
    if (available(fp, lists, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "numLists=%zu\n", lists->size());
        pr_indent(fp, indent);
        fprintf(fp, "numElements=%d\n", lists->numElements());
    }
    return indent;
}

void pr_block(FILE* fp, int indent, const char* title, const t_block* block, gmx_bool bShowNumbers)
{
    if (available(fp, block, indent, title))
    {
        indent    = pr_block_title(fp, indent, title, block);
        int start = 0;
        if (block->index[start] != 0)
        {
            fprintf(fp, "block->index[%d] should be 0\n", start);
        }
        else
        {
            for (int i = 0; i < block->nr; i++)
            {
                int end = block->index[i + 1];
                pr_indent(fp, indent);
                if (end <= start)
                {
                    fprintf(fp, "%s[%d]={}\n", title, i);
                }
                else
                {
                    fprintf(fp,
                            "%s[%d]={%d..%d}\n",
                            title,
                            bShowNumbers ? i : -1,
                            bShowNumbers ? start : -1,
                            bShowNumbers ? end - 1 : -1);
                }
                start = end;
            }
        }
    }
}

void pr_blocka(FILE* fp, int indent, const char* title, gmx::ArrayRef<const IndexGroup> blocks, gmx_bool bShowNumbers)
{
    indent = pr_blocka_title(fp, indent, title, gmx::ssize(blocks));

    for (int i = 0; i < gmx::ssize(blocks); i++)
    {
        const auto& block = blocks[i].particleIndices;
        int         size  = pr_indent(fp, indent);
        if (block.empty())
        {
            size += fprintf(fp, "%s[%d]={", title, i);
        }
        else
        {
            size += fprintf(fp, "%s[%d]={", title, bShowNumbers ? i : -1);
        }
        bool firstElement = true;
        for (const int a : block)
        {
            if (!firstElement)
            {
                size += fprintf(fp, ", ");
            }
            if ((size) > (USE_WIDTH))
            {
                fprintf(fp, "\n");
                size = pr_indent(fp, indent + INDENT);
            }
            size += fprintf(fp, "%d", a);
            firstElement = false;
        }
        fprintf(fp, "}\n");
    }
}

void pr_listoflists(FILE* fp, int indent, const char* title, const gmx::ListOfLists<int>* lists, gmx_bool bShowNumbers)
{
    if (available(fp, lists, indent, title))
    {
        indent = pr_listoflists_title(fp, indent, title, lists);
        for (gmx::Index i = 0; i < lists->ssize(); i++)
        {
            int                      size = pr_indent(fp, indent);
            gmx::ArrayRef<const int> list = (*lists)[i];
            if (list.empty())
            {
                size += fprintf(fp, "%s[%d]={", title, int(i));
            }
            else
            {
                size += fprintf(fp, "%s[%d][num=%zu]={", title, bShowNumbers ? int(i) : -1, list.size());
            }
            bool isFirst = true;
            for (const int j : list)
            {
                if (!isFirst)
                {
                    size += fprintf(fp, ", ");
                }
                if ((size) > (USE_WIDTH))
                {
                    fprintf(fp, "\n");
                    size = pr_indent(fp, indent + INDENT);
                }
                size += fprintf(fp, "%d", j);
                isFirst = false;
            }
            fprintf(fp, "}\n");
        }
    }
}
