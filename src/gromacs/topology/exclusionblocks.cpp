/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "exclusionblocks.h"

#include <algorithm>

#include "gromacs/topology/block.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

void initExclusionBlocks(ExclusionBlocks *b2, int natom)
{
    int i;

    b2->nr = natom;
    snew(b2->nra, b2->nr);
    snew(b2->a, b2->nr);
    for (i = 0; (i < b2->nr); i++)
    {
        b2->a[i] = nullptr;
    }
}

void doneExclusionBlocks(ExclusionBlocks *b2)
{
    int i;

    if (b2->nr)
    {
        for (i = 0; (i < b2->nr); i++)
        {
            sfree(b2->a[i]);
        }
        sfree(b2->a);
        sfree(b2->nra);
        b2->nr = 0;
    }
}

void blockaToExclusionBlocks(t_blocka *b, ExclusionBlocks *b2)
{
    int     i;
    int     j, a;

    for (i = 0; (i < b->nr); i++)
    {
        for (j = b->index[i]; (j < b->index[i+1]); j++)
        {
            a = b->a[j];
            srenew(b2->a[i], ++b2->nra[i]);
            b2->a[i][b2->nra[i]-1] = a;
        }
    }
}

void exclusionBlocksToBlocka(ExclusionBlocks *b2, t_blocka *b)
{
    int     i, nra;
    int     j;

    nra = 0;
    for (i = 0; (i < b2->nr); i++)
    {
        b->index[i] = nra;
        for (j = 0; (j < b2->nra[i]); j++)
        {
            b->a[nra+j] = b2->a[i][j];
        }
        nra += b2->nra[i];
    }
    /* terminate list */
    b->index[i] = nra;
}

void mergeExclusions(t_blocka *excl, ExclusionBlocks *b2)
{
    int     i, k;
    int     j;
    int     nra;

    if (!b2->nr)
    {
        return;
    }
    GMX_RELEASE_ASSERT(b2->nr == excl->nr, "Cannot merge exclusions for "
                       "blocks that do not describe the same number "
                       "of particles");

    /* First copy all entries from excl to b2 */
    blockaToExclusionBlocks(excl, b2);

    /* Count and sort the exclusions */
    nra = 0;
    for (i = 0; (i < b2->nr); i++)
    {
        if (b2->nra[i] > 0)
        {
            /* remove double entries */
            std::sort(b2->a[i], b2->a[i]+b2->nra[i]);
            k = 1;
            for (j = 1; (j < b2->nra[i]); j++)
            {
                if (b2->a[i][j] != b2->a[i][k-1])
                {
                    b2->a[i][k] = b2->a[i][j];
                    k++;
                }
            }
            b2->nra[i] = k;
            nra       += b2->nra[i];
        }
    }
    excl->nra = nra;
    srenew(excl->a, excl->nra);

    exclusionBlocksToBlocka(b2, excl);
}

} // namespace gmx
