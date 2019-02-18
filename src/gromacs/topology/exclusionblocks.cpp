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

#include "exclusionblocks.h"

#include <algorithm>

#include "gromacs/topology/block.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

std::vector<ExclusionBlock> blockaToExclusionBlocks(const t_blocka *b, gmx::ArrayRef<const ExclusionBlock> b2)
{
    std::vector<ExclusionBlock> blocks(b2.begin(), b2.end());

    for (int i = 0; (i < b->nr); i++)
    {
        for (int j = b->index[i]; (j < b->index[i+1]); j++)
        {
            blocks[i].atomNumber.push_back(b->a[j]);
        }
    }
    return blocks;
}

void exclusionBlocksToBlocka(gmx::ArrayRef<const ExclusionBlock> b2, t_blocka *b)
{
    int nra = 0;
    int i   = 0;
    for (const auto &block : b2)
    {
        b->index[i] = nra;
        int j = 0;
        for (const auto &entry : block.atomNumber)
        {
            b->a[nra+j] = entry;
            j++;
        }
        nra += block.nra();
        i++;
    }
    /* terminate list */
    b->index[i] = nra;
}

std::vector<ExclusionBlock> mergeExclusions(t_blocka *excl, gmx::ArrayRef<const ExclusionBlock> b2)
{
    std::vector<ExclusionBlock> newBlocks;

    if (b2.empty())
    {
        return newBlocks;
    }
    GMX_RELEASE_ASSERT(b2.ssize() == excl->nr, "Cannot merge exclusions for "
                       "blocks that do not describe the same number "
                       "of particles");

    /* First copy all entries from excl to b2 */
    newBlocks = blockaToExclusionBlocks(excl, b2);

    /* Count and sort the exclusions */
    int nra = 0;
    for (auto &block : newBlocks)
    {
        if (block.nra() > 0)
        {
            /* remove double entries */
            std::sort(block.atomNumber.begin(), block.atomNumber.end());
            for (auto atom = block.atomNumber.begin() + 1; atom != block.atomNumber.end(); )
            {
                auto prev = atom - 1;
                if (*prev == *atom)
                {
                    atom = block.atomNumber.erase(atom);
                }
                else
                {
                    ++atom;
                }
            }
            nra       += block.nra();
        }
    }
    excl->nra = nra;
    srenew(excl->a, excl->nra);

    exclusionBlocksToBlocka(b2, excl);
    return newBlocks;
}

} // namespace gmx
