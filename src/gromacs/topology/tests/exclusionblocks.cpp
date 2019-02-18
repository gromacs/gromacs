/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements test of exclusionblock routines
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include "gromacs/topology/exclusionblocks.h"

#include <gtest/gtest.h>

#include "gromacs/topology/block.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
namespace
{

void hackBlocka(t_blocka *b, gmx::ArrayRef<const int> index)
{
    srenew(b->index, b->nr+2);
    srenew(b->a, b->nra+index.size());
    for (int i = 0; (i < index.ssize()); i++)
    {
        b->a[b->nra++] = index[i];
    }
    b->nr++;
    b->index[b->nr] = b->nra;
}

void fillBlocka(t_blocka *b)
{
    init_blocka(b);

    std::vector<int> index = {0, 4, 7};

    hackBlocka(b, index);
    index = {1, 5, 8, 10};
    hackBlocka(b, index);
    index = {2, 6, 9, 11, 12};
    hackBlocka(b, index);
}

void initialize(ExclusionBlocks *b, int natom)
{
    initExclusionBlocks(b, natom);
    EXPECT_EQ(b->nr, natom);
    ASSERT_NE(b->nra, nullptr);
    ASSERT_NE(b->a, nullptr);
    for (int i = 0; i < natom; i++)
    {
        EXPECT_EQ(b->a[i], nullptr);
    }
}

int fillExclusionBlock(ExclusionBlocks *b)
{
    snew(b->nra, b->nr);
    snew(b->a, b->nr);
    std::vector < std::vector < int>> index = {{0, 4, 7}, {1, 5, 8, 10}, {2, 6, 9, 11, 12}};
    int nra = 0;
    for (int i = 0; i < b->nr; i++)
    {
        snew(b->a[i], index[i].size());
        for (unsigned j = 0; j < index[i].size(); j++)
        {
            b->a[i][j] = index[i][j];
        }
        b->nra[i] = index[i].size();
        nra      += b->nra[i];
    }
    return nra;
}

void compareBlocks(ExclusionBlocks *b, t_blocka *ba)
{
    for (int i = 0; i < b->nr; i++)
    {
        int index  = ba->index[i];
        for (int j = 0; j < b->nra[i]; j++)
        {
            int pos = index + j;
            EXPECT_EQ(b->a[i][j], ba->a[pos]);
        }
    }
}

void cleanUp(ExclusionBlocks *b, t_blocka *ba)
{
    doneExclusionBlocks(b);
    if (ba)
    {
        done_blocka(ba);
    }
}

TEST(ExclusionBlockTest, EmptyOnInit)
{
    ExclusionBlocks b;
    int             natom  = 3;
    initialize(&b, natom);
    cleanUp(&b, nullptr);
}

TEST(ExclusionBlockTest, ConvertBlockAToExclusionBlocks)
{
    t_blocka ba;
    fillBlocka(&ba);

    int             natom = 3;
    ExclusionBlocks b;
    initialize(&b, natom);

    ASSERT_EQ(b.nr, ba.nr);

    blockaToExclusionBlocks(&ba, &b);

    compareBlocks(&b, &ba);

    cleanUp(&b, &ba);
}

TEST(ExclusionBlockTest, ConvertExclusionBlockToBlocka)
{
    int             natom = 3;
    ExclusionBlocks b;
    initialize(&b, natom);
    int             nra = fillExclusionBlock(&b);

    t_blocka        ba;
    init_blocka(&ba);
    snew(ba.a, nra+1);
    snew(ba.index, b.nr);
    exclusionBlocksToBlocka(&b, &ba);

    compareBlocks(&b, &ba);

    cleanUp(&b, &ba);
}

TEST(ExclusionBlockTest, MergeExclusions)
{
    t_blocka ba;
    fillBlocka(&ba);

    int             natom = 3;
    ExclusionBlocks b;
    initialize(&b, natom);

    ASSERT_EQ(b.nr, ba.nr);

    mergeExclusions(&ba, &b);

    compareBlocks(&b, &ba);

    cleanUp(&b, &ba);
}

}  // namespace

}  // namespace gmx
