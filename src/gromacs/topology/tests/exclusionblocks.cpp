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

static void hackBlocka(t_blocka *b, gmx::ArrayRef<const int> index)
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

static void fillBlocka(t_blocka *b)
{
    init_blocka(b);
    // Prepare a dummy group
    std::vector<int> index = {0, 4, 7};
    // Add the group to the index
    hackBlocka(b, index);
    index = {1, 5, 8, 10};
    hackBlocka(b, index);
    index = {2, 6, 9, 11, 12};
    hackBlocka(b, index);
}

static void initialize(gmx::ArrayRef<ExclusionBlock> b, int natom)
{
    EXPECT_EQ(b.size(), natom);
    for (auto &i : b)
    {
        EXPECT_TRUE(i.atomNumber.empty());
    }
}

TEST(ExclusionBlockTest, EmptyOnInit)
{
    int                         natom  = 3;
    std::vector<ExclusionBlock> b(natom);
    initialize(b, natom);
}

TEST(ExclusionBlockTest, ConvertBlockAToExclusionBlocks)
{
    t_blocka ba;
    fillBlocka(&ba);

    int                         natom = 3;
    std::vector<ExclusionBlock> b(natom);
    initialize(b, natom);

    ASSERT_EQ(b.size(), ba.nr);

    b = blockaToExclusionBlocks(&ba, b);

    for (int i = 0; i < natom; i++)
    {
        int index  = ba.index[i];
        for (int j = 0; j < b[i].nra(); j++)
        {
            int pos = index + j;
            EXPECT_EQ(b[i].atomNumber[j], ba.a[pos]);
        }
    }
}

TEST(ExclusionBlocks, MergeExclusions)
{
    t_blocka ba;
    fillBlocka(&ba);

    int                         natom = 3;
    std::vector<ExclusionBlock> b(natom);
    initialize(b, natom);

    ASSERT_EQ(b.size(), ba.nr);

    b = blockaToExclusionBlocks(&ba, b);

    b = mergeExclusions(&ba, b);

    for (int i = 0; i < natom; i++)
    {
        int index  = ba.index[i];
        for (int j = 0; j < b[i].nra(); j++)
        {
            int pos = index + j;
            EXPECT_EQ(b[i].atomNumber[j], ba.a[pos]);
        }
    }
}

}  // namespace

}  // namespace gmx
