/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements test of exclusionblock routines
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_topology
 */
#include "gmxpre.h"

#include "gromacs/topology/exclusionblocks.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/topology/block.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/cmdlinetest.h"

namespace gmx
{
namespace testing
{
namespace
{

//! Add a new group to t_blocka
void addGroupToBlocka(t_blocka* b, gmx::ArrayRef<const int> indices)
{
    srenew(b->index, b->nr + 2);
    srenew(b->a, b->nra + indices.size());
    for (Index i = 0; i < indices.ssize(); i++)
    {
        b->a[b->nra++] = indices[i];
    }
    b->nr++;
    b->index[b->nr] = b->nra;
}

//! Fill ExclusionBlock with data.
int fillExclusionBlock(gmx::ArrayRef<ExclusionBlock> b)
{
    std::vector<std::vector<int>> indices = { { 0, 4, 7 }, { 1, 5, 8, 10 }, { 2, 6, 9, 11, 12 } };
    int                           nra     = 0;
    for (Index i = 0; i < b.ssize(); i++)
    {
        b[i].atomNumber.clear();
        for (const auto& j : indices[i])
        {
            b[i].atomNumber.push_back(j);
        }
        nra += b[i].nra();
    }
    return nra;
}

//! Fill the t_blocka with some datastructures
void makeTestBlockAData(t_blocka* ba)
{
    init_blocka(ba);

    std::vector<int> indices = { 12, 11, 9, 6, 2 };
    addGroupToBlocka(ba, indices);
    indices = { 10, 8, 5, 1 };
    addGroupToBlocka(ba, indices);
    indices = { 7, 4, 0 };
    addGroupToBlocka(ba, indices);
}

//! Return ListOfLists filled with some datastructures
ListOfLists<int> makeTestListOfLists()
{
    ListOfLists<int> list;

    std::vector<int> indices = { 12, 11, 9, 6, 2 };
    list.pushBack(indices);
    indices = { 10, 8, 5, 1 };
    list.pushBack(indices);
    indices = { 7, 4, 0 };
    list.pushBack(indices);

    return list;
}

class ExclusionBlockTest : public ::testing::Test
{
public:
    ExclusionBlockTest()
    {
        const int natom = 3;
        makeTestBlockAData(&ba_);
        list_ = makeTestListOfLists();
        b_.resize(natom);
    }
    ~ExclusionBlockTest() override { done_blocka(&ba_); }

    void compareBlocks()
    {
        for (Index i = 0; i < gmx::ssize(b_); i++)
        {
            int index = ba_.index[i];
            for (int j = 0; j < b_[i].nra(); j++)
            {
                int pos = index + j;
                EXPECT_EQ(b_[i].atomNumber[j], ba_.a[pos])
                        << "Block mismatch at " << i << " , " << j << ".";
            }
        }
    }

    void compareBlocksAndList()
    {
        GMX_RELEASE_ASSERT(gmx::ssize(b_) == list_.ssize(), "The list counts should match");
        for (Index i = 0; i < gmx::ssize(b_); i++)
        {
            gmx::ArrayRef<const int> jList = list_[i];
            ASSERT_EQ(b_[i].nra(), jList.ssize()) << "Block size mismatch at " << i << ".";
            EXPECT_THAT(b_[i].atomNumber, ::testing::Pointwise(::testing::Eq(), jList));
        }
    }

protected:
    t_blocka                    ba_;
    ListOfLists<int>            list_;
    std::vector<ExclusionBlock> b_;
};

TEST_F(ExclusionBlockTest, ConvertBlockAToExclusionBlocks)
{
    blockaToExclusionBlocks(&ba_, b_);
    compareBlocks();
}

TEST_F(ExclusionBlockTest, ConvertExclusionBlockToBlocka)
{
    int nra = fillExclusionBlock(b_);
    srenew(ba_.a, nra + 1);
    srenew(ba_.index, b_.size() + 1);
    exclusionBlocksToBlocka(b_, &ba_);
    compareBlocks();
}

TEST_F(ExclusionBlockTest, MergeExclusions)
{
    mergeExclusions(&list_, b_);
    compareBlocksAndList();
}

} // namespace

} // namespace testing

} // namespace gmx
