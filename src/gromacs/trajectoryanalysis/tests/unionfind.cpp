/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for the union-find implementation in unionfind.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/modules/unionfind.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace
{

/********************************************************************
 * Tests for UnionFinder
 */

TEST(UnionFinderTest, WorksEmpty)
{
    using ::testing::IsEmpty;
    gmx::UnionFinder finder;
    finder.init(0);
    EXPECT_THAT(finder.allSizes(), IsEmpty());
}

TEST(UnionFinderTest, BasicMerges)
{
    using ::testing::Each;
    using ::testing::UnorderedElementsAre;
    gmx::UnionFinder finder;
    finder.init(7);
    EXPECT_EQ(7U, finder.allSizes().size());
    EXPECT_THAT(finder.allSizes(), Each(1));
    finder.merge(0, 1);
    finder.merge(2, 3);
    finder.merge(0, 4);
    finder.merge(3, 5);
    EXPECT_THAT(finder.allSizes(), UnorderedElementsAre(3, 3, 1));
}

TEST(UnionFinderTest, LargerMerges)
{
    using ::testing::UnorderedElementsAre;
    gmx::UnionFinder finder;
    finder.init(7);
    finder.merge(0, 1);
    finder.merge(2, 3);
    finder.merge(4, 5);
    finder.merge(3, 5);
    finder.merge(1, 5);
    EXPECT_THAT(finder.allSizes(), UnorderedElementsAre(6, 1));
}

TEST(UnionFinderTest, LongRightMerge)
{
    using ::testing::ElementsAre;
    gmx::UnionFinder finder;
    finder.init(7);
    finder.merge(0, 1);
    finder.merge(1, 2);
    finder.merge(2, 3);
    finder.merge(3, 4);
    finder.merge(4, 5);
    finder.merge(5, 6);
    EXPECT_THAT(finder.allSizes(), ElementsAre(7));
}

TEST(UnionFinderTest, LongLeftMerge)
{
    using ::testing::ElementsAre;
    gmx::UnionFinder finder;
    finder.init(7);
    finder.merge(5, 6);
    finder.merge(4, 5);
    finder.merge(3, 4);
    finder.merge(2, 3);
    finder.merge(1, 2);
    finder.merge(0, 1);
    EXPECT_THAT(finder.allSizes(), ElementsAre(7));
}

/********************************************************************
 * Tests for MappedUnionFinder
 */

TEST(MappedUnionFinderTest, BasicMerges)
{
    using ::testing::Each;
    const int              mapping[] = {1, 1, 2, 2, 4, 6};
    gmx::MappedUnionFinder finder;
    finder.initWithGroupIndices(mapping);
    EXPECT_EQ(4U, finder.allSizes().size());
    EXPECT_THAT(finder.allSizes(), Each(1));
    finder.mergeGroups(1, 4);
    finder.mergeGroups(2, 6);
    EXPECT_EQ(2U, finder.allSizes().size());
    EXPECT_THAT(finder.allSizes(), Each(2));
}

} // namespace
