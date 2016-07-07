/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * Tests for general functionality in gmx::GroupIndexManager and
 * gmx::GroupIndex, which is only accesible through the manager.
 *
 * \author Christian Blau <cblau@gwdg.de>
 */
#include "gmxpre.h"

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "gromacs/mdlib/groupindex.h"
#include "gromacs/mdlib/groupindexmanager.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

TEST(GroupIndexManager, CanAddEmptyGroupIndex)
{
    bool                                bParallel = false;
    GroupIndexManager                   manager(bParallel);
    GroupIndexManager::GroupIndexHandle empty_group         = manager.add( 0, nullptr);
    GroupIndexManager::GroupIndexHandle another_empty_group = manager.add( 0, nullptr);
};

TEST(GroupIndexManager, CanAddandReadGroupIndexIndices)
{
    bool                                bParallel = false;
    GroupIndexManager                   manager(bParallel);
    std::vector<int>                    index = {5, 10};
    std::vector<int>                    read_index;
    GroupIndexManager::GroupIndexHandle new_group = manager.add(index.size(), index.data());

    for (const auto &i : new_group->localIndex())
    {
        read_index.push_back(i);
    }
    ;
    ASSERT_THAT(read_index, testing::ContainerEq(index));
};

TEST(GroupIndexManager, CanClean)
{
    bool                                bParallel = false;
    GroupIndexManager                   manager(bParallel);
    GroupIndexManager::GroupIndexHandle new_group = manager.add(0, nullptr);
    /* Loose the handle to the index group */
    new_group.reset();
    /* Clean up al index groups with no handle outside the manager */
    manager.clean();
    ASSERT_EQ(manager.numberOfManagedGroupIndices(), 0);
};


} // namespace test
} // namespace gmx
