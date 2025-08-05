/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Tests for the MpiComm reductions and the HierarchicalReducer sub-class
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/mpicomm.h"

#include "config.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/mpitest.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

// Tests that a reduction works
void testReduction(const MpiComm& mpiComm)
{
    std::vector<int> a(1, (1 << mpiComm.rank()));

    int sum = 0;
    for (int i = 0; i < mpiComm.size(); i++)
    {
        sum += (1 << i);
    }

    mpiComm.sumReduce(a);

    EXPECT_EQ(a[0], sum);
}

TEST(MpiComm, SingleRankNoComm)
{
    GMX_MPI_TEST(AllowAnyRankCount);

    MpiComm mpiComm(MpiComm::SingleRank{});

    EXPECT_EQ(mpiComm.size(), 1);
    EXPECT_EQ(mpiComm.isSerial(), true);
    EXPECT_EQ(mpiComm.isParallel(), false);

    // This tests that reductions also work with a single rank and without an MPI_Comm
    testReduction(mpiComm);
}

TEST(MpiComm, CommWorld)
{
    GMX_MPI_TEST(AllowAnyRankCount);

    MpiComm mpiComm(MPI_COMM_WORLD);

    EXPECT_EQ(mpiComm.isSerial(), mpiComm.size() == 1);
    EXPECT_EQ(mpiComm.isParallel(), mpiComm.size() > 1);

    testReduction(mpiComm);
}

// Creates and MpiComm for comm.
// Initializes an HierarchicalReducer using \p physicalNodeId.
// Checks whether the presence of the HierarchicalReducer matches \p expectHierarchicalReducer.
// Tests that a reduction works
void testMpiComm(MPI_Comm comm, const int physicalNodeId, const bool expectHierarchicalReducer)
{
    MpiComm mpiComm(comm);

    mpiComm.initializeHierarchicalReductions(physicalNodeId);

    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), expectHierarchicalReducer);

    testReduction(mpiComm);
}

// Each rank is on the same node: no hierarchy
TEST(MpiComm, CreatesHierarchicalReducer4Ranks1Node)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    testMpiComm(MPI_COMM_WORLD, 0, false);
}

// Each rank is on a separate node: no hierarchy
TEST(MpiComm, CreatesHierarchicalReducer4Ranks4Nodes)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    testMpiComm(MPI_COMM_WORLD, rank, false);
}

// Two ranks per node: hierarchy
TEST(MpiComm, CreatesHierarchicalReducer4Ranks2EqualNodes)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    testMpiComm(MPI_COMM_WORLD, rank / 2, true);
}

// One node with one rank, one with three ranks: hierarchy
TEST(MpiComm, CreatesHierarchicalReducer4Ranks2UnequalNodes)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    testMpiComm(MPI_COMM_WORLD, rank == 0 ? 0 : 1, true);
}

TEST(MpiComm, ConstructorThrowsOnNull)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    int isInitialized;
    MPI_Initialized(&isInitialized);
    ASSERT_EQ(isInitialized == 0, false);

    EXPECT_THROW_GMX(MpiComm(MPI_COMM_NULL), InvalidInputError);
}

// Two ranks per node: hierarchy
TEST(MpiComm, CopyConstructorWorks)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    MpiComm mpiComm(MPI_COMM_WORLD);

    mpiComm.initializeHierarchicalReductions(mpiComm.rank() / 2);

    MpiComm mpiCommCopy(mpiComm);

    EXPECT_EQ(mpiCommCopy.size(), mpiComm.size());

    EXPECT_EQ(mpiCommCopy.getHierarchicalReductionsReport().has_value(), true);
    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), true)
            << "ensure source reducer still present";
}

CLANG_DIAGNOSTIC_IGNORE("-Wself-assign-overloaded")
CLANG_DIAGNOSTIC_IGNORE("-Wself-move")
#if __GNUC__ > 12
GCC_DIAGNOSTIC_IGNORE("-Wself-move")
#endif

TEST(MpiComm, AssigmentWorks)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    MpiComm mpiComm(MPI_COMM_WORLD);

    mpiComm.initializeHierarchicalReductions(mpiComm.rank() / 2);

    mpiComm = mpiComm;

    EXPECT_EQ(mpiComm.size(), 4);

    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), true);
}

// Two ranks per node: hierarchy
TEST(MpiComm, AssigmentWorksForSelf)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    MpiComm mpiComm(MPI_COMM_WORLD);

    mpiComm.initializeHierarchicalReductions(mpiComm.rank() / 2);

    MpiComm mpiCommCopy = mpiComm;

    EXPECT_EQ(mpiCommCopy.size(), mpiComm.size());

    EXPECT_EQ(mpiCommCopy.getHierarchicalReductionsReport().has_value(), true);
    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), true)
            << "ensure source reducer still present";
}

// A pass of the static analyzer build/check complains about use of move and self move.
// Suppression with #ifndef __clang_analyzer__ wasn't effective enough,
// so we also check for GMX_CLANG_ANALYZER
#if !defined(__clang_analyzer__) && !GMX_CLANG_ANALYZER

TEST(MpiComm, MoveAssigmentWorks)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    MpiComm mpiComm(MPI_COMM_WORLD);

    mpiComm.initializeHierarchicalReductions(mpiComm.rank() / 2);

    MpiComm mpiCommCopy(MpiComm::SingleRank{});
    mpiCommCopy = std::move(mpiComm);

    EXPECT_EQ(mpiCommCopy.size(), 4);

    EXPECT_EQ(mpiCommCopy.getHierarchicalReductionsReport().has_value(), true);
    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), false);
}

TEST(MpiComm, MoveAssigmentWorksForSelf)
{
    GMX_MPI_TEST(RequireRankCount<4>);

    MpiComm mpiComm(MPI_COMM_WORLD);

    mpiComm.initializeHierarchicalReductions(mpiComm.rank() / 2);

    mpiComm = std::move(mpiComm);

    EXPECT_EQ(mpiComm.size(), 4);

    EXPECT_EQ(mpiComm.getHierarchicalReductionsReport().has_value(), true);
}

#endif // !defined(__clang_analyzer__) && !GMX_CLANG_ANALYZER

CLANG_DIAGNOSTIC_RESET
#if __GNUC__ > 12
GCC_DIAGNOSTIC_RESET
#endif

} // namespace

} // namespace test

} // namespace gmx
