/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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

#include <array>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/basenetwork.h"

#include "testutils/mpitest.h"

#include "threadaffinitytest.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::ThreadAffinityTestHelper;

//! Helper to establish test requirements on MPI ranks
class RequireEvenRankCountWithAtLeastFourRanks
{
public:
    //! Function to require even ranks with at least four ranks
    static bool conditionSatisfied(const int numRanks)
    {
        return (numRanks > 2) && (numRanks % 2 == 0);
    }
    //! Text to echo when skipping a test that does not satisfy the requirement
    inline static const char* s_skipReason = "an even rank count of at least four is required";
};

TEST(ThreadAffinityMultiRankTest, PinsWholeNode)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    ThreadAffinityTestHelper helper;
    helper.setLogicalProcessorCount(getNumberOfTestMpiRanks());
    helper.expectPinningMessage(false, 1);
    helper.expectAffinitySet(gmx_node_rank());
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, PinsWithOffsetAndStride)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    helper.setOffsetAndStride(1, 2);
    helper.setLogicalProcessorCount(2 * getNumberOfTestMpiRanks());
    helper.expectWarningMatchingRegex("Applying core pinning offset 1");
    helper.expectPinningMessage(true, 2);
    helper.expectAffinitySet(1 + 2 * gmx_node_rank());
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, PinsTwoNodes)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setPhysicalNodeId(gmx_node_rank() / 2);
    helper.setLogicalProcessorCount(2);
    helper.expectPinningMessage(false, 1);
    helper.expectAffinitySet(gmx_node_rank() % 2);
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, DoesNothingWhenDisabled)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::Off);
    helper.setLogicalProcessorCount(getNumberOfTestMpiRanks());
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, HandlesTooManyThreadsWithAuto)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    ThreadAffinityTestHelper helper;
    const int                threadsPerRank = 2;
    helper.setLogicalProcessorCount(threadsPerRank * getNumberOfTestMpiRanks() - 1);
    helper.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper.setAffinity(threadsPerRank);
}

TEST(ThreadAffinityMultiRankTest, HandlesTooManyThreadsWithForce)
{
    GMX_MPI_TEST(AllowAnyRankCount);
    ThreadAffinityTestHelper helper;
    const int                threadsPerRank = 2;
    helper.setAffinityOption(ThreadAffinity::On);
    helper.setLogicalProcessorCount(threadsPerRank * getNumberOfTestMpiRanks() - 1);
    helper.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper.setAffinity(threadsPerRank);
}

class ThreadAffinityHeterogeneousNodesTest : public ::testing::Test
{
public:
    static int  currentNode() { return gmx_node_rank() / 2; }
    static int  indexInNode() { return gmx_node_rank() % 2; }
    static bool isMain() { return gmx_node_rank() == 0; }

    static void setupNodes(ThreadAffinityTestHelper* helper, int coresOnNodeZero, int coresOnOtherNodes)
    {
        const int node = currentNode();
        helper->setPhysicalNodeId(node);
        helper->setLogicalProcessorCount(node == 0 ? coresOnNodeZero : coresOnOtherNodes);
    }
    static void expectNodeAffinitySet(ThreadAffinityTestHelper* helper, int node, int core)
    {
        if (currentNode() == node)
        {
            helper->expectAffinitySet(core);
        }
    }
};

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsOnMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    setupNodes(&helper, 2, 1);
    helper.expectWarningMatchingRegexIf("Oversubscribing available/permitted CPUs",
                                        isMain() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsOnNonMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    setupNodes(&helper, 1, 2);
    helper.expectWarningMatchingRegexIf("Oversubscribing available/permitted CPUs", currentNode() == 0);
    if (currentNode() >= 1)
    {
        helper.expectPinningMessage(false, 1);
        expectNodeAffinitySet(&helper, currentNode(), indexInNode());
    }
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesUnknownHardwareOnNonMain)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    setupNodes(&helper, 2, 0);
    helper.expectWarningMatchingRegexIf("No information on available logical cpus",
                                        isMain() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsAutomaticallyOnMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    setupNodes(&helper, 2, 1);
    helper.expectWarningMatchingRegexIf("Oversubscribing available/permitted CPUs",
                                        isMain() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsAutomaticallyOnNonMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    setupNodes(&helper, 1, 2);
    helper.expectWarningMatchingRegexIf("Oversubscribing available/permitted CPUs", currentNode() == 0);
    if (currentNode() >= 1)
    {
        helper.expectPinningMessage(false, 1);
        expectNodeAffinitySet(&helper, currentNode(), indexInNode());
    }
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesInvalidOffsetOnNonMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    helper.setOffsetAndStride(2, 0);
    setupNodes(&helper, 4, 2);
    helper.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper.expectWarningMatchingRegexIf("Requested offset too large", isMain() || currentNode() >= 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode() + 2);
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesInvalidStrideOnNonMainOnly)
{
    GMX_MPI_TEST(RequireEvenRankCountWithAtLeastFourRanks);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(ThreadAffinity::On);
    helper.setOffsetAndStride(0, 2);
    setupNodes(&helper, 4, 2);
    helper.expectWarningMatchingRegexIf("Requested stride too large", isMain() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(true, 2);
    }
    expectNodeAffinitySet(&helper, 0, 2 * indexInNode());
    helper.setAffinity(1);
}

} // namespace
} // namespace test
} // namespace gmx
