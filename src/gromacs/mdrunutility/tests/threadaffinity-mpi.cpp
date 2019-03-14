/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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

#include <array>

#include <gtest/gtest.h>

#include "gromacs/utility/basenetwork.h"

#include "testutils/mpitest.h"

#include "threadaffinitytest.h"

namespace
{

using gmx::test::ThreadAffinityTestHelper;

TEST(ThreadAffinityMultiRankTest, PinsWholeNode)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setLogicalProcessorCount(4);
    helper.expectPinningMessage(false, 1);
    helper.expectAffinitySet(gmx_node_rank());
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, PinsWithOffsetAndStride)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    helper.setOffsetAndStride(1, 2);
    helper.setLogicalProcessorCount(8);
    helper.expectWarningMatchingRegex("Applying core pinning offset 1");
    helper.expectPinningMessage(true, 2);
    helper.expectAffinitySet(1 + 2*gmx_node_rank());
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, PinsTwoNodes)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setPhysicalNodeId(gmx_node_rank()/2);
    helper.setLogicalProcessorCount(2);
    helper.expectPinningMessage(false, 1);
    helper.expectAffinitySet(gmx_node_rank()%2);
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, DoesNothingWhenDisabled)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffOFF);
    helper.setLogicalProcessorCount(4);
    helper.setAffinity(1);
}

TEST(ThreadAffinityMultiRankTest, HandlesTooManyThreadsWithAuto)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setLogicalProcessorCount(6);
    helper.expectWarningMatchingRegex("Oversubscribing the CPU");
    helper.setAffinity(2);
}

TEST(ThreadAffinityMultiRankTest, HandlesTooManyThreadsWithForce)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    helper.setLogicalProcessorCount(6);
    helper.expectWarningMatchingRegex("Oversubscribing the CPU");
    helper.setAffinity(2);
}

class ThreadAffinityHeterogeneousNodesTest : public ::testing::Test
{
    public:
        int currentNode() const { return gmx_node_rank() / 2; }
        int indexInNode() const { return gmx_node_rank() % 2; }
        bool isMaster() const { return gmx_node_rank() == 0; }

        void setupNodes(ThreadAffinityTestHelper *helper, std::array<int, 2> cores)
        {
            const int node = currentNode();
            helper->setPhysicalNodeId(node);
            helper->setLogicalProcessorCount(cores[node]);
        }
        void expectNodeAffinitySet(ThreadAffinityTestHelper *helper, int node, int core)
        {
            if (currentNode() == node)
            {
                helper->expectAffinitySet(core);
            }
        }
};

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsOnMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    setupNodes(&helper, {{2, 1}});
    helper.expectWarningMatchingRegexIf("Oversubscribing the CPU", isMaster() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsOnNonMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    setupNodes(&helper, {{1, 2}});
    helper.expectWarningMatchingRegexIf("Oversubscribing the CPU", currentNode() == 0);
    if (currentNode() == 1)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 1, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesUnknownHardwareOnNonMaster)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    setupNodes(&helper, {{2, 0}});
    helper.expectWarningMatchingRegexIf("No information on available cores", isMaster() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsAutomaticallyOnMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    setupNodes(&helper, {{2, 1}});
    helper.expectWarningMatchingRegexIf("Oversubscribing the CPU", isMaster() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, PinsAutomaticallyOnNonMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    setupNodes(&helper, {{1, 2}});
    helper.expectWarningMatchingRegexIf("Oversubscribing the CPU", currentNode() == 0);
    if (currentNode() == 1)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 1, indexInNode());
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesInvalidOffsetOnNonMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    helper.setOffsetAndStride(2, 0);
    setupNodes(&helper, {{4, 2}});
    helper.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper.expectWarningMatchingRegexIf("Requested offset too large", isMaster() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(false, 1);
    }
    expectNodeAffinitySet(&helper, 0, indexInNode()+2);
    helper.setAffinity(1);
}

TEST_F(ThreadAffinityHeterogeneousNodesTest, HandlesInvalidStrideOnNonMasterOnly)
{
    GMX_MPI_TEST(4);
    ThreadAffinityTestHelper helper;
    helper.setAffinityOption(threadaffON);
    helper.setOffsetAndStride(0, 2);
    setupNodes(&helper, {{4, 2}});
    helper.expectWarningMatchingRegexIf("Requested stride too large", isMaster() || currentNode() == 1);
    if (currentNode() == 0)
    {
        helper.expectPinningMessage(true, 2);
    }
    expectNodeAffinitySet(&helper, 0, 2*indexInNode());
    helper.setAffinity(1);
}

} // namespace
