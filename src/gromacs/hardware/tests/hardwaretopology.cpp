/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Tests for gmx::HardwareTopology
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/hardware/hardwaretopology.h"

#include "config.h"

#include <algorithm>

#include <gtest/gtest.h>

namespace
{

// There is no way we can compare to any reference data since that
// depends on the architecture, but we can at least make sure that it
// works to execute the tests and that they are self-consistent

// Although it is not strictly an error, for the very basic execution tests
// we also report if we cannot extract the hardware topology on systems
// where we expect to be able to. Since this might happen to users, we
// provide a bit more information and ask them to mail us in this case.

TEST(HardwareTopologyTest, Execute)
{
    // There is no way we can compare to any reference data since that
    // depends on the architecture, but we can at least make sure that it
    // works to execute the tests

    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    // If we cannot even find the number of logical processors we want to flag it
    EXPECT_GT(hwTop.supportLevel(), gmx::HardwareTopology::SupportLevel::None)
    << "Cannot determine number of processors. " << std::endl
    << "GROMACS might still work, but it will likely hurt your performance." << std::endl
    << "Please mail gmx-developers@gromacs.org so we can try to fix it.";
}

#if GMX_HWLOC
TEST(HardwareTopologyTest, HwlocExecute)
{
#if defined(__linux__)
    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    // On Linux with hwloc support we should be able to get at least basic information
    EXPECT_GE(hwTop.supportLevel(), gmx::HardwareTopology::SupportLevel::Basic)
    << "Cannot determine basic hardware topology from hwloc. GROMACS will still\n" << std::endl
    << "work, but it might affect your performance for large nodes." << std::endl
    << "Please mail gmx-developers@gromacs.org so we can try to fix it.";
#endif
}
#endif

TEST(HardwareTopologyTest, ProcessorSelfconsistency)
{
    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
    {
        int socketsInMachine = hwTop.machine().sockets.size();
        int coresPerSocket   = hwTop.machine().sockets[0].cores.size();
        int hwThreadsPerCore = hwTop.machine().sockets[0].cores[0].hwThreads.size();

        // Check that logical processor information is reasonable
        for (auto &l : hwTop.machine().logicalProcessors)
        {
            EXPECT_TRUE(l.socketRankInMachine >= 0 && l.socketRankInMachine < socketsInMachine);
            EXPECT_TRUE(l.coreRankInSocket >= 0 && l.coreRankInSocket < coresPerSocket);
            EXPECT_TRUE(l.hwThreadRankInCore >= 0 && l.hwThreadRankInCore < hwThreadsPerCore);
        }

        // Double-check that the tree is self-consistent with logical processor info
        for (int s = 0; s < socketsInMachine; s++)
        {
            for (int c = 0; c < coresPerSocket; c++)
            {
                for (int t = 0; t < hwThreadsPerCore; t++)
                {
                    int idx = hwTop.machine().sockets[s].cores[c].hwThreads[t].logicalProcessorId;
                    EXPECT_LT(idx, hwTop.machine().logicalProcessorCount);
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].socketRankInMachine, s);
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].coreRankInSocket, c) << "logical:" << idx;
                    EXPECT_EQ(hwTop.machine().logicalProcessors[idx].hwThreadRankInCore, t);
                }
            }
        }
    }
}

TEST(HardwareTopologyTest, NumaCacheSelfconsistency)
{
    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Full)
    {
        // Check that numa node id corresponds to rank
        for (std::size_t i = 0; i < hwTop.machine().numa.nodes.size(); i++)
        {
            EXPECT_EQ(hwTop.machine().numa.nodes[i].id, i);
        }

        // Check that the sum of numa domains is the total processor count
        int processorsinNumaNudes = 0;
        for (auto &n : hwTop.machine().numa.nodes)
        {
            processorsinNumaNudes += n.logicalProcessorId.size();
        }
        EXPECT_EQ(processorsinNumaNudes, hwTop.machine().logicalProcessorCount);

        // Check that every processor is in a numa domain (i.e., that they are unique)
        std::vector<int> v(hwTop.machine().logicalProcessorCount);
        for (auto &elem : v)
        {
            elem = 0;
        }
        for (auto &n : hwTop.machine().numa.nodes)
        {
            for (auto &idx : n.logicalProcessorId)
            {
                v[idx] = 1;
            }
        }
        int uniqueProcessorsinNumaNudes = std::count(v.begin(), v.end(), 1);
        EXPECT_EQ(uniqueProcessorsinNumaNudes, hwTop.machine().logicalProcessorCount);

        // We must have some memory in a numa node
        for (auto &n : hwTop.machine().numa.nodes)
        {
            EXPECT_GT(n.memory, 0);
        }

        // Check latency matrix size and contents
        EXPECT_GT(hwTop.machine().numa.baseLatency, 0);
        EXPECT_GT(hwTop.machine().numa.maxRelativeLatency, 0);
        // Check number of rows matches # numa nodes
        EXPECT_EQ(hwTop.machine().numa.relativeLatency.size(), hwTop.machine().numa.nodes.size());
        for (auto &v2 : hwTop.machine().numa.relativeLatency)
        {
            // Check that size of each row matches # numa nodes
            EXPECT_EQ(v2.size(), hwTop.machine().numa.nodes.size());
            for (auto &latency : v2)
            {
                // Latency values should be positive
                EXPECT_GT(latency, 0);
            }
        }

        // Check cache. The hwloc cache detection is fragile and can report
        // 0 for line size or associativity (=unknown), so we just check the size.
        for (auto &c : hwTop.machine().caches)
        {
            EXPECT_GT(c.size, 0);
        }
    }
}


} // namespace
