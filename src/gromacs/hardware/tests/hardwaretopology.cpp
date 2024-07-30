/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for gmx::HardwareTopology
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_hardware
 */
#include "gmxpre.h"

#include "gromacs/hardware/hardwaretopology.h"

#include "config.h"

#include <cstddef>

#include <algorithm>
#include <array>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

// Simple tests that execute actual hardware detection. We cannot compare
// the exact output to any reference, since it is hardware-dependent.

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
            << "Please make a post at the GROMACS development forum at" << std::endl
            << "https://gromacs.bioexcel.eu/c/gromacs-developers/10 so we can try to fix it.";
}

#if GMX_USE_HWLOC
TEST(HardwareTopologyTest, HwlocExecute)
{
#    if defined(__linux__)
    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    // On Linux with hwloc support we should be able to get at least basic information
    EXPECT_GE(hwTop.supportLevel(), gmx::HardwareTopology::SupportLevel::Basic)
            << "Cannot determine basic hardware topology from hwloc. GROMACS will still" << std::endl
            << "work, your performance on large numbers of nodes might be affected." << std::endl
            << "Please make a post at the GROMACS development forum at" << std::endl
            << "https://gromacs.bioexcel.eu/c/gromacs-developers/10 so we can try to fix it.";
#    endif
}
#endif

TEST(HardwareTopologyTest, ProcessorSelfconsistency)
{
    gmx::HardwareTopology hwTop(gmx::HardwareTopology::detect());

    if (hwTop.supportLevel() >= gmx::HardwareTopology::SupportLevel::Basic)
    {
        SCOPED_TRACE(gmx::formatString("max threads from hardware topology: %d", hwTop.maxThreads()));

        int packagesInTopology = hwTop.machine().packages.size();

        auto logicalProcessors = hwTop.machine().logicalProcessors;
        for (auto logicalProcessorIt = logicalProcessors.begin();
             logicalProcessorIt != logicalProcessors.end();
             ++logicalProcessorIt)
        {
            int packageRankInTopology = logicalProcessorIt->packageRankInTopology;
            int coreRankInPackage     = logicalProcessorIt->coreRankInPackage;
            int puRankInCore          = logicalProcessorIt->processingUnitRankInCore;

            int coresInPackage = hwTop.machine().packages[packageRankInTopology].cores.size();
            int pusInCore      = hwTop.machine()
                                    .packages[packageRankInTopology]
                                    .cores[coreRankInPackage]
                                    .processingUnits.size();

            // Check that logical processor information contains
            // reasonable values.
            SCOPED_TRACE(gmx::formatString("Socket rank in topology: %d",
                                           logicalProcessorIt->packageRankInTopology));
            SCOPED_TRACE(gmx::formatString("Core rank in socket:    %d",
                                           logicalProcessorIt->coreRankInPackage));
            SCOPED_TRACE(gmx::formatString("PU rank in core: %d",
                                           logicalProcessorIt->processingUnitRankInCore));
            EXPECT_TRUE(packageRankInTopology >= 0 && packageRankInTopology < packagesInTopology);
            EXPECT_TRUE(coreRankInPackage >= 0 && coreRankInPackage < coresInPackage);
            EXPECT_TRUE(puRankInCore >= 0 && puRankInCore < pusInCore);
            // Check that logical processor information is distinct
            // for each logical processor.

            for (auto remainingLogicalProcessorIt = logicalProcessorIt + 1;
                 remainingLogicalProcessorIt != logicalProcessors.end();
                 ++remainingLogicalProcessorIt)
            {
                SCOPED_TRACE(gmx::formatString("Other socket rank in topology: %d",
                                               remainingLogicalProcessorIt->packageRankInTopology));
                SCOPED_TRACE(gmx::formatString("Other core rank in socket:    %d",
                                               remainingLogicalProcessorIt->coreRankInPackage));
                SCOPED_TRACE(gmx::formatString("Other hw thread rank in core: %d",
                                               remainingLogicalProcessorIt->processingUnitRankInCore));
                EXPECT_TRUE((logicalProcessorIt->packageRankInTopology
                             != remainingLogicalProcessorIt->packageRankInTopology)
                            || (logicalProcessorIt->coreRankInPackage != remainingLogicalProcessorIt->coreRankInPackage)
                            || (logicalProcessorIt->processingUnitRankInCore
                                != remainingLogicalProcessorIt->processingUnitRankInCore))
                        << "This pair of logical processors have the same descriptive information, "
                           "which is an error";
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
        for (const auto& n : hwTop.machine().numa.nodes)
        {
            processorsinNumaNudes += n.processingUnits.size();
        }
        EXPECT_EQ(processorsinNumaNudes, hwTop.machine().logicalProcessors.size());

        // Check that every processor is in a numa domain (i.e., that they are unique)
        std::vector<int> v(hwTop.machine().logicalProcessors.size());
        for (auto& elem : v)
        {
            elem = 0;
        }
        for (const auto& n : hwTop.machine().numa.nodes)
        {
            for (const auto& idx : n.processingUnits)
            {
                v[idx] = 1;
            }
        }
        int uniqueProcessorsinNumaNudes = std::count(v.begin(), v.end(), 1);
        EXPECT_EQ(uniqueProcessorsinNumaNudes, hwTop.machine().logicalProcessors.size());

        // We must have some memory in a numa node
        for (const auto& n : hwTop.machine().numa.nodes)
        {
            EXPECT_GT(n.memory, 0);
        }

        // Check latency matrix size and contents
        EXPECT_GT(hwTop.machine().numa.baseLatency, 0);
        EXPECT_GT(hwTop.machine().numa.maxRelativeLatency, 0);
        // Check number of rows matches # numa nodes
        EXPECT_EQ(hwTop.machine().numa.relativeLatency.size(), hwTop.machine().numa.nodes.size());
        for (const auto& v2 : hwTop.machine().numa.relativeLatency)
        {
            // Check that size of each row matches # numa nodes
            EXPECT_EQ(v2.size(), hwTop.machine().numa.nodes.size());
            for (const auto& latency : v2)
            {
                // Latency values should be positive
                EXPECT_GT(latency, 0);
            }
        }

        // We don't check cache fields because these tests depend both
        // on whether hwloc can detect things correctly, and then
        // whether GROMACS code packages the results correctly. The
        // hwloc cache detection is fragile and can report 0 for cache
        // size, line size or associativity (=unknown), so GROMACS
        // doesn't test anything related to it.
        //
        // TODO Use proper unit tests on mock hardware to test that
        // HardwareTopology construction is doing its job, rather than
        // brittle tests that require that hwloc works correctly on
        // the user's hardware even when GROMACS is barely using the
        // values returned by hwloc.
    }
}

} // namespace
} // namespace test
} // namespace gmx
