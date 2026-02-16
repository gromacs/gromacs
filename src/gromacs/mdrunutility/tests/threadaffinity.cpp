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

#include "gromacs/mdrunutility/threadaffinity.h"

#include "config.h"

#include <numeric>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/hardware/hw_info.h"

#include "testutils/testfilemanager.h"

#include "threadaffinitytest.h"

namespace gmx
{
namespace test
{
namespace
{

class ThreadAffinityTest : public ::testing::Test
{
public:
    gmx::test::ThreadAffinityTestHelper helper_;
};

TEST_F(ThreadAffinityTest, DoesNothingWhenDisabled)
{
    helper_.setAffinityOption(ThreadAffinity::Off);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWhenNotSupported)
{
    helper_.setAffinitySupported(false);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooFewUserSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("The number of threads is not equal to the number of");
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooManyUserSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithAutoAndTooManyAutoSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.setTotNumThreadsIsAuto(true);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithUnknownHardware)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(0);
    helper_.expectWarningMatchingRegex("No information on available logical cpus");
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooManyThreads)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Oversubscribing available/permitted CPUs");
    helper_.setAffinity(8);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeOffset)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(2, 0);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper_.expectWarningMatchingRegex("Requested offset too large");
    helper_.setAffinity(3);
}

TEST_F(ThreadAffinityTest, DoesNothingWithTooLargeStride)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(0, 2);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Requested stride too large");
    helper_.setAffinity(3);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithAuto)
{
    helper_.setLogicalProcessorCount(1);
    helper_.expectAffinitySet(0);
    helper_.expectPinningMessage(false, 1);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet(0);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsSingleThreadWithOffsetWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(2, 0);
    helper_.setLogicalProcessorCount(4);
    helper_.expectWarningMatchingRegex("Applying core pinning offset 2");
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet(2);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithSingleThread)
{
    helper_.setLogicalProcessorCount(1);
    helper_.expectPinningMessage(false, 1);
    // It would be good to check for the warning, but that currently goes to stderr
    helper_.expectAffinitySetThatFails(0);
    helper_.setAffinity(1);
}


TEST_F(ThreadAffinityTest, PinsWithAffinityAndShiftedMask)
{
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    helper_.setLogicalProcessorCount(2);
    helper_.setExternalAffinitySet({ 1 });
    helper_.expectAffinitySet(1);
    helper_.setAffinity(1);
}

// TODO: If it wouldn't result in a multitude of #if's, it would be nice
// to somehow indicate in a no-OpenMP build that some tests are missing.
#if GMX_OPENMP
TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithAuto)
{
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 1);
    helper_.expectAffinitySet({ 0, 1 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsMultipleThreadsWithStrideWhenForced)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setOffsetAndStride(0, 2);
    helper_.setLogicalProcessorCount(4);
    helper_.expectPinningMessage(true, 2);
    helper_.expectAffinitySet({ 0, 2 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsWithAutoAndFewerAutoSetThreads)
{
    helper_.setLogicalProcessorCount(4);
    helper_.setTotNumThreadsIsAuto(true);
    helper_.expectPinningMessage(false, 2);
    helper_.expectAffinitySet({ 0, 2 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, HandlesPinningFailureWithOneThreadFailing)
{
    helper_.setAffinityOption(ThreadAffinity::On);
    helper_.setLogicalProcessorCount(2);
    helper_.expectPinningMessage(false, 1);
    helper_.expectWarningMatchingRegex("Affinity setting for 1/2 threads failed");
    helper_.expectGenericFailureMessage();
    helper_.expectAffinitySet(0);
    helper_.expectAffinitySetThatFails(1);
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsWithAffinityAndFullMask)
{
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    helper_.setLogicalProcessorCount(2);
    helper_.setExternalAffinitySet({ 0, 1 });
    helper_.expectAffinitySet(0);
    helper_.setAffinity(1);
}

TEST_F(ThreadAffinityTest, PinsWithAffinityAndHoleyMask)
{
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    helper_.setLogicalProcessorCount(4);
    helper_.setExternalAffinitySet({ 1, 3 });
    helper_.expectAffinitySet({ 1, 3 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, RefusesWhenAffinityMaskTooSmall)
{
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    helper_.setLogicalProcessorCount(4);
    helper_.setExternalAffinitySet({ 1 });
    helper_.expectWarningMatchingRegex(
            "The number of threads on one rank is greater than .* in the process affinity");
    helper_.setAffinity(2);
}

// Borrow test data from hardware tests to get a more realistic topology
TEST_F(ThreadAffinityTest, PinsToHalfSocketWithSmt)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 2 sockets, 8 cores each, 2 threads each
    std::vector<int> allowedCpus(32);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 4, 5, 6, 7, 20, 21, 22, 23 });
    helper_.setTopologyFromSavedMock(hwTestFile.getInputFilePath("XeonE52620v4_Cgroups2NoLimit").string(),
                                     allowedCpus,
                                     externalAffinitySet);
    helper_.expectAffinitySet({ 4, 5 }); // OS index, cores 4-5, first hw thread
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsToHalfSocketWithNoSmt)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 2 sockets, 8 cores each, 2 threads each
    std::vector<int> allowedCpus(32);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 4, 5, 6, 7 });
    helper_.setTopologyFromSavedMock(hwTestFile.getInputFilePath("XeonE52620v4_Cgroups2NoLimit").string(),
                                     allowedCpus,
                                     externalAffinitySet);
    helper_.expectAffinitySet({ 4, 5 });
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsToPCoresInHybridSystem)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 8 P-cores with 2 threads each, 8 E-cores with 1 thread each
    std::vector<int> allowedCpus(24);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 0, 1, 2, 3 }); // Both HW threads on P-cores #0-1
    helper_.setTopologyFromSavedMock(
            hwTestFile.getInputFilePath("Core12900K_Cgroups2CpuLimit650pct").string(),
            allowedCpus,
            externalAffinitySet);
    helper_.expectAffinitySet({ 0, 2 }); // First HW threads on P-cores #0 and #1
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsToSecondHwThreadsOnPCoresInHybridSystem)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 8 P-cores with 2 threads each, 8 E-cores with 1 thread each
    std::vector<int> allowedCpus(24);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 1, 3, 5, 7 }); // Second HW threads on P-cores #0-4
    helper_.setTopologyFromSavedMock(
            hwTestFile.getInputFilePath("Core12900K_Cgroups2CpuLimit650pct").string(),
            allowedCpus,
            externalAffinitySet);
    helper_.expectAffinitySet({ 1, 3 }); // Second HW threads on P-cores #0 and #1
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsToSinglePCoreInHybridSystem)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 8 P-cores with 2 threads each, 8 E-cores with 1 thread each
    std::vector<int> allowedCpus(24);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 0, 1 }); // Both HW threads on P-core #0
    helper_.setTopologyFromSavedMock(
            hwTestFile.getInputFilePath("Core12900K_Cgroups2CpuLimit650pct").string(),
            allowedCpus,
            externalAffinitySet);
    helper_.expectAffinitySet({ 0, 1 }); // Both HW threads on P-core #0
    helper_.setAffinity(2);
}

TEST_F(ThreadAffinityTest, PinsToECoresInHybridSystem)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 8 P-cores with 2 threads each, 8 E-cores with 1 thread each
    std::vector<int> allowedCpus(24);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 20, 21, 22, 23 }); // E-cores #4-7
    helper_.setTopologyFromSavedMock(
            hwTestFile.getInputFilePath("Core12900K_Cgroups2CpuLimit650pct").string(),
            allowedCpus,
            externalAffinitySet);
    helper_.expectAffinitySet({ 20, 21 }); // E-cores #4 and #5
    helper_.setAffinity(2);
}


TEST_F(ThreadAffinityTest, PinsToMixOfPAndECoresInHybridSystem)
{
    TestFileManager hwTestFile;
    hwTestFile.setInputDataDirectory(hwTestFile.getInputDataDirectory().parent_path().parent_path()
                                     / "hardware" / "tests");
    helper_.setAffinityOption(ThreadAffinity::Inherit);
    // 8 P-cores with 2 threads each, 8 E-cores with 1 thread each
    std::vector<int> allowedCpus(24);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);
    std::vector<int> externalAffinitySet({ 2, 3, 16, 17 }); // P-core #1 and E-cores #0-1
    helper_.setTopologyFromSavedMock(
            hwTestFile.getInputFilePath("Core12900K_Cgroups2CpuLimit650pct").string(),
            allowedCpus,
            externalAffinitySet);
    helper_.expectAffinitySet({ 2, 16 }); // First HW thread of P-core #1 and E-core #0
    helper_.setAffinity(2);
}

#endif

} // namespace
} // namespace test
} // namespace gmx
