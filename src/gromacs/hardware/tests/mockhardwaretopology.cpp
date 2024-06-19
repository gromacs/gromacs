/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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

#include "config.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <iterator>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include "gromacs/hardware/hardwaretopology.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

// Utility functions to check parts of the hardware topology data structure in detail.

/*! \brief Check if fields in a LogicalProcessor struct match
 *
 * \param checker  Initialized checker
 * \param cpu      Initialized LogicalProcessor structure in hardwaretopology.
 *
 * \note This routine is typically called from checkMachine().
 */
void checkLogicalProcessor(TestReferenceChecker* checker, const HardwareTopology::LogicalProcessor& cpu)
{
    TestReferenceChecker compound(checker->checkCompound("LogicalProcessorInfo", nullptr));
    compound.checkInteger(cpu.puId, "puId");
    compound.checkInteger(cpu.osId, "osId");
    compound.checkInteger(cpu.packageRankInTopology, "packageRankInTopology");
    compound.checkInteger(cpu.coreRankInPackage, "coreRankInPackage");
    compound.checkInteger(cpu.processingUnitRankInCore, "processingUnitRankInCore");
    // We don't test numa node id, since it might differ between mocking and real tests
}

/*! \brief Check if fields in a ProcessingUnit struct match
 *
 * \param checker  Initialized checker
 * \param pu       Initialized processing unit structure in hardwaretopology.
 *
 * \note This routine is typically called from checkCore().
 */
void checkProcessingUnit(TestReferenceChecker* checker, const HardwareTopology::ProcessingUnit& pu)
{
    TestReferenceChecker compound(checker->checkCompound("ProcessingUnit", nullptr));
    compound.checkInteger(pu.id, "id");
    compound.checkInteger(pu.osId, "osId");
}

/*! \brief Check if fields in a Core struct match
 *
 * \param checker  Initialized checker
 * \param core     Initialized Core structure in hardwaretopology.
 *
 * \note This routine is typically called from checkPackage().
 */
void checkCore(TestReferenceChecker* checker, const HardwareTopology::Core& core)
{
    TestReferenceChecker compound(checker->checkCompound("Core", nullptr));
    compound.checkInteger(core.id, "id");
    compound.checkSequence(
            core.processingUnits.begin(), core.processingUnits.end(), "processingUnits", checkProcessingUnit);
}

/*! \brief Check if fields in a Package struct match
 *
 * \param checker  Initialized checker
 * \param package  Initialized Package structure in hardwaretopology.
 *
 * \note This routine is typically called from checkMachine().
 */
void checkPackage(TestReferenceChecker* checker, const HardwareTopology::Package& package)
{
    TestReferenceChecker compound(checker->checkCompound("Package", nullptr));
    compound.checkInteger(package.id, "id");
    compound.checkSequence(package.cores.begin(), package.cores.end(), "processingUnits", checkCore);
}

/*! \brief Check if an entry in the osIdToPuId map matches
 *
 * \param checker  Initialized checker
 * \param mapEntry Reference to single entry in osId to puId map
 *
 * \note This routine is typically called from checkMachine().
 */
void checkMap(TestReferenceChecker* checker, const std::pair<int, int>& mapEntry)
{
    TestReferenceChecker compound(checker->checkCompound("MapEntry", nullptr));
    compound.checkInteger(mapEntry.first, "osId");
    compound.checkInteger(mapEntry.second, "puId");
}

/*! \brief Check if entire machine structure in a hardware topology matches
 *
 * \param checker  Initialized checker
 * \param machine  Initialized Machine structure from hardware topology
 */
void checkMachine(TestReferenceChecker* checker, const HardwareTopology::Machine& machine)
{
    TestReferenceChecker compound(checker->checkCompound("Machine", nullptr));
    compound.checkSequence(machine.logicalProcessors.begin(),
                           machine.logicalProcessors.end(),
                           "LogicalProcessors",
                           checkLogicalProcessor);
    compound.checkSequence(machine.osIdToPuId.begin(), machine.osIdToPuId.end(), "OsIdToPuIdMap", checkMap);
    compound.checkSequence(machine.packages.begin(), machine.packages.end(), "Packages", checkPackage);
    // We don't check caches, numa, or devices since we can't mock them.
}

void checkHardwareTopology(TestReferenceChecker* checker, const HardwareTopology& hwTop)
{
    TestReferenceChecker compound(checker->checkCompound("HardwareTopology", nullptr));
    checkMachine(&compound, hwTop.machine());
    compound.checkInteger(static_cast<int>(hwTop.supportLevel()), "SupportLevel");
    compound.checkBoolean(hwTop.isThisSystem(), "IsThisSystem");
    compound.checkFloat(hwTop.cpuLimit(), "CpuLimit");
    compound.checkInteger(hwTop.maxThreads(), "MaxThreads");
}

//! Type name for processor map
using ProcessorMap = std::map<int, std::array<int, 3>>;

//! Test parameter input struct
using HardwareDetectionTestParams = std::tuple<ProcessorMap, std::vector<int>, std::string>;


//! Hardware map for XeonE52620v4
const ProcessorMap logicalProcessorIdMapXeonE52620v4 = {
    { 0, { 0, 0, 0 } },  { 1, { 0, 1, 0 } },  { 2, { 0, 2, 0 } },  { 3, { 0, 3, 0 } },
    { 4, { 0, 4, 0 } },  { 5, { 0, 5, 0 } },  { 6, { 0, 6, 0 } },  { 7, { 0, 7, 0 } },
    { 8, { 1, 0, 0 } },  { 9, { 1, 1, 0 } },  { 10, { 1, 2, 0 } }, { 11, { 1, 3, 0 } },
    { 12, { 1, 4, 0 } }, { 13, { 1, 5, 0 } }, { 14, { 1, 6, 0 } }, { 15, { 1, 7, 0 } },
    { 16, { 0, 0, 1 } }, { 17, { 0, 1, 1 } }, { 18, { 0, 2, 1 } }, { 19, { 0, 3, 1 } },
    { 20, { 0, 4, 1 } }, { 21, { 0, 5, 1 } }, { 22, { 0, 6, 1 } }, { 23, { 0, 7, 1 } },
    { 24, { 1, 0, 1 } }, { 25, { 1, 1, 1 } }, { 26, { 1, 2, 1 } }, { 27, { 1, 3, 1 } },
    { 28, { 1, 4, 1 } }, { 29, { 1, 5, 1 } }, { 30, { 1, 6, 1 } }, { 31, { 1, 7, 1 } },
};

//! Hardware map for Xeon4116
const ProcessorMap logicalProcessorIdMapXeon4116 = {
    { 0, { 0, 0, 0 } },   { 1, { 0, 1, 0 } },   { 2, { 0, 2, 0 } },   { 3, { 0, 3, 0 } },
    { 4, { 0, 4, 0 } },   { 5, { 0, 5, 0 } },   { 6, { 0, 8, 0 } },   { 7, { 0, 9, 0 } },
    { 8, { 0, 10, 0 } },  { 9, { 0, 11, 0 } },  { 10, { 0, 12, 0 } }, { 11, { 0, 13, 0 } },
    { 12, { 1, 0, 0 } },  { 13, { 1, 1, 0 } },  { 14, { 1, 2, 0 } },  { 15, { 1, 3, 0 } },
    { 16, { 1, 4, 0 } },  { 17, { 1, 5, 0 } },  { 18, { 1, 8, 0 } },  { 19, { 1, 9, 0 } },
    { 20, { 1, 10, 0 } }, { 21, { 1, 11, 0 } }, { 22, { 1, 12, 0 } }, { 23, { 1, 13, 0 } },
    { 24, { 0, 0, 1 } },  { 25, { 0, 1, 1 } },  { 26, { 0, 2, 1 } },  { 27, { 0, 3, 1 } },
    { 28, { 0, 4, 1 } },  { 29, { 0, 5, 1 } },  { 30, { 0, 8, 1 } },  { 31, { 0, 9, 1 } },
    { 32, { 0, 10, 1 } }, { 33, { 0, 11, 1 } }, { 34, { 0, 12, 1 } }, { 35, { 0, 13, 1 } },
    { 36, { 1, 0, 1 } },  { 37, { 1, 1, 1 } },  { 38, { 1, 2, 1 } },  { 39, { 1, 3, 1 } },
    { 40, { 1, 4, 1 } },  { 41, { 1, 5, 1 } },  { 42, { 1, 8, 1 } },  { 43, { 1, 9, 1 } },
    { 44, { 1, 10, 1 } }, { 45, { 1, 11, 1 } }, { 46, { 1, 12, 1 } }, { 47, { 1, 13, 1 } },
};

//! Hardware map for Core12900K
const ProcessorMap logicalProcessorIdMapCore12900K = {
    { 0, { 0, 0, 0 } },   { 1, { 0, 0, 1 } },   { 2, { 0, 4, 0 } },   { 3, { 0, 4, 1 } },
    { 4, { 0, 8, 0 } },   { 5, { 0, 8, 1 } },   { 6, { 0, 12, 0 } },  { 7, { 0, 12, 1 } },
    { 8, { 0, 16, 0 } },  { 9, { 0, 16, 1 } },  { 10, { 0, 20, 0 } }, { 11, { 0, 20, 1 } },
    { 12, { 0, 24, 0 } }, { 13, { 0, 24, 1 } }, { 14, { 0, 28, 0 } }, { 15, { 0, 28, 1 } },
    { 16, { 0, 32, 0 } }, { 17, { 0, 33, 0 } }, { 18, { 0, 34, 0 } }, { 19, { 0, 35, 0 } },
    { 20, { 0, 36, 0 } }, { 21, { 0, 37, 0 } }, { 22, { 0, 38, 0 } }, { 23, { 0, 39, 0 } },
};

//! Empty allowed CPU vector
const std::vector<int> emptyCPUVector = {};
//! Slurm slot 1 of 4 XeonE52620v4
const std::vector<int> slurm1of4XeonE52620v4 = { 0, 1, 2, 3, 16, 17, 18, 19 };
//! Slurm slot 2 of 4 XeonE52620v4
const std::vector<int> slurm2of4XeonE52620v4 = { 4, 5, 6, 7, 20, 21, 22, 23 };
//! Slurm slot 3 of 4 XeonE52620v4
const std::vector<int> slurm3of4XeonE52620v4 = { 8, 9, 10, 11, 24, 25, 26, 27 };
//! Slurm slot 4 of 4 XeonE52620v4
const std::vector<int> slurm4of4XeonE52620v4 = { 12, 13, 14, 15, 28, 29, 30, 31 };
//! Slurm slot 1 of 4 Xeon4116
const std::vector<int> slurm1of4Xeon4116 = { 0, 1, 2, 3, 4, 5, 24, 25, 26, 27, 28, 29 };
//! Slurm slot 2 of 4 Xeon4116
const std::vector<int> slurm2of4Xeon4116 = { 6, 7, 8, 9, 10, 11, 30, 31, 32, 33, 34, 35 };
//! Slurm slot 3 of 4 Xeon4116
const std::vector<int> slurm3of4Xeon4116 = { 12, 13, 14, 15, 16, 17, 36, 37, 38, 39, 40, 41 };
//! Slurm slot 4 of 4 Xeon4116
const std::vector<int> slurm4of4Xeon4116 = { 18, 19, 20, 21, 22, 23, 42, 43, 44, 45, 46, 47 };

/*! \brief
 * Mocking test with data collected from a few systems.
 *
 * The easiest way to generate reference data
 * for these test is to temporarily replace the mocking constructor with a call to detecting the
 * topology on the hardware in question, and then only change the reference field IsThisSystem to
 * 'false'. Consult the README file in the tests directory for how to capture the mock input data to
 * use e.g. for logical processor IDs from cpuinfo or filesystem reference files.
 */
class MockHardwareTopologyTest :
    public ::testing::Test,
    public ::testing::WithParamInterface<HardwareDetectionTestParams>
{
private:
    TestReferenceData    data_;
    TestReferenceChecker checker_;

public:
    MockHardwareTopologyTest() : checker_(data_.rootChecker()) {}

    void runTest(const HardwareTopology& hwTop);
};

void MockHardwareTopologyTest::runTest(const HardwareTopology& hwTop)
{
    checkHardwareTopology(&checker_, hwTop);
}

TEST_P(MockHardwareTopologyTest, DetectsHardware)
{
    const auto& params                = GetParam();
    const auto& logicalProcessorIdMap = std::get<0>(params);
    const auto& allowedCpus           = std::get<1>(params);
    const auto& filePath              = std::get<2>(params);

    ProcessorMap allowedProcessorIdMap;
    const bool   checkCpus = !allowedCpus.empty();

    if (checkCpus)
    {
        for (const int cpu : allowedCpus)
        {
            // Use at(), since operator[] on maps returns non-const reference
            allowedProcessorIdMap[cpu] = logicalProcessorIdMap.at(cpu);
        }
    }

    HardwareTopology hwTopFromMap(checkCpus ? allowedProcessorIdMap : logicalProcessorIdMap,
                                  TestFileManager::getInputFilePath(filePath).string());
    runTest(hwTopFromMap);
}

INSTANTIATE_TEST_SUITE_P(XeonE52620,
                         MockHardwareTopologyTest,
                         ::testing::Combine(::testing::Values(logicalProcessorIdMapXeonE52620v4),
                                            ::testing::Values(emptyCPUVector,
                                                              slurm1of4XeonE52620v4,
                                                              slurm2of4XeonE52620v4,
                                                              slurm3of4XeonE52620v4,
                                                              slurm4of4XeonE52620v4),
                                            ::testing::Values("XeonE52620v4_Cgroups2NoLimit")));


INSTANTIATE_TEST_SUITE_P(Xeon4116,
                         MockHardwareTopologyTest,
                         ::testing::Combine(::testing::Values(logicalProcessorIdMapXeon4116),
                                            ::testing::Values(emptyCPUVector,
                                                              slurm1of4Xeon4116,
                                                              slurm2of4Xeon4116,
                                                              slurm3of4Xeon4116,
                                                              slurm4of4Xeon4116),
                                            ::testing::Values("Xeon4116_Cgroups1CpuLimit400pct")));


INSTANTIATE_TEST_SUITE_P(
        Core12900K,
        MockHardwareTopologyTest,
        ::testing::Combine(::testing::Values(logicalProcessorIdMapCore12900K),
                           ::testing::Values(emptyCPUVector),
                           ::testing::Values("Core12900K_Cgroups2NoLimit",
                                             "Core12900K_Cgroups2CpuLimit650pct")));

//! Test parameter struct for special systems.
using SpecialSystemHardwareParams = std::tuple<int, int, std::string>;

class MockHardwareTopologySpecialSystemTest :
    public ::testing::Test,
    public ::testing::WithParamInterface<SpecialSystemHardwareParams>
{
private:
    TestReferenceData    data_;
    TestReferenceChecker checker_;

public:
    MockHardwareTopologySpecialSystemTest() : checker_(data_.rootChecker()) {}

    void runTest(const HardwareTopology& hwTop);
};

void MockHardwareTopologySpecialSystemTest::runTest(const HardwareTopology& hwTop)
{
    checkHardwareTopology(&checker_, hwTop);
}

TEST_P(MockHardwareTopologySpecialSystemTest, DetectsHardware)
{
    const auto& params            = GetParam();
    const int   systemSize        = std::get<0>(params);
    const int   firstAvailableCpu = std::get<1>(params);
    const auto& filePath          = std::get<2>(params);

    std::vector<int> allowedCpus(systemSize);
    std::iota(allowedCpus.begin(), allowedCpus.end(), firstAvailableCpu);

    HardwareTopology hwTopFromSysFs(TestFileManager::getInputFilePath(filePath).string(), allowedCpus);
    runTest(hwTopFromSysFs);
};

INSTANTIATE_TEST_SUITE_P(
        Core12900K,
        MockHardwareTopologySpecialSystemTest,
        ::testing::Combine(::testing::Values(24),
                           ::testing::Values(0),
                           ::testing::Values("Core12900K_Cgroups2CpuLimit650pct")));

INSTANTIATE_TEST_SUITE_P(Power9,
                         MockHardwareTopologySpecialSystemTest,
                         ::testing::Combine(::testing::Values(32),
                                            ::testing::Values(0),
                                            ::testing::Values("Power9_Cgroups2NoLimit")));

INSTANTIATE_TEST_SUITE_P(A64fx,
                         MockHardwareTopologySpecialSystemTest,
                         ::testing::Combine(::testing::Values(48),
                                            ::testing::Values(12),
                                            ::testing::Values("A64fx_Cgroups1CpuLimit2400pct")));


} // namespace
} // namespace test
} // namespace gmx
