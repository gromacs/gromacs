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

#include "gromacs/hardware/hardwaretopology.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

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

// Mocking test with data collected from a few systems. The easiest way to generate reference data
// for these test is to temporarily replace the mocking constructor with a call to detecting the
// topology on the hardware in question, and then only change the reference field IsThisSystem to
// 'false'. Consult the README file in the tests directory for how to capture the mock input data to
// use e.g. for logical processor IDs from cpuinfo or filesystem reference files.

class MockHardwareTopologyXeonE52620v4Test : public ::testing::Test
{
public:
    const std::map<int, std::array<int, 3>> logicalProcessorIdMap_ = {
        { 0, { 0, 0, 0 } },  { 1, { 0, 1, 0 } },  { 2, { 0, 2, 0 } },  { 3, { 0, 3, 0 } },
        { 4, { 0, 4, 0 } },  { 5, { 0, 5, 0 } },  { 6, { 0, 6, 0 } },  { 7, { 0, 7, 0 } },
        { 8, { 1, 0, 0 } },  { 9, { 1, 1, 0 } },  { 10, { 1, 2, 0 } }, { 11, { 1, 3, 0 } },
        { 12, { 1, 4, 0 } }, { 13, { 1, 5, 0 } }, { 14, { 1, 6, 0 } }, { 15, { 1, 7, 0 } },
        { 16, { 0, 0, 1 } }, { 17, { 0, 1, 1 } }, { 18, { 0, 2, 1 } }, { 19, { 0, 3, 1 } },
        { 20, { 0, 4, 1 } }, { 21, { 0, 5, 1 } }, { 22, { 0, 6, 1 } }, { 23, { 0, 7, 1 } },
        { 24, { 1, 0, 1 } }, { 25, { 1, 1, 1 } }, { 26, { 1, 2, 1 } }, { 27, { 1, 3, 1 } },
        { 28, { 1, 4, 1 } }, { 29, { 1, 5, 1 } }, { 30, { 1, 6, 1 } }, { 31, { 1, 7, 1 } },
    };
};

TEST_F(MockHardwareTopologyXeonE52620v4Test, Cgroups2NoLimit)
{
    // Plain 2-socket Xeon E5-2620v4 system, no cpuset limitations, no cpu usage limits
    // Only the allowed CPU array differs, so we use the same mock root path for all these tests.
    HardwareTopology hwTopFromMap(
            logicalProcessorIdMap_, TestFileManager::getInputFilePath("XeonE52620v4_Cgroups2NoLimit"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

TEST_F(MockHardwareTopologyXeonE52620v4Test, SlurmSlot1of4_Cgroups2NoLimit)
{
    // Plain 2-socket Xeon E5-2620v4 system, using 1st Slurm slot out of 4 on system. No cpu usage limits.
    // Only the allowed CPU array differs, so we use the same mock root path for all these tests.
    const std::vector<int>            allowedCpus{ 0, 1, 2, 3, 16, 17, 18, 19 };
    std::map<int, std::array<int, 3>> allowedProcessorIdMap;

    for (const int cpu : allowedCpus)
    {
        // Use at(), since operator[] on maps returns non-const reference
        allowedProcessorIdMap[cpu] = logicalProcessorIdMap_.at(cpu);
    }
    HardwareTopology hwTopFromMap(
            allowedProcessorIdMap, TestFileManager::getInputFilePath("XeonE52620v4_Cgroups2NoLimit"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

TEST_F(MockHardwareTopologyXeonE52620v4Test, SlurmSlot2of4_Cgroups2NoLimit)
{
    // Plain 2-socket Xeon E5-2620v4 system, using 2nd Slurm slot out of 4 on system. No cpu usage limits.
    // Only the allowed CPU array differs, so we use the same mock root path for all these tests.
    const std::vector<int>            allowedCpus{ 4, 5, 6, 7, 20, 21, 22, 23 };
    std::map<int, std::array<int, 3>> allowedProcessorIdMap;

    for (int cpu : allowedCpus)
    {
        // Use at(), since operator[] on maps returns non-const reference
        allowedProcessorIdMap[cpu] = logicalProcessorIdMap_.at(cpu);
    }

    HardwareTopology hwTopFromMap(
            logicalProcessorIdMap_, TestFileManager::getInputFilePath("XeonE52620v4_Cgroups2NoLimit"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

TEST_F(MockHardwareTopologyXeonE52620v4Test, SlurmSlot3of4_Cgroups2NoLimit)
{
    // Plain 2-socket Xeon E5-2620v4 system, using 3rd Slurm slot out of 4 on system. No cpu usage limits.
    // Only the allowed CPU array differs, so we use the same mock root path for all these tests.
    const std::vector<int> allowedCpus{ 8, 9, 10, 11, 24, 25, 26, 27 };
    HardwareTopology       hwTopFromSysFs(
            TestFileManager::getInputFilePath("XeonE52620v4_Cgroups2NoLimit"), allowedCpus);

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromSysFs);
}

TEST_F(MockHardwareTopologyXeonE52620v4Test, SlurmSlot4of4_Cgroups2NoLimit)
{
    // Plain 2-socket Xeon E5-2620v4 system, using 4th Slurm slot out of 4 on system. No cpu usage limits.
    // Only the allowed CPU array differs, so we use the same mock root path for all these tests.
    const std::vector<int> allowedCpus{ 12, 13, 14, 15, 28, 29, 30, 31 };
    HardwareTopology       hwTopFromSysFs(
            TestFileManager::getInputFilePath("XeonE52620v4_Cgroups2NoLimit"), allowedCpus);

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromSysFs);
}

class MockHardwareTopologyXeon4116Test : public ::testing::Test
{
public:
    const std::map<int, std::array<int, 3>> logicalProcessorIdMap_ = {
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
};

TEST_F(MockHardwareTopologyXeon4116Test, Cgroups1CpuLimit400pct)
{
    // Kubernetes pod running on Xeon Silver 4116, with a cpu limit of 400%.
    HardwareTopology hwTopFromMap(
            logicalProcessorIdMap_, TestFileManager::getInputFilePath("Xeon4116_Cgroups1CpuLimit400pct"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

class MockHardwareTopologyCore12900KTest : public ::testing::Test
{
public:
    const std::map<int, std::array<int, 3>> logicalProcessorIdMap_ = {
        { 0, { 0, 0, 0 } },   { 1, { 0, 0, 1 } },   { 2, { 0, 4, 0 } },   { 3, { 0, 4, 1 } },
        { 4, { 0, 8, 0 } },   { 5, { 0, 8, 1 } },   { 6, { 0, 12, 0 } },  { 7, { 0, 12, 1 } },
        { 8, { 0, 16, 0 } },  { 9, { 0, 16, 1 } },  { 10, { 0, 20, 0 } }, { 11, { 0, 20, 1 } },
        { 12, { 0, 24, 0 } }, { 13, { 0, 24, 1 } }, { 14, { 0, 28, 0 } }, { 15, { 0, 28, 1 } },
        { 16, { 0, 32, 0 } }, { 17, { 0, 33, 0 } }, { 18, { 0, 34, 0 } }, { 19, { 0, 35, 0 } },
        { 20, { 0, 36, 0 } }, { 21, { 0, 37, 0 } }, { 22, { 0, 38, 0 } }, { 23, { 0, 39, 0 } },
    };
};

TEST_F(MockHardwareTopologyCore12900KTest, MapCgroups2NoLimit)
{
    // Alder lake (core i9-12900k) with 8 performance cores (SMT2) and 8 efficiency cores (no SMT).
    // This test has cgroups directories present, but the cpu controller is not enabled, meaning no limit.
    HardwareTopology hwTopFromMap(logicalProcessorIdMap_,
                                  TestFileManager::getInputFilePath("Core12900K_Cgroups2NoLimit"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

TEST_F(MockHardwareTopologyCore12900KTest, MapCgroups2CpuLimit650pct)
{
    // Docker container running on Core i9-12900k, with cgroups2 controller limiting cpu to 650%.
    // This test parses cpu topology from the map of processor IDs (which is derived from x86 APIC IDs)
    HardwareTopology hwTopFromMap(
            logicalProcessorIdMap_, TestFileManager::getInputFilePath("Core12900K_Cgroups2CpuLimit650pct"));

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromMap);
}

TEST_F(MockHardwareTopologyCore12900KTest, SysFsCgroups2CpuLimit650pct)
{
    // Docker container running on Core i9-12900k, with cgroups2 controller limiting cpu to 650%.
    // This test parses cpu topology from the file system.
    std::vector<int> allowedCpus(24); // 8*2+8 = 24 logical cpus.
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);

    HardwareTopology hwTopFromSysFs(
            TestFileManager::getInputFilePath("Core12900K_Cgroups2CpuLimit650pct"), allowedCpus);

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromSysFs);
}

TEST(MockHardwareTopologyPower9Test, SysFsCgroups2NoLimit)
{
    // Power9, 8-core CPU with SMT4
    std::vector<int> allowedCpus(32); // 8*4=32 logical cpus
    std::iota(allowedCpus.begin(), allowedCpus.end(), 0);

    HardwareTopology hwTopFromSysFs(TestFileManager::getInputFilePath("Power9_Cgroups2NoLimit"), allowedCpus);

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromSysFs);
}

TEST(MockHardwareTopologyA64fxTest, SysFsCgroups1Limit2400pct)
{
    // A64fx. 52 cores, with quite special numbering, and only 48 available.
    // The cgroups layout is also special, with multiple different mount points. To check that
    // we parse all such files, we have faked a cgroups limit of 2400% CPU.
    std::vector<int> allowedCpus(48);
    std::iota(allowedCpus.begin(), allowedCpus.end(), 12); // Logical CPUs 12-59 are available to users on this A64fx

    HardwareTopology hwTopFromSysFs(
            TestFileManager::getInputFilePath("A64fx_Cgroups1CpuLimit2400pct"), allowedCpus);

    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());

    checkHardwareTopology(&checker, hwTopFromSysFs);
}
} // namespace
} // namespace test
} // namespace gmx
