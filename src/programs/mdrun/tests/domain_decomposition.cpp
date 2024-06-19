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
 * Tests special cases in domain decomposition
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdlib>

#include <algorithm>
#include <array>
#include <filesystem>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"

namespace
{

// Define all the possible parameters for the test

enum class ElectrostaticsFlavor : int
{
    ReactionField,
    Pme,
    Count
};

const char* enumValueToString(const ElectrostaticsFlavor enumValue)
{
    static constexpr gmx::EnumerationArray<ElectrostaticsFlavor, const char*> s_names = { "RF",
                                                                                          "PME" };
    return s_names[enumValue];
}

enum class CouplingFlavor : int
{
    No,
    TemperatureAndPressure,
    Count
};

const char* enumValueToString(const CouplingFlavor enumValue)
{
    static constexpr gmx::EnumerationArray<CouplingFlavor, const char*> s_names = {
        "No", "TemperatureAndPressure"
    };
    return s_names[enumValue];
}

enum class OffloadFlavor : int
{
    Cpu,
    Gpu,
    Count
};

const char* enumValueToString(const OffloadFlavor enumValue)
{
    // Must be lower-case so we can use it when constructing mdrun command line
    static constexpr gmx::EnumerationArray<OffloadFlavor, const char*> s_names = { "cpu", "gpu" };
    return s_names[enumValue];
}

// Helper aliases
using NonbondedFlavor = OffloadFlavor;
using UpdateFlavor    = OffloadFlavor;
using PmeFlavor       = OffloadFlavor;

enum class SeparatePmeRankFlavor : int
{
    None = 0,
    One  = 1,
    Two  = 2,
    Count
};

//! \brief Tuple containing parameters for MDP/TPR file generation
using MdpFlavor = std::tuple<ElectrostaticsFlavor, CouplingFlavor>;
//! \brief Tuple containing parameters for mdrun command line
using RuntimeFlavor = std::tuple<NonbondedFlavor, PmeFlavor, UpdateFlavor, SeparatePmeRankFlavor>;

//! \brief Parameters for parametrized test fixture
using DomDecSpecialCasesTestParameters = std::tuple<MdpFlavor, RuntimeFlavor>;

//! \brief Verify whether the test configuration is valid and worth running
std::optional<std::string> reasonsTestIsInvalid(MdpFlavor       mdpFlavor,
                                                RuntimeFlavor   runtimeFlavor,
                                                bool gmx_unused haveCompatibleDevices)
{
    const auto electrostaticsFlavor = std::get<0>(mdpFlavor);
    const auto [nonbondedFlavor, pmeFlavor, updateFlavor, separatePmeRankFlavor] = runtimeFlavor;

    const int  numRanks       = gmx::test::getNumberOfTestMpiRanks();
    const bool haveAnyGpuWork = (pmeFlavor == PmeFlavor::Gpu || updateFlavor == UpdateFlavor::Gpu
                                 || nonbondedFlavor == NonbondedFlavor::Gpu);

    gmx::MessageStringCollector errorReasons;
    errorReasons.startContext("Test configuration is invalid:");

    // MDP-related
    errorReasons.appendIf(electrostaticsFlavor != ElectrostaticsFlavor::Pme && pmeFlavor == PmeFlavor::Gpu,
                          "Cannot offload PME to GPU when PME is not used");
    errorReasons.appendIf(electrostaticsFlavor != ElectrostaticsFlavor::Pme
                                  && separatePmeRankFlavor != SeparatePmeRankFlavor::None,
                          "Cannot have separate PME ranks when PME is not used");
    // GPU-related
    errorReasons.appendIf(haveAnyGpuWork && (GMX_GPU == 0),
                          "Cannot use GPU offload in non-GPU build");
#if GMX_GPU
    errorReasons.appendIf(haveAnyGpuWork && !haveCompatibleDevices,
                          "Cannot use GPU offload without a compatible GPU");
    errorReasons.appendIf(GMX_GPU_OPENCL && updateFlavor == UpdateFlavor::Gpu,
                          "GPU Update not supported with OpenCL");
    errorReasons.appendIf(updateFlavor == UpdateFlavor::Gpu && pmeFlavor == PmeFlavor::Cpu
                                  && separatePmeRankFlavor != SeparatePmeRankFlavor::None,
                          "Can not use GPU update and CPU PME on a separate rank");
    errorReasons.appendIf(GMX_GPU_HIP, "HIP kernels are not implemented yet");
#endif
    errorReasons.appendIf(haveAnyGpuWork && nonbondedFlavor == NonbondedFlavor::Cpu,
                          "Cannot offload PME or Update to GPU without offloading Nonbondeds");
    if ((getenv("GMX_GPU_PME_DECOMPOSITION")) == nullptr)
    {
        errorReasons.appendIf(
                pmeFlavor == PmeFlavor::Gpu && separatePmeRankFlavor == SeparatePmeRankFlavor::Two,
                "Cannot use more than one separate PME rank with GPU PME");
    }
    // Related to the number of ranks
    errorReasons.appendIf(numRanks == 1 && separatePmeRankFlavor != SeparatePmeRankFlavor::None,
                          "Cannot use separate PME rank with only one total rank");
    errorReasons.appendIf(numRanks > 1 && separatePmeRankFlavor == SeparatePmeRankFlavor::None
                                  && pmeFlavor == PmeFlavor::Gpu,
                          "Cannot use GPU PME offload with multiple PME+PP ranks");
    errorReasons.appendIf(numRanks == 2 && separatePmeRankFlavor == SeparatePmeRankFlavor::Two,
                          "Cannot use two separate PME ranks when there are only two ranks total");
    errorReasons.finishContext();
    if (errorReasons.isEmpty())
    {
        return std::nullopt;
    }
    else
    {
        return std::make_optional(errorReasons.toString());
    }
}

//! \brief Help GoogleTest name our tests
std::string nameOfTest(const testing::TestParamInfo<DomDecSpecialCasesTestParameters>& info)
{
    const auto [mdpFlavor, runtimeFlavor]                                        = info.param;
    const auto [electrostaticsFlavor, couplingFlavor]                            = mdpFlavor;
    const auto [nonbondedFlavor, pmeFlavor, updateFlavor, separatePmeRankFlavor] = runtimeFlavor;

    std::string testName = gmx::formatString("%s_%s_coupling_nb%s_pme%s_update%s_npme%d",
                                             enumValueToString(electrostaticsFlavor),
                                             enumValueToString(couplingFlavor),
                                             enumValueToString(nonbondedFlavor),
                                             enumValueToString(pmeFlavor),
                                             enumValueToString(updateFlavor),
                                             static_cast<int>(separatePmeRankFlavor));

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = gmx::replaceAll(testName, "-", "");
    testName = gmx::replaceAll(testName, " ", "_");
    return testName;
}

//! \brief Generate the contents of the MDP file
std::string buildMdpInputFileContent(MdpFlavor mdpFlavor)
{
    const auto [electrostaticsFlavor, couplingFlavor] = mdpFlavor;
    std::string mdpInputFileContent{
        "cutoff-scheme = verlet\n"
        "nsteps = 20\n"
        "nstcomm = 10\n"
        "nstcalcenergy = 5\n"
        "nstenergy = 5\n"
        "vdwtype = cut-off\n"
    };

    switch (electrostaticsFlavor)
    {
        case ElectrostaticsFlavor::ReactionField:
            mdpInputFileContent.append("coulombtype = reaction-field\n");
            break;
        case ElectrostaticsFlavor::Pme:
            mdpInputFileContent.append(
                    "coulombtype = pme\n"
                    "pme-order = 4\n");
            break;
        default: GMX_RELEASE_ASSERT(false, "Invalid value");
    }

    switch (couplingFlavor)
    {
        case CouplingFlavor::No: mdpInputFileContent.append("pcoupl = no\ntcoupl = no\n"); break;
        case CouplingFlavor::TemperatureAndPressure:
            mdpInputFileContent.append(
                    "pcoupl = c-rescale\n"
                    "ref-p = 1.0\n"
                    "compressibility = 4.5e-5\n");
            mdpInputFileContent.append(
                    "tcoupl = v-rescale\n"
                    "tc-grps = System\n"
                    "ref_t = 298\n"
                    "tau_t = 0.5\n");
            break;
        default: GMX_RELEASE_ASSERT(false, "Invalid value");
    }
    return mdpInputFileContent;
}

//! \brief Generate the mdrun command line
gmx::test::CommandLine buildMdrunCommandLine(RuntimeFlavor runtimeFlavor)
{
    const auto [nonbondedFlavor, pmeFlavor, updateFlavor, separatePmeRankFlavor] = runtimeFlavor;
    gmx::test::CommandLine mdrunCommandLine;
    mdrunCommandLine.addOption("-nb", enumValueToString(nonbondedFlavor));
    mdrunCommandLine.addOption("-pme", enumValueToString(pmeFlavor));
    mdrunCommandLine.addOption("-update", enumValueToString(updateFlavor));
    mdrunCommandLine.addOption("-npme", static_cast<int>(separatePmeRankFlavor));
    mdrunCommandLine.addOption("-notunepme");
    return mdrunCommandLine;
}

// Needed for std::unordered_map<MdpFlavor, ...>
struct MdpFlavorHash
{
    std::size_t operator()(const MdpFlavor& mdpFlavor) const
    {
        const auto [electrostaticsFlavor, couplingFlavor] = mdpFlavor;
        return static_cast<int>(electrostaticsFlavor) * static_cast<int>(CouplingFlavor::Count)
               + static_cast<int>(couplingFlavor);
    }
};

// We need to iterate over all valid MdpFlavor values when generating TPRs
constexpr std::array<MdpFlavor, 4> sc_mdpFlavors = {
    std::make_tuple(ElectrostaticsFlavor::ReactionField, CouplingFlavor::No),
    std::make_tuple(ElectrostaticsFlavor::ReactionField, CouplingFlavor::TemperatureAndPressure),
    std::make_tuple(ElectrostaticsFlavor::Pme, CouplingFlavor::No),
    std::make_tuple(ElectrostaticsFlavor::Pme, CouplingFlavor::TemperatureAndPressure),
};

//! Test fixture for domain decomposition special cases
class DomDecSpecialCasesTest :
    public gmx::test::MdrunTestFixture,
    public ::testing::WithParamInterface<DomDecSpecialCasesTestParameters>
{
public:
    //! Names of tpr files built by grompp in SetUpTestSuite to run in tests
    inline static std::unordered_map<MdpFlavor, std::string, MdpFlavorHash> s_tprFileNames;
    //! Mutex to protect creation of the TestFileManager with thread-MPI
    inline static std::mutex s_mutexForTestFileManager;
    //! Manager needed in SetUpTestSuite to handle making tpr files
    inline static std::unique_ptr<gmx::test::TestFileManager> s_testFileManager;
    //! Runs grompp to prepare the tpr files that are reused in the tests
    static void SetUpTestSuite();
    //! Cleans up the tpr files
    static void TearDownTestSuite();
};


void DomDecSpecialCasesTest::SetUpTestSuite()
{
    MdrunTestFixture::SetUpTestSuite();
    // Ensure we only make one TestFileManager per process, which
    // ensures there is no race with thread-MPI. Whichever thread gets
    // the lock first initializes s_testFileManager, and by the time
    // the rest get the lock it is already initialized.  Without
    // thread-MPI, there's only one thread, so the initialization is
    // trivial.
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        if (!s_testFileManager)
        {
            s_testFileManager = std::make_unique<gmx::test::TestFileManager>();
        }
    }

    for (MdpFlavor mdpFlavor : sc_mdpFlavors)
    {
        gmx::test::SimulationRunner runner(s_testFileManager.get());
        runner.useStringAsMdpFile(buildMdpInputFileContent(mdpFlavor));
        runner.useTopGroAndNdxFromDatabase("spc2");
        // Change the GRO file to the one where both water molecules are close by.
        runner.useGroFromDatabase("spc2_and_vacuum");

        std::string tprFileNameSuffix = gmx::formatString("%s_%s.tpr",
                                                          enumValueToString(std::get<0>(mdpFlavor)),
                                                          enumValueToString(std::get<1>(mdpFlavor)));
        std::replace(tprFileNameSuffix.begin(), tprFileNameSuffix.end(), ' ', '_');
        runner.tprFileName_ = s_testFileManager->getTemporaryFilePath(tprFileNameSuffix).string();
        // Note that only one rank actually generates a tpr file
        ASSERT_EQ(0, runner.callGrompp());
        s_tprFileNames[mdpFlavor] = runner.tprFileName_;
    }
}

void DomDecSpecialCasesTest::TearDownTestSuite()
{
    MdrunTestFixture::TearDownTestSuite();
    // Ensure we only clean up the TestFileManager once per process, avoiding races with thread-MPI.
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        s_testFileManager.reset(nullptr);
    }
}

//! When run with 2+ domains, ensures an empty cell, to make sure that zero-sized things work
TEST_P(DomDecSpecialCasesTest, EmptyDomain)
{
    gmx::test::checkTestNameLength();
    const auto [mdpFlavor, runtimeFlavor] = GetParam();

    const bool haveCompatibleDevices = !getCompatibleDevices(s_hwinfo->deviceInfoList).empty();
    if (auto reasons = reasonsTestIsInvalid(mdpFlavor, runtimeFlavor, haveCompatibleDevices);
        reasons.has_value())
    {
        /* We might wish to make sure that mdrun fails in such cases,
         * but death tests do not work well with any kind of threading. */
        GTEST_SKIP() << *reasons;
    }

    gmx::test::CommandLine mdrunCommandLine = buildMdrunCommandLine(runtimeFlavor);
    runner_.tprFileName_                    = s_tprFileNames[mdpFlavor];
    EXPECT_EQ(runner_.callMdrun(mdrunCommandLine), 0);
}

#if GMX_GPU
constexpr std::array<OffloadFlavor, 2> sc_offloadFlavors{ OffloadFlavor::Cpu, OffloadFlavor::Gpu };
#else
// We would have skippen GPU tests anyway, but why even bother instantiating them
constexpr std::array<OffloadFlavor, 1> sc_offloadFlavors{ OffloadFlavor::Cpu };
#endif

INSTANTIATE_TEST_SUITE_P(
        DomainDecomposition,
        DomDecSpecialCasesTest,
        ::testing::Combine(::testing::ValuesIn(sc_mdpFlavors),
                           ::testing::Combine(::testing::ValuesIn(sc_offloadFlavors),
                                              ::testing::ValuesIn(sc_offloadFlavors),
                                              ::testing::ValuesIn(sc_offloadFlavors),
                                              ::testing::Values(SeparatePmeRankFlavor::None,
                                                                SeparatePmeRankFlavor::One,
                                                                SeparatePmeRankFlavor::Two))

                                   ),
        nameOfTest);

} // namespace
