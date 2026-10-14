/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests that the forces written to the trajectory are complete.
 *
 * \author Andrey Alekseenko <al42and@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <cstdlib>

#include <array>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/mpitest.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "moduletest.h"
#include "trajectorycomparison.h"

namespace gmx
{
namespace test
{
namespace
{

//! Parameters for the parameterized test fixture: ( integrator; offload flags; number of PME ranks )
using ForceOutputTestParameters = std::tuple<std::string, std::string, int>;

const NameOfTestFromTuple<ForceOutputTestParameters> sc_testNamer{ std::make_tuple(
        useString,
        useString,
        [](int numPmeRanks) { return formatString("npme%d", numPmeRanks); }) };

// Every configuration has to reproduce the same forces for a given integrator.
const RefDataFilenameMaker<ForceOutputTestParameters> sc_refDataFilenameMaker{
    std::make_tuple(useString, toEmptyString<std::string>, toEmptyString<int>)
};

//! Integrators to test. Each gets its own TPR file and reference data.
const std::array<std::string, 2> c_integrators{ "md", "sd" };

class ForceOutputTest : public MdrunTestFixture, public ::testing::WithParamInterface<ForceOutputTestParameters>
{
public:
    inline static const std::string                            s_inputFile = "spc-and-methanol";
    inline static std::unordered_map<std::string, std::string> s_tprFileNames;
    inline static std::mutex                                   s_mutexForTestFileManager;
    inline static std::unique_ptr<TestFileManager>             s_testFileManager;
    static void                                                SetUpTestSuite();
    static void                                                TearDownTestSuite();
    static MessageStringCollector getSkipMessagesIfNecessary(const CommandLine& commandLine,
                                                             const std::string& integrator,
                                                             int                numPmeRanks);
};

void ForceOutputTest::SetUpTestSuite()
{
    MdrunTestFixture::SetUpTestSuite();
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        if (!s_testFileManager)
        {
            s_testFileManager = std::make_unique<TestFileManager>();
        }
    }

    for (const std::string& integrator : c_integrators)
    {
        // Stochastic dynamics needs a temperature-coupling group, and a fixed seed
        const std::string integratorMdpOptions =
                (integrator == "sd") ? "integrator              = sd\n"
                                       "tc-grps                 = System\n"
                                       "tau-t                   = 1\n"
                                       "ref-t                   = 300\n"
                                       "ld-seed                 = 1234\n"
                                     : "integrator              = " + integrator + "\n";
        SimulationRunner runner(s_testFileManager.get());
        runner.useTopGroAndNdxFromDatabase(s_inputFile);
        runner.useStringAsMdpFile(
                "dt                      = 0.002\n"
                "nsteps                  = 20\n"
                "nstxout                 = 3\n"
                "nstfout                 = 3\n"
                "nstcalcenergy           = 100\n"
                "nstenergy               = 100\n"
                "nstlog                  = 100\n"
                "nstlist                 = 10\n"
                "cutoff-scheme           = Verlet\n"
                "verlet-buffer-tolerance = -1\n"
                "rlist                   = 1.0\n"
                "coulombtype             = PME\n"
                "rcoulomb                = 0.9\n"
                "fourierspacing          = 0.12\n"
                "pme-order               = 4\n"
                "vdwtype                 = cut-off\n"
                "rvdw                    = 0.9\n"
                "tcoupl                  = no\n"
                "pcoupl                  = no\n"
                "constraints             = none\n"
                "gen-vel                 = no\n"
                "pbc                     = xyz\n"
                + integratorMdpOptions);
        runner.tprFileName_ = s_testFileManager->getTemporaryFilePath(integrator + ".tpr").string();
        ASSERT_EQ(0, runner.callGrompp());
        s_tprFileNames[integrator] = runner.tprFileName_;
    }
}

void ForceOutputTest::TearDownTestSuite()
{
    MdrunTestFixture::TearDownTestSuite();
    {
        std::lock_guard<std::mutex> testFileManagerLock(s_mutexForTestFileManager);
        s_testFileManager.reset(nullptr);
    }
}

MessageStringCollector ForceOutputTest::getSkipMessagesIfNecessary(const CommandLine& commandLine,
                                                                   const std::string& integrator,
                                                                   int                numPmeRanks)
{
    MessageStringCollector messages;
    messages.startContext("Test is being skipped because:");

    const int  numRanks           = getNumberOfTestMpiRanks();
    const bool useSeparatePmeRank = (numPmeRanks > 0);
    const bool useGpu             = (commandLine.argumentOf("-nb").value() == "gpu");
    const bool useGpuPme          = (commandLine.argumentOf("-pme").value() == "gpu");
    const bool useGpuUpdate       = (commandLine.argumentOf("-update").value() == "gpu");

    messages.appendIf(useSeparatePmeRank && numRanks == 1,
                      "it targets a separate PME rank, but the simulation is using only one rank");
    messages.appendIf(useGpu && getCompatibleDevices(s_hwinfo->deviceInfoList).empty(),
                      "it targets GPU execution, but no compatible devices were detected");
    messages.appendIf(useGpu && !GpuConfigurationCapabilities::Nonbonded,
                      "it targets GPU execution, which this build does not support");
    messages.appendIf(useGpuPme && !GpuConfigurationCapabilities::Pme,
                      "it targets PME on the GPU, which this build does not support");
    messages.appendIf(useGpuUpdate && integrator == "md" && !GpuConfigurationCapabilities::UpdateLeapfrog,
                      "it targets the update on the GPU, which this build does not support");
    messages.appendIf(
            useGpuUpdate && integrator == "sd" && !GpuConfigurationCapabilities::UpdateSD,
            "it targets the update with SD on the GPU, which this build does not support");
    messages.appendIf(useGpuPme && !useSeparatePmeRank && numRanks > 1,
                      "it targets PME decomposition, which this test does not set up");

    messages.finishContext();
    return messages;
}

TEST_P(ForceOutputTest, ReproducesForces)
{
    const auto [integrator, offloadFlags, numPmeRanks] = GetParam();
    CommandLine commandLine(splitString(offloadFlags));
    commandLine.addOption("-npme", numPmeRanks);
    commandLine.addOption("-notunepme");

    MessageStringCollector skipMessages = getSkipMessagesIfNecessary(commandLine, integrator, numPmeRanks);
    if (!skipMessages.isEmpty())
    {
        GTEST_SKIP() << skipMessages.toString();
    }

    runner_.tprFileName_ = s_tprFileNames.at(integrator);
    ASSERT_EQ(0, runner_.callMdrun(commandLine));

    // Do the checks only on the main rank
    if (gmx_node_rank() != 0)
    {
        return;
    }

    // Reference data comes from double-precision CPU build
    TrajectoryTolerances tolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    tolerances.coordinates          = relativeToleranceAsFloatingPoint(1.0, 1.0e-5);
    tolerances.forces               = relativeToleranceAsFloatingPoint(1000.0, 1.0e-3);
    TrajectoryFrameMatchSettings matchSettings;
    matchSettings.coordinatesComparison = ComparisonConditions::MustCompare;
    matchSettings.velocitiesComparison  = ComparisonConditions::NoComparison;
    matchSettings.forcesComparison      = ComparisonConditions::MustCompare;

    TestReferenceData    refData(sc_refDataFilenameMaker(GetParam()));
    TestReferenceChecker rootChecker(refData.rootChecker());
    checkTrajectoryAgainstReferenceData(runner_.fullPrecisionTrajectoryFileName_,
                                        TrajectoryComparison(matchSettings, tolerances),
                                        &rootChecker);
}


constexpr std::array<int, 2> c_numberOfPmeRanks{ 0, 1 };

const auto c_forceOutputConfigurations =
        ::testing::Combine(::testing::ValuesIn(c_integrators),
                           ::testing::Values("-nb cpu -pme cpu -update cpu",
                                             "-nb gpu -pme cpu -update cpu",
                                             "-nb gpu -pme gpu -update cpu",
                                             "-nb gpu -pme gpu -update gpu"),
                           ::testing::ValuesIn(c_numberOfPmeRanks));

INSTANTIATE_TEST_SUITE_P(TrajectoryForces, ForceOutputTest, c_forceOutputConfigurations, sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
