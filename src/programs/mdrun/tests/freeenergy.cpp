/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Tests to compare free energy simulations to reference
 *
 * The tests are parameterized over CPU and every compatible GPU. CPU and GPU
 * runs share reference data; only the test name records the hardware context.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/hardware_test_fixture.h"
#include "testutils/mpitest.h"
#include "testutils/naming.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/xvgtest.h"

#include "programs/mdrun/tests/comparison_helpers.h"
#include "programs/mdrun/tests/energycomparison.h"
#include "programs/mdrun/tests/trajectorycomparison.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx::test
{
namespace
{

using MaxNumWarnings           = int;
using ListOfInteractionsToTest = std::vector<InteractionFunction>;

//! Physical system under test. This contains what selects the reference data.
using FreeEnergyReferenceInputConfig = std::tuple<std::string, MaxNumWarnings, ListOfInteractionsToTest>;

//! No execution modes.
using FreeEnergyReferenceTestHelper =
        HardwareAndExecutionTestHelper<FreeEnergyReferenceInputConfig, std::tuple<>>;

//! Keep the historical reference-data suffix so existing files still match.
std::string formatFreeEnergySimulationName(const std::string& simulationName)
{
    return simulationName + (GMX_DOUBLE ? "_d" : "_s");
}

//! Simulation name contributes to reference data. Warnings and checked terms do not.
static const auto sc_referenceConfigFormatters = std::make_tuple(
        formatFreeEnergySimulationName,
        [](MaxNumWarnings) { return std::string{}; },
        [](const ListOfInteractionsToTest&) { return std::string{}; });
static const auto sc_executionModeFormatters = std::make_tuple();

static const NameOfTestFromTuple<FreeEnergyReferenceTestHelper::DynamicParameters> sc_referenceTestNamer =
        FreeEnergyReferenceTestHelper::testNamer(sc_referenceConfigFormatters, sc_executionModeFormatters);

const std::vector<FreeEnergyReferenceInputConfig> sc_freeEnergyReferenceConfigs = {
    { "coulandvdwsequential_coul",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda } },
    { "coulandvdwsequential_vdw",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda } },
    { "coulandvdwtogether", MaxNumWarnings(1), { InteractionFunction::dVremainingdLambda } },
    { "coulandvdwtogether-net-charge", MaxNumWarnings(2), { InteractionFunction::dVremainingdLambda } },
    { "coulandvdwtogether-decouple-counter-charge", MaxNumWarnings(2), { InteractionFunction::dVremainingdLambda } },
    { "expanded",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda } },
    // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
    { "relative",
      MaxNumWarnings(11),
      { InteractionFunction::dVremainingdLambda,
        InteractionFunction::dVCoulombdLambda,
        InteractionFunction::dVvanderWaalsdLambda,
        InteractionFunction::dVbondeddLambda } },
    // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
    { "relative-position-restraints",
      MaxNumWarnings(11),
      { InteractionFunction::dVremainingdLambda,
        InteractionFunction::dVCoulombdLambda,
        InteractionFunction::dVvanderWaalsdLambda,
        InteractionFunction::dVbondeddLambda,
        InteractionFunction::dVrestraintdLambda } },
    { "restraints", MaxNumWarnings(1), { InteractionFunction::dVrestraintdLambda } },
    { "simtemp", MaxNumWarnings(1), {} },
    { "transformAtoB", MaxNumWarnings(1), { InteractionFunction::dVremainingdLambda } },
    { "vdwalone", MaxNumWarnings(1), { InteractionFunction::dVremainingdLambda } },
    { "cmap-perturbation", // end-to-end test for CMAP FEP
      MaxNumWarnings(0),
      { InteractionFunction::dVremainingdLambda } }
};

//! Forces-only trajectory comparison for the free-energy reference test.
TrajectoryComparison freeEnergyTrajectoryComparison()
{
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ false,
                                                          false,
                                                          false,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    trajectoryTolerances.forces = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 6.0e-5 : 5.0e-4);
    return { trajectoryMatchSettings, trajectoryTolerances };
}

void addFreeEnergySkipReasons(MessageStringCollector& skipReasons)
{
    // Reproducibility checks keep a tight tolerance, so rank count stays small. See also #3741.
    constexpr int maxNumRanks       = 8;
    const int     numRanksAvailable = getNumberOfTestMpiRanks();
    skipReasons.appendIf(numRanksAvailable > maxNumRanks,
                         formatString("Rank count is %d, but these tests support at most %d ranks.",
                                      numRanksAvailable,
                                      maxNumRanks));
    bool ciEnvIsDefined = std::getenv("GITLAB_CI") != nullptr;
    skipReasons.appendIf(ciEnvIsDefined && GMX_GPU_OPENCL,
                         "Skipping tests with OpenCL in CI, as compilation there is too slow");
}

//! `-nbfe` follows the hardware context.
std::vector<SimulationOptionTuple> makeFreeEnergyMdrunOptions(bool isGpuTest)
{
    return { SimulationOptionTuple{ "-nbfe", isGpuTest ? "gpu" : "cpu" } };
}

/*! \brief Compare a free-energy simulation to reference data on one hardware context.
 *
 * \c SimulationRunner still needs the mdrun fixture's communicator and hardware
 * detection, so suite setup delegates to \c MdrunTestFixtureBase.
 */
class FreeEnergyReferenceTest : public HardwareTestFixture<FreeEnergyReferenceTestHelper>
{
public:
    static void SetUpTestSuite() { MdrunTestFixtureBase::SetUpTestSuite(); }
    static void TearDownTestSuite() { MdrunTestFixtureBase::TearDownTestSuite(); }

    ~FreeEnergyReferenceTest() override
    {
#if GMX_LIB_MPI
        MPI_Barrier(MdrunTestFixtureBase::s_communicator);
#endif
    }

protected:
    FreeEnergyReferenceTest() :
        HardwareTestFixture(sc_referenceConfigFormatters), runner_(&fileManager_)
    {
        checkTestNameLength();
    }

    void addCustomSkipReasons(MessageStringCollector& skipReasons) override
    {
        addFreeEnergySkipReasons(skipReasons);
        const auto& simulationName = std::get<0>(GetParam());
        // Expanded ensemble simulations are not implemented on GPUs yet.
        skipReasons.appendIf(isGpuTest() && simulationName == "expanded",
                             "Expanded ensemble simulations are not implemented on GPUs");
    }

    //! Manages temporary files during the test.
    TestFileManager fileManager_;
    //! Helper object to manage the preparation for and call of mdrun
    SimulationRunner runner_;
};

TEST_P(FreeEnergyReferenceTest, WithinTolerances)
{
    auto [simulationName, maxNumWarnings, interactionsList, _] = GetParam();

    SCOPED_TRACE(formatString("Comparing FEP simulation '%s' to reference on %s",
                              simulationName.c_str(),
                              hardwareContext()->description().c_str()));

    // Tolerance set to pass with identical code version and a range of different test setups for most tests
    const auto defaultEnergyTolerance = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 5e-6 : 5e-5);
    // Some simulations are significantly longer, so they need a larger tolerance
    const auto longEnergyTolerance = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 3e-5 : 2e-4);
    const bool isLongSimulation = (simulationName == "expanded");
    const auto energyTolerance  = isLongSimulation ? longEnergyTolerance : defaultEnergyTolerance;

    EnergyTermsToCompare energyTermsToCompare{
        { interaction_function[InteractionFunction::PotentialEnergy].longname, energyTolerance }
    };
    for (const auto& interaction : interactionsList)
    {
        energyTermsToCompare.emplace(interaction_function[interaction].longname, energyTolerance);
    }

    const TrajectoryComparison trajectoryComparison = freeEnergyTrajectoryComparison();

    // Set simulation file names
    auto simulationTrajectoryFileName = fileManager_.getTemporaryFilePath("trajectory.trr");
    auto simulationEdrFileName        = fileManager_.getTemporaryFilePath("energy.edr");
    auto simulationDhdlFileName       = fileManager_.getTemporaryFilePath("dhdl.xvg");

    // Run grompp
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("sim.tpr").string();
    runner_.useTopGroAndMdpFromFepTestDatabase(simulationName);
    runner_.setMaxWarn(maxNumWarnings);
    runGrompp(&runner_);

    // Do mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulationTrajectoryFileName.string();
    runner_.edrFileName_                     = simulationEdrFileName.string();
    runner_.dhdlFileName_                    = simulationDhdlFileName.string();

    runMdrun(&runner_, makeFreeEnergyMdrunOptions(isGpuTest()));

    /* Currently used tests write trajectory (x/v/f) frames every 20 steps.
     * Except for the expanded ensemble test, all tests run for 20 steps total,
     * so the trajectory has a frame at step 0 and one at the final step.
     *
     * The forces in the first frame are the initial force evaluation at step 0,
     * which is reproducible in all precisions and rank counts (and on the GPU),
     * so that frame is always checked. This also covers the FEP force path that
     * was previously tested by a separate single-step test.
     *
     * The forces in the final frame have accumulated integration differences and
     * are only reproducible in double precision using a single rank, so that
     * second frame is only checked there.
     *
     * Note that this only concerns trajectory frames; energy frames are checked
     * in all cases. */
    const bool alsoCheckFinalForces = (GMX_DOUBLE && (getNumberOfTestMpiRanks() == 1));
    const MaxNumFrames numForceFramesToCheck = alsoCheckFinalForces ? MaxNumFrames(2) : MaxNumFrames(1);

    // Compare simulation results. Hardware is omitted from the reference-data name.
    TestReferenceChecker& rootChecker = checker();
    // Check that the energies agree with the refdata within tolerance.
    checkEnergiesAgainstReferenceData(simulationEdrFileName.string(), energyTermsToCompare, &rootChecker);
    // Check that the forces agree with the refdata within tolerance.
    checkTrajectoryAgainstReferenceData(
            simulationTrajectoryFileName, trajectoryComparison, &rootChecker, numForceFramesToCheck);
    if (File::exists(simulationDhdlFileName, File::returnFalseOnError))
    {
        TextInputFile dhdlFile(simulationDhdlFileName);
        auto          settings = XvgMatchSettings();
        settings.tolerance     = defaultEnergyTolerance;
        checkXvgFile(&dhdlFile, &rootChecker, settings);
    }
}

INSTANTIATE_TEST_SUITE_P(EquivalentToReference,
                         FreeEnergyReferenceTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_freeEnergyReferenceConfigs),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::NonbondedFE))),
                                 flattenTupleWithHardwareContext<FreeEnergyReferenceInputConfig>()),
                         sc_referenceTestNamer);

} // namespace
} // namespace gmx::test
