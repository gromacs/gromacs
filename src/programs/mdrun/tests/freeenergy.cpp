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

#include <algorithm>
#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/readinp.h"
#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/fixedcapacityvector.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

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

//! Some tests can not use GPU update due to unsupported features.
using GpuUpdateSupport = YesNoType;

/*! \brief Where the update (integration) phase runs.
 *
 * This is an execution mode: it changes how the result is computed, not the
 * result itself, so it contributes to test names but not to reference-data
 * names (both update locations compare against the same reference data). The
 * GPU update path is only exercised on GPU hardware. */
enum class FreeEnergyUpdate : int
{
    //! Run the update on the CPU (mdrun '-update cpu').
    Cpu,
    //! Run the update on the GPU (mdrun '-update gpu').
    Gpu,
    //! Signals the last enum entry.
    Count
};

//! Physical system under test. This contains what selects the reference data.
using FreeEnergyReferenceInputConfig =
        std::tuple<IntegrationAlgorithm, std::string, MaxNumWarnings, ListOfInteractionsToTest>;

//! Execution modes (how, not what): whether the update runs on the CPU or the GPU.
using FreeEnergyReferenceExecutionModes = std::tuple<FreeEnergyUpdate>;

//! The update location is an execution mode, so it is skipped from reference-data names.
using FreeEnergyReferenceTestHelper =
        HardwareAndExecutionTestHelper<FreeEnergyReferenceInputConfig, FreeEnergyReferenceExecutionModes>;

//! The GPU-update execution mode contributes a "gpu" token (the update location); CPU update
//! adds nothing. The token is kept short so that the longest test names (the SD1 variant of
//! coulandvdwtogether-decouple-counter-charge) stay within the 120-character limit.
static const EnumerationArray<FreeEnergyUpdate, const char*> sc_updateNames = { { "", "gpu" } };

//! Formatter for integrator
auto formatIntegrator = [](IntegrationAlgorithm integrationAlgorithm)
{ return formatString("%s_%s", enumValueToString(integrationAlgorithm), GMX_DOUBLE ? "d" : "s"); };

//! Simulation name and integrator contribute to reference data. Warnings and checked terms do not.
static const auto sc_referenceConfigFormatters = std::make_tuple(formatIntegrator,
                                                                 useString,
                                                                 toEmptyString<MaxNumWarnings>,
                                                                 toEmptyString<ListOfInteractionsToTest>);
//! The update location contributes to test names only, not to reference-data names.
static const auto sc_executionModeFormatters = std::make_tuple(sc_updateNames);

static const NameOfTestFromTuple<FreeEnergyReferenceTestHelper::DynamicParameters> sc_referenceTestNamer =
        FreeEnergyReferenceTestHelper::testNamer(sc_referenceConfigFormatters, sc_executionModeFormatters);

/*! \brief Describes a free-energy system and the integrators to run it with.
 *
 * Some systems are intentionally with a specific integrator:
 *  - expanded ensemble and simulated tempering require the md-vv integrator;
 */
struct FreeEnergySystem
{
    //! Name of the system in the FEP test database.
    std::string name;
    //! Maximum number of warnings tolerated by grompp.
    MaxNumWarnings maxNumWarnings;
    //! Energy terms (besides potential energy) to compare against reference data.
    ListOfInteractionsToTest interactionsToTest;
    //! Which integrator variants may be run with the GPU update.
    GpuUpdateSupport gpuUpdateSupport;
    //! List of supported integrators
    std::vector<IntegrationAlgorithm> supportedIntegrators;
};

/*! \brief The free-energy systems under test.
 *
 * Depending on the test system, we use the predefined integrator or
 * test for both MD and SD1. In addition to this, we also test GPU update
 * for those cases where it is supported by the integrator and the test setup.
 */
const std::vector<FreeEnergySystem> sc_freeEnergySystems = {
    { "coulandvdwsequential_coul",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    { "coulandvdwsequential_vdw",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    { "coulandvdwtogether",
      MaxNumWarnings(1),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    { "coulandvdwtogether-net-charge",
      MaxNumWarnings(2),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    { "coulandvdwtogether-decouple-counter-charge",
      MaxNumWarnings(2),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    // Expanded ensemble requires the md-vv integrator
    { "expanded",
      MaxNumWarnings(1),
      { InteractionFunction::dVCoulombdLambda, InteractionFunction::dVvanderWaalsdLambda },
      GpuUpdateSupport::No,
      { IntegrationAlgorithm::VV } },
    // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x).
    // Perturbs masses, so GPU update is not supported.
    { "relative",
      MaxNumWarnings(11),
      { InteractionFunction::dVremainingdLambda,
        InteractionFunction::dVCoulombdLambda,
        InteractionFunction::dVvanderWaalsdLambda,
        InteractionFunction::dVbondeddLambda },
      GpuUpdateSupport::No,
      { IntegrationAlgorithm::MD } },
    // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x).
    // Perturbs masses, so GPU update is not supported.
    { "relative-position-restraints",
      MaxNumWarnings(11),
      { InteractionFunction::dVremainingdLambda,
        InteractionFunction::dVCoulombdLambda,
        InteractionFunction::dVvanderWaalsdLambda,
        InteractionFunction::dVbondeddLambda,
        InteractionFunction::dVrestraintdLambda },
      GpuUpdateSupport::No,
      { IntegrationAlgorithm::MD } },
    { "restraints",
      MaxNumWarnings(1),
      { InteractionFunction::dVrestraintdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD, IntegrationAlgorithm::SD1 } },
    // Simulated tempering requires the md-vv integrator
    { "simtemp", MaxNumWarnings(1), {}, GpuUpdateSupport::No, { IntegrationAlgorithm::VV } },
    // Perturbs masses and uses md-vv, so GPU update is not supported.
    { "transformAtoB",
      MaxNumWarnings(1),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::No,
      { IntegrationAlgorithm::VV } },
    // Uses md-vv
    { "vdwalone",
      MaxNumWarnings(1),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::No,
      { IntegrationAlgorithm::VV } },
    // End-to-end test for CMAP FEP. No temperature coupling means no SD integrator
    { "cmap-perturbation",
      MaxNumWarnings(0),
      { InteractionFunction::dVremainingdLambda },
      GpuUpdateSupport::Yes,
      { IntegrationAlgorithm::MD } }
};

//! Look up the descriptor for a named free-energy system.
const FreeEnergySystem& findFreeEnergySystem(const std::string& name)
{
    const auto it = std::find_if(sc_freeEnergySystems.begin(),
                                 sc_freeEnergySystems.end(),
                                 [&](const FreeEnergySystem& system) { return system.name == name; });
    GMX_RELEASE_ASSERT(it != sc_freeEnergySystems.end(),
                       "Encountered an unknown free-energy system name");
    return *it;
}

//! Whether the GPU update path can be used for this system and integrator combination.
bool freeEnergyGpuUpdateIsSupported(const std::string& simulationName, const IntegrationAlgorithm testIntegrator)
{
    if (testIntegrator == IntegrationAlgorithm::VV || testIntegrator == IntegrationAlgorithm::SD1)
    {
        // no support for md-vv or sd integrator, change once SD1 supports GPU update
        return false;
    }
    return findFreeEnergySystem(simulationName).gpuUpdateSupport == GpuUpdateSupport::Yes;
}

const std::vector<FreeEnergyReferenceInputConfig> sc_freeEnergyReferenceConfigs = []()
{
    std::vector<FreeEnergyReferenceInputConfig> configs;
    for (const FreeEnergySystem& system : sc_freeEnergySystems)
    {
        for (const IntegrationAlgorithm& testIntegrator : system.supportedIntegrators)
        {
            configs.emplace_back(
                    testIntegrator, system.name, system.maxNumWarnings, system.interactionsToTest);
        }
    }
    return configs;
}();

//! Update locations to test. GPU update is skipped at run time on CPU hardware.
const std::vector<FreeEnergyUpdate> sc_freeEnergyUpdateModes = { FreeEnergyUpdate::Cpu,
                                                                 FreeEnergyUpdate::Gpu };

/*! \brief Modify the \c mdpContents in situ by setting the chosen integrator*/
void addIntegratorToMdpContents(std::string* mdpContents, const IntegrationAlgorithm integrator)
{
    mdpContents->append(formatString("integrator = %s\n", enumValueToString(integrator)));
    if (integrator == IntegrationAlgorithm::SD1)
    {
        // Fix the Langevin-dynamics seed so the stochastic run is reproducible.
        mdpContents->append("ld-seed = 1993\n");
    }
}

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

//! `-nbfe` follows the hardware context; `-update` follows the execution mode.
std::vector<SimulationOptionTuple> makeFreeEnergyMdrunOptions(bool isGpuTest, FreeEnergyUpdate updateMode)
{
    return { SimulationOptionTuple{ "-nbfe", isGpuTest ? "gpu" : "cpu" },
             SimulationOptionTuple{ "-update", updateMode == FreeEnergyUpdate::Gpu ? "gpu" : "cpu" } };
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
        const IntegrationAlgorithm testIntegrator = std::get<0>(GetParam());
        const auto&                simulationName = std::get<1>(GetParam());
        const auto                 updateMode     = std::get<4>(GetParam());
        // Expanded ensemble simulations are not implemented on GPUs yet.
        skipReasons.appendIf(isGpuTest() && simulationName == "expanded",
                             "Expanded ensemble simulations are not implemented on GPUs");
        // The GPU update path can only be exercised on GPU hardware.
        skipReasons.appendIf(updateMode == FreeEnergyUpdate::Gpu && !isGpuTest(),
                             "GPU update requires GPU hardware");
        // GPU update supports only the md integrator and does not support
        // free-energy perturbation of masses or constraints.
        skipReasons.appendIf(updateMode == FreeEnergyUpdate::Gpu && isGpuTest()
                                     && !freeEnergyGpuUpdateIsSupported(simulationName, testIntegrator),
                             "GPU update is not supported for this system and integrator");
    }

    //! Manages temporary files during the test.
    TestFileManager fileManager_;
    //! Helper object to manage the preparation for and call of mdrun
    SimulationRunner runner_;
};

TEST_P(FreeEnergyReferenceTest, WithinTolerances)
{
    auto [integrator, simulationName, maxNumWarnings, interactionsList, updateMode, _] = GetParam();

    SCOPED_TRACE(
            formatString("Comparing FEP simulation '%s' using the %s integrator with %s update to "
                         "reference on %s",
                         simulationName.c_str(),
                         enumValueToString(integrator),
                         updateMode == FreeEnergyUpdate::Gpu ? "GPU" : "CPU",
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
    // We manually add the integrator to the mdp file
    std::string mdpFile = TextReader::readFileToString(runner_.mdpFileName_);
    addIntegratorToMdpContents(&mdpFile, integrator);
    const auto fullMdpFileName = fileManager_.getTemporaryFilePath("final.mdp");
    TextWriter::writeFileFromString(fullMdpFileName, mdpFile);
    runner_.mdpFileName_ = fullMdpFileName.string();
    runner_.setMaxWarn(maxNumWarnings);
    runGrompp(&runner_);

    // Do mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulationTrajectoryFileName.string();
    runner_.edrFileName_                     = simulationEdrFileName.string();
    runner_.dhdlFileName_                    = simulationDhdlFileName.string();

    runMdrun(&runner_, makeFreeEnergyMdrunOptions(isGpuTest(), updateMode));

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
        // The dhdl derivatives diverge like the energies for the SD1 GPU update.
        settings.tolerance = defaultEnergyTolerance;
        checkXvgFile(&dhdlFile, &rootChecker, settings);
    }
}

//! Flatten tuple<InputConfig, ExecutionMode, hardware context> into the flat parameter tuple.
auto flattenFreeEnergyTupleWithHardwareContext()
{
    return [](const std::tuple<FreeEnergyReferenceInputConfig, FreeEnergyUpdate, const TestHardwareContext*>& nested)
    {
        return std::tuple_cat(std::get<0>(nested),
                              std::make_tuple(std::get<1>(nested), std::get<2>(nested)));
    };
}

INSTANTIATE_TEST_SUITE_P(EquivalentToReference,
                         FreeEnergyReferenceTest,
                         ::testing::ConvertGenerator(
                                 ::testing::Combine(::testing::ValuesIn(sc_freeEnergyReferenceConfigs),
                                                    ::testing::ValuesIn(sc_freeEnergyUpdateModes),
                                                    ::testing::ValuesIn(getHardwareContextsWithCapability(
                                                            GpuConfigurationCapabilities::NonbondedFE))),
                                 flattenFreeEnergyTupleWithHardwareContext()),
                         sc_referenceTestNamer);

} // namespace
} // namespace gmx::test
