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
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"
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

/*! \brief Test fixture base for free energy calculations
 *
 * This test ensures that selected free energy perturbation calculations produce
 * results identical to an earlier version. The results of this earlier version
 * have been verified manually to ensure physical correctness.
 */
using MaxNumWarnings                = int;
using ListOfInteractionsToTest      = std::vector<int>;
using FreeEnergyReferenceTestParams = std::tuple<std::string, MaxNumWarnings, ListOfInteractionsToTest>;
class FreeEnergyReferenceTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<FreeEnergyReferenceTestParams>
{
public:
    struct PrintParametersToString
    {
        template<class ParamType>
        std::string operator()(const testing::TestParamInfo<ParamType>& parameter) const
        {
            auto simulationName = std::get<0>(parameter.param);
            std::replace(simulationName.begin(), simulationName.end(), '-', '_');
            return simulationName + (GMX_DOUBLE ? "_d" : "_s");
        }
    };
};

TEST_P(FreeEnergyReferenceTest, WithinTolerances)
{
    checkTestNameLength();
    const auto& simulationName   = std::get<0>(GetParam());
    const auto  maxNumWarnings   = std::get<1>(GetParam());
    const auto& interactionsList = std::get<2>(GetParam());

    // As these tests check reproducibility, we restrict the maximum number
    // of ranks to allow us to keep the tolerances tight. See also #3741.
    const int     numRanksAvailable = getNumberOfTestMpiRanks();
    constexpr int maxNumRanks       = 8;
    if (numRanksAvailable > maxNumRanks)
    {
        fprintf(stdout,
                "The FEP tests cannot run with %d ranks.\n"
                "The maximum number of ranks supported is %d.",
                numRanksAvailable,
                maxNumRanks);
        return;
    }

    SCOPED_TRACE(formatString("Comparing FEP simulation '%s' to reference", simulationName.c_str()));

    // Tolerance set to pass with identical code version and a range of different test setups for most tests
    const auto defaultEnergyTolerance = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 5e-6 : 5e-5);
    // Some simulations are significantly longer, so they need a larger tolerance
    const auto longEnergyTolerance = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 2e-5 : 2e-4);
    const bool isLongSimulation    = (simulationName == "expanded");
    const auto energyTolerance = isLongSimulation ? longEnergyTolerance : defaultEnergyTolerance;

    EnergyTermsToCompare energyTermsToCompare{ { interaction_function[F_EPOT].longname, energyTolerance } };
    for (const auto& interaction : interactionsList)
    {
        energyTermsToCompare.emplace(interaction_function[interaction].longname, energyTolerance);
    }

    // Specify how trajectory frame matching must work (only testing forces).
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ false,
                                                          false,
                                                          false,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    trajectoryTolerances.forces = relativeToleranceAsFloatingPoint(100.0, GMX_DOUBLE ? 5.0e-5 : 5.0e-4);

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

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
    runMdrun(&runner_);

    /* Currently used tests write trajectory (x/v/f) frames every 20 steps.
     * Except for the expanded ensemble test, all tests run for 20 steps total.
     * As the tolerances are relatively strict, we need to restrict the number of
     * force frames we can expect to match.
     * Testing more than the first force frame is only feasible in double precision
     * using a single rank.
     * Testing one force frame is only feasible in double precision.
     * Note that this only concerns trajectory frames, energy frames are checked
     * in all cases. */
    const bool testTwoTrajectoryFrames = (GMX_DOUBLE && (getNumberOfTestMpiRanks() == 1));
    const bool testOneTrajectoryFrame  = GMX_DOUBLE;

    // Compare simulation results
    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());
    // Check that the energies agree with the refdata within tolerance.
    checkEnergiesAgainstReferenceData(simulationEdrFileName.string(), energyTermsToCompare, &rootChecker);
    // Check that the trajectories agree with the refdata within tolerance.
    if (testTwoTrajectoryFrames)
    {
        checkTrajectoryAgainstReferenceData(
                simulationTrajectoryFileName.string(), trajectoryComparison, &rootChecker, MaxNumFrames(2));
    }
    else if (testOneTrajectoryFrame)
    {
        checkTrajectoryAgainstReferenceData(
                simulationTrajectoryFileName.string(), trajectoryComparison, &rootChecker, MaxNumFrames(1));
    }
    else
    {
        checkTrajectoryAgainstReferenceData(
                simulationTrajectoryFileName.string(), trajectoryComparison, &rootChecker, MaxNumFrames(0));
    }
    if (File::exists(simulationDhdlFileName, File::returnFalseOnError))
    {
        TextInputFile dhdlFile(simulationDhdlFileName);
        auto          settings = XvgMatchSettings();
        settings.tolerance     = defaultEnergyTolerance;
        checkXvgFile(&dhdlFile, &rootChecker, settings);
    }
}

// TODO: The time for OpenCL kernel compilation means these tests time
//       out. Once that compilation is cached for the whole process, these
//       tests can run in such configurations.
#if !GMX_GPU_OPENCL
INSTANTIATE_TEST_SUITE_P(
        EquivalentToReference,
        FreeEnergyReferenceTest,
        ::testing::Values(
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_coul",
                                               MaxNumWarnings(1),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_vdw",
                                               MaxNumWarnings(1),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether", MaxNumWarnings(1), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether-net-charge", MaxNumWarnings(2), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether-decouple-counter-charge",
                                               MaxNumWarnings(1),
                                               { F_DVDL } },
                FreeEnergyReferenceTestParams{ "expanded", MaxNumWarnings(1), { F_DVDL_COUL, F_DVDL_VDW } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{ "relative",
                                               MaxNumWarnings(11),
                                               { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{
                        "relative-position-restraints",
                        MaxNumWarnings(11),
                        { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "restraints", MaxNumWarnings(1), { F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "simtemp", MaxNumWarnings(1), {} },
                FreeEnergyReferenceTestParams{ "transformAtoB", MaxNumWarnings(1), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "vdwalone", MaxNumWarnings(1), { F_DVDL } }),
        FreeEnergyReferenceTest::PrintParametersToString());
#else
INSTANTIATE_TEST_SUITE_P(
        DISABLED_EquivalentToReference,
        FreeEnergyReferenceTest,
        ::testing::Values(
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_coul",
                                               MaxNumWarnings(1),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_vdw",
                                               MaxNumWarnings(1),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether", MaxNumWarnings(1), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether-net-charge", MaxNumWarnings(2), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether-decouple-counter-charge",
                                               MaxNumWarnings(1),
                                               { F_DVDL } },
                FreeEnergyReferenceTestParams{ "expanded", MaxNumWarnings(1), { F_DVDL_COUL, F_DVDL_VDW } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{ "relative",
                                               MaxNumWarnings(11),
                                               { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{
                        "relative-position-restraints",
                        MaxNumWarnings(11),
                        { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "restraints", MaxNumWarnings(1), { F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "simtemp", MaxNumWarnings(1), {} },
                FreeEnergyReferenceTestParams{ "transformAtoB", MaxNumWarnings(1), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "vdwalone", MaxNumWarnings(1), { F_DVDL } }),
        FreeEnergyReferenceTest::PrintParametersToString());
#endif

} // namespace
} // namespace gmx::test
