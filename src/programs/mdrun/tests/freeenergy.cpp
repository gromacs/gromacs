/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * Tests to compare free energy simulations to reference
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/refdata.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"
#include "testutils/xvgtest.h"

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
    const auto& simulationName   = std::get<0>(GetParam());
    const auto  maxNumWarnings   = std::get<1>(GetParam());
    const auto& interactionsList = std::get<2>(GetParam());

    // TODO In similar tests, we are checking if the tests
    //      can be run with the number of MPI ranks available

    SCOPED_TRACE(formatString("Comparing FEP simulation '%s' to reference", simulationName.c_str()));

    // Tolerance set to pass with identical code version and a range of different test setups
    const auto defaultEnergyTolerance = relativeToleranceAsFloatingPoint(50.0, GMX_DOUBLE ? 1e-5 : 1e-4);

    EnergyTermsToCompare energyTermsToCompare{ { interaction_function[F_EPOT].longname,
                                                 defaultEnergyTolerance } };
    for (const auto& interaction : interactionsList)
    {
        energyTermsToCompare.emplace(interaction_function[interaction].longname, defaultEnergyTolerance);
    }

    // Specify how trajectory frame matching must work (only testing forces).
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ false,
                                                          false,
                                                          false,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

    // Set simulation file names
    auto simulationTrajectoryFileName = fileManager_.getTemporaryFilePath("trajectory.trr");
    auto simulationEdrFileName        = fileManager_.getTemporaryFilePath("energy.edr");
    auto simulationDhdlFileName       = fileManager_.getTemporaryFilePath("dhdl.xvg");

    // Run grompp
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("sim.tpr");
    runner_.useTopGroAndMdpFromFepTestDatabase(simulationName);
    runGrompp(&runner_, { SimulationOptionTuple("-maxwarn", std::to_string(maxNumWarnings)) });

    // Do mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulationTrajectoryFileName;
    runner_.edrFileName_                     = simulationEdrFileName;
    runner_.dhdlFileName_                    = simulationDhdlFileName;
    runMdrun(&runner_);

    // Currently used tests write trajectory (x/v/f) frames every 20 steps.
    // Testing more than the first force frame is only feasible in double precision
    // using a single rank.
    // Note that this only concerns trajectory frames, energy frames are checked
    // in all cases.
    const bool testAllTrajectoryFrames = (GMX_DOUBLE && (getNumberOfTestMpiRanks() == 1));

    // Compare simulation results
    TestReferenceData    refData;
    TestReferenceChecker rootChecker(refData.rootChecker());
    // Check that the energies agree with the refdata within tolerance.
    checkEnergiesAgainstReferenceData(simulationEdrFileName, energyTermsToCompare, &rootChecker);
    // Check that the trajectories agree with the refdata within tolerance.
    if (testAllTrajectoryFrames)
    {
        checkTrajectoryAgainstReferenceData(simulationTrajectoryFileName, trajectoryComparison, &rootChecker);
    }
    else
    {
        checkTrajectoryAgainstReferenceData(simulationTrajectoryFileName, trajectoryComparison,
                                            &rootChecker, MaxNumFrames(1));
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
INSTANTIATE_TEST_CASE_P(
        FreeEnergyCalculationsAreEquivalentToReference,
        FreeEnergyReferenceTest,
        ::testing::Values(
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_coul",
                                               MaxNumWarnings(0),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_vdw",
                                               MaxNumWarnings(0),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether", MaxNumWarnings(0), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "expanded", MaxNumWarnings(0), { F_DVDL_COUL, F_DVDL_VDW } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{ "relative",
                                               MaxNumWarnings(10),
                                               { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{
                        "relative-position-restraints",
                        MaxNumWarnings(10),
                        { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "restraints", MaxNumWarnings(0), { F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "simtemp", MaxNumWarnings(0), {} },
                FreeEnergyReferenceTestParams{ "transformAtoB", MaxNumWarnings(0), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "vdwalone", MaxNumWarnings(0), { F_DVDL } }),
        FreeEnergyReferenceTest::PrintParametersToString());
#else
INSTANTIATE_TEST_CASE_P(
        DISABLED_FreeEnergyCalculationsAreEquivalentToReference,
        FreeEnergyReferenceTest,
        ::testing::Values(
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_coul",
                                               MaxNumWarnings(0),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwsequential_vdw",
                                               MaxNumWarnings(0),
                                               { F_DVDL_COUL, F_DVDL_VDW } },
                FreeEnergyReferenceTestParams{ "coulandvdwtogether", MaxNumWarnings(0), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "expanded", MaxNumWarnings(0), { F_DVDL_COUL, F_DVDL_VDW } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{ "relative",
                                               MaxNumWarnings(10),
                                               { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED } },
                // Tolerated warnings: No default bonded interaction types for perturbed atoms (10x)
                FreeEnergyReferenceTestParams{
                        "relative-position-restraints",
                        MaxNumWarnings(10),
                        { F_DVDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "restraints", MaxNumWarnings(0), { F_DVDL_RESTRAINT } },
                FreeEnergyReferenceTestParams{ "simtemp", MaxNumWarnings(0), {} },
                FreeEnergyReferenceTestParams{ "transformAtoB", MaxNumWarnings(0), { F_DVDL } },
                FreeEnergyReferenceTestParams{ "vdwalone", MaxNumWarnings(0), { F_DVDL } }),
        FreeEnergyReferenceTest::PrintParametersToString());
#endif

} // namespace
} // namespace gmx::test
