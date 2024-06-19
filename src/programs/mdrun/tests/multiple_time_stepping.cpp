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
 * Tests to compare that multiple time stepping is (nearly) identical to normal integration.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <filesystem>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "programs/mdrun/tests/comparison_helpers.h"
#include "programs/mdrun/tests/energycomparison.h"
#include "programs/mdrun/tests/trajectorycomparison.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief Test fixture base for two integration schemes
 *
 * This test ensures that integration with(out) different multiple time stepping
 * schemes (called via different mdp options) yield near identical energies,
 * forces and virial at step 0 and similar energies and virial after 4 steps.
 */
using MtsComparisonTestParams = std::tuple<std::string, std::string>;
class MtsComparisonTest : public MdrunTestFixture, public ::testing::WithParamInterface<MtsComparisonTestParams>
{
};

//! Returns set of energy terms to compare with associated tolerances
EnergyTermsToCompare energyTermsToCompare(const real energyTol, const real virialTol)
{
    return EnergyTermsToCompare{ { { interaction_function[F_EPOT].longname,
                                     relativeToleranceAsFloatingPoint(100.0, energyTol) },
                                   { "Vir-XX", relativeToleranceAsFloatingPoint(30.0, virialTol) },
                                   { "Vir-YY", relativeToleranceAsFloatingPoint(30.0, virialTol) },
                                   { "Vir-ZZ", relativeToleranceAsFloatingPoint(30.0, virialTol) } } };
}

TEST_P(MtsComparisonTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto mtsScheme      = std::get<1>(params);

    // Note that there should be no relevant limitation on MPI ranks and OpenMP threads
    SCOPED_TRACE(formatString(
            "Comparing for '%s' no MTS with MTS scheme '%s'", simulationName.c_str(), mtsScheme.c_str()));

    const bool isPullTest = (mtsScheme.find("pull") != std::string::npos);

    const int numSteps         = 4;
    auto      sharedMdpOptions = gmx::formatString(
            "integrator   = md\n"
            "dt           = 0.001\n"
            "nsteps       = %d\n"
            "verlet-buffer-tolerance = -1\n"
            "rlist        = 1.0\n"
            "coulomb-type = %s\n"
            "vdw-type     = cut-off\n"
            "rcoulomb     = 0.9\n"
            "rvdw         = 0.9\n"
            "constraints  = h-bonds\n",
            numSteps,
            isPullTest ? "reaction-field" : "PME");

    if (isPullTest)
    {
        sharedMdpOptions +=
                "pull                 = yes\n"
                "pull-ngroups         = 2\n"
                "pull-group1-name     = FirstWaterMolecule\n"
                "pull-group2-name     = SecondWaterMolecule\n"
                "pull-ncoords         = 1\n"
                "pull-coord1-type     = umbrella\n"
                "pull-coord1-geometry = distance\n"
                "pull-coord1-groups   = 1 2\n"
                "pull-coord1-init     = 1\n"
                "pull-coord1-k        = 10000\n";
    }

    // set nstfout to > numSteps so we only write forces at step 0
    const int nstfout       = 2 * numSteps;
    auto      refMdpOptions = sharedMdpOptions
                         + gmx::formatString(
                                 "mts       = no\n"
                                 "nstcalcenergy = %d\n"
                                 "nstenergy = %d\n"
                                 "nstxout   = 0\n"
                                 "nstvout   = 0\n"
                                 "nstfout   = %d\n",
                                 numSteps,
                                 numSteps,
                                 nstfout);

    auto mtsMdpOptions = sharedMdpOptions
                         + gmx::formatString(
                                 "mts        = yes\n"
                                 "mts-levels = 2\n"
                                 "mts-level2-forces = %s\n"
                                 "mts-level2-factor = 2\n"
                                 "nstcalcenergy = %d\n"
                                 "nstenergy  = %d\n"
                                 "nstxout    = 0\n"
                                 "nstvout    = 0\n"
                                 "nstfout    = %d\n",
                                 mtsScheme.c_str(),
                                 numSteps,
                                 numSteps,
                                 nstfout);

    // At step 0 the energy and virial should only differ due to rounding errors
    EnergyTermsToCompare energyTermsToCompareStep0 = energyTermsToCompare(0.001, 0.01);
    EnergyTermsToCompare energyTermsToCompareAllSteps =
            energyTermsToCompare(mtsScheme == "pme" ? 0.015 : 0.04, mtsScheme == "pme" ? 0.1 : 0.2);

    // Specify how trajectory frame matching must work.
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,
                                                          true,
                                                          true,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::NoComparison,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    // Tolerances for force comparison needs to be somewhat higher than the current default
    // when comparing between forces CPU and GPU, due to summation order differences
    trajectoryTolerances.forces = relativeToleranceAsFloatingPoint(1000.0, GMX_DOUBLE ? 1.0e-7 : 1.0e-4);

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

    // Set file names
    auto simulator1TrajectoryFileName = fileManager_.getTemporaryFilePath("sim1.trr");
    auto simulator1EdrFileName        = fileManager_.getTemporaryFilePath("sim1.edr");
    auto simulator2TrajectoryFileName = fileManager_.getTemporaryFilePath("sim2.trr");
    auto simulator2EdrFileName        = fileManager_.getTemporaryFilePath("sim2.edr");

    // Run grompp
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("sim.tpr").string();
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(refMdpOptions);
    runGrompp(&runner_);

    // Do first mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator1TrajectoryFileName.string();
    runner_.edrFileName_                     = simulator1EdrFileName.string();
    runMdrun(&runner_);

    runner_.useStringAsMdpFile(mtsMdpOptions);
    runGrompp(&runner_);

    // Do second mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator2TrajectoryFileName.string();
    runner_.edrFileName_                     = simulator2EdrFileName.string();
    runMdrun(&runner_);

    // Compare simulation results at step 0, which should be identical
    compareEnergies(simulator1EdrFileName.string(),
                    simulator2EdrFileName.string(),
                    energyTermsToCompareStep0,
                    MaxNumFrames(1));
    compareTrajectories(simulator1TrajectoryFileName.string(),
                        simulator2TrajectoryFileName.string(),
                        trajectoryComparison);

    // Compare energies at the last step (and step 0 again) with lower tolerance
    compareEnergies(simulator1EdrFileName.string(),
                    simulator2EdrFileName.string(),
                    energyTermsToCompareAllSteps,
                    MaxNumFrames::compareAllFrames());
}

INSTANTIATE_TEST_SUITE_P(
        MultipleTimeSteppingIsNearSingleTimeStepping,
        MtsComparisonTest,
        ::testing::Combine(::testing::Values("ala"),
                           ::testing::Values("longrange-nonbonded",
                                             "longrange-nonbonded nonbonded pair dihedral")));

INSTANTIATE_TEST_SUITE_P(MultipleTimeSteppingIsNearSingleTimeSteppingPull,
                         MtsComparisonTest,
                         ::testing::Combine(::testing::Values("spc2"), ::testing::Values("pull")));

} // namespace
} // namespace test
} // namespace gmx
