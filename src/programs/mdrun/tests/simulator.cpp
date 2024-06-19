/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests to compare two simulators which are expected to be identical
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdio>
#include <cstdlib>

#include <filesystem>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

#include <gtest/gtest.h>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/gmxassert.h"
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

/*! \brief Test fixture base for two equivalent simulators
 *
 * This test ensures that two simulator code paths (called via different mdp
 * options and/or environment variables) yield identical coordinate, velocity,
 * box, force and energy trajectories, up to some (arbitrary) precision.
 *
 * These tests are useful to check that re-implementations of existing simulators
 * are correct, and that different code paths expected to yield identical results
 * are equivalent.
 */
using SimulatorComparisonTestParams =
        std::tuple<std::tuple<std::string, std::string, std::string, std::string, MdpParameterDatabase>, std::string>;
class SimulatorComparisonTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<SimulatorComparisonTestParams>
{
};

TEST_P(SimulatorComparisonTest, WithinTolerances)
{
    const auto& params              = GetParam();
    const auto& mdpParams           = std::get<0>(params);
    const auto& simulationName      = std::get<0>(mdpParams);
    const auto& integrator          = std::get<1>(mdpParams);
    const auto& tcoupling           = std::get<2>(mdpParams);
    const auto& pcoupling           = std::get<3>(mdpParams);
    const auto& additionalParameter = std::get<4>(mdpParams);
    const auto& environmentVariable = std::get<1>(params);

    int maxNumWarnings = 0;

    const bool hasConstraints = (simulationName != "argon12");

    // TODO At some point we should also test PME-only ranks.
    const int numRanksAvailable = getNumberOfTestMpiRanks();
    if (!isNumberOfPpRanksSupported(simulationName, numRanksAvailable))
    {
        fprintf(stdout,
                "Test system '%s' cannot run with %d ranks.\n"
                "The supported numbers are: %s\n",
                simulationName.c_str(),
                numRanksAvailable,
                reportNumbersOfPpRanksSupported(simulationName).c_str());
        return;
    }

    if (integrator == "md-vv" && pcoupling == "Parrinello-Rahman")
    {
        // do_md calls this MTTK, requires Nose-Hoover, and
        // does not work with constraints or anisotropically
        return;
    }

    if (tcoupling == "berendsen")
    {
        // "Using Berendsen temperature coupling invalidates the true ensemble
        maxNumWarnings++;
    }
    if (pcoupling == "berendsen")
    {
        // "Using Berendsen pressure coupling invalidates the true ensemble for the thermostat"
        maxNumWarnings++;
    }
    if (tcoupling == "andersen" && hasConstraints)
    {
        // Constraints are not allowed with non-massive Andersen
        return;
    }

    if (tcoupling == "nose-hoover" && pcoupling == "berendsen")
    {
        if (integrator == "md-vv")
        {
            // Combination not allowed by legacy do_md.
            return;
        }
    }

    const bool systemHasConstraints = (simulationName != "argon12");
    if (pcoupling == "mttk" && (tcoupling != "nose-hoover" || systemHasConstraints))
    {
        // Legacy mttk works only with Nose-Hoover and without constraints
        return;
    }

    // We should reenable C-rescale here when it supports NPH
    if (pcoupling == "c-rescale" && tcoupling == "no")
    {
        return;
    }

    const std::string envVariableModSimOn  = "GMX_USE_MODULAR_SIMULATOR";
    const std::string envVariableModSimOff = "GMX_DISABLE_MODULAR_SIMULATOR";

    GMX_RELEASE_ASSERT(
            environmentVariable == envVariableModSimOn || environmentVariable == envVariableModSimOff,
            ("Expected tested environment variable to be " + envVariableModSimOn + " or " + envVariableModSimOff)
                    .c_str());

    const auto hasConservedField = !(tcoupling == "no" && pcoupling == "no")
                                   && !(tcoupling == "andersen-massive" || tcoupling == "andersen");

    SCOPED_TRACE(formatString(
            "Comparing two simulations of '%s' "
            "with integrator '%s', '%s' temperature coupling, and '%s' pressure coupling "
            "switching environment variable '%s'",
            simulationName.c_str(),
            integrator.c_str(),
            tcoupling.c_str(),
            pcoupling.c_str(),
            environmentVariable.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(
            simulationName.c_str(), integrator.c_str(), tcoupling.c_str(), pcoupling.c_str(), additionalParameter);
    if (tcoupling == "andersen")
    {
        // Fixes error "nstcomm must be 1, not 4 for Andersen, as velocities of
        //              atoms in coupled groups are randomized every time step"
        mdpFieldValues["nstcomm"]       = "1";
        mdpFieldValues["nstcalcenergy"] = "1";
    }

    if (pcoupling == "mttk")
    {
        // Standard parameters use compressibility of 5e-5
        // Increasing compressibility makes this test significantly more sensitive
        mdpFieldValues["compressibility"] = "1";
    }

    EnergyTermsToCompare energyTermsToCompare{ {
            { interaction_function[F_EPOT].longname, relativeToleranceAsPrecisionDependentUlp(60.0, 200, 160) },
            { interaction_function[F_EKIN].longname, relativeToleranceAsPrecisionDependentUlp(60.0, 200, 160) },
            { interaction_function[F_PRES].longname,
              relativeToleranceAsPrecisionDependentFloatingPoint(10.0, 0.01, 0.001) },
    } };
    if (hasConservedField)
    {
        energyTermsToCompare.emplace(interaction_function[F_ECONSERVED].longname,
                                     relativeToleranceAsPrecisionDependentUlp(50.0, 100, 80));
    }

    if (simulationName == "argon12")
    {
        // Without constraints, we can be more strict
        energyTermsToCompare = { {
                { interaction_function[F_EPOT].longname,
                  relativeToleranceAsPrecisionDependentUlp(10.0, 24, 80) },
                { interaction_function[F_EKIN].longname,
                  relativeToleranceAsPrecisionDependentUlp(10.0, 24, 80) },
                { interaction_function[F_PRES].longname,
                  relativeToleranceAsPrecisionDependentFloatingPoint(10.0, 0.001, 0.0001) },
        } };
        if (hasConservedField)
        {
            energyTermsToCompare.emplace(interaction_function[F_ECONSERVED].longname,
                                         relativeToleranceAsPrecisionDependentUlp(10.0, 24, 80));
        }
    }

    // Specify how trajectory frame matching must work.
    const TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,
                                                                true,
                                                                true,
                                                                ComparisonConditions::MustCompare,
                                                                ComparisonConditions::MustCompare,
                                                                ComparisonConditions::MustCompare,
                                                                MaxNumFrames::compareAllFrames() };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    if (simulationName != "argon12")
    {
        trajectoryTolerances.velocities = trajectoryTolerances.coordinates;
    }

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    const TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

    // Set file names
    const auto simulator1TrajectoryFileName = fileManager_.getTemporaryFilePath("sim1.trr");
    const auto simulator1EdrFileName        = fileManager_.getTemporaryFilePath("sim1.edr");
    const auto simulator2TrajectoryFileName = fileManager_.getTemporaryFilePath("sim2.trr");
    const auto simulator2EdrFileName        = fileManager_.getTemporaryFilePath("sim2.edr");

    // Run grompp
    runner_.tprFileName_ = fileManager_.getTemporaryFilePath("sim.tpr").string();
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    runner_.setMaxWarn(maxNumWarnings);
    runGrompp(&runner_);

    // Backup current state of both environment variables and unset them
    const char* environmentVariableBackupOn  = getenv(envVariableModSimOn.c_str());
    const char* environmentVariableBackupOff = getenv(envVariableModSimOff.c_str());
    gmxUnsetenv(envVariableModSimOn.c_str());
    gmxUnsetenv(envVariableModSimOff.c_str());

    // Do first mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator1TrajectoryFileName.string();
    runner_.edrFileName_                     = simulator1EdrFileName.string();
    runMdrun(&runner_);

    // Set tested environment variable
    const int overWriteEnvironmentVariable = 1;
    gmxSetenv(environmentVariable.c_str(), "ON", overWriteEnvironmentVariable);

    // Do second mdrun
    runner_.fullPrecisionTrajectoryFileName_ = simulator2TrajectoryFileName.string();
    runner_.edrFileName_                     = simulator2EdrFileName.string();
    runMdrun(&runner_);

    // Unset tested environment variable
    gmxUnsetenv(environmentVariable.c_str());
    // Reset both environment variables to leave further tests undisturbed
    if (environmentVariableBackupOn != nullptr)
    {
        gmxSetenv(environmentVariable.c_str(), environmentVariableBackupOn, overWriteEnvironmentVariable);
    }
    if (environmentVariableBackupOff != nullptr)
    {
        gmxSetenv(environmentVariable.c_str(), environmentVariableBackupOff, overWriteEnvironmentVariable);
    }

    // Compare simulation results
    compareEnergies(simulator1EdrFileName.string(), simulator2EdrFileName.string(), energyTermsToCompare);
    compareTrajectories(simulator1TrajectoryFileName.string(),
                        simulator2TrajectoryFileName.string(),
                        trajectoryComparison);
}

// TODO: The time for OpenCL kernel compilation means these tests time
//       out. Once that compilation is cached for the whole process, these
//       tests can run in such configurations.
// These tests are very sensitive, so we only run them in double precision.
// As we change call ordering, they might actually become too strict to be useful.
#if !GMX_GPU_OPENCL && GMX_DOUBLE
INSTANTIATE_TEST_SUITE_P(
        SimulatorsAreEquivalentDefaultModular,
        SimulatorComparisonTest,
        ::testing::Combine(
                ::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                   ::testing::Values("md-vv"),
                                   ::testing::Values("no",
                                                     "v-rescale",
                                                     "berendsen",
                                                     "nose-hoover",
                                                     "andersen-massive",
                                                     "andersen"),
                                   ::testing::Values("mttk", "no", "berendsen", "c-rescale", "mttk"),
                                   ::testing::Values(MdpParameterDatabase::Default)),
                ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(
        SimulatorsAreEquivalentDefaultLegacy,
        SimulatorComparisonTest,
        ::testing::Combine(
                ::testing::Combine(
                        ::testing::Values("argon12", "tip3p5"),
                        ::testing::Values("md"),
                        ::testing::Values("no", "v-rescale", "berendsen", "nose-hoover"),
                        ::testing::Values("no", "Parrinello-Rahman", "berendsen", "c-rescale"),
                        ::testing::Values(MdpParameterDatabase::Default)),
                ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(SimulatorsAreEquivalentDefaultModularPull,
                         SimulatorComparisonTest,
                         ::testing::Combine(::testing::Combine(::testing::Values("spc2"),
                                                               ::testing::Values("md-vv"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values(MdpParameterDatabase::Pull)),
                                            ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(SimulatorsAreEquivalentDefaultLegacyPull,
                         SimulatorComparisonTest,
                         ::testing::Combine(::testing::Combine(::testing::Values("spc2"),
                                                               ::testing::Values("md"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values(MdpParameterDatabase::Pull)),
                                            ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
#else
INSTANTIATE_TEST_SUITE_P(
        DISABLED_SimulatorsAreEquivalentDefaultModular,
        SimulatorComparisonTest,
        ::testing::Combine(
                ::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                   ::testing::Values("md-vv"),
                                   ::testing::Values("no",
                                                     "v-rescale",
                                                     "berendsen",
                                                     "andersen-massive",
                                                     "andersen",
                                                     "nose-hoover"),
                                   ::testing::Values("no", "berendsen", "c-rescale", "mttk"),
                                   ::testing::Values(MdpParameterDatabase::Default)),
                ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(
        DISABLED_SimulatorsAreEquivalentDefaultLegacy,
        SimulatorComparisonTest,
        ::testing::Combine(
                ::testing::Combine(
                        ::testing::Values("argon12", "tip3p5"),
                        ::testing::Values("md"),
                        ::testing::Values("no", "v-rescale", "berendsen", "nose-hoover"),
                        ::testing::Values("no", "Parrinello-Rahman", "berendsen", "c-rescale"),
                        ::testing::Values(MdpParameterDatabase::Default)),
                ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(DISABLED_SimulatorsAreEquivalentDefaultModularPull,
                         SimulatorComparisonTest,
                         ::testing::Combine(::testing::Combine(::testing::Values("spc2"),
                                                               ::testing::Values("md-vv"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values(MdpParameterDatabase::Pull)),
                                            ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_SUITE_P(DISABLED_SimulatorsAreEquivalentDefaultLegacyPull,
                         SimulatorComparisonTest,
                         ::testing::Combine(::testing::Combine(::testing::Values("spc2"),
                                                               ::testing::Values("md"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values("no"),
                                                               ::testing::Values(MdpParameterDatabase::Pull)),
                                            ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
#endif

} // namespace
} // namespace test
} // namespace gmx
