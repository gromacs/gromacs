/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests to compare two simulators which are expected to be identical
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/stringutil.h"

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
        std::tuple<std::tuple<std::string, std::string, std::string, std::string>, std::string>;
class SimulatorComparisonTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<SimulatorComparisonTestParams>
{
};

TEST_P(SimulatorComparisonTest, WithinTolerances)
{
    auto params         = GetParam();
    auto mdpParams      = std::get<0>(params);
    auto simulationName = std::get<0>(mdpParams);
    auto integrator     = std::get<1>(mdpParams);
    auto tcoupling      = std::get<2>(mdpParams);
    auto pcoupling      = std::get<3>(mdpParams);
    auto envVariable    = std::get<1>(params);

    if (integrator == "md-vv" && pcoupling == "Parrinello-Rahman")
    {
        // do_md calls this MTTK, requires Nose-Hoover, and
        // does not work with constraints or anisotropically
        return;
    }

    SCOPED_TRACE(formatString(
            "Comparing two simulations of '%s' "
            "with integrator '%s' and '%s' temperature coupling, "
            "switching environment variable '%s'",
            simulationName.c_str(), integrator.c_str(), tcoupling.c_str(), envVariable.c_str()));

    auto mdpFieldValues = prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(),
                                                tcoupling.c_str(), pcoupling.c_str());

    EnergyTermsToCompare energyTermsToCompare{ {
            { interaction_function[F_EPOT].longname, relativeToleranceAsPrecisionDependentUlp(10.0, 100, 80) },
            { interaction_function[F_EKIN].longname, relativeToleranceAsPrecisionDependentUlp(60.0, 100, 80) },
            { interaction_function[F_PRES].longname,
              relativeToleranceAsPrecisionDependentFloatingPoint(10.0, 0.01, 0.001) },
    } };

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
    }

    // Specify how trajectory frame matching must work.
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,
                                                          true,
                                                          true,
                                                          ComparisonConditions::MustCompare,
                                                          ComparisonConditions::MustCompare,
                                                          ComparisonConditions::MustCompare };
    TrajectoryTolerances trajectoryTolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    if (simulationName != "argon12")
    {
        trajectoryTolerances.velocities = trajectoryTolerances.coordinates;
    }

    // Build the functor that will compare reference and test
    // trajectory frames in the chosen way.
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings, trajectoryTolerances };

    int numWarningsToTolerate = 0;
    executeSimulatorComparisonTest(envVariable, &fileManager_, &runner_, simulationName,
                                   numWarningsToTolerate, mdpFieldValues, energyTermsToCompare,
                                   trajectoryComparison);
}

// TODO: The time for OpenCL kernel compilation means these tests time
//       out. Once that compilation is cached for the whole process, these
//       tests can run in such configurations.
// These tests are very sensitive, so we only run them in double precision.
// As we change call ordering, they might actually become too strict to be useful.
#if GMX_GPU != GMX_GPU_OPENCL && GMX_DOUBLE
INSTANTIATE_TEST_CASE_P(SimulatorsAreEquivalentDefaultModular,
                        SimulatorComparisonTest,
                        ::testing::Combine(::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                                              ::testing::Values("md-vv"),
                                                              ::testing::Values("no", "v-rescale"),
                                                              ::testing::Values("no")),
                                           ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_CASE_P(
        SimulatorsAreEquivalentDefaultLegacy,
        SimulatorComparisonTest,
        ::testing::Combine(::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                              ::testing::Values("md"),
                                              ::testing::Values("no", "v-rescale"),
                                              ::testing::Values("no", "Parrinello-Rahman")),
                           ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_SimulatorsAreEquivalentDefaultModular,
                        SimulatorComparisonTest,
                        ::testing::Combine(::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                                              ::testing::Values("md-vv"),
                                                              ::testing::Values("no", "v-rescale"),
                                                              ::testing::Values("no")),
                                           ::testing::Values("GMX_DISABLE_MODULAR_SIMULATOR")));
INSTANTIATE_TEST_CASE_P(
        DISABLED_SimulatorsAreEquivalentDefaultLegacy,
        SimulatorComparisonTest,
        ::testing::Combine(::testing::Combine(::testing::Values("argon12", "tip3p5"),
                                              ::testing::Values("md"),
                                              ::testing::Values("no", "v-rescale"),
                                              ::testing::Values("no", "Parrinello-Rahman")),
                           ::testing::Values("GMX_USE_MODULAR_SIMULATOR")));
#endif

} // namespace
} // namespace test
} // namespace gmx
