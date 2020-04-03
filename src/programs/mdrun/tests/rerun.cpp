/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * Tests for the mdrun -rerun functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
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
/*! \brief Test fixture base for mdrun -rerun
 *
 * This test ensures mdrun can run a simulation, writing a trajectory
 * and matching energies, and reproduce the same energies from a rerun
 * to within a tight tolerance. It says nothing about whether a rerun
 * can reproduce energies from a trajectory generated with older code,
 * since that is not a useful property. Whether mdrun produced correct
 * energies then and now needs different kinds of testing, but if
 * true, this test ensures the rerun has the expected property.
 *
 * The limitations of mdrun and its output means that reproducing the
 * same energies is currently only meaningful for integration without
 * thermostats or barostats, however the present form of the test
 * infrastructure has in-principle support for such, if that is ever
 * needed/useful.
 *
 * We should also not compare pressure, because with constraints the
 * non-search steps need a much larger tolerance, and per Issue #1868
 * we should stop computing pressure in reruns anyway.
 *
 * Similarly, per 1868, in the present implementation the kinetic
 * energy quantities are not generally reproducible, either.
 *
 * The choices for tolerance are arbitrary but sufficient.  Rerun does
 * pair search every frame, so it cannot in general exactly reproduce
 * quantities from a normal run, because the accumulation order
 * differs. (Nor does it reproduce pair-search frames exactly,
 * either). */
class MdrunRerunTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string>>
{
public:
    //! Trajectory components to compare
    static const TrajectoryFrameMatchSettings trajectoryMatchSettings;
};

// Compare box, positions and forces, but not velocities
// (velocities are ignored in reruns)
const TrajectoryFrameMatchSettings MdrunRerunTest::trajectoryMatchSettings = {
    true,
    true,
    true,
    ComparisonConditions::MustCompare,
    ComparisonConditions::NoComparison,
    ComparisonConditions::MustCompare
};

TEST_P(MdrunRerunTest, WithinTolerances)
{
    auto params         = GetParam();
    auto simulationName = std::get<0>(params);
    auto integrator     = std::get<1>(params);
    SCOPED_TRACE(
            formatString("Comparing normal and rerun of simulation '%s' "
                         "with integrator '%s'",
                         simulationName.c_str(), integrator.c_str()));

    auto mdpFieldValues =
            prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no");

    // bd is much less reproducible in a rerun than the other integrators
    const int            toleranceScaleFactor = (integrator == "bd") ? 2 : 1;
    EnergyTermsToCompare energyTermsToCompare{ {
            { interaction_function[F_EPOT].longname,
              relativeToleranceAsPrecisionDependentUlp(10.0, 24 * toleranceScaleFactor,
                                                       40 * toleranceScaleFactor) },
    } };

    // Specify how trajectory frame matching must work
    TrajectoryComparison trajectoryComparison{ trajectoryMatchSettings,
                                               TrajectoryComparison::s_defaultTrajectoryTolerances };

    int numWarningsToTolerate = 0;
    executeRerunTest(&fileManager_, &runner_, simulationName, numWarningsToTolerate, mdpFieldValues,
                     energyTermsToCompare, trajectoryComparison);
}

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if GMX_GPU != GMX_GPU_OPENCL
INSTANTIATE_TEST_CASE_P(
        NormalMdrunIsReproduced,
        MdrunRerunTest,
        ::testing::Combine(::testing::Values("argon12", "tip3p5", "alanine_vsite_vacuo"),
                           ::testing::Values("md", "md-vv", "bd", "sd")));
#else
INSTANTIATE_TEST_CASE_P(
        DISABLED_NormalMdrunIsReproduced,
        MdrunRerunTest,
        ::testing::Combine(::testing::Values("argon12", "tip3p5", "alanine_vsite_vacuo"),
                           ::testing::Values("md", "md-vv", "bd", "sd")));
#endif

class MdrunRerunFreeEnergyTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string, int>>
{
};

TEST_P(MdrunRerunFreeEnergyTest, WithinTolerances)
{
    auto params          = GetParam();
    auto simulationName  = std::get<0>(params);
    auto integrator      = std::get<1>(params);
    auto initLambdaState = std::get<2>(params);
    SCOPED_TRACE(
            formatString("Comparing normal and rerun of simulation '%s' "
                         "with integrator '%s' for initial lambda state %d",
                         simulationName.c_str(), integrator.c_str(), initLambdaState));

    auto mdpFieldValues =
            prepareMdpFieldValues(simulationName.c_str(), integrator.c_str(), "no", "no");
    mdpFieldValues["other"] += formatString("\ninit-lambda-state = %d", initLambdaState);

    EnergyTermsToCompare energyTermsToCompare{
        { { interaction_function[F_EPOT].longname, relativeToleranceAsPrecisionDependentUlp(10.0, 24, 32) },
          { interaction_function[F_DVDL_COUL].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8) },
          { interaction_function[F_DVDL_VDW].longname, relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8) },
          { interaction_function[F_DVDL_BONDED].longname,
            relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8) },
          { interaction_function[F_DVDL_RESTRAINT].longname,
            relativeToleranceAsPrecisionDependentUlp(1.0, 8, 8) } }
    };

    // Specify how trajectory frame matching must work
    TrajectoryComparison trajectoryComparison{ MdrunRerunTest::trajectoryMatchSettings,
                                               TrajectoryComparison::s_defaultTrajectoryTolerances };

    // The md integrator triggers a warning for nearly decoupled
    // states, which we need to suppress. TODO sometimes?
    int numWarningsToTolerate = (integrator == "md") ? 1 : 0;
    executeRerunTest(&fileManager_, &runner_, simulationName, numWarningsToTolerate, mdpFieldValues,
                     energyTermsToCompare, trajectoryComparison);
}

// TODO The time for OpenCL kernel compilation means these tests time
// out. Once that compilation is cached for the whole process, these
// tests can run in such configurations.
#if GMX_GPU != GMX_GPU_OPENCL
INSTANTIATE_TEST_CASE_P(MdrunIsReproduced,
                        MdrunRerunFreeEnergyTest,
                        ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                           ::testing::Values("md", "md-vv", "sd"),
                                           ::testing::Range(0, 11)));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_MdrunIsReproduced,
                        MdrunRerunFreeEnergyTest,
                        ::testing::Combine(::testing::Values("nonanol_vacuo"),
                                           ::testing::Values("md", "md-vv", "sd"),
                                           ::testing::Range(0, 11)));
#endif

} // namespace
} // namespace test
} // namespace gmx
