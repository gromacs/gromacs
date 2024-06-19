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
 * \brief Tests for checkpoint writing sanity checks
 *
 * Checks that final checkpoint is equal to final trajectory output.
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

#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/tests/comparison_helpers.h"
#include "programs/mdrun/tests/trajectorycomparison.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx::test
{
namespace
{

class CheckpointCoordinatesSanityChecks :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<std::tuple<std::string, std::string, std::string, std::string>>
{
public:
    void runSimulation(MdpFieldValues mdpFieldValues, int numSteps)
    {
        mdpFieldValues["nsteps"] = toString(numSteps);
        // Trajectories have the initial and the last frame
        mdpFieldValues["nstxout"] = toString(numSteps);
        mdpFieldValues["nstvout"] = toString(numSteps);
        mdpFieldValues["nstfout"] = toString(0);

        // Run grompp
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        runGrompp(&runner_);

        // Do first mdrun
        runMdrun(&runner_, {});
    }

    static void compareCptAndTrr(const std::string&          trrFileName,
                                 const std::string&          cptFileName,
                                 const TrajectoryComparison& trajectoryComparison)
    {
        TrajectoryFrameReader trrReader(trrFileName);
        TrajectoryFrameReader cptReader(cptFileName);
        // Checkpoint has at least one frame
        EXPECT_TRUE(cptReader.readNextFrame());
        // Trajectory has at least two frames
        EXPECT_TRUE(trrReader.readNextFrame());
        EXPECT_NO_THROW(trrReader.frame());
        EXPECT_TRUE(trrReader.readNextFrame());

        // Now compare frames
        trajectoryComparison(cptReader.frame(), trrReader.frame());

        // Files had exactly 1 / 2 frames
        EXPECT_FALSE(cptReader.readNextFrame());
        EXPECT_FALSE(trrReader.readNextFrame());
    }
};

TEST_P(CheckpointCoordinatesSanityChecks, WithinTolerances)
{
    const auto& params              = GetParam();
    const auto& simulationName      = std::get<0>(params);
    const auto& integrator          = std::get<1>(params);
    const auto& temperatureCoupling = std::get<2>(params);
    const auto& pressureCoupling    = std::get<3>(params);

    // Specify how trajectory frame matching must work.
    TrajectoryFrameMatchSettings trajectoryMatchSettings{ true,
                                                          true,
                                                          true,
                                                          ComparisonConditions::MustCompare,
                                                          ComparisonConditions::MustCompare,
                                                          ComparisonConditions::NoComparison,
                                                          MaxNumFrames::compareAllFrames() };
    if (integrator == "md-vv")
    {
        // When using md-vv and modular simulator, the velocities are expected to be off by
        // 1/2 dt between checkpoint (top of the loop) and trajectory (full time step state)
        trajectoryMatchSettings.velocitiesComparison = ComparisonConditions::NoComparison;
    }
    const TrajectoryTolerances trajectoryTolerances{
        defaultRealTolerance(), defaultRealTolerance(), defaultRealTolerance(), defaultRealTolerance()
    };

    const auto mdpFieldValues =
            prepareMdpFieldValues(simulationName, integrator, temperatureCoupling, pressureCoupling);
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    // Set file names
    const auto trrFileName = fileManager_.getTemporaryFilePath(".trr");

    SCOPED_TRACE(formatString(
            "Checking the sanity of the checkpointed coordinates using system '%s' "
            "with integrator '%s', '%s' temperature coupling, and '%s' pressure coupling ",
            simulationName.c_str(),
            integrator.c_str(),
            temperatureCoupling.c_str(),
            pressureCoupling.c_str()));

    SCOPED_TRACE("End of trajectory sanity");
    // Running a few steps - we expect the checkpoint to be equal
    // to the final configuration
    runSimulation(mdpFieldValues, 16);
    compareCptAndTrr(runner_.fullPrecisionTrajectoryFileName_,
                     runner_.cptOutputFileName_,
                     { trajectoryMatchSettings, trajectoryTolerances });
}

#if !GMX_GPU_OPENCL
INSTANTIATE_TEST_SUITE_P(CheckpointCoordinatesAreSane,
                         CheckpointCoordinatesSanityChecks,
                         ::testing::Combine(::testing::Values("spc2"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("no"),
                                            ::testing::Values("no")));
#else
INSTANTIATE_TEST_SUITE_P(DISABLED_CheckpointCoordinatesAreSane,
                         CheckpointCoordinatesSanityChecks,
                         ::testing::Combine(::testing::Values("spc2"),
                                            ::testing::Values("md", "md-vv"),
                                            ::testing::Values("no"),
                                            ::testing::Values("no")));
#endif

} // namespace
} // namespace gmx::test
