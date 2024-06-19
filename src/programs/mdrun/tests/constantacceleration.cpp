/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief End-to-end tests checking sanity of results of simulations
 *        containing constant acceleration groups
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstddef>

#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"
#include "testutils/trajectoryreader.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx::test
{
namespace
{
/*! \brief Test fixture checking the velocities of aomts
 *
 * This tests that velocities of non-accelerated atoms are zero
 * and that velocities of accelarated atoms match acceleration*time.
 * This code assumes the first half the atoms are non-accelerated
 * and the second half are accelerated
 */
using AccelerationGroupTestParams = std::tuple<std::string, std::string>;
class AccelerationGroupTest :
    public MdrunTestFixture,
    public ::testing::WithParamInterface<AccelerationGroupTestParams>
{
public:
    //! Check that the velocities are zero or accelerated
    static void checkVelocities(const std::string&            trajectoryName,
                                const RVec                    acceleration,
                                const FloatingPointTolerance& tolerance)
    {
        const size_t c_groupSize = 3;

        const std::vector<RVec> zeroVelocities(c_groupSize, RVec{ 0, 0, 0 });

        TrajectoryFrameReader trajectoryFrameReader(trajectoryName);
        while (trajectoryFrameReader.readNextFrame())
        {
            const auto frame = trajectoryFrameReader.frame();
            GMX_RELEASE_ASSERT(frame.v().size() == 2 * c_groupSize,
                               "Expect velocities for both atom groups");

            const RVec              referenceVelocity = real(frame.time()) * acceleration;
            const std::vector<RVec> referenceVelocities(c_groupSize, referenceVelocity);

            SCOPED_TRACE("Checking velocities of non-accelerated atoms");
            ArrayRef<const RVec> nonAcceleratedVelocities = frame.v().subArray(0, c_groupSize);
            EXPECT_THAT(nonAcceleratedVelocities, Pointwise(RVecEq(tolerance), zeroVelocities));

            SCOPED_TRACE("Checking velocities of accelerated atoms");
            ArrayRef<const RVec> acceleratedVelocities = frame.v().subArray(c_groupSize, c_groupSize);
            EXPECT_THAT(acceleratedVelocities, Pointwise(RVecEq(tolerance), referenceVelocities));
        }
    }
};

TEST_P(AccelerationGroupTest, WithinTolerances)
{
    const auto& params         = GetParam();
    const auto& integrator     = std::get<0>(params);
    const auto& tcoupling      = std::get<1>(params);
    const auto& simulationName = "spc2";

    // Prepare mdp input
    auto mdpFieldValues       = prepareMdpFieldValues(simulationName, integrator, tcoupling, "no");
    mdpFieldValues["nsteps"]  = "8";
    mdpFieldValues["dt"]      = "0.002";
    mdpFieldValues["nstxout"] = "0";
    mdpFieldValues["nstvout"] = "8";
    mdpFieldValues["nstfout"] = "0";
    mdpFieldValues["comm-mode"] = "none";
    // The two groups will not see each other when the cut-off is 0.9 nm
    mdpFieldValues["coulombtype"]             = "reaction-field";
    mdpFieldValues["rcoulomb"]                = "0.8";
    mdpFieldValues["rvdw"]                    = "0.8";
    mdpFieldValues["verlet-buffer-tolerance"] = "-1";
    mdpFieldValues["rlist"]                   = "0.9";
    // Couple the (non-)accelerated groups separately, so their velocties are independent
    mdpFieldValues["tc-grps"] = "firstWaterMolecule secondWaterMolecule";
    mdpFieldValues["ref-t"]   = "0.001 0.001";
    // Use weak coupling so we can check vecolities of atoms with a tight tolerance
    mdpFieldValues["tau-t"]      = "10.0 10.0";
    const RVec c_acceleration    = { 2.0, 3.0, 1.5 };
    mdpFieldValues["acc-grps"]   = "secondWaterMolecule";
    mdpFieldValues["accelerate"] = "2.0 3.0 1.5";
    // Set initial velocities to zero
    mdpFieldValues["gen-vel"]  = "yes";
    mdpFieldValues["gen-temp"] = "0";

    // Run grompp
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    runGrompp(&runner_);
    // Run mdrun
    runMdrun(&runner_);

    // T-coupling causes changes in the velocities up to 1e-4
    auto tolerance = absoluteTolerance((GMX_DOUBLE && tcoupling == "no") ? 1e-10 : 1e-4);

    // Check the velocities of the non-accelerated and accelerated groups
    checkVelocities(runner_.fullPrecisionTrajectoryFileName_, c_acceleration, tolerance);
}

// The v-rescale case check that Ekin computation and temperature coupling
// can be performed independently for atoms groups, so the accelerations
// are not affected. This can be useful in practice.
INSTANTIATE_TEST_SUITE_P(AccelerationWorks,
                         AccelerationGroupTest,
                         ::testing::Combine(::testing::Values("md", "md-vv"),
                                            ::testing::Values("no", "v-rescale")));
} // namespace
} // namespace gmx::test
