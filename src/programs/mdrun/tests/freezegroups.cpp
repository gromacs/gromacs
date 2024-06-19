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
 *        containing freeze groups
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <array>
#include <string>
#include <tuple>
#include <utility>
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
/*! \brief Test fixture checking sanity of freeze group results
 *
 * This tests the sanity of simulation results containing fully and partially
 * frozen atoms. For fully frozen atoms, it checks that their reported position
 * is identical for all steps, and that their velocity is zero. For partially
 * frozen atoms (for simplicity only in z-direction), it checks that their
 * position is identical in the frozen dimension for all steps, and that their
 * velocity is zero in the frozen dimension.
 */
using FreezeGroupTestParams = std::tuple<std::string, std::string, std::string>;
class FreezeGroupTest : public MdrunTestFixture, public ::testing::WithParamInterface<FreezeGroupTestParams>
{
public:
    //! Check that the frozen positions don't change and velocities are zero
    static void checkFreezeGroups(const std::string&            trajectoryName,
                                  ArrayRef<const unsigned int>  fullyFrozenAtoms,
                                  ArrayRef<const unsigned int>  partiallyFrozenAtomsDimZ,
                                  const FloatingPointTolerance& positionsTolerance,
                                  const FloatingPointTolerance& velocitiesTolerance)
    {
        auto [fullyFrozenPositions, fullyFrozenVelocities] =
                getFrozenPositionsAndVelocities(trajectoryName, fullyFrozenAtoms);
        auto [partiallyFrozenPositions, partiallyFrozenVelocities] =
                getFrozenPositionsAndVelocities(trajectoryName, partiallyFrozenAtomsDimZ);
        GMX_RELEASE_ASSERT(fullyFrozenPositions.size() == fullyFrozenVelocities.size(),
                           "Position and velocity trajectory don't have the same length.");
        GMX_RELEASE_ASSERT(partiallyFrozenPositions.size() == partiallyFrozenVelocities.size(),
                           "Position and velocity trajectory don't have the same length.");
        GMX_RELEASE_ASSERT(fullyFrozenPositions.size() == partiallyFrozenPositions.size(),
                           "Fully and partially frozen trajectory don't have the same length.");
        const auto trajectorySize = fullyFrozenPositions.size();

        for (auto frameIdx = decltype(trajectorySize){ 0 }; frameIdx < trajectorySize; frameIdx++)
        {
            SCOPED_TRACE(formatString("Checking frame %lu", frameIdx + 1));
            if (frameIdx > 0)
            {
                checkFullyFrozenPositions(fullyFrozenPositions[frameIdx],
                                          fullyFrozenPositions[frameIdx - 1],
                                          positionsTolerance);
                checkZDimFrozenPositions(partiallyFrozenPositions[frameIdx],
                                         partiallyFrozenPositions[frameIdx - 1],
                                         positionsTolerance);
            }
            checkFullyFrozenVelocities(fullyFrozenVelocities[frameIdx], velocitiesTolerance);
            checkZDimFrozenVelocities(partiallyFrozenVelocities[frameIdx], velocitiesTolerance);
        }
    }

    //! Check that fully frozen frame velocities are zero
    static void checkFullyFrozenVelocities(ArrayRef<const RVec>          velocities,
                                           const FloatingPointTolerance& velocitiesTolerance)
    {
        SCOPED_TRACE("Checking fully frozen velocity frame");
        std::vector<RVec> zeroVelocities(velocities.size(), RVec{ 0, 0, 0 });
        EXPECT_THAT(zeroVelocities, Pointwise(RVecEq(velocitiesTolerance), velocities));
    }
    //! Check that z-dimension frozen frame velocities are zero
    static void checkZDimFrozenVelocities(ArrayRef<const RVec>          velocities,
                                          const FloatingPointTolerance& velocitiesTolerance)
    {
        SCOPED_TRACE("Checking z-dimension frozen velocity frame");
        std::vector<real> zVelocities;
        for (const auto& v : velocities)
        {
            zVelocities.emplace_back(v[ZZ]);
        }
        std::vector<real> zeroVelocities(zVelocities.size(), 0);
        EXPECT_THAT(zeroVelocities, Pointwise(RealEq(velocitiesTolerance), zVelocities));
    }
    //! Check that fully frozen frame positions are static
    static void checkFullyFrozenPositions(ArrayRef<const RVec>          positions,
                                          ArrayRef<const RVec>          previousPositions,
                                          const FloatingPointTolerance& positionsTolerance)
    {
        SCOPED_TRACE("Checking fully frozen position frame");
        EXPECT_THAT(previousPositions, Pointwise(RVecEq(positionsTolerance), positions));
    }
    //! Check that z-dimension frozen frame positions are zero
    static void checkZDimFrozenPositions(ArrayRef<const RVec>          positions,
                                         ArrayRef<const RVec>          previousPositions,
                                         const FloatingPointTolerance& positionsTolerance)
    {
        SCOPED_TRACE("Checking z-dimension frozen position frame");
        std::vector<real> zPositions;
        for (const auto& p : positions)
        {
            zPositions.emplace_back(p[ZZ]);
        }
        std::vector<real> zPrevPositions;
        for (const auto& p : previousPositions)
        {
            zPrevPositions.emplace_back(p[ZZ]);
        }
        EXPECT_THAT(zPrevPositions, Pointwise(RealEq(positionsTolerance), zPositions));
    }

    static std::tuple<std::vector<std::vector<RVec>>, std::vector<std::vector<RVec>>>
    getFrozenPositionsAndVelocities(const std::string& trajectoryName, ArrayRef<const unsigned int> frozenAtoms)
    {
        std::vector<std::vector<RVec>> positions;
        std::vector<std::vector<RVec>> velocities;

        TrajectoryFrameReader trajectoryFrameReader(trajectoryName);
        while (trajectoryFrameReader.readNextFrame())
        {
            const auto frame = trajectoryFrameReader.frame();
            positions.emplace_back();
            velocities.emplace_back();
            for (const auto& index : frozenAtoms)
            {
                positions.back().emplace_back(frame.x().at(index));
                velocities.back().emplace_back(frame.v().at(index));
            }
        }

        return { std::move(positions), std::move(velocities) };
    }
};

TEST_P(FreezeGroupTest, WithinTolerances)
{
    const auto& params         = GetParam();
    const auto& integrator     = std::get<0>(params);
    const auto& tcoupling      = std::get<1>(params);
    const auto& pcoupling      = std::get<2>(params);
    const auto& simulationName = "alanine_vacuo";

    constexpr std::array<unsigned int, 5> backbone   = { 4, 6, 8, 14, 16 };
    constexpr std::array<unsigned int, 3> sideChainH = { 0, 10, 18 };

    if (integrator == "md-vv" && pcoupling == "parrinello-rahman")
    {
        GTEST_SKIP() << "Parrinello-Rahman is not implemented in md-vv.";
    }

    // Prepare mdp input
    auto mdpFieldValues = prepareMdpFieldValues(simulationName, integrator, tcoupling, pcoupling);
    mdpFieldValues["nsteps"]      = "8";
    mdpFieldValues["nstxout"]     = "4";
    mdpFieldValues["nstvout"]     = "4";
    mdpFieldValues["freezegrps"]  = "Backbone SideChain-H";
    mdpFieldValues["freezedim"]   = "Y Y Y N N Y";
    mdpFieldValues["constraints"] = "h-bonds";

    if (pcoupling == "c-rescale" && tcoupling == "no" && integrator != "sd" && integrator != "bd")
    {
        mdpFieldValues["ensemble-temperature-setting"] = "constant";
        mdpFieldValues["ensemble-temperature"]         = "0";
    }

    // Run grompp
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    // Allow one warning for COMM removal + partially frozen atoms, and one for berendsen t-coupling
    int maxWarn = 1;
    if (pcoupling == "berendsen")
    {
        maxWarn++;
    }
    runner_.setMaxWarn(maxWarn);
    runGrompp(&runner_);
    // Run mdrun
    runMdrun(&runner_);

    // Check frozen atoms
    FloatingPointTolerance positionsTolerance  = ulpTolerance(0);
    FloatingPointTolerance velocitiesTolerance = ulpTolerance(0);
    checkFreezeGroups(
            runner_.fullPrecisionTrajectoryFileName_, backbone, sideChainH, positionsTolerance, velocitiesTolerance);
}

INSTANTIATE_TEST_SUITE_P(
        FreezeWorks,
        FreezeGroupTest,
        ::testing::Combine(::testing::Values("md", "md-vv", "sd", "bd"),
                           ::testing::Values("no"),
                           ::testing::Values("no", "c-rescale", "parrinello-rahman")));
} // namespace
} // namespace gmx::test
