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
 * This implements SimulationState tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include <vector>

#include "gromacs/nblib/box.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/topology.h"

#include "testutils/testasserts.h"

#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

void compareValues(const std::vector<gmx::RVec>& ref, const std::vector<gmx::RVec>& test)
{
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            EXPECT_EQ(ref[i][j], test.at(i)[j]);
        }
    }
}

TEST(NBlibTest, CanConstructSimulationState)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    EXPECT_NO_THROW(argonSimulationStateBuilder.setupSimulationState());
}

TEST(NBlibTest, SimulationStateThrowsCoordinateNAN)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    argonSimulationStateBuilder.setCoordinate(2, 0, NAN);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsCoordinateINF)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    argonSimulationStateBuilder.setCoordinate(2, 0, INFINITY);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsVelocityNAN)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    argonSimulationStateBuilder.setVelocity(2, 0, NAN);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsVelocityINF)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    argonSimulationStateBuilder.setVelocity(2, 0, INFINITY);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateCanMove)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    EXPECT_NO_THROW(SimulationState movedSimState = std::move(simState));
}

TEST(NBlibTest, SimulationStateCanAssign)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    EXPECT_NO_THROW(const SimulationState& gmx_unused AssignedSimState = simState);
}

TEST(NBlibTest, SimulationStateHasBox)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    const Box&                  testBox  = simState.box();
    const Box&                  refBox   = argonSimulationStateBuilder.box();
    // GTEST does not like the comparison operator in a different namespace
    const bool compare = (refBox == testBox);
    EXPECT_TRUE(compare);
}

TEST(NBlibTest, SimulationStateHasCorrectCoordinates)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    std::vector<gmx::RVec>      test     = simState.coordinates();
    std::vector<gmx::RVec>      ref      = argonSimulationStateBuilder.coordinates();
    compareValues(ref, test);
}

TEST(NBlibTest, SimulationStateHasCorrectVelocities)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder;
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    std::vector<gmx::RVec>      test     = simState.velocities();
    std::vector<gmx::RVec>      ref      = argonSimulationStateBuilder.velocities();
    compareValues(ref, test);
}

} // namespace
} // namespace test
} // namespace nblib