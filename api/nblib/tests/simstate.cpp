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
 * This implements SimulationState tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <vector>

#include "testutils/testasserts.h"

#include "nblib/box.h"
#include "nblib/exception.h"
#include "nblib/simulationstate.h"
#include "nblib/topology.h"

#include "simulationstateimpl.h"
#include "testhelpers.h"
#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

//! Utility function to compare 2 std::vectors of gmx::RVec used to compare cartesians
void compareValues(const std::vector<Vec3>& ref, const std::vector<Vec3>& test)
{
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (int j = 0; j < dimSize; j++)
        {
            EXPECT_EQ(ref[i][j], test.at(i)[j]);
        }
    }
}

TEST(NBlibTest, CanConstructSimulationState)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    EXPECT_NO_THROW(argonSimulationStateBuilder.setupSimulationState());
}

TEST(NBlibTest, SimulationStateThrowsCoordinateNAN)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    argonSimulationStateBuilder.setCoordinate(2, 0, NAN);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), InputException);
}

TEST(NBlibTest, SimulationStateThrowsCoordinateINF)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    argonSimulationStateBuilder.setCoordinate(2, 0, INFINITY);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), InputException);
}

TEST(NBlibTest, SimulationStateThrowsVelocityNAN)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    argonSimulationStateBuilder.setVelocity(2, 0, NAN);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), InputException);
}

TEST(NBlibTest, SimulationStateThrowsVelocityINF)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    argonSimulationStateBuilder.setVelocity(2, 0, INFINITY);
    EXPECT_THROW(argonSimulationStateBuilder.setupSimulationState(), InputException);
}

TEST(NBlibTest, SimulationStateCanMove)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    EXPECT_NO_THROW(SimulationState movedSimState = std::move(simState));
}

TEST(NBlibTest, SimulationStateCanAssign)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    EXPECT_NO_THROW(const SimulationState& gmx_unused AssignedSimState = simState);
}

TEST(NBlibTest, SimulationStateHasBox)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    const Box&                  testBox  = simState.box();
    const Box&                  refBox   = argonSimulationStateBuilder.box();
    // GTEST does not like the comparison operator in a different namespace
    const bool compare = (refBox == testBox);
    EXPECT_TRUE(compare);
}

TEST(NBlibTest, SimulationStateHasCorrectCoordinates)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    std::vector<Vec3>           test     = simState.coordinates();
    std::vector<Vec3>           ref      = argonSimulationStateBuilder.coordinates();
    compareValues(ref, test);
}

TEST(NBlibTest, SimulationStateHasCorrectVelocities)
{
    ArgonSimulationStateBuilder argonSimulationStateBuilder(fftypes::GROMOS43A1);
    SimulationState             simState = argonSimulationStateBuilder.setupSimulationState();
    std::vector<Vec3>           test     = simState.velocities();
    std::vector<Vec3>           ref      = argonSimulationStateBuilder.velocities();
    compareValues(ref, test);
}

} // namespace
} // namespace test
} // namespace nblib
