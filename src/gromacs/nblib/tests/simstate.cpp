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
 * This implements SimulationState tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include <vector>
#include <stdio.h>

#include "gromacs/nblib/atomtype.h"
#include "gromacs/nblib/box.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/nbkerneloptions.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/topology.h"

#include "testutils/testasserts.h"

#include "testsystems.h"

namespace nblib
{
namespace test
{
namespace
{

static void compareValues(const std::vector<gmx::RVec>& ref, const std::vector<gmx::RVec>& test)
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
    SimulationStateTester simulationStateTester;
    EXPECT_NO_THROW(simulationStateTester.setupSimulationState());
}

TEST(NBlibTest, SimulationStateThrowsCoordinateNAN)
{
    SimulationStateTester simulationStateTester;
    simulationStateTester.setCoordinate(2, 0, NAN);
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsCoordinateINF)
{
    SimulationStateTester simulationStateTester;
    simulationStateTester.setCoordinate(2, 0, INFINITY);
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsVelocityNAN)
{
    SimulationStateTester simulationStateTester;
    simulationStateTester.setVelocity(2, 0, NAN);
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateThrowsVelocityINF)
{
    SimulationStateTester simulationStateTester;
    simulationStateTester.setVelocity(2, 0, INFINITY);
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateCanMove)
{
    SimulationStateTester simulationStateTester;
    SimulationState       simState = simulationStateTester.setupSimulationState();
    EXPECT_NO_THROW(SimulationState movedSimState = std::move(simState));
}

TEST(NBlibTest, SimulationStateCanAssign)
{
    SimulationStateTester simulationStateTester;
    SimulationState       simState = simulationStateTester.setupSimulationState();
    EXPECT_NO_THROW(const SimulationState& gmx_unused AssignedSimState = simState);
}

TEST(NBlibTest, SimulationStateHasTopology)
{
    SimulationStateTester simulationStateTester;
    SimulationState       simState = simulationStateTester.setupSimulationState();
    const Topology&       testTop  = simState.topology();
    const Topology&       refTop   = simulationStateTester.topology();
    EXPECT_EQ(refTop, testTop);
}

TEST(NBlibTest, SimulationStateHasBox)
{
    SimulationStateTester simulationStateTester;
    SimulationState       simState = simulationStateTester.setupSimulationState();
    const Box&            testBox  = simState.box();
    const Box&            refBox   = simulationStateTester.box();
    EXPECT_EQ(refBox, testBox);
}

TEST(NBlibTest, SimulationStateHasCorrectCoordinates)
{
    SimulationStateTester  simulationStateTester;
    SimulationState        simState = simulationStateTester.setupSimulationState();
    std::vector<gmx::RVec> test     = simState.coordinates();
    std::vector<gmx::RVec> ref      = simulationStateTester.coordinates();
    compareValues(ref, test);
}

TEST(NBlibTest, SimulationStateHasCorrectVelocities)
{
    SimulationStateTester  simulationStateTester;
    SimulationState        simState = simulationStateTester.setupSimulationState();
    std::vector<gmx::RVec> test     = simState.velocities();
    std::vector<gmx::RVec> ref      = simulationStateTester.velocities();
    compareValues(ref, test);
}

} // namespace
} // namespace test
} // namespace nblib
