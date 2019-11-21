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
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 */

#include <vector>

#include "gmxpre.h"

#include "gromacs/nblib/box.h"
#include "gromacs/nblib/setup.h"
#include "gromacs/nblib/topology.h"
#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/simulationstate.h"

#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

class SimulationStateTester
{
public:
    std::vector<gmx::RVec> coords_;
    std::vector<gmx::RVec> vel_;

    Box box_;
    TopologyBuilder topologyBuilder_;

    SimulationStateTester() : box_(2.05449)
    {
        constexpr int NumArgonAtoms = 3;

        AtomType argonAtom("AR", 39.94800, 0.0, 0.0062647225, 9.847044e-06);

        MoleculeType argonMolecule("AR");
        argonMolecule.addAtom("AR", argonAtom);

        TopologyBuilder topBuilder;
        topBuilder.addMolecule(argonMolecule, NumArgonAtoms);

        Topology top = topBuilder.buildTopology();

        Box box(7.73950);

        std::vector<gmx::RVec> coords = {
            { 5.158, 6.923, 3.413 },
            { 2.891, 6.634, 0.759 },
            { 4.356, 2.932, 1.414 },
        };

        SimulationState simulationState(coords, box, top);

        EXPECT_EQ(top.getMasses().size(), NumArgonAtoms);
    }

    void setupSimulationState()
    {
        auto topology = topologyBuilder_.buildTopology();
        SimulationState(coords_, box_, topology, vel_);
    }
};

TEST(NBlibTest, SimulationStateArgonBox)
{
    SimulationStateTester simulationStateTester;
    EXPECT_NO_THROW(simulationStateTester.setupSimulationState());
}

TEST(NBlibTest, SimulationStateArgonBoxCoordThrowNAN)
{
    SimulationStateTester simulationStateTester;
    simulationStateTester.coords_[2][0] = NAN;
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

TEST(NBlibTest, SimulationStateArgonBoxCoordThrowINF)
{

    SimulationStateTester simulationStateTester;
    simulationStateTester.coords_[2][0] = INFINITY;
    EXPECT_THROW(simulationStateTester.setupSimulationState(), gmx::InvalidInputError);
}

}  // namespace
}  // namespace test
}  // namespace nblib
