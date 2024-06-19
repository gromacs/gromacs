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
 * This implements molecule setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/integrator.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/testasserts.h"

#include "nblib/interactions.h"
#include "nblib/molecules.h"
#include "nblib/particletype.h"
#include "nblib/simulationstate.h"
#include "nblib/topology.h"
#include "nblib/util/util.hpp"
#include "nblib/vector.h"

namespace nblib
{
namespace test
{
namespace
{

TEST(NBlibTest, IntegratorWorks)
{
    int  numAtoms = 1;
    int  numSteps = 100;
    real dt       = 0.001;

    ParticleType particleType(ParticleTypeName("H"), Mass(1.0));
    Molecule     molecule(MoleculeName("SomeMolecule"));
    molecule.addParticle(ParticleName("SomeAtom"), particleType);

    ParticleTypesInteractions interactions;
    interactions.add(particleType.name(), C6{ 0 }, C12{ 0 });

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(molecule, numAtoms);
    topologyBuilder.addParticleTypesInteractions(interactions);
    Topology topology = topologyBuilder.buildTopology();

    std::vector<Vec3> x(numAtoms, { 0.0, 0.0, 0.0 });
    std::vector<Vec3> v(numAtoms, { 0.0, 0.0, 0.0 });
    std::vector<Vec3> f(numAtoms, { 1.0, 2.0, 0.0 });

    Box box(100);

    std::vector<Vec3> x0(x);
    std::vector<Vec3> v0(v);

    SimulationState simulationState(x, v, f, box, topology);
    put_atoms_in_box(PbcType::Xyz, box.legacyMatrix(), x0);

    LeapFrog integrator(simulationState.topology(), simulationState.box());

    gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(numSteps * 0.000005);
    for (int step = 0; step < numSteps; step++)
    {
        real totalTime = step * dt;

        Vec3 xAnalytical;
        Vec3 vAnalytical;

        for (int i = 0; i < numAtoms; i++)
        {
            for (int d = 0; d < dimSize; d++)
            {
                // Analytical solution for constant-force particle movement
                int  typeIndex = simulationState.topology().getParticleTypeIdOfAllParticles()[i];
                real im = 1.0 / simulationState.topology().getParticleTypes()[typeIndex].mass();
                xAnalytical[d] =
                        x0[i][d] + v0[i][d] * totalTime + 0.5 * f[i][d] * totalTime * totalTime * im;
                vAnalytical[d] = v0[i][d] + f[i][d] * totalTime * im;

                EXPECT_REAL_EQ_TOL(xAnalytical[d], simulationState.coordinates()[i][d], tolerance)
                        << formatString(
                                   "Coordinate {} of atom {} is different from analytical solution "
                                   "at step {}.",
                                   d,
                                   i,
                                   step);

                EXPECT_REAL_EQ_TOL(vAnalytical[d], simulationState.velocities()[i][d], tolerance)
                        << formatString(
                                   "Velocity component {} of atom {} is different from analytical "
                                   "solution at step {}.",
                                   d,
                                   i,
                                   step);
            }
            integrator.integrate(dt,
                                 simulationState.coordinates(),
                                 simulationState.velocities(),
                                 simulationState.forces());
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
