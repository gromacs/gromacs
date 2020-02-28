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
 * This implements molecule setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/nblib/integrator.h"

#include "gromacs/nblib/molecules.h"
#include "gromacs/nblib/particletype.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

#include "testsystems.h"


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

    ParticleType particleType(ParticleName("H"), Mass(1.0));
    Molecule     molecule("SomeMolecule");
    molecule.addParticle(ParticleName("SomeAtom"), particleType);

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(molecule, numAtoms);
    Topology topology = topologyBuilder.buildTopology();

    std::vector<gmx::RVec> x(numAtoms);
    std::vector<gmx::RVec> v(numAtoms);
    std::vector<gmx::RVec> f(numAtoms);

    f[0][XX] = 1.0;
    f[0][YY] = 2.0;
    f[0][ZZ] = 0.0;

    Box box(100);

    std::vector<gmx::RVec> x0(x);
    std::vector<gmx::RVec> v0(v);

    SimulationState simulationState(x, v, f, box, topology);

    LeapFrog integrator(simulationState);

    gmx::test::FloatingPointTolerance tolerance = gmx::test::absoluteTolerance(numSteps * 0.000005);
    for (int step = 0; step < numSteps; step++)
    {
        real totalTime = step * dt;

        gmx::RVec xAnalytical;
        gmx::RVec vAnalytical;

        for (int i = 0; i < numAtoms; i++)
        {
            for (int d = 0; d < DIM; d++)
            {
                // Analytical solution for constant-force particle movement
                int  typeIndex = simulationState.topology().getParticleTypeIdOfAllParticles()[i];
                real im = 1.0 / simulationState.topology().getParticleTypes()[typeIndex].mass();
                xAnalytical[d] =
                        x0[i][d] + v0[i][d] * totalTime + 0.5 * f[i][d] * totalTime * totalTime * im;
                vAnalytical[d] = v0[i][d] + f[i][d] * totalTime * im;

                EXPECT_REAL_EQ_TOL(xAnalytical[d], simulationState.coordinates()[i][d], tolerance)
                        << gmx::formatString(
                                   "Coordinate %d of atom %d is different from analytical solution "
                                   "at step %d.",
                                   d, i, step);

                EXPECT_REAL_EQ_TOL(vAnalytical[d], simulationState.velocities()[i][d], tolerance)
                        << gmx::formatString(
                                   "Velocity component %d of atom %d is different from analytical "
                                   "solution at step %d.",
                                   d, i, step);
            }
            integrator.integrate(dt);
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
