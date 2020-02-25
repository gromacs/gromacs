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
 * This implements topology setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gromacs/nblib/forcecalculator.h"

#include "gromacs/nblib/particletype.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/topology.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

// This is defined in src/gromacs/mdtypes/forcerec.h but there is also a
// legacy C6 macro defined there that conflicts with the nblib C6 type.
// Todo: Once that C6 has been refactored into a regular function, this
//       file can just include forcerec.h
#define SET_CGINFO_HAS_VDW(cgi) (cgi) = ((cgi) | (1 << 23))

class KernelSystemTester
{
public:
    std::vector<gmx::RVec> coordinates;
    std::vector<gmx::RVec> velocities;

    Box             box;
    TopologyBuilder topologyBuilder;

    KernelSystemTester() : box(7.25449)
    {
        constexpr int numWaters = 2;

        //! Define Particle Type
        ParticleType Ow(ParticleName("Ow"), Mass(16), C6(6.), C12(12.));
        ParticleType Hw(ParticleName("Hw"), Mass(1), C6(0.6), C12(0.12));

        //! Define Molecule
        Molecule water("water");

        //! Add the particles
        water.addParticle(ParticleName("Oxygen"), Charge(-0.6), Ow);
        water.addParticle(ParticleName("H1"), Charge(+0.3), Hw);
        water.addParticle(ParticleName("H2"), Charge(+0.3), Hw);

        //! Add the exclusions
        water.addExclusion("Oxygen", "H1");
        water.addExclusion("Oxygen", "H2");
        water.addExclusion("H1", "H2");

        // Todo: Add bonds functionality so this can be used/tested
        // water.addHarmonicBond(HarmonicType{1, 2, "H1", "Oxygen"});
        // water.addHarmonicBond(HarmonicType{1, 2, "H2", "Oxygen"});

        //! Add some molecules to the topology
        topologyBuilder.addMolecule(water, numWaters);

        coordinates = {
            { 0.569, 1.275, 1.165 }, { 0.476, 1.268, 1.128 }, { 0.580, 1.364, 1.209 },
            { 1.555, 1.511, 0.703 }, { 1.498, 1.495, 0.784 }, { 1.496, 1.521, 0.623 },
        };

        velocities = {
            { 0.569, 1.215, 1.965 }, { 0.669, 1.225, 1.865 }, { 0.769, 1.235, 1.765 },
            { 0.869, 1.245, 1.665 }, { 0.169, 0.275, 1.565 }, { 0.269, 2.275, 1.465 },
        };
    }

    NBKernelSystem setupKernelSystem()
    {
        Topology topology       = topologyBuilder.buildTopology();
        auto     simState       = SimulationState(coordinates, box, topology, velocities);
        auto     nbKernelSystem = NBKernelSystem(simState);
        return nbKernelSystem;
    }
};

/*
TEST(NBlibTest, KernelSystemHasNumParticles)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const int test = kernelSystem.numParticles;
    const int ref  = 6;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasNonbondedParameters)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test = kernelSystem.nonbondedParameters;
    std::vector<real> ref;
    ref.resize(kernelSystem.numParticles*kernelSystem.numParticles*2, 0);
    ref[0] = 6;
    ref[1] = 12;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasParticlesTypes)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<int> test = kernelSystem.particleTypes;
    std::vector<int> ref;
    ref.resize(kernelSystem.numParticles, 0);
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasCharges)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test = kernelSystem.charges;
    const std::vector<real> ref  = { -0.6, 0.3, 0.3, -0.6, 0.3, 0.3 };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasMasses)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test = kernelSystem.masses;
    const std::vector<real> ref  = { 16., 1., 1., 16., 1., 1. };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasParticleInfoAllVdw)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<int> test = kernelSystem.particleInfoAllVdw;
    std::vector<int> ref;
    ref.resize(kernelSystem.numParticles);
    for (size_t particleI = 0; particleI < ref.size(); particleI++)
    {
        SET_CGINFO_HAS_VDW(ref[particleI]);
    }
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasCoordinates)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<gmx::RVec> test = kernelSystem.coordinates;
    const std::vector<gmx::RVec> ref  = kernelSystemTester.coordinates;
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        EXPECT_EQ(ref[i][j], test[i][j]);
    }
}

TEST(NBlibTest, KernelSystemHasVelocities)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<gmx::RVec> test = kernelSystem.velocities;
    const std::vector<gmx::RVec> ref  = kernelSystemTester.velocities;
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        EXPECT_EQ(ref[i][j], test[i][j]);
    }
}

TEST(NBlibTest, KernelSystemHasExclusions)
{
    KernelSystemTester kernelSystemTester;
    auto kernelSystem = kernelSystemTester.setupKernelSystem();
    t_blocka testBlocka     = kernelSystem.excls;
    std::vector<gmx::ExclusionBlock> testExclusionBlocks;

    //! Setting t_blocka.nr is needed for conversion to ExclusionBlock
    testBlocka.nr = kernelSystem.numParticles;
    testExclusionBlocks.resize(kernelSystem.numParticles);
    blockaToExclusionBlocks(&testBlocka, testExclusionBlocks);

    std::vector<std::vector<int>> refExclusionBlocks = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
                                                         { 3, 4, 5 }, { 3, 4, 5 }, { 3, 4, 5 } };
    for (size_t particle = 0; particle < refExclusionBlocks.size(); particle++)
    {
        for (size_t exclusion = 0; exclusion < refExclusionBlocks[particle].size(); exclusion++)
        {
            EXPECT_EQ(refExclusionBlocks[particle][exclusion], testExclusionBlocks[particle].particleNumber[exclusion]);
        }
    }
}
*/

} // namespace
} // namespace test
} // namespace nblib
