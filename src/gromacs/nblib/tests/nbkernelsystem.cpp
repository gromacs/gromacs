/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
#include <gromacs/nblib/forcecalculator.h>
#include "gmxpre.h"

#include "gromacs/nblib/nbkernelsystem.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/nblib/atomtype.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nblib/topology.h"
#include "gromacs/topology/exclusionblocks.h"

#include "testutils/testasserts.h"

namespace nblib
{
namespace test
{
namespace
{

using ::testing::Eq;
using ::testing::Pointwise;

//! Compares all element between two lists of lists
//! Todo: unify this with the identical function in nbkernelsystem test make this a method
//!       of ListOfLists<>
template<typename T>
void compareLists(const gmx::ListOfLists<T>& list, const std::vector<std::vector<T>>& v)
{
    ASSERT_EQ(list.size(), v.size());
    for (std::size_t i = 0; i < list.size(); i++)
    {
        ASSERT_EQ(list[i].size(), v[i].size());
        EXPECT_THAT(list[i], Pointwise(Eq(), v[i]));
    }
}

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

        //! Define Atom Type
        AtomType Ow(AtomName("Ow"), Mass(16), C6(6.), C12(12.));
        AtomType Hw(AtomName("Hw"), Mass(1), C6(0.6), C12(0.12));

        //! Define Molecule
        Molecule water("water");

        //! Add the atoms
        water.addAtom(AtomName("Oxygen"), Charge(-0.6), Ow);
        water.addAtom(AtomName("H1"), Charge(+0.3), Hw);
        water.addAtom(AtomName("H2"), Charge(+0.3), Hw);

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

    ForceCalculator setupForceCalculator()
    {
        Topology topology       = topologyBuilder.buildTopology();
        auto     simState       = SimulationState(coordinates, box, topology, velocities);
        auto     options        = NBKernelOptions();
        options.nbnxmSimd       = BenchMarkKernels::SimdNo;
        auto     forceCalculator = ForceCalculator(simState, options);
        return forceCalculator;
    }
};

TEST(NBlibTest, canRunForceCalculator)
{
    KernelSystemTester kernelSystemTester;
    auto               kernelSystem = kernelSystemTester.setupKernelSystem();
    const int          test         = kernelSystem.numAtoms;
    const int          ref          = 6;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasNonbondedParameters)
{
    KernelSystemTester      kernelSystemTester;
    auto                    kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test         = kernelSystem.nonbondedParameters;
    std::vector<real>       ref;
    ref.resize(kernelSystem.numAtoms * kernelSystem.numAtoms * 2, 0);
    ref[0] = 6;
    ref[1] = 12;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasAtomTypes)
{
    KernelSystemTester     kernelSystemTester;
    auto                   kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<int> test         = kernelSystem.atomTypes;
    std::vector<int>       ref;
    ref.resize(kernelSystem.numAtoms, 0);
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasCharges)
{
    KernelSystemTester      kernelSystemTester;
    auto                    kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test         = kernelSystem.charges;
    const std::vector<real> ref          = { -0.6, 0.3, 0.3, -0.6, 0.3, 0.3 };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasMasses)
{
    KernelSystemTester      kernelSystemTester;
    auto                    kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<real> test         = kernelSystem.masses;
    const std::vector<real> ref          = { 16., 1., 1., 16., 1., 1. };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasAtomInfoAllVdw)
{
    KernelSystemTester     kernelSystemTester;
    auto                   kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<int> test         = kernelSystem.atomInfoAllVdw;
    std::vector<int>       ref;
    ref.resize(kernelSystem.numAtoms);
    for (size_t atomI = 0; atomI < ref.size(); atomI++)
    {
        SET_CGINFO_HAS_VDW(ref[atomI]);
    }
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, KernelSystemHasCoordinates)
{
    KernelSystemTester           kernelSystemTester;
    auto                         kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<gmx::RVec> test         = kernelSystem.coordinates;
    const std::vector<gmx::RVec> ref          = kernelSystemTester.coordinates;
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
            EXPECT_EQ(ref[i][j], test[i][j]);
    }
}

TEST(NBlibTest, KernelSystemHasVelocities)
{
    KernelSystemTester           kernelSystemTester;
    auto                         kernelSystem = kernelSystemTester.setupKernelSystem();
    const std::vector<gmx::RVec> test         = kernelSystem.velocities;
    const std::vector<gmx::RVec> ref          = kernelSystemTester.velocities;
    for (size_t i = 0; i < ref.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
            EXPECT_EQ(ref[i][j], test[i][j]);
    }
}

TEST(NBlibTest, TopologyHasExclusions)
{
    KernelSystemTester    kernelSystemTester;
    auto                  kernelSystem   = kernelSystemTester.setupKernelSystem();
    gmx::ListOfLists<int> testExclusions = kernelSystem.excls;

    const std::vector<std::vector<int>> refExclusions = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
                                                          { 3, 4, 5 }, { 3, 4, 5 }, { 3, 4, 5 } };

    compareLists(testExclusions, refExclusions);
}

} // namespace
} // namespace test
} // namespace nblib
