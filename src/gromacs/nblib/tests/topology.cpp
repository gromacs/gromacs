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
#include "gmxpre.h"

#include "gromacs/nblib/topology.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/nblib/atomtype.h"
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

class TwoWaterMolecules
{
public:
    TwoWaterMolecules() :
        Ow(AtomName("Ow"), Mass(16), C6(6.), C12(12.)),
        Hw(AtomName("Hw"), Mass(1), C6(0.6), C12(0.12))
    {
    }

    Topology buildTopology()
    {
        //! Manually Create Molecule (Water)

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

        //! Build the topology
        TopologyBuilder topologyBuilder;

        //! Add some molecules to the topology
        topologyBuilder.addMolecule(water, 2);
        Topology topology = topologyBuilder.buildTopology();
        return topology;
    }

    int numAtoms = 6;

    //! Define Atom Type
    AtomType Ow;
    AtomType Hw;
};

TEST(NBlibTest, TopologyHasNumAtoms)
{
    TwoWaterMolecules waters;
    Topology          watersTopology = waters.buildTopology();
    const int         test           = watersTopology.numAtoms();
    const int         ref            = 6;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasCharges)
{
    TwoWaterMolecules       waters;
    Topology                watersTopology = waters.buildTopology();
    const std::vector<real> test           = watersTopology.getCharges();
    const std::vector<real> ref            = { -0.6, 0.3, 0.3, -0.6, 0.3, 0.3 };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasMasses)
{
    TwoWaterMolecules       waters;
    Topology                watersTopology = waters.buildTopology();

    const std::vector<real> test = expandQuantity(watersTopology, &AtomType::mass);
    const std::vector<real> ref            = { 16., 1., 1., 16., 1., 1. };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasAtomTypes)
{
    TwoWaterMolecules           waters;
    Topology                    watersTopology = waters.buildTopology();
    const std::vector<AtomType> test           = watersTopology.getAtomTypes();
    const std::vector<AtomType> ref            = { waters.Ow, waters.Hw };
    const std::vector<AtomType> ref2           = { waters.Hw, waters.Ow };
    EXPECT_TRUE(ref == test || ref2 == test);
}

TEST(NBlibTest, TopologyHasAtomTypeIds)
{
    TwoWaterMolecules waters;
    Topology          watersTopology = waters.buildTopology();

    const std::vector<int>      testIds   = watersTopology.getAtomTypeIdOfAllAtoms();
    const std::vector<AtomType> testTypes = watersTopology.getAtomTypes();

    std::vector<AtomType> testTypesExpanded;
    for (int i : testIds)
    {
        testTypesExpanded.push_back(testTypes[i]);
    }

    const std::vector<AtomType> ref = { waters.Ow, waters.Hw, waters.Hw,
                                        waters.Ow, waters.Hw, waters.Hw };

    EXPECT_TRUE(ref == testTypesExpanded);
}

TEST(NBlibTest, TopologyThrowsIdenticalAtomType)
{
    //! User error: Two different AtomTypes with the same name
    AtomType U235(AtomName("Uranium"), Mass(235), C6(6.), C12(12.));
    AtomType U238(AtomName("Uranium"), Mass(238), C6(6.), C12(12.));

    Molecule ud235("UraniumDimer235");
    ud235.addAtom(AtomName("U1"), U235);
    ud235.addAtom(AtomName("U2"), U235);

    Molecule ud238("UraniumDimer238");
    ud238.addAtom(AtomName("U1"), U238);
    ud238.addAtom(AtomName("U2"), U238);

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(ud235, 1);
    EXPECT_THROW(topologyBuilder.addMolecule(ud238, 1), gmx::InvalidInputError);
}

TEST(NBlibTest, TopologyHasExclusions)
{
    TwoWaterMolecules     waters;
    Topology              watersTopology = waters.buildTopology();
    gmx::ListOfLists<int> testExclusions = watersTopology.getGmxExclusions();

    const std::vector<std::vector<int>> refExclusions = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
                                                          { 3, 4, 5 }, { 3, 4, 5 }, { 3, 4, 5 } };

    compareLists(testExclusions, refExclusions);
}

TEST(NBlibTest, TopologyHasNonbondedParameters)
{
    TwoWaterMolecules                         waters;
    Topology                                  watersTopology = waters.buildTopology();

    const std::vector<real> testC6 = expandQuantity(watersTopology, &AtomType::c6);
    const std::vector<real> testC12 = expandQuantity(watersTopology, &AtomType::c12);

    const std::vector<real> refC6 = {6, 0.6, 0.6, 6, 0.6, 0.6};
    const std::vector<real> refC12 = {12, 0.12, 0.12, 12, 0.12, 0.12};

    EXPECT_EQ(refC6, testC6);
    EXPECT_EQ(refC12, testC12);
}

//! Todo: this belongs to ForceCalculator
//TEST(NBlibTest, TopologyHasAtomInfoAllVdw)
//{
//    TwoWaterMolecules      waters;
//    Topology               watersTopology = waters.buildTopology();
//    const std::vector<int> test           = watersTopology.getAtomInfoAllVdw();
//    std::vector<int>       ref;
//    ref.resize(watersTopology.numAtoms());
//    for (size_t atomI = 0; atomI < ref.size(); atomI++)
//    {
//        SET_CGINFO_HAS_VDW(ref[atomI]);
//    }
//    EXPECT_EQ(ref, test);
//}

TEST(NBlibTest, toGmxExclusionBlockWorks)
{
    std::vector<std::tuple<int, int>> testInput{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 },
                                                 { 1, 2 }, { 2, 0 }, { 2, 1 }, { 2, 2 } };

    std::vector<gmx::ExclusionBlock> reference;

    gmx::ExclusionBlock localBlock;
    localBlock.atomNumber.push_back(0);
    localBlock.atomNumber.push_back(1);
    localBlock.atomNumber.push_back(2);

    reference.push_back(localBlock);
    reference.push_back(localBlock);
    reference.push_back(localBlock);

    std::vector<gmx::ExclusionBlock> probe = detail::toGmxExclusionBlock(testInput);

    ASSERT_EQ(reference.size(), probe.size());
    for (size_t i = 0; i < reference.size(); ++i)
    {
        ASSERT_EQ(reference[i].atomNumber.size(), probe[i].atomNumber.size());
        for (size_t j = 0; j < reference[i].atomNumber.size(); ++j)
        {
            EXPECT_EQ(reference[i].atomNumber[j], probe[i].atomNumber[j]);
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
