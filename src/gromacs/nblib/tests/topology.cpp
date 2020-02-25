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

#include "gromacs/nblib/topology.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/nblib/particletype.h"
#include "gromacs/topology/exclusionblocks.h"

#include "testutils/testasserts.h"

#include "testsystems.h"

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

TEST(NBlibTest, TopologyHasNumParticles)
{
    WaterTopology waters;
    Topology      watersTopology = waters.buildTopology(2);
    const int     test           = watersTopology.numParticles();
    const int     ref            = 6;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasCharges)
{
    WaterTopology           waters;
    Topology                watersTopology = waters.buildTopology(2);
    const std::vector<real> test           = watersTopology.getCharges();
    const std::vector<real> ref            = { Charges.at("Ow"), Charges.at("Hw"), Charges.at("Hw"),
                                    Charges.at("Ow"), Charges.at("Hw"), Charges.at("Hw") };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasMasses)
{
    WaterTopology waters;
    Topology      watersTopology = waters.buildTopology(2);

    const Mass              refOwMass = waters.water().at("Ow").mass();
    const Mass              refHwMass = waters.water().at("H").mass();
    const std::vector<real> ref = { refOwMass, refHwMass, refHwMass, refOwMass, refHwMass, refHwMass };
    const std::vector<real> test = expandQuantity(watersTopology, &ParticleType::mass);
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasParticleTypes)
{
    WaterTopology                   waters;
    Topology                        watersTopology = waters.buildTopology(2);
    const std::vector<ParticleType> test           = watersTopology.getParticleTypes();
    const ParticleType              refOw          = waters.water().at("Ow");
    const ParticleType              refHw          = waters.water().at("H");
    const std::vector<ParticleType> ref            = { refOw, refHw };
    const std::vector<ParticleType> ref2           = { refHw, refOw };
    EXPECT_TRUE(ref == test || ref2 == test);
}

TEST(NBlibTest, TopologyHasParticleTypeIds)
{
    WaterTopology waters;
    Topology      watersTopology = waters.buildTopology(2);

    const std::vector<int>          testIds   = watersTopology.getParticleTypeIdOfAllParticles();
    const std::vector<ParticleType> testTypes = watersTopology.getParticleTypes();

    std::vector<ParticleType> testTypesExpanded;
    for (int i : testIds)
    {
        testTypesExpanded.push_back(testTypes[i]);
    }

    const ParticleType              refOw = waters.water().at("Ow");
    const ParticleType              refHw = waters.water().at("H");
    const std::vector<ParticleType> ref   = { refOw, refHw, refHw, refOw, refHw, refHw };

    EXPECT_TRUE(ref == testTypesExpanded);
}

TEST(NBlibTest, TopologyThrowsIdenticalParticleType)
{
    //! User error: Two different ParticleTypes with the same name
    ParticleType U235(ParticleName("Uranium"), Mass(235), C6(6.), C12(12.));
    ParticleType U238(ParticleName("Uranium"), Mass(238), C6(6.), C12(12.));

    Molecule ud235("UraniumDimer235");
    ud235.addParticle(ParticleName("U1"), U235);
    ud235.addParticle(ParticleName("U2"), U235);

    Molecule ud238("UraniumDimer238");
    ud238.addParticle(ParticleName("U1"), U238);
    ud238.addParticle(ParticleName("U2"), U238);

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(ud235, 1);
    EXPECT_THROW(topologyBuilder.addMolecule(ud238, 1), gmx::InvalidInputError);
}

TEST(NBlibTest, TopologyHasExclusions)
{
    WaterTopology         waters;
    Topology              watersTopology = waters.buildTopology(2);
    gmx::ListOfLists<int> testExclusions = watersTopology.getGmxExclusions();

    const std::vector<std::vector<int>> refExclusions = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
                                                          { 3, 4, 5 }, { 3, 4, 5 }, { 3, 4, 5 } };

    compareLists(testExclusions, refExclusions);
}

TEST(NBlibTest, TopologyHasNonbondedParameters)
{
    WaterTopology waters;
    Topology      watersTopology = waters.buildTopology(2);

    const Mass              refOwC6  = waters.water().at("Ow").c6();
    const Mass              refOwC12 = waters.water().at("Ow").c12();
    const Mass              refHwC6  = waters.water().at("H").c6();
    const Mass              refHwC12 = waters.water().at("H").c12();
    const std::vector<real> refC6    = { refOwC6, refHwC6, refHwC6, refOwC6, refHwC6, refHwC6 };
    const std::vector<real> refC12 = { refOwC12, refHwC12, refHwC12, refOwC12, refHwC12, refHwC12 };
    const std::vector<real> testC6 = expandQuantity(watersTopology, &ParticleType::c6);
    const std::vector<real> testC12 = expandQuantity(watersTopology, &ParticleType::c12);

    EXPECT_EQ(refC6, testC6);
    EXPECT_EQ(refC12, testC12);
}

//! Todo: this belongs to ForceCalculator
// TEST(NBlibTest, TopologyHasParticlesInfoAllVdw)
//{
//    TwoWaterMolecules      waters;
//    Topology               watersTopology = waters.buildTopology();
//    const std::vector<int> test           = watersTopology.getParticleInfoAllVdw();
//    std::vector<int>       ref;
//    ref.resize(watersTopology.numParticles());
//    for (size_t particleI = 0; particleI < ref.size(); particleI++)
//    {
//        SET_CGINFO_HAS_VDW(ref[particleI]);
//    }
//    EXPECT_EQ(ref, test);
//}

TEST(NBlibTest, toGmxExclusionBlockWorks)
{
    std::vector<std::tuple<int, int>> testInput{ { 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 },
                                                 { 1, 2 }, { 2, 0 }, { 2, 1 }, { 2, 2 } };

    std::vector<gmx::ExclusionBlock> reference;

    gmx::ExclusionBlock localBlock;
    localBlock.particleNumber.push_back(0);
    localBlock.particleNumber.push_back(1);
    localBlock.particleNumber.push_back(2);

    reference.push_back(localBlock);
    reference.push_back(localBlock);
    reference.push_back(localBlock);

    std::vector<gmx::ExclusionBlock> probe = detail::toGmxExclusionBlock(testInput);

    ASSERT_EQ(reference.size(), probe.size());
    for (size_t i = 0; i < reference.size(); ++i)
    {
        ASSERT_EQ(reference[i].particleNumber.size(), probe[i].particleNumber.size());
        for (size_t j = 0; j < reference[i].particleNumber.size(); ++j)
        {
            EXPECT_EQ(reference[i].particleNumber[j], probe[i].particleNumber[j]);
        }
    }
}

} // namespace
} // namespace test
} // namespace nblib
