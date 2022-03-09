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
 * This implements topology setup tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/topology.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/utility/listoflists.h"

#include "nblib/exception.h"
#include "nblib/particlesequencer.h"
#include "nblib/particletype.h"
#include "nblib/sequencing.hpp"
#include "nblib/tests/testsystems.h"
#include "nblib/topologyhelpers.h"

namespace nblib
{
namespace test
{
namespace
{

using ::testing::Eq;
using ::testing::Pointwise;

//! \brief Compares all element between two lists of lists
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

TEST(NBlibTest, TopologyHasNumParticles)
{
    WaterTopologyBuilder waters;
    Topology             watersTopology = waters.buildTopology(2);
    const int            test           = watersTopology.numParticles();
    const int            ref            = 6;
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasCharges)
{
    WaterTopologyBuilder     waters;
    Topology                 watersTopology = waters.buildTopology(2);
    const std::vector<real>& test           = watersTopology.getCharges();
    const std::vector<real>& ref = { Charges.at("Ow"), Charges.at("Hw"), Charges.at("Hw"),
                                     Charges.at("Ow"), Charges.at("Hw"), Charges.at("Hw") };
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasMasses)
{
    WaterTopologyBuilder waters;
    Topology             watersTopology = waters.buildTopology(2);

    const Mass              refOwMass = waters.water().at("Ow").mass();
    const Mass              refHwMass = waters.water().at("H").mass();
    const std::vector<Mass> ref = { refOwMass, refHwMass, refHwMass, refOwMass, refHwMass, refHwMass };
    const std::vector<Mass> test = expandQuantity(watersTopology, &ParticleType::mass);
    EXPECT_EQ(ref, test);
}

TEST(NBlibTest, TopologyHasParticleTypes)
{
    WaterTopologyBuilder             waters;
    Topology                         watersTopology = waters.buildTopology(2);
    const std::vector<ParticleType>& test           = watersTopology.getParticleTypes();
    const ParticleType               refOw          = waters.water().at("Ow");
    const ParticleType               refHw          = waters.water().at("H");
    const std::vector<ParticleType>& ref            = { refOw, refHw };
    const std::vector<ParticleType>& ref2           = { refHw, refOw };
    EXPECT_TRUE(ref == test || ref2 == test);
}

TEST(NBlibTest, TopologyHasParticleTypeIds)
{
    WaterTopologyBuilder waters;
    Topology             watersTopology = waters.buildTopology(2);

    const std::vector<int>&          testIds   = watersTopology.getParticleTypeIdOfAllParticles();
    const std::vector<ParticleType>& testTypes = watersTopology.getParticleTypes();

    std::vector<ParticleType> testTypesExpanded;
    testTypesExpanded.reserve(testTypes.size());
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
    ParticleType U235(ParticleTypeName("Uranium"), Mass(235));
    ParticleType U238(ParticleTypeName("Uranium"), Mass(238));

    Molecule ud235(MoleculeName("UraniumDimer235"));
    ud235.addParticle(ParticleName("U1"), U235);
    ud235.addParticle(ParticleName("U2"), U235);

    Molecule ud238(MoleculeName("UraniumDimer238"));
    ud238.addParticle(ParticleName("U1"), U238);
    ud238.addParticle(ParticleName("U2"), U238);

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(ud235, 1);
    EXPECT_THROW(topologyBuilder.addMolecule(ud238, 1), InputException);
}

TEST(NBlibTest, TopologyHasExclusions)
{
    WaterTopologyBuilder        waters;
    Topology                    watersTopology = waters.buildTopology(2);
    ExclusionLists<int>         exclusionLists = watersTopology.exclusionLists();
    const gmx::ListOfLists<int> testExclusions(std::move(exclusionLists.ListRanges),
                                               std::move(exclusionLists.ListElements));

    const std::vector<std::vector<int>>& refExclusions = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 },
                                                           { 3, 4, 5 }, { 3, 4, 5 }, { 3, 4, 5 } };

    compareLists(testExclusions, refExclusions);
}

TEST(NBlibTest, TopologyHasSequencing)
{
    WaterTopologyBuilder waters;
    Topology             watersTopology = waters.buildTopology(2);

    EXPECT_EQ(0, watersTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("Oxygen")));
    EXPECT_EQ(1, watersTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H1")));
    EXPECT_EQ(2, watersTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H2")));
    EXPECT_EQ(3, watersTopology.sequenceID(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("Oxygen")));
    EXPECT_EQ(4, watersTopology.sequenceID(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("H1")));
    EXPECT_EQ(5, watersTopology.sequenceID(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("H2")));
}

TEST(NBlibTest, TopologyCanAggregateBonds)
{
    Molecule water    = WaterMoleculeBuilder{}.waterMolecule();
    Molecule methanol = MethanolMoleculeBuilder{}.methanolMolecule();

    std::vector<std::tuple<Molecule, int>> molecules{ std::make_tuple(water, 2),
                                                      std::make_tuple(methanol, 1) };
    std::vector<HarmonicBondType>          bonds;
    std::vector<size_t>                    bondsExpansion;
    std::tie(bondsExpansion, bonds) = detail::collectInteractions<HarmonicBondType>(molecules);

    std::vector<HarmonicBondType> bondsTest;
    // use the expansionArray (bondsExpansion) to expand to the full list if bonds
    std::transform(begin(bondsExpansion),
                   end(bondsExpansion),
                   std::back_inserter(bondsTest),
                   [&bonds](size_t i) { return bonds[i]; });

    std::vector<HarmonicBondType> waterBonds =
            pickType<HarmonicBondType>(water.interactionData()).interactionTypes_;
    std::vector<HarmonicBondType> methanolBonds =
            pickType<HarmonicBondType>(methanol.interactionData()).interactionTypes_;

    std::vector<HarmonicBondType> bondsReference;
    std::copy(begin(waterBonds), end(waterBonds), std::back_inserter(bondsReference));
    std::copy(begin(waterBonds), end(waterBonds), std::back_inserter(bondsReference));
    std::copy(begin(methanolBonds), end(methanolBonds), std::back_inserter(bondsReference));

    EXPECT_EQ(bondsTest, bondsReference);
}

TEST(NBlibTest, TopologyCanSequencePairIDs)
{
    Molecule water    = WaterMoleculeBuilder{}.waterMolecule();
    Molecule methanol = MethanolMoleculeBuilder{}.methanolMolecule();

    std::vector<std::tuple<Molecule, int>> molecules{ std::make_tuple(water, 2),
                                                      std::make_tuple(methanol, 1) };
    ParticleSequencer                      particleSequencer;
    particleSequencer.build(molecules);
    auto pairs = detail::sequenceIDs<HarmonicBondType>(molecules, particleSequencer);

    int Ow1 = particleSequencer(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("Oxygen"));
    int H11 = particleSequencer(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H1"));
    int H12 = particleSequencer(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H2"));
    int Ow2 = particleSequencer(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("Oxygen"));
    int H21 = particleSequencer(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("H1"));
    int H22 = particleSequencer(MoleculeName("SOL"), 1, ResidueName("SOL"), ParticleName("H2"));

    int Me  = particleSequencer(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("Me1"));
    int MeO = particleSequencer(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("O2"));
    int MeH = particleSequencer(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("H3"));

    /// \cond DO_NOT_DOCUMENT
#define SORT(i, j) (i < j) ? i : j, (i < j) ? j : i

    std::vector<CoordinateIndex<HarmonicBondType>> pairs_reference{
        { SORT(Ow1, H11) }, { SORT(Ow1, H12) }, { SORT(Ow2, H21) },
        { SORT(Ow2, H22) }, { SORT(MeO, MeH) }, { SORT(MeO, Me) }
    };
#undef SORT
    /// \endcond

    EXPECT_EQ(pairs, pairs_reference);
}

TEST(NBlibTest, TopologySequenceIdThrows)
{
    Molecule water    = WaterMoleculeBuilder{}.waterMolecule();
    Molecule methanol = MethanolMoleculeBuilder{}.methanolMolecule();

    std::vector<std::tuple<Molecule, int>> molecules{ std::make_tuple(water, 2),
                                                      std::make_tuple(methanol, 1) };
    ParticleSequencer                      particleSequencer;
    particleSequencer.build(molecules);
    auto pairs = detail::sequenceIDs<HarmonicBondType>(molecules, particleSequencer);

    // Input error: no particle called O-Atom in molecule "water"
    EXPECT_THROW(particleSequencer(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("O-Atom")),
                 InputException);
}

TEST(NBlibTest, TopologyCanEliminateDuplicateBonds)
{
    HarmonicBondType b1(1.0, 2.0);
    HarmonicBondType b2(1.1, 2.1);
    HarmonicBondType b3(1.2, 2.2);

    // can be compressed to {b1,b2,b3} + {1,1,2,0,1,0,2,2}
    std::vector<HarmonicBondType> bonds{ b2, b2, b3, b1, b2, b1, b3, b3 };

    // expected output
    std::vector<HarmonicBondType> uniqueBondsReference{ b1, b2, b3 };
    std::vector<size_t>           indicesReference{ 1, 1, 2, 0, 1, 0, 2, 2 };

    std::tuple<std::vector<size_t>, std::vector<HarmonicBondType>> bondData =
            detail::eliminateDuplicateInteractions(bonds);

    auto indices     = std::get<0>(bondData);
    auto uniqueBonds = std::get<1>(bondData);

    EXPECT_EQ(uniqueBondsReference, uniqueBonds);
    EXPECT_EQ(indicesReference, indices);
}

TEST(NBlibTest, TopologyListedInteractions)
{
    // Todo: add angles to SpcMethanolTopologyBuilder and extend interactions_reference below to
    // Todo: include angles

    Topology spcTopology = SpcMethanolTopologyBuilder{}.buildTopology(1, 2);

    auto  interactionData = spcTopology.getInteractionData();
    auto& harmonicBonds   = pickType<HarmonicBondType>(interactionData);

    auto& indices = harmonicBonds.indices;
    auto& bonds   = harmonicBonds.parameters;

    std::map<std::tuple<int, int>, HarmonicBondType> interactions_test;
    for (auto& ituple : indices)
    {
        interactions_test[std::make_tuple(std::get<0>(ituple), std::get<1>(ituple))] =
                bonds[std::get<2>(ituple)];
    }

    // there should be 3 unique HarmonicBondType instances
    EXPECT_EQ(bonds.size(), 3);
    // and 6 interaction pairs (bonds)
    EXPECT_EQ(indices.size(), 6);

    HarmonicBondType ohBond(1., 1.);
    HarmonicBondType ohBondMethanol(1.01, 1.02);
    HarmonicBondType ometBond(1.1, 1.2);

    std::map<std::tuple<int, int>, HarmonicBondType> interactions_reference;

    int Ow = spcTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("Oxygen"));
    int H1 = spcTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H1"));
    int H2 = spcTopology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H2"));

    int Me1 = spcTopology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("Me1"));
    int MeO1 = spcTopology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("O2"));
    int MeH1 = spcTopology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("H3"));

    int Me2 = spcTopology.sequenceID(MoleculeName("MeOH"), 1, ResidueName("MeOH"), ParticleName("Me1"));
    int MeO2 = spcTopology.sequenceID(MoleculeName("MeOH"), 1, ResidueName("MeOH"), ParticleName("O2"));
    int MeH2 = spcTopology.sequenceID(MoleculeName("MeOH"), 1, ResidueName("MeOH"), ParticleName("H3"));

    /// \cond DO_NOT_DOCUMENT
#define SORT(i, j) (i < j) ? i : j, (i < j) ? j : i
    interactions_reference[std::make_tuple(SORT(Ow, H1))]     = ohBond;
    interactions_reference[std::make_tuple(SORT(Ow, H2))]     = ohBond;
    interactions_reference[std::make_tuple(SORT(MeO1, MeH1))] = ohBondMethanol;
    interactions_reference[std::make_tuple(SORT(MeO1, Me1))]  = ometBond;
    interactions_reference[std::make_tuple(SORT(MeO2, MeH2))] = ohBondMethanol;
    interactions_reference[std::make_tuple(SORT(MeO2, Me2))]  = ometBond;
#undef SORT
    /// \endcond

    EXPECT_TRUE(std::equal(
            begin(interactions_reference), end(interactions_reference), begin(interactions_test)));
}

TEST(NBlibTest, TopologyListedInteractionsMultipleTypes)
{
    Molecule water    = WaterMoleculeBuilder{}.waterMolecule();
    Molecule methanol = MethanolMoleculeBuilder{}.methanolMolecule();

    CubicBondType testBond(1., 1., 1.);
    HarmonicAngle testAngle(1, Degrees(1));

    water.addInteraction(ParticleName("H1"), ParticleName("H2"), testBond);
    water.addInteraction(ParticleName("H1"), ParticleName("Oxygen"), ParticleName("H2"), testAngle);

    ParticleTypesInteractions nbInteractions;
    std::vector<std::string>  particleTypeNames = { "Ow", "H", "OMet", "CMet" };
    for (const auto& name : particleTypeNames)
    {
        nbInteractions.add(ParticleTypeName(name), C6(0), C12(0));
    }

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(water, 1);
    topologyBuilder.addMolecule(methanol, 1);
    topologyBuilder.addParticleTypesInteractions(nbInteractions);

    Topology topology = topologyBuilder.buildTopology();

    auto  interactionData   = topology.getInteractionData();
    auto& harmonicBonds     = pickType<HarmonicBondType>(interactionData);
    auto& cubicBonds        = pickType<CubicBondType>(interactionData);
    auto& angleInteractions = pickType<HarmonicAngle>(interactionData);

    HarmonicBondType              ohBond(1., 1.);
    HarmonicBondType              ohBondMethanol(1.01, 1.02);
    HarmonicBondType              ometBond(1.1, 1.2);
    std::vector<HarmonicBondType> harmonicBondsReference{ ohBond, ohBondMethanol, ometBond };

    EXPECT_EQ(harmonicBonds.parameters, harmonicBondsReference);

    int H1 = topology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H1"));
    int H2 = topology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("H2"));
    int Ow = topology.sequenceID(MoleculeName("SOL"), 0, ResidueName("SOL"), ParticleName("Oxygen"));

    int Me1 = topology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("Me1"));
    int MeO1 = topology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("O2"));
    int MeH1 = topology.sequenceID(MoleculeName("MeOH"), 0, ResidueName("MeOH"), ParticleName("H3"));

    std::vector<CubicBondType>                   cubicBondsReference{ testBond };
    std::vector<InteractionIndex<CubicBondType>> cubicIndicesReference{
        { std::min(H1, H2), std::max(H1, H2), 0 }
    };
    EXPECT_EQ(cubicBondsReference, cubicBonds.parameters);
    EXPECT_EQ(cubicIndicesReference, cubicBonds.indices);

    HarmonicAngle                                methanolAngle(397.5, Degrees(108.52));
    std::vector<HarmonicAngle>                   angleReference{ testAngle, methanolAngle };
    std::vector<InteractionIndex<HarmonicAngle>> angleIndicesReference{
        { std::min(H1, H2), Ow, std::max(H1, H2), 0 }, { std::min(MeH1, MeO1), Me1, std::max(MeO1, MeH1), 1 }
    };
    EXPECT_EQ(angleReference, angleInteractions.parameters);
    EXPECT_EQ(angleIndicesReference, angleInteractions.indices);
}

TEST(NBlibTest, TopologyInvalidParticleInInteractionThrows)
{
    Molecule water    = WaterMoleculeBuilder{}.waterMolecule();
    Molecule methanol = MethanolMoleculeBuilder{}.methanolMolecule();

    HarmonicBondType testBond(1., 1.);

    // Invalid input: no particle named "Iron" in molecule water
    water.addInteraction(ParticleName("H1"), ParticleName("Iron"), testBond);

    ParticleTypesInteractions nbInteractions;
    std::vector<std::string>  particleTypeNames = { "Ow", "H", "OMet", "CMet" };
    for (const auto& name : particleTypeNames)
    {
        nbInteractions.add(ParticleTypeName(name), C6(0), C12(0));
    }

    TopologyBuilder topologyBuilder;
    topologyBuilder.addMolecule(water, 1);
    topologyBuilder.addMolecule(methanol, 1);
    topologyBuilder.addParticleTypesInteractions(nbInteractions);

    EXPECT_THROW(topologyBuilder.buildTopology(), InputException);
}


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

    std::vector<gmx::ExclusionBlock> probe = toGmxExclusionBlock(testInput);

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
