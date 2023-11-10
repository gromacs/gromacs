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
 * This implements basic nblib utility tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/listed_forces/aggregate_transformations.hpp"

#include <numeric>

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "linear_chain_input.hpp"
#include "poly_dimethyl_butene_input.hpp"


namespace nblib
{
namespace test
{
namespace
{

ListedInteractionData someBonds()
{
    ListedInteractionData         interactions;
    HarmonicBondType              bond{ 10, 0.1 };
    std::vector<HarmonicBondType> bonds{ bond };
    pickType<HarmonicBondType>(interactions).parametersA = bonds;

    std::vector<InteractionIndex<HarmonicBondType>> bondIndices{
        { 0, 1, 0 }, { 1, 2, 0 }, { 2, 3, 0 }, { 3, 4, 0 }
    };
    pickType<HarmonicBondType>(interactions).indices = std::move(bondIndices);

    return interactions;
}

TEST(ListedAggregates, DeleteInteractions)
{
    ListedInteractionData interactions = someBonds();

    std::vector<int> deleteList{ 0, 2, 3 };

    deleteInteractions(pickType<HarmonicBondType>(interactions), deleteList);

    std::vector<InteractionIndex<HarmonicBondType>> reference{ { 1, 2, 0 } };
    EXPECT_EQ(pickType<HarmonicBondType>(interactions).indices, reference);
}

ListedInteractionData someBondsAndAngles()
{
    ListedInteractionData         interactions;
    HarmonicBondType              bond1{ 10, 0.1 };
    HarmonicBondType              bond2{ 20, 0.2 };
    std::vector<HarmonicBondType> bonds{ bond1, bond2 };
    pickType<HarmonicBondType>(interactions).parametersA = bonds;

    HarmonicAngle              angle1(100, Degrees(100));
    HarmonicAngle              angle2(200, Degrees(101));
    std::vector<HarmonicAngle> angles{ angle1, angle2 };
    pickType<HarmonicAngle>(interactions).parametersA = angles;

    std::vector<InteractionIndex<HarmonicBondType>> bondIndices{ { 0, 1, 0 }, { 1, 2, 1 }, { 2, 3, 1 } };
    pickType<HarmonicBondType>(interactions).indices = std::move(bondIndices);

    std::vector<InteractionIndex<HarmonicAngle>> angleIndices{ { 0, 1, 2, 0 }, { 1, 2, 3, 1 } };
    pickType<HarmonicAngle>(interactions).indices = std::move(angleIndices);

    return interactions;
}

TEST(ListedAggregates, MigrateTwoToThree)
{
    using AggregateType = ThreeCenterAggregate<HarmonicBondType, HarmonicAngle>;

    ListedInteractionData interactions = someBondsAndAngles();
    auto                  bonds        = pickType<HarmonicBondType>(interactions);
    auto                  angles       = pickType<HarmonicAngle>(interactions);

    ListedTypeData<AggregateType> aggregateData;
    migrateTwoToThree(pickType<HarmonicBondType>(interactions),
                      pickType<HarmonicAngle>(interactions),
                      aggregateData);

    EXPECT_EQ(aggregateData.indices.size(), 2);
    EXPECT_EQ(pickType<HarmonicBondType>(interactions).indices.size(), 1);

    EXPECT_EQ(aggregateData.parametersA[0].twoCenter(), bonds.parametersA[0]);
    EXPECT_EQ(aggregateData.parametersA[1].twoCenter(), bonds.parametersA[1]);
    EXPECT_EQ(aggregateData.parametersA[0].carrier(), angles.parametersA[0]);
    EXPECT_EQ(aggregateData.parametersA[1].carrier(), angles.parametersA[1]);

    EXPECT_TRUE(aggregateData.parametersA[0].manifest & AggregateType::bond_ij);
    EXPECT_TRUE(aggregateData.parametersA[1].manifest & AggregateType::bond_ij);

    std::vector<InteractionIndex<AggregateType>> aggregateIndices{ { 0, 1, 2, 0 }, { 1, 2, 3, 1 } };
    EXPECT_EQ(aggregateData.indices, aggregateIndices);
}

TEST(ListedAggregates, MigrateLinearChain)
{
    using AggregateType = ThreeCenterAggregate<HarmonicBondType, HarmonicAngle>;

    int                           length = 10;
    LinearChainData               linearChain(length, 0.0);
    ListedInteractionData         interactions = linearChain.interactions;
    auto&                         bonds        = pickType<HarmonicBondType>(interactions);
    auto&                         angles       = pickType<HarmonicAngle>(interactions);
    ListedTypeData<AggregateType> aggregates;

    EXPECT_EQ(length - 1, bonds.indices.size());

    migrateTwoToThree(bonds, angles, aggregates);

    // all except the last bond get migrated
    EXPECT_EQ(1, bonds.indices.size());
}

TEST(ListedAggregates, MigratePolymerBonds)
{
    using AggregateType = ThreeCenterAggregate<HarmonicBondType, HarmonicAngle>;

    int                           length = 1;
    PolyDimethylButene            polymer(length);
    ListedInteractionData         interactions = polymer.interactions;
    auto&                         bonds        = pickType<HarmonicBondType>(interactions);
    auto&                         angles       = pickType<HarmonicAngle>(interactions);
    ListedTypeData<AggregateType> aggregates;

    EXPECT_EQ(6 * length, bonds.indices.size());

    migrateTwoToThree(bonds, angles, aggregates);

    std::cout << bonds.indices[0][0] << " " << bonds.indices[0][1] << std::endl;
    // all except the last bond get migrated
    EXPECT_EQ(0, bonds.indices.size());
}

TEST(ListedAggregates, MigratePolymerAngles)
{
    using AggregateType = TypeListElement_t<0, AggregateTypes>;

    int                   length = 2;
    PolyDimethylButene    polymer(length);
    ListedInteractionData interactions = polymer.interactions;
    auto& bonds      = pickType<typename AggregateType::TwoCenterAggregateType>(interactions);
    auto& angles     = pickType<typename AggregateType::ThreeCenterAggregateType>(interactions);
    auto& dihedrals  = pickType<typename AggregateType::CarrierType>(interactions);
    auto& pairs      = pickType<typename AggregateType::PairAggregateType>(interactions);
    auto& aggregates = pickType<AggregateType>(interactions);

    EXPECT_EQ(12 * length, angles.indices.size());
    EXPECT_EQ(18 * length, dihedrals.indices.size());

    migrateParameters(dihedrals, aggregates);
    migrateTwoToFour(bonds, aggregates);
    migrateThreeToFour(angles, aggregates);
    migratePairs(pairs, aggregates);

    EXPECT_EQ(18 * length, aggregates.indices.size());

    // all except the last angle get migrated
    EXPECT_EQ(2 + 2 * length, angles.indices.size());
}

TEST(ListedAggregates, HashTest)
{
    std::unordered_map<util::array<int, 4>, int, ArrayHasher<4>> map;

    for (int i = 0; i < 10; ++i)
    {
        IndexArray<4> idx{ i, i + 1, i + 2, i + 3 };
        map[idx] = i;
    }

    IndexArray<4> test{ 2, 3, 4, 5 };

    EXPECT_EQ(map[test], 2);
}

TEST(ListedAggregates, BondConfig)
{
    int                 config = 0;
    util::array<int, 2> grouping{ 1, 3 };

    config = setBondIndices(grouping, config);

    auto probe = getBondIndices(config);

    EXPECT_EQ(probe[0], 1);
    EXPECT_EQ(probe[1], 3);
}

TEST(ListedAggregates, PairConfig)
{
    int                 config = 0;
    util::array<int, 2> grouping{ 1, 3 };

    config = setPairIndices(grouping, config);

    auto probe = getPairIndices(config);

    EXPECT_EQ(probe[0], 1);
    EXPECT_EQ(probe[1], 3);
}

TEST(ListedAggregates, HashIntegrate)
{
    int                   length = 2;
    PolyDimethylButene    polymer(length);
    ListedInteractionData interactions = polymer.interactions;

    auto& bonds     = pickType<HarmonicBondType>(interactions);
    auto& angles    = pickType<HarmonicAngle>(interactions);
    auto& dihedrals = pickType<ProperDihedral>(interactions);
    auto& pairs     = pickType<PairLJType>(interactions);

    pairs.indices.push_back({ 4, 0, 40 });
    pairs.indices.push_back({ 3, 5, 41 });
    pairs.indices.push_back({ 4, 12, 42 });

    auto bondIndexCopy  = bonds.indices;
    auto angleIndexCopy = angles.indices;
    auto pairIndexCopy  = pairs.indices;

    // for (auto d : dihedrals.indices)
    //{
    //    printf("%d %d %d %d\n", d[0], d[1], d[2], d[3]);
    //}

    std::vector<IndexArray<4>> aggregates(dihedrals.indices.size());

    integrateBonds(bonds.indices, dihedrals.indices, aggregates);
    integrateAngles(angles.indices, dihedrals.indices, aggregates);
    integratePairs(pairs.indices, dihedrals.indices, aggregates);

    for (size_t i = 0; i < aggregates.size(); ++i)
    {
        auto host = dihedrals.indices[i];

        if (hasBond(aggregates[i][0]))
        {
            auto config = getBondIndices(aggregates[i][0]);
            int  bond_i = host[config[0]];
            int  bond_j = host[config[1]];
            auto it = std::find_if(begin(bondIndexCopy), end(bondIndexCopy), [bond_i, bond_j](auto i) {
                return i[0] == bond_i && i[1] == bond_j;
            });
            ASSERT_NE(it, bondIndexCopy.end());
            // check that the referenced bond parameter index matches the original
            EXPECT_EQ((*it)[2], aggregates[i][1]);
        }
        if (hasAngle(aggregates[i][0]))
        {
            auto config  = getAngleIndices(aggregates[i][0]);
            int  angle_i = host[config[0]];
            int  angle_j = host[config[1]];
            int  angle_k = host[config[2]];
            auto it      = std::find_if(
                    begin(angleIndexCopy), end(angleIndexCopy), [angle_i, angle_j, angle_k](auto i) {
                        return i[0] == angle_i && i[1] == angle_j && i[2] == angle_k;
                    });
            ASSERT_NE(it, angleIndexCopy.end());
            // check that the referenced bond parameter index matches the original
            EXPECT_EQ((*it)[3], aggregates[i][2]);
        }
        if (hasPair(aggregates[i][0]))
        {
            auto config = getPairIndices(aggregates[i][0]);
            int  bond_i = host[config[0]];
            int  bond_j = host[config[1]];
            auto it = std::find_if(begin(pairIndexCopy), end(pairIndexCopy), [bond_i, bond_j](auto i) {
                return i[0] == bond_i && i[1] == bond_j;
            });
            ASSERT_NE(it, pairIndexCopy.end());
            // check that the referenced bond parameter index matches the original
            EXPECT_EQ((*it)[2], aggregates[i][3]);
        }
    }

    int numIntegratedBonds =
            std::count_if(begin(aggregates), end(aggregates), [](auto i) { return hasBond(i[0]); });
    int numIntegratedAngles =
            std::count_if(begin(aggregates), end(aggregates), [](auto i) { return hasAngle(i[0]); });
    int numIntegratedPairs =
            std::count_if(begin(aggregates), end(aggregates), [](auto i) { return hasPair(i[0]); });

    EXPECT_EQ(numIntegratedPairs, 3);

    EXPECT_EQ(bonds.indices.size() + numIntegratedBonds, bondIndexCopy.size());
    EXPECT_EQ(angles.indices.size() + numIntegratedAngles, angleIndexCopy.size());
    EXPECT_EQ(pairs.indices.size() + numIntegratedPairs, pairIndexCopy.size());
}

} // namespace
} // namespace test
} // namespace nblib
