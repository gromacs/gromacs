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
#include <array>

#include <gtest/gtest.h>

#include "listed_forces/conversionscommon.h"

#include "gromacs/topology/forcefieldparameters.h"

#include "testutils/testasserts.h"

#include "nblib/box.h"

#include "listedtesthelpers.h"
#include "testhelpers.h"

namespace nblib
{
namespace test
{
namespace
{

ListedInteractionData someBondsAndAngles()
{
    ListedInteractionData         interactions;
    HarmonicBondType              bond1{ 10, 0.1 };
    HarmonicBondType              bond2{ 20, 0.2 };
    std::vector<HarmonicBondType> bonds{ bond1, bond2 };
    pickType<HarmonicBondType>(interactions).parameters = bonds;

    HarmonicAngle              angle1(100, Degrees(100));
    HarmonicAngle              angle2(200, Degrees(101));
    std::vector<HarmonicAngle> angles{ angle1, angle2 };
    pickType<HarmonicAngle>(interactions).parameters = angles;

    std::vector<InteractionIndex<HarmonicBondType>> bondIndices{ { 0, 1, 0 }, { 1, 2, 0 }, { 2, 3, 1 } };
    pickType<HarmonicBondType>(interactions).indices = std::move(bondIndices);

    std::vector<InteractionIndex<HarmonicAngle>> angleIndices{ { 0, 1, 2, 0 }, { 1, 2, 3, 1 } };
    pickType<HarmonicAngle>(interactions).indices = std::move(angleIndices);

    return interactions;
}

TEST(ListedShims, ParameterConversion)
{
    ListedInteractionData interactions = someBondsAndAngles();

    auto [idef, gmx_params] = convertToGmxInteractions(interactions);

    EXPECT_EQ(gmx_params->iparams.size(), 4);
    EXPECT_EQ(gmx_params->iparams[0].harmonic.rA,
              pickType<HarmonicBondType>(interactions).parameters[0].equilConstant());
    EXPECT_REAL_EQ_TOL(gmx_params->iparams[2].harmonic.rA,
                       pickType<HarmonicAngle>(interactions).parameters[0].equilConstant() / DEG2RAD,
                       gmx::test::defaultRealTolerance());

    EXPECT_EQ(idef->il[F_BONDS].iatoms.size(), 9);
    std::vector<int> bondIatoms{ 0, 0, 1, 0, 1, 2, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_BONDS].iatoms, bondIatoms);
    std::vector<int> angleIatoms{ 2, 0, 1, 2, 3, 1, 2, 3 };
    EXPECT_EQ(idef->il[F_ANGLES].iatoms, angleIatoms);
    idef->clear();
}

// Dummy no-op function to ignore unused variables
template<class T>
void ignore(const T&)
{
}

TEST(ListedShims, GmxToNblibConversion)
{
    ListedInteractionData interactions = someBondsAndAngles();

    auto [idef, gmx_params] = convertToGmxInteractions(interactions);
    ignore(gmx_params);

    ListedInteractionData convertedData = convertToNblibInteractions(*idef);

    // compare some bond parameter data
    size_t numHarmonicBondsOrig = pickType<HarmonicBondType>(interactions).parameters.size();
    size_t numHarmonicBondsConv = pickType<HarmonicBondType>(convertedData).parameters.size();
    EXPECT_EQ(numHarmonicBondsOrig, numHarmonicBondsConv);

    size_t numHarmonicAnglesOrig = pickType<HarmonicAngle>(interactions).parameters.size();
    size_t numHarmonicAnglesConv = pickType<HarmonicAngle>(convertedData).parameters.size();
    EXPECT_EQ(numHarmonicAnglesOrig, numHarmonicAnglesConv);

    HarmonicAngle harmonicAngleOrig0 = pickType<HarmonicAngle>(interactions).parameters[0];
    HarmonicAngle harmonicAngleConv0 = pickType<HarmonicAngle>(convertedData).parameters[0];
    EXPECT_TRUE(harmonicAngleOrig0 == harmonicAngleConv0);

    // compare some index data
    auto haIndicesOrig = pickType<HarmonicAngle>(interactions).indices[0];
    auto haIndicesConv = pickType<HarmonicAngle>(convertedData).indices[0];
    EXPECT_EQ(haIndicesOrig, haIndicesConv);
}

std::vector<gmx::RVec> twoCenterCoordinates = { { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 } };

std::vector<gmx::RVec> threeCenterCoordinates = { { 1.382, 1.573, 1.482 },
                                                  { 1.281, 1.559, 1.596 },
                                                  { 1.292, 1.422, 1.663 } };

std::vector<gmx::RVec> fourCentercoordinates = { { 1.382, 1.573, 1.482 },
                                                 { 1.281, 1.559, 1.596 },
                                                 { 1.292, 1.422, 1.663 },
                                                 { 1.189, 1.407, 1.775 } };

std::array<std::vector<gmx::RVec>, 3> NCenterCoordinates{ twoCenterCoordinates,
                                                          threeCenterCoordinates,
                                                          fourCentercoordinates };

template<class Interaction>
struct TypeInput
{
    typedef Interaction         type; // needed for pickType
    ListedTypeData<Interaction> interactionData;
    static constexpr int        MinimumCenterCount = 2;
    // assign coordinates depending on the number of centers in the interaction type from the array above
    std::vector<gmx::RVec> coordinates =
            std::get<NCenter<Interaction>{} - MinimumCenterCount>(NCenterCoordinates);
};

std::tuple TestInput{
    // Two Center Types
    TypeInput<HarmonicBondType>{
            { { { HarmonicBondType(500.0, 0.15) } }, { indexVector<HarmonicBondType>() } } },
    TypeInput<G96BondType>{ { { { G96BondType(50.0, 0.15) } }, { indexVector<G96BondType>() } } },
    TypeInput<CubicBondType>{ { { { CubicBondType(50.0, 2.0, 0.16) } }, { indexVector<CubicBondType>() } } },
    TypeInput<MorseBondType>{ { { { MorseBondType(30.0, 2.7, 0.15) } }, { indexVector<MorseBondType>() } } },
    TypeInput<FENEBondType>{ { { { FENEBondType(5.0, 0.4) } }, { indexVector<FENEBondType>() } } },
    // Three Center Types
    TypeInput<HarmonicAngle>{
            { { { HarmonicAngle(2.2, Degrees(91.0)) } }, { indexVector<HarmonicAngle>() } } },
    TypeInput<G96Angle>{ { { { G96Angle(50.0, Degrees(100)) } }, { indexVector<G96Angle>() } } },
    TypeInput<RestrictedAngle>{
            { { { RestrictedAngle(50.0, Degrees(100)) } }, { indexVector<RestrictedAngle>() } } },
    TypeInput<LinearAngle>{ { { { LinearAngle(50.0, 0.4) } }, { indexVector<LinearAngle>() } } },
    TypeInput<QuarticAngle>{ { { { QuarticAngle(1.1, 2.3, 4.6, 7.8, 9.2, Degrees(87)) } },
                               { indexVector<QuarticAngle>() } } },
    TypeInput<CrossBondBond>{ { { { CrossBondBond(45.0, 0.8, 0.7) } }, { indexVector<CrossBondBond>() } } },
    TypeInput<CrossBondAngle>{
            { { { CrossBondAngle(45.0, 0.8, 0.7, 0.3) } }, { indexVector<CrossBondAngle>() } } },
    // Four Center Types
    TypeInput<ProperDihedral>{
            { { { ProperDihedral(Degrees(45), 2.3, 1) } }, { indexVector<ProperDihedral>() } } }
};

template<class... Ts>
ListedInteractionData combineTestInput(std::tuple<Ts...> testInput)
{
    ListedInteractionData interactionData;
    // transfer all elements of testInput into the returned ListedInteractionData
    // use a lambda + for_each_tuple
    auto copyParamsOneType = [&interactionData](const auto& typeInput) {
        for (size_t i = 0; i < typeInput.interactionData.parameters.size(); i++)
        {
            auto interactionParams = typeInput.interactionData.parameters[i];
            using InteractionType  = decltype(interactionParams);
            pickType<InteractionType>(interactionData).parameters.push_back(interactionParams);

            auto indices = typeInput.interactionData.indices[i];
            pickType<InteractionType>(interactionData).indices.push_back(indices);
        }
    };
    for_each_tuple(copyParamsOneType, testInput);

    return interactionData;
}

TEST(NBlibTest, GmxToNblibConversionAllTypes)
{
    ListedInteractionData originalData = combineTestInput(TestInput);

    auto [idef, gmx_params] = convertToGmxInteractions(originalData);
    ignore(gmx_params);

    ListedInteractionData convertedData = convertToNblibInteractions(*idef);

    auto compareParamsAndIndices = [&convertedData](auto& original) {
        if (!original.parameters.empty())
        {
            using InteractionType = typename std::decay_t<decltype(original)>::type;
            const auto converted  = pickType<InteractionType>(convertedData);

            // compare parameters
            // comparing the two cosine angles with lower tolerance as this test introduces
            // numerical errors via the x != arccos(cos(x)) comparison
            if constexpr (std::is_same<InteractionType, G96Angle>::value
                          || std::is_same<InteractionType, RestrictedAngle>::value)
            {
                // There is only one interaction parameter checked per loop through the parameter
                // list, so it is correct to only check the first element in the parameters list
                // when we can't guaranteed equality as in the other conversions
                EXPECT_EQ(original.parameters[0].forceConstant(), converted.parameters[0].forceConstant());
                const real theta_orig = original.parameters[0].equilConstant();
                const real theta_conv = converted.parameters[0].equilConstant();
                EXPECT_FLOAT_DOUBLE_EQ_TOL(theta_conv,
                                           theta_orig,
                                           theta_orig,
                                           gmx::test::relativeToleranceAsFloatingPoint(theta_orig, 1e-6));
            }
            else
            {
                EXPECT_TRUE(original.parameters == converted.parameters);
            }

            // compare indices
            EXPECT_TRUE(original.indices == converted.indices);
        }
    };
    for_each_tuple(compareParamsAndIndices, originalData);
}

TEST(NBlibTest, EndToEndListedComparison)
{
    int                   numParticles    = 4;
    ListedInteractionData interactionData = combineTestInput(TestInput);

    Box                    box(1.0);
    std::vector<gmx::RVec> coordinates = {
        { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 }, { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 }
    };

    compareNblibAndGmxListedImplementations(interactionData, coordinates, numParticles, 1, box, 1e-2);
}

template<typename Interaction>
class NblibGmxListed : public testing::Test
{
public:
    void compareNblibAndGmx()
    {
        ListedInteractionData  interactionData;
        TypeInput<Interaction> typeInput       = pickType<Interaction>(TestInput);
        pickType<Interaction>(interactionData) = typeInput.interactionData;
        compareNblibAndGmxListedImplementations(
                interactionData, typeInput.coordinates, typeInput.coordinates.size(), 1, Box(1.0), 1e-4);
    }
};

// Extract a typelist interaction types from the TestInput tuple
template<class TestInputofInteractionType>
using ExtractType = typename TestInputofInteractionType::type;
using ListedTypes = Map<ExtractType, decltype(TestInput)>;
using TestTypes   = Reduce<::testing::Types, ListedTypes>;

TYPED_TEST_SUITE_P(NblibGmxListed);

TYPED_TEST_P(NblibGmxListed, SameForcesOnBoth)
{
    this->compareNblibAndGmx();
}

REGISTER_TYPED_TEST_SUITE_P(NblibGmxListed, SameForcesOnBoth);

INSTANTIATE_TYPED_TEST_SUITE_P(CompareEachTypeInNblibAndGmx, NblibGmxListed, TestTypes);

} // namespace
} // namespace test
} // namespace nblib
