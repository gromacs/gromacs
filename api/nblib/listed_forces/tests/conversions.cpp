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

#include "gromacs/topology/forcefieldparameters.h"

#include "testutils/testasserts.h"

#include "nblib/box.h"
#include "nblib/listed_forces/conversionscommon.h"
#include "nblib/tests/testhelpers.h"

#include "listedtesthelpers.h"
#include "listedtypeinput.hpp"

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
    pickType<HarmonicBondType>(interactions).parametersA = bonds;

    HarmonicAngle              angle1(100, Degrees(100));
    HarmonicAngle              angle2(200, Degrees(101));
    std::vector<HarmonicAngle> angles{ angle1, angle2 };
    pickType<HarmonicAngle>(interactions).parametersA = angles;

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
              pickType<HarmonicBondType>(interactions).parametersA[0].equilConstant());
    EXPECT_REAL_EQ_TOL(gmx_params->iparams[2].harmonic.rA,
                       pickType<HarmonicAngle>(interactions).parametersA[0].equilConstant() / DEG2RAD,
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
    size_t numHarmonicBondsOrig = pickType<HarmonicBondType>(interactions).parametersA.size();
    size_t numHarmonicBondsConv = pickType<HarmonicBondType>(convertedData).parametersA.size();
    EXPECT_EQ(numHarmonicBondsOrig, numHarmonicBondsConv);

    size_t numHarmonicAnglesOrig = pickType<HarmonicAngle>(interactions).parametersA.size();
    size_t numHarmonicAnglesConv = pickType<HarmonicAngle>(convertedData).parametersA.size();
    EXPECT_EQ(numHarmonicAnglesOrig, numHarmonicAnglesConv);

    HarmonicAngle harmonicAngleOrig0 = pickType<HarmonicAngle>(interactions).parametersA[0];
    HarmonicAngle harmonicAngleConv0 = pickType<HarmonicAngle>(convertedData).parametersA[0];
    EXPECT_TRUE(harmonicAngleOrig0 == harmonicAngleConv0);

    // compare some index data
    auto haIndicesOrig = pickType<HarmonicAngle>(interactions).indices[0];
    auto haIndicesConv = pickType<HarmonicAngle>(convertedData).indices[0];
    EXPECT_EQ(haIndicesOrig, haIndicesConv);
}

TEST(ListedShims, UBInteractionConversion)
{
    gmx_ffparams_t ffparams;
    t_iparams      param;
    param.u_b.kUBA    = 50.0;
    param.u_b.r13A    = 0.1;
    param.u_b.kthetaA = 100.0;
    param.u_b.thetaA  = 0.4;
    ffparams.iparams.push_back(param);
    InteractionDefinitions interactionDefinitions(ffparams);

    constexpr size_t gmxListedID = FindIndex<UreyBradley, GmxToNblibMapping>::value;
    interactionDefinitions.il[gmxListedID].iatoms = { 0, 0, 1, 2 };

    ListedInteractionData convertedData = convertToNblibInteractions(interactionDefinitions);

    size_t numHarmonicBondsConv = pickType<HarmonicBondType>(convertedData).parametersA.size();
    EXPECT_EQ(numHarmonicBondsConv, 1);

    size_t numHarmonicAnglesConv = pickType<HarmonicAngle>(convertedData).parametersA.size();
    EXPECT_EQ(numHarmonicAnglesConv, 1);

    HarmonicBondType harmonicBondConv  = pickType<HarmonicBondType>(convertedData).parametersA[0];
    HarmonicAngle    harmonicAngleConv = pickType<HarmonicAngle>(convertedData).parametersA[0];

    EXPECT_EQ(harmonicBondConv.forceConstant(), param.u_b.kUBA);
    EXPECT_EQ(harmonicBondConv.equilConstant(), param.u_b.r13A);

    EXPECT_EQ(harmonicAngleConv.forceConstant(), param.u_b.kthetaA);
    // Need to convert back to radians
    EXPECT_EQ(harmonicAngleConv.equilConstant() / DEG2RAD, param.u_b.thetaA);

    auto hbIndices = pickType<HarmonicBondType>(convertedData).indices[0];
    auto haIndices = pickType<HarmonicAngle>(convertedData).indices[0];

    EXPECT_EQ(hbIndices[0], 0);
    EXPECT_EQ(hbIndices[1], 2);

    EXPECT_EQ(haIndices[0], 0);
    EXPECT_EQ(haIndices[1], 1);
    EXPECT_EQ(haIndices[2], 2);
}

TEST(NBlibTest, GmxToNblibConversionAllTypes)
{
    ListedInteractionData originalData = combineTestInput(TestInput);

    auto [idef, gmx_params] = convertToGmxInteractions(originalData);
    ignore(gmx_params);

    ListedInteractionData convertedData = convertToNblibInteractions(*idef);

    auto compareParamsAndIndices = [&convertedData](auto& original) {
        if (!original.parametersA.empty())
        {
            using InteractionType = typename std::decay_t<decltype(original)>::type;

            auto converted = pickType<InteractionType>(convertedData);

            // compare parameters
            // comparing the two cosine angles with lower tolerance as this test introduces
            // numerical errors via the x != arccos(cos(x)) comparison
            if constexpr (std::is_same<InteractionType, G96Angle>::value
                          || std::is_same<InteractionType, RestrictedAngle>::value)
            {
                EXPECT_EQ(original.parametersA[0].forceConstant(),
                          converted.parametersA[0].forceConstant());
                real theta_orig = original.parametersA[0].equilConstant();
                real theta_conv = converted.parametersA[0].equilConstant();
                EXPECT_FLOAT_DOUBLE_EQ_TOL(theta_conv,
                                           theta_orig,
                                           theta_orig,
                                           gmx::test::relativeToleranceAsFloatingPoint(theta_orig, 1e-6));
            }
            else
            {
                EXPECT_TRUE(original.parametersA[0] == converted.parametersA[0]);
            }

            // compare indices
            EXPECT_TRUE(original.indices[0] == converted.indices[0]);
        }
    };
    for_each_tuple(compareParamsAndIndices, originalData);
}

//! \brief test forces of all listed types simultaneously
TEST(NBlibTest, EndToEndListedComparison)
{
    int  numParticles    = 4;
    auto interactionData = combineTestInput(TestInput);

    Box                    box(1.0);
    Vec3                   centerOfMass = { 0.5, 0.6, 0.0 };
    std::vector<gmx::RVec> coordinates  = {
        { 1.382, 1.573, 1.482 }, { 1.281, 1.559, 1.596 }, { 1.292, 1.422, 1.663 }, { 1.189, 1.407, 1.775 }
    };

    compareNblibAndGmxListedImplementations(
            interactionData, coordinates, charges, numParticles, 1, box, centerOfMass, 1e-2);
}

template<typename Interaction>
class NblibGmxListed : public testing::Test
{
public:
    void compareNblibAndGmx()
    {
        ListedInteractionData interactionData;

        TypeInput<Interaction> typeInput    = pickType<Interaction>(TestInput);
        size_t                 numParticles = typeInput.coordinates.size();

        pickType<Interaction>(interactionData) = typeInput.interactionData;

        // adjusting the charges vector size to the number of centers in the test case (1 to 4)
        std::vector<real> chargesToSize = subsetVector(charges, numParticles);
        compareNblibAndGmxListedImplementations(interactionData,
                                                typeInput.coordinates,
                                                chargesToSize,
                                                typeInput.coordinates.size(),
                                                1,
                                                Box(1.0),
                                                { 0.5, 0.6, 0.0 },
                                                1e-4);
    }
};

//! \brief test the listed force calculator output for each listed type one-by-one
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
