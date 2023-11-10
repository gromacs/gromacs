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
 * This implements basic nblib box tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <gtest/gtest.h>

#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/dataflowpolarization.hpp"
#include "nblib/tests/testhelpers.h"
#include "nblib/util/array.hpp"

#include "listedtesthelpers.h"
#include "listedtypeinput.hpp"

namespace nblib
{

/*
// This test uses charges and coordinates from the GROMACS test and not from the listedtypeinput.hpp
static void compareForcesAndEnergy(TypeInput<PairLJType> input)
{
    gmx::ArrayRef<const PairLJType> paramsA(input.interactionData.parametersA);

    PbcHolder            pbcHolder(PbcType::Xyz, Box(1.0));
    test::RefDataChecker refDataChecker(1e-4, "CompareTypes_ChargedLJPair.xml");

    KernelEnergy<real> energy;

    std::vector<real> charge = { 1.0, -0.5, -0.5 };

    std::vector<util::array<real, 3>> chargedPairCoordinates = { { 0.0, 0.0, 0.0 },
                                                                 { 1.0, 1.0, 1.0 },
                                                                 { 1.1, 1.2, 1.3 } };

    std::vector<util::array<real, 3>> forces(chargedPairCoordinates.size(),
                                             util::array<real, 3>{ 0, 0, 0 });

    std::vector<InteractionIndex<PairLJType>> chargedPairIndices = { { 1, 2, 0 }, { 0, 2, 0 } };

    energy = computePolarizationForces(
            gmx::ArrayRef<const InteractionIndex<PairLJType>>(chargedPairIndices),
            paramsA,
            paramsA,
            gmx::ArrayRef<const util::array<real, 3>>(chargedPairCoordinates),
            gmx::ArrayRef<real>(charge),
            NoFepLambdaType{},
            &forces,
            gmx::ArrayRef<std::nullptr_t>{},
            pbcHolder);

    refDataChecker.testReal(energy.eVdw(), "EVdw");
    refDataChecker.testReal(energy.eCoul(), "ECoul");
    refDataChecker.test3DVectors<util::array<real, 3>>(forces, "forces");
}
 */

template<class Interaction>
void compareForcesAndEnergy(TypeInput<Interaction> input)
{
    gmx::ArrayRef<const InteractionIndex<Interaction>> indices(input.interactionData.indices);
    gmx::ArrayRef<const Interaction>                   paramsA(input.interactionData.parametersA);
    gmx::ArrayRef<const gmx::RVec>                     x(input.coordinates);

    PbcHolder                         pbcHolder(PbcType::Xyz, Box(1.5));
    test::RefDataChecker              refDataChecker(1e-4, ("CompareTypes_" + input.name + ".xml"));
    std::vector<util::array<real, 3>> forces(input.coordinates.size(), util::array<real, 3>{ 0, 0, 0 });

    KernelEnergy<real> energy;

    if constexpr (Contains<Interaction, BasicListedTypes>{})
    {
        energy = computeForces(indices, paramsA, x, &forces, pbcHolder);
    }
    if constexpr (Contains<Interaction, PolarizationTypes>{})
    {
        std::vector<real> charge = { 1.5, -2.0 };
        energy                   = computePolarizationForces(
                indices, paramsA, paramsA, x, charge, NoFepLambdaType{}, &forces, gmx::ArrayRef<std::nullptr_t>{}, pbcHolder);
    }

    refDataChecker.testReal(energy.potentialEnergy(), "Epot");
    refDataChecker.test3DVectors<util::array<real, 3>>(forces, "forces");
}

// Tests for position restraints implemented in positionrestraints.cpp
static void compareForcesAndEnergy(TypeInput<PositionRestraints> input)
{
    std::ignore = (input);
}

template<class Interaction>
class TypeTests : public testing::Test
{
public:
    void compareWithRefData()
    {
        TypeInput<Interaction> typeInput = pickType<Interaction>(TestInput);
        compareForcesAndEnergy(typeInput);
    }
};

//! \brief test the listed force calculator output for each listed type one-by-one
TYPED_TEST_SUITE_P(TypeTests);

TYPED_TEST_P(TypeTests, RefComp)
{
    this->compareWithRefData();
}

REGISTER_TYPED_TEST_SUITE_P(TypeTests, RefComp);

INSTANTIATE_TYPED_TEST_SUITE_P(CompareTypes, TypeTests, TestTypes);

} // namespace nblib
