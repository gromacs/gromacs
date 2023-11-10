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
 * Input data and processing for the tests that compare nblib's listed
 * forces implementation. For use with both GROMACS implementation and
 * ref data comparisons.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_LISTEDFORCES_LISTEDTYPEINPUT_HPP
#define NBLIB_LISTEDFORCES_LISTEDTYPEINPUT_HPP

#include "nblib/listed_forces/definitions.h"
#include "nblib/vector.h"

namespace nblib
{

static std::vector<gmx::RVec> testCoordinates = { { 1.382, 1.573, 1.482 },
                                                  { 1.281, 1.559, 1.596 },
                                                  { 1.292, 1.422, 1.663 },
                                                  { 1.189, 1.407, 1.775 } };

static std::vector<real> charges = { 1.5, -2.0, 1.5, -1.0 };

template<class Interaction>
struct TypeInput
{
    typedef Interaction         type; // needed for pickType
    ListedTypeData<Interaction> interactionData;
    std::string                 name;
    // assign  coordinates depending on the number of centers in the interaction type from the array above
    std::vector<gmx::RVec> coordinates = subsetVector(testCoordinates, NCenter<Interaction>{});
};

static std::tuple TestInput{
    // One center Types
    TypeInput<PositionRestraints>{ { { { PositionRestraints(0.5, 0.6, 0.0, 0, 200, 400) } },
                                     {},
                                     { indexVector<PositionRestraints>() } },
                                   "PositionRestraints" },
    // Two Center Types
    TypeInput<HarmonicBondType>{
            { { { HarmonicBondType(500.0, 0.15) } }, {}, { indexVector<HarmonicBondType>() } },
            "HarmonicBond" },
    TypeInput<G96BondType>{ { { { G96BondType(50.0, 0.15) } }, {}, { indexVector<G96BondType>() } },
                            "G96Bond" },
    TypeInput<CubicBondType>{ { { { CubicBondType(50.0, 2.0, 0.16) } }, {}, { indexVector<CubicBondType>() } },
                              "CubicBond" },
    TypeInput<MorseBondType>{ { { { MorseBondType(30.0, 2.7, 0.15) } }, {}, { indexVector<MorseBondType>() } },
                              "MorseBond" },
    TypeInput<FENEBondType>{ { { { FENEBondType(5.0, 0.4) } }, {}, { indexVector<FENEBondType>() } },
                             "FENEBond" },
    // Polarization Types
    TypeInput<SimplePolarization>{
            { { { SimplePolarization(0.12) } }, {}, { indexVector<SimplePolarization>() } },
            "SimplePolarization" },
    // TypeInput<PairLJType>{
    //    { { { PairLJType(C6(0.001458), C12(1.0062882e-6)) } }, {}, indexVector<PairLJType>() }, "ChargedLJPair" },
    // Three Center Types
    TypeInput<HarmonicAngle>{
            { { { HarmonicAngle(50.0, Degrees(100)) } }, {}, { indexVector<HarmonicAngle>() } },
            "HarmonicAngle" },
    TypeInput<G96Angle>{ { { { G96Angle(50.0, Degrees(100)) } }, {}, { indexVector<G96Angle>() } },
                         "CosineAngle" },
    TypeInput<RestrictedAngle>{
            { { { RestrictedAngle(50.0, Degrees(100)) } }, {}, { indexVector<RestrictedAngle>() } },
            "RestrictedAngle" },
    TypeInput<LinearAngle>{ { { { LinearAngle(50.0, 0.4) } }, {}, { indexVector<LinearAngle>() } },
                            "LinearAngle" },
    TypeInput<QuarticAngle>{ { { { QuarticAngle(1.1, 2.3, 4.6, 7.8, 9.2, Degrees(87)) } },
                               {},
                               { indexVector<QuarticAngle>() } },
                             "QuarticAngle" },
    TypeInput<CrossBondBond>{ { { { CrossBondBond(45.0, 0.8, 0.7) } }, {}, { indexVector<CrossBondBond>() } },
                              "CrossBondBond" },
    TypeInput<CrossBondAngle>{
            { { { CrossBondAngle(45.0, 0.8, 0.7, 0.3) } }, {}, { indexVector<CrossBondAngle>() } },
            "CrossBondAngle" },
    // Four Center Types
    TypeInput<ProperDihedral>{
            { { { ProperDihedral(Degrees(-105), 15.0, 2) } }, {}, { indexVector<ProperDihedral>() } },
            "ProperDihedral" },
    TypeInput<ImproperDihedral>{
            { { { ImproperDihedral(Degrees(100), 50) } }, {}, { indexVector<ImproperDihedral>() } },
            "ImproperDihedral" },
    TypeInput<RyckaertBellemanDihedral>{ { { { RyckaertBellemanDihedral(-7.35, 13.6, 8.4, -16.7, 1.3, 12.4) } },
                                           {},
                                           { indexVector<RyckaertBellemanDihedral>() } },
                                         "RBDihedral" }
};

//! \brief Converts the input data tuple above into a ListedInteractionData package that can be used as input for ListedForceCalculator
template<class... Ts>
ListedInteractionData combineTestInput(std::tuple<Ts...> testInput)
{
    ListedInteractionData interactionData;
    // transfer all elements of testInput into the returned ListedInteractionData
    // use a lambda + for_each_tuple
    auto copyParamsOneType = [&interactionData](const auto& typeInput) {
        for (size_t i = 0; i < typeInput.interactionData.parametersA.size(); i++)
        {
            auto interactionParams = typeInput.interactionData.parametersA[i];
            using InteractionType  = decltype(interactionParams);

            pickType<InteractionType>(interactionData).parametersA.push_back(interactionParams);
            pickType<InteractionType>(interactionData).parametersB.push_back(interactionParams);

            auto indices = typeInput.interactionData.indices[i];
            pickType<InteractionType>(interactionData).indices.push_back(indices);
        }
    };
    for_each_tuple(copyParamsOneType, testInput);

    return interactionData;
}

// Meta-function to return the ::type member of the argument
template<class TestInputofInteractionType>
using ExtractType = typename TestInputofInteractionType::type;

// a TypeList with all the listed interaction types that occur in TestInput
using ListedTypes = Map<ExtractType, decltype(TestInput)>;
// put the listed types into a gtest type-list
using TestTypes = Reduce<::testing::Types, ListedTypes>;

} // namespace nblib

#endif // NBLIB_LISTEDFORCES_LISTEDTYPEINPUT_HPP
