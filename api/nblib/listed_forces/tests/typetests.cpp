/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * This implements basic nblib box tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gromacs/utility/arrayref.h"
#include "nblib/listed_forces/dataflow.hpp"
#include "nblib/listed_forces/tests/listedtesthelpers.h"
#include "nblib/tests/testhelpers.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace nblib
{

//! Coordinates for testing
static const std::vector<gmx::RVec> c_coordinatesForDihTests = { { 0.0, 0.0, 0.0 },
                                                                 { 0.0, 0.0, 0.2 },
                                                                 { 0.005, 0.0, 0.1 },
                                                                 { -0.001, 0.1, 0.0 } };

//! Coordinates for testing angles
static const std::vector<gmx::RVec> c_coordinatesForAngleTests = { { 1.382, 1.573, 1.482 },
                                                                   { 1.281, 1.559, 1.596 },
                                                                   { 1.292, 1.422, 1.663 } };

//! Coordinates for testing bonds
static const std::vector<gmx::RVec> c_coordinatesForBondTests = { { 1.382, 1.573, 1.482 },
                                                                  { 1.281, 1.559, 1.596 } };

//! Function types for testing Harmonic bonds
static const std::vector<HarmonicBondType> c_InputHarmonicBonds = { { HarmonicBondType(500.0, 0.15) } };

//! Function types for testing G96 bonds
static const std::vector<G96BondType> c_InputG96Bonds = { { G96BondType(50.0, 0.15) } };

//! Function types for testing cubic bonds
static const std::vector<CubicBondType> c_InputCubicBonds = { { CubicBondType(50.0, 2.0, 0.16) } };

//! Function types for testing FENE bonds
static const std::vector<FENEBondType> c_InputFeneBonds = { { FENEBondType(5.0, 0.4) } };

//! Function types for testing Harmonic angles
static const std::vector<HarmonicAngle> c_InputHarmonicAngles = { { HarmonicAngle(50.0, Degrees(100)) } };

//! Function types for testing dihedrals
static const std::vector<ProperDihedral> c_InputDihs = { { ProperDihedral(Degrees(-105.0), 15.0, 2) } };

template<class Interaction, std::enable_if_t<Contains<Interaction, SupportedListedTypes>{}>* = nullptr>
void checkForcesAndEnergiesWithRefData(std::vector<Interaction> input, gmx::ArrayRef<const gmx::RVec> x)
{
    auto                   indices = indexVector<Interaction>();
    PbcHolder              pbcHolder(PbcType::Xyz, Box(1.5));
    test::RefDataChecker   refDataChecker(1e-4);
    std::vector<gmx::RVec> forces(x.size(), gmx::RVec{ 0, 0, 0 });

    auto energy = computeForces(gmx::ArrayRef<const InteractionIndex<Interaction>>(indices),
                                gmx::ArrayRef<const Interaction>(input),
                                x,
                                &forces,
                                pbcHolder);

    refDataChecker.testReal(energy, "Epot");
    refDataChecker.testArrays<gmx::RVec>(forces, "forces");
}

TEST(FourCenter, ListedForcesProperDihedralTest)
{
    checkForcesAndEnergiesWithRefData(c_InputDihs, c_coordinatesForDihTests);
}

TEST(ThreeCenter, ListedForcesHarmonicAngleTest)
{
    checkForcesAndEnergiesWithRefData(c_InputHarmonicAngles, c_coordinatesForAngleTests);
}

TEST(TwoCenter, ListedForcesHarmonicBondTest)
{
    checkForcesAndEnergiesWithRefData(c_InputHarmonicBonds, c_coordinatesForBondTests);
}

TEST(TwoCenter, ListedForcesG96BondTest)
{
    checkForcesAndEnergiesWithRefData(c_InputG96Bonds, c_coordinatesForBondTests);
}

TEST(TwoCenter, ListedForcesCubicBondTest)
{
    checkForcesAndEnergiesWithRefData(c_InputCubicBonds, c_coordinatesForBondTests);
}

TEST(TwoCenter, ListedForcesFeneBondTest)
{
    checkForcesAndEnergiesWithRefData(c_InputFeneBonds, c_coordinatesForBondTests);
}

} // namespace nblib
