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
 * Tests COM handling code.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/com.h"

#include <iterator>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbcenums.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

/*! \brief Populates a molectype for generate a graph
 *
 * This moleculetype has the first and last atom not connected to the other atoms
 * so we test optimizations in the shift code that generates the actual graph only
 * for the atom range from the first to the last connected atoms.
 *
 * Defines three residues to test residue and molecule COM handling.
 *
 * \todo Return moltype once t_atoms is a proper class.
 */
void populateMoleculeType(gmx_moltype_t* moltype)
{
    constexpr int atomNumber        = 5;
    constexpr int residueNumber     = 3;
    moltype->atoms.nr               = atomNumber;
    moltype->ilist[F_CONSTR].iatoms = { 0, 1, 2 };
    moltype->ilist[F_ANGLES].iatoms = { 1, 2, 1, 3 };

    snew(moltype->atoms.atom, atomNumber);
    snew(moltype->atoms.resinfo, residueNumber);
    moltype->atoms.atom[0].resind = 0;
    moltype->atoms.atom[1].resind = 1;
    moltype->atoms.atom[2].resind = 1;
    moltype->atoms.atom[3].resind = 1;
    moltype->atoms.atom[4].resind = 2;

    moltype->atoms.atom[0].m = 42;
    moltype->atoms.atom[1].m = 13;
    moltype->atoms.atom[2].m = 7;
    moltype->atoms.atom[3].m = 23;
    moltype->atoms.atom[4].m = 2;
}

//! Set up initial coordinates.
std::vector<RVec> initialCoordinates()
{
    std::vector<RVec> coordinates;
    coordinates.emplace_back(-1, 0, 3);
    coordinates.emplace_back(1.5, 1.5, 4.5);
    coordinates.emplace_back(1.6, 1.5, 4.5);
    coordinates.emplace_back(1.6, 1.7, 4.5);
    coordinates.emplace_back(1, 4, 2);

    coordinates.emplace_back(1, 0, -3);
    coordinates.emplace_back(-1.5, -1.5, -4.5);
    coordinates.emplace_back(-1.6, -1.5, -4.5);
    coordinates.emplace_back(-1.6, -1.7, -4.5);
    coordinates.emplace_back(-1, -4, -2);

    return coordinates;
}

using COMInPlaceTestParams = std::tuple<UnitCellType, CenteringType, PbcType>;
using test::FloatingPointTolerance;

/*! \brief
 * Test fixture for checking correct molecule COM treatment.
 */
class COMInPlaceTest : public ::testing::Test, public ::testing::WithParamInterface<COMInPlaceTestParams>
{
public:
    COMInPlaceTest();

    //! Run the test with the given input for molecule COM.
    void runTestMolecule(const matrix box);
    //! Run the test with the given input for residue COM.
    void runTestResidue(const matrix box);

private:
    //! Coordinates to use.
    std::vector<RVec> testCoordinates_;
    //! Dummy topology to use.
    gmx_mtop_t testTopology_;
    //! Reference data holder.
    TestReferenceData data_;
    //! Checker for data.
    TestReferenceChecker checker_;
};

COMInPlaceTest::COMInPlaceTest() :
    testCoordinates_(initialCoordinates()), checker_(data_.rootChecker())
{
    auto& moltype = testTopology_.moltype.emplace_back();
    populateMoleculeType(&moltype);
    auto& molblock       = testTopology_.molblock.emplace_back();
    molblock.nmol        = 2;
    molblock.type        = 0;
    testTopology_.natoms = moltype.atoms.nr * molblock.nmol;
    testTopology_.finalize();
    FloatingPointTolerance tolerance(
            FloatingPointTolerance(1.0e-6, 1.0e-6, 1.0e-8, 1.0e-12, 10000, 100, false));
    checker_.setDefaultTolerance(tolerance);
}

void COMInPlaceTest::runTestMolecule(const matrix box)
{
    auto params   = GetParam();
    auto unitcell = std::get<0>(params);
    auto center   = std::get<1>(params);
    auto pbcType  = std::get<2>(params);
    placeCoordinatesWithCOMInBox(
            pbcType, unitcell, center, box, testCoordinates_, testTopology_, COMShiftType::Molecule);
    std::string testString = "Molecule " + std::string(unitCellTypeNames(unitcell))
                             + std::string(centerTypeNames(center)) + c_pbcTypeNames[pbcType]
                             + c_pbcTypeNames[guessPbcType(box)];
    checker_.checkSequence(std::begin(testCoordinates_), std::end(testCoordinates_), testString.c_str());
}

void COMInPlaceTest::runTestResidue(const matrix box)
{
    auto params   = GetParam();
    auto unitcell = std::get<0>(params);
    auto center   = std::get<1>(params);
    auto pbcType  = std::get<2>(params);
    placeCoordinatesWithCOMInBox(
            pbcType, unitcell, center, box, testCoordinates_, testTopology_, COMShiftType::Residue);
    std::string testString = "Residue " + std::string(unitCellTypeNames(unitcell))
                             + std::string(centerTypeNames(center)) + c_pbcTypeNames[pbcType]
                             + c_pbcTypeNames[guessPbcType(box)];
    checker_.checkSequence(std::begin(testCoordinates_), std::end(testCoordinates_), testString.c_str());
}

TEST(ShiftTest, CoordinateShiftWorks)
{
    auto       coords = initialCoordinates();
    const RVec shift(1.5, -2.2, 0);

    shiftAtoms(shift, coords);
    EXPECT_FLOAT_EQ(coords[4][0], 2.5);
    EXPECT_FLOAT_EQ(coords[3][1], -0.5);
    EXPECT_FLOAT_EQ(coords[9][2], -2);
}


TEST_P(COMInPlaceTest, MatrixDefault)
{
    const matrix box = { { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 } };

    runTestMolecule(box);
    runTestResidue(box);
}

INSTANTIATE_TEST_SUITE_P(CorrectCoordinates,
                         COMInPlaceTest,
                         testing::Combine(::testing::Values(UnitCellType::Compact,
                                                            UnitCellType::Rectangular,
                                                            UnitCellType::Triclinic),
                                          ::testing::Values(CenteringType::Rectangular,
                                                            CenteringType::Triclinic,
                                                            CenteringType::Zero),
                                          ::testing::Values(PbcType::No, PbcType::Xyz, PbcType::XY)));

// TODO add PbcType::Screw once it is fully supported.

} // namespace

} // namespace test

} // namespace gmx
