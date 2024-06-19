/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests PBC code
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/pbc.h"

#include "config.h"

#include <cfenv>
#include <cstddef>

#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

/*! \brief Small value used for testing points close to edges
 *
 * The value is defined a bit arbitrarily. It should be small enough to test edge cases
 * but not so small that the test fails because of rounding errors.
 */
static real returnSmallValue(real L)
{
    return 100 * std::numeric_limits<real>::epsilon() * L;
}

//! Relative tolerance for the expected floating point rounding error when wrapping large numbers
constexpr real relativeTolerance = 0.5 * std::numeric_limits<real>::epsilon();

TEST(PbcTest, CalcShiftsWorks)
{
    // Choose box vector entries whose magnitudes will lead to unique
    // shift vector values when the largest box shift in any dimension
    // is two.
    std::vector<gmx::RVec> shiftVectors(c_numShiftVectors);

    const matrix box = { { 0.01, 1, -100 }, { 300, -0.03, 3 }, { -6, -600, 0.06 } };
    calc_shifts(box, shiftVectors);

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    checker.checkSequence(std::begin(shiftVectors), std::end(shiftVectors), "ShiftVectors");
}

//! Testing if points that are already inside a cubic box are mapped back to themselves
TEST(PbcTest, PutAtomsInCubicBoxAlreadyInBox)
{
    // Testing cubic box with edge length 3.5 nm
    const real                   smallValue = returnSmallValue(3.5);
    const FloatingPointTolerance scTolerance = relativeToleranceAsFloatingPoint(3.5, relativeTolerance);
    const matrix                 box         = { { 3.5, 0, 0 }, { 0, 3.5, 0 }, { 0, 0, 3.5 } };
    const real                   coordCloseToEdge = 3.5 - smallValue;
    std::vector<gmx::RVec>       pointsPrePbc{ { 1.0, 2.0, 3.0 },
                                         { smallValue, 1.0, 1.0 },
                                         { coordCloseToEdge, coordCloseToEdge, coordCloseToEdge } };
    std::vector<gmx::RVec>       pointsPostPbc(pointsPrePbc);
    put_atoms_in_box(PbcType::Xyz, box, pointsPrePbc);
    for (std::size_t i = 0; i < pointsPrePbc.size(); ++i)
    {
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][0], pointsPostPbc[i][0], scTolerance);
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][1], pointsPostPbc[i][1], scTolerance);
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][2], pointsPostPbc[i][2], scTolerance);
    }
}

//! Testing if points outside a cubic box are mapped to their periodic image inside the box
TEST(PbcTest, PutAtomsInCubicBoxFromOutsideBox)
{
    // Testing cubic box with edge length 3.5 nm
    const real                   smallValue = returnSmallValue(3.5);
    const FloatingPointTolerance scTolerance = relativeToleranceAsFloatingPoint(3.5, relativeTolerance);
    const matrix                 box         = { { 3.5, 0, 0 }, { 0, 3.5, 0 }, { 0, 0, 3.5 } };
    const real                   outsideCloseToEdge = 3.5 + smallValue;
    const real                   insideCloseToEdge  = 3.5 - smallValue;
    std::vector<gmx::RVec>       pointsPrePbc{ { -0.5, 1.0, 4.5 },
                                         { 10000.0, 20000.0, -40000.0 },
                                         { -smallValue, 1.0, 1.0 },
                                         { outsideCloseToEdge, outsideCloseToEdge, outsideCloseToEdge } };
    std::vector<gmx::RVec>       pointsPostPbc{ { 3.0, 1.0, 1.0 },
                                          { 0.5, 1.0, 1.5 },
                                          { insideCloseToEdge, 1.0, 1.0 },
                                          { smallValue, smallValue, smallValue } };
    put_atoms_in_box(PbcType::Xyz, box, pointsPrePbc);
    for (std::size_t i = 0; i < pointsPrePbc.size(); ++i)
    {
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][0], pointsPostPbc[i][0], scTolerance);
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][1], pointsPostPbc[i][1], scTolerance);
        EXPECT_REAL_EQ_TOL(pointsPrePbc[i][2], pointsPostPbc[i][2], scTolerance);
    }
}

// TODO: test edge cases for triclinic boxes
//! Testing if a point ouside a triclinic box is mapped to its periodic image inside the box
TEST(PbcTest, PutAtomsInTriclinicBoxFromOutsideBox)
{
    std::vector<gmx::RVec>       pointOutside{ { -1.0, 1.0, 4.0 } };
    const matrix                 box         = { { 3.5, 0, 0 }, { 0, 3.5, 0 }, { 1.5, 0, 3.5 } };
    const FloatingPointTolerance scTolerance = relativeToleranceAsFloatingPoint(3.5, relativeTolerance);
    put_atoms_in_box(PbcType::Xyz, box, pointOutside);
    EXPECT_REAL_EQ_TOL(pointOutside[0][0], 1.0, scTolerance);
    EXPECT_REAL_EQ_TOL(pointOutside[0][1], 1.0, scTolerance);
    EXPECT_REAL_EQ_TOL(pointOutside[0][2], 0.5, scTolerance);
}

#if HAVE_FEDISABLEEXCEPT && !defined(__riscv)
/*! \brief Testing if put_atoms_in_box enters into an infinite loop
 *
 * This test checks that put_atoms_in_box does not enter into an infinite loop in case floating
 * point exceptions are not enabled. The test is not run in case floating point exceptions are
 * enabled (FE_INVALID). Untested on riscv architectures.
 */
TEST(PbcTest, PutAtomsInBoxHandlesInf)
{
    if (fegetexcept() & FE_INVALID)
    {
        GTEST_SKIP() << "Cannot test if floating point exceptions are enabled (FE_INVALID)";
    }
    else
    {
        const matrix           box           = { { 3.5, 0, 0 }, { 0, 3.5, 0 }, { 0, 0, 3.5 } };
        real                   infiniteValue = std::numeric_limits<real>::infinity();
        std::vector<gmx::RVec> points{ { infiniteValue, infiniteValue, infiniteValue } };
        put_atoms_in_box(PbcType::Xyz, box, points);
    }
}
#endif

} // namespace test

} // namespace gmx
