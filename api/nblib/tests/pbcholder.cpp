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
 * This PbcHolder tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include <iostream>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

#include "nblib/pbc.hpp"
#include "nblib/util/array.hpp"

using gmx::test::defaultRealTolerance;

namespace nblib
{

TEST(NBlibTest, PbcHolderWorks)
{
    Box box(10, 10, 10);

    PbcHolder pbcHolder(PbcType::Xyz, box);

    gmx::RVec x1{ 1.0, 1.1, 0.9 }, x2{ 9, 8.9, 9.1 };
    gmx::RVec dx;

    pbcHolder.dxAiuc(x1, x2, dx);
    gmx::RVec ref{ 2, 2.2, 1.8 };

    EXPECT_REAL_EQ_TOL(ref[0], dx[0], gmx::test::relativeToleranceAsFloatingPoint(ref[0], 1e-6));
    EXPECT_REAL_EQ_TOL(ref[1], dx[1], gmx::test::relativeToleranceAsFloatingPoint(ref[0], 1e-6));
    EXPECT_REAL_EQ_TOL(ref[2], dx[2], gmx::test::relativeToleranceAsFloatingPoint(ref[0], 1e-6));
}

inline util::array<int, 3> shiftIndexToXyz(int s)
{
    int ix = s % gmx::detail::c_nBoxX;
    int iy = ((s - ix) % (gmx::detail::c_nBoxY * gmx::detail::c_nBoxX)) / gmx::detail::c_nBoxX;
    int iz = (s - ix - iy) / (gmx::detail::c_nBoxY * gmx::detail::c_nBoxX);

    return { ix - gmx::c_dBoxX, iy - gmx::c_dBoxY, iz - gmx::c_dBoxZ };
}

TEST(NBlibTest, PbcShiftIndex)
{
    Box       box(10, 10, 10);
    PbcHolder pbc(PbcType::Xyz, box);

    gmx::RVec x1{ 1, 0, 0 }, x2{ 9, 0, 0 };
    gmx::RVec dx;

    int s1 = pbc.dxAiuc(x1, x2, dx);
    int s2 = pbc.dxAiuc(x2, x1, dx);

    auto xyz      = shiftIndexToXyz(s1);
    int  s2decode = gmx::xyzToShiftIndex(-xyz[0], -xyz[1], -xyz[2]);
    EXPECT_EQ(s2, s2decode);

    std::vector<util::array<int, 3>> ivecs{
        { 0, 0, 0 },   { 1, 0, 0 }, { -1, 0, 0 },  { 0, 1, 0 },  { 0, -1, 0 },  { 0, 0, 0 },
        { 0, 0, -1 },  { 2, 0, 0 }, { -2, 0, 0 },  { 2, 0, 1 },  { -2, 0, 1 },  { 2, 1, 0 },
        { -2, 1, 0 },  { 2, 1, 1 }, { -2, 1, 1 },  { 2, -1, 0 }, { -2, -1, 0 }, { 2, -1, 1 },
        { -2, -1, 1 }, { 1, 1, 1 }, { -1, -1, -1 }
    };

    for (auto iv : ivecs)
    {
        int sIdx = gmx::xyzToShiftIndex(iv[0], iv[1], iv[2]);

        auto decoded = shiftIndexToXyz(sIdx);
        EXPECT_EQ(iv, decoded);
    }
}

} // namespace nblib
