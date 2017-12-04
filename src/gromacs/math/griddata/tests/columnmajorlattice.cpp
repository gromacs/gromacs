/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests column major lattice indexing routines
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/griddata/columnmajorlattice.h"

#include <vector>
#include <string>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace internal
{

namespace
{

TEST(ColumnMajorLatticeTest, canConstruct)
{
    ColumnMajorLattice<1> minimalLAttice {{{
                                               1
                                           }}};
    ColumnMajorLattice<2> twoDLattice {{{
                                            4, 2
                                        }}};
    ColumnMajorLattice<4> fourDLattice {{{
                                             4, 2, 78, 1
                                         }}};
}

TEST(ColumnMajorLatticeTest, numberLatticePoints)
{
    //TODO: implement and test overflow behaviour
    ASSERT_EQ(ColumnMajorLattice<4>({{{ 4, 2, 78, 5}}}).getNumLatticePoints(), 3120);
}

TEST(ColumnMajorLatticeTest, lineariseIsInverseOfVectorise)
{
    ColumnMajorLattice<1> onedLattice {{{
                                            1
                                        }}};

    int oneDLinearIndex {
        0
    };
    ASSERT_EQ(oneDLinearIndex, onedLattice.lineariseVectorIndex(onedLattice.vectoriseLinearIndex(oneDLinearIndex)));

    // Test all indices in three-dimensional grid
    ColumnMajorLattice<3> threeDimensionalLattice {{{
                                                        2, 3, 4
                                                    }}};

    std::vector<int> threeDLinearIndices(threeDimensionalLattice.getNumLatticePoints());
    std::iota (std::begin(threeDLinearIndices), std::end(threeDLinearIndices), 0);
    for (const auto &i : threeDLinearIndices)
    {
        ASSERT_EQ(i, threeDimensionalLattice.lineariseVectorIndex(threeDimensionalLattice.vectoriseLinearIndex(i)));
    }
}

TEST(ColumnMajorLatticeTest, vectoriseIsInverseOfLinearise)
{
    // Test all indices in three-dimensional grid
    ColumnMajorLattice<3> threeDimensionalLattice {{{
                                                        2, 3, 4
                                                    }}};

    std::vector<ColumnMajorLattice<3>::MultiIndex> threeDVectorIndices;
    std::vector<int> x_indices = {0, 1};
    std::vector<int> y_indices = {0, 1, 2};
    std::vector<int> z_indices = {0, 1, 2, 3};

    for (const auto &ix : x_indices)
    {
        for (const auto &iy : y_indices)
        {
            for (const auto &iz : z_indices)
            {
                threeDVectorIndices.push_back({{ix, iy, iz}});
            }
        }
    }
    for (const auto &i : threeDVectorIndices)
    {
        ASSERT_EQ(i, threeDimensionalLattice.vectoriseLinearIndex(threeDimensionalLattice.lineariseVectorIndex(i)));
    }
}

TEST(ColumnMajorLatticeTest, throwsWhenOutOfBounds)
{
    ColumnMajorLattice<3> threeDimensionalLattice {{{
                                                        2, 3, 4
                                                    }}};
    EXPECT_THROW_GMX(threeDimensionalLattice.lineariseVectorIndex({{-1, 0, 0}}), gmx::RangeError);
    EXPECT_THROW_GMX(threeDimensionalLattice.lineariseVectorIndex({{3, 1, 1}}), gmx::RangeError);
    EXPECT_THROW_GMX(threeDimensionalLattice.vectoriseLinearIndex(-1), gmx::RangeError);
    EXPECT_THROW_GMX(threeDimensionalLattice.vectoriseLinearIndex(24), gmx::RangeError);

}


} // namespace

} // internal

} // test

} // gmx
