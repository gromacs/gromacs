/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Tests matrices
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/matrix.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

namespace
{

class MatrixTest : public ::testing::Test
{
    public:
        MatrixTest()
        {
            std::fill(begin(matrix_), end(matrix_), testNumber_ - 1);
        }
    protected:
        Matrix3x3 matrix_;
        real      testNumber_ = 42;
};


TEST_F(MatrixTest, canConstructAndFill)
{
    for (const auto &x : matrix_)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MatrixTest, canSetValues)
{
    matrix_(1, 1) = testNumber_;
    EXPECT_EQ(testNumber_, matrix_(1, 1) );
}

TEST_F(MatrixTest, canCopyAssign)
{
    Matrix3x3 other;
    other = matrix_;
    auto      twoDArrayIt = begin(matrix_);
    for (const auto &x : other)
    {
        EXPECT_EQ(*twoDArrayIt, x);
        ++twoDArrayIt;
    }
}

TEST_F(MatrixTest, canSwap)
{
    Matrix3x3 other;
    other.swap(matrix_);
    for (const auto &x : other)
    {
        EXPECT_EQ(testNumber_ - 1, x);
    }
}

TEST_F(MatrixTest, staticMultiDimArrayExtent)
{
    EXPECT_EQ(matrix_.extent(0), DIM);
    EXPECT_EQ(matrix_.extent(1), DIM);
}

} // namespace

} // namespace test

} // namespace gmx
