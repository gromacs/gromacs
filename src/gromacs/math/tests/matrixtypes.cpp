/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Tests gmx::Matrix3x3 implementation.
 *
 * The main point of these tests is to check that all different constructs
 * using gmx::Matrix3x3 compile, and that some of the non-trivial conversions
 * to/from matrix work as intended.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/matrixtypes.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"

namespace
{

using gmx::Matrix3x3;
using gmx::Matrix3x3Upper;
using gmx::Matrix3x3Lower;

TEST(Matrix3x3Test, CanBeStoredInVector)
{
    std::vector<Matrix3x3> v;
    v.emplace_back(1, 2, 3, 4, 5, 6, 7, 8, 9);
    v.resize(2);
    EXPECT_EQ(1, v[0][XX][XX]);
    EXPECT_EQ(2, v[0][XX][YY]);
    EXPECT_EQ(3, v[0][XX][ZZ]);
    EXPECT_EQ(4, v[0][YY][XX]);
    EXPECT_EQ(5, v[0][YY][YY]);
    EXPECT_EQ(6, v[0][YY][ZZ]);
    EXPECT_EQ(7, v[0][ZZ][XX]);
    EXPECT_EQ(8, v[0][ZZ][YY]);
    EXPECT_EQ(9, v[0][ZZ][ZZ]);
}

TEST(Matrix3x3Test, ConvertsImplicitlyFrom_matrix)
{
    std::vector<Matrix3x3> v;
    matrix                 x = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    v.emplace_back(x);
    EXPECT_EQ(1, v[0][XX][XX]);
    EXPECT_EQ(2, v[0][XX][YY]);
    EXPECT_EQ(3, v[0][XX][ZZ]);
    EXPECT_EQ(4, v[0][YY][XX]);
    EXPECT_EQ(5, v[0][YY][YY]);
    EXPECT_EQ(6, v[0][YY][ZZ]);
    EXPECT_EQ(7, v[0][ZZ][XX]);
    EXPECT_EQ(8, v[0][ZZ][YY]);
    EXPECT_EQ(9, v[0][ZZ][ZZ]);
}

TEST(Matrix3x3Test, ConvertsImplicitlyTo_matrix)
{
    std::vector<Matrix3x3> v;
    v.emplace_back(1, 2, 3, 4, 5, 6, 7, 8, 9);
    matrix                 x;
    copy_mat(v[0], x);
    EXPECT_EQ(1, x[XX][XX]);
    EXPECT_EQ(2, x[XX][YY]);
    EXPECT_EQ(3, x[XX][ZZ]);
    EXPECT_EQ(4, x[YY][XX]);
    EXPECT_EQ(5, x[YY][YY]);
    EXPECT_EQ(6, x[YY][ZZ]);
    EXPECT_EQ(7, x[ZZ][XX]);
    EXPECT_EQ(8, x[ZZ][YY]);
    EXPECT_EQ(9, x[ZZ][ZZ]);
}

TEST(Matrix3x3Test, WorksAsMutable_matrix)
{
    std::vector<Matrix3x3> v;
    v.emplace_back(1, 2, 3, 4, 5, 6, 7, 8, 9);
    matrix                 x = {{4, 5, 6}, {7, 8, 9}, {10, 11, 12}};
    copy_mat(x, v[0]);
    EXPECT_EQ(4, v[0][XX][XX]);
    EXPECT_EQ(5, v[0][XX][YY]);
    EXPECT_EQ(6, v[0][XX][ZZ]);
    EXPECT_EQ(7, v[0][YY][XX]);
    EXPECT_EQ(8, v[0][YY][YY]);
    EXPECT_EQ(9, v[0][YY][ZZ]);
    EXPECT_EQ(10, v[0][ZZ][XX]);
    EXPECT_EQ(11, v[0][ZZ][YY]);
    EXPECT_EQ(12, v[0][ZZ][ZZ]);
}

TEST(Matrix3x3UpperTest, ConvertsImplicitlyFrom_matrix)
{
    std::vector<Matrix3x3Upper> v;
    matrix                      x = {{1, 2, 3}, {0, 5, 6}, {0, 0.0, 9}};
    v.emplace_back(x);
    EXPECT_EQ(1, v[0][XX][XX]);
    EXPECT_EQ(2, v[0][XX][YY]);
    EXPECT_EQ(3, v[0][XX][ZZ]);
    EXPECT_EQ(0, v[0][YY][XX]);
    EXPECT_EQ(5, v[0][YY][YY]);
    EXPECT_EQ(6, v[0][YY][ZZ]);
    EXPECT_EQ(0, v[0][ZZ][XX]);
    EXPECT_EQ(0, v[0][ZZ][YY]);
    EXPECT_EQ(9, v[0][ZZ][ZZ]);
}

TEST(Matrix3x3UpperTest, UpholdsTheConstraintOnConstruction)
{
    std::vector<Matrix3x3Upper> v;
    matrix                      x = {{1, 2, 3}, {0, 5, 6}, {7, 0.0, 9}};
    EXPECT_DEATH(v.emplace_back(x), "The matrix constraint was violated");
}

TEST(Matrix3x3UpperTest, UpholdsTheConstraintOnAssignment)
{
    std::vector<Matrix3x3Upper> v;
    matrix                      x = {{1, 0, 0}, {0, 5, 0}, {0, 0.0, 9}};
    v.emplace_back(x);
    EXPECT_DEATH(v[0].setValue(YY, XX, 7), "The matrix constraint was violated");
}

//! Test Matrix3x3 function for testing type conversion. Could be a lazy copy_mat() replacement.
void testMatrixCopy(const Matrix3x3 &src, Matrix3x3 &dst)
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            dst.setValue(i, j, src[i][j]);
        }
    }
}

TEST(Matrix3x3UpperTest, ConvertsImplicitlyToMatrix3x3)
{
    // This tests that Matrix3x3Upper can be passed as Matrix3x3 &.
    // The opposite (passing Matrix3x3 as Matrix3x3Upper &) will not compile.
    std::vector<Matrix3x3>      v2;
    matrix                      x = {{1, 2, 3}, {0, 5, 6}, {0, 0.0, 9}};
    v2.emplace_back(x);
    Matrix3x3Upper              v;
    testMatrixCopy(v2[0], v);
    EXPECT_EQ(1, v[XX][XX]);
    EXPECT_EQ(2, v[XX][YY]);
    EXPECT_EQ(3, v[XX][ZZ]);
    EXPECT_EQ(0, v[YY][XX]);
    EXPECT_EQ(5, v[YY][YY]);
    EXPECT_EQ(6, v[YY][ZZ]);
    EXPECT_EQ(0, v[ZZ][XX]);
    EXPECT_EQ(0, v[ZZ][YY]);
    EXPECT_EQ(9, v[ZZ][ZZ]);
}

TEST(Matrix3x3LowerTest, ConvertsImplicitlyFrom_matrix)
{
    std::vector<Matrix3x3Lower> v;
    matrix                      x = {{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};
    v.emplace_back(x);
    EXPECT_EQ(1, v[0][XX][XX]);
    EXPECT_EQ(0, v[0][XX][YY]);
    EXPECT_EQ(0, v[0][XX][ZZ]);
    EXPECT_EQ(4, v[0][YY][XX]);
    EXPECT_EQ(5, v[0][YY][YY]);
    EXPECT_EQ(0, v[0][YY][ZZ]);
    EXPECT_EQ(7, v[0][ZZ][XX]);
    EXPECT_EQ(8, v[0][ZZ][YY]);
    EXPECT_EQ(9, v[0][ZZ][ZZ]);
}

TEST(Matrix3x3LowerTest, UpholdsTheConstraintOnConstruction)
{
    std::vector<Matrix3x3Lower> v;
    matrix                      x = {{0, 0, 3}, {0, 5, 0}, {7, 0.0, 9}};
    EXPECT_DEATH(v.emplace_back(x), "The matrix constraint was violated");
}

TEST(Matrix3x3LowerTest, UpholdsTheConstraintOnAssignment)
{
    std::vector<Matrix3x3Lower> v;
    matrix                      x = {{0, 0, 0}, {0, 5, 0}, {7, 0.0, 9}};
    v.emplace_back(x);
    EXPECT_DEATH(v[0].setValue(XX, YY, 7), "The matrix constraint was violated");
}

TEST(Matrix3x3LowerTest, ConvertsImplicitlyToMatrix3x3)
{
    // This tests that Matrix3x3Lower can be passed as Matrix3x3 &.
    // The opposite (passing Matrix3x3 as Matrix3x3Lower &) will not compile.
    std::vector<Matrix3x3>      v2;
    matrix                      x = {{1, 0, 0}, {4, 5, 0}, {7, 8, 9}};
    v2.emplace_back(x);
    Matrix3x3Lower              v;
    testMatrixCopy(v2[0], v);
    EXPECT_EQ(1, v[XX][XX]);
    EXPECT_EQ(0, v[XX][YY]);
    EXPECT_EQ(0, v[XX][ZZ]);
    EXPECT_EQ(4, v[YY][XX]);
    EXPECT_EQ(5, v[YY][YY]);
    EXPECT_EQ(0, v[YY][ZZ]);
    EXPECT_EQ(7, v[ZZ][XX]);
    EXPECT_EQ(8, v[ZZ][YY]);
    EXPECT_EQ(9, v[ZZ][ZZ]);
}

} // namespace
