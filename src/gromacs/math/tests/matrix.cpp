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
 * Tests matrices
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \author Anatolii Titov <Wapuk-cobaka@yandex.ru>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/matrix.h"

#include <array>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

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
    MatrixTest() :
        matrix_{ { testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1,
                   testNumber_ - 1 } }
    {
    }

protected:
    Matrix3x3                                      matrix_;
    static constexpr real                          testNumber_     = 42;
    static const size_t                            matrixSize_     = DIM * DIM;
    static constexpr std::array<real, matrixSize_> array123456789_ = { 1., 2., 3., 4., 5.,
                                                                       6., 7., 8., 9. };
    static constexpr std::array<real, matrixSize_> array987654321_ = { 9., 8., 7., 6., 5.,
                                                                       4., 3., 2., 1. };
};

TEST_F(MatrixTest, supportsConstexpr)
{
    constexpr Matrix3x3 newMatrix(array123456789_);
    EXPECT_EQ(newMatrix(XX, XX), 1.);
    EXPECT_EQ(newMatrix(XX, YY), 2.);
    EXPECT_EQ(newMatrix(XX, ZZ), 3.);
    EXPECT_EQ(newMatrix(YY, XX), 4.);
    EXPECT_EQ(newMatrix(YY, YY), 5.);
    EXPECT_EQ(newMatrix(YY, ZZ), 6.);
    EXPECT_EQ(newMatrix(ZZ, XX), 7.);
    EXPECT_EQ(newMatrix(ZZ, YY), 8.);
    EXPECT_EQ(newMatrix(ZZ, ZZ), 9.);
}

TEST_F(MatrixTest, legacyAccessRead)
{
    Matrix3x3 newMatrix(array123456789_);
    EXPECT_EQ(newMatrix[XX][XX], 1.);
    EXPECT_EQ(newMatrix[XX][YY], 2.);
    EXPECT_EQ(newMatrix[XX][ZZ], 3.);
    EXPECT_EQ(newMatrix[YY][XX], 4.);
    EXPECT_EQ(newMatrix[YY][YY], 5.);
    EXPECT_EQ(newMatrix[YY][ZZ], 6.);
    EXPECT_EQ(newMatrix[ZZ][XX], 7.);
    EXPECT_EQ(newMatrix[ZZ][YY], 8.);
    EXPECT_EQ(newMatrix[ZZ][ZZ], 9.);
}

TEST_F(MatrixTest, legacyAccessReadConst)
{
    const Matrix3x3 newMatrix(array123456789_);
    EXPECT_EQ(newMatrix[XX][XX], 1.);
    EXPECT_EQ(newMatrix[XX][YY], 2.);
    EXPECT_EQ(newMatrix[XX][ZZ], 3.);
    EXPECT_EQ(newMatrix[YY][XX], 4.);
    EXPECT_EQ(newMatrix[YY][YY], 5.);
    EXPECT_EQ(newMatrix[YY][ZZ], 6.);
    EXPECT_EQ(newMatrix[ZZ][XX], 7.);
    EXPECT_EQ(newMatrix[ZZ][YY], 8.);
    EXPECT_EQ(newMatrix[ZZ][ZZ], 9.);
}

TEST_F(MatrixTest, canGetRVecFromRow)
{
    Matrix3x3 newMatrix(array123456789_);
    RVec      row0 = newMatrix[0];
    EXPECT_EQ(row0[XX], 1.);
    EXPECT_EQ(row0[YY], 2.);
    EXPECT_EQ(row0[ZZ], 3.);
}

TEST_F(MatrixTest, constTests)
{
    std::array<real, matrixSize_> arr = array123456789_;
    Matrix3x3                     newMatrix(arr);
    const Matrix3x3               constNewMatrix1(arr);
    const Matrix3x3               constNewMatrix2(array123456789_);
    const Matrix3x3               constNewMatrix3 = newMatrix;
    const Matrix3x3               constNewMatrix4 = constNewMatrix3;
    EXPECT_EQ(newMatrix, constNewMatrix1);
    EXPECT_EQ(newMatrix, constNewMatrix2);
    EXPECT_EQ(newMatrix, constNewMatrix3);
    EXPECT_EQ(constNewMatrix3, constNewMatrix4);
    RVec       row0_1      = newMatrix[XX];
    RVec       row0_2      = constNewMatrix1[XX];
    const RVec constRow0_1 = constNewMatrix1[XX];
    const RVec constRow0_2 = newMatrix[XX];
    EXPECT_EQ(row0_1, row0_2);
    EXPECT_EQ(row0_1, constRow0_1);
    EXPECT_EQ(row0_1, constRow0_2);
    EXPECT_EQ(constRow0_1[XX], 1.);
    EXPECT_EQ(constRow0_1[YY], 2.);
    EXPECT_EQ(constRow0_1[ZZ], 3.);
}

TEST_F(MatrixTest, canSetFromArray)
{
    Matrix3x3 newMatrix(array123456789_);
    EXPECT_EQ(newMatrix(XX, XX), 1.);
    EXPECT_EQ(newMatrix(XX, YY), 2.);
    EXPECT_EQ(newMatrix(XX, ZZ), 3.);
    EXPECT_EQ(newMatrix(YY, XX), 4.);
    EXPECT_EQ(newMatrix(YY, YY), 5.);
    EXPECT_EQ(newMatrix(YY, ZZ), 6.);
    EXPECT_EQ(newMatrix(ZZ, XX), 7.);
    EXPECT_EQ(newMatrix(ZZ, YY), 8.);
    EXPECT_EQ(newMatrix(ZZ, ZZ), 9.);
}

TEST_F(MatrixTest, canSetStaticallyFromList)
{
    Matrix3x3 newMatrix = { { 1., 2., 3., 4., 5., 6., 7., 8., 9. } };
    EXPECT_EQ(newMatrix(XX, XX), 1.);
    EXPECT_EQ(newMatrix(XX, YY), 2.);
    EXPECT_EQ(newMatrix(XX, ZZ), 3.);
    EXPECT_EQ(newMatrix(YY, XX), 4.);
    EXPECT_EQ(newMatrix(YY, YY), 5.);
    EXPECT_EQ(newMatrix(YY, ZZ), 6.);
    EXPECT_EQ(newMatrix(ZZ, XX), 7.);
    EXPECT_EQ(newMatrix(ZZ, YY), 8.);
    EXPECT_EQ(newMatrix(ZZ, ZZ), 9.);
}

TEST_F(MatrixTest, canConstructAndFill)
{
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            EXPECT_EQ(testNumber_ - 1, matrix_[i][j]);
        }
    }
}

TEST_F(MatrixTest, canSetValues)
{
    matrix_(YY, YY) = testNumber_;
    EXPECT_EQ(testNumber_, matrix_(YY, YY));
}

TEST_F(MatrixTest, legacySetValues)
{
    matrix_[YY][YY] = testNumber_;
    EXPECT_EQ(testNumber_, matrix_(YY, YY));
}

TEST_F(MatrixTest, canCopyAssign)
{
    Matrix3x3 other;
    other = matrix_;
    using ::testing::Eq;
    using ::testing::Pointwise;
    EXPECT_THAT(other.toArrayRef(), Pointwise(Eq(), matrix_.toArrayRef()));
}

TEST_F(MatrixTest, canSwap)
{
    Matrix3x3 other;
    matrix_(XX, XX) = testNumber_;
    other.swap(matrix_);
    EXPECT_EQ(testNumber_, other[XX][XX]);
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            if (i == XX && j == XX)
            {
                continue;
            }
            EXPECT_EQ(testNumber_ - 1., other[i][j]);
        }
    }
}

TEST_F(MatrixTest, canClear)
{
    Matrix3x3 newMatrix = array123456789_;
    EXPECT_EQ(newMatrix(XX, XX), 1.);
    EXPECT_EQ(newMatrix(XX, YY), 2.);
    EXPECT_EQ(newMatrix(XX, ZZ), 3.);
    EXPECT_EQ(newMatrix(YY, XX), 4.);
    EXPECT_EQ(newMatrix(YY, YY), 5.);
    EXPECT_EQ(newMatrix(YY, ZZ), 6.);
    EXPECT_EQ(newMatrix(ZZ, XX), 7.);
    EXPECT_EQ(newMatrix(ZZ, YY), 8.);
    EXPECT_EQ(newMatrix(ZZ, ZZ), 9.);
    newMatrix.clear();
    EXPECT_EQ(newMatrix(XX, XX), 0);
    EXPECT_EQ(newMatrix(XX, YY), 0);
    EXPECT_EQ(newMatrix(XX, ZZ), 0);
    EXPECT_EQ(newMatrix(YY, XX), 0);
    EXPECT_EQ(newMatrix(YY, YY), 0);
    EXPECT_EQ(newMatrix(YY, ZZ), 0);
    EXPECT_EQ(newMatrix(ZZ, XX), 0);
    EXPECT_EQ(newMatrix(ZZ, YY), 0);
    EXPECT_EQ(newMatrix(ZZ, ZZ), 0);
}


TEST_F(MatrixTest, canAddMatrix)
{
    const Matrix3x3 matrixA = array123456789_;
    const Matrix3x3 matrixB = array987654321_;
    Matrix3x3       sum;
    sum = matrixA + matrixB;
    EXPECT_EQ(sum(XX, XX), 10.);
    EXPECT_EQ(sum(XX, YY), 10.);
    EXPECT_EQ(sum(XX, ZZ), 10.);
    EXPECT_EQ(sum(YY, XX), 10.);
    EXPECT_EQ(sum(YY, YY), 10.);
    EXPECT_EQ(sum(YY, ZZ), 10.);
    EXPECT_EQ(sum(ZZ, XX), 10.);
    EXPECT_EQ(sum(ZZ, YY), 10.);
    EXPECT_EQ(sum(ZZ, ZZ), 10.);
}

TEST_F(MatrixTest, canAddAssignMatrix)
{
    Matrix3x3       matrixA = array123456789_;
    const Matrix3x3 matrixB = array987654321_;
    matrixA += matrixB;
    EXPECT_EQ(matrixA(XX, XX), 10.);
    EXPECT_EQ(matrixA(XX, YY), 10.);
    EXPECT_EQ(matrixA(XX, ZZ), 10.);
    EXPECT_EQ(matrixA(YY, XX), 10.);
    EXPECT_EQ(matrixA(YY, YY), 10.);
    EXPECT_EQ(matrixA(YY, ZZ), 10.);
    EXPECT_EQ(matrixA(ZZ, XX), 10.);
    EXPECT_EQ(matrixA(ZZ, YY), 10.);
    EXPECT_EQ(matrixA(ZZ, ZZ), 10.);
}

TEST_F(MatrixTest, canSubtractMatrix)
{
    const Matrix3x3 matrixA = array123456789_;
    const Matrix3x3 matrixB = array987654321_;
    Matrix3x3       difference;
    difference = matrixA - matrixB;
    EXPECT_EQ(difference(XX, XX), -8.);
    EXPECT_EQ(difference(XX, YY), -6.);
    EXPECT_EQ(difference(XX, ZZ), -4.);
    EXPECT_EQ(difference(YY, XX), -2.);
    EXPECT_EQ(difference(YY, YY), 0.);
    EXPECT_EQ(difference(YY, ZZ), 2.);
    EXPECT_EQ(difference(ZZ, XX), 4.);
    EXPECT_EQ(difference(ZZ, YY), 6.);
    EXPECT_EQ(difference(ZZ, ZZ), 8.);
}

TEST_F(MatrixTest, canSubtractAssignMatrix)
{
    Matrix3x3       matrixA = array123456789_;
    const Matrix3x3 matrixB = array987654321_;
    matrixA -= matrixB;
    EXPECT_EQ(matrixA(XX, XX), -8.);
    EXPECT_EQ(matrixA(XX, YY), -6.);
    EXPECT_EQ(matrixA(XX, ZZ), -4.);
    EXPECT_EQ(matrixA(YY, XX), -2.);
    EXPECT_EQ(matrixA(YY, YY), 0.);
    EXPECT_EQ(matrixA(YY, ZZ), 2.);
    EXPECT_EQ(matrixA(ZZ, XX), 4.);
    EXPECT_EQ(matrixA(ZZ, YY), 6.);
    EXPECT_EQ(matrixA(ZZ, ZZ), 8.);
}

TEST_F(MatrixTest, canNegateMatrix)
{
    Matrix3x3 matrixA = array123456789_;
    Matrix3x3 matrixB;
    matrixB = -matrixA;
    EXPECT_EQ(matrixB(XX, XX), -1.);
    EXPECT_EQ(matrixB(XX, YY), -2.);
    EXPECT_EQ(matrixB(XX, ZZ), -3.);
    EXPECT_EQ(matrixB(YY, XX), -4.);
    EXPECT_EQ(matrixB(YY, YY), -5.);
    EXPECT_EQ(matrixB(YY, ZZ), -6.);
    EXPECT_EQ(matrixB(ZZ, XX), -7.);
    EXPECT_EQ(matrixB(ZZ, YY), -8.);
    EXPECT_EQ(matrixB(ZZ, ZZ), -9.);
}

TEST_F(MatrixTest, MatrixScalarMultiplication)
{
    const Matrix3x3 matrix = { { 1., 2., 3., 3., 2., 1., 1., 1., 1. } };
    real            scalar = 2.;
    Matrix3x3       scaledMatrix;
    scaledMatrix = scalar * matrix;
    EXPECT_EQ(scaledMatrix(XX, XX), 2.);
    EXPECT_EQ(scaledMatrix(XX, YY), 4.);
    EXPECT_EQ(scaledMatrix(XX, ZZ), 6.);
    EXPECT_EQ(scaledMatrix(YY, XX), 6.);
    EXPECT_EQ(scaledMatrix(YY, YY), 4.);
    EXPECT_EQ(scaledMatrix(YY, ZZ), 2.);
    EXPECT_EQ(scaledMatrix(ZZ, XX), 2.);
    EXPECT_EQ(scaledMatrix(ZZ, YY), 2.);
    EXPECT_EQ(scaledMatrix(ZZ, ZZ), 2.);
}

TEST_F(MatrixTest, MatrixScalarDivision)
{
    const Matrix3x3 matrix = { { 1., 2., 3., 3., 2., 1., 1., 1., 1. } };
    real            scalar = 2.;
    Matrix3x3       scaledMatrix;
    scaledMatrix = matrix / scalar;
    EXPECT_EQ(scaledMatrix(XX, XX), 0.5);
    EXPECT_EQ(scaledMatrix(XX, YY), 1.);
    EXPECT_EQ(scaledMatrix(XX, ZZ), 1.5);
    EXPECT_EQ(scaledMatrix(YY, XX), 1.5);
    EXPECT_EQ(scaledMatrix(YY, YY), 1.);
    EXPECT_EQ(scaledMatrix(YY, ZZ), 0.5);
    EXPECT_EQ(scaledMatrix(ZZ, XX), 0.5);
    EXPECT_EQ(scaledMatrix(ZZ, YY), 0.5);
    EXPECT_EQ(scaledMatrix(ZZ, ZZ), 0.5);
}

TEST_F(MatrixTest, MatrixVectorMultiplicationOperator)
{
    const Matrix3x3 matrix = { { 0.1, 1., 0.1, 0.4, 1., 0.6, 0.7, 0.8, 0.9 } };
    const RVec      vector(1., 2., 3.);
    RVec            result = matrix * vector;
    EXPECT_REAL_EQ(2.4, result[XX]);
    EXPECT_REAL_EQ(4.2, result[YY]);
    EXPECT_REAL_EQ(5.0, result[ZZ]);
}

TEST_F(MatrixTest, MatrixMatrixInnerProduct)
{
    const Matrix3x3 matrixA(array123456789_);
    const Matrix3x3 matrixB(array987654321_);
    Matrix3x3       matrixC = inner(matrixA, matrixB);
    EXPECT_REAL_EQ(matrixC(XX, XX), 30.);
    EXPECT_REAL_EQ(matrixC(XX, YY), 24.);
    EXPECT_REAL_EQ(matrixC(XX, ZZ), 18.);
    EXPECT_REAL_EQ(matrixC(YY, XX), 84.);
    EXPECT_REAL_EQ(matrixC(YY, YY), 69.);
    EXPECT_REAL_EQ(matrixC(YY, ZZ), 54.);
    EXPECT_REAL_EQ(matrixC(ZZ, XX), 138.);
    EXPECT_REAL_EQ(matrixC(ZZ, YY), 114.);
    EXPECT_REAL_EQ(matrixC(ZZ, ZZ), 90.);
}

TEST_F(MatrixTest, MatrixMatrixMultiplication)
{
    const Matrix3x3 matrixA(array123456789_);
    const Matrix3x3 matrixB(array987654321_);
    Matrix3x3       matrixC = matrixA * matrixB;
    EXPECT_REAL_EQ(matrixC(XX, XX), 30.);
    EXPECT_REAL_EQ(matrixC(XX, YY), 24.);
    EXPECT_REAL_EQ(matrixC(XX, ZZ), 18.);
    EXPECT_REAL_EQ(matrixC(YY, XX), 84.);
    EXPECT_REAL_EQ(matrixC(YY, YY), 69.);
    EXPECT_REAL_EQ(matrixC(YY, ZZ), 54.);
    EXPECT_REAL_EQ(matrixC(ZZ, XX), 138.);
    EXPECT_REAL_EQ(matrixC(ZZ, YY), 114.);
    EXPECT_REAL_EQ(matrixC(ZZ, ZZ), 90.);
}

TEST_F(MatrixTest, determinantWorks)
{
    const Matrix3x3 mat = { { 1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0 } };
    EXPECT_EQ(determinant(mat), 1.);
}

TEST_F(MatrixTest, noninvertableDeterminantIsZero)
{
    const Matrix3x3 mat = { { 1., 0., 0., 0., 1., 0., 0., 0., 0. } };
    EXPECT_EQ(determinant(mat), 0.);
}

TEST_F(MatrixTest, determinantOfDiagonalMatrix)
{
    const Matrix3x3 mat = { { 2., 0., 0., 0., 3., 0., 0., 0., 4. } };
    EXPECT_EQ(determinant(mat), 24.);
}

TEST_F(MatrixTest, traceWorks)
{
    const Matrix3x3 mat = { { 1.5, 9., 9., 9., 2.0, 9., 9., 9., 0.25 } };
    EXPECT_EQ(trace(mat), 3.75);
}

TEST_F(MatrixTest, transposeWorks)
{
    const Matrix3x3 asymmetricMat(array123456789_);

    const Matrix3x3 transposedAsymmetricMat = transpose(asymmetricMat);
    EXPECT_EQ(asymmetricMat(XX, XX), transposedAsymmetricMat(XX, XX));
    EXPECT_EQ(asymmetricMat(XX, YY), transposedAsymmetricMat(YY, XX));
    EXPECT_EQ(asymmetricMat(XX, ZZ), transposedAsymmetricMat(ZZ, XX));
    EXPECT_EQ(asymmetricMat(YY, XX), transposedAsymmetricMat(XX, YY));
    EXPECT_EQ(asymmetricMat(YY, YY), transposedAsymmetricMat(YY, YY));
    EXPECT_EQ(asymmetricMat(YY, ZZ), transposedAsymmetricMat(ZZ, YY));
    EXPECT_EQ(asymmetricMat(ZZ, XX), transposedAsymmetricMat(XX, ZZ));
    EXPECT_EQ(asymmetricMat(ZZ, YY), transposedAsymmetricMat(YY, ZZ));
    EXPECT_EQ(asymmetricMat(ZZ, ZZ), transposedAsymmetricMat(ZZ, ZZ));
}

TEST_F(MatrixTest, transposeOfSymmetricMatrix)
{
    const Matrix3x3 symmetricMat           = { { 1., 2., 3., 2., 5., 6., 3., 6., 9. } };
    const Matrix3x3 transposedSymmetricMat = transpose(symmetricMat);
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            EXPECT_EQ(symmetricMat(i, j), transposedSymmetricMat(i, j));
        }
    }
}

TEST_F(MatrixTest, canCreateFromLegacyMatrix)
{
    matrix          legacyMatrix = { { 1., 2., 3. }, { 4., 5., 6. }, { 7., 8., 9. } };
    const Matrix3x3 fromLegacy   = createMatrix3x3FromLegacyMatrix(legacyMatrix);
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            EXPECT_EQ(fromLegacy(i, j), legacyMatrix[i][j]);
        }
    }
}

TEST_F(MatrixTest, canFillLegacyMatrix)
{
    matrix legacyMatrix = { { -2. } };
    fillLegacyMatrix(matrix_, legacyMatrix);
    for (size_t i = XX; i < DIM; ++i)
    {
        for (size_t j = XX; j < DIM; ++j)
        {
            EXPECT_EQ(legacyMatrix[i][j], matrix_(i, j));
        }
    }
}

TEST_F(MatrixTest, IdentityMatrix)
{
    const Matrix3x3 realIdMatrix = identityMatrix<real>();
    EXPECT_REAL_EQ(realIdMatrix(XX, XX), 1.);
    EXPECT_REAL_EQ(realIdMatrix(XX, YY), 0.);
    EXPECT_REAL_EQ(realIdMatrix(XX, ZZ), 0.);
    EXPECT_REAL_EQ(realIdMatrix(YY, XX), 0.);
    EXPECT_REAL_EQ(realIdMatrix(YY, YY), 1.);
    EXPECT_REAL_EQ(realIdMatrix(YY, ZZ), 0.);
    EXPECT_REAL_EQ(realIdMatrix(ZZ, XX), 0.);
    EXPECT_REAL_EQ(realIdMatrix(ZZ, YY), 0.);
    EXPECT_REAL_EQ(realIdMatrix(ZZ, ZZ), 1.);
}

TEST_F(MatrixTest, DiagonalMatrix)
{
    const Matrix3x3 realDMatrix = diagonalMatrix<real>(10);
    EXPECT_REAL_EQ(realDMatrix(XX, XX), 10.);
    EXPECT_REAL_EQ(realDMatrix(XX, YY), 0.);
    EXPECT_REAL_EQ(realDMatrix(XX, ZZ), 0.);
    EXPECT_REAL_EQ(realDMatrix(YY, XX), 0.);
    EXPECT_REAL_EQ(realDMatrix(YY, YY), 10.);
    EXPECT_REAL_EQ(realDMatrix(YY, ZZ), 0.);
    EXPECT_REAL_EQ(realDMatrix(ZZ, XX), 0.);
    EXPECT_REAL_EQ(realDMatrix(ZZ, YY), 0.);
    EXPECT_REAL_EQ(realDMatrix(ZZ, ZZ), 10.);
}

} // namespace

} // namespace test

} // namespace gmx
