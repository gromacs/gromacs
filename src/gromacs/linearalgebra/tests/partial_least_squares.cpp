/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests partial least squares implementation.
 *
 * Both the asssociated data types (FMatrix and FVector) and
 * the function pls_denham() are tested
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 */
#include "gmxpre.h"

#include "gromacs/linearalgebra/partial_least_squares.h"

#include <vector>

#include <gtest/gtest.h>

#include "testutils/testasserts.h"

namespace
{

TEST(PLSTest, EmptyFVectorIsEmpty)
{
    FVector v = FVector();
    EXPECT_EQ(0, v.Length());
}

TEST(PLSTest, FVectorConvertsToFORTRAN)
{
    FVector v = FVector(3);
    v[0] = 1;
    v[1] = 2;
    v[2] = 3;
    real *fv = v.toF();
    EXPECT_EQ(1.0, fv[0]);
    EXPECT_EQ(2.0, fv[1]);
    EXPECT_EQ(3.0, fv[2]);
}

TEST(PLSTest, EmptyFMatrixIsEmpty)
{
    FMatrix m = FMatrix();
    EXPECT_EQ(0, m.NRows());
    EXPECT_EQ(0, m.NCols());
}

TEST(PLSTest, FMatrixConvertsToFORTRAN)
{
    int     i, j;

    FMatrix m = FMatrix(3, 3);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            m(i, j) = 10*(i+1) + (j+1);
        }
    }
    real *fm = m.toF();
    EXPECT_EQ(11.0, fm[0]);
    EXPECT_EQ(21.0, fm[1]);
    EXPECT_EQ(31.0, fm[2]);
    EXPECT_EQ(12.0, fm[3]);
    EXPECT_EQ(22.0, fm[4]);
    EXPECT_EQ(32.0, fm[5]);
    EXPECT_EQ(13.0, fm[6]);
    EXPECT_EQ(23.0, fm[7]);
    EXPECT_EQ(33.0, fm[8]);
}

TEST(PLSTest, PLS_Denham)
{
    int     i, k;
    int     nFrames_  = 5;
    int     nAtoms_   = 5;
    int     nDim      = 2;

    FMatrix fframes_ = FMatrix(nFrames_, nAtoms_);
    FVector avg_     = FVector(nAtoms_);
    FMatrix w_       =  FMatrix(nAtoms_, nDim);
    FVector q_       = FVector(nDim);
    FVector yTrain_  = FVector(nFrames_);
    FVector PLSvec(nAtoms_);

    // Matrix elements randomly generated once (and centered)
    fframes_(0, 0) = -2.7;
    fframes_(0, 1) = 3.0;
    fframes_(0, 2) = 2.8;
    fframes_(0, 3) = -1.6;
    fframes_(0, 4) = -1.7;
    fframes_(1, 0) = -1.0;
    fframes_(1, 1) = 1.9;
    fframes_(1, 2) = -3.1;
    fframes_(1, 3) = -2.8;
    fframes_(1, 4) = -1.0;
    fframes_(2, 0) = -0.2;
    fframes_(2, 1) = -2.7;
    fframes_(2, 2) = -3.2;
    fframes_(2, 3) = 3.3;
    fframes_(2, 4) = 3.3;
    fframes_(3, 0) = 1.6;
    fframes_(3, 1) = 1.8;
    fframes_(3, 2) = 0.5;
    fframes_(3, 3) = -1.7;
    fframes_(3, 4) = 0.5;
    fframes_(4, 0) = 2.4;
    fframes_(4, 1) = -4.2;
    fframes_(4, 2) = 3.0;
    fframes_(4, 3) = 2.7;
    fframes_(4, 4) = -1.3;

    // Vector elements randomly generated once (and centered)
    yTrain_[0] = 0.8;
    yTrain_[1] = -2.6;
    yTrain_[2] = 3.7;
    yTrain_[3] = 0.1;
    yTrain_[4] = -2.1;

    pls_denham(&fframes_, &yTrain_, nFrames_, nDim, &w_, &q_);

    for (i = 0; i < nAtoms_; i++)
    {
        for (k = 0; k < nDim; k++)
        {
            PLSvec[i] += q_[k] * w_(i, k);
        }
    }

    // these values have been calculated with a working version
    // of pls_denham. As the results have been given only with
    // single digit precision, the tolerance is rather high
    EXPECT_REAL_EQ_TOL(235.6, PLSvec[0], gmx::test::absoluteTolerance(0.1));
    EXPECT_REAL_EQ_TOL(-1287.2, PLSvec[1], gmx::test::absoluteTolerance(0.1));
    EXPECT_REAL_EQ_TOL(-1138.4, PLSvec[2], gmx::test::absoluteTolerance(0.1));
    EXPECT_REAL_EQ_TOL(1315.1, PLSvec[3], gmx::test::absoluteTolerance(0.1));
    EXPECT_REAL_EQ_TOL(1144.4, PLSvec[4], gmx::test::absoluteTolerance(0.1));
}


} // namespace
