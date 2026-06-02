/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * Tests utilities for generating trajectory frame data.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/generate_frame_data.h"

#include <array>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/matrix.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"


namespace gmx
{
namespace test
{
namespace
{

template<typename ValueType>
class TrajectoryFrameGenerationTest : public ::testing::Test
{
};

using VecTypes = ::testing::Types<float[3], double[3], BasicVector<float>, BasicVector<double>>;
TYPED_TEST_SUITE(TrajectoryFrameGenerationTest, VecTypes);

TYPED_TEST(TrajectoryFrameGenerationTest, SequentialPositive)
{
    constexpr int numAtoms  = 3;
    constexpr int numFrames = 5;

    constexpr double startOffset = 15.0;
    constexpr double frameOffset = 3.0;
    constexpr double valueOffset = 2.0;
    constexpr double dimOffset   = 0.1;

    TrajectoryFrameDataGenerator frameDataGenerator(
            TrajectoryFrameMode::Positive, { startOffset, frameOffset, valueOffset, dimOffset });

    for (int frame = 0; frame < numFrames; ++frame)
    {
        std::array<TypeParam, numAtoms> frameValues;
        frameDataGenerator(frameValues);

        for (int i = 0; i < numAtoms; ++i)
        {
            for (int d = 0; d < DIM; ++d)
            {
                EXPECT_FLOAT_EQ(frameValues[i][d],
                                startOffset + (frame * frameOffset) + (i * valueOffset)
                                        + ((d + 1) * dimOffset));
            }
        }
    }
}

TYPED_TEST(TrajectoryFrameGenerationTest, SequentialNegative)
{
    constexpr int numAtoms  = 3;
    constexpr int numFrames = 5;

    constexpr double startOffset = 15.0;
    constexpr double frameOffset = 3.0;
    constexpr double valueOffset = 2.0;
    constexpr double dimOffset   = 0.1;

    TrajectoryFrameDataGenerator frameDataGenerator(
            TrajectoryFrameMode::Negative, { startOffset, frameOffset, valueOffset, dimOffset });

    for (int frame = 0; frame < numFrames; ++frame)
    {
        std::array<TypeParam, numAtoms> frameValues;
        frameDataGenerator(frameValues);

        for (int i = 0; i < numAtoms; ++i)
        {
            for (int d = 0; d < DIM; ++d)
            {
                EXPECT_FLOAT_EQ(frameValues[i][d],
                                -(startOffset + (frame * frameOffset) + (i * valueOffset)
                                  + ((d + 1) * dimOffset)));
            }
        }
    }
}

TYPED_TEST(TrajectoryFrameGenerationTest, SequentialAlternating)
{
    constexpr int numAtoms  = 3;
    constexpr int numFrames = 5;

    constexpr double startOffset = 15.0;
    constexpr double frameOffset = 3.0;
    constexpr double valueOffset = 2.0;
    constexpr double dimOffset   = 0.1;

    TrajectoryFrameDataGenerator frameDataGenerator(
            TrajectoryFrameMode::Alternating, { startOffset, frameOffset, valueOffset, dimOffset });

    for (int frame = 0; frame < numFrames; ++frame)
    {
        std::array<TypeParam, numAtoms> frameValues;
        frameDataGenerator(frameValues);
        double sign = -1.0;

        for (int i = 0; i < numAtoms; ++i)
        {
            for (int d = 0; d < DIM; ++d)
            {
                EXPECT_FLOAT_EQ(frameValues[i][d],
                                sign
                                        * (startOffset + (frame * frameOffset) + (i * valueOffset)
                                           + ((d + 1) * dimOffset)));
            }
            sign *= -1.0;
        }
    }
}

TEST(TrajectoryFrameGenerationTest, DefaultGenerator)
{
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1.0, 1e-6));

    constexpr int numAtoms      = 3;
    constexpr int numIterations = 3;

    std::array<RVec, numAtoms>   values;
    TrajectoryFrameDataGenerator defaultGenerator;
    for (int i = 0; i < numIterations; ++i)
    {
        defaultGenerator(values);
        checker.checkSequence(makeConstArrayRef(values), formatString("frame %d", i).c_str());
    }
}

template<typename ValueType>
class MatrixFrameGenerationTest : public ::testing::Test
{
};

using MatrixTypes =
        ::testing::Types<float[3][3], double[3][3], BasicMatrix3x3<float>, BasicMatrix3x3<double>>;
TYPED_TEST_SUITE(MatrixFrameGenerationTest, MatrixTypes);

TYPED_TEST(MatrixFrameGenerationTest, Sequential)
{
    constexpr int numFrames = 5;

    constexpr double startOffset = 15.0;
    constexpr double frameOffset = 3.0;
    constexpr double valueOffset = 2.0;

    MatrixFrameDataGenerator matrixDataGenerator({ startOffset, frameOffset, valueOffset });

    for (int frame = 0; frame < numFrames; ++frame)
    {
        TypeParam box;
        matrixDataGenerator(box);

        for (int i = 0; i < DIM; ++i)
        {
            for (int j = 0; j < DIM; ++j)
            {
                EXPECT_FLOAT_EQ(box[i][j],
                                startOffset + (frame * frameOffset) + (valueOffset * ((i * DIM) + j + 1)));
            }
        }
    }
}

TEST(MatrixFrameGenerationTest, DefaultGenerator)
{
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(relativeToleranceAsFloatingPoint(1.0, 1e-6));

    constexpr int numIterations = 3;

    Matrix3x3                box;
    MatrixFrameDataGenerator defaultGenerator;
    for (int i = 0; i < numIterations; ++i)
    {
        defaultGenerator(box);
        checker.checkSequence(makeConstArrayRef(box.toConstArrayRef()),
                              formatString("frame %d", i).c_str());
    }
}

} // namespace
} // namespace test
} // namespace gmx
