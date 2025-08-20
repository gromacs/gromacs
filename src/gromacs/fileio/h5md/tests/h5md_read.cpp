/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Tests for reading data from H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_read.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_dataset.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for an open H5md file
template<typename ValueType>
class H5mdReadNumericPrimitiveTest : public H5mdTestBase
{
};

//! \brief Test fixture for an open H5md file
template<typename ValueType>
class H5mdReadBasicVectorListTest : public H5mdTestBase
{
};

//! \brief List of numeric primitives for which to create tests
using NumericPrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;
TYPED_TEST_SUITE(H5mdReadNumericPrimitiveTest, NumericPrimitiveList);

//! \brief List of numeric primitives for which to create real-valued tests
using RealPrimitiveList = ::testing::Types<float, double>;
TYPED_TEST_SUITE(H5mdReadBasicVectorListTest, RealPrimitiveList);

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueFrom1dSetWorks)
{
    const std::vector<TypeParam>     values = { 9, 3, 1, 0, 6 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(values.size()) // enable writing the values above into the data set as frames
                    .build();

    // Write values to data set in bulk for efficiency
    H5Dwrite(dataSet.id(), dataSet.dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

    TypeParam value;
    for (hsize_t frameIndex = 0; frameIndex < values.size(); ++frameIndex)
    {
        readFrame(dataSet, frameIndex, value);
        EXPECT_FLOAT_EQ(value, values[frameIndex]);
    }
}

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueOutsideOfSetBoundsThrows)
{
    const std::vector<TypeParam>     values = { 9, 3, 1, 0, 6 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(values.size()) // enable writing the values above into the data set as frames
                    .build();

    // Write values to data set in bulk for efficiency
    H5Dwrite(dataSet.id(), dataSet.dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

    TypeParam value;
    EXPECT_NO_THROW(readFrame(dataSet, values.size() - 1, value))
            << "Sanity check failed: reading last value must work";
    EXPECT_THROW(readFrame(dataSet, values.size(), value), gmx::FileIOError)
            << "Must throw when reading out of bounds";
}

TYPED_TEST(H5mdReadBasicVectorListTest, ReadFrameWorks)
{
    constexpr int                                      numAtoms = 5;
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build();

    // Write 5 RVecs to the data set (this matches its complete size: 1 frame * 5 atoms * 3 values)
    const std::array<gmx::BasicVector<TypeParam>, numAtoms> rvecListToWrite = {
        BasicVector<TypeParam>{ 0.0, 1.0, 2.0 },
        BasicVector<TypeParam>{ 3.0, 4.0, 5.0 },
        BasicVector<TypeParam>{ 6.0, 7.0, 8.0 },
        BasicVector<TypeParam>{ 9.0, 10.0, 11.0 },
        BasicVector<TypeParam>{ 12.0, 13.0, 14.0 }
    };
    H5Dwrite(dataSet.id(), dataSet.dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, rvecListToWrite.data());

    // Read back the frame data (5 atoms * 3 values at index = 0) and compare to the original
    std::array<gmx::BasicVector<TypeParam>, numAtoms> readBuffer;
    readBuffer.fill(BasicVector<TypeParam>{ -1.0, -1.0, -1.0 }); // initialize buffer to sentinel values
    readFrame(dataSet, 0, makeArrayRef(readBuffer));

    EXPECT_THAT(readBuffer, ::testing::Pointwise(::testing::Eq(), rvecListToWrite));
}

TYPED_TEST(H5mdReadBasicVectorListTest, ReadingFrameFromNonSequentialIndexWorks)
{
    constexpr int                                      numAtoms  = 5;
    constexpr hsize_t                                  numFrames = 5;
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(numFrames)
                    .build();

    // Prepare a list of numFrames rvec-lists and write it to the data set in bulk
    // Total number of values matches that of the data set: 5 frames * 5 atoms * 3 values
    std::array<std::array<gmx::BasicVector<TypeParam>, numAtoms>, numFrames> perFrameRvecLists;
    for (int frameIndex = 0; frameIndex < gmx::ssize(perFrameRvecLists); ++frameIndex)
    {
        const TypeParam frameValue = 1000 * frameIndex; // assign unique per-frame values to list;
        for (int atomIndex = 0; atomIndex < gmx::ssize(perFrameRvecLists[frameIndex]); ++atomIndex)
        {
            perFrameRvecLists[frameIndex][atomIndex][XX] = frameValue + 100 * atomIndex;
            perFrameRvecLists[frameIndex][atomIndex][YY] = frameValue + 10 * atomIndex;
            perFrameRvecLists[frameIndex][atomIndex][ZZ] = frameValue + atomIndex;
        }
    }
    H5Dwrite(dataSet.id(), dataSet.dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, perFrameRvecLists.data());


    // Check that reading all 5 frames in non-sequential order works fine
    for (const int i : { 3, 1, 4, 2, 0 })
    {
        std::array<gmx::BasicVector<TypeParam>, numAtoms> readBuffer;
        readBuffer.fill(BasicVector<TypeParam>{ -1.0, -1.0, -1.0 });
        readFrame(dataSet, i, makeArrayRef(readBuffer));

        EXPECT_THAT(readBuffer, ::testing::Pointwise(::testing::Eq(), perFrameRvecLists[i]));
    }
}

TYPED_TEST(H5mdReadBasicVectorListTest, OutputBufferWithNonMatchingDimensionsThrows)
{
    constexpr int                                      numAtoms = 3;
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build();

    std::array<BasicVector<TypeParam>, numAtoms> readBufferCorrect;
    ASSERT_NO_THROW(readFrame(dataSet, 0, makeArrayRef(readBufferCorrect)))
            << "Sanity check failed: reading with a correctly sized buffer must work";

    std::array<BasicVector<TypeParam>, numAtoms - 1> smallReadBuffer;
    EXPECT_THROW(readFrame(dataSet, 0, makeArrayRef(smallReadBuffer)), gmx::FileIOError);
    std::array<BasicVector<TypeParam>, numAtoms + 1> largeReadBuffer;
    EXPECT_THROW(readFrame(dataSet, 0, makeArrayRef(largeReadBuffer)), gmx::FileIOError);
}

TYPED_TEST(H5mdReadBasicVectorListTest, CorrectRowOrderIsUsedForReadingArrayData)
{
    constexpr int                                      numAtoms = 5;
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build();

    // Write a 1d array with values 0..(numAtoms * DIM) to the data set
    std::array<TypeParam, DIM * numAtoms> writeBuffer1d;
    for (int i = 0; i < DIM * numAtoms; ++i)
    {
        writeBuffer1d[i] = static_cast<TypeParam>(i);
    }
    H5Dwrite(dataSet.id(), dataSet.dataType(), H5S_ALL, H5S_ALL, H5P_DEFAULT, writeBuffer1d.data());

    // Ensure that we are consistent when reading it back as a 2d [numAtoms, DIM] array
    std::array<BasicVector<TypeParam>, numAtoms> readBuffer;
    readFrame(dataSet, 0, makeArrayRef(readBuffer));

    int index1d = 0;
    for (int i = 0; i < numAtoms; ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            EXPECT_REAL_EQ(readBuffer[i][d], writeBuffer1d[index1d]);
            index1d++;
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
