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
 * \brief Tests for writing data to H5MD (HDF5) files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_write.h"

#include <hdf5.h>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/h5md_dataset.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_read.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for an open H5md file
template<typename ValueType>
class H5mdWriteNumericPrimitiveTest : public H5mdTestBase
{
};

//! \brief Test fixture for an open H5md file
template<typename ValueType>
class H5mdWriteBasicVectorListTest : public H5mdTestBase
{
};

//! \brief List of numeric primitives for which to create tests
using NumericPrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;
TYPED_TEST_SUITE(H5mdWriteNumericPrimitiveTest, NumericPrimitiveList);

//! \brief List of numeric primitives for which to create real-valued tests
using RealPrimitiveList = ::testing::Types<float, double>;
TYPED_TEST_SUITE(H5mdWriteBasicVectorListTest, RealPrimitiveList);

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueTo1dSetWritesToCorrectIndex)
{
    // List of values to write into data set one-by-one with our tested function
    const std::vector<TypeParam> valuesToWrite = { 9, 3, 1, 0, 6 };
    const auto [dataSet, dataSetGuard]         = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(valuesToWrite.size()) // enable reading and writing with index = 0
                    .build());

    for (hsize_t writeIndex = 0; writeIndex < valuesToWrite.size(); ++writeIndex)
    {
        writeFrame(dataSet, writeIndex, valuesToWrite[writeIndex]);
    }

    // Read back the values in bulk, then compare to those which were written
    std::vector<TypeParam> readBuffer(valuesToWrite.size());
    H5Dread(dataSet, hdf5DataTypeFor<TypeParam>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, readBuffer.data());
    EXPECT_THAT(readBuffer, ::testing::Pointwise(::testing::Eq(), valuesToWrite));
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueOutsideOfSetBoundsThrows)
{
    constexpr hsize_t numFrames = 5;
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                                         .withNumFrames(numFrames)
                                         .build());

    constexpr TypeParam value = 0;
    ASSERT_NO_THROW(writeFrame(dataSet, numFrames - 1, value))
            << "Sanity check failed: write last value must work";
    EXPECT_THROW(writeFrame(dataSet, numFrames, value), gmx::FileIOError)
            << "Must throw when writing out of bounds";
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueOfNonMatchingTypeThrows)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                                         .withNumFrames(1) // enable reading and writing with index = 0
                                         .build());

    if (!std::is_same_v<TypeParam, float>)
    {
        constexpr float value = 0;
        EXPECT_THROW(writeFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, double>)
    {
        constexpr double value = 0;
        EXPECT_THROW(writeFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, int32_t>)
    {
        constexpr int32_t value = 0;
        EXPECT_THROW(writeFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, int64_t>)
    {
        constexpr int64_t value = 0;
        EXPECT_THROW(writeFrame(dataSet, 0, value), gmx::FileIOError);
    }
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueToNonValidDataSetThrows)
{
    hid_t               invalidDataSet = H5I_INVALID_HID;
    constexpr TypeParam value          = 10;
    {
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                        .withNumFrames(1) // enable reading and writing with index = 0
                        .build());

        invalidDataSet = dataSet;
        ASSERT_NO_THROW(writeFrame(invalidDataSet, 0, value))
                << "Sanity check failed: writing value must work before the guard closes the data "
                   "set";
    }

    EXPECT_THROW(writeFrame(invalidDataSet, 0, value), gmx::FileIOError)
            << "We must throw when trying to write to a closed data set";
    EXPECT_THROW(writeFrame(H5I_INVALID_HID, 0, value), gmx::FileIOError)
            << "We must throw when trying to write to a non-existing data set";
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueForReadOnlyFileThrows)
{
    // For this test we use a locally scoped file to first open it in write-mode
    // and create a data set, then reopen the file as read-only mode and try to
    // write into the existing data set which should result in a throw.
    TestFileManager       fileManager;
    std::filesystem::path fileName      = fileManager.getTemporaryFilePath("ref.h5md");
    constexpr char        dataSetName[] = "testDataSet";

    {
        SCOPED_TRACE("Create data set in write-mode file");
        H5md file(fileName, H5mdFileMode::Write);
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(H5mdFrameDataSetBuilder<TypeParam>(file.fileid(), dataSetName)
                                             .withNumFrames(1) // enable reading and writing with index = 0
                                             .build());

        constexpr TypeParam value = 10;
        ASSERT_NO_THROW(writeFrame(dataSet, 0, value))
                << "Sanity check failed: writing value to write-mode file must work";
    }
    {
        SCOPED_TRACE("Open file in read-only mode and write a value to its data set");
        H5md file(fileName, H5mdFileMode::Read);
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(openDataSet(file.fileid(), dataSetName));

        constexpr TypeParam value = 10;
        EXPECT_THROW(writeFrame(dataSet, 0, value), gmx::FileIOError)
                << "Writing value to read-only mode file must throw";
    }
}

TYPED_TEST(H5mdWriteBasicVectorListTest, WriteFrameWorks)
{
    constexpr int numAtoms             = 5;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build());
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    const std::array<gmx::BasicVector<TypeParam>, numAtoms> rvecList = {
        BasicVector<TypeParam>{ 0.0, 1.0, 2.0 },
        BasicVector<TypeParam>{ 3.0, 4.0, 5.0 },
        BasicVector<TypeParam>{ 6.0, 7.0, 8.0 },
        BasicVector<TypeParam>{ 9.0, 10.0, 11.0 },
        BasicVector<TypeParam>{ 12.0, 13.0, 14.0 }
    };
    writeFrame(dataSet, 0, makeConstArrayRef(rvecList));

    std::array<gmx::BasicVector<TypeParam>, numAtoms> readBuffer;
    readBuffer.fill(BasicVector<TypeParam>{ -1.0, -1.0, -1.0 }); // initialize buffer to fill value to avoid memory reuse
    H5Dread(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, readBuffer.data());
    EXPECT_THAT(readBuffer, ::testing::Pointwise(::testing::Eq(), rvecList));
}

TYPED_TEST(H5mdWriteBasicVectorListTest, WriteFrameToIndexWorks)
{
    constexpr int     numAtoms         = 5;
    constexpr hsize_t numFrames        = 5;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(numFrames)
                    .build());
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    // Write 5 frames to the data set in order
    std::array<std::array<gmx::BasicVector<TypeParam>, numAtoms>, numFrames> perFrameRvecLists;
    for (hsize_t frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        std::array<gmx::BasicVector<TypeParam>, numAtoms>& rvecList = perFrameRvecLists[frameIndex];
        const TypeParam frameValue = 1000 * frameIndex; // assign unique per-frame values to list
        for (int atomIndex = 0; atomIndex < gmx::ssize(rvecList); ++atomIndex)
        {
            rvecList[atomIndex][XX] = frameValue + 100 * atomIndex;
            rvecList[atomIndex][YY] = frameValue + 10 * atomIndex;
            rvecList[atomIndex][ZZ] = frameValue + atomIndex;
        }

        writeFrame(dataSet, frameIndex, makeConstArrayRef(rvecList));
    }

    // Read back all written rvec-lists into a buffer with the same data type dimensions:
    // 5 frames * 5 atoms * 3 values, then compare it to the rvec-lists which were
    // written one at a time above
    std::array<std::array<gmx::BasicVector<TypeParam>, numAtoms>, numFrames> readBufferAllFrames;
    H5Dread(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, readBufferAllFrames.data());

    for (int i = 0; i < gmx::ssize(perFrameRvecLists); ++i)
    {
        EXPECT_THAT(readBufferAllFrames[i], ::testing::Pointwise(::testing::Eq(), perFrameRvecLists[i]));
    }
}

TYPED_TEST(H5mdWriteBasicVectorListTest, WritingFrameToNonSequentialIndexWorks)
{
    constexpr int     numAtoms         = 1;
    constexpr hsize_t numFrames        = 3;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(numFrames)
                    .build());

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

    // write to the middle index (1) last to ensure it does not write into 0 or 2
    for (const int i : { 0, 2, 1 })
    {
        writeFrame(dataSet, i, makeConstArrayRef(perFrameRvecLists[i]));
    }

    for (int i = 0; i < gmx::ssize(perFrameRvecLists); ++i)
    {
        std::array<gmx::BasicVector<TypeParam>, numAtoms> readBuffer;
        readFrame(dataSet, i, makeArrayRef(readBuffer));
        EXPECT_THAT(readBuffer, ::testing::Pointwise(::testing::Eq(), perFrameRvecLists[i]));
    }
}

TYPED_TEST(H5mdWriteBasicVectorListTest, InputBufferWithNonMatchingDimensionsThrows)
{
    constexpr int numAtoms             = 3;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build());

    std::array<BasicVector<TypeParam>, numAtoms> readBufferCorrect;
    ASSERT_NO_THROW(writeFrame(dataSet, 0, makeConstArrayRef(readBufferCorrect)))
            << "Sanity check failed: reading with a correctly sized buffer must work";

    std::array<BasicVector<TypeParam>, numAtoms - 1> smallReadBuffer;
    EXPECT_THROW(writeFrame(dataSet, 0, makeConstArrayRef(smallReadBuffer)), gmx::FileIOError);
    std::array<BasicVector<TypeParam>, numAtoms + 1> largeReadBuffer;
    EXPECT_THROW(writeFrame(dataSet, 0, makeConstArrayRef(largeReadBuffer)), gmx::FileIOError);
}

TYPED_TEST(H5mdWriteBasicVectorListTest, CorrectRowOrderIsUsedForWritingArrayData)
{
    constexpr int numAtoms             = 5;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build());
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    // Write a 2d [numAtoms, DIM] array with values in row-major increasing order as a single frame
    const std::array<gmx::BasicVector<TypeParam>, numAtoms> rvecList = {
        BasicVector<TypeParam>{ 0.0, 1.0, 2.0 },
        BasicVector<TypeParam>{ 3.0, 4.0, 5.0 },
        BasicVector<TypeParam>{ 6.0, 7.0, 8.0 },
        BasicVector<TypeParam>{ 9.0, 10.0, 11.0 },
        BasicVector<TypeParam>{ 12.0, 13.0, 14.0 }
    };
    writeFrame(dataSet, 0, makeConstArrayRef(rvecList));

    // Ensure that it was written in row-major order by reading it back as a 1d array
    std::array<TypeParam, DIM * rvecList.size()> readBuffer1d;
    readBuffer1d.fill(-1.0);
    H5Dread(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, readBuffer1d.data());

    int index1d = 0;
    for (int i = 0; i < gmx::ssize(rvecList); ++i)
    {
        for (int d = 0; d < DIM; ++d)
        {
            EXPECT_REAL_EQ(rvecList[i][d], readBuffer1d[index1d]);
            index1d++;
        }
    }
}

TYPED_TEST(H5mdWriteBasicVectorListTest, WriteAfterResizingDataSetWorks)
{
    constexpr int numAtoms             = 1;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build());

    const std::array<gmx::BasicVector<TypeParam>, numAtoms> writeBuffer;
    ASSERT_NO_THROW(writeFrame(dataSet, 0, makeConstArrayRef(writeBuffer)))
            << "Sanity check failed: must be able to write to index 0 before resize";
    ASSERT_THROW(writeFrame(dataSet, 1, makeConstArrayRef(writeBuffer)), gmx::FileIOError)
            << "Sanity check failed: must throw when trying to write into index 1 before resize";

    setNumFrames(dataSet, 2);
    ASSERT_NO_THROW(writeFrame(dataSet, 1, makeConstArrayRef(writeBuffer)))
            << "Must be able to write to index 1 after resize";
}

TYPED_TEST(H5mdWriteBasicVectorListTest, WriteSingleValueToBasicVectorListDataSetThrows)
{
    const auto [rvecListDataSet, rvecListDataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ 1 })
                    .withNumFrames(1) // enable reading and writing with index = 0
                    .build());

    const TypeParam writeBuffer = 1.0;
    EXPECT_THROW(writeFrame<TypeParam>(rvecListDataSet, 0, writeBuffer), gmx::FileIOError);
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteBasicVectorListToSingleValueDataSetThrows)
{
    // Only compile when reading ArrayRef<BasicVector<real>> from a data set of <real> type:
    // reading functions for other ArrayRef<BasicVector<T>> types are not implemented.
    if constexpr (std::is_same_v<TypeParam, float> || std::is_same_v<TypeParam, double>)
    {
        const auto [primitiveTypeDataSet, primitiveDataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                        .withNumFrames(1) // enable reading and writing with index = 0
                        .build());

        const std::array<BasicVector<TypeParam>, 1> writeBuffer;
        EXPECT_THROW(writeFrame(primitiveTypeDataSet, 0, makeConstArrayRef(writeBuffer)), gmx::FileIOError);
    }
}

} // namespace
} // namespace test
} // namespace gmx
