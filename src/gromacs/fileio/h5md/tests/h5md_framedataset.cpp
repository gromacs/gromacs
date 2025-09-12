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
 * \brief Tests for H5MD frame data set routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_framedataset.h"

#include <hdf5.h>

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

/**@{
 * \brief Test suite for reading values from data sets with a parametrized frame dimension.
 *
 * Our test coverage should assure that data set reading is consistent for a suitable
 * set of frame dimensions and templated types. These tests are parametrized over a set
 * of frame dimensions, ranging from scalar-array data sets (1d) to multidimensional
 * sets.
 *
 * Individual test cases cover templated types which are primitive and compound types
 * to ensure that our supported types are tested. Types which are tested below are:
 *
 *  - int32_t:            Covers all basic primitives: int32_t, int64_t, float, double
 *  - BasicVector<float>  Covers BasicVector<T> for T = float, double
 */
struct TestFrameDimensions
{
    TestFrameDimensions(std::initializer_list<hsize_t> initList) : frameDims_(initList)
    {
        for (hsize_t d : frameDims_)
        {
            numValuesPerFrame_ *= d;
        }
    }

    DataSetDims frameDims_;

    int numValuesPerFrame_ = 1;
};

//! \brief Helper function for GTest to print dimension parameters.
void PrintTo(const TestFrameDimensions& info, std::ostream* os)
{
    const auto toString = [](const hsize_t value) -> std::string { return std::to_string(value); };
    *os << "FrameDims {" << gmx::formatAndJoin(info.frameDims_, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfTest(const ::testing::TestParamInfo<TestFrameDimensions>& info)
{
    const auto toString = [](const hsize_t value) -> std::string { return std::to_string(value); };

    std::string testName = "FrameDims_";
    if (info.param.frameDims_.empty())
    {
        testName.append("Empty");
    }
    else
    {
        testName.append(gmx::formatAndJoin(info.param.frameDims_, "_", toString));
    }

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");

    return testName;
}

//! Test fixture inheriting parametrization capabilities.
class WithFrameDims : public H5mdTestBase, public ::testing::WithParamInterface<TestFrameDimensions>
{
};

/**@{
 * \brief Test H5mdFrameDataSet<T>.readFrame() for int32_t and BasicVector<float>
 */
TEST_P(WithFrameDims, WriteAndReadFrameWorksForPrimitiveDataSets)
{
    using ValueType                       = int32_t;
    const std::vector<hsize_t>  frameDims = GetParam().frameDims_;
    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    // write values 0, 1, 2, ... into the data set as frames
    std::vector<ValueType> valuesFrame0(GetParam().numValuesPerFrame_, 0);
    std::iota(valuesFrame0.begin(), valuesFrame0.end(), 0);
    std::vector<ValueType> valuesFrame1(GetParam().numValuesPerFrame_, 0);
    std::iota(valuesFrame0.begin(), valuesFrame0.end(), GetParam().numValuesPerFrame_);

    ASSERT_EQ(dataSet.numFrames(), 0) << "Sanity check failed: numFrames() should be 0";
    dataSet.writeNextFrame(valuesFrame0);
    EXPECT_EQ(dataSet.numFrames(), 1) << "Error writing frame 0: numFrames() should be 1";
    dataSet.writeNextFrame(valuesFrame1);
    EXPECT_EQ(dataSet.numFrames(), 2) << "Error writing frame 1: numFrames() should be 2";

    std::vector<ValueType> readFrameBuffer(GetParam().numValuesPerFrame_, 0);
    dataSet.readFrame(0, readFrameBuffer);
    EXPECT_EQ(readFrameBuffer, valuesFrame0) << "Error reading frame 0";
    dataSet.readFrame(1, readFrameBuffer);
    EXPECT_EQ(readFrameBuffer, valuesFrame1) << "Error reading frame 1";
}
TEST_P(WithFrameDims, WriteAndReadFrameWorksForVectorDataSets)
{
    using ValueType                      = gmx::BasicVector<float>;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    // write values 0, 1, 2, ... into the data set as two frames
    std::vector<ValueType> valuesFrame0;
    valuesFrame0.reserve(GetParam().numValuesPerFrame_);
    for (int i = 0; i < GetParam().numValuesPerFrame_; ++i)
    {
        valuesFrame0.push_back(
                { static_cast<float>(i), static_cast<float>(100 * i), static_cast<float>(10000 * i) });
    }

    std::vector<ValueType> valuesFrame1;
    valuesFrame1.reserve(GetParam().numValuesPerFrame_);
    for (int i = 0; i < GetParam().numValuesPerFrame_; ++i)
    {
        // Multiply values by 2 to separate them from frame 0
        valuesFrame1.push_back({ static_cast<float>(2 * i),
                                 static_cast<float>(2 * 100 * i),
                                 static_cast<float>(2 * 10000 * i) });
    }

    ASSERT_EQ(dataSet.numFrames(), 0) << "Sanity check failed: numFrames() should be 0";
    dataSet.writeNextFrame(valuesFrame0);
    EXPECT_EQ(dataSet.numFrames(), 1) << "Error writing frame 0: numFrames() should be 1";
    dataSet.writeNextFrame(valuesFrame1);
    EXPECT_EQ(dataSet.numFrames(), 2) << "Error writing frame 1: numFrames() should be 2";

    std::vector<ValueType> readFrameBuffer;
    readFrameBuffer.resize(GetParam().numValuesPerFrame_);
    dataSet.readFrame(0, readFrameBuffer);
    EXPECT_THAT(readFrameBuffer, ::testing::Pointwise(::testing::Eq(), valuesFrame0))
            << "Error reading frame 0";
    dataSet.readFrame(1, readFrameBuffer);
    EXPECT_THAT(readFrameBuffer, ::testing::Pointwise(::testing::Eq(), valuesFrame1))
            << "Error reading frame 1";
}
/**@}*/

/**@{
 * \brief Test non-sequential read indices with H5mdFrameDataSet<T>.readFrame() for int32_t and BasicVector<float>
 */
TEST_P(WithFrameDims, ReadFrameWorksForNonSequentialIndicesPrimitive)
{
    using ValueType                      = int32_t;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;
    constexpr int              numFrames = 4;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    std::vector<std::vector<ValueType>> valuesPerFrame;
    for (int frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        std::vector<ValueType> values(GetParam().numValuesPerFrame_, 0);
        std::iota(values.begin(), values.end(), 1000 * frameIndex);
        valuesPerFrame.push_back(values);
        dataSet.writeNextFrame(values);
    }

    for (const int frameIndex : { 3, 1, 0, 2 })
    {
        ASSERT_TRUE(frameIndex < numFrames) << "Sanity check failed: too large frame index";
        std::vector<ValueType> readBuffer(GetParam().numValuesPerFrame_, 0);
        dataSet.readFrame(frameIndex, readBuffer);
        EXPECT_EQ(readBuffer, valuesPerFrame[frameIndex]);
    }
}
TEST_P(WithFrameDims, ReadFrameWorksForNonSequentialIndicesVector)
{
    using ValueType                      = gmx::BasicVector<float>;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;
    constexpr int              numFrames = 4;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    std::vector<std::vector<ValueType>> valuesPerFrame;
    for (int frameIndex = 0; frameIndex < numFrames; ++frameIndex)
    {
        std::vector<ValueType> values;
        values.reserve(GetParam().numValuesPerFrame_);
        for (int i = 0; i < GetParam().numValuesPerFrame_; ++i)
        {
            values.push_back({ static_cast<float>((1000 * frameIndex) + (100 * i)),
                               static_cast<float>((1000 * frameIndex) + (10 * i)),
                               static_cast<float>((1000 * frameIndex) + i) });
        }

        valuesPerFrame.push_back(values);
        dataSet.writeNextFrame(values);
    }

    for (const int frameIndex : { 3, 1, 0, 2 })
    {
        ASSERT_TRUE(frameIndex < numFrames) << "Sanity check failed: too large frame index";
        std::vector<ValueType> readBuffer(GetParam().numValuesPerFrame_, { 0.0, 0.0, 0.0 });
        dataSet.readFrame(frameIndex, readBuffer);
        EXPECT_EQ(readBuffer, valuesPerFrame[frameIndex]);
    }
}
/**@}*/

/**@{
 * \brief Test incorrect read buffer size with H5mdFrameDataSet<T>.readFrame() for int32_t and BasicVector<float>
 */
TEST_P(WithFrameDims, ReadFrameThrowsForNonMatchingInputBufferSizePrimitive)
{
    using ValueType                      = int32_t;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    std::vector<ValueType> values(GetParam().numValuesPerFrame_, 0);
    std::iota(values.begin(), values.end(), 0);
    dataSet.writeNextFrame(values);

    std::vector<ValueType> readBufferEmpty;
    EXPECT_THROW(dataSet.readFrame(0, readBufferEmpty), gmx::FileIOError);
    std::vector<ValueType> readBufferTooSmall(GetParam().numValuesPerFrame_ - 1, 0);
    EXPECT_THROW(dataSet.readFrame(0, readBufferTooSmall), gmx::FileIOError);
    std::vector<ValueType> readBufferTooLarge(GetParam().numValuesPerFrame_ + 1, 0);
    EXPECT_THROW(dataSet.readFrame(0, readBufferTooLarge), gmx::FileIOError);
    std::vector<ValueType> readBufferJustRight(GetParam().numValuesPerFrame_, 0);
    EXPECT_NO_THROW(dataSet.readFrame(0, readBufferJustRight));
}
TEST_P(WithFrameDims, ReadFrameThrowsForNonMatchingInputBufferSizeVector)
{
    using ValueType                      = gmx::BasicVector<float>;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .build();

    std::vector<ValueType> values;
    values.reserve(GetParam().numValuesPerFrame_);
    for (int i = 0; i < GetParam().numValuesPerFrame_; ++i)
    {
        values.push_back(
                { static_cast<float>((100 * i)), static_cast<float>((10 * i)), static_cast<float>(i) });
    }
    dataSet.writeNextFrame(values);

    std::vector<ValueType> readBufferEmpty;
    EXPECT_THROW(dataSet.readFrame(0, readBufferEmpty), gmx::FileIOError);
    std::vector<ValueType> readBufferTooSmall(GetParam().numValuesPerFrame_ - 1, { 0.0, 0.0, 0.0 });
    EXPECT_THROW(dataSet.readFrame(0, readBufferTooSmall), gmx::FileIOError);
    std::vector<ValueType> readBufferTooLarge(GetParam().numValuesPerFrame_ + 1, { 0.0, 0.0, 0.0 });
    EXPECT_THROW(dataSet.readFrame(0, readBufferTooLarge), gmx::FileIOError);
    std::vector<ValueType> readBufferJustRight(GetParam().numValuesPerFrame_, { 0.0, 0.0, 0.0 });
    EXPECT_NO_THROW(dataSet.readFrame(0, readBufferJustRight));
}
/**@}*/

/**@{
 * \brief Test incorrect frame indices with H5mdFrameDataSet<T>.readFrame() for int32_t and BasicVector<float>
 */
TEST_P(WithFrameDims, ReadFrameThrowsForInvalidFrameIndexPrimitive)
{
    using ValueType                      = int32_t;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;
    constexpr int              numFrames = 5;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .withNumFrames(numFrames)
                    .build();

    std::vector<ValueType> readBuffer(GetParam().numValuesPerFrame_, 0);
    EXPECT_NO_THROW(dataSet.readFrame(numFrames - 1, readBuffer));
    EXPECT_THROW(dataSet.readFrame(numFrames, readBuffer), gmx::FileIOError);
    EXPECT_THROW(dataSet.readFrame(numFrames + 1, readBuffer), gmx::FileIOError);
}
TEST_P(WithFrameDims, ReadFrameThrowsForInvalidFrameIndexVector)
{
    using ValueType                      = gmx::BasicVector<float>;
    const std::vector<hsize_t> frameDims = GetParam().frameDims_;
    constexpr int              numFrames = 5;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDims)
                    .withNumFrames(numFrames)
                    .build();

    std::vector<ValueType> readBuffer(GetParam().numValuesPerFrame_, { 0.0, 0.0, 0.0 });
    EXPECT_NO_THROW(dataSet.readFrame(numFrames - 1, readBuffer));
    EXPECT_THROW(dataSet.readFrame(numFrames, readBuffer), gmx::FileIOError);
    EXPECT_THROW(dataSet.readFrame(numFrames + 1, readBuffer), gmx::FileIOError);
}
/**@}*/

//! \brief Set of frame dimension parameters to instantiate test suite for.
const TestFrameDimensions g_testFrameDims[] = { {},    // scalar data set: no frame dimensions
                                                { 5 }, // single frame dimension, e.g. T[numAtoms]
                                                { 5, 3 } }; // multi-dim frame dimension, e.g. T[numAtoms][DIM]
INSTANTIATE_TEST_SUITE_P(H5mdReadTest, WithFrameDims, ::testing::ValuesIn(g_testFrameDims), nameOfTest);
/**@}*/

//! Test fixture for non-parametrized frame data set tests.
using H5mdFrameDataSetTest = H5mdTestBase;

TEST_F(H5mdFrameDataSetTest, WriteNextFrameIndexingBeginsAtZeroForNewDataSets)
{
    using ValueType = int32_t;

    constexpr int               numFrames = 5;
    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .build();

    ValueType readBuffer;

    EXPECT_THROW(dataSet.readFrame(numFrames, arrayRefFromArray(&readBuffer, 1)), gmx::FileIOError)
            << "Frame at created index should not exist before call to writeNextFrame()";

    ValueType valueToWrite = 7;
    dataSet.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1));

    EXPECT_NO_THROW(dataSet.readFrame(numFrames, arrayRefFromArray(&readBuffer, 1)))
            << "Read from created index should work after call to writeNextFrame()";
    EXPECT_EQ(readBuffer, valueToWrite);
}

TEST_F(H5mdFrameDataSetTest, WriteNextFrameIndexingWorksAfterOpen)
{
    using ValueType         = int32_t;
    constexpr int numFrames = 5;

    {
        SCOPED_TRACE("Create a data set, write numFrames values to it and close at end of scope");
        H5mdFrameDataSet<ValueType> dataSet =
                H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet").build();

        for (ValueType i = 0; i < numFrames; ++i)
        {
            dataSet.writeNextFrame(constArrayRefFromArray(&i, 1));
        }
    }
    {
        SCOPED_TRACE("Open data set and assert that writeNextIndex write into index = numFrames");
        H5mdFrameDataSet<ValueType> dataSet(this->fileid(), "testDataSet");

        ValueType value = 103;
        dataSet.writeNextFrame(constArrayRefFromArray(&value, 1));

        ValueType readBuffer;
        dataSet.readFrame(numFrames, arrayRefFromArray(&readBuffer, 1));
        EXPECT_EQ(readBuffer, value);
    }
}

TEST_F(H5mdFrameDataSetTest, WriteAfterMaxNumFramesThrows)
{
    using ValueType              = int32_t;
    constexpr int   maxNumFrames = 2;
    const ValueType valueToWrite = 50;

    H5mdFrameDataSet<ValueType> dataSet =
            H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                    .withMaxNumFrames(maxNumFrames)
                    .build();
    for (ValueType i = 0; i < maxNumFrames; ++i)
    {
        dataSet.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1));
    }

    EXPECT_THROW(dataSet.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1)), gmx::FileIOError)
            << "Must throw when writing more than maxNumFrames to a data set";
}

TEST_F(H5mdFrameDataSetTest, WriteAfterMaxNumFramesThrowsAfterReopen)
{
    using ValueType              = int32_t;
    constexpr int   maxNumFrames = 2;
    const ValueType valueToWrite = 50;

    {
        SCOPED_TRACE("Create data set and write maxNumFrames to it, then close it at end of scope");
        H5mdFrameDataSet<ValueType> dataSet =
                H5mdFrameDataSetBuilder<ValueType>(this->fileid(), "testDataSet")
                        .withMaxNumFrames(maxNumFrames)
                        .build();
        for (ValueType i = 0; i < maxNumFrames; ++i)
        {
            dataSet.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1));
        }
    }
    {
        SCOPED_TRACE("Reopen the data set and try to write another frame to it");
        H5mdFrameDataSet<ValueType> dataSet(fileid(), "testDataSet");
        EXPECT_THROW(dataSet.writeNextFrame(constArrayRefFromArray(&valueToWrite, 1)), gmx::FileIOError)
                << "Must throw when writing more than maxNumFrames to a data set";
    }
}


} // namespace
} // namespace test
} // namespace gmx
