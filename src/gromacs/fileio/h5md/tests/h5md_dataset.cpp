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
 * \brief Tests for H5MD file data set creation, opening and manipulation.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_dataset.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for all primitive types
template<typename ValueType>
class H5mdDataSetTest : public ::testing::Test
{
};

//! \brief List of primitives for which to create tests
using PrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;

// Set up suites for testing of all relevant types
TYPED_TEST_SUITE(H5mdDataSetTest, PrimitiveList);

TEST(H5mdDataSetTest, Create1dFrameDataSetCreatesEmpty1dSet)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"));

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    const int numDims                      = H5Sget_simple_extent_ndims(dataSpace);
    EXPECT_EQ(numDims, 1);

    hsize_t dataSetDims; // numDims = 1 so we only need to allocate memory for 1 dimension value
    H5Sget_simple_extent_dims(dataSpace, &dataSetDims, nullptr);
    EXPECT_EQ(dataSetDims, 0);
}

TYPED_TEST(H5mdDataSetTest, Create1dFrameDataSetsCorrectDataType)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(file.fileid(), "testDataSet"));
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    const TypeParam                 value = 1;
    const ArrayRef<const TypeParam> array = arrayRefFromArray(&value, 1);

    EXPECT_TRUE(valueTypeIsDataType(dataType, array));
}

TEST(H5mdDataSetTest, Create1dFrameDataSetSetsUnlimitedMaxSize)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"));

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    DataSetDims<1> dataSetMaxDims;
    H5Sget_simple_extent_dims(dataSpace, nullptr, dataSetMaxDims.data());

    EXPECT_EQ(dataSetMaxDims.at(0), H5S_UNLIMITED);
}

TEST(H5mdDataSetTest, Create1dFrameDataSetThrowsForNonExistingContainer)
{
    EXPECT_THROW(create1dFrameDataSet<int32_t>(H5I_INVALID_HID, "testDataSet"), gmx::FileIOError);
}

TEST(H5mdDataSetTest, Create1dFrameDataSetThrowsForEmptyName)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    EXPECT_THROW(create1dFrameDataSet<int32_t>(file.fileid(), ""), gmx::FileIOError);
}

TEST(H5mdDataSetTest, Create1dFrameDataSetThrowsForDuplicateName)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    constexpr char dataSetName[] = "testDataSet";
    create1dFrameDataSet<int32_t>(file.fileid(), dataSetName);
    EXPECT_THROW(create1dFrameDataSet<int32_t>(file.fileid(), dataSetName), gmx::FileIOError);
}

TEST(H5mdDataSetTest, Create1dFrameDataSetThrowsForReadOnlyFile)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    {
        SCOPED_TRACE("Create file.");
        H5md file(fileName, H5mdFileMode::Write);
    }
    {
        SCOPED_TRACE("Open file for reading and try to create data set.");
        H5md file(fileName, H5mdFileMode::Read);
        EXPECT_THROW(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"), gmx::FileIOError);
    }
}

TEST(H5mdDataSetTest, OpenDataSetWorksForWriteModeFiles)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md           file(fileName, H5mdFileMode::Write);
    constexpr char dataSetName[] = "testDataSet";

    {
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), dataSetName));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        ASSERT_GT(H5Tequal(dataType, H5T_NATIVE_INT32), 0)
                << "Sanity check failed: data types should match after creation";
    }
    {
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(openDataSet(file.fileid(), dataSetName));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        EXPECT_GT(H5Tequal(dataType, H5T_NATIVE_INT32), 0) << "Data types must match after opening";
    }
}

TEST(H5mdDataSetTest, OpenDataSetWorksForReadOnlyFiles)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    constexpr char dataSetName[] = "testDataSet";

    {
        H5md file(fileName, H5mdFileMode::Write);
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), dataSetName));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        ASSERT_GT(H5Tequal(dataType, H5T_NATIVE_INT32), 0)
                << "Sanity check failed: data types should match after creation";
    }
    {
        H5md file(fileName, H5mdFileMode::Read);
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(openDataSet(file.fileid(), dataSetName));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        EXPECT_GT(H5Tequal(dataType, H5T_NATIVE_INT32), 0) << "Data types must match after opening";
    }
}

TEST(H5mdDataSetTest, OpenDataSetThrowsForInvalidContainer)
{
    EXPECT_THROW(openDataSet(H5I_INVALID_HID, "testDataSet"), gmx::FileIOError);
}

TEST(H5mdDataSetTest, OpenDataSetThrowsForInvalidSetName)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    {
        // Create a data set to ensure that there is something in the file which we cannot read with bad names
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"));
        ASSERT_NE(dataSet, H5I_INVALID_HID) << "Sanity check failed: data set should be created";
    }

    EXPECT_THROW(openDataSet(file.fileid(), ""), gmx::FileIOError) << "Should throw for empty name";
    EXPECT_THROW(openDataSet(file.fileid(), "aBadIdea"), gmx::FileIOError)
            << "Should throw for bad name";
}

TEST(H5mdDataSetTest, GetNumFramesFor1dDataSetsReturnsSize)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"));

    EXPECT_EQ(getNumFrames<1>(dataSet), 0) << "1d frame data sets are created as empty";

    const hsize_t newSize = 6;
    H5Dset_extent(dataSet, &newSize);
    EXPECT_EQ(getNumFrames<1>(dataSet), newSize) << "Number of frames should match new size";
}

TEST(H5mdDataSetTest, GetNumFramesThrowsForInvalidDataSetHandle)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    H5md file(fileName, H5mdFileMode::Write);

    // We use a scope guard to create a data set and store its handle here; as the scope exits
    // the data set is closed and this turns invalid
    hid_t dataSetToTest = H5I_INVALID_HID;

    {
        // Create a data set with some data in the file to ensure that there is something
        // there which cannot be read with bad data set handles.
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(create1dFrameDataSet<int32_t>(file.fileid(), "testDataSet"));

        dataSetToTest = dataSet;
        ASSERT_NE(dataSetToTest, H5I_INVALID_HID)
                << "Sanity check failed: data set handle must be valid before exiting scope";
    }

    EXPECT_THROW(getNumFrames<1>(dataSetToTest), gmx::FileIOError);
    EXPECT_THROW(getNumFrames<1>(H5I_INVALID_HID), gmx::FileIOError);
}

} // namespace
} // namespace test
} // namespace gmx
