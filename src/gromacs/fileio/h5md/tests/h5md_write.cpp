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
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"

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

//! \brief List of numeric primitives for which to create tests
using PrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;

// Set up suites for testing of all relevant types
// Right now these are just numeric primitives but we will also have strings,
// gmx::RVecs, gmx::Matrix3x3s, etc. Some of these cannot run the same
// suite, so more may be added.
TYPED_TEST_SUITE(H5mdWriteNumericPrimitiveTest, PrimitiveList);

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueTo1dSetWritesToCorrectIndex)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));

    // List of values to write into data set one-by-one with our tested function
    const std::vector<TypeParam> valuesToWrite = { 9, 3, 1, 0, 6 };
    const hsize_t                numFrames     = valuesToWrite.size();
    H5Dset_extent(dataSet, &numFrames);

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
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));
    constexpr hsize_t numFrames = 5;
    H5Dset_extent(dataSet, &numFrames);

    constexpr TypeParam value = 0;
    ASSERT_NO_THROW(writeFrame(dataSet, numFrames - 1, value))
            << "Sanity check failed: write last value must work";
    EXPECT_THROW(writeFrame(dataSet, numFrames, value), gmx::FileIOError)
            << "Must throw when writing out of bounds";
}

TYPED_TEST(H5mdWriteNumericPrimitiveTest, WriteValueOfNonMatchingTypeThrows)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));
    constexpr hsize_t numFrames = 1; // enable writing into index = 0
    H5Dset_extent(dataSet, &numFrames);

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
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));
        constexpr hsize_t numFrames = 1; // enable writing into index = 0
        H5Dset_extent(dataSet, &numFrames);

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
                makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(file.fileid(), dataSetName));

        constexpr hsize_t numFrames = 1;
        H5Dset_extent(dataSet, &numFrames);
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

} // namespace
} // namespace test
} // namespace gmx
