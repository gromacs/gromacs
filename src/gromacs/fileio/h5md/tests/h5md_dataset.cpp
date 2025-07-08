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
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/math/vectypes.h"
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
class H5mdPrimitiveDataSetTest : public H5mdTestBase
{
};

template<typename ValueType>
class H5mdBasicVectorListDataSetTest : public H5mdTestBase
{
};

//! \brief List of primitives for which to create tests
using NumericPrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;
TYPED_TEST_SUITE(H5mdPrimitiveDataSetTest, NumericPrimitiveList);

//! \brief List of numeric primitives for which to create real-valued tests
using RealPrimitiveList = ::testing::Types<float, double>;
TYPED_TEST_SUITE(H5mdBasicVectorListDataSetTest, RealPrimitiveList);

TYPED_TEST(H5mdPrimitiveDataSetTest, OpenDataSetWorksForWriteModeFiles)
{
    constexpr char dataSetName[] = "testDataSet";

    {
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), dataSetName).build());
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        ASSERT_GT(H5Tequal(dataType, hdf5DataTypeFor<TypeParam>()), 0)
                << "Sanity check failed: data types should match after creation";
    }
    {
        const auto [dataSet, dataSetGuard] =
                makeH5mdDataSetGuard(openDataSet(this->fileid(), dataSetName));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        EXPECT_GT(H5Tequal(dataType, hdf5DataTypeFor<TypeParam>()), 0)
                << "Data types must match after opening";
    }
}

TEST(H5mdDataSetTest, OpenDataSetWorksForReadOnlyFiles)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    constexpr char dataSetName[] = "testDataSet";

    {
        H5md file(fileName, H5mdFileMode::Write);
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<int32_t>(file.fileid(), dataSetName).build());
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

TYPED_TEST(H5mdPrimitiveDataSetTest, OpenDataSetThrowsForInvalidSetName)
{
    {
        // Create a data set to ensure that there is something in the file which we cannot read with bad names
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build());
        ASSERT_TRUE(handleIsValid(dataSet))
                << "Sanity check failed: data set handle should be valid";
    }

    EXPECT_THROW(openDataSet(this->fileid(), ""), gmx::FileIOError)
            << "Should throw for empty name";
    EXPECT_THROW(openDataSet(this->fileid(), "aBadIdea"), gmx::FileIOError)
            << "Should throw for bad name";
}

TYPED_TEST(H5mdPrimitiveDataSetTest, GetNumFramesFor1dDataSetsReturnsSize)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build());

    const hsize_t newSize = 6;
    H5Dset_extent(dataSet, &newSize);
    EXPECT_EQ(getNumFrames(dataSet), newSize) << "Number of frames should match new size";
}

TYPED_TEST(H5mdPrimitiveDataSetTest, GetNumFramesFor3dDataSetsReturnsDim0Value)
{
    const std::vector<hsize_t> frameDimensions = { 5, 2 };
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                                         .withFrameDimension(frameDimensions)
                                         .build());

    const hsize_t     newSize = 6;
    const DataSetDims newDims = { newSize, frameDimensions[0], frameDimensions[1] };
    H5Dset_extent(dataSet, newDims.data());
    EXPECT_EQ(getNumFrames(dataSet), newSize)
            << "Number of frames should match new size along dim 0";
}

TYPED_TEST(H5mdPrimitiveDataSetTest, GetNumFramesThrowsForInvalidDataSetHandle)
{
    // We use a scope guard to create a data set and store its handle here; as the scope exits
    // the data set is closed and this turns invalid
    hid_t dataSetToTest = H5I_INVALID_HID;

    {
        // Create a data set with some data in the file to ensure that there is something
        // there which cannot be read with bad data set handles.
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build());

        dataSetToTest = dataSet;
        ASSERT_TRUE(handleIsValid(dataSetToTest))
                << "Sanity check failed: data set handle must be valid before exiting scope";
    }

    EXPECT_THROW(getNumFrames(dataSetToTest), gmx::FileIOError);
    EXPECT_THROW(getNumFrames(H5I_INVALID_HID), gmx::FileIOError);
}

TYPED_TEST(H5mdPrimitiveDataSetTest, SetNumFramesWorks)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build());

    {
        SCOPED_TRACE("Validate data set size before setting number of frames");
        const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
        DataSetDims dataSetDims(1, 0);
        H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr);
        ASSERT_EQ(dataSetDims, (DataSetDims{ 0 }))
                << "Sanity check failed: unexpected initial size of data set";
    }

    constexpr hsize_t newNumFrames = 2;
    setNumFrames(dataSet, newNumFrames);
    // The data space is invalidated after changing the data set size, so we cannot
    // reuse the one created for checking the initial size above
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    DataSetDims dataSetDims(1, 0);
    H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr);
    EXPECT_EQ(dataSetDims, (DataSetDims{ newNumFrames }))
            << "Incorrect data set size after setNumFrames";
}

TYPED_TEST(H5mdPrimitiveDataSetTest, SetNumFramesCanShrinkDataset)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build());

    setNumFrames(dataSet, 2);
    setNumFrames(dataSet, 1);

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    DataSetDims dataSetDims(1, 0);
    H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr);
    EXPECT_EQ(dataSetDims, (DataSetDims{ 1 })) << "Incorrect data set size after setNumFrames";
}

TYPED_TEST(H5mdBasicVectorListDataSetTest, SetNumFramesWorksForBasicVector)
{
    constexpr int numAtoms             = 5;
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withFrameDimension({ numAtoms })
                    .build());

    {
        SCOPED_TRACE("Validate data set size before setting number of frames");
        const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
        DataSetDims dataSetDims(3, 0);
        H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr);
        ASSERT_EQ(dataSetDims, (DataSetDims{ 0, numAtoms, DIM }))
                << "Sanity check failed: unexpected initial size of data set";
    }

    constexpr hsize_t newNumFrames = 2;
    setNumFrames(dataSet, newNumFrames);
    // The data space is invalidated after changing the data set size, so we cannot
    // reuse the one created for checking the initial size above
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    DataSetDims dataSetDims(3, 0);
    H5Sget_simple_extent_dims(dataSpace, dataSetDims.data(), nullptr);
    EXPECT_EQ(dataSetDims, (DataSetDims{ newNumFrames, numAtoms, DIM }))
            << "Incorrect data set size after setNumFrames";
}

} // namespace
} // namespace test
} // namespace gmx
