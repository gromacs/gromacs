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
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"

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

//! \brief List of numeric primitives for which to create tests
using PrimitiveList = ::testing::Types<int32_t, int64_t, float, double>;

// Set up suites for testing of all relevant types
// Right now these are just numeric primitives but we will also have strings,
// gmx::RVecs, gmx::Matrix3x3s, etc. Some of these cannot run the same
// suite, so more may be added.
TYPED_TEST_SUITE(H5mdReadNumericPrimitiveTest, PrimitiveList);

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueFrom1dSetWorks)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));

    // Write values to data set in bulk for efficiency
    const std::vector<TypeParam> values    = { 9, 3, 1, 0, 6 };
    const hsize_t                numFrames = values.size();
    H5Dset_extent(dataSet, &numFrames);
    H5Dwrite(dataSet, hdf5DataTypeFor<TypeParam>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

    TypeParam value;
    for (hsize_t frameIndex = 0; frameIndex < values.size(); ++frameIndex)
    {
        readFrame(dataSet, frameIndex, value);
        EXPECT_FLOAT_EQ(value, values[frameIndex]);
    }
}

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueOutsideOfSetBoundsThrows)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));

    // Write values to data set in bulk for efficiency
    const std::vector<TypeParam> values    = { 9, 3, 1, 0, 6 };
    const hsize_t                numFrames = values.size();
    H5Dset_extent(dataSet, &numFrames);
    H5Dwrite(dataSet, hdf5DataTypeFor<TypeParam>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

    TypeParam value;
    EXPECT_NO_THROW(readFrame(dataSet, values.size() - 1, value))
            << "Sanity check failed: reading last value must work";
    EXPECT_THROW(readFrame(dataSet, values.size(), value), gmx::FileIOError)
            << "Must throw when reading out of bounds";
}

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueOfNonMatchingTypeThrows)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));

    // Write values to data set in bulk for efficiency
    const std::vector<TypeParam> values    = { 9, 3, 1, 0, 6 };
    const hsize_t                numFrames = values.size();
    H5Dset_extent(dataSet, &numFrames);
    H5Dwrite(dataSet, hdf5DataTypeFor<TypeParam>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

    if (!std::is_same_v<TypeParam, float>)
    {
        float value;
        EXPECT_THROW(readFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, double>)
    {
        double value;
        EXPECT_THROW(readFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, int32_t>)
    {
        int32_t value;
        EXPECT_THROW(readFrame(dataSet, 0, value), gmx::FileIOError);
    }
    if (!std::is_same_v<TypeParam, int64_t>)
    {
        int64_t value;
        EXPECT_THROW(readFrame(dataSet, 0, value), gmx::FileIOError);
    }
}

TYPED_TEST(H5mdReadNumericPrimitiveTest, ReadValueFromNonValidDataSetThrows)
{
    const auto [dataSet, dataSetGuard] =
            makeH5mdDataSetGuard(create1dFrameDataSet<TypeParam>(this->fileid(), "testDataSet"));

    TypeParam value;
    {
        // Write values to data set in bulk for efficiency
        const std::vector<TypeParam> values    = { 9, 3, 1, 0, 6 };
        const hsize_t                numFrames = values.size();
        H5Dset_extent(dataSet, &numFrames);
        H5Dwrite(dataSet, hdf5DataTypeFor<TypeParam>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());

        ASSERT_NO_THROW(readFrame(dataSet, 0, value))
                << "Sanity check failed: reading value from open set must work";
        H5Dclose(dataSet);
    }

    EXPECT_THROW(readFrame(dataSet, 0, value), gmx::FileIOError)
            << "We must throw when trying to read from a closed data set";
    EXPECT_THROW(readFrame(H5I_INVALID_HID, 0, value), gmx::FileIOError)
            << "We must throw when trying to read from a non-existing data set";
}

} // namespace
} // namespace test
} // namespace gmx
