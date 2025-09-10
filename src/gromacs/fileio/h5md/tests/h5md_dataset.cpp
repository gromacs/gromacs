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

#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for all relevant data set types
template<typename ValueType>
class H5mdDataSetTest : public H5mdTestBase
{
};

//! \brief List of all data types to create tests for
using DataTypesToTest =
        ::testing::Types<int32_t, int64_t, float, double, gmx::BasicVector<float>, gmx::BasicVector<double>>;
TYPED_TEST_SUITE(H5mdDataSetTest, DataTypesToTest);

TYPED_TEST(H5mdDataSetTest, NumFramesFor1dDataSetsReturnsSize)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();

    const hsize_t newSize = 6;
    setNumFrames(dataSet, newSize);
    EXPECT_EQ(getNumFrames(dataSet), newSize) << "Number of frames should match new size";
}

TYPED_TEST(H5mdDataSetTest, NumFramesFor3dDataSetsReturnsDim0Value)
{
    const std::vector<hsize_t>       frameDimensions = { 5, 2 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDimensions)
                    .build();

    const hsize_t newSize = 6;
    setNumFrames(dataSet, newSize);
    EXPECT_EQ(getNumFrames(dataSet), newSize)
            << "Number of frames should match new size along dim 0";
}

TYPED_TEST(H5mdDataSetTest, SetNumFramesWorks)
{
    H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();

    DataSetDims expectedDims = dataSet.dims();
    ASSERT_EQ(expectedDims[0], 0) << "Sanity check failed: unexpected initial size of data set";

    constexpr hsize_t newNumFrames = 6;
    setNumFrames(dataSet, newNumFrames);
    expectedDims[0] = newNumFrames;
    EXPECT_EQ(expectedDims, dataSet.dims()) << "Incorrect data set size after setNumFrames";
}

TYPED_TEST(H5mdDataSetTest, SetNumFramesFor3dDataSetsWorks)
{
    const std::vector<hsize_t> frameDimensions = { 5, 2 };
    H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDimensions)
                    .build();

    DataSetDims expectedDims = dataSet.dims();
    ASSERT_EQ(expectedDims[0], 0) << "Sanity check failed: unexpected initial size of data set";

    constexpr hsize_t newNumFrames = 6;
    setNumFrames(dataSet, newNumFrames);
    expectedDims[0] = newNumFrames;
    EXPECT_EQ(expectedDims, dataSet.dims()) << "Incorrect data set size after setNumFrames";
}

TYPED_TEST(H5mdDataSetTest, SetNumFramesCanShrinkDataset)
{
    H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    DataSetDims expectedDims = dataSet.dims();

    setNumFrames(dataSet, 2);
    setNumFrames(dataSet, 1);

    expectedDims[0] = 1;
    EXPECT_EQ(dataSet.dims(), expectedDims) << "Incorrect data set size after setNumFrames";
}

} // namespace
} // namespace test
} // namespace gmx
