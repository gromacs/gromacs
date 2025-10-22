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
 * \brief Tests for H5MD frame data set builder routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

hsize_t dataSetSize(const hid_t hdf, const char* datasetName)
{
    const auto dataSet                     = H5mdDataSetBase<std::string>(hdf, datasetName);
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    std::vector<hsize_t> dims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    return dims[0];
}

hsize_t dataSetCapacity(const hid_t hdf, const char* datasetName)
{
    const auto dataSet                     = H5mdDataSetBase<std::string>(hdf, datasetName);
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    std::vector<hsize_t> maxDims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, nullptr, maxDims.data());
    return maxDims[0];
}

template<typename ValueType>
class H5mdFrameDataSetBuilderTest : public H5mdTestBase
{
};

template<typename ValueType>
class H5mdNumericPrimitiveFrameDataSetBuilderTest : public H5mdTestBase
{
};

template<typename ValueType>
class H5mdBasicVectorFrameDataSetBuilderTest : public H5mdTestBase
{
};

using DataTypes =
        ::testing::Types<int32_t, int64_t, float, double, gmx::BasicVector<float>, gmx::BasicVector<double>>;
TYPED_TEST_SUITE(H5mdFrameDataSetBuilderTest, DataTypes);

using NumericPrimitives = ::testing::Types<int32_t, int64_t, float, double>;
TYPED_TEST_SUITE(H5mdNumericPrimitiveFrameDataSetBuilderTest, NumericPrimitives);

using RealPrimitives = ::testing::Types<float, double>;
TYPED_TEST_SUITE(H5mdBasicVectorFrameDataSetBuilderTest, RealPrimitives);

TYPED_TEST(H5mdFrameDataSetBuilderTest, DefaultNumFramesIsZero)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    std::vector<hsize_t> dims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], 0) << "Incorrect number of frames in data set";
}

TYPED_TEST(H5mdFrameDataSetBuilderTest, SetNumFramesWorks)
{
    constexpr int                    numFrames = 5;
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    std::vector<hsize_t> dims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames) << "Incorrect number of frames in data set";
}

TYPED_TEST(H5mdFrameDataSetBuilderTest, DefaultMaxNumFramesIsUnlimited)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    std::vector<hsize_t> maxDims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, nullptr, maxDims.data());
    EXPECT_EQ(maxDims[0], H5S_UNLIMITED) << "Incorrect maximum number of frames in data set";
}

TYPED_TEST(H5mdFrameDataSetBuilderTest, SetMaxNumFramesWorks)
{
    constexpr int                    maxNumFrames = 5;
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withMaxNumFrames(maxNumFrames)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    std::vector<hsize_t> maxDims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, nullptr, maxDims.data());
    EXPECT_EQ(maxDims[0], maxNumFrames) << "Incorrect maximum number of frames in data set";
}

using H5mdStringFrameDataSetBuilderTest = H5mdTestBase;

TEST_F(H5mdStringFrameDataSetBuilderTest, StringDSetWithDefaultFrameNumber)
{
    {
        SCOPED_TRACE("Check variable length string data set");

        EXPECT_NO_THROW(
                H5mdFrameDataSetBuilder<std::string>(this->fileid(), "DataSetWithVariableLength")
                        .withVariableStringLength()
                        .build());

        const auto dataSet =
                H5mdDataSetBase<std::string>(this->fileid(), "DataSetWithVariableLength");

        EXPECT_EQ(H5Tget_class(dataSet.dataType()), H5T_STRING);
        EXPECT_TRUE(H5Tis_variable_str(dataSet.dataType()));
        EXPECT_EQ(dataSetSize(this->fileid(), "DataSetWithVariableLength"), 0)
                << "Default number of frames should be 0";
        EXPECT_EQ(dataSetCapacity(this->fileid(), "DataSetWithVariableLength"), H5S_UNLIMITED)
                << "Default maximum number of frames should be unlimited";
    }

    {
        SCOPED_TRACE("Check fixed length string data set");

        const auto dataSet =
                H5mdFrameDataSetBuilder<std::string>(this->fileid(), "DataSetWithMaxLength")
                        .withMaxStringLength(10)
                        .build();

        EXPECT_EQ(H5Tget_class(dataSet.dataType()), H5T_STRING);
        EXPECT_EQ(H5Tget_size(dataSet.dataType()), 10);
        EXPECT_EQ(dataSetSize(this->fileid(), "DataSetWithMaxLength"), 0)
                << "Default number of frames should be 0";
        EXPECT_EQ(dataSetCapacity(this->fileid(), "DataSetWithMaxLength"), H5S_UNLIMITED)
                << "Default maximum number of frames should be unlimited";
    }
}


TEST_F(H5mdStringFrameDataSetBuilderTest, StringDataSetWithMaxFrames)
{
    const int numFrames = 5;
    {
        SCOPED_TRACE("Check variable length string data set");

        const auto dataSet =
                H5mdFrameDataSetBuilder<std::string>(this->fileid(), "DataSetWithMaxLength")
                        .withMaxStringLength(10)
                        .withNumFrames(numFrames)
                        .withMaxNumFrames(numFrames)
                        .build();

        EXPECT_EQ(H5Tget_class(dataSet.dataType()), H5T_STRING);
        EXPECT_EQ(H5Tget_size(dataSet.dataType()), 10);
        EXPECT_EQ(dataSetSize(this->fileid(), "DataSetWithMaxLength"), numFrames)
                << "Default number of frames should be 0";
        EXPECT_EQ(dataSetCapacity(this->fileid(), "DataSetWithMaxLength"), numFrames)
                << "Default maximum number of frames should be unlimited";
    }

    {
        SCOPED_TRACE("Check fixed length string data set");

        const auto dataSet =
                H5mdFrameDataSetBuilder<std::string>(this->fileid(), "DataSetWithVariableLength")
                        .withVariableStringLength()
                        .withNumFrames(numFrames)
                        .withMaxNumFrames(numFrames)
                        .build();

        EXPECT_EQ(H5Tget_class(dataSet.dataType()), H5T_STRING);
        EXPECT_TRUE(H5Tis_variable_str(dataSet.dataType()));
        EXPECT_EQ(dataSetSize(this->fileid(), "DataSetWithVariableLength"), numFrames)
                << "Default number of frames should be 0";
        EXPECT_EQ(dataSetCapacity(this->fileid(), "DataSetWithVariableLength"), numFrames)
                << "Default maximum number of frames should be unlimited";
    }
}

TYPED_TEST(H5mdNumericPrimitiveFrameDataSetBuilderTest, DataTypeIsCorrect)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    EXPECT_TRUE(valueTypeIsDataType<TypeParam>(dataSet.dataType()));
}

TYPED_TEST(H5mdNumericPrimitiveFrameDataSetBuilderTest, DefaultDimsIs1d)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    EXPECT_EQ(H5Sget_simple_extent_ndims(dataSpace), 1)
            << "Primitive data sets should be 1d by default";
}

TYPED_TEST(H5mdNumericPrimitiveFrameDataSetBuilderTest, SetFrameDimensionAddsDimension1)
{
    constexpr int                    numFrames      = 3;
    const std::vector<hsize_t>       frameDimension = { 5 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .withFrameDimension(frameDimension)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    ASSERT_EQ(H5Sget_simple_extent_ndims(dataSpace), 2)
            << "Adding one dimension should result in 2d";
    std::vector<hsize_t> dims(2, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames);
    EXPECT_EQ(dims[1], frameDimension[0]);
}

TYPED_TEST(H5mdNumericPrimitiveFrameDataSetBuilderTest, SetFrameDimensionAddsDimension2)
{
    constexpr int                    numFrames      = 3;
    const std::vector<hsize_t>       frameDimension = { 5, 7 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .withFrameDimension(frameDimension)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    ASSERT_EQ(H5Sget_simple_extent_ndims(dataSpace), 3)
            << "Adding two dimensions should result in 3d";
    std::vector<hsize_t> dims(3, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames);
    EXPECT_EQ(dims[1], frameDimension[0]);
    EXPECT_EQ(dims[2], frameDimension[1]);
}

TYPED_TEST(H5mdBasicVectorFrameDataSetBuilderTest, DataTypeIsCorrect)
{
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet").build();
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet.id()));
    EXPECT_TRUE(valueTypeIsDataType<TypeParam>(dataType));
}

TYPED_TEST(H5mdBasicVectorFrameDataSetBuilderTest, DefaultDimsIs2d)
{
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet").build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    EXPECT_EQ(H5Sget_simple_extent_ndims(dataSpace), 2)
            << "BasicVector data sets should be 2d by default";
}

TYPED_TEST(H5mdBasicVectorFrameDataSetBuilderTest, BasicVectorDimIsInnerDimension)
{
    constexpr int                                      numFrames = 5;
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    std::vector<hsize_t> dims(2, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames) << "Outer dimension should be numFrames for row major order";
    EXPECT_EQ(dims[1], DIM) << "Inner dimension should be DIM for row major order";
}

TYPED_TEST(H5mdBasicVectorFrameDataSetBuilderTest, SetFrameDimensionAddsDimensionInCenter1)
{
    constexpr int                                      numFrames      = 3;
    const std::vector<hsize_t>                         frameDimension = { 5 };
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .withFrameDimension(frameDimension)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    ASSERT_EQ(H5Sget_simple_extent_ndims(dataSpace), 3)
            << "Adding one dimension should result in 3d";
    std::vector<hsize_t> dims(3, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames) << "Outer dimension should be numFrames for row major order";
    EXPECT_EQ(dims[1], frameDimension[0])
            << "Frame dimension should be in center for row major order";
    EXPECT_EQ(dims[2], DIM) << "Inner dimension should be DIM for row major order";
}

TYPED_TEST(H5mdBasicVectorFrameDataSetBuilderTest, SetFrameDimensionAddsDimensionInCenter2)
{
    constexpr int                                      numFrames      = 3;
    const std::vector<hsize_t>                         frameDimension = { 5, 7 };
    const H5mdDataSetBase<gmx::BasicVector<TypeParam>> dataSet =
            H5mdFrameDataSetBuilder<gmx::BasicVector<TypeParam>>(this->fileid(), "testDataSet")
                    .withNumFrames(numFrames)
                    .withFrameDimension(frameDimension)
                    .build();
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));

    ASSERT_EQ(H5Sget_simple_extent_ndims(dataSpace), 4)
            << "Adding two dimensions should result in 4d";
    std::vector<hsize_t> dims(4, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    EXPECT_EQ(dims[0], numFrames) << "Outer dimension should be numFrames for row major order";
    EXPECT_EQ(dims[1], frameDimension[0])
            << "Frame dimension 0 should be in center for row major order";
    EXPECT_EQ(dims[2], frameDimension[1])
            << "Frame dimension 1 should be in center for row major order";
    EXPECT_EQ(dims[3], DIM) << "Inner dimension should be DIM for row major order";
}

TYPED_TEST(H5mdFrameDataSetBuilderTest, UnitAttributeNotSetByDefault)
{
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();

    EXPECT_FALSE(getAttribute<std::string>(dataSet.id(), "unit").has_value());
}

TYPED_TEST(H5mdFrameDataSetBuilderTest, WithUnitAttribute)
{
    constexpr char                   unit[] = "cm+2 s-1";
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withUnit("unused unit") // ensure that only the last .withUnit() value is used
                    .withUnit(unit)
                    .build();

    const std::optional<std::string> unitAttribute = getAttribute<std::string>(dataSet.id(), "unit");
    ASSERT_TRUE(unitAttribute.has_value()) << "Unit attribute was not set";
    EXPECT_EQ(unitAttribute.value(), unit) << "Incorrect unit attribute";
}

} // namespace
} // namespace test
} // namespace gmx
