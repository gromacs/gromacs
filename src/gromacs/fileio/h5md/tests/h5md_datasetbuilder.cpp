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
 * \brief Tests for H5MD data set builder routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"

#include <hdf5.h>

#include <iostream>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Helper struct to parametrize tests for combinations of dimensions, max dimensions and chunk dimensions.
struct DataSetDimensions
{
    DataSetDimensions(const std::vector<hsize_t>& dims) : dims_{ dims }, expectedResult_{ true } {}

    DataSetDimensions(const std::vector<hsize_t>& dims,
                      const std::vector<hsize_t>& maxDims,
                      const std::vector<hsize_t>& chunkDims,
                      const bool                  expectedResult) :
        dims_{ dims }, maxDims_{ maxDims }, chunkDims_{ chunkDims }, expectedResult_{ expectedResult }
    {
    }

    //!< Per-dimension size of data set.
    std::vector<hsize_t> dims_;

    //!< Per-dimension maximum size of data set.
    std::vector<hsize_t> maxDims_;

    //!< Chunk dimensions for data set.
    std::vector<hsize_t> chunkDims_;

    //!< Whether the above combination of dimension values is valid.
    bool expectedResult_;
};

//! \brief Helper function for GTest to print dimension parameters.
void PrintTo(const DataSetDimensions& info, std::ostream* os)
{
    const auto toString = [&](hsize_t value) -> std::string
    {
        if (value == H5S_UNLIMITED)
        {
            return "Inf";
        }
        else
        {
            return std::to_string(value);
        }
    };

    *os << "Dims {" << gmx::formatAndJoin(info.dims_, ", ", toString) << "}" << ", "
        << "MaxDims {" << gmx::formatAndJoin(info.maxDims_, ", ", toString) << "}" << ", "
        << "ChunkDims {" << gmx::formatAndJoin(info.chunkDims_, ", ", toString) << "}" << ", "
        << "ExpectedResult = " << (info.expectedResult_ ? "true" : "false");
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfTest(const ::testing::TestParamInfo<DataSetDimensions>& info)
{
    std::vector<std::string> partNames;

    const auto toString = [](const hsize_t value) -> std::string
    {
        if (value == H5S_UNLIMITED)
        {
            return "Inf";
        }
        else
        {
            return std::to_string(value);
        }
    };

    if (!info.param.dims_.empty())
    {
        std::string partName;
        partName.append("Dims");
        for (const int d : info.param.dims_)
        {
            partName.append(toString(d));
        }
        partNames.push_back(partName);
    }
    if (!info.param.maxDims_.empty())
    {
        std::string partName;
        partName.append("MaxDims");
        for (const int d : info.param.maxDims_)
        {
            partName.append(toString(d));
        }
        partNames.push_back(partName);
    }
    if (!info.param.chunkDims_.empty())
    {
        std::string partName;
        partName.append("chunkDims");
        for (const int d : info.param.chunkDims_)
        {
            partName.append(toString(d));
        }
        partNames.push_back(partName);
    }

    std::string testName = gmx::joinStrings(partNames.cbegin(), partNames.cend(), "_");
    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");

    if (testName.empty())
    {
        return "Empty";
    }
    else
    {
        return testName;
    }
}

//! Test fixture for data set builders, inheriting parametrization capabilities.
template<typename ValueType>
class H5mdDataSetBuilderTestBase : public H5mdTestBase, public ::testing::WithParamInterface<DataSetDimensions>
{
};

//! Parametrized test suite for data sets of simple primitives
using Primitive = H5mdDataSetBuilderTestBase<int32_t>;

//! Parametrized test suite for data sets of BasicVector<T>
using BasicVector = H5mdDataSetBuilderTestBase<gmx::BasicVector<float>>;

TEST_P(Primitive, Works)
{
    const DataSetDimensions& testParam = GetParam();

    H5mdDataSetBuilder<int32_t> builder(this->fileid(), "testDataSet");
    builder.withDimension(testParam.dims_);
    builder.withMaxDimension(testParam.maxDims_);
    builder.withChunkDimension(testParam.chunkDims_);

    if (testParam.expectedResult_)
    {
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(builder.build());

        // Check in order:
        // data type, data set dimensions, max dimensions and chunk dimensions
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        EXPECT_TRUE(valueTypeIsDataType<int32_t>(dataType)) << "Incorrect data type in data set";

        const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
        EXPECT_EQ(H5Sget_simple_extent_ndims(dataSpace), testParam.dims_.size())
                << "Incorrect number of dimensions";
        std::vector<hsize_t> dims(testParam.dims_.size(), 0);
        std::vector<hsize_t> maxDims(testParam.dims_.size(), 0);
        H5Sget_simple_extent_dims(dataSpace, dims.data(), maxDims.data());
        EXPECT_EQ(dims, testParam.dims_) << "Incorrect dimensions";

        if (testParam.maxDims_.empty())
        {
            EXPECT_EQ(maxDims, testParam.dims_)
                    << "Incorrect max dimensions: should be set to dims if not specified";
        }
        else
        {
            EXPECT_EQ(maxDims, testParam.maxDims_) << "Incorrect max dimensions";
        }

        const auto [propertyList, propertyListGuard] =
                makeH5mdPropertyListGuard(H5Dget_create_plist(dataSet));
        std::vector<hsize_t> chunkDims(testParam.dims_.size(), 0);
        H5Pget_chunk(propertyList, testParam.dims_.size(), chunkDims.data());
        if (testParam.chunkDims_.empty())
        {
            std::vector<hsize_t> expectedChunkDims = testParam.dims_;
            for (hsize_t& d : expectedChunkDims)
            {
                if (d == 0)
                {
                    d = 1;
                }
            }
            EXPECT_EQ(chunkDims, expectedChunkDims) << "Incorrect chunk dimensions: should be set "
                                                       "to dims (min. 1) if not specified";
        }
        else
        {
            EXPECT_EQ(chunkDims, testParam.chunkDims_) << "Incorrect chunk dimensions";
        }
    }
    else
    {
        EXPECT_THROW(builder.build(), gmx::FileIOError);
    }
}

TEST_P(BasicVector, Works)
{
    const DataSetDimensions& testParam = GetParam();

    // For a data set of type BasicVector<T> an extra dimension is added as an inner
    // row to contain DIM extra values of type T. This test is identical to that for
    // primitive types but the expected dims is the input testParam.dims_ + [3],
    // and so on for maxDims_ and chunkDims_.

    H5mdDataSetBuilder<gmx::BasicVector<float>> builder(this->fileid(), "testDataSet");
    builder.withDimension(testParam.dims_);
    builder.withMaxDimension(testParam.maxDims_);
    builder.withChunkDimension(testParam.chunkDims_);

    if (testParam.expectedResult_)
    {
        const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(builder.build());

        // Check in order:
        // data type, data set dimensions, max dimensions and chunk dimensions
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));
        EXPECT_TRUE(valueTypeIsDataType<float>(dataType)) << "Incorrect data type in data set";

        const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
        // Add the extra vector dimension here!
        const size_t expectedNumDims = testParam.dims_.size() + 1;
        EXPECT_EQ(H5Sget_simple_extent_ndims(dataSpace), expectedNumDims)
                << "Incorrect number of dimensions: should be input dims + 1 for vector";

        // And here!
        std::vector<hsize_t> dims(expectedNumDims, 0);
        std::vector<hsize_t> maxDims(expectedNumDims, 0);
        H5Sget_simple_extent_dims(dataSpace, dims.data(), maxDims.data());

        // And here!
        std::vector<hsize_t> expectedDims = testParam.dims_;
        expectedDims.push_back(DIM);
        EXPECT_EQ(dims, expectedDims) << "Incorrect dimensions";

        if (testParam.maxDims_.empty())
        {
            EXPECT_EQ(maxDims, expectedDims)
                    << "Incorrect max dimensions: should be set to dims if not specified";
        }
        else
        {
            // And here!
            std::vector<hsize_t> expectedMaxDims = testParam.maxDims_;
            expectedMaxDims.push_back(DIM);
            EXPECT_EQ(maxDims, expectedMaxDims) << "Incorrect max dimensions";
        }

        const auto [propertyList, propertyListGuard] =
                makeH5mdPropertyListGuard(H5Dget_create_plist(dataSet));
        // And here!
        std::vector<hsize_t> chunkDims(expectedNumDims, 0);
        H5Pget_chunk(propertyList, chunkDims.size(), chunkDims.data());
        if (testParam.chunkDims_.empty())
        {
            std::vector<hsize_t> expectedChunkDims = expectedDims;
            for (hsize_t& d : expectedChunkDims)
            {
                if (d == 0)
                {
                    d = 1;
                }
            }
            EXPECT_EQ(chunkDims, expectedChunkDims) << "Incorrect chunk dimensions: should be set "
                                                       "to dims (min. 1) if not specified";
        }
        else
        {
            // And here!
            std::vector<hsize_t> expectedChunkDims = testParam.chunkDims_;
            expectedChunkDims.push_back(DIM);
            EXPECT_EQ(chunkDims, expectedChunkDims) << "Incorrect chunk dimensions";
        }
    }
    else
    {
        EXPECT_THROW(builder.build(), gmx::FileIOError);
    }
}

//! \brief Data set dimensions which should work when constructing data sets.
const DataSetDimensions g_validDataSetDimensions[] = {
    // Size = 0 along all dimensions
    { { 0 } },
    { { 0, 0 } },
    { { 0, 0, 0 } },
    // Size = 1 along all dimensions
    { { 1 } },
    { { 1, 1 } },
    { { 1, 1, 1 } },
    // Mixed sizes along all dimensions
    { { 5 } },
    { { 5, 1, 2 } },
    { { 0, 5, 2 } },
    { { 0, 1, 5 } }
};

//! \brief Combinations which should work when constructing data sets.
const DataSetDimensions g_validDimensionCombinations[] = {
    // Empty maxDims and chunkDims are valid
    { { 0 }, {}, {}, true },
    { { 1 }, {}, {}, true },
    // maxDims or chunkDims same size as dims are valid
    { { 0 }, { 0 }, {}, true },
    { { 0 }, {}, { 1 }, true },
    { { 0 }, { 0 }, { 1 }, true },
    { { 0, 0 }, { 0, 0 }, { 1, 1 }, true },
    { { 1, 2 }, { 1, 2 }, {}, true },
    { { 1, 2 }, {}, { 1, 2 }, true },
    { { 1, 2 }, { 1, 2 }, { 1, 2 }, true },
    // maxDims can be larger than dims
    { { 0, 1 }, { 1, 1 }, {}, true },
    { { 0, 1 }, { 0, 2 }, {}, true },
    { { 0, 1 }, { 1, 2 }, {}, true },
    // maxDims can be unlimited along any dimension
    { { 0, 1 }, { H5S_UNLIMITED, 1 }, {}, true },
    { { 0, 1 }, { 0, H5S_UNLIMITED }, {}, true },
    { { 0, 1 }, { H5S_UNLIMITED, H5S_UNLIMITED }, {}, true }
};

//! \brief Combinations which should result in a throw when constructing data sets.
const DataSetDimensions g_invalidDimensionCombinations[] = {
    // Zero dimensions is not valid
    { {}, {}, {}, false },
    // Unlimited data set dims along any dimension is not valid
    { { H5S_UNLIMITED }, {}, {}, false },
    { { H5S_UNLIMITED, 0 }, {}, {}, false },
    { { 0, H5S_UNLIMITED }, {}, {}, false },
    // maxDims or chunkDims of different dimensions than dims is not valid (if set)
    { { 0, 0 }, { 0 }, {}, false },
    { { 0, 0 }, { 0, 0, 0 }, {}, false },
    { { 0, 0 }, {}, { 0 }, false },
    { { 0, 0 }, {}, { 0, 0, 0 }, false },
    // maxDims smaller than dims along any dimension is not valid
    { { 1, 1 }, { 0, 1 }, {}, false },
    { { 1, 1 }, { 1, 0 }, {}, false },
    // chunkDims equal to zero along any dimension is not valid
    { { 1, 1 }, {}, { 0, 1 }, false },
    { { 1, 1 }, {}, { 1, 0 }, false },
};

/*! \brief Instantiate parametrized test suites for int32_t data sets.
 *
 * These cover the creation of data sets with a primitive type for a large set
 * of input dimension combinations.
 */
//! {
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderValidDims,
                         Primitive,
                         ::testing::ValuesIn(g_validDataSetDimensions),
                         nameOfTest);
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderValidDimCombinations,
                         Primitive,
                         ::testing::ValuesIn(g_validDimensionCombinations),
                         nameOfTest);
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderInvalidDimCombinations,
                         Primitive,
                         ::testing::ValuesIn(g_invalidDimensionCombinations),
                         nameOfTest);
//! }

/*! \brief Instantiate parametrized test suites for BasicVector<float> data sets.
 *
 * These cover the creation of data sets with a dimensional type for a large set
 * of input dimension combinations.
 */
//! {
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderValidDims,
                         BasicVector,
                         ::testing::ValuesIn(g_validDataSetDimensions),
                         nameOfTest);
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderValidDimCombinations,
                         BasicVector,
                         ::testing::ValuesIn(g_validDimensionCombinations),
                         nameOfTest);
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderInvalidDimCombinations,
                         BasicVector,
                         ::testing::ValuesIn(g_invalidDimensionCombinations),
                         nameOfTest);
//! }

//! \brief Test fixture for data set builder.
using H5mdDataSetBuilderTest = H5mdTestBase;

TEST_F(H5mdDataSetBuilderTest, ThrowsForUnsetDimension)
{
    EXPECT_THROW(H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").build(), gmx::FileIOError);
    EXPECT_NO_THROW(H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension({ 0 }).build());
}

TEST_F(H5mdDataSetBuilderTest, ThrowsForDuplicateNameInGroup)
{
    const std::string sharedName = "testDataSet";
    ASSERT_NO_THROW(H5mdDataSetBuilder<int32_t>(fileid(), sharedName).withDimension({ 0 }).build());
    EXPECT_THROW(H5mdDataSetBuilder<int32_t>(fileid(), sharedName).withDimension({ 0 }).build(),
                 gmx::FileIOError);
}

TEST_F(H5mdDataSetBuilderTest, ThrowsForEmptyName)
{
    EXPECT_THROW(H5mdDataSetBuilder<int32_t>(fileid(), "").withDimension({ 0 }).build(), gmx::FileIOError);
}

TEST_F(H5mdDataSetBuilderTest, ThrowsForInvalidContainer)
{
    EXPECT_THROW(H5mdDataSetBuilder<int32_t>(H5I_INVALID_HID, "testDataSet").withDimension({ 0 }).build(),
                 gmx::FileIOError);
}

} // namespace
} // namespace test
} // namespace gmx
