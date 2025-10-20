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
 * \brief Tests for H5MD fixed data set class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_fixeddataset.h"

#include <hdf5.h>

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Helper struct to parametrize tests for combinations of dimensions, max dimensions and chunk dimensions.
struct TestDimensions
{
    TestDimensions(const std::vector<hsize_t>& dims) : dims_{ dims } {}

    //!< Data set dimensions.
    std::vector<hsize_t> dims_;
};

//! \brief Helper function for GTest to print dimension parameters.
void PrintTo(const TestDimensions& info, std::ostream* os)
{
    const auto toString = [&](hsize_t value) -> std::string { return std::to_string(value); };

    *os << "Dims {" << gmx::formatAndJoin(info.dims_, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfTest(const ::testing::TestParamInfo<TestDimensions>& info)
{
    const auto  toString = [](const hsize_t value) -> std::string { return std::to_string(value); };
    std::string testName = gmx::formatAndJoin(info.param.dims_, "_", toString);

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

//! \brief Parametrized test fixture for tests with different data set dimensions for BasicVector<T> sets.
class WithDims : public H5mdTestBase, public ::testing::WithParamInterface<TestDimensions>
{
};

const TestDimensions g_dataSetDimsToTest[] = { { { 0 } }, // edge case: 0 values in data set
                                               { { 0, 0 } },
                                               { { 3 } },
                                               { { 0, 3 } },
                                               { { 5, 2, 3 } } };

TEST_P(WithDims, DimsIsConsistent)
{
    const DataSetDims         dims = GetParam().dims_;
    H5mdFixedDataSet<int32_t> dataSet =
            H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension(dims).build();

    EXPECT_EQ(dataSet.dims(), dims);
}

TEST_P(WithDims, NumValuesIsConsistent)
{
    const DataSetDims         dims = GetParam().dims_;
    H5mdFixedDataSet<int32_t> dataSet =
            H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension(dims).build();

    size_t numValues = 1;
    for (const auto d : dims)
    {
        numValues *= d;
    }

    EXPECT_EQ(dataSet.numValues(), numValues);
}

TEST_P(WithDims, ReadAndWriteDataWorksForPrimitives)
{
    const DataSetDims         dims = GetParam().dims_;
    H5mdFixedDataSet<int32_t> dataSet =
            H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension(dims).build();

    std::vector<int32_t> valuesToWrite(dataSet.numValues());
    std::iota(valuesToWrite.begin(), valuesToWrite.end(), 0);

    dataSet.writeData(valuesToWrite);
    std::vector<int32_t> readBuffer(valuesToWrite.size());
    dataSet.readData(readBuffer);
    EXPECT_EQ(readBuffer, valuesToWrite);
}

TEST_P(WithDims, ReadAndWriteDataWorksForBasicVector)
{
    const DataSetDims                    dims = GetParam().dims_;
    H5mdFixedDataSet<BasicVector<float>> dataSet =
            H5mdDataSetBuilder<BasicVector<float>>(fileid(), "testDataSet").withDimension(dims).build();

    EXPECT_EQ(dataSet.dims(), dims);

    std::vector<BasicVector<float>> valuesToWrite(dataSet.numValues());
    std::generate(valuesToWrite.begin(),
                  valuesToWrite.end(),
                  []()
                  {
                      static float             value = 0.0;
                      const BasicVector<float> ret   = { value, value + 0.1f, value + 0.01f };
                      ++value;
                      return ret;
                  });

    dataSet.writeData(valuesToWrite);
    std::vector<BasicVector<float>> readBuffer(valuesToWrite.size());
    dataSet.readData(readBuffer);
    EXPECT_EQ(readBuffer, valuesToWrite);
}

TEST_P(WithDims, OpenDataSetIsConsistent)
{
    const DataSetDims dims = GetParam().dims_;
    {
        SCOPED_TRACE("Create a data set and close it at the end of scope");
        H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension(dims).build();
    }

    H5mdFixedDataSet<int32_t> dataSet(fileid(), "testDataSet");
    EXPECT_EQ(dataSet.dims(), dims);

    std::vector<int32_t> valuesToWrite(dataSet.numValues());
    std::iota(valuesToWrite.begin(), valuesToWrite.end(), 0);

    dataSet.writeData(valuesToWrite);
    std::vector<int32_t> readBuffer(valuesToWrite.size());
    dataSet.readData(readBuffer);
    EXPECT_EQ(readBuffer, valuesToWrite);
}

TEST_P(WithDims, OpenDataSetThrowsForWrongType)
{
    const DataSetDims dims = GetParam().dims_;
    {
        SCOPED_TRACE("Create a data set and close it at the end of scope");
        H5mdDataSetBuilder<int32_t>(fileid(), "int32").withDimension(dims).build();
    }

    ASSERT_NO_THROW(H5mdFixedDataSet<int32_t>(fileid(), "int32"))
            << "Sanity check failed: Open with same type must work";
    EXPECT_THROW(H5mdFixedDataSet<float>(fileid(), "int32"), gmx::FileIOError)
            << "Must throw when opening as float";
}

TEST_P(WithDims, ThrowsForReadAndWriteWithWrongSizeBuffers)
{
    const DataSetDims         dims = GetParam().dims_;
    H5mdFixedDataSet<int32_t> dataSet =
            H5mdDataSetBuilder<int32_t>(fileid(), "testDataSet").withDimension(dims).build();


    if (dataSet.numValues() > 0)
    {
        std::vector<int32_t> bufferTooSmall(dataSet.numValues() - 1);
        EXPECT_THROW(dataSet.writeData(bufferTooSmall), gmx::FileIOError);
        EXPECT_THROW(dataSet.readData(bufferTooSmall), gmx::FileIOError);
    }

    std::vector<int32_t> bufferTooLarge(dataSet.numValues() + 1);
    EXPECT_THROW(dataSet.writeData(bufferTooLarge), gmx::FileIOError);
    EXPECT_THROW(dataSet.readData(bufferTooLarge), gmx::FileIOError);
}

INSTANTIATE_TEST_SUITE_P(H5mdFixedDataSetTest, WithDims, ::testing::ValuesIn(g_dataSetDimsToTest), nameOfTest);

} // namespace
} // namespace test
} // namespace gmx
