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
 * \brief Tests for H5MD data set base class.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_datasetbase.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture for all relevant data set types
template<typename ValueType>
class H5mdDataSetBaseTest : public H5mdTestBase
{
};

//! \brief List of all data types to create tests for
using DataTypesToTest =
        ::testing::Types<int32_t, int64_t, float, double, BasicVector<float>, BasicVector<double>>;
TYPED_TEST_SUITE(H5mdDataSetBaseTest, DataTypesToTest);

TYPED_TEST(H5mdDataSetBaseTest, DataTypesAreCorrect)
{
    H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();

    if constexpr (std::is_same_v<TypeParam, BasicVector<float>>)
    {
        EXPECT_TRUE(valueTypeIsDataType<float>(dataSet.storedDataType()));
    }
    else if constexpr (std::is_same_v<TypeParam, BasicVector<double>>)
    {
        EXPECT_TRUE(valueTypeIsDataType<double>(dataSet.storedDataType()));
    }
    else
    {
        EXPECT_TRUE(valueTypeIsDataType<TypeParam>(dataSet.storedDataType()));
    }

    // Since we are creating the data type in this test the native data type
    // will always match the stored data type, so this will always return >0.
    EXPECT_GT(H5Tequal(dataSet.nativeDataType(), dataSet.storedDataType()), 0)
            << "Native data type does not match data type";
}

TYPED_TEST(H5mdDataSetBaseTest, DestructorClosesHandles)
{
    hid_t dataSetHandle;
    hid_t dataTypeHandle;
    hid_t nativeDataTypeHandle;

    {
        H5mdDataSetBase<TypeParam> dataSet =
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
        dataSetHandle        = dataSet.id();
        dataTypeHandle       = dataSet.storedDataType();
        nativeDataTypeHandle = dataSet.nativeDataType();

        ASSERT_TRUE(handleIsValid(dataSetHandle))
                << "Sanity check failed: data set handle must be valid before end of scope";
        ASSERT_TRUE(handleIsValid(dataTypeHandle))
                << "Sanity check failed: data type handle must be valid before end of scope";
        ASSERT_TRUE(handleIsValid(nativeDataTypeHandle))
                << "Sanity check failed: native data type handle must be valid before end of scope";
    }

    EXPECT_FALSE(handleIsValid(dataSetHandle)) << "Data set handle was not closed at end of scope";
    EXPECT_FALSE(handleIsValid(dataTypeHandle))
            << "Data type handle was not closed at end of scope";
    EXPECT_FALSE(handleIsValid(nativeDataTypeHandle))
            << "Native data type handle was not closed at end of scope";
}

TYPED_TEST(H5mdDataSetBaseTest, OpenDataSetWorksForWriteModeFiles)
{
    constexpr char dataSetName[] = "testDataSet";

    {
        const H5mdDataSetBase<TypeParam> dataSet =
                H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), dataSetName).build();
    }
    {
        const H5mdDataSetBase<TypeParam> dataSet(this->fileid(), dataSetName);
        if constexpr (std::is_same_v<TypeParam, BasicVector<float>>)
        {
            EXPECT_TRUE(valueTypeIsDataType<float>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
        else if constexpr (std::is_same_v<TypeParam, BasicVector<double>>)
        {
            EXPECT_TRUE(valueTypeIsDataType<double>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
        else
        {
            EXPECT_TRUE(valueTypeIsDataType<TypeParam>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
    }
}

TYPED_TEST(H5mdDataSetBaseTest, OpenDataSetWorksForReadOnlyFiles)
{
    TestFileManager       fileManager;
    std::filesystem::path fileName = fileManager.getTemporaryFilePath("ref.h5md");

    constexpr char dataSetName[] = "testDataSet";

    {
        H5md                             file(fileName, H5mdFileMode::Write);
        const H5mdDataSetBase<TypeParam> dataSet =
                H5mdFrameDataSetBuilder<TypeParam>(file.fileid(), dataSetName).build();
    }
    {
        H5md                             file(fileName, H5mdFileMode::Read);
        const H5mdDataSetBase<TypeParam> dataSet(file.fileid(), dataSetName);
        if constexpr (std::is_same_v<TypeParam, BasicVector<float>>)
        {
            EXPECT_TRUE(valueTypeIsDataType<float>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
        else if constexpr (std::is_same_v<TypeParam, BasicVector<double>>)
        {
            EXPECT_TRUE(valueTypeIsDataType<double>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
        else
        {
            EXPECT_TRUE(valueTypeIsDataType<TypeParam>(dataSet.storedDataType()))
                    << "Data types must match after opening";
        }
    }
}

//! \brief Test function for actual test below
//
// Open a data set as compiled type and assert that it works if the type matches
// that of the original, and throws if it does not. The data set name is inferred
// from the original type parameter: the actual test sets these up.
//
// \tparam AsType Type to open the data set as
// \tparam OriginalType Original data set type, used to infer the data set name using its typeid
// \param[in] container Parent container in which the data set is created
template<typename AsType, typename OriginalType>
void testOpenDataSetAsType(const hid_t container)
{
    if constexpr (std::is_same_v<AsType, OriginalType>)
    {
        EXPECT_NO_THROW_GMX(H5mdDataSetBase<AsType>(container, typeid(OriginalType).name()))
                << "Sanity check failed: could not open data set of same type";
    }
    else
    {
        EXPECT_THROW_GMX(H5mdDataSetBase<AsType>(container, typeid(OriginalType).name()), H5mdError)
                << formatString(
                           "Must throw if compiled data set type does not match original type %s",
                           typeid(OriginalType).name());
    }
}

TYPED_TEST(H5mdDataSetBaseTest, OpenDataSetThrowsIfOriginalTypeDoesNotMatchTemplate)
{
    if constexpr (std::is_same_v<TypeParam, BasicVector<float>>
                  || std::is_same_v<TypeParam, BasicVector<double>>)
    {
        // Non-primitive types require more advanced testing: their base types should
        // be matched against the native type, and the data set dimensions checked.
        return;
    }

    H5mdFrameDataSetBuilder<int32_t>(this->fileid(), typeid(int32_t).name()).build();
    H5mdFrameDataSetBuilder<int64_t>(this->fileid(), typeid(int64_t).name()).build();
    H5mdFrameDataSetBuilder<float>(this->fileid(), typeid(float).name()).build();
    H5mdFrameDataSetBuilder<double>(this->fileid(), typeid(double).name()).build();

    testOpenDataSetAsType<TypeParam, int32_t>(this->fileid());
    testOpenDataSetAsType<TypeParam, int64_t>(this->fileid());
    testOpenDataSetAsType<TypeParam, float>(this->fileid());
    testOpenDataSetAsType<TypeParam, double>(this->fileid());
}

TYPED_TEST(H5mdDataSetBaseTest, OpenDataSetThrowsForInvalidContainer)
{
    EXPECT_THROW_GMX(H5mdDataSetBase<TypeParam>(H5I_INVALID_HID, "testDataSet"), H5mdError);
}

TYPED_TEST(H5mdDataSetBaseTest, OpenDataSetThrowsForInvalidSetName)
{
    {
        // Create a data set to ensure that there is something in the file which we cannot read with bad names
        H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet").build();
    }

    EXPECT_THROW_GMX(H5mdDataSetBase<TypeParam>(this->fileid(), ""), H5mdError)
            << "Should throw for empty name";
    EXPECT_THROW_GMX(H5mdDataSetBase<TypeParam>(this->fileid(), "aBadIdea"), H5mdError)
            << "Should throw for bad name";
}

TYPED_TEST(H5mdDataSetBaseTest, DimsReturnsDataSetDimensions)
{
    const std::vector<hsize_t>       frameDimensions = { 5, 2 };
    const H5mdDataSetBase<TypeParam> dataSet =
            H5mdFrameDataSetBuilder<TypeParam>(this->fileid(), "testDataSet")
                    .withFrameDimension(frameDimensions)
                    .build();

    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    DataSetDims expectedDims(H5Sget_simple_extent_ndims(dataSpace), 0);
    H5Sget_simple_extent_dims(dataSpace, expectedDims.data(), nullptr);
    EXPECT_EQ(dataSet.dims(), expectedDims) << "dims() does not return the data set dimensions";
}

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

    *os << "Dims {" << formatAndJoin(info.dims_, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfTest(const ::testing::TestParamInfo<TestDimensions>& info)
{
    const auto  toString = [](const hsize_t value) -> std::string { return std::to_string(value); };
    std::string testName = formatAndJoin(info.param.dims_, "_", toString);

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
class H5mdDataSetBaseBasicVectorTest : public H5mdTestBase, public ::testing::WithParamInterface<TestDimensions>
{
};

const TestDimensions g_dataSetDimsToTest[] = { { { 3 } }, { { 0, 3 } }, { { 0, 0, 3 } }, { { 5, 2, 3 } } };

TEST_P(H5mdDataSetBaseBasicVectorTest, A)
{
    DataSetDims dims = GetParam().dims_;

    {
        SCOPED_TRACE("Construct float and double type data sets with inner dimensions 2, 3 and 4");

        dims.back() = 2;
        H5mdDataSetBuilder<float>(this->fileid(), "float_2").withDimension(dims).build();
        H5mdDataSetBuilder<double>(this->fileid(), "double_2").withDimension(dims).build();
        dims.back() = 3;
        H5mdDataSetBuilder<float>(this->fileid(), "float_3").withDimension(dims).build();
        H5mdDataSetBuilder<double>(this->fileid(), "double_3").withDimension(dims).build();
        dims.back() = 4;
        H5mdDataSetBuilder<float>(this->fileid(), "float_4").withDimension(dims).build();
        H5mdDataSetBuilder<double>(this->fileid(), "double_4").withDimension(dims).build();
    }
    {
        SCOPED_TRACE("Open the data sets as BasicVector<float/double>");
        ASSERT_NO_THROW_GMX(H5mdDataSetBase<BasicVector<float>>(this->fileid(), "float_3"))
                << "Sanity check failed: should not throw for inner dimension = 3";
        ASSERT_NO_THROW_GMX(H5mdDataSetBase<BasicVector<double>>(this->fileid(), "double_3"))
                << "Sanity check failed: should not throw for inner dimension = 3";
        EXPECT_THROW_GMX(H5mdDataSetBase<BasicVector<float>>(this->fileid(), "float_2"), H5mdError)
                << "Must throw for inner dimension = 2 (BasicVector<float>)";
        EXPECT_THROW_GMX(H5mdDataSetBase<BasicVector<double>>(this->fileid(), "double_2"), H5mdError)
                << "Must throw for inner dimension = 2 (BasicVector<double>)";
        EXPECT_THROW_GMX(H5mdDataSetBase<BasicVector<float>>(this->fileid(), "float_4"), H5mdError)
                << "Must throw for inner dimension = 4 (BasicVector<float>)";
        EXPECT_THROW_GMX(H5mdDataSetBase<BasicVector<double>>(this->fileid(), "double_4"), H5mdError)
                << "Must throw for inner dimension = 4 (BasicVector<double>)";
    }
}

INSTANTIATE_TEST_SUITE_P(H5mdDataSetBuilderValidDims,
                         H5mdDataSetBaseBasicVectorTest,
                         ::testing::ValuesIn(g_dataSetDimsToTest),
                         nameOfTest);

/******************************************************************************
 * TEST SUITE FOR READING AND WRITING STRINGS                                 *
 *                                                                            *
 * This suite is parametrized to run over both variable and fixed-size string *
 * data sets, ensuring that both kinds pass string-specific tests.            *
 ******************************************************************************/

// Maximum length for variable-size strings in these tests
static constexpr size_t sc_maxVariableStringLength = 4096;

//! \brief Type of string data set with fixed max string size, or nullopt for variable-string data sets.
struct StringDataSetType
{
    //!< Fixed-size of string (includes terminating '\0')
    std::optional<size_t> maxStringLength;
};

//! \brief Test fixture for an open H5md file with string data sets
class H5mdDataSetBaseStringTest : public H5mdTestBase, public ::testing::WithParamInterface<StringDataSetType>
{
public:
    //! \brief Create and return a data set for the current string test parameter.
    H5mdDataSetBase<std::string> createDataSet(const DataSetDims& dims)
    {
        H5mdDataSetBuilder<std::string> builder = H5mdDataSetBuilder<std::string>(this->fileid(), name_);
        builder.withDimension(dims);

        isFixedStringDataSet_ = GetParam().maxStringLength.has_value();
        if (isFixedStringDataSet_)
        {
            const size_t maxStringLength = GetParam().maxStringLength.value();
            builder.withMaxStringLength(maxStringLength);
        }

        return builder.build();
    }

    bool isFixedStringDataSet_;

    const char* name_ = "testDataSet";
};

//! \brief Helper function for GTest to print size parameter.
void PrintTo(const StringDataSetType& param, std::ostream* os)
{
    if (param.maxStringLength.has_value())
    {
        *os << "FixedLength(" << param.maxStringLength.value() << ")";
    }
    else
    {
        *os << "VariableLength";
    }
}

//! \brief Helper function for GTest to construct test names.
static std::string nameOfStringTest(const ::testing::TestParamInfo<StringDataSetType>& info)
{
    if (info.param.maxStringLength.has_value())
    {
        return formatString("FixedLengthString%lu", info.param.maxStringLength.value());
    }
    else
    {
        return "VariableLengthString";
    }
}

TEST_P(H5mdDataSetBaseStringTest, FullDataSetIo)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringLength = GetParam().maxStringLength.value_or(sc_maxVariableStringLength);

    constexpr int                nx = 3, ny = 2, nz = 4;
    const DataSetDims            dims    = { nx, ny, nz };
    H5mdDataSetBase<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToWrite;
    stringsToWrite.reserve(nx * ny * nz);
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                stringsToWrite.push_back(formatString("%d.%d.%d", i, j, k));
            }
        }
    }

    ASSERT_TRUE(dataSet.write(stringsToWrite, H5S_ALL, H5S_ALL)) << "Could not write string data";

    std::vector<std::string> readBuffer(stringsToWrite.size());
    ASSERT_TRUE(dataSet.read(readBuffer, H5S_ALL, H5S_ALL)) << "Could not read string data";

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        // Compare testMaxStringLength - 1 characters to account for the terminating '\0'
        ASSERT_EQ(readBuffer[i], stringsToWrite[i].substr(0, testMaxStringLength - 1));
    }
}

TEST_P(H5mdDataSetBaseStringTest, RandomAccessRead)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringLength = GetParam().maxStringLength.value_or(sc_maxVariableStringLength);

    constexpr int                nx = 3, ny = 2, nz = 4;
    const DataSetDims            dims    = { nx, ny, nz };
    H5mdDataSetBase<std::string> dataSet = createDataSet(dims);

    // Write easily identifiable strings to each point in the data set
    std::vector<std::string> stringsToWrite;
    stringsToWrite.reserve(nx * ny * nz);
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                stringsToWrite.push_back(formatString("%d.%d.%d", i, j, k));
            }
        }
    }
    dataSet.write(stringsToWrite, H5S_ALL, H5S_ALL);

    // Define a set of points to read
    constexpr hsize_t                             numPoints = 4;
    std::array<std::array<hsize_t, 3>, numPoints> coords    = {
        { { 0, 0, 0 }, { 2, 1, 3 }, { 1, 1, 1 }, { 1, 1, 2 } }
    };
    const auto [memoryDataSpace, memoryDataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, &numPoints, nullptr));
    // Define in-file selection for points
    const auto [fileDataSpace, fileDataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    H5Sselect_elements(
            fileDataSpace, H5S_SELECT_SET, numPoints, reinterpret_cast<const hsize_t*>(coords.data()));

    std::vector<std::string> readBuffer(numPoints);
    ASSERT_TRUE(dataSet.read(readBuffer, memoryDataSpace, fileDataSpace))
            << "Could not read string data from points";
    // Compare testMaxStringLength - 1 characters to account for the terminating '\0'
    EXPECT_EQ(readBuffer[0], stringsToWrite[0].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[1], stringsToWrite[23].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[2], stringsToWrite[13].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[3], stringsToWrite[14].substr(0, testMaxStringLength - 1));
}

TEST_P(H5mdDataSetBaseStringTest, RandomAccessWrite)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringLength = GetParam().maxStringLength.value_or(sc_maxVariableStringLength);

    const DataSetDims            dims    = { 5, 4 };
    H5mdDataSetBase<std::string> dataSet = createDataSet(dims);

    constexpr hsize_t                             numPoints = 4;
    std::array<std::array<hsize_t, 2>, numPoints> coords    = {
        { { 0, 0 }, { 4, 3 }, { 2, 1 }, { 2, 2 } },
    };
    std::array<std::string, numPoints> stringsToWrite = { "0.0", "4.3", "2.1", "2.2" };
    const auto [memoryDataSpace, memoryDataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, &numPoints, nullptr));

    // Define in-file selection for points
    const auto [fileDataSpace, fileDataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet.id()));
    H5Sselect_elements(
            fileDataSpace, H5S_SELECT_SET, numPoints, reinterpret_cast<const hsize_t*>(coords.data()));

    ASSERT_TRUE(dataSet.write(stringsToWrite, memoryDataSpace, fileDataSpace))
            << "Could not write string data to points";
    std::vector<std::string> readBuffer(numPoints);
    dataSet.read(readBuffer, memoryDataSpace, fileDataSpace);
    // Compare testMaxStringLength - 1 characters to account for the terminating '\0'
    EXPECT_EQ(readBuffer[0], stringsToWrite[0].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[1], stringsToWrite[1].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[2], stringsToWrite[2].substr(0, testMaxStringLength - 1));
    EXPECT_EQ(readBuffer[3], stringsToWrite[3].substr(0, testMaxStringLength - 1));
}

const StringDataSetType g_stringDataSetTypes[] = {
    { std::nullopt }, // variable-size string
    // { 0 }, HDF5 cannot have fixed-size 0 strings!
    { 1 }, // fixed-size 1 string (empty, since we need space for '\0')
    { 2 }, // fixed-size 2 string (single character)
    { 16 } // fixed-size 16 string
};
INSTANTIATE_TEST_SUITE_P(H5mdDataSetBaseTest,
                         H5mdDataSetBaseStringTest,
                         ::testing::ValuesIn(g_stringDataSetTypes),
                         nameOfStringTest);

} // namespace
} // namespace test
} // namespace gmx
