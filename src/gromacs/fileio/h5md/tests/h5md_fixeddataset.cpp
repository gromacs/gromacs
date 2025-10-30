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

#include <array>
#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
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

//! \brief Helper function to read fixed-size strings from a data set.
//
// Added for testing. Remove once read methods are properly implemented
// in \c H5mdFixedDataSet.
std::vector<std::string> readFixedSizeStringData(const hid_t container, const char* name)
{
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(H5Dopen(container, name, H5P_DEFAULT));
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet));

    // Get the maximum string size (including the terminating '\0')
    const size_t maxStringSize = H5Tget_size(dataType);

    // Get the total number of values in data set
    const auto [dataSpace, dataSpaceGuard] = makeH5mdDataSpaceGuard(H5Dget_space(dataSet));
    const int            numDims           = H5Sget_simple_extent_ndims(dataSpace);
    std::vector<hsize_t> dims(numDims, 0);
    H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    int numValues = 1;
    for (const hsize_t d : dims)
    {
        numValues *= d;
    }

    std::vector<char> readBuffer(numValues * maxStringSize);

    throwUponH5mdError(H5Dread(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, readBuffer.data()) < 0,
                       "Error writing data.");

    std::vector<std::string> stringValues;
    stringValues.reserve(numValues);
    auto startOfNextString = readBuffer.cbegin();
    for (int i = 0; i < numValues; ++i)
    {
        stringValues.emplace_back(startOfNextString,
                                  std::find(startOfNextString, startOfNextString + maxStringSize, '\0'));
        startOfNextString += maxStringSize;
    }
    return stringValues;
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

TEST_P(WithDims, WriteDataWorksForFixedSizeStrings)
{
    const DataSetDims dims          = GetParam().dims_;
    constexpr int     maxStringSize = 3; // without terminating '\0'

    H5mdFixedDataSet<std::string> dataSet =
            H5mdDataSetBuilder<std::string>(fileid(), "testDataSet")
                    .withMaxStringLength(maxStringSize + 1) // with terminating '\0'
                    .withDimension(dims)
                    .build();

    std::vector<std::string> stringsToWrite;
    stringsToWrite.reserve(dataSet.numValues());
    for (hsize_t i = 0; i < dataSet.numValues(); ++i)
    {
        std::string uniqueString;
        uniqueString.resize(maxStringSize);

        // Fill each string with unique values, indexing within 'a'..='z' repeating
        std::generate(uniqueString.begin(),
                      uniqueString.end(),
                      []()
                      {
                          constexpr char offset = 'a';
                          constexpr char range  = 'z' - 'a';
                          static int     j      = 0;
                          return offset + (j++ % range);
                      });
        stringsToWrite.push_back(uniqueString);
    }

    dataSet.writeData(stringsToWrite);

    // TODO: Once fixed-size string reading is implemented this test should
    // be updated to use it.
    const std::vector<std::string> readStrings = readFixedSizeStringData(fileid(), "testDataSet");
    for (hsize_t i = 0; i < dataSet.numValues(); ++i)
    {
        EXPECT_EQ(readStrings[i], stringsToWrite[i]);
    }
}

TEST_P(WithDims, ReadAndWriteDataWorksForVariableSizeStrings)
{
    const DataSetDims dims = GetParam().dims_;

    H5mdFixedDataSet<std::string> dataSet = H5mdDataSetBuilder<std::string>(fileid(), "testDataSet")
                                                    .withVariableStringLength()
                                                    .withDimension(dims)
                                                    .build();

    // We want this to test "arbitrary" string lengths, so define a table
    // of different lengths to use for each string written to and read from
    // the data set in this test.
    constexpr std::array<int, 8> variableStringLengths = { 5, 1, 2, 13, 0, 3, 8, 5 };

    std::vector<std::string> stringsToWrite;
    stringsToWrite.reserve(dataSet.numValues());
    for (hsize_t i = 0; i < dataSet.numValues(); ++i)
    {
        const int   stringSize = variableStringLengths.at(i % variableStringLengths.size());
        std::string uniqueString;
        uniqueString.resize(stringSize);

        // Fill each string with unique values, indexing within 'a'..='z' repeating
        std::generate(uniqueString.begin(),
                      uniqueString.end(),
                      []()
                      {
                          constexpr char offset = 'a';
                          constexpr char range  = 'z' - 'a';
                          static int     j      = 0;
                          return offset + (j++ % range);
                      });
        stringsToWrite.push_back(uniqueString);
    }

    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readStringsBuffer(dataSet.numValues());
    dataSet.readData(readStringsBuffer);
    for (hsize_t i = 0; i < dataSet.numValues(); ++i)
    {
        EXPECT_EQ(readStringsBuffer[i], stringsToWrite[i]);
    }
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

/******************************************************************************
 * TEST SUITE FOR READING AND WRITING STRINGS                                 *
 *                                                                            *
 * This suite is parametrized to run over both variable and fixed-size string *
 * data sets, ensuring that both kinds pass string-specific tests.            *
 ******************************************************************************/

// Maximum length for variable-size strings in these tests: used as the upper bound
// when checking expected results, which requires a .substr() call when testing
// fixed-size strings. Using this as an upper limit for variable strings simplifies
// the tests.
constexpr size_t c_maxVariableStringLength = 4096;

//! \brief Type of string data set with fixed max string size, or nullopt for variable-string data sets.
struct StringDataSetType
{
    //!< Fixed-size of string (includes terminating '\0')
    std::optional<size_t> maxStringSize;
};

//! \brief Test fixture for an open H5md file with string data sets
class StringTypes : public H5mdTestBase, public ::testing::WithParamInterface<StringDataSetType>
{
public:
    //! \brief Create and return a data set for the current string test parameter.
    H5mdFixedDataSet<std::string> createDataSet(const DataSetDims& dims)
    {
        H5mdDataSetBuilder<std::string> builder = H5mdDataSetBuilder<std::string>(this->fileid(), name_);
        builder.withDimension(dims);

        isFixedStringDataSet_ = GetParam().maxStringSize.has_value();
        if (isFixedStringDataSet_)
        {
            const size_t maxStringSize = GetParam().maxStringSize.value();
            builder.withMaxStringLength(maxStringSize);
        }

        return builder.build();
    }

    bool isFixedStringDataSet_;

    const char* name_ = "testDataSet";
};

//! \brief Helper function for GTest to print size parameter.
void PrintTo(const StringDataSetType& param, std::ostream* os)
{
    if (param.maxStringSize.has_value())
    {
        *os << "FixedSize(" << param.maxStringSize.value() << ")";
    }
    else
    {
        *os << "VariableSize";
    }
}

//! \brief Helper function for GTest to construct test names.
static std::string nameOfStringTest(const ::testing::TestParamInfo<StringDataSetType>& info)
{
    if (info.param.maxStringSize.has_value())
    {
        return gmx::formatString("FixedSizeString%lu", info.param.maxStringSize.value());
    }
    else
    {
        return "VariableSizeString";
    }
}

TEST_P(StringTypes, Size1StringsWork)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringSize = GetParam().maxStringSize.value_or(c_maxVariableStringLength);

    const DataSetDims             dims    = { 3 };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToWrite(dataSet.numValues(), "A");
    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        // For all strings except fixed-size == 1 we can compare the full string above,
        // if the size is 1 we know that the created string is empty since the written
        // string includes the '\0' terminator
        EXPECT_EQ(readBuffer[i], testMaxStringSize == 1 ? "" : stringsToWrite[i]);
    }
}

TEST_P(StringTypes, EmptyStringsWork)
{
    const DataSetDims             dims    = { 3 };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToWrite(dataSet.numValues(), "");
    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        EXPECT_EQ(readBuffer[i], "");
    }
}

TEST_P(StringTypes, ExactFixedSizeLengthWorks)
{
    // Variable-size strings have no fixed size to test for, so just return
    if (!GetParam().maxStringSize.has_value())
    {
        return;
    }

    // Maximum string size (including terminating '\0')
    const size_t testMaxStringSize = GetParam().maxStringSize.value_or(c_maxVariableStringLength);

    const DataSetDims             dims    = { 3 };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToWrite;
    stringsToWrite.reserve(dataSet.numValues());
    for (hsize_t i = 0; i < dataSet.numValues(); ++i)
    {
        std::string uniqueString;
        uniqueString.resize(testMaxStringSize - 1); // -1 to account for the '\0' terminator

        // Fill each string with unique values, indexing within 'a'..='z' repeating
        std::generate(uniqueString.begin(),
                      uniqueString.end(),
                      []()
                      {
                          constexpr char offset = 'a';
                          constexpr char range  = 'z' - 'a';
                          static int     j      = 0;
                          return offset + (j++ % range);
                      });
        stringsToWrite.push_back(uniqueString);
    }
    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        EXPECT_EQ(readBuffer[i], stringsToWrite[i]);
    }
}

TEST_P(StringTypes, LongStringsAreTrimmedToMaxSize)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringSize = GetParam().maxStringSize.value_or(c_maxVariableStringLength);

    const DataSetDims             dims    = { 3 };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToWrite(
            dataSet.numValues(), "A long string which should be trimmed when written as fixed-size");
    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        // Compare testMaxStringSize - 1 characters to account for the terminating '\0'
        EXPECT_EQ(readBuffer[i], stringsToWrite[i].substr(0, testMaxStringSize - 1));
    }
}

TEST_P(StringTypes, WritingOverwritesOldData)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringSize = GetParam().maxStringSize.value_or(c_maxVariableStringLength);

    const DataSetDims             dims    = { 3 };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    std::vector<std::string> stringsToOverwrite(dataSet.numValues(),
                                                "A long string which to overwrite");
    dataSet.writeData(stringsToOverwrite);

    std::vector<std::string> finalStrings(dataSet.numValues(), "A");
    dataSet.writeData(finalStrings);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(finalStrings); ++i)
    {
        // Compare testMaxStringSize - 1 characters to account for the terminating '\0'
        EXPECT_EQ(readBuffer[i], finalStrings[i].substr(0, testMaxStringSize - 1));
    }
}

TEST_P(StringTypes, DifferentCharacterSets)
{
    // Maximum string size (including terminating '\0')
    const size_t testMaxStringSize = GetParam().maxStringSize.value_or(c_maxVariableStringLength);

    const std::vector<std::string> stringsToWrite = {
        "\ttabs\nand\bback\band\tnew\nlines", // control characters
        "100%&#@//\\and.,counting!~<>",       // mildly unusual symbols
        "Hallå!",                             // UTF-8: Swedish
        "你好",                               // UTF-8: Chinese
        "مرحبا",                              // UTF-8: Arabic
        "Привет",                             // UTF-8: Russian
        "こんにちは",                         // UTF-8: Japanese
        "안녕하세요",                         // UTF-8: Korean
        "😀😃😄😁🔥",                         // UTF-8: Emojis
    };

    const DataSetDims             dims    = { stringsToWrite.size() };
    H5mdFixedDataSet<std::string> dataSet = createDataSet(dims);

    dataSet.writeData(stringsToWrite);

    std::vector<std::string> readBuffer;
    if (isFixedStringDataSet_)
    {
        // TODO: Once fixed-size string reading is implemented this test should
        // be updated to use it.
        readBuffer = readFixedSizeStringData(fileid(), name_);
    }
    else
    {
        readBuffer.resize(dataSet.numValues());
        dataSet.readData(readBuffer);
    }

    for (int i = 0; i < gmx::ssize(stringsToWrite); ++i)
    {
        // Compare testMaxStringSize - 1 characters to account for the terminating '\0'
        EXPECT_EQ(readBuffer[i], stringsToWrite[i].substr(0, testMaxStringSize - 1));
    }
}

const StringDataSetType g_stringDataSetTypes[] = {
    { std::nullopt }, // variable-size string
    // { 0 }, HDF5 cannot have fixed-size 0 strings! Tested explicitly below.
    { 1 }, // fixed-size 1 string (empty, since we need space for '\0')
    { 2 }, // fixed-size 2 string (single character)
    { 16 } // fixed-size 16 string
};

INSTANTIATE_TEST_SUITE_P(H5mdFixedDataSetTest, StringTypes, ::testing::ValuesIn(g_stringDataSetTypes), nameOfStringTest);

TEST_F(StringTypes, ThrowsWhenConstructingFixedSize0String)
{
    // H5mdFixedDataSet constructs from H5mdDataSetBase, so ensure that
    // we cannot construct those with fixed-size 0 (or smaller) strings
    EXPECT_THROW(H5mdDataSetBuilder<std::string>(fileid(), "negativeSize").withMaxStringLength(-1), gmx::FileIOError)
            << "Must throw when constructing data set for fixed-size -1 strings";
    EXPECT_THROW(H5mdDataSetBuilder<std::string>(fileid(), "size0").withMaxStringLength(0), gmx::FileIOError)
            << "Must throw when constructing data set for fixed-size 0 strings";
    EXPECT_NO_THROW(H5mdDataSetBuilder<std::string>(fileid(), "size1").withMaxStringLength(1))
            << "Sanity check: Must not throw for fixed-size 1 strings";
}

} // namespace
} // namespace test
} // namespace gmx
