/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Tests for H5MD compound data set i/o.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_compound_type.h"

#include <hdf5.h>

#include <array>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_framedataset.h"
#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_topologycompound.h"
#include "gromacs/fileio/h5md/h5md_type.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

using H5mdCompoundTypeTest = H5mdTestBase;

// A compound type is recognised by its static constructHdf5DataType(), and a type with
// fixed-length string members that can be stored with a given size by the sized overload
static_assert(valueTypeIsHdf5CompoundType<MoleculeBlock>());
static_assert(!valueTypeIsHdf5CompoundType<int64_t>());
static_assert(!valueTypeIsHdf5CompoundType<std::string>());
static_assert(valueTypeHasSizedStrings<MoleculeBlock>());
// A type gets a units attribute only if it declares sc_units, which MoleculeBlock does not
static_assert(!valueTypeHasUnits<MoleculeBlock>());

//! \brief Create the compound data set \p name of one record with the given members in \p container.
template<typename RecordType>
void createRawDataSet(const hid_t container, const char* name, const hid_t dataType, const RecordType& record)
{
    const hsize_t numValues = 1;
    const auto [dataSpace, dataSpaceGuard] =
            makeH5mdDataSpaceGuard(H5Screate_simple(1, &numValues, nullptr));
    const auto [dataSet, dataSetGuard] = makeH5mdDataSetGuard(
            H5Dcreate(container, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_GE(dataSet, 0);
    ASSERT_GE(H5Dwrite(dataSet, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &record), 0);
}

TEST_F(H5mdCompoundTypeTest, FixedDataSetWorks)
{
    const DataSetDims dims      = { 5, 3 };
    const int         numValues = std::reduce(dims.cbegin(), dims.cend(), 1, std::multiplies<>());
    H5mdFixedDataSet<MoleculeBlock> dataSet =
            H5mdDataSetBuilder<MoleculeBlock>(fileid(), "testDataSet").withDimension(dims).build();

    std::vector<MoleculeBlock> values;
    values.reserve(numValues);
    for (int i = 0; i < numValues; ++i)
    {
        values.emplace_back("block" + std::to_string(i), 2 * i + 1);
    }

    dataSet.writeData(values);

    std::vector<MoleculeBlock> readBuffer(numValues);
    dataSet.readData(readBuffer);

    for (int i = 0; i < numValues; ++i)
    {
        EXPECT_STREQ(values[i].name.data(), readBuffer[i].name.data());
        EXPECT_EQ(values[i].count, readBuffer[i].count);
    }
}

TEST_F(H5mdCompoundTypeTest, StoresNamesWithGivenSize)
{
    const std::vector<std::string> names = { "Protein", "SOL", "NA" };
    std::vector<MoleculeBlock>     values;
    for (Index i = 0; i < gmx::ssize(names); ++i)
    {
        values.emplace_back(names[i], static_cast<int32_t>(10 * i + 1));
    }

    // The longest name plus its terminating null
    const int nameSize = 8;
    {
        H5mdFixedDataSet<MoleculeBlock> dataSet =
                H5mdDataSetBuilder<MoleculeBlock>(fileid(), "sized")
                        .withDimension({ values.size() })
                        .withMaxStringLength(nameSize)
                        .build();
        dataSet.writeData(values);
    }

    const H5mdDataSetBase<MoleculeBlock> storedDataSet(fileid(), "sized");
    const hid_t                          storedDataType = storedDataSet.storedDataType();
    const auto [nameDataType, nameDataTypeGuard] =
            makeH5mdTypeGuard(H5Tget_member_type(storedDataType, 0));
    EXPECT_EQ(H5Tget_size(nameDataType), static_cast<size_t>(nameSize));
    EXPECT_EQ(H5Tget_member_class(storedDataType, 1), H5T_INTEGER);
    EXPECT_EQ(H5Tget_size(storedDataType), nameSize + sizeof(int32_t)) << "Members must be packed";
    EXPECT_EQ(H5Aexists(storedDataSet.id(), "units"), 0) << "A type without sc_units gets no units";

    H5mdFixedDataSet<MoleculeBlock> dataSet(fileid(), "sized");
    std::vector<MoleculeBlock>      readBuffer(values.size());
    dataSet.readData(readBuffer);
    for (Index i = 0; i < gmx::ssize(names); ++i)
    {
        EXPECT_EQ(readBuffer[i].name.data(), names[i]);
        EXPECT_EQ(readBuffer[i].count, values[i].count);
    }
}

TEST_F(H5mdCompoundTypeTest, ThrowsForNamesThatDoNotFit)
{
    const std::string tooLong(MoleculeBlock::sc_nameCapacity, 'A');
    EXPECT_THROW_GMX(MoleculeBlock(tooLong, 1), InvalidInputError);
    EXPECT_NO_THROW_GMX(MoleculeBlock(tooLong.substr(1), 1))
            << "Must hold the longest name grompp reads";

    EXPECT_THROW_GMX(H5Tclose(MoleculeBlock::constructHdf5DataType(0)), H5mdError);
    EXPECT_THROW_GMX(H5Tclose(MoleculeBlock::constructHdf5DataType(MoleculeBlock::sc_nameCapacity + 1)),
                     H5mdError);
}

TEST_F(H5mdCompoundTypeTest, OpenDataSetThrowsForWrongType)
{
    //! \brief A molecule block record as another writer may store it.
    struct StoredBlock
    {
        std::array<char, 8> name;
        int32_t             count;
    };
    const StoredBlock record = { { 'S', 'O', 'L' }, 3 };
    {
        SCOPED_TRACE("Create a data set of a primitive type and close it at the end of scope");
        H5mdDataSetBuilder<int64_t>(fileid(), "int64").withDimension({ 3 }).build();
    }
    {
        SCOPED_TRACE("Create a compound data set with another name for one member");
        const auto [nameDataType, nameDataTypeGuard] =
                makeH5mdTypeGuard(hdf5DataTypeForFixedLengthString(record.name.size()));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(createCompoundDataType(
                sizeof(StoredBlock),
                { { "molecule_name", HOFFSET(StoredBlock, name), nameDataType },
                  { "count", HOFFSET(StoredBlock, count), hdf5DataTypeFor<int32_t>() } }));
        createRawDataSet(fileid(), "otherMemberName", dataType, record);
    }
    {
        SCOPED_TRACE("Create a compound data set with a variable-length name");
        const auto [nameDataType, nameDataTypeGuard] =
                makeH5mdTypeGuard(hdf5DataTypeForVariableLengthString());
        struct
        {
            const char* name;
            int32_t     count;
        } const variableLengthRecord         = { "SOL", 3 };
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(createCompoundDataType(
                sizeof(variableLengthRecord),
                { { "name", HOFFSET(decltype(variableLengthRecord), name), nameDataType },
                  { "count", HOFFSET(decltype(variableLengthRecord), count), hdf5DataTypeFor<int32_t>() } }));
        createRawDataSet(fileid(), "variableLengthName", dataType, variableLengthRecord);
    }

    EXPECT_THROW_GMX(H5mdFixedDataSet<MoleculeBlock>(fileid(), "int64"), H5mdError)
            << "Must throw when opening a primitive data set";
    EXPECT_THROW_GMX(H5mdFixedDataSet<MoleculeBlock>(fileid(), "otherMemberName"), H5mdError)
            << "Must throw when a member name differs";
    EXPECT_THROW_GMX(H5mdFixedDataSet<MoleculeBlock>(fileid(), "variableLengthName"), H5mdError)
            << "Must throw when the name is not a fixed-length string";
}

TEST_F(H5mdCompoundTypeTest, OpenDataSetThrowsForNamesLongerThanCapacity)
{
    // Reading would cut the names, so such a data set is not opened
    struct LongNameBlock
    {
        std::array<char, MoleculeBlock::sc_nameCapacity + 1> name;
        int32_t                                              count;
    };
    const LongNameBlock record = { {}, 1 };
    const auto [nameDataType, nameDataTypeGuard] =
            makeH5mdTypeGuard(hdf5DataTypeForFixedLengthString(record.name.size()));
    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(createCompoundDataType(
            sizeof(LongNameBlock),
            { { "name", HOFFSET(LongNameBlock, name), nameDataType },
              { "count", HOFFSET(LongNameBlock, count), hdf5DataTypeFor<int32_t>() } }));
    createRawDataSet(fileid(), "longNames", dataType, record);

    EXPECT_THROW_GMX(H5mdFixedDataSet<MoleculeBlock>(fileid(), "longNames"), H5mdError);
}

TEST_F(H5mdCompoundTypeTest, OpensDataSetWithOtherWidthOfMember)
{
    // Another writer may store the count as a 64-bit integer, which the specification allows
    struct StoredBlockInt64
    {
        std::array<char, 8> name;
        int64_t             count;
    };
    const StoredBlockInt64 record = { { 'S', 'O', 'L' }, 1000 };
    {
        const auto [nameDataType, nameDataTypeGuard] =
                makeH5mdTypeGuard(hdf5DataTypeForFixedLengthString(record.name.size()));
        const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(createCompoundDataType(
                sizeof(StoredBlockInt64),
                { { "name", HOFFSET(StoredBlockInt64, name), nameDataType },
                  { "count", HOFFSET(StoredBlockInt64, count), hdf5DataTypeFor<int64_t>() } }));
        createRawDataSet(fileid(), "int64Count", dataType, record);
    }

    H5mdFixedDataSet<MoleculeBlock> dataSet(fileid(), "int64Count");
    std::vector<MoleculeBlock>      readBuffer(1);
    dataSet.readData(readBuffer);
    EXPECT_STREQ(readBuffer[0].name.data(), "SOL");
    EXPECT_EQ(readBuffer[0].count, 1000);
}

TEST_F(H5mdCompoundTypeTest, FrameDataSetRoundTrip)
{
    const std::vector<std::string> names  = { "one", "two", "three" };
    const std::vector<int32_t>     counts = { 12, 5, 30 };

    {
        H5mdFrameDataSet<MoleculeBlock> dataSet =
                H5mdFrameDataSetBuilder<MoleculeBlock>(fileid(), "testDataSet")
                        .withMaxNumFrames(names.size())
                        .withMaxStringLength(6)
                        .build();
        for (Index i = 0; i < gmx::ssize(names); ++i)
        {
            const MoleculeBlock block(names[i], counts[i]);
            dataSet.writeNextFrame(constArrayRefFromArray(&block, 1));
        }
    }

    H5mdFrameDataSet<MoleculeBlock> dataSet(fileid(), "testDataSet");
    ASSERT_EQ(dataSet.numFrames(), names.size());
    // A longer name held before reading is overwritten completely
    MoleculeBlock block("held before reading", 0);
    for (hsize_t i = 0; i < dataSet.numFrames(); ++i)
    {
        dataSet.readFrame(i, arrayRefFromArray(&block, 1));
        EXPECT_STREQ(block.name.data(), names[i].c_str());
        EXPECT_EQ(block.count, counts[i]);
    }
}

} // namespace
} // namespace test
} // namespace gmx
