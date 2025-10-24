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
 * \brief Tests for H5MD utility functions.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_util.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_framedatasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"

namespace gmx
{
namespace test
{
namespace
{

//! \brief Test fixture alias for functions checking whether objects exist.
using H5mdObjectExistsTest = H5mdTestBase;

//! \brief Test fixture alias for functions asserting HDF5 handles.
using H5mdHandleIsValidTest = H5mdTestBase;

//! \brief Test fixture alias for handle name utilities.
using H5mdPathsOfHandleTest = H5mdTestBase;

TEST_F(H5mdObjectExistsTest, FindsObjectsInFileRoot)
{
    EXPECT_FALSE(objectExists(fileid(), "testGroup"));
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    EXPECT_TRUE(objectExists(fileid(), "testGroup"));

    EXPECT_FALSE(objectExists(fileid(), "testDataSet"));
    H5mdFrameDataSetBuilder<float>(fileid(), "testDataSet").build();
    EXPECT_TRUE(objectExists(fileid(), "testDataSet"));
}

TEST_F(H5mdObjectExistsTest, FindsObjectsInGroups)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));

    EXPECT_FALSE(objectExists(group, "testSubGroup"));
    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    EXPECT_TRUE(objectExists(group, "testSubGroup"));

    EXPECT_FALSE(objectExists(group, "testDataSet"));
    H5mdFrameDataSetBuilder<float>(group, "testDataSet").build();
    EXPECT_TRUE(objectExists(group, "testDataSet"));
}

TEST_F(H5mdObjectExistsTest, DoesNotSearchForObjectInsideSubGroups)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));

    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    EXPECT_FALSE(objectExists(fileid(), "testSubGroup")) << "Must not find sub group in root";
    EXPECT_TRUE(objectExists(group, "testSubGroup")) << "Must find sub group inside parent group";

    H5mdFrameDataSetBuilder<float>(group, "testDataSet").build();
    EXPECT_FALSE(objectExists(fileid(), "testDataSet")) << "Must not find data set in root";
    EXPECT_TRUE(objectExists(group, "testDataSet")) << "Must find data set inside parent group";
}

TEST_F(H5mdObjectExistsTest, SearchesExplicitMultiLevelPath)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "/path/to/testGroup"));
    const auto [subGroup, subGroupGuard] = makeH5mdGroupGuard(createGroup(group, "testSubGroup"));
    H5mdFrameDataSetBuilder<float>(subGroup, "testDataSet").build();

    EXPECT_TRUE(objectExists(subGroup, "testDataSet"));
    EXPECT_TRUE(objectExists(group, "testSubGroup/testDataSet"));
    EXPECT_TRUE(objectExists(fileid(), "/path/to/testGroup/testSubGroup/testDataSet"));
}

TEST_F(H5mdHandleIsValidTest, ReturnTrueForValidHandlesToObjects)
{
    EXPECT_TRUE(handleIsValid(fileid()));

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    EXPECT_TRUE(handleIsValid(group));

    const H5mdDataSetBase<float> dataSet =
            H5mdFrameDataSetBuilder<float>(fileid(), "testDataSet").build();
    EXPECT_TRUE(handleIsValid(dataSet.id()));

    const auto [dataType, dataTypeGuard] = makeH5mdTypeGuard(H5Dget_type(dataSet.id()));
    EXPECT_TRUE(handleIsValid(dataType));
}

TEST_F(H5mdHandleIsValidTest, ReturnFalseForInvalidHandles)
{
    EXPECT_FALSE(handleIsValid(H5I_INVALID_HID));
}

TEST_F(H5mdHandleIsValidTest, ReturnFalseForHandlesToClosedObjects)
{
    hid_t handleToClose = H5I_INVALID_HID;
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
        handleToClose                  = group;
        ASSERT_TRUE(handleIsValid(handleToClose))
                << "Sanity check: must be valid before scope guard exits scope";
    }
    EXPECT_FALSE(handleIsValid(handleToClose));
}

TEST_F(H5mdPathsOfHandleTest, HandlePathOfGroups)
{
    {
        const auto [handleTest, guard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
        const auto path                = getHandlePath(handleTest);
        ASSERT_TRUE(path.has_value());
        EXPECT_EQ(path.value(), "/testGroup");
    }

    {
        SCOPED_TRACE("Test long path");
        const auto [handleTest, guard] = makeH5mdGroupGuard(
                createGroup(fileid(), "/This/is/a/extremely/long/path/to/a/testgp"));
        const auto path = getHandlePath(handleTest);
        ASSERT_TRUE(path.has_value());
        EXPECT_EQ(path.value(), "/This/is/a/extremely/long/path/to/a/testgp");
    }

    {
        SCOPED_TRACE("Test path with trailing slash");
        const auto [handleTest, guard] =
                makeH5mdGroupGuard(createGroup(fileid(), "slash/in/the/end/"));
        const auto path = getHandlePath(handleTest);
        ASSERT_TRUE(path.has_value());
        EXPECT_EQ(path.value(), "/slash/in/the/end");
    }
}

TEST_F(H5mdPathsOfHandleTest, BaseNameOfGroups)
{
    {
        const auto [handleTest, guard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
        EXPECT_EQ(getHandleBaseName(handleTest).value_or("NoAName"), "testGroup");
    }

    {
        SCOPED_TRACE("Test long path");
        const auto [handleTest, guard] = makeH5mdGroupGuard(
                createGroup(fileid(), "/This/is/a/extremely/long/path/to/a/testgp"));
        EXPECT_EQ(getHandleBaseName(handleTest).value_or("NoAName"), "testgp");
    }

    {
        SCOPED_TRACE("Test path with trailing slash");
        const auto [handleTest, guard] =
                makeH5mdGroupGuard(createGroup(fileid(), "slash/in/the/end/"));
        EXPECT_EQ(getHandleBaseName(handleTest).value_or("NoAName"), "end");
    }
}

TEST_F(H5mdPathsOfHandleTest, HandleNameOfDatasets)
{
    const auto [grpHandle, guard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    const H5mdDataSetBase<float> dataSet =
            H5mdFrameDataSetBuilder<float>(grpHandle, "testDataSet").build();
    const auto datasetPath = getHandlePath(dataSet.id());
    ASSERT_TRUE(datasetPath.has_value());
    EXPECT_EQ(datasetPath.value(), "/testGroup/testDataSet");
}

TEST_F(H5mdPathsOfHandleTest, BaseNameOfDatasets)
{
    const auto [grpHandle, guard] = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
    const H5mdDataSetBase<float> dataSet =
            H5mdFrameDataSetBuilder<float>(grpHandle, "testDataSet").build();
    EXPECT_EQ(getHandleBaseName(dataSet.id()).value_or("NoAName"), "testDataSet");
}

TEST_F(H5mdPathsOfHandleTest, ReturnEmptyForInvalidHandle)
{
    {
        SCOPED_TRACE("Test invalid HDF5 handle");
        EXPECT_EQ(getHandlePath(H5I_INVALID_HID), std::nullopt);
        EXPECT_EQ(getHandleBaseName(H5I_INVALID_HID), std::nullopt);
    }

    {
        SCOPED_TRACE("Test uninitialized HDF5 handle");
        hid_t handleTest = 1;
        EXPECT_EQ(getHandlePath(handleTest), std::nullopt);
        EXPECT_EQ(getHandleBaseName(handleTest), std::nullopt);
    }

    {
        SCOPED_TRACE("Test closed HDF5 handle");
        hid_t handleTest;
        {
            // Within the Scope guard, the handle is valid
            auto resource = makeH5mdGroupGuard(createGroup(fileid(), "testGroup"));
            handleTest    = resource.first;
            ASSERT_EQ(getHandlePath(handleTest).value_or("NotAPath"), "/testGroup")
                    << "Sanity check: must be valid before scope guard exits scope";
        }
        EXPECT_EQ(getHandlePath(handleTest), std::nullopt);
        EXPECT_EQ(getHandleBaseName(handleTest), std::nullopt);
    }
}

TEST_F(H5mdPathsOfHandleTest, NameOfRoot)
{
    {
        SCOPED_TRACE("The full name of the root group should be /");
        EXPECT_EQ(getHandlePath(fileid()).value_or("NotAPath"), "/");
    }

    {
        SCOPED_TRACE("The base name of the root group should be an empty string");
        EXPECT_EQ(getHandleBaseName(fileid()).value_or("NoAName"), "");
    }
}


TEST(H5mdStringBuffer, EstimateSizeOfStringBuffer)
{
    {
        SCOPED_TRACE("Empty string vector");

        std::vector<std::string> empty_strings = { "", "", "", "" };
        auto [maxLength, strCount]             = estimateBufferSize(empty_strings.begin(),
                                                        empty_strings.end(),
                                                        [](const auto& it) { return it->c_str(); });
        EXPECT_EQ(maxLength, 0) << "Max length should be 0 for empty strings";
        EXPECT_EQ(strCount, empty_strings.size()) << "String count should match sequence size";
    }

    {
        SCOPED_TRACE("String vector with same length");

        std::vector<std::string> sequence = { "ALA", "SER", "ARG", "MET", "ACE", "NH2" };
        auto [maxLength, strCount]        = estimateBufferSize(
                sequence.begin(), sequence.end(), [](const auto& it) { return it->c_str(); });
        EXPECT_EQ(maxLength, 3) << "Max length should be 3 for uniform strings";
        EXPECT_EQ(strCount, sequence.size()) << "String count should match sequence size";
    }

    {
        SCOPED_TRACE("String vector with varied length");

        std::vector<std::string> strings = { "one",      "two",      "three",   "four",
                                             "five",     "six",      "seven",   "eight",
                                             "nine",     "ten",      "eleven",  "twelve",
                                             "thirteen", "fourteen", "fifteen", "sixteen" };
        auto [maxLength, strCount]       = estimateBufferSize(
                strings.begin(), strings.end(), [](const auto& it) { return it->c_str(); });
        // Maximum length 8 ("thirteen" and "fourteen")
        EXPECT_EQ(maxLength, 8) << "Max length should be 8 for varied strings";
        EXPECT_EQ(strCount, strings.size()) << "String count should match sequence size";
    }
}

TEST(H5mdStringBuffer, PackStringsInSameLength)
{
    std::vector<std::string> sequence = { "ALA", "SER", "ARG", "MET", "ACE", "NH2" };

    auto [maxLength, strCount] = estimateBufferSize(
            sequence.begin(), sequence.end(), [](const auto& it) { return it->c_str(); });
    std::vector<char> buffer(sequence.size() * (maxLength + 1), '\0');
    buffer = packBufferViaIterator(sequence.begin(),
                                   sequence.end(),
                                   std::move(buffer),
                                   maxLength,
                                   [](const auto& it) { return it->c_str(); });

    EXPECT_EQ(buffer.size(), strCount * (maxLength + 1))
            << "Buffer size should match the strCount x (maxLength + 1)";
    for (size_t i = 0; i < sequence.size(); ++i)
    {
        EXPECT_EQ(std::string(buffer.data() + i * (maxLength + 1)), sequence[i])
                << "Buffer content should match the original string";
    }
}

TEST(H5mdStringBuffer, PackStringsInVariedLength)
{
    std::vector<std::string> strings = { "one",      "two",      "three",   "four",
                                         "five",     "six",      "seven",   "eight",
                                         "nine",     "ten",      "eleven",  "twelve",
                                         "thirteen", "fourteen", "fifteen", "sixteen" };

    auto [maxLength, strCount] = estimateBufferSize(
            strings.begin(), strings.end(), [](const auto& it) { return it->c_str(); });
    std::vector<char> buffer(strings.size() * (maxLength + 1), '\0');
    buffer = packBufferViaIterator(strings.begin(),
                                   strings.end(),
                                   std::move(buffer),
                                   maxLength,
                                   [](const auto& it) { return it->c_str(); });

    EXPECT_EQ(buffer.size(), strCount * (maxLength + 1))
            << "Buffer size should match the strCount x (maxLength + 1)";
    for (size_t i = 0; i < strings.size(); ++i)
    {
        EXPECT_EQ(std::string(buffer.data() + i * (maxLength + 1)), strings[i])
                << "Buffer content should match the original string";
    }
}

} // namespace
} // namespace test
} // namespace gmx
