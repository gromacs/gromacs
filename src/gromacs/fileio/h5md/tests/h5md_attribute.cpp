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
 * \brief
 * Tests for reading and writing of HDF5 attributes in H5MD files.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_attribute.h"

#include <optional>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

using H5mdAttributeTest = H5mdTestBase;

/*! \brief Test the reading and writing of scalar attributes.
 */
TEST_F(H5mdAttributeTest, ScalarAttribute)
{
    // Set variables to test
    std::string       groupName             = "/h5md/testattrs/";
    const int32_t     referenceCreationYear = 2019;
    const int64_t     referenceMembers      = 290;
    const uint32_t    referenceAge          = 21;
    const uint64_t    referenceBirthday     = 11111201;
    const float       referenceHeight       = 162.0;
    const double      referenceScore        = 85.456;
    const std::string referenceName         = "Waai Fu";
    const char*       referenceAddress      = "Lee's Detective Agency, Lungmen, Terra";

    {
        SCOPED_TRACE("Testing H5MD writing of scalar attributes.");
        // Numeric attributes
        setAttribute<int32_t>(fileid(), groupName + "creation_year", referenceCreationYear);
        setAttribute<int64_t>(fileid(), groupName + "members", referenceMembers);
        setAttribute<uint32_t>(fileid(), groupName + "age", referenceAge);
        setAttribute<uint64_t>(fileid(), groupName + "birthday", referenceBirthday);
        setAttribute<float>(fileid(), groupName + "height", referenceHeight);
        setAttribute<double>(fileid(), groupName + "score", referenceScore);

        // String-like attributes
        setAttribute<std::string>(fileid(), groupName + "name", referenceName);
        setAttribute<const char*>(fileid(), groupName + "address", referenceAddress);
    }

    {
        EXPECT_EQ(getAttribute<int32_t>(fileid(), groupName + "creation_year"), referenceCreationYear);
        EXPECT_EQ(getAttribute<int64_t>(fileid(), groupName + "members"), referenceMembers);
        EXPECT_EQ(getAttribute<uint32_t>(fileid(), groupName + "age"), referenceAge);
        EXPECT_EQ(getAttribute<uint64_t>(fileid(), groupName + "birthday"), referenceBirthday);
        EXPECT_EQ(getAttribute<float>(fileid(), groupName + "height"), referenceHeight);
        EXPECT_EQ(getAttribute<double>(fileid(), groupName + "score"), referenceScore);
        EXPECT_EQ(getAttribute<std::string>(fileid(), groupName + "name"), referenceName);
        EXPECT_EQ(getAttribute<std::string>(fileid(), groupName + "address"), referenceAddress);
    }
}

/*! \brief Test the reading and writing of array-like (std::vector) attributes.
 */
TEST_F(H5mdAttributeTest, ArrayAttribute)
{
    std::string                 groupName              = "/h5md/testattrs/";
    const std::vector<int32_t>  referenceResIDs        = { 24, 25, 26, 27, 28, 29 };
    const std::vector<int64_t>  referenceAtomicNumbers = { 7, 1, 1, 1, 6, 1 };
    const std::vector<uint32_t> referenceAtomIDs = { 9999, 10000, 10001, 10002, 10003, 10004 };
    const std::vector<uint64_t> referenceMelt    = { 10000000, 2000000, 3000000,
                                                     4000000,  5000000, 6000000 };
    const std::vector<float>    referenceMasses  = { 14.0027, 1.008, 1.008, 1.008, 12.011, 1.008 };
    const std::vector<double>   referencePositions           = { 0.07334561, 5.038612, 6.275699,
                                                                 0.05687081, 5.024728, 6.373351,
                                                                 0.03213923, 4.963619, 6.223948 };
    const std::vector<std::string> referenceResidueNames     = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                                 "CYS", "GLU", "LEU", "ALA", "ALA" };
    const std::vector<const char*> referenceResidueNamesCStr = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                                 "CYS", "GLU", "LEU", "ALA", "ALA" };

    {
        SCOPED_TRACE("Testing H5MD writing of array-like attributes.");
        // Numeric attributes
        setAttributeVector<int32_t>(fileid(), groupName + "index", referenceResIDs);
        setAttributeVector<int64_t>(fileid(), groupName + "atomic_numbers", referenceAtomicNumbers);
        setAttributeVector<uint32_t>(fileid(), groupName + "atom_ids", referenceAtomIDs);
        setAttributeVector<uint64_t>(fileid(), groupName + "melt", referenceMelt);
        setAttributeVector<float>(fileid(), groupName + "masses", referenceMasses);
        setAttributeVector<double>(fileid(), groupName + "positions", referencePositions);

        // String-like attributes
        setAttributeVector<const char*>(fileid(), groupName + "residue_names_c", referenceResidueNamesCStr);
        setAttributeVector<std::string>(fileid(), groupName + "residue_names", referenceResidueNames);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of array-like attributes.");

        EXPECT_EQ(getAttributeVector<int32_t>(fileid(), groupName + "index"), referenceResIDs);
        EXPECT_EQ(getAttributeVector<int64_t>(fileid(), groupName + "atomic_numbers"),
                  referenceAtomicNumbers);
        EXPECT_EQ(getAttributeVector<uint32_t>(fileid(), groupName + "atom_ids"), referenceAtomIDs);
        EXPECT_EQ(getAttributeVector<uint64_t>(fileid(), groupName + "melt"), referenceMelt);
        EXPECT_EQ(getAttributeVector<float>(fileid(), groupName + "masses"), referenceMasses);
        EXPECT_EQ(getAttributeVector<double>(fileid(), groupName + "positions"), referencePositions);
        EXPECT_EQ(getAttributeVector<std::string>(fileid(), groupName + "residue_names"),
                  referenceResidueNames);
        EXPECT_EQ(getAttributeVector<std::string>(fileid(), groupName + "residue_names_c"),
                  referenceResidueNames);
    }
}

/*! \brief Test edge cases for reading and writing of scalar attributes.
 *
 * The int32_t and std::string are tested.
 * This test checks the behavior of the H5MD attribute functions when dealing with
 * getting with a mismatched data types,
 * getting from an empty path,
 * getting from an invalid container,
 * and writing to an existing attribute.
 */
TEST_F(H5mdAttributeTest, ScalarAttributeEdgeCases)
{
    std::string groupName = "/h5md/testattrs/";

    const int32_t     refNumber = 171;
    const std::string refName   = "Margaret Nearl";
    setAttribute<int32_t>(fileid(), groupName + "height", refNumber);
    setAttribute<std::string>(fileid(), groupName + "name", refName);

    {
        SCOPED_TRACE("Testing H5MD reading from a mismatched data type.");
        // Supposed to fail to match the required type and the attribute type
        EXPECT_THROW(getAttribute<int64_t>(fileid(), groupName + "height"), gmx::FileIOError);
        EXPECT_THROW(getAttribute<int64_t>(fileid(), groupName + "name"), gmx::FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD write to an existing attribute.");
        // Supposed to fail to rewrite to an existing attribute
        EXPECT_THROW(setAttribute<std::string>(fileid(), groupName + "name", "Maria Nearl"),
                     gmx::FileIOError);
        EXPECT_THROW(setAttribute<int32_t>(fileid(), groupName + "height", 165), gmx::FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of invalid container or path.");
        // Return nullopt upon reading a non-existing attribute
        EXPECT_EQ(getAttribute<int32_t>(fileid(), "NotAnAttribute"), std::nullopt);
        // Return nullopt upon empty attribute name
        EXPECT_EQ(getAttribute<std::string>(fileid(), ""), std::nullopt);
        // Return nullopt upon invalid container
        EXPECT_EQ(getAttribute<int32_t>(H5I_INVALID_HID, groupName + "height"), std::nullopt);
    }
}

/*! \brief Test edge cases for reading and writing of vector-like attributes.
 *
 * The int32_t and std::string are tested.
 * This test checks the behavior of the H5MD attribute functions when dealing with
 * getting with a mismatched data types,
 * getting from an empty path,
 * getting from an invalid container,
 * getting with an empty vector,
 * writing of non-canonical strings,
 * and writing to an existing attribute.
 */
TEST_F(H5mdAttributeTest, VectorAttributeEdgeCases)
{
    std::string                    groupName           = "/h5md/testattrs/";
    const std::vector<int32_t>     referenceResIDs     = { 24, 25, 26, 27, 28, 29 };
    const std::vector<std::string> refNames            = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                           "CYS", "GLU", "LEU", "ALA", "ALA" };
    const std::vector<std::string> nonCanonicalStrings = {
        "\ttabs\nand\bback\band\tnew\nlines", // control characters
        "100%&#@//\\and.,counting!~<>",       // mildly unusual symbols
        "Hallå!",                             // UTF-8: Swedish
        "你好",                               // UTF-8: Chinese
        "مرحبا",                              // UTF-8: Arabic
        "Привет",                             // UTF-8: Russian
        "こんにちは",                         // UTF-8: Japanese
        "안녕하세요",                         // UTF-8: Korean
        "😀😃😄😁🔥",                         // UTF-8: Emojis
        "This is a very long string that exceeds the typical length of a residue name, "
        "but it is still a valid string for testing purposes." // long string
    };

    setAttributeVector<int32_t>(fileid(), groupName + "index", referenceResIDs);
    setAttributeVector<std::string>(fileid(), groupName + "residue_names", refNames);
    setAttributeVector<std::string>(fileid(), groupName + "non_canonical_strings", nonCanonicalStrings);

    {
        SCOPED_TRACE("Testing H5MD reading from a mismatched data type.");
        EXPECT_THROW(getAttributeVector<int64_t>(fileid(), groupName + "index"), gmx::FileIOError);
        EXPECT_THROW(getAttributeVector<int64_t>(fileid(), groupName + "residue_names"), gmx::FileIOError);
    }

    {
        // NOTE: Throw if the attribute already exists (!5205)
        SCOPED_TRACE("Testing H5MD write to an existing attribute.");
        const std::vector<int32_t>     resID2   = { 1, 2, 3, 4 };
        const std::vector<std::string> resName2 = { "ALA", "GLY", "PHE", "VAL" };
        EXPECT_THROW(setAttributeVector<int32_t>(fileid(), groupName + "index", resID2), gmx::FileIOError);
        EXPECT_THROW(setAttributeVector<std::string>(fileid(), groupName + "residue_names", resName2),
                     gmx::FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD writing an empty vector.");
        std::vector<int32_t> emptyIntegerVector;
        EXPECT_THROW(setAttributeVector<int32_t>(fileid(), groupName + "empty_int", emptyIntegerVector),
                     gmx::FileIOError);
        std::vector<std::string> emptyStringVector;
        EXPECT_THROW(setAttributeVector<std::string>(fileid(), groupName + "empty_str", emptyStringVector),
                     gmx::FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD reading from an invalid container or path.");
        // Return nullopt upon reading a non-existing attribute
        EXPECT_EQ(getAttributeVector<int32_t>(fileid(), "NotAnAttribute"), std::nullopt);
        // Return nullopt upon empty attribute name
        EXPECT_EQ(getAttributeVector<std::string>(fileid(), ""), std::nullopt);
        // Return nullopt upon invalid container
        EXPECT_EQ(getAttributeVector<int32_t>(H5I_INVALID_HID, groupName + "index"), std::nullopt);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of previously written attributes.");
        auto readNonCanonicalStrings =
                getAttributeVector<std::string>(fileid(), groupName + "non_canonical_strings");
        ASSERT_TRUE(readNonCanonicalStrings.has_value());
        EXPECT_EQ(*readNonCanonicalStrings, nonCanonicalStrings);
        auto readResIDs = getAttributeVector<int32_t>(fileid(), groupName + "index");
        ASSERT_TRUE(readResIDs.has_value());
        EXPECT_EQ(*readResIDs, referenceResIDs);
        auto readResNames = getAttributeVector<std::string>(fileid(), groupName + "residue_names");
        ASSERT_TRUE(readResNames.has_value());
        EXPECT_EQ(*readResNames, refNames);
    }
}


} // namespace
} // namespace test
} // namespace gmx
