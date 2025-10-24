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
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"

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
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
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
        setAttribute(group, "creation_year", referenceCreationYear);
        setAttribute(group, "members", referenceMembers);
        setAttribute(group, "age", referenceAge);
        setAttribute(group, "birthday", referenceBirthday);
        setAttribute(group, "height", referenceHeight);
        setAttribute(group, "score", referenceScore);

        // String-like attributes
        setAttribute(group, "name", referenceName);
        setAttribute(group, "address", referenceAddress);
        setAttribute(group, "hobby", "Kungfu");
    }

    {
        EXPECT_EQ(getAttribute<int32_t>(group, "creation_year"), referenceCreationYear);
        EXPECT_EQ(getAttribute<int64_t>(group, "members"), referenceMembers);
        EXPECT_EQ(getAttribute<uint32_t>(group, "age"), referenceAge);
        EXPECT_EQ(getAttribute<uint64_t>(group, "birthday"), referenceBirthday);
        EXPECT_EQ(getAttribute<float>(group, "height"), referenceHeight);
        EXPECT_EQ(getAttribute<double>(group, "score"), referenceScore);
        EXPECT_EQ(getAttribute<std::string>(group, "name"), referenceName);
        EXPECT_EQ(getAttribute<std::string>(group, "address"), referenceAddress);
    }
}

/*! \brief Test the reading and writing of array-like (std::vector) attributes.
 *
 * Implicitly convert std::vector to ArrayRef when writing attributes.
 */
TEST_F(H5mdAttributeTest, NumericAttributeViaVector)
{
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const std::vector<int32_t>  referenceResIDs        = { 24, 25, 26, 27, 28, 29 };
    const std::vector<int64_t>  referenceAtomicNumbers = { 7, 1, 1, 1, 6, 1 };
    const std::vector<uint32_t> referenceAtomIDs = { 9999, 10000, 10001, 10002, 10003, 10004 };
    const std::vector<uint64_t> referenceMelt    = { 10000000, 2000000, 3000000,
                                                     4000000,  5000000, 6000000 };
    const std::vector<float>    referenceMasses  = { 14.0027, 1.008, 1.008, 1.008, 12.011, 1.008 };
    const std::vector<double>   referencePositions = { 0.07334561, 5.038612, 6.275699,
                                                       0.05687081, 5.024728, 6.373351,
                                                       0.03213923, 4.963619, 6.223948 };

    {
        SCOPED_TRACE("Testing H5MD writing of array-like attributes.");
        // Numeric attributes
        setAttributeVector<int32_t>(group, "index", referenceResIDs);
        setAttributeVector<int64_t>(group, "atomic_numbers", referenceAtomicNumbers);
        setAttributeVector<uint32_t>(group, "atom_ids", referenceAtomIDs);
        setAttributeVector<uint64_t>(group, "melt", referenceMelt);
        setAttributeVector<float>(group, "masses", referenceMasses);
        setAttributeVector<double>(group, "positions", referencePositions);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of array-like attributes.");

        EXPECT_EQ(getAttributeVector<int32_t>(group, "index"), referenceResIDs);
        EXPECT_EQ(getAttributeVector<int64_t>(group, "atomic_numbers"), referenceAtomicNumbers);
        EXPECT_EQ(getAttributeVector<uint32_t>(group, "atom_ids"), referenceAtomIDs);
        EXPECT_EQ(getAttributeVector<uint64_t>(group, "melt"), referenceMelt);
        EXPECT_EQ(getAttributeVector<float>(group, "masses"), referenceMasses);
        EXPECT_EQ(getAttributeVector<double>(group, "positions"), referencePositions);
    }
}

/*! \brief Test the reading and writing of array-like (gmx::ArrayRef) attributes.
 */
TEST_F(H5mdAttributeTest, NumericAttributeViaArrayRef)
{
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const std::vector<int32_t>  referenceResIDs        = { 24, 25, 26, 27, 28, 29 };
    const std::vector<int64_t>  referenceAtomicNumbers = { 7, 1, 1, 1, 6, 1 };
    const std::vector<uint32_t> referenceAtomIDs = { 9999, 10000, 10001, 10002, 10003, 10004 };
    const std::vector<uint64_t> referenceMelt    = { 10000000, 2000000, 3000000,
                                                     4000000,  5000000, 6000000 };
    const std::vector<float>    referenceMasses  = { 14.0027, 1.008, 1.008, 1.008, 12.011, 1.008 };
    const std::vector<double>   referencePositions = { 0.07334561, 5.038612, 6.275699,
                                                       0.05687081, 5.024728, 6.373351,
                                                       0.03213923, 4.963619, 6.223948 };

    {
        SCOPED_TRACE("Testing H5MD writing of array-like attributes.");

        setAttributeVector<int32_t>(group, "index", referenceResIDs);
        setAttributeVector<int64_t>(group, "atomic_numbers", referenceAtomicNumbers);
        setAttributeVector<uint32_t>(group, "atom_ids", referenceAtomIDs);
        setAttributeVector<uint64_t>(group, "melt", referenceMelt);
        setAttributeVector<float>(group, "masses", referenceMasses);
        setAttributeVector<double>(group, "positions", referencePositions);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of array-like attributes.");

        EXPECT_EQ(getAttributeVector<int32_t>(group, "index"), referenceResIDs);
        EXPECT_EQ(getAttributeVector<int64_t>(group, "atomic_numbers"), referenceAtomicNumbers);
        EXPECT_EQ(getAttributeVector<uint32_t>(group, "atom_ids"), referenceAtomIDs);
        EXPECT_EQ(getAttributeVector<uint64_t>(group, "melt"), referenceMelt);
        EXPECT_EQ(getAttributeVector<float>(group, "masses"), referenceMasses);
        EXPECT_EQ(getAttributeVector<double>(group, "positions"), referencePositions);
    }
}

TEST_F(H5mdAttributeTest, StringAttributeViaArrayRef)
{
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const std::vector<std::string> refSeq     = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                  "CYS", "GLU", "LEU", "ALA", "ALA" };
    const std::vector<const char*> refSeqCStr = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                  "CYS", "GLU", "LEU", "ALA", "ALA" };

    auto refSeq_ref     = makeArrayRef(refSeq);
    auto refSeqCStr_ref = makeArrayRef(refSeqCStr);
    {
        SCOPED_TRACE("Writing string attributes with an empty buffer.");

        auto buffer1 =
                setAttributeStringVector(group, "sequence", {}, refSeq_ref.begin(), refSeq_ref.end());

        auto buffer2 = setAttributeStringVector(
                group, "sequence_c", {}, refSeqCStr_ref.begin(), refSeqCStr_ref.end());

        EXPECT_EQ(buffer1.size(), 4 * refSeq.size()) << "Buffer size is incorrect.";
        EXPECT_EQ(buffer2.size(), 4 * refSeqCStr.size()) << "Buffer size is incorrect.";
    }

    {
        const auto ret1 = getAttributeVector<std::string>(group, "sequence");
        const auto ret2 = getAttributeVector<std::string>(group, "sequence_c");
        EXPECT_TRUE(ret1.has_value());
        EXPECT_TRUE(ret2.has_value());
        EXPECT_EQ(ret1.value(), refSeq);
        EXPECT_EQ(ret2.value(), refSeq);
    }
}

TEST_F(H5mdAttributeTest, StringBufferReuse)
{
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const std::vector<std::string> refSeq     = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                  "CYS", "GLU", "LEU", "ALA", "ALA" };
    const std::vector<const char*> refSeqCStr = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                  "CYS", "GLU", "LEU", "ALA", "ALA" };

    auto refSeq_ref     = makeArrayRef(refSeq);
    auto refSeqCStr_ref = makeArrayRef(refSeqCStr);
    {
        SCOPED_TRACE("Testing r-value to initialize the buffer.");

        auto buffer =
                setAttributeStringVector(group, "sequence0", {}, refSeq_ref.begin(), refSeq_ref.end());
        auto bufferCStr = setAttributeStringVector(
                group, "sequence_c0", std::move(buffer), refSeqCStr_ref.begin(), refSeqCStr_ref.end());
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence0"), refSeq);
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence_c0"), refSeq);
    }

    {
        SCOPED_TRACE("Testing moving an empty buffer to set the attribute.");

        std::vector<char> buffer;
        buffer = setAttributeStringVector(
                group, "sequence1", std::move(buffer), refSeq_ref.begin(), refSeq_ref.end());
        buffer = setAttributeStringVector(
                group, "sequence_c1", std::move(buffer), refSeqCStr_ref.begin(), refSeqCStr_ref.end());

        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence1"), refSeq);
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence_c1"), refSeq);
    }

    {
        SCOPED_TRACE("Testing moving an known-size buffer to set the attribute.");
        std::vector<char> buffer = std::vector<char>(4 * refSeq_ref.size(), '\0');
        buffer                   = setAttributeStringVector(
                group, "sequence2", std::move(buffer), refSeq_ref.begin(), refSeq_ref.end());
        buffer = setAttributeStringVector(
                group, "sequence_c2", std::move(buffer), refSeqCStr_ref.begin(), refSeqCStr_ref.end());

        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence2"), refSeq);
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence_c2"), refSeq);
    }

    {
        SCOPED_TRACE("Testing moving an preallocated large-enough buffer to set the attribute.");

        std::vector<char> buffer = std::vector<char>(256 * refSeq_ref.size(), '\0');
        buffer                   = setAttributeStringVector(
                group, "sequence3", std::move(buffer), refSeq_ref.begin(), refSeq_ref.end());
        buffer = setAttributeStringVector(
                group, "sequence_c3", std::move(buffer), refSeqCStr_ref.begin(), refSeqCStr_ref.end());

        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence3"), refSeq);
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence_c3"), refSeq);
    }

    {
        SCOPED_TRACE("Testing allocating a small buffer and resize it on-the-fly.");
        std::vector<char> buffer = std::vector<char>(2 * refSeq_ref.size(), '\0');
        EXPECT_EQ(buffer.size(), 2 * refSeq_ref.size());
        buffer = setAttributeStringVector(
                group, "sequence4", std::move(buffer), refSeq_ref.begin(), refSeq_ref.end());
        EXPECT_EQ(buffer.size(), 4 * refSeq_ref.size());
        buffer = setAttributeStringVector(
                group, "sequence_c4", std::move(buffer), refSeqCStr_ref.begin(), refSeqCStr_ref.end());
        EXPECT_EQ(buffer.size(), 4 * refSeq_ref.size());

        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence4"), refSeq);
        EXPECT_EQ(getAttributeVector<std::string>(group, "sequence_c4"), refSeq);
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
    std::string groupName       = "/h5md/testattrs/";
    auto [group, groupGuard]    = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const int32_t     refNumber = 171;
    const std::string refName   = "Margaret Nearl";
    setAttribute(group, "height", refNumber);
    setAttribute(group, "name", refName);

    {
        SCOPED_TRACE("Testing H5MD reading from a mismatched data type.");
        // Supposed to fail to match the required type and the attribute type
        EXPECT_THROW(getAttribute<int64_t>(group, "height"), FileIOError);
        EXPECT_THROW(getAttribute<int64_t>(group, "name"), FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD write to an existing attribute.");
        // Supposed to fail to rewrite to an existing attribute
        EXPECT_THROW(setAttribute(group, "name", "Maria Nearl"), FileIOError);
        EXPECT_THROW(setAttribute(group, "height", 165), FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of invalid container or path.");
        // Return nullopt upon reading a non-existing attribute
        EXPECT_EQ(getAttribute<int32_t>(group, "NotAnAttribute"), std::nullopt);
        // Return nullopt upon empty attribute name
        EXPECT_EQ(getAttribute<std::string>(group, ""), std::nullopt);
        // Return nullopt upon invalid container
        EXPECT_EQ(getAttribute<int32_t>(H5I_INVALID_HID, "height"), std::nullopt);
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
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    const std::vector<int32_t>     referenceResIDs     = { 24, 25, 26, 27, 28, 29 };
    const std::vector<std::string> refNames            = { "LYS", "VAL", "PHE", "GLY", "ARG",
                                                           "CYS", "GLU", "LEU", "ALA", "ALA" };
    const std::vector<std::string> nonCanonicalStrings = {
        // control characters
        "\ttabs\nand\bback\band\tnew\nlines",
        // mildly unusual symbols
        "100%&#@//\\and.,counting!~<>",
        // UTF-8: Swedish
        "Hallå!",
        // UTF-8: Chinese
        "你好",
        // UTF-8: Arabic
        "مرحبا",
        // UTF-8: Russian
        "Привет",
        // UTF-8: Japanese
        "こんにちは",
        // UTF-8: Korean
        "안녕하세요",
        // UTF-8: Emojis
        "😀😃😄😁🔥",
        // long string
        "This is a very long string that exceeds the typical length of a residue name, "
        "but it is still a valid string for testing purposes."
    };

    setAttributeVector<int32_t>(group, "index", referenceResIDs);

    // Use an empty buffer to initialize the reusable buffer
    auto refNames_ref = makeArrayRef(refNames);
    auto NCString_ref = makeArrayRef(nonCanonicalStrings);
    auto buffer       = setAttributeStringVector(
            group, "residue_names", {}, refNames_ref.begin(), refNames_ref.end());
    buffer = setAttributeStringVector(
            group, "non_canonical_strings", {}, NCString_ref.begin(), NCString_ref.end());

    {
        SCOPED_TRACE("Testing H5MD reading from a mismatched data type.");
        EXPECT_THROW(getAttributeVector<int64_t>(group, "index"), FileIOError);
        EXPECT_THROW(getAttributeVector<int64_t>(group, "residue_names"), FileIOError);
    }

    {
        // NOTE: Throw if the attribute already exists (!5205)
        SCOPED_TRACE("Testing H5MD write to an existing attribute.");
        const std::vector<int32_t>     resID2       = { 1, 2, 3, 4 };
        const std::vector<std::string> resName2     = { "ALA", "GLY", "PHE", "VAL" };
        auto                           resName2_ref = makeArrayRef(resName2);
        EXPECT_THROW(setAttributeVector<int32_t>(group, "index", resID2), FileIOError);
        EXPECT_THROW(buffer = setAttributeStringVector(group,
                                                       "residue_names",
                                                       std::move(buffer),
                                                       resName2_ref.begin(),
                                                       resName2_ref.end()),
                     FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD writing an empty vector.");
        // NOTE: Empty vector throws in H5Awrite
        std::vector<int32_t> emptyIntegerVector{};
        EXPECT_THROW(setAttributeVector<int32_t>(group, "empty_int", emptyIntegerVector), FileIOError);

        // NOTE: Empty buffer throws as expected similar to numerical attributes
        std::vector<std::string>    emptyStringVector{};
        ArrayRef<const std::string> emptyStringRef = makeArrayRef(emptyStringVector);
        EXPECT_THROW(buffer = setAttributeStringVector(
                             group, "empty_string", {}, emptyStringRef.begin(), emptyStringRef.end()),
                     FileIOError);
    }

    {
        SCOPED_TRACE("Testing H5MD reading from an invalid container or path.");
        // Return nullopt upon reading a non-existing attribute
        EXPECT_EQ(getAttributeVector<int32_t>(group, "NotAnAttribute"), std::nullopt);
        // Return nullopt upon empty attribute name
        EXPECT_EQ(getAttributeVector<std::string>(group, ""), std::nullopt);
        // Return nullopt upon invalid container
        EXPECT_EQ(getAttributeVector<int32_t>(H5I_INVALID_HID, "index"), std::nullopt);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of previously written attributes.");
        auto readNonCanonicalStrings =
                getAttributeVector<std::string>(group, "non_canonical_strings");
        ASSERT_TRUE(readNonCanonicalStrings.has_value());
        EXPECT_EQ(*readNonCanonicalStrings, nonCanonicalStrings);
        auto readResIDs = getAttributeVector<int32_t>(group, "index");
        ASSERT_TRUE(readResIDs.has_value());
        EXPECT_EQ(*readResIDs, referenceResIDs);
        auto readResNames = getAttributeVector<std::string>(group, "residue_names");
        ASSERT_TRUE(readResNames.has_value());
        EXPECT_EQ(*readResNames, refNames);
    }
}


/*! \brief
 * Creates dummy topology with two differently sized residues.
 *
 * Residue begin and end are set to allow checking routines
 * that make use of the boundaries.
 *
 * \returns The residue ranges.
 */
void createTwoResidueTopology(gmx_mtop_t* mtop, const ArrayRef<char*> atomNames)
{
    auto& moltype        = mtop->moltype.emplace_back();
    int   residueOneSize = 5;
    int   residueTwoSize = 4;
    moltype.atoms.nr     = residueOneSize + residueTwoSize;
    snew(moltype.atoms.atom, residueOneSize + residueTwoSize);
    snew(moltype.atoms.atomname, residueOneSize + residueTwoSize);
    auto atomNameIterator = atomNames.begin();
    for (int i = 0; i < residueOneSize; i++, atomNameIterator++)
    {
        moltype.atoms.atom[i].resind = 0;
        moltype.atoms.atomname[i]    = atomNameIterator.data();
    }
    for (int i = residueOneSize; i < residueOneSize + residueTwoSize; i++, atomNameIterator++)
    {
        moltype.atoms.atom[i].resind = 1;
        moltype.atoms.atomname[i]    = atomNameIterator.data();
    }

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;
    mtop->natoms           = moltype.atoms.nr * mtop->molblock[0].nmol;
}

/*! \brief Test the reading and writing of topology information via iterating over atoms.
 */
TEST_F(H5mdAttributeTest, AtomNameAttributes)
{
    // Mark is assuming that we actually do want to write atom names
    // to the H5md file - but we should really have a working
    // prototype in mdrun to which we are adding the things that we
    // have in the plan, whereever that is.
    std::vector<std::string> atomNames = { "AA", "BBB", "C", "DD", "EEE", "F", "GGG", "HH", "III" };
    std::vector<char*>       handlesToAtomNames;
    for (const auto& atomName : atomNames)
    {
        handlesToAtomNames.push_back(const_cast<char*>(atomName.data()));
    }
    gmx_mtop_t mtop;
    createTwoResidueTopology(&mtop, handlesToAtomNames);
    mtop.finalize();
    AtomRange   atoms(mtop);
    std::string groupName    = "/h5md/testattrs/";
    auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), groupName.c_str()));
    {
        SCOPED_TRACE("Testing H5MD writing of array-like attributes.");
        // Caller may know this already (e.g. residue names and atom
        // names have standard maximum lengths), but if not, it's easy
        // for them to compute and hard to provide a good general helper
        // function.
        auto [maxLength, strCount] = estimateBufferSize(atoms.begin(),
                                                        atoms.end(),
                                                        [](const AtomIterator& atomIt) -> const char*
                                                        { return atomIt->atomName(); });
        EXPECT_EQ(maxLength, 3);

        // First, prepare storage
        std::vector<char> packingBuffer;
        // The caller already knows the atom count and typically the max
        // length of an atom name, so reserve storage appropriately
        packingBuffer.resize(mtop.natoms * (maxLength + 1));

        // Use utility function to do the thing that is specific to
        // how H5md wants to pack these string-vector attributes
        packingBuffer = packBufferViaIterator(atoms.begin(),
                                              atoms.end(),
                                              std::move(packingBuffer),
                                              maxLength,
                                              [](const AtomIterator& atomIt) -> const char*
                                              { return atomIt->atomName(); });
        // Do the low-level H5md things with the packed buffer
        setStringAttributeByBuffer(group, "atom_names", mtop.natoms, maxLength, packingBuffer);
    }

    {
        SCOPED_TRACE("Testing H5MD reading of array-like attributes.");

        EXPECT_EQ(getAttributeVector<std::string>(group, "atom_names"), atomNames);
    }
}


} // namespace
} // namespace test
} // namespace gmx
