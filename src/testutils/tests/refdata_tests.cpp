/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests utilities for test reference data.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/refdata.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/variant.h"

#include "testutils/testasserts.h"
#include "testutils/testexceptions.h"

namespace
{

using gmx::test::TestReferenceData;
using gmx::test::TestReferenceChecker;

TEST(ReferenceDataTest, HandlesSimpleData)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkBoolean(true, "bool");
        checker.checkInteger(1, "int");
        checker.checkInt64(1ULL<<42, "int64");
        checker.checkUInt64(1ULL<<42, "uint64");
        checker.checkDouble(0.5, "real");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkBoolean(true, "bool");
        checker.checkInteger(1, "int");
        checker.checkInt64(1ULL<<42, "int64");
        checker.checkUInt64(1ULL<<42, "uint64");
        checker.checkDouble(0.5, "real");
        checker.checkString("Test", "string");
    }
}

TEST(ReferenceDataTest, HandlesFloatingPointData)
{
    const float  floatValue  = 4.0f/3.0f;
    const double doubleValue = 4.0/3.0;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkDouble(doubleValue, "double");
        checker.checkReal(doubleValue, "real");
        checker.checkFloat(floatValue, "float");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkDouble(doubleValue, "double");
        checker.checkReal(floatValue, "real");
        checker.checkReal(doubleValue, "real");
        checker.checkFloat(floatValue, "float");
    }
}

TEST(ReferenceDataTest, HandlesPresenceChecks)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_TRUE(checker.checkPresent(true, "present"));
        checker.checkInteger(1, "present");
        EXPECT_FALSE(checker.checkPresent(false, "absent"));
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        // Assigned to avoid warnings about potentially uninitialized value.
        bool                 bRet = true;
        EXPECT_TRUE(checker.checkPresent(true, "present"));
        checker.checkInteger(1, "present");
        EXPECT_NONFATAL_FAILURE(bRet = checker.checkPresent(false, "present"), "");
        EXPECT_FALSE(bRet);
        EXPECT_NONFATAL_FAILURE(bRet = checker.checkPresent(true, "absent"), "");
        EXPECT_FALSE(bRet);
        EXPECT_FALSE(checker.checkPresent(false, "absent"));
    }
}


TEST(ReferenceDataTest, HandlesStringBlockData)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkTextBlock("Line1\nLine2\n", "block");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkTextBlock("Line1\nLine2\n", "block");
        checker.checkString("Line1\nLine2\n", "block");
        checker.checkTextBlock("Test", "string");
    }
}


TEST(ReferenceDataTest, HandlesVectorData)
{
    int    veci[3] = { -1, 3, 5 };
    float  vecf[3] = { -2.3f, 1.43f, 2.5f };
    double vecd[3] = { -2.3, 1.43, 2.5 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVector(veci, "ivec");
        checker.checkVector(vecf, "fvec");
        checker.checkVector(vecd, "dvec");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVector(veci, "ivec");
        checker.checkVector(vecf, "fvec");
        checker.checkVector(vecd, "dvec");
    }
}


TEST(ReferenceDataTest, HandlesSequenceData)
{
    const int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkSequenceArray(5, seq, "seq");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkSequenceArray(5, seq, "seq");
    }
}

//! Helper typedef
typedef double dvec[3];
//! Helper function for HandlesSequenceOfCustomData
void checkCustomVector(TestReferenceChecker *checker, const dvec &value)
{
    checker->checkVector(value, nullptr);
}

TEST(ReferenceDataTest, HandlesSequenceOfCustomData)
{
    const dvec seq[] = { {-3, 4, 5}, {-2.3, 5, 0} };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkSequence(std::begin(seq), std::end(seq), "seq", checkCustomVector);
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkSequence(std::begin(seq), std::end(seq), "seq", checkCustomVector);
    }
}


TEST(ReferenceDataTest, HandlesIncorrectData)
{
    int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkDouble(0.5, "real");
        checker.checkString("Test", "string");
        checker.checkSequenceArray(5, seq, "seq");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkInteger(2, "int"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkDouble(0.3, "real"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test2", "string"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(4, seq, "seq"), "");
        seq[0] = 2;
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(5, seq, "seq"), "");
    }
}

TEST(ReferenceDataTest, HandlesIncorrectDataType)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkCompound("Compound", "compound");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkString("1", "int"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkCompound("OtherCompound", "compound"), "");
    }
}


TEST(ReferenceDataTest, HandlesMissingData)
{
    const int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkSequenceArray(5, seq, "seq");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkInteger(1, "missing"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(5, seq, "missing"), "");
        // Needed to not make the test fail because of unused "int" and "seq".
        EXPECT_NONFATAL_FAILURE(checker.checkUnusedEntries(), "");
    }
}


TEST(ReferenceDataTest, HandlesUncheckedData)
{
    const int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkSequenceArray(5, seq, "seq");
        checker.checkUnusedEntries();
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        EXPECT_NONFATAL_FAILURE(checker.checkUnusedEntries(), "");
    }
}


TEST(ReferenceDataTest, HandlesUncheckedDataInSequence)
{
    const int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkSequenceArray(5, seq, "seq");
        checker.checkUnusedEntries();
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(3, seq, "seq"), "");
        // It might be nicer to not report the unused sequence entries also
        // here, but both behaviors are quite OK.
        EXPECT_NONFATAL_FAILURE(checker.checkUnusedEntries(), "");
    }
}


TEST(ReferenceDataTest, HandlesUncheckedDataInCompound)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        TestReferenceChecker compound(checker.checkCompound("Compound", "Compound"));
        compound.checkInteger(1, "int1");
        compound.checkInteger(2, "int2");
        compound.checkUnusedEntries();
        checker.checkInteger(1, "int");
        checker.checkUnusedEntries();
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        TestReferenceChecker compound(checker.checkCompound("Compound", "Compound"));
        compound.checkInteger(1, "int1");
        EXPECT_NONFATAL_FAILURE(compound.checkUnusedEntries(), "");
        checker.checkInteger(1, "int");
        checker.checkUnusedEntries();
    }
}


TEST(ReferenceDataTest, HandlesVariants)
{
    using gmx::Variant;
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVariant(Variant::create<bool>(true), "bool");
        checker.checkVariant(Variant::create<int>(1), "int");
        checker.checkVariant(Variant::create<double>(3.5), "real");
        checker.checkVariant(Variant::create<std::string>("foo"), "str");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVariant(Variant::create<bool>(true), "bool");
        checker.checkVariant(Variant::create<int>(1), "int");
        checker.checkVariant(Variant::create<double>(3.5), "real");
        checker.checkVariant(Variant::create<std::string>("foo"), "str");
    }
}

//! Helper for building a KeyValueTree for testing.
gmx::KeyValueTreeObject buildKeyValueTree(bool full)
{
    gmx::KeyValueTreeBuilder builder;
    auto                     root = builder.rootObject();
    auto                     obj  = root.addObject("o");
    obj.addValue<int>("i", 1);
    if (full)
    {
        obj.addValue<std::string>("s", "x");
    }
    auto arr  = root.addUniformArray<int>("a");
    arr.addValue(2);
    arr.addValue(3);
    root.addValue<std::string>("s", "y");
    return builder.build();
}


TEST(ReferenceDataTest, HandlesKeyValueTree)
{
    gmx::KeyValueTreeObject tree = buildKeyValueTree(true);
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(tree, "tree");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(tree, "tree");
    }
}


TEST(ReferenceDataTest, HandlesKeyValueTreeExtraKey)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(buildKeyValueTree(false), "tree");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkKeyValueTreeObject(buildKeyValueTree(true), "tree"), "");
    }
}


TEST(ReferenceDataTest, HandlesKeyValueTreeMissingKey)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkKeyValueTreeObject(buildKeyValueTree(true), "tree");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkKeyValueTreeObject(buildKeyValueTree(false), "tree"), "");
    }
}


TEST(ReferenceDataTest, HandlesVariantsWithIncorrectValue)
{
    using gmx::Variant;
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVariant(Variant::create<bool>(true), "bool");
        checker.checkVariant(Variant::create<int>(1), "int");
        checker.checkVariant(Variant::create<double>(3.5), "real");
        checker.checkVariant(Variant::create<std::string>("foo"), "str");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<bool>(false), "bool"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<int>(2), "int"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<double>(2.5), "real"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<std::string>("bar"), "str"), "");
    }
}


TEST(ReferenceDataTest, HandlesVariantsWithIncorrectType)
{
    using gmx::Variant;
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkVariant(Variant::create<bool>(true), "bool");
        checker.checkVariant(Variant::create<int>(1), "int");
        checker.checkVariant(Variant::create<double>(3.5), "real");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<int>(1), "bool"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<bool>(true), "int"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkVariant(Variant::create<int>(2), "real"), "");
    }
}


TEST(ReferenceDataTest, HandlesMissingReferenceDataFile)
{
    const int seq[5] = { -1, 3, 5, 2, 4 };

    EXPECT_NONFATAL_FAILURE({
                                TestReferenceData data(gmx::test::erefdataCompare);
                                TestReferenceChecker checker(data.rootChecker());
                                checker.checkInteger(1, "int");
                                checker.checkDouble(0.5, "real");
                                checker.checkString("Test", "string");
                                checker.checkSequenceArray(5, seq, "seq");
                            }, "");
}


TEST(ReferenceDataTest, HandlesSpecialCharactersInStrings)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        // Note that '\r' is not handled correctly in string or
        // stringblock (see the TODO in createElementContents), so
        // don't try to test it
        checker.checkString("\"<'>\n&\\/;", "string");
        checker.checkTextBlock("\"<'>\n&\\/;", "stringblock");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("\"<'>\n&\\/;", "string");
        checker.checkTextBlock("\"<'>\n&\\/;", "stringblock");
    }
}

TEST(ReferenceDataTest, HandlesStringsWithTextAndWhitespace)
{
    const char *strings[] = { "  test", "test  ", "  test  ", "the test", "\ntest", "\n\ntest", "test\n", "test\n\n" };
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        for (const auto &s : strings)
        {
            checker.checkString(s, nullptr);
            checker.checkTextBlock(s, nullptr);
        }
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        for (const auto &s : strings)
        {
            checker.checkString(s, nullptr);
            checker.checkTextBlock(s, nullptr);
        }
    }
}

TEST(ReferenceDataTest, HandlesEmptyStrings)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("", "Empty");
        // GROMACS cannot use an empty line in a reference data String
        // until https://github.com/leethomason/tinyxml2/issues/432 is
        // resolved.
        EXPECT_THROW_GMX(checker.checkString("\n", "EmptyLine"), gmx::test::TestException);
        checker.checkTextBlock("", "EmptyBlock");
        checker.checkTextBlock("\n", "EmptyLineBlock");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("", "Empty");
        EXPECT_THROW_GMX(checker.checkString("\n", "EmptyLine"), gmx::test::TestException);
        checker.checkTextBlock("", "EmptyBlock");
        checker.checkTextBlock("\n", "EmptyLineBlock");
    }
}


TEST(ReferenceDataTest, HandlesEmbeddedCdataEndTagInTextBlock)
{
    /* stringblocks are implemented as CDATA fields, and the first
       appearance of "]]>" always terminates the CDATA field. If a
       string to be stored in a stringblock would contain such text,
       then a quality XML writer would escape the text somehow, and
       read it back in a matching way. This test verifies that the
       overall implementation copes with this issue. (GROMACS tests
       don't actually depend on this behaviour, but it might be nice
       to have / know about.) */
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkTextBlock(" ]]> ", "stringblock");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkTextBlock(" ]]> ", "stringblock");
    }
}


TEST(ReferenceDataTest, HandlesSequenceItemIndices)
{
    int seq[5] = { -1, 3, 5, 2, 4 };

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkSequenceArray(5, seq, "seq");
        checker.checkSequenceArray(5, seq, "seq2");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        seq[0] = 2;
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(5, seq, "seq"), "seq/[0]");
        seq[0] = -1;
        seq[3] = 0;
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(5, seq, "seq"), "seq/[3]");
        EXPECT_NONFATAL_FAILURE(checker.checkSequenceArray(5, seq, "seq2"), "seq2/[3]");
    }
}


TEST(ReferenceDataTest, HandlesMultipleChecksAgainstSameData)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", "string");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test2", "string"), "");
        checker.checkTextBlock("TestString", "stringblock");
        EXPECT_NONFATAL_FAILURE(checker.checkTextBlock("TestString2", "stringblock"), "");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", "string");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test2", "string"), "");
        checker.checkTextBlock("TestString", "stringblock");
        EXPECT_NONFATAL_FAILURE(checker.checkTextBlock("TestString2", "stringblock"), "");
    }
}


TEST(ReferenceDataTest, HandlesMultipleNullIds)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", nullptr);
        checker.checkString("Test2", nullptr);
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", nullptr);
        checker.checkString("Test2", nullptr);
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test", nullptr), "");
    }
}


TEST(ReferenceDataTest, HandlesMultipleComparisonsAgainstNullIds)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int1");
        checker.checkString("Test", nullptr);
        checker.checkString("Test2", nullptr);
        checker.checkInteger(2, "int2");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test3", nullptr), "");
        checker.checkString("Test2", nullptr);
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int1");
        checker.checkString("Test", nullptr);
        checker.checkString("Test2", nullptr);
        checker.checkInteger(2, "int2");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test3", nullptr), "");
        checker.checkInteger(1, "int1");
        checker.checkString("Test", nullptr);
        checker.checkString("Test2", nullptr);
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test", nullptr), "");
        checker.checkInteger(2, "int2");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test3", nullptr), "");
    }
}


TEST(ReferenceDataTest, HandlesReadingValues)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkUChar('A', "char");
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        EXPECT_EQ('A', checker.readUChar("char"));
        EXPECT_EQ(1, checker.readInteger("int"));
        EXPECT_EQ("Test", checker.readString("string"));
    }
}


TEST(ReferenceDataTest, HandlesUpdateChangedWithoutChanges)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
        TestReferenceChecker compound(checker.checkCompound("Compound", "Compound"));
        compound.checkInteger(2, "int");
    }
    {
        TestReferenceData    data(gmx::test::erefdataUpdateChanged);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
        TestReferenceChecker compound(checker.checkCompound("Compound", "Compound"));
        compound.checkInteger(2, "int");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
        TestReferenceChecker compound(checker.checkCompound("Compound", "Compound"));
        compound.checkInteger(2, "int");
    }
}

TEST(ReferenceDataTest, HandlesUpdateChangedWithValueChanges)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataUpdateChanged);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(2, "int");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(2, "int");
        checker.checkString("Test", "string");
    }
}

TEST(ReferenceDataTest, HandlesUpdateChangedWithTypeChanges)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "foo");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataUpdateChanged);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("foo", "foo");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("foo", "foo");
        checker.checkString("Test", "string");
    }
}

TEST(ReferenceDataTest, HandlesUpdateChangedWithCompoundChanges)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "1");
        TestReferenceChecker compound(checker.checkCompound("Compound", "2"));
        compound.checkInteger(2, "int");
    }
    {
        TestReferenceData    data(gmx::test::erefdataUpdateChanged);
        TestReferenceChecker checker(data.rootChecker());
        TestReferenceChecker compound(checker.checkCompound("Compound", "1"));
        compound.checkInteger(2, "int");
        checker.checkString("Test", "2");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        TestReferenceChecker compound(checker.checkCompound("Compound", "1"));
        compound.checkInteger(2, "int");
        checker.checkString("Test", "2");
    }
}

TEST(ReferenceDataTest, HandlesUpdateChangedWithRemovedEntries)
{
    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(1, "int");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataUpdateChanged);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(2, "int");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkInteger(2, "int");
    }
}

} // namespace
