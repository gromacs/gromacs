/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

#include <vector>

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

namespace
{

TEST(ReferenceDataTest, HandlesSimpleData)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkBoolean(true, "int");
        checker.checkInteger(1, "int");
        checker.checkInt64(1ULL<<42, "int64");
        checker.checkUInt64(1ULL<<42, "uint64");
        checker.checkDouble(0.5, "real");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkBoolean(true, "int");
        checker.checkInteger(1, "int");
        checker.checkInt64(1ULL<<42, "int64");
        checker.checkUInt64(1ULL<<42, "uint64");
        checker.checkDouble(0.5, "real");
        checker.checkString("Test", "string");
    }
}

TEST(ReferenceDataTest, HandlesFloatingPointData)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

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
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkStringBlock("Line1\nLine2\n", "block");
        checker.checkString("Test", "string");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkStringBlock("Line1\nLine2\n", "block");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Line1\nLine2\n", "block"), "");
        EXPECT_NONFATAL_FAILURE(checker.checkStringBlock("Test", "string"), "");
    }
}


TEST(ReferenceDataTest, HandlesVectorData)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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


TEST(ReferenceDataTest, HandlesIncorrectData)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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


TEST(ReferenceDataTest, HandlesMissingData)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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
    }
}


TEST(ReferenceDataTest, HandlesMissingReferenceDataFile)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("\"<'>\n \r &\\/;", "string");
        // \r is not handled correctly
        checker.checkStringBlock("\"<'>\n ]]> &\\/;", "stringblock");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("\"<'>\n \r &\\/;", "string");
        checker.checkStringBlock("\"<'>\n ]]> &\\/;", "stringblock");
    }
}


TEST(ReferenceDataTest, HandlesSequenceItemIndices)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;
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
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", "string");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test2", "string"), "");
        checker.checkStringBlock("TestString", "stringblock");
        EXPECT_NONFATAL_FAILURE(checker.checkStringBlock("TestString2", "stringblock"), "");
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", "string");
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test2", "string"), "");
        checker.checkStringBlock("TestString", "stringblock");
        EXPECT_NONFATAL_FAILURE(checker.checkStringBlock("TestString2", "stringblock"), "");
    }
}


TEST(ReferenceDataTest, HandlesMultipleNullIds)
{
    using gmx::test::TestReferenceData;
    using gmx::test::TestReferenceChecker;

    {
        TestReferenceData    data(gmx::test::erefdataUpdateAll);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", NULL);
        checker.checkString("Test2", NULL);
    }
    {
        TestReferenceData    data(gmx::test::erefdataCompare);
        TestReferenceChecker checker(data.rootChecker());
        checker.checkString("Test", NULL);
        checker.checkString("Test2", NULL);
        EXPECT_NONFATAL_FAILURE(checker.checkString("Test", NULL), "");
    }
}

} // namespace
