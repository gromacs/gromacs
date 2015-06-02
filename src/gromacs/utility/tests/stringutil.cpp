/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
 * Tests for string utility functions and classes.
 *
 * For development, the tests can be run with a '-stdout' command-line option
 * to print out the help to stdout instead of using the XML reference
 * framework.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/stringutil.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"

namespace
{

/********************************************************************
 * Tests for simple string utilities
 */

TEST(StringUtilityTest, StartsWith)
{
    EXPECT_TRUE(gmx::startsWith("foobar", "foo"));
    EXPECT_TRUE(gmx::startsWith("foobar", ""));
    EXPECT_TRUE(gmx::startsWith("", ""));
    EXPECT_FALSE(gmx::startsWith("", "foobar"));
    EXPECT_FALSE(gmx::startsWith("foo", "foobar"));
    EXPECT_FALSE(gmx::startsWith("foobar", "oob"));
    EXPECT_TRUE(gmx::startsWith(std::string("foobar"), "foo"));
    EXPECT_TRUE(gmx::startsWith(std::string("foobar"), ""));
    EXPECT_TRUE(gmx::startsWith(std::string(""), ""));
    EXPECT_FALSE(gmx::startsWith(std::string(""), "foobar"));
    EXPECT_FALSE(gmx::startsWith(std::string("foo"), "foobar"));
    EXPECT_FALSE(gmx::startsWith(std::string("foobar"), "oob"));
}

TEST(StringUtilityTest, EndsWith)
{
    EXPECT_TRUE(gmx::endsWith("foobar", "bar"));
    EXPECT_TRUE(gmx::endsWith("foobar", NULL));
    EXPECT_TRUE(gmx::endsWith("foobar", ""));
    EXPECT_TRUE(gmx::endsWith("", ""));
    EXPECT_FALSE(gmx::endsWith("", "foobar"));
    EXPECT_FALSE(gmx::endsWith("foobar", "bbar"));
    EXPECT_FALSE(gmx::endsWith("foobar", "barr"));
    EXPECT_FALSE(gmx::endsWith("foobar", "foofoobar"));
}

TEST(StringUtilityTest, StripSuffixIfPresent)
{
    EXPECT_EQ("foo", gmx::stripSuffixIfPresent("foobar", "bar"));
    EXPECT_EQ("foobar", gmx::stripSuffixIfPresent("foobar", NULL));
    EXPECT_EQ("foobar", gmx::stripSuffixIfPresent("foobar", ""));
    EXPECT_EQ("foobar", gmx::stripSuffixIfPresent("foobar", "bbar"));
    EXPECT_EQ("foobar", gmx::stripSuffixIfPresent("foobar", "barr"));
    EXPECT_EQ("foobar", gmx::stripSuffixIfPresent("foobar", "foofoobar"));
}

TEST(StringUtilityTest, StripString)
{
    EXPECT_EQ("", gmx::stripString(""));
    EXPECT_EQ("foo", gmx::stripString("foo"));
    EXPECT_EQ("foo", gmx::stripString("  foo"));
    EXPECT_EQ("foo", gmx::stripString("foo "));
    EXPECT_EQ("f o o", gmx::stripString(" f o o  "));
}

TEST(StringUtilityTest, SplitString)
{
    using ::testing::ElementsAre;
    using ::testing::IsEmpty;
    using ::testing::Matcher;
    Matcher<std::vector<std::string> > matcher = ElementsAre("foo", "bar");
    EXPECT_THAT(gmx::splitString("foo bar"), matcher);
    EXPECT_THAT(gmx::splitString("  foo bar"), matcher);
    EXPECT_THAT(gmx::splitString("foo bar  "), matcher);
    EXPECT_THAT(gmx::splitString(" foo \t bar  "), matcher);
    EXPECT_THAT(gmx::splitString(""), IsEmpty());
    EXPECT_THAT(gmx::splitString("   "), IsEmpty());
}

/********************************************************************
 * Tests for formatString()
 */

TEST(FormatStringTest, HandlesBasicFormatting)
{
    EXPECT_EQ("12 abc", gmx::formatString("%d %s", 12, "abc"));
}

TEST(FormatStringTest, HandlesLongStrings)
{
    std::string longString = gmx::formatString("%*c%d", 2000, 'x', 10);
    EXPECT_EQ(2002U, longString.length());
    EXPECT_EQ("x10", longString.substr(1999));
}

/********************************************************************
 * Tests for StringFormatter
 */

TEST(StringFormatterTest, HandlesBasicFormatting)
{
    int value = 103;
    EXPECT_EQ("103", gmx::StringFormatter("%d") (value));
    EXPECT_EQ("null", gmx::StringFormatter("null") (value));
}

/********************************************************************
 * Tests for formatAndJoin
 */

TEST(formatAndJoinTest, Works)
{
    const char * const words[] = { "The", "quick", "brown", "fox" };
    EXPECT_EQ("The       .quick     .brown     .fox       ",
              gmx::formatAndJoin(gmx::ConstArrayRef<const char *>(words), ".",
                                 gmx::StringFormatter("%-10s")));

    const int values[] = { 0, 1, 4 };
    EXPECT_EQ("0,1,4", gmx::formatAndJoin(gmx::ConstArrayRef<int>(values), ",",
                                          gmx::StringFormatter("%d")));
}

/********************************************************************
 * Tests for joinStrings
 */

TEST(JoinStringsTest, Works)
{
    const char * const               words[] = { "The", "quick", "brown", "fox" };
    gmx::ConstArrayRef<const char *> refToWords(words);
    EXPECT_EQ("The; quick; brown; fox", gmx::joinStrings(refToWords.begin(), refToWords.end(), "; "));
    EXPECT_EQ("The-quick-brown-fox", gmx::joinStrings(refToWords, "-"));
    EXPECT_EQ("The-quick-brown-fox", gmx::joinStrings(words, "-"));
}

/********************************************************************
 * Tests for replaceAll() and replaceAllWords()
 */

TEST(ReplaceAllTest, HandlesEmptyStrings)
{
    EXPECT_EQ("", gmx::replaceAll("", "aaa", "bbbb"));
    EXPECT_EQ("", gmx::replaceAllWords("", "aaa", "bbbb"));
}

TEST(ReplaceAllTest, HandlesNoMatches)
{
    const std::string text("Text with no matches");
    EXPECT_EQ(text, gmx::replaceAll(text, "aaa", "bbbb"));
    EXPECT_EQ(text, gmx::replaceAllWords(text, "aaa", "bbbb"));
}

TEST(ReplaceAllTest, HandlesMatchesAtEnds)
{
    EXPECT_EQ("bbbbtext", gmx::replaceAll("aaatext", "aaa", "bbbb"));
    EXPECT_EQ("textbbbb", gmx::replaceAll("textaaa", "aaa", "bbbb"));
    EXPECT_EQ("bbbb text", gmx::replaceAllWords("aaa text", "aaa", "bbbb"));
    EXPECT_EQ("text bbbb", gmx::replaceAllWords("text aaa", "aaa", "bbbb"));
}

TEST(ReplaceAllTest, HandlesMultipleMatches)
{
    const std::string text("Text aaa with multiple aaa matches");
    EXPECT_EQ("Text bbbb with multiple bbbb matches",
              gmx::replaceAll(text, "aaa", "bbbb"));
    EXPECT_EQ("Text bbbb with multiple bbbb matches",
              gmx::replaceAllWords(text, "aaa", "bbbb"));
}

TEST(ReplaceAllTest, HandlesWordBoundaries)
{
    const std::string text("Text aaax with one word aaa match");
    EXPECT_EQ("Text aaax with one word bbbb match",
              gmx::replaceAllWords(text, "aaa", "bbbb"));
}

TEST(ReplaceAllTest, HandlesPossibleRecursiveMatches)
{
    const std::string text("Text with recursive aaabbbbbb matches");
    EXPECT_EQ("Text with recursive aaaaaabbb matches",
              gmx::replaceAll(text, "aaabbb", "aaaaaa"));
}

/********************************************************************
 * Tests for TextLineWrapper
 */

//! Simple test string for wrapping.
const char g_wrapText[] = "A quick brown fox jumps over the lazy dog";
//! Test string for wrapping with embedded line breaks.
const char g_wrapText2[] = "A quick brown fox jumps\nover the lazy dog";
//! Test string for wrapping with embedded line breaks and an empty line.
const char g_wrapText3[] = "A quick brown fox jumps\n\nover the lazy dog";
//! Test string for wrapping with a long word.
const char g_wrapTextLongWord[]
    = "A quick brown fox jumps awordthatoverflowsaline over the lazy dog";
//! Test string for wrapping with extra whitespace.
const char g_wrapTextWhitespace[] = " A quick brown   fox jumps  \n over the lazy dog";

//! Test fixture for gmx::TextLineWrapper.
typedef gmx::test::StringTestBase TextLineWrapperTest;

TEST_F(TextLineWrapperTest, HandlesEmptyStrings)
{
    gmx::TextLineWrapper wrapper;

    EXPECT_EQ("", wrapper.wrapToString(""));
    EXPECT_EQ("", wrapper.wrapToString("   "));
    EXPECT_TRUE(wrapper.wrapToVector("").empty());
    EXPECT_TRUE(wrapper.wrapToString("   ").empty());
}

TEST_F(TextLineWrapperTest, HandlesTrailingNewlines)
{
    gmx::TextLineWrapper wrapper;

    EXPECT_EQ("line", wrapper.wrapToString("line"));
    EXPECT_EQ("line\n", wrapper.wrapToString("line\n"));
    EXPECT_EQ("line\n\n", wrapper.wrapToString("line\n\n"));
    EXPECT_EQ("\n", wrapper.wrapToString("\n"));
    EXPECT_EQ("\n\n", wrapper.wrapToString("\n\n"));
    {
        std::vector<std::string> wrapped(wrapper.wrapToVector("line"));
        ASSERT_EQ(1U, wrapped.size());
        EXPECT_EQ("line", wrapped[0]);
    }
    {
        std::vector<std::string> wrapped(wrapper.wrapToVector("line\n"));
        ASSERT_EQ(1U, wrapped.size());
        EXPECT_EQ("line", wrapped[0]);
    }
    {
        std::vector<std::string> wrapped(wrapper.wrapToVector("line\n\n"));
        ASSERT_EQ(2U, wrapped.size());
        EXPECT_EQ("line", wrapped[0]);
        EXPECT_EQ("", wrapped[1]);
    }
    {
        std::vector<std::string> wrapped(wrapper.wrapToVector("\n"));
        ASSERT_EQ(1U, wrapped.size());
        EXPECT_EQ("", wrapped[0]);
    }
    {
        std::vector<std::string> wrapped(wrapper.wrapToVector("\n\n"));
        ASSERT_EQ(2U, wrapped.size());
        EXPECT_EQ("", wrapped[0]);
        EXPECT_EQ("", wrapped[1]);
    }
}

TEST_F(TextLineWrapperTest, WrapsCorrectly)
{
    gmx::TextLineWrapper wrapper;

    wrapper.settings().setLineLength(10);
    checkText(wrapper.wrapToString(g_wrapText), "WrappedAt10");
    std::vector<std::string> wrapped(wrapper.wrapToVector(g_wrapText));
    checker().checkSequence(wrapped.begin(), wrapped.end(), "WrappedToVector");
    wrapper.settings().setLineLength(13);
    checkText(wrapper.wrapToString(g_wrapText), "WrappedAt13");
    wrapper.settings().setLineLength(14);
    checkText(wrapper.wrapToString(g_wrapText), "WrappedAt14");
    checkText(wrapper.wrapToString(g_wrapTextLongWord), "WrappedWithLongWord");
}

TEST_F(TextLineWrapperTest, WrapsCorrectlyWithExistingBreaks)
{
    gmx::TextLineWrapper wrapper;

    checkText(wrapper.wrapToString(g_wrapText2), "WrappedWithNoLimit");
    wrapper.settings().setLineLength(10);
    checkText(wrapper.wrapToString(g_wrapText2), "WrappedAt10");
    wrapper.settings().setLineLength(14);
    checkText(wrapper.wrapToString(g_wrapText2), "WrappedAt14");
}

TEST_F(TextLineWrapperTest, HandlesIndent)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setIndent(2);

    checkText(wrapper.wrapToString(g_wrapText2), "WrappedWithNoLimit");
    wrapper.settings().setLineLength(16);
    checkText(wrapper.wrapToString(g_wrapText2), "WrappedAt14");
}

TEST_F(TextLineWrapperTest, HandlesIndentWithEmptyLines)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setIndent(2);

    checkText(wrapper.wrapToString(g_wrapText3), "WrappedWithNoLimit");
    wrapper.settings().setLineLength(16);
    checkText(wrapper.wrapToString(g_wrapText3), "WrappedAt14");
}

TEST_F(TextLineWrapperTest, HandlesHangingIndent)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setFirstLineIndent(2);
    wrapper.settings().setIndent(4);

    checkText(wrapper.wrapToString(g_wrapText2), "WrappedWithNoLimit");
    wrapper.settings().setLineLength(16);
    checkText(wrapper.wrapToString(g_wrapText2), "WrappedAt14/12");
}

TEST_F(TextLineWrapperTest, HandlesContinuationCharacter)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setFirstLineIndent(2);
    wrapper.settings().setIndent(4);
    wrapper.settings().setContinuationChar('\\');

    wrapper.settings().setLineLength(16);
    checkText(wrapper.wrapToString(g_wrapText2), "WrappedAt14/12");
}

TEST_F(TextLineWrapperTest, WrapsCorrectlyWithExtraWhitespace)
{
    gmx::TextLineWrapper wrapper;

    wrapper.settings().setLineLength(14);
    wrapper.settings().setStripLeadingWhitespace(true);
    checkText(wrapper.wrapToString(g_wrapTextWhitespace),
              "WrappedAt14StripLeading");
    wrapper.settings().setStripLeadingWhitespace(false);
    checkText(wrapper.wrapToString(g_wrapTextWhitespace),
              "WrappedAt14PreserveLeading");
}

} // namespace
