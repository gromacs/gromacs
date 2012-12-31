/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Tests for string utility functions and classes.
 *
 * For development, the tests can be run with a '-stdout' command-line option
 * to print out the help to stdout instead of using the XML reference
 * framework.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"

namespace
{

/********************************************************************
 * Tests for simple string utilities
 */

TEST(StringUtilityTest, StartsWithWorks)
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

TEST(StringUtilityTest, EndsWithWorks)
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
 * Tests for concatenateStrings()
 */

//! Test fixture for gmx::concatenateStrings().
typedef gmx::test::StringTestBase ConcatenateStringsTest;

TEST_F(ConcatenateStringsTest, HandlesDifferentStringEndings)
{
    static const char * const strings[] =
    {
        "First string",
        "Second string ",
        "Third string\n",
        "Fourth string",
        ""
    };
    checkText(gmx::concatenateStrings(strings), "CombinedStrings");
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
