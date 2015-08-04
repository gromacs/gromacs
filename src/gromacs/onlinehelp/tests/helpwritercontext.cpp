/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Tests for help text formatting.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "gromacs/onlinehelp/helpwritercontext.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/stringtest.h"

namespace
{

/********************************************************************
 * Tests for HelpWriterContext
 */

class HelpWriterContextTest : public gmx::test::StringTestBase
{
    public:
        void testFormatting(const std::string     &text,
                            gmx::HelpOutputFormat  format,
                            const char            *id)
        {
            gmx::HelpWriterContext context(NULL, format);
            std::string            result
                = context.substituteMarkupAndWrapToString(settings_, text);
            if (id == NULL)
            {
                switch (format)
                {
                    case gmx::eHelpOutputFormat_Console:
                        id = "Console";
                        break;
                    case gmx::eHelpOutputFormat_Rst:
                        id = "reStructuredText";
                        break;
                    default:
                        GMX_RELEASE_ASSERT(false, "Help output format testing not implemented");
                }
            }
            checkText(result, id);
        }
        void testFormatting(const gmx::ConstArrayRef<const char *> &text)
        {
            std::string testText = gmx::joinStrings(text, "\n");
            testFormatting(testText, gmx::eHelpOutputFormat_Console, NULL);
            testFormatting(testText, gmx::eHelpOutputFormat_Rst, NULL);
        }

        gmx::TextLineWrapperSettings settings_;
};

TEST_F(HelpWriterContextTest, FormatsParagraphs)
{
    const char *const text[] = {
        "A quick",
        "brown fox",
        "jumps over a lazy dog.[PAR]",
        "A quick brown fox",
        "jumps over",
        "a lazy dog."
    };
    settings_.setLineLength(13);
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsRstStyleParagraphs)
{
    const char *const text[] = {
        "A quick",
        "brown fox",
        "jumps over a lazy dog.",
        "",
        "A quick brown fox",
        "jumps over",
        "a lazy dog."
    };
    settings_.setLineLength(13);
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, CleansUpExtraWhitespace)
{
    const char *const text[] = {
        "",
        "A quick ",
        "brown fox ",
        "jumps over a lazy dog.",
        "[PAR]",
        "A quick brown fox",
        "jumps over",
        "a lazy dog.[PAR]"
    };
    settings_.setLineLength(13);
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsLiteralText)
{
    const char *const text[] = {
        "Sample paragraph::",
        "",
        "    literal block",
        "    another line",
        "",
        "Normal paragraph",
        "with wrapping"
    };
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsLiteralTextAtBeginning)
{
    const char *const text[] = {
        "::",
        "",
        "    literal block",
        "    another line",
        "",
        "Normal paragraph"
    };
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsLiteralTextWithIndentation)
{
    const char *const text[] = {
        "Sample paragraph::",
        "",
        "    literal block",
        "      another indented line",
        "",
        "Normal paragraph",
        "with wrapping"
    };
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsBulletList)
{
    const char *const text[] = {
        "Sample list:",
        "",
        "* first item",
        "* second item that",
        "  spans multiple lines",
        "* third item that has a single long line",
        "",
        "Normal paragraph"
    };
    // Wrapping to rst with a fixed line length does not currently work
    // correctly, but it is not used, either.
    settings_.setLineLength(15);
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsEnumeratedList)
{
    const char *const text[] = {
        "Sample list:",
        "",
        "1. first item",
        "2. second item that",
        "   spans multiple lines",
        "3. third item that has a single long line",
        "",
        "9.  Item with extra indentation",
        "10. Double digit item",
        "",
        "Normal paragraph"
    };
    // Wrapping to rst with a fixed line length does not currently work
    // correctly, but it is not used, either.
    settings_.setLineLength(15);
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsSimpleTable)
{
    const char *const text[] = {
        "Simple table:",
        "",
        "============  =============",
        "First column  Second header",
        "============  =============",
        "text          text",
        "============  =============",
        "",
        "Normal paragraph",
        "again."
    };
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsGridTable)
{
    const char *const text[] = {
        "Grid table:",
        "",
        "+--------------+---------------+",
        "| First column | Second header |",
        "+--------------+---------------+",
        "| text         | text          |",
        "+--------------+---------------+",
        "",
        "Normal paragraph",
        "again."
    };
    testFormatting(text);
}

TEST_F(HelpWriterContextTest, FormatsTitles)
{
    const char *const text[] = {
        "Title",
        "=====",
        "Some text without spacing",
        "",
        "Subtitle",
        "++++++++",
        "",
        "More text",
    };
    testFormatting(text);
}

} // namespace
