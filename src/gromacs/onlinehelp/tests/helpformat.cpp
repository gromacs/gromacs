/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Tests for help string formatting functions and classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "gromacs/onlinehelp/helpformat.h"

#include <gtest/gtest.h>

#include "testutils/stringtest.h"

namespace
{

//! Simple test string for wrapping.
const char g_wrapText[] = "A quick brown fox jumps over the lazy dog";
//! Test string for wrapping with embedded line breaks.
const char g_wrapText2[] = "A quick brown fox jumps\nover the lazy dog";

/********************************************************************
 * Tests for TextTableFormatter
 */

class TextTableFormatterTest : public gmx::test::StringTestBase
{
    public:
        void setupStandardColumns();

        gmx::TextTableFormatter formatter_;
};

void TextTableFormatterTest::setupStandardColumns()
{
    formatter_.addColumn("Col1", 4, false);
    formatter_.addColumn("Col2", 4, false);
    formatter_.addColumn("Col3Wrap", 14, true);
    formatter_.addColumn("Col4Wrap", 14, true);
}

TEST_F(TextTableFormatterTest, HandlesBasicCase)
{
    setupStandardColumns();
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedTable");
}

TEST_F(TextTableFormatterTest, HandlesEmptyColumnTitles)
{
    formatter_.addColumn(NULL, 4, false);
    formatter_.addColumn("", 4, false);
    formatter_.addColumn(NULL, 14, true);
    formatter_.addColumn("", 14, true);

    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedTable");
}

TEST_F(TextTableFormatterTest, HandlesIndentation)
{
    setupStandardColumns();
    formatter_.setFirstColumnIndent(4);
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedTable");
}

TEST_F(TextTableFormatterTest, HandlesOverflowingLines)
{
    setupStandardColumns();
    formatter_.clear();
    formatter_.addColumnLine(0, "foobar");
    formatter_.addColumnLine(1, "barfoo");
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedTable");
    formatter_.clear();
    formatter_.addColumnLine(0, "foobar");
    formatter_.addColumnLine(1, "barfoo");
    formatter_.setColumnFirstLineOffset(1, 1);
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedRow2");
    formatter_.clear();
    formatter_.addColumnLine(0, "foobar");
    formatter_.addColumnLine(1, "barfoo");
    formatter_.setColumnFirstLineOffset(1, 1);
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.setColumnFirstLineOffset(2, 2);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedRow3");
    // Test a case where a column value overflows even the next column.
    formatter_.addColumnLine(0, "foobarfoobar");
    formatter_.addColumnLine(1, "barfoobarfoo");
    formatter_.addColumnLine(2, g_wrapText);
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedRow4");
}

TEST_F(TextTableFormatterTest, HandlesLastColumnFolding)
{
    setupStandardColumns();
    formatter_.setFoldLastColumnToNextLine(4);
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, "text");
    formatter_.addColumnLine(3, g_wrapText);
    checkText(formatter_.formatRow(), "FormattedTable");
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, "text");
    formatter_.addColumnLine(3, g_wrapText2);
    checkText(formatter_.formatRow(), "FormattedRow2");
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, "text");
    formatter_.addColumnLine(3, "A quick brown fox jumps");
    checkText(formatter_.formatRow(), "FormattedRow3");
}

TEST_F(TextTableFormatterTest, HandlesEmptyColumns)
{
    setupStandardColumns();
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(3, "Text");
    checkText(formatter_.formatRow(), "FormattedTable");
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.setColumnFirstLineOffset(2, 1);
    formatter_.addColumnLine(3, "Text");
    checkText(formatter_.formatRow(), "FormattedRow2");
    formatter_.clear();
    formatter_.addColumnLine(0, "foo");
    formatter_.addColumnLine(1, "bar");
    formatter_.addColumnLine(2, "");
    formatter_.setColumnFirstLineOffset(2, 1);
    formatter_.addColumnLine(3, "Text");
    checkText(formatter_.formatRow(), "FormattedRow3");
}

} // namespace
