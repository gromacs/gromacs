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
 * Tests for help string formatting functions and classes.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_onlinehelp
 */
#include <gtest/gtest.h>

#include "gromacs/onlinehelp/helpformat.h"

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
