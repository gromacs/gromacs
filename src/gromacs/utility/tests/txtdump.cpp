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
 * Tests for text dumping utilities
 *
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/txtdump.h"

#include <cstdio>

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/posixmemstream.h"
#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{
namespace
{

using DumpingTextTest = StringTestBase;

TEST_F(DumpingTextTest, LegacyIndentingWorks)
{
    {
        PosixMemstream stream;
        const int      amountToIndent = 0;
        EXPECT_EQ(pr_indent(stream.stream(), amountToIndent), amountToIndent);
        if (stream.canCheckBufferContents())
        {
            checkText(stream.toString(), "indent 0");
        }
    }
    {
        PosixMemstream stream;
        const int      amountToIndent = 2;
        EXPECT_EQ(pr_indent(stream.stream(), amountToIndent), amountToIndent);
        if (stream.canCheckBufferContents())
        {
            checkText(stream.toString(), "indent 2");
        }
    }
}

TEST_F(DumpingTextTest, ReportsWhenNotAvailableIdentically)
{
    const int   amountToIndent = 2;
    const char* title          = "The title";
    const char* test           = nullptr;
    {
        PosixMemstream stream;
        EXPECT_EQ(available(stream.stream(), test, amountToIndent, title), false);
        if (stream.canCheckBufferContents())
        {
            checkText(stream.toString(), "output");
        }
    }
    {
        StringOutputStream stringStream;
        TextWriter         writer(&stringStream);
        ScopedIndenter     indenter = writer.addScopedIndentation(amountToIndent);
        EXPECT_EQ(checkIfAvailable(&writer, title, test), false);
        checkText(stringStream.toString(), "output");
    }
}

TEST_F(DumpingTextTest, QuietWhenAvailableIdentically)
{
    const int   amountToIndent = 2;
    const char* title          = "The title";
    const char* test           = "test";
    {
        PosixMemstream stream;
        EXPECT_EQ(available(stream.stream(), test, amountToIndent, title), true);
        if (stream.canCheckBufferContents())
        {
            checkText(stream.toString(), "output");
        }
    }
    {
        StringOutputStream stringStream;
        TextWriter         writer(&stringStream);
        ScopedIndenter     indenter = writer.addScopedIndentation(amountToIndent);
        EXPECT_EQ(checkIfAvailable(&writer, title, test), true);
        checkText(stringStream.toString(), "output");
    }
}

TEST_F(DumpingTextTest, PrintsTitleAndSizeIdentically)
{
    const int   amountToIndent = 2;
    const char* title          = "The title";
    const int   numberToPrint  = 4;
    {
        PosixMemstream stream;
        EXPECT_EQ(pr_title_n(stream.stream(), amountToIndent, title, numberToPrint), amountToIndent + 3);
        if (stream.canCheckBufferContents())
        {
            checkText(stream.toString(), "output");
        }
    }
    {
        StringOutputStream stringStream;
        TextWriter         writer(&stringStream);
        ScopedIndenter     indenter = writer.addScopedIndentation(amountToIndent);
        printTitleAndSize(&writer, title, numberToPrint);
        checkText(stringStream.toString(), "output");
    }
}

} // namespace
} // namespace test
} // namespace gmx
