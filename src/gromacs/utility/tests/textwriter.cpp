/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Tests for gmx::TextWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/textwriter.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/stringtest.h"

namespace gmx
{
namespace test
{
namespace
{

class TextWriterTest : public gmx::test::StringTestBase
{
public:
    TextWriterTest() : writer_(&stream_) {}

    void checkOutput() { checkText(stream_.toString(), "Output"); }

    gmx::StringOutputStream stream_;
    gmx::TextWriter         writer_;
};

TEST_F(TextWriterTest, WritesLines)
{
    writer_.writeLine("Explicit newline\n");
    writer_.writeLine("Implicit newline");
    writer_.writeLine(std::string("Explicit newline\n"));
    writer_.writeLine(std::string("Implicit newline"));
    writer_.writeLine();
    checkOutput();
}

TEST_F(TextWriterTest, WritesLinesInParts)
{
    writer_.writeString("Partial ");
    writer_.writeString("spaced");
    writer_.writeString(" line");
    writer_.writeLine();
    writer_.writeString(std::string("Partial "));
    writer_.writeString(std::string("spaced"));
    writer_.writeString(std::string(" line"));
    writer_.writeLine();
    checkOutput();
}

TEST_F(TextWriterTest, WritesWrappedLines)
{
    writer_.wrapperSettings().setIndent(2);
    writer_.wrapperSettings().setLineLength(15);
    writer_.writeLine("Wrapped and indented text");
    writer_.writeLine(std::string("Wrapped and indented text"));
    writer_.writeLine();
    checkOutput();
}

TEST_F(TextWriterTest, WritesLinesInPartsWithWrapper)
{
    writer_.wrapperSettings().setLineLength(50);
    writer_.writeString("Partial ");
    writer_.writeString("spaced");
    writer_.writeString(" line");
    writer_.writeLine();
    writer_.writeString(std::string("Partial "));
    writer_.writeString(std::string("spaced"));
    writer_.writeString(std::string(" line"));
    writer_.writeLine();
    checkOutput();
}

TEST_F(TextWriterTest, TracksNewlines)
{
    writer_.ensureLineBreak();
    writer_.ensureEmptyLine();
    writer_.writeString("First line");
    writer_.ensureLineBreak();
    writer_.ensureLineBreak();
    writer_.writeString("Second line");
    writer_.ensureEmptyLine();
    writer_.writeLine("Third line");
    writer_.ensureEmptyLine();
    writer_.writeString(std::string("Fourth line\n"));
    writer_.ensureLineBreak();
    writer_.writeString(std::string("Fifth line\n\n"));
    writer_.ensureEmptyLine();
    writer_.writeString(std::string("Sixth line"));
    writer_.ensureEmptyLine();
    checkOutput();
}

TEST_F(TextWriterTest, PreservesTrailingWhitespace)
{
    writer_.writeString("Line   ");
    writer_.writeLine();
    writer_.writeString(std::string("Line   "));
    writer_.writeLine();
    writer_.writeLine("Line   ");
    writer_.writeLine(std::string("Line   "));
    writer_.writeString("Line   \n");
    writer_.writeString(std::string("Line   \n"));
    checkOutput();
}

} // namespace
} // namespace test
} // namespace gmx
