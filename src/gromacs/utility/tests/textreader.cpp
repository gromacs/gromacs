/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for gmx::TextReader.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/textreader.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/stringstream.h"
//#include "gromacs/utility/stringutil.h"

namespace
{

TEST(TextReaderTest, ReadsTrimmedLines)
{
    std::vector<std::string> lines =
    {
        "expected text",
        "expected text \t "
    };
    gmx::StringInputStream   stream(lines);
    gmx::TextReader          reader(&stream);
    std::string              line;
    for (std::size_t i = 0; i != lines.size(); ++i)
    {
        EXPECT_TRUE(reader.readLineTrimmed(&line));
        EXPECT_STREQ("expected text", line.c_str());
    }
    EXPECT_FALSE(reader.readLineTrimmed(&line));
}

TEST(TextReaderTest, ReadsLinesWithoutCommentsAndTrimmed)
{
    std::vector<std::string> lines =
    {
        "expected text",
        "expected text \t ",
        "expected text#",
        "expected text\t #",
        "expected text   # not expected"
    };
    gmx::StringInputStream   stream(lines);
    gmx::TextReader          reader(&stream);
    reader.setCommentChar('#');
    std::string              line;
    for (std::size_t i = 0; i != lines.size(); ++i)
    {
        EXPECT_TRUE(reader.readLineWithoutCommentsAndTrimmed(&line));
        EXPECT_STREQ("expected text", line.c_str());
    }
    EXPECT_FALSE(reader.readLineWithoutCommentsAndTrimmed(&line));
}

TEST(TextReaderTest, ReadsLinesWithoutCommentsAndTrimmedThatParseEmpty)
{
    std::vector<std::string> lines =
    {
        "",
        " \t ",
        "#",
        "\t #",
        "   # not expected"
    };
    gmx::StringInputStream   stream(lines);
    gmx::TextReader          reader(&stream);
    reader.setCommentChar('#');
    std::string              line;
    for (std::size_t i = 0; i != lines.size(); ++i)
    {
        EXPECT_TRUE(reader.readLineWithoutCommentsAndTrimmed(&line));
        EXPECT_STREQ("", line.c_str());
    }
    EXPECT_FALSE(reader.readLineWithoutCommentsAndTrimmed(&line));
}

} // namespace
