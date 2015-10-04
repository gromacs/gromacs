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
 * Tests for reading xvg files
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/xvgreader.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/*! \brief Helper typedef for the exception most often thrown
 *
 * This is useful as a central point change so one can inspect the
 * text of the reason contained in the exception by hacking the test
 * to fail.
 *
 * TODO Find a nice way for a developer to do this without having to
 * change the code to inspect the reason message. */
typedef InvalidInputError ExpectedException;

//! Array of lines to fake an .xvg file
const char * xvgLines[] = {
    "# This file was created Sun Oct  4 19:42:39 2015",
    "#",
    "@    title \"dH/d\\xl\\f{} and \\xD\\f{}H\"",
    "@    xaxis  label \"Time (ps)\"",
    "@    yaxis  label \"dH/d\\xl\\f{} and \\xD\\f{}H (kJ/mol [\\xl\\f{}]\\S-1\\N)\"",
    "@TYPE xy",
    "@ subtitle \"T = 298 (K) \\xl\\f{} state 12: (fep-lambda, coul-lambda, vdw-lambda, bonded-lambda) = (0.0000, 1.0000, 0.6000, 0.6000)\"",
    "@ view 0.15, 0.15, 0.75, 0.85",
    "@ legend on",
    "@ s0 legend \"dH/d\\xl\\f{} fep-lambda = 0.0000\"",
    "@ s1 legend \"dH/d\\xl\\f{} coul-lambda = 1.0000\"",
    "@ s2 legend \"dH/d\\xl\\f{} vdw-lambda = 0.6000\"",
    "@ s3 legend \"dH/d\\xl\\f{} bonded-lambda = 0.6000\"",
    "@ s4 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.5000, 0.5000)\"",
    "@ s5 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.6000, 0.6000)\"",
    "@ s6 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.6500, 0.6500)\"",
    "@ s7 legend \"pV (kJ/mol)\"",
    "@ some other line",
    "0.0000 0.0000000 -4.9941921 53.961510 -31.388611 -1.8704289 0.0000000 1.2094297 1.6674712",
    "0.0030 0.0000000 -4.9344740 53.162354 -31.474159 -1.7668755 0.0000000 1.1680595 1.6700362",
    "0.0060 0.0000000 -4.8577471 53.050510 -31.679737 -1.7318789 0.0000000 1.1528042 1.6698903",
    "&",
    "skipped"
};

//! Helper function to write an actual file to read in
void makeXvgFileToRead(const std::string &filename)
{
    std::string fileContents;
    std::for_each(std::begin(xvgLines), std::end(xvgLines), [&](const char *s) { fileContents += s; fileContents += "\n"; });
    gmx::TextWriter::writeFileFromString(filename, fileContents);
}

TEST(XvgTest, read_xvg_legendCanReadFromFile)
{
    double                   **data;
    int                        numColumns;
    char                      *subtitle;
    char                     **legend;
    gmx::test::TestFileManager fileManager;
    std::string                filename = fileManager.getTemporaryFilePath("sample.xvg");
    makeXvgFileToRead(filename);
    int                        numRows = read_xvg_legend(filename.c_str(), &data, &numColumns, &subtitle, &legend);
    EXPECT_EQ(9, numColumns);
    EXPECT_EQ(3, numRows);
    EXPECT_EQ(-4.9941921, data[2][0]);
    // TODO The next test (or something similar) should pass, so is a
    // bug. But we'll replace read_xvg_legend (only used in gmx bar)
    // with something like readXvgTable, and that'll be OK (probably
    // with matching fix to gmx bar).
    EXPECT_NE("Time (ps)", legend[0]);
    EXPECT_STREQ("dH/d\\xl\\f{} bonded-lambda = 0.6000", legend[3]);
    EXPECT_STREQ("T = 298 (K) \\xl\\f{} state 12: (fep-lambda, coul-lambda, vdw-lambda, bonded-lambda) = (0.0000, 1.0000, 0.6000, 0.6000)", subtitle);
    for (int i = 0; i < numColumns; ++i)
    {
        sfree(data[i]);
    }
    sfree(data);
    sfree(legend);
    sfree(subtitle);
}

TEST(XvgTest, CanReadFromStream)
{
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    gmx::XvgTable          table = gmx::readXvgTable(&reader);

    EXPECT_EQ("T = 298 (K) \\xl\\f{} state 12: (fep-lambda, coul-lambda, vdw-lambda, bonded-lambda) = (0.0000, 1.0000, 0.6000, 0.6000)", table.subtitle_);
    EXPECT_EQ("dH/d\\xl\\f{} and \\xD\\f{}H (kJ/mol [\\xl\\f{}]\\S-1\\N)", table.yAxisLabel_);
    EXPECT_EQ("Time (ps)", table.ordinate_.legend_);
    EXPECT_EQ(8, table.columns_.size());
    EXPECT_EQ(3, table.columns_[0].rows_.size());
    EXPECT_EQ(-4.9941921, table.columns_[1].rows_[0]);
    EXPECT_EQ("dH/d\\xl\\f{} bonded-lambda = 0.6000", table.columns_[3].legend_);
    EXPECT_EQ("dH/d\\xl\\f{} and \\xD\\f{}H", table.title_);
}

TEST(XvgTest, ThrowsIfLineHasTooFewColumns)
{
    const char            *xvgLines[] = {
        "3 4 5",
        "5 4",
        "6 4 5",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_THROW_GMX(readXvgTable(&reader), ExpectedException);
}

TEST(XvgTest, ThrowsIfLineHasTooManyColumns)
{
    const char            *xvgLines[] = {
        "3 4 5",
        "5 4 7 8",
        "6 4 5",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_THROW_GMX(readXvgTable(&reader), ExpectedException);
}

TEST(XvgTest, ThrowsIfLineIsEmpty)
{
    const char            *xvgLines[] = {
        "3e-5 4.2 -5",
        "",
        "6 4 5",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_THROW_GMX(readXvgTable(&reader), ExpectedException);
}

TEST(XvgTest, ToleratesUnusualWhitespace)
{
    const char            *xvgLines[] = {
        "  3\t4.2 -5",
        "2  3.1   1e-10\t",
        "6\t\t4E+10 5   ",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_NO_THROW_GMX(readXvgTable(&reader));
}

TEST(XvgTest, ThrowsIfDataIsNotNumerical)
{
    const char            *xvgLines[] = {
        "3 4 5 1 2",
        "2 the quick brown fox",
        "6 4 5 2 23",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_THROW_GMX(readXvgTable(&reader), ExpectedException);
}

TEST(XvgTest, ReadsPlainDataFile)
{
    const char            *xvgLines[] = {
        "3 4 5",
        "6 4 5"
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_NO_THROW_GMX(readXvgTable(&reader));
}

TEST(XvgTest, DoesNotThrowIfSetOfLegendIsIncomplete)
{
    const char            *xvgLines[] = {
        "@legend string 0 \"first\"",
        "@legend string 2 \"second\"",
        "3 4 5",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_NO_THROW_GMX(readXvgTable(&reader));
}

TEST(XvgTest, ThrowsIfLegendStringIsDuplicated)
{
    const char            *xvgLines[] = {
        "@legend string 0 \"first\"",
        "@legend string 0 \"second\"",
        "3 4 5",
    };
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    EXPECT_THROW_GMX(readXvgTable(&reader), ExpectedException);
}

} // namespace
} // namespace
