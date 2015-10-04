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

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/* TODO The contents of the string from which we make the .xvg file
   are duplicated below, but we'll kill read_xvg_legend and its test
   shortly */

//! Helper function to write a file to read in
void makeXvgFileToRead(const std::string &filename)
{
    gmx::TextWriter::writeFileFromString
        (filename,
        "# This file was created Sun Oct  4 19:42:39 2015\n"
        "#\n"
        "@    title \"dH/d\\xl\\f{} and \\xD\\f{}H\"\n"
        "@    xaxis  label \"Time (ps)\"\n"
        "@    yaxis  label \"dH/d\\xl\\f{} and \\xD\\f{}H (kJ/mol [\\xl\\f{}]\\S-1\\N)\"\n"
        "@TYPE xy\n"
        "@ subtitle \"T = 298 (K) \\xl\\f{} state 12: (fep-lambda, coul-lambda, vdw-lambda, bonded-lambda) = (0.0000, 1.0000, 0.6000, 0.6000)\"\n"
        "@ view 0.15, 0.15, 0.75, 0.85\n"
        "@ legend on\n"
        "@ s0 legend \"dH/d\\xl\\f{} fep-lambda = 0.0000\"\n"
        "@ s1 legend \"dH/d\\xl\\f{} coul-lambda = 1.0000\"\n"
        "@ s2 legend \"dH/d\\xl\\f{} vdw-lambda = 0.6000\"\n"
        "@ s3 legend \"dH/d\\xl\\f{} bonded-lambda = 0.6000\"\n"
        "@ s4 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.5000, 0.5000)\"\n"
        "@ s5 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.6000, 0.6000)\"\n"
        "@ s6 legend \"\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.6500, 0.6500)\"\n"
        "@ s7 legend \"pV (kJ/mol)\"\n"
        "@ some other line\n"
        "0.0000 0.0000000 -4.9941921 53.961510 -31.388611 -1.8704289 0.0000000 1.2094297 1.6674712\n"
        "0.0030 0.0000000 -4.9344740 53.162354 -31.474159 -1.7668755 0.0000000 1.1680595 1.6700362\n"
        "0.0060 0.0000000 -4.8577471 53.050510 -31.679737 -1.7318789 0.0000000 1.1528042 1.6698903\n"
        "&\n"
        "skipped");
}

TEST(XvgTest, CanReadFromFile)
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
    // TODO A test like this should pass, so is a bug. But we'll
    // replace read_xvg_legend (only used in gmx bar) with
    // readXvgTable, and that'll be OK (probably with matching fix to
    // gmx bar).
    EXPECT_NE("Time (ps)", legend[0]);
    // TODO this is the legend for .xvg data set s3, but internally
    // this should be is column with index 4, because the abscissa is
    // a column and we don't always read it. Fixed in the new version.
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

//! Array of lines to fake a file
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

TEST(XvgTest, CanReadFromStream)
{
    gmx::StringInputStream stream(xvgLines);
    gmx::TextReader        reader(&stream);
    gmx::XvgTable          table = gmx::readXvgTable(&reader);

    EXPECT_EQ(9, table.columns_.size());
    EXPECT_EQ(3, table.columns_[0].rows_.size());
    EXPECT_EQ(-4.9941921, table.columns_[2].rows_[0]);
    EXPECT_EQ("Time (ps)", table.xAxisLabel_);
    EXPECT_EQ("Time (ps)", table.columns_[0].legend_);
    // This is the legend for s3, which is correctly in the column
    // with index 4
    EXPECT_EQ("dH/d\\xl\\f{} bonded-lambda = 0.6000", table.columns_[4].legend_);
    EXPECT_EQ("T = 298 (K) \\xl\\f{} state 12: (fep-lambda, coul-lambda, vdw-lambda, bonded-lambda) = (0.0000, 1.0000, 0.6000, 0.6000)", table.subtitle_);
}

/* TODO There are many ways to break an .xvg file such that
   readXvgTable will throw exceptions. Add tests for these behaviours
   when removing calls to the deprecated read_xvg*, hopefully having
   found out more about what the expected behaviour of readXvgTable
   should be. */

} // namespace
} // namespace
