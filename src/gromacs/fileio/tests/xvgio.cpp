/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * Implements tests for xvg file operations
 *
 * \author Joe Jordan <ejjordan@kth.se>
 */
#include "gmxpre.h"

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{

class XvgioTest : public ::testing::Test
{
public:
    XvgioTest() { referenceFilename_ = fileManager_.getTemporaryFilePath("ref.xvg"); }

    const std::string& referenceFilename() const { return referenceFilename_; }

    const std::string& referenceContents() const { return referenceContents_; }

    void useStringAsXvgFile(const std::string& xvgString) { referenceContents_ = xvgString; }

    void writeXvgFile()
    {
        gmx::TextWriter::writeFileFromString(referenceFilename(), referenceContents());
    }

    static void compareValues(basic_mdspan<const double, dynamicExtents2D> ref,
                              basic_mdspan<const double, dynamicExtents2D> test)
    {
        // The xvg reading routines use a column-major layout, while we would
        // like to enforce row major behaviour everywhere else. This requires
        // this test to swap the orders between reference data and test data.
        // Hence, we compare extent(0) with extent(1) and [i][j] with [j][i].
        EXPECT_EQ(ref.extent(0), test.extent(1));
        EXPECT_EQ(ref.extent(1), test.extent(0));

        for (std::ptrdiff_t i = 0; i < ref.extent(0); i++)
        {
            for (std::ptrdiff_t j = 0; j < ref.extent(1); j++)
            {
                EXPECT_DOUBLE_EQ(ref[i][j], test[j][i]);
            }
        }
    }

private:
    gmx::test::TestFileManager fileManager_;
    std::string                referenceFilename_;
    std::string                referenceContents_;
};

TEST_F(XvgioTest, readXvgIntWorks)
{
    useStringAsXvgFile(
            "1 2 3\n"
            "4 5 6\n");
    writeXvgFile();
    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgTestData = readXvgData(referenceFilename());

    const int                                            numRows    = 2;
    const int                                            numColumns = 3;
    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgRefData(numRows, numColumns);
    std::iota(begin(xvgRefData), end(xvgRefData), 1);

    compareValues(xvgRefData.asConstView(), xvgTestData.asConstView());
}

TEST_F(XvgioTest, readXvgRealWorks)
{
    useStringAsXvgFile(
            "1.1 2.2\n"
            "3.3 4.4\n"
            "5.5 6.6\n");
    writeXvgFile();
    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgTestData = readXvgData(referenceFilename());

    const int                                            numRows    = 3;
    const int                                            numColumns = 2;
    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgRefData(numRows, numColumns);
    std::generate(begin(xvgRefData), end(xvgRefData), [n = 0.0]() mutable {
        n += 1.1;
        return n;
    });
    compareValues(xvgRefData.asConstView(), xvgTestData.asConstView());
}

TEST_F(XvgioTest, readXvgIgnoreCommentLineWorks)
{
    useStringAsXvgFile(
            "1 2 3\n"
            "#comment\n"
            "4 5 6\n");
    writeXvgFile();

    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgTestData = readXvgData(referenceFilename());

    const int                                            numRows    = 2;
    const int                                            numColumns = 3;
    MultiDimArray<std::vector<double>, dynamicExtents2D> xvgRefData(numRows, numColumns);
    std::iota(begin(xvgRefData), end(xvgRefData), 1);

    compareValues(xvgRefData.asConstView(), xvgTestData.asConstView());
}

// TODO Remove this test once all calls to read_xvg have been ported to readXvgData
TEST_F(XvgioTest, readXvgDeprecatedWorks)
{
    useStringAsXvgFile(
            "1 2 3\n"
            "4 5 6\n");
    writeXvgFile();
    std::vector<std::vector<double>> xvgData = { { 1, 4 }, { 2, 5 }, { 3, 6 } };

    double** xvgTestData = nullptr;
    int      testNumColumns;
    int      testNumRows = read_xvg(referenceFilename().c_str(), &xvgTestData, &testNumColumns);

    double** xvgRefData    = nullptr;
    int      refNumColumns = 3;
    int      refNumRows    = 2;

    EXPECT_EQ(refNumColumns, testNumColumns);
    EXPECT_EQ(refNumRows, testNumRows);

    // Set the reference data
    snew(xvgRefData, refNumColumns);
    for (int column = 0; column < refNumColumns; column++)
    {
        snew(xvgRefData[column], refNumRows);
        for (int row = 0; row < refNumRows; row++)
        {
            xvgRefData[column][row] = xvgData[column][row];
        }
    }

    // Check that the reference and test data match
    for (int column = 0; column < refNumColumns; column++)
    {
        for (int row = 0; row < refNumRows; row++)
        {
            EXPECT_EQ(xvgRefData[column][row], xvgTestData[column][row]);
        }
    }

    // Free the reference and test data memory
    for (int column = 0; column < refNumColumns; column++)
    {
        sfree(xvgRefData[column]);
        sfree(xvgTestData[column]);
    }
    sfree(xvgRefData);
    sfree(xvgTestData);
}

} // namespace test
} // namespace gmx
