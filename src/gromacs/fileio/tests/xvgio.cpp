/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Implements tests for xvg file operations
 *
 * \author Joe Jordan <ejjordan@kth.se>
 */
#include "gmxpre.h"

#include <cstddef>

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{
namespace
{

void compareValues(basic_mdspan<const double, dynamicExtents2D> ref,
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

} // namespace

/*! \brief
 * Convienience type for testing read_xvg_time.
 *
 * Fields are: haveStartTime, haveEndTime
 */
using XvgrTimeReadingParams = std::tuple<bool, bool>;

class XvgioTest : public ::testing::Test, public ::testing::WithParamInterface<XvgrTimeReadingParams>
{
public:
    XvgioTest() { referenceFilename_ = fileManager_.getTemporaryFilePath("ref.xvg").string(); }

    const std::string& referenceFilename() const { return referenceFilename_; }

    const std::string& referenceContents() const { return referenceContents_; }

    void useStringAsXvgFile(const std::string& xvgString) { referenceContents_ = xvgString; }

    void writeXvgFile() const
    {
        gmx::TextWriter::writeFileFromString(referenceFilename(), referenceContents());
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

namespace
{
//! Helper for checking file contents against reference file
void checkMatrix(TestReferenceChecker* checker, basic_mdspan<const double, dynamicExtents2D> input)
{
    TestReferenceChecker compound(checker->checkCompound("XvgMatrix", nullptr));

    std::vector<real> values;
    for (int i = 0; i < input.extent(0); ++i)
    {
        std::vector<real> row(input.extent(1));
        for (int j = 0; j < input.extent(1); ++j)
        {
            row[j] = input[i][j];
        }
        values.insert(values.end(), row.begin(), row.end());
    }
    compound.checkSequence(values.begin(), values.end(), "ColumnValues");
}

} // namespace

TEST_P(XvgioTest, readXvgTimeSeriesWorks)
{
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    useStringAsXvgFile(
            "0.2 1 2 3\n"
            "0.4 1 3 2\n"
            "0.8 3 2 1\n"
            "1.0 3 1 2\n"
            "1.4 1 2 3\n");
    writeXvgFile();

    auto params       = GetParam();
    bool useStartTime = std::get<0>(params);
    bool useEndTime   = std::get<1>(params);

    auto timeSeriesData = readXvgTimeSeries(referenceFilename(),
                                            useStartTime ? std::make_optional(0.3) : std::nullopt,
                                            useEndTime ? std::make_optional(1.2) : std::nullopt);

    checkMatrix(&checker, timeSeriesData);
}

INSTANTIATE_TEST_SUITE_P(XvgReadTimeSeries,
                         XvgioTest,
                         ::testing::Combine(

                                 ::testing::Values(true, false),
                                 ::testing::Values(true, false)));


} // namespace test
} // namespace gmx
