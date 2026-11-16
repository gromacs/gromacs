/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Tests for compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_compression.h"

#include <hdf5.h>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/binary_information.h"
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

//! Compression scheme and the HDF5 filters that are expected to be set for it
struct TestCompressionData
{
    H5mdCompression           compression;
    std::vector<H5Z_filter_t> expectedHdf5Filters;
};

//! \brief Helper function for GTest to print test parameters
void PrintTo(const TestCompressionData& info, std::ostream* os)
{
    *os << enumValueToString(info.compression);
}

class H5mdCompressionFilterTest : public H5mdTestBase, public ::testing::WithParamInterface<TestCompressionData>
{
};

//! \brief Helper function for GTest to construct test names.
std::string nameOfCompressionTest(const ::testing::TestParamInfo<TestCompressionData>& info)
{
    std::string testName(enumValueToString(info.param.compression));

    // Note that the returned names must be unique and may use only
    // alphanumeric ASCII characters. It's not supposed to contain
    // underscores (see the GoogleTest FAQ
    // why-should-test-suite-names-and-test-names-not-contain-underscore),
    // but doing so works for now, is likely to remain so, and makes
    // such test names much more readable.
    testName = replaceAll(testName, "-", "_");
    testName = replaceAll(testName, ".", "_");
    testName = replaceAll(testName, " ", "_");

    return testName;
}

TEST_P(H5mdCompressionFilterTest, FiltersAreRegistered)
{
    if (!sz3CompressionIsAvailable() && GetParam().compression == H5mdCompression::LossySz3BioMD)
    {
        GTEST_SKIP() << "GROMACS is not built with SZ3 support";
    }

    const auto [propertyList, propertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));
    setCompressionFilters(propertyList, GetParam().compression, 5e-4);

    ASSERT_EQ(H5Pget_nfilters(propertyList), GetParam().expectedHdf5Filters.size());
    for (const H5Z_filter_t filter : GetParam().expectedHdf5Filters)
    {
        // H5Pget_filter_by_id2 returns: >=0 if the filter (second argument) is set, else <0
        EXPECT_GE(H5Pget_filter_by_id2(propertyList, filter, nullptr, nullptr, nullptr, 0, nullptr, nullptr), 0)
                << "Expected filter was not registered";
    }
}

TEST_P(H5mdCompressionFilterTest, GetSz3AbsoluteErrorWorks)
{
    if (!sz3CompressionIsAvailable())
    {
        GTEST_SKIP() << "GROMACS is not built with SZ3 support";
    }

    constexpr double absoluteError = 5e-4;

    H5mdDataSetBase<RVec> dataSet = H5mdDataSetBuilder<RVec>(fileid(), "testDataSet")
                                            .withCompression(GetParam().compression, absoluteError)
                                            .withDimension({ 0, 0 })
                                            .build();

    if (GetParam().compression == H5mdCompression::LossySz3BioMD)
    {
        EXPECT_THAT(getSz3AbsoluteError(dataSet.id()), ::testing::Optional(absoluteError));
    }
    else
    {
        EXPECT_FALSE(getSz3AbsoluteError(dataSet.id()).has_value());
    }
}

TEST_P(H5mdCompressionFilterTest, IsCompressedWithSz3Works)
{
    if (!sz3CompressionIsAvailable())
    {
        GTEST_SKIP() << "GROMACS is not built with SZ3 support";
    }

    H5mdDataSetBase<BasicVector<float>> dataSet =
            H5mdDataSetBuilder<BasicVector<float>>(this->fileid(), "testDataSet")
                    .withDimension({ 1, 1 })
                    .withCompression(GetParam().compression, 5e-4)
                    .build();

    if (GetParam().compression == H5mdCompression::LossySz3BioMD)
    {
        EXPECT_TRUE(isCompressedWithSz3(dataSet.id()));
    }
    else
    {
        EXPECT_FALSE(isCompressedWithSz3(dataSet.id()));
    }
}

static const TestCompressionData sc_testCompressionData[] = {
    { H5mdCompression::Uncompressed, {} },
    { H5mdCompression::LosslessShuffle, { H5Z_FILTER_DEFLATE, H5Z_FILTER_SHUFFLE } },
    { H5mdCompression::LosslessNoShuffle, { H5Z_FILTER_DEFLATE } },
    { H5mdCompression::LossySz3BioMD, { c_sz3FilterId } },
};

INSTANTIATE_TEST_SUITE_P(Works,
                         H5mdCompressionFilterTest,
                         ::testing::ValuesIn(sc_testCompressionData),
                         nameOfCompressionTest);

} // namespace
} // namespace test
} // namespace gmx
