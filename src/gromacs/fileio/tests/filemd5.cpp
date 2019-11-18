/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Tests for routines for computing MD5 sums on files.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include <cstdio>

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

class FileMD5Test : public ::testing::Test
{
public:
    void prepareFile(int lengthInBytes)
    {
        // Fill some memory with some arbitrary bits.
        std::vector<char> data(lengthInBytes);
        std::iota(data.begin(), data.end(), 1);
        // Binary mode ensures it works the same on all OS
        FILE* fp = fopen(filename_.c_str(), "wb");
        fwrite(data.data(), sizeof(char), data.size(), fp);
        fclose(fp);
    }
    ~FileMD5Test() override
    {
        if (file_)
        {
            gmx_fio_close(file_);
        }
    }
    TestFileManager fileManager_;
    // Make sure the file extension is one that gmx_fio_open will
    // recognize to open as binary.
    std::string filename_ = fileManager_.getTemporaryFilePath("data.edr");
    t_fileio*   file_     = nullptr;
};

TEST_F(FileMD5Test, CanComputeMD5)
{
    prepareFile(1000);
    file_ = gmx_fio_open(filename_.c_str(), "r+");

    std::array<unsigned char, 16> digest = { 0 };
    // Chosen to be less than the full file length
    gmx_off_t offset             = 64;
    gmx_off_t expectedLength     = 64;
    gmx_off_t lengthActuallyRead = gmx_fio_get_file_md5(file_, offset, &digest);

    EXPECT_EQ(expectedLength, lengthActuallyRead);
    // Did we compute an actual reproducible checksum?
    auto total = std::accumulate(digest.begin(), digest.end(), 0);
    EXPECT_EQ(2111, total);
}

TEST_F(FileMD5Test, ReturnsErrorIfFileModeIsWrong)
{
    prepareFile(1000);
    file_ = gmx_fio_open(filename_.c_str(), "r");

    std::array<unsigned char, 16> digest;
    gmx_off_t                     offset             = 100;
    gmx_off_t                     lengthActuallyRead = gmx_fio_get_file_md5(file_, offset, &digest);
    EXPECT_EQ(-1, lengthActuallyRead);
}

} // namespace
} // namespace test
} // namespace gmx
