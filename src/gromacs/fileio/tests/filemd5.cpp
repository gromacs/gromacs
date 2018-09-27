/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "gromacs/fileio/gmxfio.h"

#include <string>

#include <gtest/gtest.h>

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
        FileMD5Test()
        {
            // Pick an arbitrary file already present during the run
            // of the test.
            auto filename = TestFileManager::getInputFilePath("CMakeLists.txt");
            file_ = gmx_fio_open(filename.c_str(), "r+");
        }
        ~FileMD5Test()
        {
            gmx_fio_close(file_);
        }
        gmx::test::TestReferenceData     data_;
        gmx::test::TestReferenceChecker  checker_ = data_.rootChecker();
        t_fileio *file_;
};

TEST_F(FileMD5Test, CanComputeMD5)
{
    unsigned char           digest[16];
    ArrayRef<unsigned char> digestView(digest);
    // Chosen to be less than the full file length
    gmx_off_t               expectedLength     = 64;
    gmx_off_t               lengthActuallyRead = gmx_fio_get_file_md5(file_, expectedLength, digest);

    checker_.checkInteger(lengthActuallyRead, "Length Read");
    checker_.checkSequence(digestView.begin(), digestView.end(), "Digest");
}

TEST_F(FileMD5Test, WillNotComputeMD5PastEndOfFile)
{
    unsigned char           digest[16];
    ArrayRef<unsigned char> digestView(digest);
    // Chosen to be greater than the full file length
    gmx_off_t               expectedLength     = 10000;
    gmx_off_t               lengthActuallyRead = gmx_fio_get_file_md5(file_, expectedLength, digest);

    checker_.checkInteger(lengthActuallyRead, "Length Read");
    checker_.checkSequence(digestView.begin(), digestView.end(), "Digest");
}

} // namespace
} // namespace test
} // namespace gmx
