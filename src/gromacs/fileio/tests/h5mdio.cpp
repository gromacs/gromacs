/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for H5MD file I/O routines
 *
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/h5md_io.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"


// #include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"


namespace
{

class H5mdIoTest : public ::testing::Test
{
public:
    H5mdIoTest()
    {
        clear_mat(testBox_);
        clear_mat(refBox_);
        referenceFilename_ = fileManager_.getTemporaryFilePath(getFileSuffix("ref")).u8string();
        testFilename_      = fileManager_.getTemporaryFilePath(getFileSuffix("test")).u8string();
    }
    void openReferenceFile(const char mode)
    {
        referenceH5mdIo_.openFile(referenceFilename_, mode);
    }
    void closeReferenceFile()
    {
        referenceH5mdIo_.closeFile();
    }
    bool isReferenceFileOpen()
    {
        return referenceH5mdIo_.isFileOpen();
    }
    void openTestFile(const char mode)
    {
        testH5mdIo_.openFile(testFilename_, mode);
    }
    void closeTestFile()
    {
        testH5mdIo_.closeFile();
    }
    bool isTestFileOpen()
    {
        return testH5mdIo_.isFileOpen();
    }

private:
    static std::string getFileSuffix(const char* type)
    {
        return std::string(type) + "." + ftp2ext(efH5MD);
    }

    gmx::test::TestFileManager fileManager_;
    std::string                referenceFilename_;
    std::string                testFilename_;
    GmxH5mdIo                  referenceH5mdIo_;
    GmxH5mdIo                  testH5mdIo_;
    std::vector<gmx::RVec>     refX_;
    matrix                     refBox_;
    std::vector<gmx::RVec>     testX_;
    matrix                     testBox_;
};

TEST_F(H5mdIoTest, CanCreateAndCloseH5mdFile)
{
    EXPECT_THROW_GMX(openReferenceFile('r'), gmx::FileIOError);
    EXPECT_FALSE(isReferenceFileOpen());
    openReferenceFile('w');
    EXPECT_TRUE(isReferenceFileOpen());
    closeReferenceFile();
    EXPECT_FALSE(isReferenceFileOpen());
    openReferenceFile('r');
    EXPECT_TRUE(isReferenceFileOpen());
    closeReferenceFile();
    EXPECT_FALSE(isReferenceFileOpen());
}

} // namespace
