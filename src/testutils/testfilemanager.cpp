/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements gmx::test::TestFileManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_testutils
 */
#include "testfilemanager.h"

#include <cstdio>

#include <algorithm>
#include <set>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * TestFileManager::Impl
 */

/*! \internal \brief
 * Private implementation class for TestFileManager.
 *
 * \ingroup module_testutils
 */
class TestFileManager::Impl
{
    public:
        //! Global test input data path set with setDataInputDirectory().
        static const char *s_inputDirectory;

        //! Container type for names of temporary files.
        typedef std::set<std::string> FileNameList;

        /*! \brief
         * Try to remove all temporary files.
         *
         * Does not throw; errors (e.g., missing files) are silently ignored.
         */
        void removeFiles();

        //! List of unique paths returned by getTemporaryFilePath().
        FileNameList            files_;
};

const char *TestFileManager::Impl::s_inputDirectory = NULL;

void TestFileManager::Impl::removeFiles()
{
    FileNameList::const_iterator i;
    for (i = files_.begin(); i != files_.end(); ++i)
    {
        std::remove(i->c_str());
    }
    files_.clear();
}

/********************************************************************
 * TestFileManager
 */

TestFileManager::TestFileManager()
    : impl_(new Impl)
{
}

TestFileManager::~TestFileManager()
{
    impl_->removeFiles();
}

std::string TestFileManager::getTemporaryFilePath(const char *suffix)
{
    // TODO: Add the path of the test binary
    std::string filename = getTestSpecificFileName(suffix);
    impl_->files_.insert(filename);
    return filename;
}

std::string TestFileManager::getTestSpecificFileName(const char *suffix)
{
    const ::testing::TestInfo *test_info =
            ::testing::UnitTest::GetInstance()->current_test_info();
    std::string                filename = std::string(test_info->test_case_name())
        + "_" + test_info->name();
    if (suffix[0] != '.')
    {
        filename.append("_");
    }
    filename.append(suffix);
    std::replace(filename.begin(), filename.end(), '/', '_');
    return filename;
}

std::string TestFileManager::getInputFilePath(const char *filename)
{
    return Path::join(getInputDataDirectory(), filename);
}

const char *TestFileManager::getInputDataDirectory()
{
    GMX_RELEASE_ASSERT(Impl::s_inputDirectory != NULL, "Test data path not set");
    return Impl::s_inputDirectory;
}

void TestFileManager::setInputDataDirectory(const char *path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(Directory::exists(path),
                       "Test data directory does not exist");
    Impl::s_inputDirectory = path;
}

} // namespace test
} // namespace gmx
