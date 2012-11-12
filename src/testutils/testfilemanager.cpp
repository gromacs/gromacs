/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
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
        FileNameList files_;
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
    std::string filename = std::string(test_info->test_case_name()) +
        "_" + test_info->name();
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
