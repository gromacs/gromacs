/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * Implements gmx::test::TestFileManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/testfilemanager.h"

#include <cstdio>

#include <algorithm>
#include <set>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"

#include "testutils/testoptions.h"

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
    //! Global test input data path set with setInputDataDirectory().
    static std::string s_inputDirectory;

    //! Global path to simulation input database set with setTestSimulationDataBaseDirectory().
    static std::string s_simulationDatabaseDirectory;

    //! Global temporary output directory for tests, set with setGlobalOutputTempDirectory().
    static const char* s_globalOutputTempDirectory;

    //! Container type for names of temporary files.
    typedef std::set<std::string> FileNameList;

    /*! \brief Constructor
     *
     * \param path Value for the outputTempDirectory, typically
     * set by default from s_globalOutputTempDirectory */
    explicit Impl(const char* path) : outputTempDirectory_(path)
    {
        GMX_RELEASE_ASSERT(Directory::exists(outputTempDirectory_),
                           "Directory for tests' temporary files does not exist");
    }

    /*! \brief
     * Try to remove all temporary files.
     *
     * Does not throw; errors (e.g., missing files) are silently ignored.
     */
    void removeFiles();

    //! List of unique paths returned by getTemporaryFilePath().
    FileNameList files_;

    /*! \brief Temporary output directory local to the current
     * test, set by a test with setOutputTempDirectory() if the
     * global default is inappropriate. */
    std::string outputTempDirectory_;
};

std::string TestFileManager::Impl::s_inputDirectory;
std::string TestFileManager::Impl::s_simulationDatabaseDirectory;
const char* TestFileManager::Impl::s_globalOutputTempDirectory = nullptr;
/** Controls whether TestFileManager should delete temporary files
    after the test finishes. */
static bool g_bDeleteFilesAfterTest = true;

//! \cond
GMX_TEST_OPTIONS(TestFileManagerOptions, options)
{
    options->addOption(
            BooleanOption("delete-temporary-files")
                    .store(&g_bDeleteFilesAfterTest)
                    .description(
                            "At the end of each test case, delete temporary and output files"));
}
//! \endcond

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

TestFileManager::TestFileManager() : impl_(new Impl(Impl::s_globalOutputTempDirectory)) {}

TestFileManager::~TestFileManager()
{
    if (g_bDeleteFilesAfterTest)
    {
        impl_->removeFiles();
    }
}

std::string TestFileManager::getTemporaryFilePath(const char* suffix)
{
    /* Configure a temporary directory from CMake, so that temporary
     * output from a test goes to a location relevant to that
     * test. Currently, files whose names are returned by this method
     * get cleaned up (by default) at the end of all tests.
     */
    std::string filename = Path::join(getOutputTempDirectory(), getTestSpecificFileName(suffix));
    impl_->files_.insert(filename);
    return filename;
}

std::string TestFileManager::getTemporaryFilePath(const std::string& suffix)
{
    return getTemporaryFilePath(suffix.c_str());
}

void TestFileManager::manageGeneratedOutputFile(const char* filename)
{
    manageGeneratedOutputFile(std::string(filename));
}
void TestFileManager::manageGeneratedOutputFile(std::string&& filename)
{
    impl_->files_.insert(std::move(filename));
}

std::string TestFileManager::getTestSpecificFileNameRoot()
{
    const ::testing::TestInfo* test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string                filenameRoot;
    if (test_info)
    {
        filenameRoot = std::string(test_info->test_suite_name()) + "_" + test_info->name();
    }
    else
    {
        const ::testing::TestSuite* test_suite_info =
                ::testing::UnitTest::GetInstance()->current_test_suite();
        filenameRoot = std::string(test_suite_info->name());
    }
    std::replace(filenameRoot.begin(), filenameRoot.end(), '/', '_');
    return filenameRoot;
}

std::string TestFileManager::getTestSpecificFileName(const char* suffix)
{
    std::string filename = getTestSpecificFileNameRoot();
    if (suffix[0] != '.')
    {
        filename.append("_");
    }
    filename.append(suffix);
    return filename;
}

std::string TestFileManager::getInputFilePath(const char* filename)
{
    // Check if file is present in local directory.
    if (File::exists(Path::join(getInputDataDirectory(), filename), File::returnFalseOnError))
    {
        return Path::join(getInputDataDirectory(), filename);
    }
    else if (File::exists(Path::join(getTestSimulationDatabaseDirectory(), filename), File::returnFalseOnError))
    {
        // Assume file is in global directory for simulation input files.
        return Path::join(getTestSimulationDatabaseDirectory(), filename);
    }
    // Assume file is present locally without full name (e.g. extension).
    return Path::join(getInputDataDirectory(), filename);
}

std::string TestFileManager::getInputFilePath(const std::string& filename)
{
    return getInputFilePath(filename.c_str());
}

const char* TestFileManager::getInputDataDirectory()
{
    GMX_RELEASE_ASSERT(!Impl::s_inputDirectory.empty(), "Path for test input files is not set");
    return Impl::s_inputDirectory.c_str();
}

const char* TestFileManager::getGlobalOutputTempDirectory()
{
    GMX_RELEASE_ASSERT(Impl::s_globalOutputTempDirectory != nullptr,
                       "Global path for temporary output files from tests is not set");
    return Impl::s_globalOutputTempDirectory;
}

const char* TestFileManager::getOutputTempDirectory() const
{
    return impl_->outputTempDirectory_.c_str();
}

const char* TestFileManager::getTestSimulationDatabaseDirectory()
{
    GMX_RELEASE_ASSERT(!Impl::s_simulationDatabaseDirectory.empty(),
                       "Path for simulation input database directory is not set");
    return Impl::s_simulationDatabaseDirectory.c_str();
}

void TestFileManager::setInputDataDirectory(const std::string& path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(Directory::exists(path), "Test data directory does not exist");
    Impl::s_inputDirectory = path;
}

void TestFileManager::setTestSimulationDatabaseDirectory(const std::string& path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(Directory::exists(path), "Simulation database directory does not exist");
    Impl::s_simulationDatabaseDirectory = path;
}

void TestFileManager::setGlobalOutputTempDirectory(const char* path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(Directory::exists(path),
                       "Directory for tests' temporary files does not exist");
    Impl::s_globalOutputTempDirectory = path;
}

void TestFileManager::setOutputTempDirectory(const std::string& path)
{
    // There could be a need to protect this with a mutex, since it is
    // intended to be used in test fixtures, not just during setup.
    GMX_RELEASE_ASSERT(Directory::exists(path),
                       "Directory for tests' temporary files does not exist");
    impl_->outputTempDirectory_ = path;
}

} // namespace test
} // namespace gmx
