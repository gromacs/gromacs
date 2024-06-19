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
#include <filesystem>
#include <memory>
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
    static std::filesystem::path s_inputDirectory;

    //! Global path to simulation input database set with setTestSimulationDataBaseDirectory().
    static std::filesystem::path s_simulationDatabaseDirectory;

    //! Global temporary output directory for tests, set with setGlobalOutputTempDirectory().
    static std::filesystem::path s_globalOutputTempDirectory;

    //! Container type for names of temporary files.
    typedef std::set<std::filesystem::path> FileNameList;

    /*! \brief Constructor
     *
     * \param path Value for the outputTempDirectory, typically
     * set by default from s_globalOutputTempDirectory */
    explicit Impl(const std::filesystem::path& path) : outputTempDirectory_(path)
    {
        GMX_RELEASE_ASSERT(std::filesystem::exists(outputTempDirectory_),
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
    std::filesystem::path outputTempDirectory_;
};

std::filesystem::path TestFileManager::Impl::s_inputDirectory;
std::filesystem::path TestFileManager::Impl::s_simulationDatabaseDirectory;
std::filesystem::path TestFileManager::Impl::s_globalOutputTempDirectory;
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
    for (const auto& file : files_)
    {
        std::filesystem::remove(file);
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

std::filesystem::path TestFileManager::getTemporaryFilePath(const std::filesystem::path& suffix)
{
    /* Configure a temporary directory from CMake, so that temporary
     * output from a test goes to a location relevant to that
     * test. Currently, files whose names are returned by this method
     * get cleaned up (by default) at the end of all tests.
     */
    auto filename = getOutputTempDirectory() / getTestSpecificFileName(suffix);
    impl_->files_.insert(filename);
    return filename;
}

void TestFileManager::manageGeneratedOutputFile(const std::filesystem::path& filename)
{
    impl_->files_.insert(filename);
}

std::filesystem::path TestFileManager::getTestSpecificFileNameRoot()
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

std::filesystem::path TestFileManager::getTestSpecificFileName(const std::filesystem::path& suffix)
{
    std::string filename = getTestSpecificFileNameRoot().string();
    if (suffix.string()[0] != '.')
    {
        filename.append("_");
    }
    filename.append(suffix.string());
    return filename;
}

std::filesystem::path TestFileManager::getInputFilePath(const std::filesystem::path& filename)
{
    const std::string& name = filename.string();
    // Check if file is present in local directory.
    std::filesystem::path file;
    return File::exists(file = getInputDataDirectory() / name, File::returnFalseOnError) ? file
           : File::exists(file = getTestSimulationDatabaseDirectory() / name, File::returnFalseOnError)
                   ? file
                   : getInputDataDirectory() / name;
}

std::filesystem::path TestFileManager::getInputDataDirectory()
{
    GMX_RELEASE_ASSERT(!Impl::s_inputDirectory.empty(), "Path for test input files is not set");
    return Impl::s_inputDirectory;
}

std::filesystem::path TestFileManager::getGlobalOutputTempDirectory()
{
    GMX_RELEASE_ASSERT(!Impl::s_globalOutputTempDirectory.empty(),
                       "Global path for temporary output files from tests is not set");
    return Impl::s_globalOutputTempDirectory;
}

std::filesystem::path TestFileManager::getOutputTempDirectory() const
{
    return impl_->outputTempDirectory_;
}

std::filesystem::path TestFileManager::getTestSimulationDatabaseDirectory()
{
    GMX_RELEASE_ASSERT(!Impl::s_simulationDatabaseDirectory.empty(),
                       "Path for simulation input database directory is not set");
    return Impl::s_simulationDatabaseDirectory;
}

void TestFileManager::setInputDataDirectory(const std::filesystem::path& path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(std::filesystem::exists(path) && std::filesystem::is_directory(path),
                       "Test data directory does not exist");
    Impl::s_inputDirectory = path;
}

void TestFileManager::setTestSimulationDatabaseDirectory(const std::filesystem::path& path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(std::filesystem::exists(path) && std::filesystem::is_directory(path),
                       "Simulation database directory does not exist");
    Impl::s_simulationDatabaseDirectory = path;
}

void TestFileManager::setGlobalOutputTempDirectory(const std::filesystem::path& path)
{
    // There is no need to protect this by a mutex, as this is called in early
    // initialization of the tests.
    GMX_RELEASE_ASSERT(std::filesystem::exists(path) && std::filesystem::is_directory(path),
                       "Directory for tests' temporary files does not exist");
    Impl::s_globalOutputTempDirectory = path;
}

void TestFileManager::setOutputTempDirectory(const std::filesystem::path& path)
{
    // There could be a need to protect this with a mutex, since it is
    // intended to be used in test fixtures, not just during setup.
    GMX_RELEASE_ASSERT(std::filesystem::exists(path) && std::filesystem::is_directory(path),
                       "Directory for tests' temporary files does not exist");
    impl_->outputTempDirectory_ = path;
}

} // namespace test
} // namespace gmx
