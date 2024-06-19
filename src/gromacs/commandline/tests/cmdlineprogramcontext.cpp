/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests for gmx::CommandLineProgramContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/cmdlineprogramcontext.h"

#include "config.h"

#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/path.h"

#include "testutils/cmdlinetest.h"

#include "buildinfo.h"

#if GMX_NATIVE_WINDOWS || GMX_CYGWIN
//! Extension for executable files on the platform.
#    define EXECUTABLE_EXTENSION ".exe"
#else
//! Extension for executable files on the platform.
#    define EXECUTABLE_EXTENSION ""
//! Defined if the platform supports symlinks and those can be tested.
#    define TEST_SYMLINKS
#endif

namespace
{

class TestExecutableEnvironment : public gmx::IExecutableEnvironment
{
public:
    TestExecutableEnvironment() :
        workingDirectory_(CMAKE_BINARY_DIR "/src/gromacs/commandline/tests/test-bin")
    {
    }

    std::filesystem::path getWorkingDirectory() const override { return workingDirectory_; }
    std::vector<std::filesystem::path> getExecutablePaths() const override { return path_; }

    std::filesystem::path              workingDirectory_;
    std::vector<std::filesystem::path> path_;

    GMX_DISALLOW_COPY_AND_ASSIGN(TestExecutableEnvironment);
};

//! Shorthand for a smart pointer to TestExecutableEnvironment.
typedef std::unique_ptr<TestExecutableEnvironment> TestExecutableEnvironmentPointer;

class CommandLineProgramContextTest : public ::testing::Test
{
public:
    CommandLineProgramContextTest() : env_(new TestExecutableEnvironment())
    {
        expectedExecutable_ = std::filesystem::path(env_->getWorkingDirectory())
                                      .append("bin/test-exe" EXECUTABLE_EXTENSION)
                                      .make_preferred()
                                      .string();
    }

    void testBinaryPathSearch(const char* argv0)
    {
        ASSERT_TRUE(env_.get() != nullptr);
        gmx::CommandLineProgramContext info(1, &argv0, std::move(env_));
        EXPECT_EQ(expectedExecutable_, info.fullBinaryPath());
    }
    void testBinaryPathSearch(const std::string& argv0) { testBinaryPathSearch(argv0.c_str()); }

    std::string                      expectedExecutable_;
    TestExecutableEnvironmentPointer env_;
};

TEST_F(CommandLineProgramContextTest, FindsBinaryWithAbsolutePath)
{
    testBinaryPathSearch(std::filesystem::path(env_->getWorkingDirectory()).append("bin/test-exe").string());
}

TEST_F(CommandLineProgramContextTest, FindsBinaryWithRelativePath)
{
    testBinaryPathSearch("bin/test-exe");
}

TEST_F(CommandLineProgramContextTest, FindsBinaryFromPath)
{
    env_->path_.push_back(std::filesystem::path(env_->getWorkingDirectory()).append("bin"));
    testBinaryPathSearch("test-exe");
}

TEST_F(CommandLineProgramContextTest, FindsBinaryFromCurrentDirectory)
{
    env_->workingDirectory_ = std::filesystem::path(env_->getWorkingDirectory()).append("bin");
    env_->path_.emplace_back("");
    testBinaryPathSearch("test-exe");
}

#ifdef TEST_SYMLINKS
TEST_F(CommandLineProgramContextTest, FindsBinaryFromAbsoluteSymLink)
{
    testBinaryPathSearch("bin/test-abs-link");
}

TEST_F(CommandLineProgramContextTest, FindsBinaryFromRelativeSymLink)
{
    testBinaryPathSearch("bin/test-rel-link");
}
#endif

} // namespace
