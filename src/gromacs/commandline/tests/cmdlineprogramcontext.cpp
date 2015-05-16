/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Tests for gmx::CommandLineProgramContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "gromacs/commandline/cmdlineprogramcontext.h"

#include "config.h"

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>

#include "buildinfo.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/path.h"

#include "testutils/cmdlinetest.h"

using gmx::test::CommandLine;
using gmx::Path;

#if (defined GMX_NATIVE_WINDOWS || defined GMX_CYGWIN)
//! Extension for executable files on the platform.
#define EXECUTABLE_EXTENSION ".exe"
#else
//! Extension for executable files on the platform.
#define EXECUTABLE_EXTENSION ""
//! Defined if the platform supports symlinks and those can be tested.
#define TEST_SYMLINKS
#endif

namespace
{

class TestExecutableEnvironment : public gmx::ExecutableEnvironmentInterface
{
    public:
        TestExecutableEnvironment()
            : workingDirectory_(CMAKE_BINARY_DIR "/src/gromacs/commandline/tests/test-bin")
        {
        }

        virtual std::string getWorkingDirectory() const
        {
            return workingDirectory_;
        }
        virtual std::vector<std::string> getExecutablePaths() const
        {
            return path_;
        }

        std::string               workingDirectory_;
        std::vector<std::string>  path_;

        GMX_DISALLOW_COPY_AND_ASSIGN(TestExecutableEnvironment);
};

//! Shorthand for a smart pointer to TestExecutableEnvironment.
typedef boost::shared_ptr<TestExecutableEnvironment>
    TestExecutableEnvironmentPointer;

class CommandLineProgramContextTest : public ::testing::Test
{
    public:
        CommandLineProgramContextTest()
            : env_(new TestExecutableEnvironment())
        {
            expectedExecutable_ =
                Path::normalize(
                        Path::join(env_->getWorkingDirectory(),
                                   "bin/test-exe" EXECUTABLE_EXTENSION));
        }

        void testBinaryPathSearch(const char *argv0)
        {
            ASSERT_TRUE(env_.get() != NULL);
            gmx::CommandLineProgramContext  info(1, &argv0, env_);
            EXPECT_EQ(expectedExecutable_, info.fullBinaryPath());
        }
        void testBinaryPathSearch(const std::string &argv0)
        {
            testBinaryPathSearch(argv0.c_str());
        }

        std::string                      expectedExecutable_;
        TestExecutableEnvironmentPointer env_;
};

TEST_F(CommandLineProgramContextTest, FindsBinaryWithAbsolutePath)
{
    testBinaryPathSearch(Path::join(env_->getWorkingDirectory(), "bin/test-exe"));
}

TEST_F(CommandLineProgramContextTest, FindsBinaryWithRelativePath)
{
    testBinaryPathSearch("bin/test-exe");
}

TEST_F(CommandLineProgramContextTest, FindsBinaryFromPath)
{
    env_->path_.push_back(Path::join(env_->getWorkingDirectory(), "bin"));
    testBinaryPathSearch("test-exe");
}

TEST_F(CommandLineProgramContextTest, FindsBinaryFromCurrentDirectory)
{
    env_->workingDirectory_ = Path::join(env_->getWorkingDirectory(), "bin");
    env_->path_.push_back("");
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
