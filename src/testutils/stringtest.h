/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::test::StringTestBase.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_STRINGTEST_H
#define GMX_TESTUTILS_STRINGTEST_H

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class FileOutputRedirectorInterface;

namespace test
{

class TestFileManager;
class TestReferenceChecker;

/*! \libinternal \brief
 * Test fixture for tests that check string formatting.
 *
 * For development, tests that use this fixture as their base can be run with a
 * '-stdout' command-line option to print out the tested strings to stdout.
 * If this flag is not given, they check the strings using the XML reference
 * framework (see TestReferenceData).
 *
 * \inlibraryapi
 * \ingroup module_testutils
 */
class StringTestBase : public ::testing::Test
{
    public:
        StringTestBase();
        ~StringTestBase();

        /*! \brief
         * Creates a redirector that directs all output to temporary files.
         *
         * \param[in] fileManager  File manager to use for temporary files.
         *
         * Can only be called once in a test.
         *
         * \see checkRedirectedOutputFiles()
         */
        FileOutputRedirectorInterface &
        initOutputRedirector(TestFileManager *fileManager);

        /*! \brief
         * Returns the root checker for this test's reference data.
         *
         * Can be used to perform custom checks against reference data (e.g.,
         * if the test needs to check some other values than plain strings.
         */
        TestReferenceChecker &checker();

        /*! \brief
         * Checks a string.
         *
         * \param[in] text  String to check.
         * \param[in] id    Unique (within a single test) id for the string.
         */
        void checkText(const std::string &text, const char *id);
        /*! \brief
         * Checks contents of a file as a single string.
         *
         * \param[in] filename  Name of the file to check.
         * \param[in] id        Unique (within a single test) id for the string.
         *
         * Provided for convenience.  Reads the contents of \p filename into a
         * single string and calls checkText().
         */
        void checkFileContents(const std::string &filename, const char *id);
        /*! \brief
         * Checks contents of all files redirected with initOutputRedirector().
         *
         * Uses the same logic as checkFileContents() to check each file
         * (including `stdout`) that has been created using the redirector
         * returned by initOutputRedirector().
         *
         * initOutputRedirector() must have been called.
         * This method should not be called if the redirector will still be
         * used for further output in the test.  Behavior is not designed for
         * checking in the middle of the test, although that could potentially
         * be changed if necessary.
         */
        void checkRedirectedOutputFiles();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif
