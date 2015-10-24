/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
   \brief
   Tests for the the verbosity_management module

   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/cstdio_ostream.h"

#include <iomanip>  // ::std::setw, ::std::setfill
#include <ios>      // ::std::left, ::std::rightusing std::cout;
#include <iostream> // ::std::cout
#include <ostream>  // ::std::ostream (<<)
#include <stdexcept>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"

namespace
{

//! flag that enables (true) or disables (false) debug output
constexpr bool dbg_bool = false;

using std::endl;


class CstdioOstreamTest : public ::testing::Test
{
    public:
        CstdioOstreamTest()
            : stdout_buf_((size_t)255, stdout), stderr_buf_((size_t)255, stderr), cnull_buf_(),
              save_cout_buf_pt_(std::cout.rdbuf()),
              save_cerr_buf_pt_(std::cerr.rdbuf()),
              debug(&cnull_buf_), error(&cnull_buf_)
        {
        }
        void SetUp()
        {
            // redirect std::cout/cerr to use fprintf for cstdio compatibility
            // such that the google testing framework can provide information on why a test failed
            std::cout.rdbuf(&stdout_buf_);
            std::cerr.rdbuf(&stderr_buf_);
            // trigger debug output
            if (dbg_bool)
            {
                activateDebugOut();
                activateErrorOut();
            }
        }
        void TearDown()
        {
            // redirect std::cout/cerr to their original target buffers from before the test
            std::cout.rdbuf(save_cout_buf_pt_);
            std::cerr.rdbuf(save_cerr_buf_pt_);
        }

        void activateDebugOut()
        {
            debug.rdbuf(&stdout_buf_);
        }
        void activateErrorOut()
        {
            error.rdbuf(&stderr_buf_);
        }
        void deactivateDebugOut()
        {
            debug.rdbuf(&cnull_buf_);
        }
        void deactivateErrorOut()
        {
            error.rdbuf(&cnull_buf_);
        }

        gmx::test::TestFileManager   fileManager_;
    private:
        // stream buffer that uses fprintf for output
        gmx::test::CstdioStreamBuffer      stdout_buf_;
        // stream buffer that uses fprintf for output
        gmx::test::CstdioStreamBuffer      stderr_buf_;
        // no-effect stream buffer that, silent sink for output
        gmx::test::CstdioStreamBuffer      cnull_buf_;
        // save a pointer to the original output buffer of std::cout
        std::streambuf                    *save_cout_buf_pt_;
        // save a pointer to the original output buffer of std::cerr
        std::streambuf                    *save_cerr_buf_pt_;
    public:
        // debugging output stream, silent by default
        gmx::test::TestableOstream         debug;
        // error output stream, silent by default
        gmx::test::TestableOstream         error;
};

/*********************************************************************/

TEST_F(CstdioOstreamTest, SetupDebugOutput)
{

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running basic test: cstdio ostream functionality"                                << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test triggering of error output:"                        << endl;

    error << "This error message should appear." << endl;
    activateErrorOut();
    EXPECT_TRUE(error)  << "activation   of error output failed.";
    deactivateErrorOut();
    error << "This error message should not appear." << endl;
    EXPECT_FALSE(error) << "deactivation of error output failed.";

    // restore the original setting
    if (dbg_bool)
    {
        activateErrorOut();
    }
    else
    {
        deactivateErrorOut();
    }

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test triggering of debug output:"                        << endl;

    debug << "This debug message should appear." << endl;
    activateDebugOut();
    EXPECT_TRUE(debug)   << "activation   of debug output failed.";
    deactivateDebugOut();
    debug << "This debug message should not appear." << endl;
    EXPECT_FALSE(debug)  << "deactivation of debug output failed.";

    // restore the original setting
    if (dbg_bool)
    {
        activateDebugOut();
    }
    else
    {
        deactivateDebugOut();
    }

    debug << "-----------------------------------------------------<<" << endl;
}


} // namespace
