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

   \ingroup module_verbosity_management
 */
#include "gmxpre.h"

#include "gromacs/utility/verbosity_management.h"

#include <fstream>  // ::std::ostream (<<)
#include <iomanip>  // ::std::setw, ::std::setfill
#include <ios>      // ::std::left, ::std::rightusing std::cout;
#include <iostream> // ::std::cout
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

//! stream handles that can be used by outside code as shown below
//! these ostreams use cstdio internally and thus do not interfere
//! with output via cstdio's (f)printf
using gmx::verbosity_management::debug;
using gmx::verbosity_management::error;
using gmx::verbosity_management::blab0;
using gmx::verbosity_management::blab1;
using gmx::verbosity_management::blab2;
using gmx::verbosity_management::blab3;


class VerbosityManagementTest : public ::testing::Test
{
    public:
        VerbosityManagementTest()
            : stdout_buf_(stdout), stderr_buf_(stderr),
              save_cout_buf_pt_(std::cout.rdbuf()),
              save_cerr_buf_pt_(std::cerr.rdbuf())
        {
        }
        void SetUp ()
        {
            //! redirect std::cout/cerr to use fprintf for cstdio compatibility
            //! such that the google testing framework can provide information on why a test failed
            std::cout.rdbuf(&stdout_buf_);
            std::cerr.rdbuf(&stderr_buf_);
        }
        void TearDown ()
        {
            //! redirect std::cout/cerr to their original target buffers from before the test
            std::cout.rdbuf(save_cout_buf_pt_);
            std::cerr.rdbuf(save_cerr_buf_pt_);
        }
        gmx::test::TestFileManager      fileManager_;
    private:
        //! stream buffer that uses fprintf for output
        gmx::CstdioStreamBuffer<255> stdout_buf_;
        //! stream buffer that uses fprintf for output
        gmx::CstdioStreamBuffer<255> stderr_buf_;
        //! save a pointer to the original output buffer of std::cout
        std::streambuf              *save_cout_buf_pt_ = std::cout.rdbuf();
        //! save a pointer to the original output buffer of std::cerr
        std::streambuf              *save_cerr_buf_pt_ = std::cout.rdbuf();
};

/*********************************************************************/

TEST_F(VerbosityManagementTest, Globals)
{
    // set verbosity
    gmx::verbosity_management::setDefaults();
    if (dbg_bool)
    {
        // default verbosity
        gmx::verbosity_management::setDebugOut(true);
        gmx::verbosity_management::setErrorOut(true);
    }

    debug << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << "Running basic test 1: globals (verbosity management)"                            << endl;
    debug << "-------------------------------------------------------------------------------" << endl;
    debug << endl;

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test triggering of debug output:"                        << endl;

    gmx::verbosity_management::setDebugOut(true);
    EXPECT_TRUE(debug) << "Switching on debug output failed";

    gmx::verbosity_management::setDebugOut();
    gmx::verbosity_management::setDebugOut(false);
    EXPECT_FALSE(debug) << "Switching off debug output failed";

    gmx::verbosity_management::setDebugOut(false);
    gmx::verbosity_management::setDebugOut();
    EXPECT_TRUE(debug) << "Switching on debug output again after switching it off failed";
    gmx::verbosity_management::setDebugOut(dbg_bool);

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test triggering of error output:"                        << endl;

    gmx::verbosity_management::setErrorOut(false);
    gmx::verbosity_management::setErrorOut(true);
    EXPECT_TRUE(error) << "Switching on error output failed";

    gmx::verbosity_management::setErrorOut();
    gmx::verbosity_management::setErrorOut(false);
    EXPECT_FALSE(error) << "Switching off error output failed";
    gmx::verbosity_management::setErrorOut(true);

    debug << ">>-----------------------------------------------------" << endl;
    debug << "Test triggering of error output:"                        << endl;

    std::vector<gmx::StreamAdapter*> osv(4);

    osv[0] = &gmx::verbosity_management::blab0;
    osv[1] = &gmx::verbosity_management::blab1;
    osv[2] = &gmx::verbosity_management::blab2;
    osv[3] = &gmx::verbosity_management::blab3;
    for (unsigned int i = 0; i < 4; ++i)
    {
        debug << "testing verbosity/blab level i = " << i << " ..." << endl;
        gmx::verbosity_management::setBlabLevel(i);
        // all ostreams blabn with n <= level i should be active ...
        for (unsigned int j = 0; (j <= i && j <= 3); ++j)
        {
            if (dbg_bool)
            {
                (*osv[j]) << "blab" << j << " output" << endl;
            }
            EXPECT_TRUE((*osv[j])) << "verbosity level j = " << j << "   <=   i = " << i << " should have been activated";
        }
        // ... the remaining ostreams blabn with n > level i should not be active ...
        for (unsigned int j = i+1; j <= 3; ++j)
        {
            if (dbg_bool)
            {
                (*osv[j]) << "blab" << j << " output" << endl;
            }
            EXPECT_FALSE((*osv[j])) << "verbosity level j = " << j << "   >   i = " << i << " should not have been activated";
        }
    }

    // restore default verbosity
    gmx::verbosity_management::setDefaults();
    if (dbg_bool)
    {
        // default verbosity
        gmx::verbosity_management::setDebugOut(true);
        gmx::verbosity_management::setErrorOut(true);
    }
    else
    {
        // default verbosity
        gmx::verbosity_management::setDebugOut(false);
        gmx::verbosity_management::setErrorOut(true);
        gmx::verbosity_management::setBlabLevel(0);
    }

/*  // this shows how to redirect output to files and back again to sdout/stderr
    gmx::verbosity_management::setRedirectStdOutToFile(true);
    gmx::verbosity_management::setRedirectErrOutToFile(true);
    error << "This output, normally destined to stderr, should be redirected to a file" << endl;
    debug << "This output, normally destined to stdout, should be redirected to a file" << endl;

    gmx::verbosity_management::setRedirectStdOutToFile(false);
    gmx::verbosity_management::setRedirectErrOutToFile(false);
    error << "This output, normally destined to stderr, should again be redirected to stderr" << endl;
    debug << "This output, normally destined to stdout, should again be redirected to stdout" << endl;
 */
    debug << "-----------------------------------------------------<<" << endl;
}


} // namespace
