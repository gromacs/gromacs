/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Test routines that parse mdp fields from grompp input and writes
 * mdp back out.
 *
 * In particular these will provide test coverage as we refactor to
 * use a new Options-based key-value-style mdp implementation to
 * support a more modular mdrun.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/readir.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/warninp.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

class GetIrTest : public ::testing::Test
{
    public:
        GetIrTest() : fileManager_(), data_(), checker_(data_.rootChecker()),
                      ir_(), mdModules_(), opts_(),
                      wi_(init_warning(FALSE, 0)), wiGuard_(wi_)

        {
            snew(opts_.include, STRLEN);
            snew(opts_.define, STRLEN);
        }
        ~GetIrTest()
        {
            done_inputrec_strings();
            sfree(opts_.include);
            sfree(opts_.define);
        }
        /*! \brief Test mdp reading and writing
         *
         * \todo Modernize read_inp and write_inp to use streams,
         * which will make these tests run faster, because they don't
         * use disk files. */
        void runTest(const std::string &inputMdpFileContents)
        {
            auto inputMdpFilename  = fileManager_.getTemporaryFilePath("input.mdp");
            auto outputMdpFilename = fileManager_.getTemporaryFilePath("output.mdp");

            TextWriter::writeFileFromString(inputMdpFilename, inputMdpFileContents);

            get_ir(inputMdpFilename.c_str(), outputMdpFilename.c_str(),
                   &mdModules_, &ir_, &opts_, WriteMdpHeader::no, wi_);
            bool failure = warning_errors_exist(wi_);
            checker_.checkBoolean(failure, "Error parsing mdp file");
            warning_reset(wi_);

            auto outputMdpContents = TextReader::readFileToString(outputMdpFilename);
            checker_.checkString(outputMdpContents, "OutputMdpFile");
        }

        TestFileManager                    fileManager_;
        TestReferenceData                  data_;
        TestReferenceChecker               checker_;
        t_inputrec                         ir_;
        MDModules                          mdModules_;
        t_gromppopts                       opts_;
        warninp_t                          wi_;
        unique_cptr<warninp, free_warning> wiGuard_;
};

TEST_F(GetIrTest, HandlesDifferentKindsOfMdpLines)
{
    const char *inputMdpFile[] = {
        "; File to run my simulation",
        "title = simulation",
        ";",
        "xtc_grps = System ; was Protein",
        "include = -I/home/me/stuff",
        "",
        "tau-t = 0.1 0.3",
        "tinit = 0.3",
        "init_step = 0",
        "nstcomm = 100",
        "integrator = steep"
    };
    runTest(joinStrings(inputMdpFile, "\n"));
}

// This case is used often by SimulationRunner::useEmptyMdpFile (see
// comments there for explanation). When we remove the group scheme,
// that usage will have to disappear, and so can this test.
TEST_F(GetIrTest, HandlesOnlyCutoffScheme)
{
    const char *inputMdpFile = "cutoff-scheme = Group\n";
    runTest(inputMdpFile);
}

// TODO Stop accepting any of these
TEST_F(GetIrTest, UserErrorsSilentlyTolerated)
{
    const char *inputMdpFile[] = {
        "title simulation",
        "xtc_grps = ",
        "= -I/home/me/stuff",
        "="
    };
    runTest(joinStrings(inputMdpFile, "\n"));
}

TEST_F(GetIrTest, EmptyInputWorks)
{
    const char *inputMdpFile = "";
    runTest(inputMdpFile);
}

// These tests observe how the electric-field keys behave, since they
// are currently the only ones using the new Options-style handling.
TEST_F(GetIrTest, ProducesOutputFromElectricField)
{
    const char *inputMdpFile = "electric-field-x = 1.2 0 0 0";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, ProducesOutputFromElectricFieldPulsed)
{
    const char *inputMdpFile = "electric-field-y = 3.7 2.0 6.5 1.0";
    runTest(inputMdpFile);
}

TEST_F(GetIrTest, ProducesOutputFromElectricFieldOscillating)
{
    const char *inputMdpFile = "electric-field-z = 3.7 7.5 0 0";
    runTest(inputMdpFile);
}

} // namespace
} // namespace
