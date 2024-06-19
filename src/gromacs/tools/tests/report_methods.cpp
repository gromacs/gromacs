/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Tests for functionality of the "report" tool to write system information.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/tools/report_methods.h"

#include <filesystem>
#include <functional>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{

namespace test
{

class ReportMethodsTest : public ::testing::Test
{
protected:
    // TODO this is changed in newer googletest versions
    //! Prepare shared resources.
    static void SetUpTestSuite() { s_tprFileHandle = new TprAndFileManager("lysozyme"); }
    //! Clean up shared resources.
    static void TearDownTestSuite()
    {
        delete s_tprFileHandle;
        s_tprFileHandle = nullptr;
    }
    //! Storage for opened file handles.
    static TprAndFileManager* s_tprFileHandle;
};

TprAndFileManager* ReportMethodsTest::s_tprFileHandle = nullptr;

/*! \brief
 * Reads a tpr for the test.
 *
 * Reads a tpr to have access to the system information for print out.
 *
 * \param[in] tprHandle Handle to the tpr to red in.
 * \param[in] mtop Pointer to topology datastructure to populate.
 * \param[in] ir Pointer to inputrec to populate.
 */
static void readTprInput(const TprAndFileManager* tprHandle, gmx_mtop_t* mtop, t_inputrec* ir)
{
    // read tpr into variables needed for output
    {
        t_state state;
        read_tpx_state(tprHandle->tprName(), ir, &state, mtop);
    }
}

TEST_F(ReportMethodsTest, WritesCorrectHeadersFormated)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         test(&stream);
    writeHeader(&test, "Test", "Test", true);

    std::string referenceString = "\\Test{Test}\n";

    EXPECT_EQ(stream.toString(), referenceString);
}
TEST_F(ReportMethodsTest, WritesCorrectHeadersUnformatted)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         test(&stream);
    writeHeader(&test, "Test", "Test", false);

    std::string referenceString = "Test: Test\n";

    EXPECT_EQ(stream.toString(), referenceString);
}

TEST_F(ReportMethodsTest, WritesCorrectInformation)
{
    gmx_mtop_t top;
    t_inputrec ir;
    EXPECT_NO_THROW(readTprInput(s_tprFileHandle, &top, &ir));

    // Test both formatted and unformatted output in the same test
    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeSystemInformation(&test, top, true);

        std::string referenceString =
                "\\subsection{Simulation system}\n"
                "A system of 1 molecules (156 atoms) was simulated.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }

    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeSystemInformation(&test, top, false);

        std::string referenceString =
                "subsection: Simulation system\n"
                "A system of 1 molecules (156 atoms) was simulated.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }

    // Test for correct parameter information as well
    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeParameterInformation(&test, ir, true);

        std::string referenceString =
                "\\subsection{Simulation settings}\n"
                "A total of 0 ns were simulated with a time step of 1 fs.\n"
                "Neighbor searching was performed every 10 steps.\n"
                "The Cut-off algorithm was used for electrostatic interactions.\n"
                "with a cut-off of 1 nm.\nA single cut-off of 1.1 nm was used "
                "for Van der Waals interactions.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }
}

TEST_F(ReportMethodsTest, ToolEndToEndTest)
{
    const char* const command[] = { "report-methods", "-s", s_tprFileHandle->tprName().c_str() };
    CommandLine       cmdline(command);
    EXPECT_EQ(0, gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ReportMethodsInfo::create, &cmdline));
}

} // namespace test

} // namespace gmx
