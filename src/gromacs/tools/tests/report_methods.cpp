/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * Tests for functionality of the "report" tool to write system information.
 *
 * \author
 */
#include "gmxpre.h"

#include "gromacs/tools/report_methods.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{

namespace test
{

namespace
{

/*! \brief
 * Generates a tpr for the test.
 *
 * Generates the tpr from a sample pdb file using grompp,and returns the path
 * to the file as std::string for reading it in later.
 *
 * \param[in] fileManager Filemanager to keep track of the input file.
 * \param[in] filename Basename of the input and output files.
 */
std::string generateTprInput(TestFileManager *fileManager, const std::string &filename)
{
// generate temporary tpr file from test system
    const std::string mdpInputFileName = fileManager->getTemporaryFilePath(filename + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, "");
    std::string       tprName = fileManager->getTemporaryFilePath(filename + ".tpr");
    {
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(filename + ".top"));
        caller.addOption("-c", TestFileManager::getInputFilePath(filename + ".pdb"));
        caller.addOption("-o", tprName);
        EXPECT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }
    return tprName;
}

/*! \brief
 * Generates and reads a tpr for the test.
 *
 * Generates a tpr from a sample pdb file using grompp, and reads it in to
 * have access to the system information for print out.
 *
 * \param[in] filename Basename of the input and output files.
 * \param[in] mtop Pointer to topology datastructure to populate.
 * \param[in] ir Pointer to inputrec to populate.
 */
void generateAndReadTprInput(const std::string &filename, gmx_mtop_t *mtop, t_inputrec *ir)
{
    TestFileManager   fileManager;
    std::string       tprName = generateTprInput(&fileManager, filename);
// read tpr into variables needed for output
    {
        t_state     state;
        read_tpx_state(tprName.c_str(), ir, &state, mtop);
    }
}

TEST(ReportMethodsTest, WritesCorrectHeadersFormated)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         test(&stream);
    writeHeader(&test, "Test", "Test", true);

    std::string referenceString = "\\Test{Test}\n";

    EXPECT_EQ(stream.toString(), referenceString);
}
TEST(ReportMethodsTest, WritesCorrectHeadersUnformatted)
{
    gmx::StringOutputStream stream;
    gmx::TextWriter         test(&stream);
    writeHeader(&test, "Test", "Test", false);

    std::string referenceString = "Test: Test\n";

    EXPECT_EQ(stream.toString(), referenceString);
}

TEST(ReportMethodsTest, WritesCorrectInformation)
{
    gmx_mtop_t top;
    t_inputrec ir;
    EXPECT_NO_THROW(generateAndReadTprInput("lysozyme", &top, &ir));

    // Test both formatted and unformatted output in the same test
    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeSystemInformation(&test, top, true);

        std::string referenceString = "\\subsection{Simulation system}\n"
            "A system of 1 molecules (156 atoms) was simulated.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }

    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeSystemInformation(&test, top, false);

        std::string referenceString = "subsection: Simulation system\n"
            "A system of 1 molecules (156 atoms) was simulated.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }

    // Test for correct parameter information as well
    {
        gmx::StringOutputStream stream;
        gmx::TextWriter         test(&stream);

        writeParameterInformation(&test, ir, true);

        std::string referenceString = "\\subsection{Simulation settings}\n"
            "A total of 0 ns were simulated with a time step of 1 fs.\n"
            "Neighbor searching was performed every 10 steps.\n"
            "The Cut-off algorithm was used for electrostatic interactions.\n"
            "with a cut-off of 1 nm.\nA single cut-off of 1.1 nm was used "
            "for Van der Waals interactions.\n";

        EXPECT_EQ(stream.toString(), referenceString);
    }
}

TEST(ReportMethodsTest, ToolEndToEndTest)
{
    TestFileManager   fileManager;
    std::string       tprName   = generateTprInput(&fileManager, "lysozyme");
    const char *const command[] = {
        "report-methods", "-s", tprName.c_str()
    };
    CommandLine       cmdline(command);
    EXPECT_EQ(0, gmx::test::CommandLineTestHelper::runModuleFactory(
                      &gmx::ReportMethodsInfo::create, &cmdline));

}

} // namespace

} // namespace test

} // namespace gmx
