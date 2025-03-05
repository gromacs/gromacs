/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for grompp directives parsing
 *
 * \author Eliane Briand <eliane@br.iand.fr>
 */

#include "gmxpre.h"

#include <filesystem>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/conftest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

using gmx::test::CommandLine;
using gmx::test::TestFileManager;

enum class ExpectedResult
{
    Success,
    Death
};

class GromppDirectiveTest :
    public ::testing::TestWithParam<std::tuple<std::string, ExpectedResult, std::string>>
{
public:
    GromppDirectiveTest() = default;

protected:
    gmx::test::TestFileManager fileManager_;
    std::string                mdpContentString_ =
            "title                   = Directive edge case test \n"
            "integrator              = md \n"
            "nsteps                  = 1 \n"
            "dt                      = 0.002 \n"
            "vdwtype                 = cutoff \n"
            "coulombtype             = cutoff \n"
            "tcoupl                  = no \n"
            "pcoupl                  = no \n"
            "pbc                     = xyz \n"
            "gen_vel                 = yes \n";
};

TEST_F(GromppDirectiveTest, edgeCaseAtomTypeNames)
{
    CommandLine cmdline;
    cmdline.addOption("grompp");

    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath("directives.mdp").string();
    gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContentString_);
    cmdline.addOption("-f", mdpInputFileName);


    cmdline.addOption("-c", TestFileManager::getInputFilePath("directives.gro").string());
    cmdline.addOption("-p", TestFileManager::getInputFilePath("directives.top").string());

    std::string outTprFilename = fileManager_.getTemporaryFilePath("directives.tpr").string();
    cmdline.addOption("-o", outTprFilename);

    ASSERT_EQ(0, gmx_grompp(cmdline.argc(), cmdline.argv()));
    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state;
        read_tpx_state(outTprFilename, &ir_after, &state, &top_after);

        int indexInMoltype = top_after.molblock[0].type;

        // Check atomic numbers (or lack thereof coded as -1)
        ASSERT_EQ(top_after.moltype[indexInMoltype].atoms.nr, 4);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[0].atomnumber, -1);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[1].atomnumber, 6);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[2].atomnumber, 7);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[3].atomnumber, -1);
    }
}

TEST_F(GromppDirectiveTest, NoteOnDihedralNotSumToZero)
{
    CommandLine cmdline;
    cmdline.addOption("grompp");

    std::string mdpString = mdpContentString_;
    mdpString += "define = -DDIHEDRAL_SUM_NOT_ZERO";

    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath("directives.mdp").string();
    gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpString);
    cmdline.addOption("-f", mdpInputFileName);


    cmdline.addOption("-c", TestFileManager::getInputFilePath("directives.gro").string());
    cmdline.addOption("-p", TestFileManager::getInputFilePath("directives.top").string());

    std::string outTprFilename = fileManager_.getTemporaryFilePath("directives.tpr").string();
    cmdline.addOption("-o", outTprFilename);

    // We cannot directly check printing of a note, but we at least check that it terminates
    // successfully.
    EXPECT_EQ(gmx_grompp(cmdline.argc(), cmdline.argv()), 0);
}

TEST_F(GromppDirectiveTest, WarnOnDihedralSumDifferentForFreeEnergy)
{
    CommandLine cmdline;
    cmdline.addOption("grompp");

    std::string mdpString = mdpContentString_;
    mdpString +=
            "define = -DDIHEDRAL_SUM_DIFFERENT_STATEA_STATEB\n"
            "free-energy = yes\n"
            "init-lambda = 0.5";

    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath("directives.mdp").string();
    gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpString);
    cmdline.addOption("-f", mdpInputFileName);


    cmdline.addOption("-c", TestFileManager::getInputFilePath("directives.gro").string());
    cmdline.addOption("-p", TestFileManager::getInputFilePath("directives.top").string());

    std::string outTprFilename = fileManager_.getTemporaryFilePath("directives.tpr").string();
    cmdline.addOption("-o", outTprFilename);

    GMX_EXPECT_DEATH_IF_SUPPORTED(gmx_grompp(cmdline.argc(), cmdline.argv()),
                                  "undesired offset in dHdl values");
}

TEST_P(GromppDirectiveTest, AcceptValidAndErrorOnInvalidCMAP)
{
    auto testParam = GetParam();

    CommandLine cmdline;
    cmdline.addOption("grompp");

    std::string mdpString = mdpContentString_;
    mdpString += std::get<0>(testParam);

    const std::string mdpInputFileName =
            fileManager_.getTemporaryFilePath("directives-cmap.mdp").string();
    gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpString);
    cmdline.addOption("-f", mdpInputFileName);


    cmdline.addOption("-c", TestFileManager::getInputFilePath("directives-cmap.gro").string());
    cmdline.addOption("-p", TestFileManager::getInputFilePath("directives-cmap.top").string());

    std::string outTprFilename = fileManager_.getTemporaryFilePath("directives-cmap.tpr").string();
    cmdline.addOption("-o", outTprFilename);

    switch (std::get<1>(testParam))
    {
        case ExpectedResult::Success:
            EXPECT_EQ(gmx_grompp(cmdline.argc(), cmdline.argv()), 0);
            break;
        case ExpectedResult::Death:
            GMX_EXPECT_DEATH_IF_SUPPORTED(gmx_grompp(cmdline.argc(), cmdline.argv()),
                                          std::get<2>(testParam));
            break;
        default: FAIL();
    }
}

std::vector<std::tuple<std::string, ExpectedResult, std::string>> cmapValidInputOutput = {
    { "", ExpectedResult::Death, "Unknown cmap torsion between atoms 1 2 3 4 5" },
    { "define = -DNOT_A_CMAPTYPE",
      ExpectedResult::Death,
      "Unknown atomtype 1 found at position 5 in cmap type" },
    { "define = -DMATCHING_CMAPTYPE", ExpectedResult::Success, "" },
    { "define = -DMATCHING_CMAPTYPE_DOUBLESPACED", ExpectedResult::Success, "" },
    { "define = -DMATCHING_CMAPTYPE_PADDED", ExpectedResult::Success, "" },
    { "define = -DMATCHING_CMAPTYPE_TABBED", ExpectedResult::Success, "" },
    { "define = -DUNKNOWN_ATOMTYPE_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Unknown bond_atomtype for Z in cmap atomtypes X Y X X Z" },
    { "define = -DTOO_MANY_ATOMTYPES_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Invalid function type for cmap type: must be a number, found Y" },
    { "define = -DTOO_FEW_ATOMTYPES_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Invalid function type for cmap type: must be 1" },
    { "define = -DINVALID_FUNCTYPE_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Invalid function type for cmap type: must be 1" },
    { "define = -DRECTANGULAR_GRID_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Not the same grid spacing in x and y for cmap grid: x=2, y=3" },
    { "define = -DUNREAL_GRID_SIZE_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Invalid cmap type grid spacings in x and y dimensions: must be numbers,\n  found Tarydium" },
    { "define = -DTOO_FEW_GRID_PARAMETERS_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Error in reading cmap parameter for atomtypes X Y X X Y: found 3,\n  expected 4" },
    { "define = -DTOO_MANY_GRID_PARAMETERS_IN_CMAPTYPE",
      ExpectedResult::Death,
      "One or more unread cmap parameters exist for atomtypes X Y X X Y" },
    { "define = -DUNREAL_GRID_PARAMETER_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Invalid cmap parameters for atomtypes X Y X X Y: must be real numbers,\n  found Tarydium" },
    { "define = -DSOME_RESIDUE_NAMES_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Incorrect format for cmap atomtypes X Y X X Y, residuetypes are required\n  for all 5 "
      "atomtypes or none" },
    { "define = -DMATCHING_RESIDUE_STARS_IN_CMAPTYPE", ExpectedResult::Success, "" },
    { "define = -DMATCHING_RESIDUE_NAMES_IN_CMAPTYPE", ExpectedResult::Success, "" },
    { "define = -DNONMATCHING_RESIDUE_NAMES_IN_CMAPTYPE",
      ExpectedResult::Death,
      "Unknown cmap torsion between atoms 1 2 3 4 5" },
    { "define = -DNOT_A_CMAP_TORSION", ExpectedResult::Death, "Too few parameters on line" },
    { "define = -DINVALID_FUNCTYPE_IN_CMAP_TORSION",
      ExpectedResult::Death,
      "Invalid function type for cmap torsion: must be 1" },
    { "define = -DUSER_SPECIFIED_CMAPTYPE", ExpectedResult::Success, "" },
    { "define = -DUSER_SPECIFIED_CMAPTYPE_OUT_OF_BOUNDS",
      ExpectedResult::Death,
      "Unable to assign a cmap type to torsion 1 2 3 4 and 5" }
};

INSTANTIATE_TEST_SUITE_P(CMAPDefinesAndErrors, GromppDirectiveTest, testing::ValuesIn(cmapValidInputOutput));

} // namespace
} // namespace test
} // namespace gmx
