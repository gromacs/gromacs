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
#include "testutils/testfilemanager.h"
#include "testutils/textblockmatchers.h"

namespace
{

using gmx::test::CommandLine;
using gmx::test::TestFileManager;

class GromppDirectiveTest : public ::testing::Test
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

    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath("directives.mdp");
    gmx::TextWriter::writeFileFromString(mdpInputFileName, mdpContentString_);
    cmdline.addOption("-f", mdpInputFileName);


    cmdline.addOption("-c", TestFileManager::getInputFilePath("directives.gro"));
    cmdline.addOption("-p", TestFileManager::getInputFilePath("directives.top"));

    std::string outTprFilename = fileManager_.getTemporaryFilePath("directives.tpr");
    cmdline.addOption("-o", outTprFilename);

    ASSERT_EQ(0, gmx_grompp(cmdline.argc(), cmdline.argv()));
    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state;
        read_tpx_state(outTprFilename.c_str(), &ir_after, &state, &top_after);

        int indexInMoltype = top_after.molblock[0].type;

        // Check atomic numbers (or lack thereof coded as -1)
        ASSERT_EQ(top_after.moltype[indexInMoltype].atoms.nr, 4);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[0].atomnumber, -1);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[1].atomnumber, 6);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[2].atomnumber, 7);
        EXPECT_EQ(top_after.moltype[indexInMoltype].atoms.atom[3].atomnumber, -1);
    }
}

} // namespace
