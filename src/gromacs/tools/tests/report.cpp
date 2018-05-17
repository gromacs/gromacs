/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "gromacs/tools/report.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/exceptions.h"
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

void generateAndReadTprInput(std::string filename, gmx_mtop_t *mtop, t_inputrec *ir)
{
    TestFileManager   fileManager;
// generate temporary tpr file from test system
    const std::string mdpInputFileName = fileManager.getTemporaryFilePath(filename + ".mdp");
    TextWriter::writeFileFromString(mdpInputFileName, "");
    std::string       tprName = fileManager.getTemporaryFilePath(filename + ".tpr");
    {
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(filename));
        caller.addOption("-c", TestFileManager::getInputFilePath(filename + ".pdb"));
        caller.addOption("-o", tprName);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }
// read tpr into variables needed for output
    {
        t_state     state;
        std::string tprName = fileManager.getTemporaryFilePath(filename + ".tpr");
        read_tpx_state(tprName.c_str(), ir, &state, mtop);
    }
}


TEST(ReportTest, WritesCorrectInformation)
{
    gmx_mtop_t top;
    t_inputrec ir;
    EXPECT_NO_THROW(generateAndReadTprInput("lysozyme", &top, &ir));
}


} // namespace

} // namespace test

} // namespace gmx
