/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for gmx convert-tpr that cause memory leaks
 *
 * Some legacy features/functions invoked in convert-tpr have memory leaks,
 * which are picked up in AddressSanitizer runs.
 * To avoid blanket suppression, the tests causing the leaks are isolated here.
 *
 * In particular, the following leaf-level functions are known to leak in
 * these tests: atomcat (mtop_util.cpp), do_atoms (tpxio.cpp), rd_groups
 * (topology/index.cpp).
 *
 * These functions are called and/or cause leaks when slicing the topology
 * according to an index group.
 *
 * \author Eliane Briand <eliane@br.iand.fr>
 */
#include "gmxpre.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/topology/topology.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stdiohelper.h"
#include "testutils/tprfilegenerator.h"

#include "convert-tpr-fixture.h"

namespace gmx
{
namespace test
{
namespace
{

TEST_F(ConvertTprTest, selectIndexTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    ASSERT_EQ(top.natoms, 156);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const std::string           command[]      = {
        "convert-tpr", "-s", tprFileHandle.tprName(), "-o", outTprFilename.string()
    };
    CommandLine cmdline(command);

    StdioTestHelper stdioHelper(&fileManager);
    stdioHelper.redirectStringToStdin("Protein-H\n");

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename, &ir_after, &state_after, &top_after);
        EXPECT_EQ(top_after.natoms, 75);
    }
}

TEST_F(ConvertTprNoVelocityTest, selectIndexTestWithoutVelocity)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    ASSERT_EQ(top.natoms, 156);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const std::string           command[]      = {
        "convert-tpr", "-s", tprFileHandle.tprName(), "-o", outTprFilename.string()
    };
    CommandLine cmdline(command);

    StdioTestHelper stdioHelper(&fileManager);
    stdioHelper.redirectStringToStdin("Protein-H\n");

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename, &ir_after, &state_after, &top_after);
        EXPECT_EQ(top_after.natoms, 75);
    }
}

} // namespace
} // namespace test
} // namespace gmx
