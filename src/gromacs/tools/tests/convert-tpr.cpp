/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * Tests for gmx convert-tpr.
 *
 * \author Eliane Briand <eliane@br.iand.fr>
 */
#include "gmxpre.h"

#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/topology/topology.h"
#include "gromacs/fileio/tpxio.h"

#include "testutils/cmdlinetest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{
namespace test
{
namespace
{

class ConvertTprTest : public ::testing::Test
{
protected:
    ConvertTprTest() : tprFileHandle("lysozyme") {}

    //! Storage for opened file handles.
    TprAndFileManager tprFileHandle;
};


void readTprInput(const char* filename, gmx_mtop_t* mtop, t_inputrec* ir)
{
    // read tpr into variables needed for output
    {
        t_state state;
        read_tpx_state(filename, ir, &state, mtop);
    }
}

TEST_F(ConvertTprTest, ExtendRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    readTprInput(tprFileHandle.tprName().c_str(), &top, &ir);

    const int64_t originalNStep = ir.nsteps;

    const int64_t extendByPs     = 100;
    std::string   extendByString = std::to_string(extendByPs);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const char* const command[]      = {
        "convert-tpr",          "-s",      tprFileHandle.tprName().c_str(), "-o",
        outTprFilename.c_str(), "-extend", extendByString.c_str()
    };
    CommandLine cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        readTprInput(outTprFilename.c_str(), &top_after, &ir_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(extendByPs / ir.delta_t));
    }


    // Extending again (tests nsteps not zero initially

    std::string anotherOutTprFilename = fileManager.getTemporaryFilePath("extended_again.tpr");

    const char* const secondCommand[] = { "convert-tpr",
                                          "-s",
                                          outTprFilename.c_str(),
                                          "-o",
                                          anotherOutTprFilename.c_str(),
                                          "-extend",
                                          extendByString.c_str() };
    CommandLine       secondCmdline(secondCommand);
    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &secondCmdline);


    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        readTprInput(anotherOutTprFilename.c_str(), &top_after, &ir_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(2 * extendByPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, UntilRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    readTprInput(tprFileHandle.tprName().c_str(), &top, &ir);

    const int64_t originalNStep = ir.nsteps;

    const int64_t untilPs       = 100;
    std::string   untilPsString = std::to_string(untilPs);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const char* const command[]      = { "convert-tpr",
                                    "-s",
                                    tprFileHandle.tprName().c_str(),
                                    "-o",
                                    outTprFilename.c_str(),
                                    "-until",
                                    untilPsString.data() };
    CommandLine       cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        readTprInput(outTprFilename.c_str(), &top_after, &ir_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(untilPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, nstepRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    readTprInput(tprFileHandle.tprName().c_str(), &top, &ir);

    const int64_t originalNStep = ir.nsteps;

    const int64_t nsteps    = 102;
    std::string   nstepsStr = std::to_string(nsteps);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const char* const command[]      = { "convert-tpr",
                                    "-s",
                                    tprFileHandle.tprName().c_str(),
                                    "-o",
                                    outTprFilename.c_str(),
                                    "-nsteps",
                                    nstepsStr.data() };
    CommandLine       cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        readTprInput(outTprFilename.c_str(), &top_after, &ir_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + nsteps);
    }
}

} // namespace
} // namespace test
} // namespace gmx
