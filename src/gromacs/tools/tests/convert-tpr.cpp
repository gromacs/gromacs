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
 * Tests for gmx convert-tpr.
 *
 * \author Eliane Briand <eliane@br.iand.fr>
 */
#include "gmxpre.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/topology/topology.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
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


TEST_F(ConvertTprTest, ExtendRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName().c_str(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t extendByPs     = 100;
    std::string   extendByString = std::to_string(extendByPs);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr").u8string();
    const char* const command[]      = {
        "convert-tpr",          "-s",      tprFileHandle.tprName().c_str(), "-o",
        outTprFilename.c_str(), "-extend", extendByString.c_str()
    };
    CommandLine cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename.c_str(), &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(extendByPs / ir.delta_t));
    }


    // Extending again (tests nsteps not zero initially

    std::string anotherOutTprFilename = fileManager.getTemporaryFilePath("extended_again.tpr").u8string();

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
        t_state    state_after;
        read_tpx_state(anotherOutTprFilename.c_str(), &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(2 * extendByPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, UntilRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName().c_str(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t untilPs       = 100;
    std::string   untilPsString = std::to_string(untilPs);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr").u8string();
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
        t_state    state;
        read_tpx_state(outTprFilename.c_str(), &ir_after, &state, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(untilPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, nstepRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName().c_str(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t nsteps    = 102;
    std::string   nstepsStr = std::to_string(nsteps);

    TestFileManager   fileManager;
    std::string       outTprFilename = fileManager.getTemporaryFilePath("extended.tpr").u8string();
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
        t_state    state_after;
        read_tpx_state(outTprFilename.c_str(), &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + nsteps);
    }
}

TEST_F(ConvertTprTest, generateVelocitiesTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName().c_str(), &ir, &state, &top);

    TestFileManager fileManager;
    std::string outTprFilename  = fileManager.getTemporaryFilePath("new_velocities.tpr").u8string();
    const char* const command[] = { "convert-tpr",
                                    "-s",
                                    tprFileHandle.tprName().c_str(),
                                    "-o",
                                    outTprFilename.c_str(),
                                    "-generate_velocities",
                                    "-velocity_temp",
                                    "300",
                                    "-velocity_seed",
                                    "12345" };
    CommandLine       cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename.c_str(), &ir_after, &state_after, &top_after);

        gmx::test::TestReferenceData    data;
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        std::vector<real>               result;
        // Check that the X coordinates did NOT change
        for (int i = 0; i < state.x.size(); i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                EXPECT_EQ(state.x[i][j], state_after.x[i][j]);
                result.push_back(state_after.v[i][j]);
            }
        }
        checker.checkSequence(result.begin(), result.end(), "ConvertTprTestgenerateVelocitiesTestV");
    }
}

} // namespace
} // namespace test
} // namespace gmx
