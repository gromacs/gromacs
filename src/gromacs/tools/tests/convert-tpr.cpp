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

#include <cstdint>

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/simulationdatabase.h"
#include "testutils/stdiohelper.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

#include "convert-tpr-fixture.h"

namespace gmx
{
namespace test
{
namespace
{

TEST_F(ConvertTprTest, ExtendRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t extendByPs     = 100;
    std::string   extendByString = std::to_string(extendByPs);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const std::string           command[]      = { "convert-tpr",           "-s",
                                    tprFileHandle.tprName(), "-o",
                                    outTprFilename.string(), "-extend",
                                    extendByString };
    CommandLine                 cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename, &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(extendByPs / ir.delta_t));
    }


    // Extending again (tests nsteps not zero initially

    const std::filesystem::path anotherOutTprFilename =
            fileManager.getTemporaryFilePath("extended_again.tpr");

    const std::string secondCommand[] = {
        "convert-tpr", "-s",          outTprFilename.string(), "-o", anotherOutTprFilename.string(),
        "-extend",     extendByString
    };
    CommandLine secondCmdline(secondCommand);
    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &secondCmdline);


    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(anotherOutTprFilename, &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(2 * extendByPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, UntilRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t untilPs       = 100;
    std::string   untilPsString = std::to_string(untilPs);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const std::string           command[]      = {
        "convert-tpr", "-s",         tprFileHandle.tprName(), "-o", outTprFilename.string(),
        "-until",      untilPsString
    };
    CommandLine cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state;
        read_tpx_state(outTprFilename, &ir_after, &state, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + gmx::roundToInt64(untilPs / ir.delta_t));
    }
}

TEST_F(ConvertTprTest, nstepRuntimeExtensionTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    const int64_t originalNStep = ir.nsteps;

    const int64_t nsteps    = 102;
    std::string   nstepsStr = std::to_string(nsteps);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename = fileManager.getTemporaryFilePath("extended.tpr");
    const std::string           command[]      = {
        "convert-tpr", "-s",     tprFileHandle.tprName(), "-o", outTprFilename.string(),
        "-nsteps",     nstepsStr
    };
    CommandLine cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename, &ir_after, &state_after, &top_after);

        EXPECT_EQ(ir_after.nsteps, originalNStep + nsteps);
    }
}

TEST_F(ConvertTprTest, generateVelocitiesTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename =
            fileManager.getTemporaryFilePath("new_velocities.tpr");
    const std::string command[] = { "convert-tpr",           "-s",
                                    tprFileHandle.tprName(), "-o",
                                    outTprFilename.string(), "-generate_velocities",
                                    "-velocity_temp",        "300",
                                    "-velocity_seed",        "12345" };
    CommandLine       cmdline(command);

    gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline);

    {
        gmx_mtop_t top_after;
        t_inputrec ir_after;
        t_state    state_after;
        read_tpx_state(outTprFilename, &ir_after, &state_after, &top_after);

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

TEST_F(ConvertTprNoVelocityTest, refuseToGenerateVelocitiesWhenTprDidNotHaveVelocitiesInitiallyTest)
{
    gmx_mtop_t top;
    t_inputrec ir;
    t_state    state;
    read_tpx_state(tprFileHandle.tprName(), &ir, &state, &top);

    TestFileManager             fileManager;
    const std::filesystem::path outTprFilename =
            fileManager.getTemporaryFilePath("new_velocities.tpr");
    const std::string command[] = { "convert-tpr",           "-s",
                                    tprFileHandle.tprName(), "-o",
                                    outTprFilename.string(), "-generate_velocities",
                                    "-velocity_temp",        "300",
                                    "-velocity_seed",        "12345" };
    CommandLine       cmdline(command);

    GMX_EXPECT_DEATH_IF_SUPPORTED(
            gmx::test::CommandLineTestHelper::runModuleFactory(&gmx::ConvertTprInfo::create, &cmdline),
            "does not contain velocities");
}

} // namespace
} // namespace test
} // namespace gmx
