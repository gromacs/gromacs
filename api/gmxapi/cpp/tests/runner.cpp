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
#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/trajectory/trajectoryframe.h"

#include "testutils/testfilemanager.h"
#include "testutils/trajectoryreader.h"

#include "programs/mdrun/tests/moduletest.h"

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "testingconfiguration.h"

namespace gmxapi
{

namespace testing
{

namespace
{

/*!
 * \brief Check that we can run a basic simulation from a simple client.
 */
TEST_F(GmxApiTest, RunnerBasicMD)
{
    makeTprFile(2);
    auto system = gmxapi::fromTprFile(runner_.tprFileName_);

    // Check our input validation.
    {
        auto           context = std::make_shared<gmxapi::Context>(gmxapi::createContext());
        gmxapi::MDArgs args    = makeMdArgs();
        args.emplace_back("-s");
        args.emplace_back("dummyfilename.tpr");
        context->setMDArgs(args);
        EXPECT_THROW(system.launch(context), gmxapi::UsageError);
    }

    {
        auto           context = std::make_shared<gmxapi::Context>(gmxapi::createContext());
        gmxapi::MDArgs args    = makeMdArgs();

        context->setMDArgs(args);
        auto session = system.launch(context);
        EXPECT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        status = session->close();
        EXPECT_TRUE(status.success());
    }
}

/*!
 * \brief Test our ability to reinitialize the libgromacs environment between simulations.
 */
TEST_F(GmxApiTest, RunnerReinitialize)
{
    auto           context = std::make_shared<gmxapi::Context>(gmxapi::createContext());
    gmxapi::MDArgs args    = makeMdArgs();

    makeTprFile(20);

    {
        context->setMDArgs(args);
        auto system  = gmxapi::fromTprFile(runner_.tprFileName_);
        auto session = system.launch(context);

        // Try to simulate an interrupt signal to catch.
        gmx_set_stop_condition(StopCondition::NextNS);

        auto status = session->run();
        EXPECT_FALSE(status.success());

        // If this assertion fails, it is not an error, but it indicates expected behavior has
        // changed and we need to consider the impact of whatever changes caused this.
        EXPECT_NE(gmx_get_stop_condition(), StopCondition::None);

        session->close();
    } // allow system and session to be destroyed.

    {
        context->setMDArgs(args);
        auto system = gmxapi::fromTprFile(runner_.tprFileName_);

        // If this assertion fails, it is not an error, but it indicates expected behavior has
        // changed and we need to consider the impact of whatever changes caused this.
        // We are expecting that the libgromacs state has retained the stop condition from the
        // previously issued SIGINT
        EXPECT_NE(gmx_get_stop_condition(), StopCondition::None);

        auto session = system.launch(context);

        // Launching a session should clear the stop condition.
        EXPECT_EQ(gmx_get_stop_condition(), StopCondition::None);

        auto status = session->run();
        EXPECT_TRUE(status.success());

        // Stop condition should still be clear.
        EXPECT_EQ(gmx_get_stop_condition(), StopCondition::None);

        session->close();
    }
}

/*!
 * \brief Test chained trajectory segments.
 *
 * Use the result of one simulation as the starting point for an exact continuation.
 */
TEST_F(GmxApiTest, RunnerChainedMD)
{
    const int segment1Steps     = 2;
    const int segment2Steps     = 2;
    const int segment2FinalStep = segment1Steps + segment2Steps;
    makeTprFile(segment1Steps);

    auto context = std::make_shared<gmxapi::Context>(gmxapi::createContext());

    {
        auto system = gmxapi::fromTprFile(runner_.tprFileName_);

        // In the absence of a SimulationResult or SimulationOutput abstraction, we have to use the
        // input parameters to infer the outputs.
        // Note: makeMdArgs may add checkpoint file options that we would want to read from.

        EXPECT_TRUE(context != nullptr);
        gmxapi::MDArgs args = makeMdArgs();
        // makeMdArgs automatically sets output trajectory.
        ASSERT_TRUE(std::any_of(
                args.cbegin(), args.cend(), [](const std::string& arg) { return arg == "-x"; }));
        ASSERT_TRUE(!runner_.reducedPrecisionTrajectoryFileName_.empty());

        context->setMDArgs(args);
        auto session = system.launch(context);
        EXPECT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        ASSERT_NO_THROW(status = session->close());
        EXPECT_TRUE(status.success());
    }

    {
        // Check that we ran the expected number of steps.
        auto reader = gmx::test::TrajectoryFrameReader(runner_.fullPrecisionTrajectoryFileName_);
        int  currentStep = 0;
        while (currentStep < segment1Steps)
        {
            currentStep = reader.frame().step();
        }
        EXPECT_FALSE(reader.readNextFrame());
    }

    // Run a new segment extending the previous trajectory.
    {
        // Re-use the name of the previous TPR file, simplifying
        // clean-up of test files under the gmx::test::MdrunTestFixture scheme.
        runner_.changeTprNsteps(segment2FinalStep);

        auto system = gmxapi::fromTprFile(runner_.tprFileName_);

        gmxapi::MDArgs args = makeMdArgs();

        // Get the checkpoint file from the previous simulation.
        ASSERT_TRUE(std::none_of(
                args.cbegin(), args.cend(), [](const std::string& arg) { return arg == "-cpi"; }));
        ASSERT_TRUE(!runner_.cptOutputFileName_.empty());
        args.emplace_back("-cpi");
        args.emplace_back(runner_.cptOutputFileName_);

        // Provide distinct output.
        args.emplace_back("-noappend");
        // Need to override the value for `-x` set by makeMdArgs.
        auto trajectory_arg = std::find(args.begin(), args.end(), "-x");
        ASSERT_TRUE(trajectory_arg != args.cend());
        ++trajectory_arg;
        // As of this writing, the TestFileManager does not actually interact
        // with the automatic GROMACS simulation file management, and output
        // file names may not be exactly as specified by command line options.
        // The suffix here works around these caveats for a particular structure
        // of this test code, but the workaround is brittle.
        *trajectory_arg = fileManager_.getTemporaryFilePath(".part0002.xtc");

        context->setMDArgs(args);
        auto session = system.launch(context);
        EXPECT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        ASSERT_NO_THROW(status = session->close());
        EXPECT_TRUE(status.success());

        // Check that the new trajectory has the expected steps.
        auto reader      = gmx::test::TrajectoryFrameReader(*trajectory_arg);
        int  currentStep = reader.frame().step();
        EXPECT_GE(currentStep, segment1Steps);
        while (currentStep < segment2FinalStep)
        {
            currentStep = reader.frame().step();
        }
        EXPECT_EQ(currentStep, segment2FinalStep);
        EXPECT_FALSE(reader.readNextFrame());
    }
}

} // end anonymous namespace

} // end namespace testing

} // end namespace gmxapi
