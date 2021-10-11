/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2020,2021, by the GROMACS development team, led by
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
#include <memory>

#include "testingconfiguration.h"
#include "gmxapi/context.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/iforceprovider.h"

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

    {
        auto           context = std::make_shared<gmxapi::Context>(gmxapi::createContext());
        gmxapi::MDArgs args    = makeMdArgs();
        // TODO the command line arguments should be set through the
        // usual command line options settings for the tests

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
 * \brief Test simulation continuation.
 *
 * Run a simulation, then extend the target number of steps and continue the simulation.
 */
TEST_F(GmxApiTest, RunnerContinuedMD)
{
    // Run a simulation, then extend the target number of steps and continue the simulation
    makeTprFile(10);
    auto system = gmxapi::fromTprFile(runner_.tprFileName_);

    {
        auto context = std::make_shared<gmxapi::Context>(gmxapi::createContext());

        {
            EXPECT_TRUE(context != nullptr);
            gmxapi::MDArgs args = makeMdArgs();

            context->setMDArgs(args);
            auto session = system.launch(context);
            EXPECT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            EXPECT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            EXPECT_TRUE(status.success());
        }

        // Reuse the context. Add MD parameters. Run a new session extending the previous trajectory.
        {
            gmxapi::MDArgs args = makeMdArgs();
            // TODO This needs to be changed to make a new tpr with convert-tpr instead.
            args.emplace_back("-nsteps");
            args.emplace_back("20");

            context->setMDArgs(args);
            auto session = system.launch(context);
            EXPECT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            EXPECT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            EXPECT_TRUE(status.success());
        }
    }
}

} // end anonymous namespace

} // end namespace testing

} // end namespace gmxapi
