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
#include <memory>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/session/resources.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md/mdsignals.h"

#include "testingconfiguration.h"

namespace
{

//! Input file for testing is built by CMake script and filename is compiled into in testingconfiguration binary.
const auto &filename = gmxapi::testing::sample_tprfilename;

/*!
 * \brief Restraint that can optionally issue an immediate stop signal.
 */
class StopSignalIssuer : public gmx::IRestraintPotential
{
    public:
        /*!
         * \brief Construct a restraint that does nothing.
         */
        StopSignalIssuer() :
            StopSignalIssuer(false)
        {}

        /*!
         * \brief Choose whether or not to issue stop signal when called.
         *
         * \param sendStopSignal If true, issue stop signal at every opportunity.
         */
        explicit StopSignalIssuer(bool sendStopSignal) :
            sendStopSignal_ {sendStopSignal}
        {}

        /*! \cond Implement IRestraintPotential */
        gmx::PotentialPointData evaluate(gmx::Vector r_site,
                                         gmx::Vector r_ref,
                                         double      t) override
        {
            (void)r_site;
            (void)r_ref;
            (void)t;
            // Note that evaluate gets called once for each site,
            // which is twice per time step for a pair restraint.
            // The following initialization logic is not atomic, but it is sufficient.
            if (!isInitialized_)
            {
                // Force is also calculated for initial step.
                simulationStartTime_ = t;
                isInitialized_       = true;
            }
            lastSimulationTime_ = t;

            if (sendStopSignal_)
            {
                auto signalSender = gmxapi::getMdrunnerSignal(resources_,
                                                              gmxapi::md::signals::STOP);
                signalSender();
            }

            return {{0., 0., 0.}, 0.};
        }

        std::vector<int> sites() const override
        {
            return {{0, 1}};
        }

        void bindSession(gmxapi::SessionResources* resources) override
        {
            resources_ = resources;
        }
        //! \endcond

        /*!
         * \brief Note simulation start time when called on the zeroeth step.
         */
        double simulationStartTime_ = 0.;

        /*!
         * \brief Record the simulation time at the last step active.
         */
        double lastSimulationTime_ = 0.;

    private:
        /*!
         * \brief Whether to consider current simulation time to be the start time.
         */
        bool isInitialized_ = false;

        /*!
         * \brief Handle through which to get signalling resources.
         */
        gmxapi::SessionResources* resources_ = nullptr;

        /*!
         * \brief Whether to issue stop signal when called.
         */
        bool sendStopSignal_ = false;
};

/*!
 * \brief Wrap a StopSignalIssuer for testing purposes.
 */
class SimpleSignalingClient : public gmxapi::MDModule
{
    public:
        /*! \cond
         * Implement gmxapi::MDModule interface.
         */
        SimpleSignalingClient() :
            restraint_(std::make_shared<StopSignalIssuer>())
        {}

        explicit SimpleSignalingClient(bool sendStopSignal) :
            restraint_(std::make_shared<StopSignalIssuer>(sendStopSignal))
        {}

        const char *name() const override
        {
            return "SimpleSignalingClient";
        }

        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            return restraint_;
        }
        //! \endcond

        /*!
         * \brief Number of steps in which this restraint was active.
         *
         * \return Number of MD time steps.
         */
        int numberOfTimesCalled() const
        {
            const auto timeElapsed =
                restraint_->lastSimulationTime_ - restraint_->simulationStartTime_;

            const auto numSteps    = timeElapsed / gmxapi::testing::testingTimestep;
            return gmx::roundToInt(numSteps);
        }

    private:
        //! restraint to provide to client or MD simulator.
        std::shared_ptr<StopSignalIssuer> restraint_;
};

/*!
 * \brief Check that we can bind to and use the stop signaler.
 */
TEST(ApiRunner, StopSignalClient)
{

    auto system  = gmxapi::fromTprFile(filename);
    auto context = std::make_shared<gmxapi::Context>();

    // Check assumptions about basic simulation behavior.
    {
        gmxapi::MDArgs args    = gmxapi::testing::mdArgs;
        args.emplace_back("-nsteps");
        args.emplace_back("3");
        args.emplace_back("-nstlist");
        args.emplace_back("1");
        // Work around unclean working directory.
        args.emplace_back("-noappend");

        context->setMDArgs(args);

        auto           restraint = std::make_shared<SimpleSignalingClient>();

        auto           session = system.launch(context);
        EXPECT_TRUE(session);

        gmxapi::addSessionRestraint(session.get(), restraint);
        EXPECT_EQ(0, restraint->numberOfTimesCalled());

        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        EXPECT_EQ(3, restraint->numberOfTimesCalled());

        status = session->close();
        EXPECT_TRUE(status.success());
    }

    // Make sure that stop signal shortens simulation.
    {
        gmxapi::MDArgs args    = gmxapi::testing::mdArgs;
        args.emplace_back("-nsteps");
        args.emplace_back("3");
        args.emplace_back("-nstlist");
        args.emplace_back("1");
        // Work around unclean working directory.
        args.emplace_back("-noappend");

        context->setMDArgs(args);

        const bool issueImmediateStopSignal = true;
        auto       restraint                = std::make_shared<SimpleSignalingClient>(issueImmediateStopSignal);

        auto       session = system.launch(context);
        EXPECT_TRUE(session);

        gmxapi::addSessionRestraint(session.get(), restraint);
        EXPECT_EQ(0, restraint->numberOfTimesCalled());

        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        EXPECT_EQ(1, restraint->numberOfTimesCalled());

        status = session->close();
        EXPECT_TRUE(status.success());
    }
}

} // end anonymous namespace
