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
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/exceptions.h"

#include "programs/mdrun/tests/moduletest.h"

#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md/mdsignals.h"
#include "gmxapi/session.h"
#include "gmxapi/session/resources.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "testingconfiguration.h"

namespace gmxapi
{
class SessionResources;

namespace testing
{

namespace
{

/*!
 * \brief Restraint that can optionally issue an immediate stop signal.
 */
class StopSignalIssuer : public gmx::IRestraintPotential
{
public:
    /*!
     * \brief Construct a restraint that does nothing.
     */
    StopSignalIssuer() : StopSignalIssuer(false) {}

    /*!
     * \brief Choose whether or not to issue stop signal when called.
     *
     * \param sendStopSignal If true, issue stop signal at every opportunity.
     */
    explicit StopSignalIssuer(bool sendStopSignal) : sendStopSignal_{ sendStopSignal } {}

    /*! \cond Implement IRestraintPotential */
    gmx::PotentialPointData evaluate(gmx::Vector /* r_site */, gmx::Vector /*  r_ref */, double t) override
    {
        // Note that evaluate gets called once for each site,
        // which is twice per time step for a pair restraint.
        // The following initialization logic is not atomic, but it is sufficient.
        if (!isInitialized())
        {
            // Force is also calculated for initial step.
            simulationStartTime_ = t;
        }
        lastSimulationTime_ = t;

        if (sendStopSignal_)
        {
            auto signalSender = gmxapi::getMdrunnerSignal(resources_, gmxapi::md::signals::STOP);
            signalSender();
        }

        return { { 0., 0., 0. }, 0. };
    }

    std::vector<int> sites() const override { return { { 0, 1 } }; }

    void bindSession(gmxapi::SessionResources* resources) override { resources_ = resources; }
    //! \endcond

    /*!
     * \brief Note simulation start time when called on the zeroeth step.
     */
    double simulationStartTime_ = 0.;

    /*!
     * \brief Record the simulation time at the last step active.
     */
    std::optional<double> lastSimulationTime_;

    /*!
     * \brief Whether restraint was ever used
     */
    bool isInitialized() const { return lastSimulationTime_.has_value(); }

private:
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
    SimpleSignalingClient() : restraint_(std::make_shared<StopSignalIssuer>()) {}

    explicit SimpleSignalingClient(bool sendStopSignal) :
        restraint_(std::make_shared<StopSignalIssuer>(sendStopSignal))
    {
    }

    const char* name() const override { return "SimpleSignalingClient"; }

    std::shared_ptr<gmx::IRestraintPotential> getRestraint() override { return restraint_; }
    //! \endcond

    /*!
     * \brief Last time this restraint was active, minus the simulation start time.
     *
     * \return Time elapsed since start.
     */
    double timeElapsedSinceStart() const
    {
        if (!restraint_->isInitialized())
        {
            GMX_THROW(gmx::InternalError("timeElapsedSinceStart called before restraint was used"));
        }
        return restraint_->lastSimulationTime_.value() - restraint_->simulationStartTime_;
    }

private:
    //! restraint to provide to client or MD simulator.
    std::shared_ptr<StopSignalIssuer> restraint_;
};

/*!
 * \brief Check that we can bind to and use the stop signaler.
 */
TEST_F(GmxApiTest, ApiRunnerStopSignalClient)
{
    const int nsteps = 4;
    makeTprFile(nsteps);

    // Check assumptions about basic simulation behavior.
    {
        const int nstlist = 1;

        auto system  = gmxapi::fromTprFile(runner_.tprFileName_);
        auto context = std::make_shared<gmxapi::Context>(gmxapi::createContext());

        gmxapi::MDArgs args = makeMdArgs();
        args.emplace_back("-nstlist");
        args.emplace_back(std::to_string(nstlist));

        context->setMDArgs(args);

        auto restraint = std::make_shared<SimpleSignalingClient>();

        auto session = system.launch(context);
        EXPECT_TRUE(session);

        gmxapi::addSessionRestraint(session.get(), restraint);
        EXPECT_THROW(restraint->timeElapsedSinceStart(), gmx::InternalError);

        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        EXPECT_EQ(nsteps, gmx::roundToInt(restraint->timeElapsedSinceStart() / getTestStepSize()));

        status = session->close();
        EXPECT_TRUE(status.success());
    }

    // Make sure that stop signal shortens simulation.
    {
        /* StopHandler promises to stop a simulation at the next NS step after the signal got communicated.
         * We don't know the communication interval, but we know that it is at most nstlist. We cannot assume
         * that the signal gets communicated on the step it is set, even if that step is a communication step.
         * As the signal is set on the first step, we know that the restraint will be called at
         * most 2*nstlist + 1 times.
         * Since the time elapsed after the first step is 0, however, we expect the elapsed time
         * divided by the step size to be at most 2*nstlist.
         */

        const int           nstlist  = 1;
        constexpr const int maxsteps = nstlist * 2 + 1;
        // This test is meaningless if the the simulation ends early without a signal.
        static_assert(
                maxsteps < nsteps,
                "Simulation is already scheduled to end before it can receive a stop signal.");

        auto system  = gmxapi::fromTprFile(runner_.tprFileName_);
        auto context = std::make_shared<gmxapi::Context>(gmxapi::createContext());

        gmxapi::MDArgs args = makeMdArgs();
        args.emplace_back("-nstlist");
        args.emplace_back(std::to_string(nstlist));
        // TODO (Ref #3256) use api functionality to extend simulation instead
        args.emplace_back("-nsteps");
        args.emplace_back(std::to_string(nsteps));

        context->setMDArgs(args);

        const bool issueImmediateStopSignal = true;
        auto       restraint = std::make_shared<SimpleSignalingClient>(issueImmediateStopSignal);

        auto session = system.launch(context);
        EXPECT_TRUE(session);

        gmxapi::addSessionRestraint(session.get(), restraint);
        EXPECT_THROW(restraint->timeElapsedSinceStart(), gmx::InternalError);

        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());

        const int steps_just_run =
                gmx::roundToInt(restraint->timeElapsedSinceStart() / getTestStepSize());
        EXPECT_LT(steps_just_run, maxsteps);

        status = session->close();
        EXPECT_TRUE(status.success());
    }
}

} // end anonymous namespace

} // end namespace testing

} // end namespace gmxapi
