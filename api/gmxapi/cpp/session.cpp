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

#include "gmxapi/session.h"

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>

#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/init.h"

#include "gmxapi/context.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/md/mdmodule.h"
#include "gmxapi/md/mdsignals.h"
#include "gmxapi/status.h"

#include "createsession.h"
#include "mdsignals.h"
#include "session_impl.h"
#include "sessionresources.h"

namespace gmxapi
{
class ContextImpl;

SignalManager::SignalManager(gmx::StopHandlerBuilder* stopHandlerBuilder) :
    state_(std::make_shared<gmx::StopSignal>(gmx::StopSignal::noSignal))
{

    /*!
     * \brief Signal issuer managed by this object.
     *
     * Created and bound when the runner is built. Subsequently, client
     * stop signals are proxied to the simulation through the session
     * resources. The MD internal signal handler calls this functor
     * during the MD loop to see if the simulation should be stopped.
     * Thus, this should execute within a very small fraction of an MD
     * step and not require any synchronization.
     */
    auto currentState     = state_;
    auto stopSignalIssuer = [currentState]() { return *currentState; };
    stopHandlerBuilder->registerStopCondition(stopSignalIssuer);
}

//! \cond
SignalManager::~SignalManager() = default;
//! \endcond

bool SessionImpl::isOpen() const noexcept
{
    // Currently, an active session is equivalent to an active Mdrunner.
    return bool(runner_);
}

Status SessionImpl::close()
{
    // Assume unsuccessful until proven otherwise.
    auto successful = Status(false);

    // When the Session is closed, we need to know that the MD output has been finalized.
    runner_.reset();
    logFilePtr_.reset();

    // \todo provide meaningful result.
    // We should be checking that resources were properly shut down, but
    // there isn't currently a way to do that.
    successful = true;
    return successful;
}

Status SessionImpl::run() noexcept
{
    // Status is failure until proven otherwise.
    auto successful = Status(false);
    GMX_ASSERT(runner_, "SessionImpl invariant implies valid Mdrunner handle.");
    auto rc = runner_->mdrunner();
    if (rc == 0)
    {
        successful = true;
    }
    return successful;
}

std::unique_ptr<SessionImpl> SessionImpl::create(std::shared_ptr<ContextImpl> context,
                                                 gmx::MdrunnerBuilder&&       runnerBuilder,
                                                 gmx::SimulationContext&&     simulationContext,
                                                 gmx::LogFilePtr              logFilehandle)
{
    // We should be able to get a communicator (or subcommunicator) through the
    // Context.
    return std::make_unique<SessionImpl>(std::move(context),
                                         std::move(runnerBuilder),
                                         std::move(simulationContext),
                                         std::move(logFilehandle));
}

SessionImpl::SessionImpl(std::shared_ptr<ContextImpl> context,
                         gmx::MdrunnerBuilder&&       runnerBuilder,
                         gmx::SimulationContext&&     simulationContext,
                         gmx::LogFilePtr              fplog) :
    context_(std::move(context)),
    simulationContext_(std::move(simulationContext)),
    logFilePtr_(std::move(fplog))
{
    GMX_ASSERT(context_, "SessionImpl invariant implies valid ContextImpl handle.");

    // \todo Session objects can have logic specialized for the runtime environment.

    auto stopHandlerBuilder = std::make_unique<gmx::StopHandlerBuilder>();
    signalManager_          = std::make_unique<SignalManager>(stopHandlerBuilder.get());
    GMX_ASSERT(signalManager_, "SessionImpl invariant includes a valid SignalManager.");

    runnerBuilder.addStopHandlerBuilder(std::move(stopHandlerBuilder));
    runner_ = std::make_unique<gmx::Mdrunner>(runnerBuilder.build());
    GMX_ASSERT(runner_, "SessionImpl invariant implies valid Mdrunner handle.");

    // For the libgromacs context, a session should explicitly reset global variables that could
    // have been set in a previous simulation during the same process.
    gmx_reset_stop_condition();
}

std::shared_ptr<Session> createSession(std::shared_ptr<ContextImpl> context,
                                       gmx::MdrunnerBuilder&&       runnerBuilder,
                                       gmx::SimulationContext&&     simulationContext,
                                       gmx::LogFilePtr              logFilehandle)
{
    auto newSession      = SessionImpl::create(std::move(context),
                                          std::move(runnerBuilder),
                                          std::move(simulationContext),
                                          std::move(logFilehandle));
    auto launchedSession = std::make_shared<Session>(std::move(newSession));
    return launchedSession;
}

Status SessionImpl::addRestraint(std::shared_ptr<gmxapi::MDModule> module)
{
    GMX_ASSERT(runner_, "SessionImpl invariant implies valid Mdrunner handle.");
    Status status{ false };

    if (module != nullptr)
    {
        const auto& name = module->name();
        if (restraints_.find(name) == restraints_.end())
        {
            auto restraint = module->getRestraint();
            if (restraint != nullptr)
            {
                restraints_.emplace(name, restraint);
                auto sessionResources = createResources(module);
                if (!sessionResources)
                {
                    status = false;
                }
                else
                {
                    runner_->addPotential(restraint, module->name());
                    status = true;
                }
            }
        }
    }
    return status;
}


SignalManager* SessionImpl::getSignalManager()
{
    SignalManager* ptr = nullptr;
    if (isOpen())
    {
        ptr = signalManager_.get();
    }
    return ptr;
}

gmx::Mdrunner* SessionImpl::getRunner()
{
    gmx::Mdrunner* runner = nullptr;
    if (runner_)
    {
        runner = runner_.get();
    }
    return runner;
}

gmxapi::SessionResources* SessionImpl::getResources(const std::string& name) const noexcept
{
    gmxapi::SessionResources* resources = nullptr;
    try
    {
        resources = resources_.at(name).get();
    }
    catch (const std::out_of_range& e)
    {
        // named operation does not have any resources registered.
    }

    return resources;
}

gmxapi::SessionResources* SessionImpl::createResources(std::shared_ptr<gmxapi::MDModule> module) noexcept
{
    // check if resources already exist for this module
    // If not, create resources and return handle.
    // Return nullptr for any failure.
    gmxapi::SessionResources* resources = nullptr;
    const auto&               name      = module->name();
    if (resources_.find(name) == resources_.end())
    {
        auto resourcesInstance = std::make_unique<SessionResources>(this, name);
        resources_.emplace(name, std::move(resourcesInstance));
        resources = resources_.at(name).get();
        // To do: This should be more dynamic.
        getSignalManager()->addSignaller(name);
        if (restraints_.find(name) != restraints_.end())
        {
            auto restraintRef = restraints_.at(name);
            auto restraint    = restraintRef.lock();
            if (restraint)
            {
                restraint->bindSession(resources);
            }
        }
    }
    return resources;
}

Session::Session(std::unique_ptr<SessionImpl> impl) noexcept
{
    if (impl != nullptr)
    {
        impl_ = std::move(impl);
    }
    GMX_ASSERT(impl_->isOpen(), "SessionImpl invariant implies valid Mdrunner handle.");
}

Status Session::run() noexcept
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");

    const Status status = impl_->run();
    return status;
}

Status Session::close()
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");

    auto status = Status(false);
    if (isOpen())
    {
        // \todo catch exceptions when we know what they might be
        status = impl_->close();
    }

    return status;
}

Session::~Session()
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");
    if (isOpen())
    {
        try
        {
            impl_->close();
        }
        catch (const std::exception&)
        {
            // \todo find some exception-safe things to do with this via the Context interface.
        }
    }
}

bool Session::isOpen() const noexcept
{
    GMX_ASSERT(impl_, "Session invariant implies valid implementation object handle.");
    const auto result = impl_->isOpen();
    return result;
}

Status addSessionRestraint(Session* session, std::shared_ptr<gmxapi::MDModule> restraint)
{
    auto status = gmxapi::Status(false);

    if (session != nullptr && restraint != nullptr)
    {
        // \todo Improve the external / library API facets
        // so the public API does not need to offer raw pointers.
        auto sessionImpl = session->getRaw();

        GMX_RELEASE_ASSERT(sessionImpl,
                           "Session invariant implies valid implementation object handle.");
        // GMX_ASSERT alone is not strong enough to convince linters not to warn of possible nullptr.
        if (sessionImpl)
        {
            status = sessionImpl->addRestraint(std::move(restraint));
        }
    }
    return status;
}

//! \cond internal
SessionImpl* Session::getRaw() const noexcept
{
    return impl_.get();
}
//! \endcond

std::shared_ptr<Session> launchSession(Context* context, const Workflow& work) noexcept
{
    auto session = context->launch(work);
    return session;
}

SessionImpl::~SessionImpl() = default;

SessionResources::SessionResources(gmxapi::SessionImpl* session, std::string name) :
    sessionImpl_(session), name_(std::move(name))
{
}

SessionResources::~SessionResources() = default;

std::string SessionResources::name() const
{
    return name_;
}

Signal SessionResources::getMdrunnerSignal(md::signals signal)
{
    //// while there is only one choice...
    if (signal != md::signals::STOP)
    {
        throw gmxapi::MissingImplementationError("This signaller only handles stop signals.");
    }

    // Get a signalling proxy for the caller.
    auto signalManager = sessionImpl_->getSignalManager();
    if (signalManager == nullptr)
    {
        throw gmxapi::ProtocolError(
                "Client requested access to a signaller that is not available.");
    }
    auto functor = signalManager->getSignal(name_, signal);

    return functor;
}

} // end namespace gmxapi
