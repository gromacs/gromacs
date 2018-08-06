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
/*! \file
 * \brief Implementation details of gmxapi::Context
 *
 */

#include "gmxapi/context.h"

#include <memory>
#include <utility>

#include <cassert>
#include <vector>
#include "programs/mdrun/runner.h"
#include "gromacs/mdtypes/tpxstate.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"

#include "workflow.h"
#include "workflow-impl.h"
#include "session-impl.h"
#include "gromacs/compat/make_unique.h"

namespace gmxapi
{

/*! \brief Temporary shim until proper messaging.
 *
 */
class warn
{
    public:
        /*! \brief Create a warning message.
         *
         * \param message must outlive warn instance...
         */
        explicit warn(const char* message) :
            message_ {message}
        { };
        /// pointer to string managed somewhere else.
        const char* message_;
};

using gmxapi::MDArgs;

/*!
 * \brief Context implementation base class.
 *
 * Execution contexts have a uniform interface specified by the API. Implementations for
 * particular execution environments can specialize / derive from this base.
 *
 * \todo Separate interface and implementation.
 */
class ContextImpl final : public std::enable_shared_from_this<ContextImpl>
{
    public:
        static std::shared_ptr<gmxapi::ContextImpl> create();

        /*!
         * \brief Default constructor.
         *
         * Don't use this. Use create() to get a shared pointer right away.
         * Otherwise, shared_from_this() is potentially dangerous.
         */
        ContextImpl();

        /*!
         * \brief Get a reference to the current status object.
         *
         * \return shared ownership of the current status.
         */
        std::shared_ptr<const Status> status() const noexcept;

        /*!
         * \brief Translate the workflow to the execution context and launch.
         *
         * \param work workflow graph
         * \return ownership of a new session
         *
         * \todo This probably makes more sense as a free function.
         */
        std::shared_ptr<Session> launch(const Workflow &work);

        /*!
         * \brief Status of the last operation in the local context.
         *
         * This pointer should always be valid while the context handle exists and
         * client code can extend the life of the object. We use a shared_ptr because
         * it may be expensive or dangerous to copy the object when it is most needed.
         */
        std::shared_ptr<Status> status_;
        std::weak_ptr<Session>  session_;

        MDArgs                  mdArgs_;
};

ContextImpl::ContextImpl() :
    status_ {std::make_shared<Status>(true)},
session_ {}
{
    assert(status_->success());
    assert(session_.expired());
}

std::shared_ptr<gmxapi::ContextImpl> ContextImpl::create()
{
    auto impl = std::make_shared<gmxapi::ContextImpl>();
    return impl;
}

std::shared_ptr<const Status> ContextImpl::status() const noexcept
{
    return status_;
}

std::shared_ptr<Session> ContextImpl::launch(const Workflow &work)
{
    // Assume failure until proven otherwise.
    assert(status_ != nullptr);
    *status_ = false;

    // Much of this implementation is not easily testable: we need tools to inspect simulation results and to modify
    // simulation inputs.

    // As default behavior, automatically extend trajectories from the checkpoint file.
    // In the future, our API for objects used to initialize a simulation needs to address the fact that currently a
    // microstate requires data from both the TPR and checkpoint file to be fully specified. Put another way, current
    // GROMACS simulations can take a "configuration" as input that does not constitute a complete microstate in terms
    // of hidden degrees of freedom (integrator/thermostat/barostat/PRNG state), but we want a clear notion of a
    // microstate for gmxapi interfaces.
    mdArgs_.emplace_back("-cpi");
    mdArgs_.emplace_back("state.cpt");

    std::shared_ptr<Session> session {
        nullptr
    };

    // This implementation can only run one workflow at a time.
    if (session_.expired())
    {
        // Check workflow spec, build graph for current context, launch and return new session.
        // \todo This is specific to the session implementation...
        auto        mdNode = work.getNode("MD");
        std::string filename {};
        if (mdNode != nullptr)
        {
            filename = mdNode->params();
        }
        auto newMdRunner = gmx::compat::make_unique<gmx::Mdrunner>();
        if (!filename.empty())
        {
            auto tpxState = gmx::TpxState::initializeFromFile(filename);
            newMdRunner->setTpx(std::move(tpxState));
            newMdRunner->initFromAPI(mdArgs_);
        }

        {
            auto newSession = SessionImpl::create(shared_from_this(),
                                                  std::move(newMdRunner));
            session = std::make_shared<Session>(std::move(newSession));
        }


//        for (auto&& node : work)
//        {
//            // If we can do something with the node, do it. If the spec is bad, error.
//        }
    }
    else
    {
        throw gmxapi::ProtocolError("Tried to launch a session while a session is still active.");
    }

    if (session != nullptr)
    {
        // Update weak reference.
        session_ = session;
        // Set successful status.
        *status_ = true;
    }
    return session;
}

// In 0.0.3 there is only one Context type
Context::Context() :
    Context {std::make_shared<ContextImpl>()}
{
    assert(impl_ != nullptr);
}

std::shared_ptr<Session> Context::launch(const Workflow &work)
{
    return impl_->launch(work);
}

Context::Context(std::shared_ptr<ContextImpl> &&impl) :
    impl_ {std::move(impl)}
{
    assert(impl_ != nullptr);
}

void Context::setMDArgs(const MDArgs &mdArgs)
{
    impl_->mdArgs_ = mdArgs;
}

Context::~Context() = default;

std::unique_ptr<Context> defaultContext()
{
    auto impl    = gmx::compat::make_unique<ContextImpl>();
    auto context = gmx::compat::make_unique<Context>(std::move(impl));
    return context;
}

} // end namespace gmxapi
