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
 * \todo Share mdrun input handling implementation via modernized modular options framework.
 * Initial implementation of `launch` relies on borrowed code from the mdrun command
 * line input processing.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 * \ingroup gmxapi
 */

#include "gmxapi/context.h"

#include <cstring>

#include <memory>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrun/logging.h"
#include "gromacs/mdrun/multisim.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/session.h"
#include "gmxapi/version.h"

#include "context-impl.h"
#include "session-impl.h"

namespace gmxapi
{

ContextImpl::ContextImpl()
{
    GMX_ASSERT(session_.expired(), "This implementation assumes an expired weak_ptr at initialization.");
}

std::shared_ptr<gmxapi::ContextImpl> ContextImpl::create()
{
    std::shared_ptr<ContextImpl> impl = std::make_shared<ContextImpl>();
    return impl;
}

std::shared_ptr<Session> ContextImpl::launch(std::string filename)
{
    using namespace gmx;
    // Much of this implementation is not easily testable: we need tools to inspect simulation results and to modify
    // simulation inputs.

    std::shared_ptr<Session> session = nullptr;

    // This implementation can only run one workflow at a time.
    if (session_.expired())
    {
        // Check workflow spec, build graph for current context, launch and return new session.

        /* As default behavior, automatically extend trajectories from the checkpoint file.
         * In the future, our API for objects used to initialize a simulation needs to address the fact that currently a
         * microstate requires data from both the TPR and checkpoint file to be fully specified. Put another way,
         * current
         * GROMACS simulations can take a "configuration" as input that does not constitute a complete microstate in
         * terms of hidden degrees of freedom (integrator/thermostat/barostat/PRNG state), but we want a clear notion of
         * a microstate for gmxapi interfaces.
         */

        // Note: these output options normalize the file names, but not their
        // paths. gmxapi 0.0.7 changes working directory for each session, so the
        // relative paths are appropriate, but a near-future version will avoid
        // changing directories once the process starts and manage file paths explicitly.
        using gmxapi::c_majorVersion;
        using gmxapi::c_minorVersion;
        using gmxapi::c_patchVersion;
        static_assert(!(c_majorVersion != 0 || c_minorVersion != 0 || c_patchVersion > 7),
                      "Developer notice: check assumptions about working directory and relative file paths for this "
                      "software version.");

        // Set input TPR name
        mdArgs_.emplace_back("-s");
        mdArgs_.emplace_back(filename);
        // Set checkpoint file name
        mdArgs_.emplace_back("-cpi");
        mdArgs_.emplace_back("state.cpt");

        // Create a mock argv. Note that argv[0] is expected to hold
        // the program name. We should not try to use
        // std::string.data() because it might not be a valid C
        // pointer, e.g. if small-string optimizations are
        // implemented.
        std::vector < std::vector < char>> argvHolder(1);
        // argv[0] is ignored, but should be a valid string (e.g. null terminated array of char)
        argvHolder[0] = {'\0'};
        for (const auto &element : mdArgs_)
        {
            argvHolder.emplace_back(element.begin(), element.end());
        }
        std::vector<const char *> argv;
        for (const auto &element : argvHolder)
        {
            argv.emplace_back(const_cast<char *>(element.data()));
        }

        auto newMdRunner = compat::make_unique<gmx::Mdrunner>();
        if (!newMdRunner->prepare(static_cast<int>(argv.size()), const_cast<char **>(argv.data())))
        {
            return nullptr;
        }
        SimulationContext simulationContext(newMdRunner->cr);
        {
            auto newSession = SessionImpl::create(shared_from_this(),
                                                  std::move(newMdRunner),
                                                  std::move(simulationContext),
                                                  &newMdRunner->fplog,
                                                  newMdRunner->ms);
            // TO DO: capture fplog and ms
            // TO DO: check builder docs to see if we need to manage anything else.
            session = std::make_shared<Session>(std::move(newSession));
        }

    }
    else
    {
        throw gmxapi::ProtocolError("Tried to launch a session while a session is still active.");
    }

    if (session != nullptr)
    {
        // Update weak reference.
        session_ = session;
    }
    return session;
}

// As of gmxapi 0.0.3 there is only one Context type
Context::Context() :
    Context {ContextImpl::create()}
{
    GMX_ASSERT(impl_, "Context requires a non-null implementation member.");
}

std::shared_ptr<Session> Context::launch(std::string filename)
{
    return impl_->launch(std::move(filename));
}

Context::Context(std::shared_ptr<ContextImpl> impl) :
    impl_ {std::move(impl)}
{
    GMX_ASSERT(impl_, "Context requires a non-null implementation member.");
}

void Context::setMDArgs(const MDArgs &mdArgs)
{
    impl_->mdArgs_ = mdArgs;
}

Context::~Context() = default;

} // end namespace gmxapi
