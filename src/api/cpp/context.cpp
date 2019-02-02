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
#include <utility>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/stophandler.h"
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
#include "gmxapi/status.h"
#include "gmxapi/version.h"

#include "context_impl.h"
#include "createsession.h"
#include "session_impl.h"
#include "workflow.h"

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

std::shared_ptr<Session> ContextImpl::launch(const Workflow &work)
{
    using namespace gmx;
    // Much of this implementation is not easily testable: we need tools to inspect simulation results and to modify
    // simulation inputs.

    std::shared_ptr<Session> launchedSession = nullptr;

    // This implementation can only run one workflow at a time.
    // Check whether we are already aware of an active session.
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

        /* As default behavior, automatically extend trajectories from the checkpoint file.
         * In the future, our API for objects used to initialize a simulation needs to address the fact that currently a
         * microstate requires data from both the TPR and checkpoint file to be fully specified. Put another way,
         * current
         * GROMACS simulations can take a "configuration" as input that does not constitute a complete microstate in
         * terms of hidden degrees of freedom (integrator/thermostat/barostat/PRNG state), but we want a clear notion of
         * a microstate for gmxapi interfaces.
         */

        // Set input TPR name
        mdArgs_.emplace_back("-s");
        mdArgs_.emplace_back(filename);

        // Set checkpoint file name
        mdArgs_.emplace_back("-cpi");
        mdArgs_.emplace_back("state.cpt");
        /* Note: we normalize the checkpoint file name, but not its full path.
         * Through version 0.0.8, gmxapi clients change working directory
         * for each session, so relative path(s) below are appropriate.
         * A future gmxapi version should avoid changing directories once the
         * process starts and instead manage files (paths) in an absolute and
         * immutable way, with abstraction provided through the Context chain-of-responsibility.
         * TODO: API abstractions for initializing simulations that may be new or partially complete.
         * Reference gmxapi milestone 13 at https://redmine.gromacs.org/issues/2585
         */

        // Create a mock argv. Note that argv[0] is expected to hold the program name.
        const int  offset = 1;
        const auto argc   = static_cast<size_t>(mdArgs_.size() + offset);
        auto       argv   = std::vector<char *>(argc, nullptr);
        // argv[0] is ignored, but should be a valid string (e.g. null terminated array of char)
        argv[0]  = new char[1];
        *argv[0] = '\0';
        for (size_t argvIndex = offset; argvIndex < argc; ++argvIndex)
        {
            const auto &mdArg = mdArgs_[argvIndex - offset];
            argv[argvIndex] = new char[mdArg.length() + 1];
            strcpy(argv[argvIndex], mdArg.c_str());
        }

        // pointer-to-t_commrec is the de facto handle type for communications record.
        // Complete shared / borrowed ownership requires a reference to this stack variable
        // (or pointer-to-pointer-to-t_commrec) since borrowing code may update the pointer.
        // \todo Define the ownership and lifetime management semantics for a communication record, handle or value type.

        // Note: this communications record initialization acts directly on
        // MPI_COMM_WORLD and is incompatible with MPI environment sharing in
        // gmxapi through 0.0.7, at least.
        options_.cr = init_commrec();

        const char *desc[]  = {"gmxapi placeholder text"};
        if (options_.updateFromCommandLine(argc, argv.data(), desc) == 0)
        {
            return nullptr;
        }

        if (MASTER(options_.cr))
        {
            options_.logFileGuard = openLogFile(ftp2fn(efLOG,
                                                       options_.filenames.size(),
                                                       options_.filenames.data()),
                                                options_.mdrunOptions.continuationOptions.appendFiles);
        }

        auto simulationContext = createSimulationContext(options_.cr);

        auto builder = MdrunnerBuilder(compat::not_null<decltype( &simulationContext)>(&simulationContext));
        builder.addSimulationMethod(options_.mdrunOptions, options_.pforce);
        builder.addDomainDecomposition(options_.domdecOptions);
        // \todo pass by value
        builder.addNonBonded(options_.nbpu_opt_choices[0]);
        // \todo pass by value
        builder.addElectrostatics(options_.pme_opt_choices[0], options_.pme_fft_opt_choices[0]);
        builder.addBondedTaskAssignment(options_.bonded_opt_choices[0]);
        builder.addNeighborList(options_.nstlist_cmdline);
        builder.addReplicaExchange(options_.replExParams);
        // \todo take ownership of multisim resources (ms)
        builder.addMultiSim(options_.ms);
        // \todo Provide parallelism resources through SimulationContext.
        // Need to establish run-time values from various inputs to provide a resource handle to Mdrunner
        builder.addHardwareOptions(options_.hw_opt);
        // \todo File names are parameters that should be managed modularly through further factoring.
        builder.addFilenames(options_.filenames);
        // Note: The gmx_output_env_t life time is not managed after the call to parse_common_args.
        // \todo Implement lifetime management for gmx_output_env_t.
        // \todo Output environment should be configured outside of Mdrunner and provided as a resource.
        builder.addOutputEnvironment(options_.oenv);
        builder.addLogFile(options_.logFileGuard.get());

        // Note, creation is not mature enough to be exposed in the external API yet.
        launchedSession = createSession(shared_from_this(),
                                        std::move(builder),
                                        simulationContext,
                                        std::move(options_.logFileGuard),
                                        options_.ms);

        // Clean up argv once builder is no longer in use
        for (auto && string : argv)
        {
            if (string != nullptr)
            {
                delete[] string;
                string = nullptr;
            }
        }

    }
    else
    {
        throw gmxapi::ProtocolError("Tried to launch a session while a session is still active.");
    }

    if (launchedSession != nullptr)
    {
        // Update weak reference.
        session_ = launchedSession;
    }
    return launchedSession;
}

// As of gmxapi 0.0.3 there is only one Context type
Context::Context() :
    Context {ContextImpl::create()}
{
    GMX_ASSERT(impl_, "Context requires a non-null implementation member.");
}

std::shared_ptr<Session> Context::launch(const Workflow &work)
{
    return impl_->launch(work);
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
