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

#include <algorithm>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/legacymdrunoptions.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrun/simulationinputhandle.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/unique_cptr.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/mpi/resourceassignment.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/version.h"

#include "context_impl.h"
#include "createsession.h"
#include "session_impl.h"
#include "workflow.h"

namespace gmxapi
{

// Support some tag dispatch. Warning: These are just aliases (not strong types).
/*!
 * \brief Logical helper alias to convert preprocessor constant to type.
 *
 * \tparam Value Provide the GMX\_LIB\_MPI macro.
 */
template<bool Value>
using hasLibraryMpi = std::bool_constant<Value>;
/* Note that a no-MPI build still uses the tMPI headers to define MPI_Comm for the
 * gmx::SimulationContext definition. The dispatching in this file accounts for
 * these two definitions of SimulationContext. gmxThreadMpi here does not imply
 * that the library was necessarily compiled with thread-MPI enabled.
 */
using gmxThreadMpi = hasLibraryMpi<false>;
using gmxLibMpi    = hasLibraryMpi<true>;
using MpiType      = std::conditional_t<GMX_LIB_MPI, gmxLibMpi, gmxThreadMpi>;

using MpiContextInitializationError = BasicException<struct MpiContextInitialization>;


/*!
 * \brief Helpers to evaluate the correct precondition for the library build.
 *
 * TODO: (#3650) Consider distinct MpiContextManager types for clearer definition of preconditions.
 */
namespace
{

[[maybe_unused]] MPI_Comm validCommunicator(const MPI_Comm& communicator, const gmxThreadMpi&)
{
    if (communicator != MPI_COMM_NULL)
    {
        throw MpiContextInitializationError(
                "Provided communicator must be MPI_COMM_NULL for GROMACS built without MPI "
                "library.");
    }
    return communicator;
}

[[maybe_unused]] MPI_Comm validCommunicator(const MPI_Comm& communicator, const gmxLibMpi&)
{
    if (communicator == MPI_COMM_NULL)
    {
        throw MpiContextInitializationError("MPI-enabled GROMACS requires a valid communicator.");
    }
    return communicator;
}

/*!
 * \brief Return the communicator if it is appropriate for the environment.
 *
 * \throws MpiContextInitializationError if communicator does not match the
 *  MpiContextManager precondition for the current library configuration.
 */
MPI_Comm validCommunicator(const MPI_Comm& communicator)
{
    return validCommunicator(communicator, MpiType());
}

//! \brief Provide a reasonable default value.
MPI_Comm validCommunicator()
{
    return GMX_LIB_MPI ? MPI_COMM_WORLD : MPI_COMM_NULL;
}

} // anonymous namespace

MpiContextManager::MpiContextManager(MPI_Comm communicator) :
    communicator_(std::make_unique<MPI_Comm>(validCommunicator(communicator)))
{
    // Safely increments the GROMACS MPI initialization counter after checking
    // whether the MPI library is already initialized. After this call, MPI_Init
    // or MPI_Init_thread has been called exactly once.
    gmx::init(nullptr, nullptr);
    GMX_RELEASE_ASSERT(!GMX_LIB_MPI || gmx_mpi_initialized(),
                       "MPI should be initialized before reaching this point.");
    if (this->communicator() != MPI_COMM_NULL)
    {
        // Synchronise at the point of acquiring a MpiContextManager.
        gmx_barrier(this->communicator());
    }
}

MpiContextManager::~MpiContextManager()
{
    if (communicator_)
    {
        // This is always safe to call. It is a no-op if
        // thread-MPI, and if the constructor completed then the
        // MPI library is initialized with reference counting.
        gmx::finalize();
    }
}

MpiContextManager::MpiContextManager() : MpiContextManager(validCommunicator()) {}

MPI_Comm MpiContextManager::communicator() const
{
    if (!communicator_)
    {
        throw UsageError("Invalid MpiContextManager. Accessed after `move`?");
    }
    return *communicator_;
}

ContextImpl::~ContextImpl() = default;

[[maybe_unused]] static Context createContext(const ResourceAssignment& resources, const gmxLibMpi&)
{
    CommHandle handle;
    resources.applyCommunicator(&handle);
    if (handle.communicator == MPI_COMM_NULL)
    {
        throw UsageError("MPI-enabled Simulator contexts require a valid communicator.");
    }
    auto contextmanager = MpiContextManager(handle.communicator);
    auto impl           = ContextImpl::create(std::move(contextmanager));
    GMX_ASSERT(impl, "ContextImpl creation method should not be able to return null.");
    auto context = Context(impl);
    return context;
}

[[maybe_unused]] static Context createContext(const ResourceAssignment& resources, const gmxThreadMpi&)
{
    if (resources.size() > 1)
    {
        throw UsageError("Only one thread-MPI Simulation per Context is supported.");
    }
    // Thread-MPI Context does not yet have a need for user-provided resources.
    // However, see #3650.
    return createContext();
}

Context createContext(const ResourceAssignment& resources)
{
    return createContext(resources, hasLibraryMpi<GMX_LIB_MPI>());
}

Context createContext()
{
    MpiContextManager contextmanager;
    auto              impl = ContextImpl::create(std::move(contextmanager));
    GMX_ASSERT(impl, "ContextImpl creation method should not be able to return null.");
    auto context = Context(impl);
    return context;
}

ContextImpl::ContextImpl(MpiContextManager&& mpi) noexcept(std::is_nothrow_constructible_v<gmx::LegacyMdrunOptions>) :
    mpi_(std::move(mpi)),
    hardwareInformation_(gmx_detect_hardware(
            gmx::PhysicalNodeCommunicator(mpi_.communicator(), gmx_physicalnode_id_hash()),
            mpi_.communicator()))
{
    // Confirm our understanding of the MpiContextManager invariant.
    GMX_ASSERT(mpi_.communicator() == MPI_COMM_NULL ? !GMX_LIB_MPI : GMX_LIB_MPI,
               "Precondition violated: inappropriate communicator for the library environment.");
    // Make sure we didn't change the data members and overlook implementation details.
    GMX_ASSERT(session_.expired(),
               "This implementation assumes an expired weak_ptr at initialization.");
}

std::shared_ptr<ContextImpl> ContextImpl::create(MpiContextManager&& mpi)
{
    std::shared_ptr<ContextImpl> impl;
    impl.reset(new ContextImpl(std::move(mpi)));
    return impl;
}

std::shared_ptr<Session> ContextImpl::launch(const Workflow& work)
{
    using namespace gmx;
    // Much of this implementation is not easily testable: we need tools to inspect simulation
    // results and to modify simulation inputs.

    std::shared_ptr<Session> launchedSession = nullptr;

    // This implementation can only run one workflow at a time.
    // Check whether we are already aware of an active session.
    if (session_.expired())
    {
        // Check workflow spec, build graph for current context, launch and return new session.
        // \todo This is specific to the session implementation...
        auto        mdNode = work.getNode("MD");
        std::string filename{};
        if (mdNode != nullptr)
        {
            filename = mdNode->params();
        }

        /* Mock up the argv interface used by option processing infrastructure.
         *
         * As default behavior, automatically extend trajectories from the checkpoint file.
         * In the future, our API for objects used to initialize a simulation needs to address the fact that currently a
         * microstate requires data from both the TPR and checkpoint file to be fully specified. Put another way,
         * current
         * GROMACS simulations can take a "configuration" as input that does not constitute a complete microstate in
         * terms of hidden degrees of freedom (integrator/thermostat/barostat/PRNG state), but we want a clear notion of
         * a microstate for gmxapi interfaces.
         *
         * TODO: Remove `-s` and `-cpi` arguments.
         *       Ref: https://gitlab.com/gromacs/gromacs/-/issues/3652
         */

        // Set input TPR name
        if (std::any_of(mdArgs_.cbegin(), mdArgs_.cend(), [](const std::string& arg) {
                return arg == "-s";
            }))
        {
            throw UsageError("gmxapi must control the simulation input, but caller provided '-s'.");
        }
        mdArgs_.emplace_back("-s");
        mdArgs_.emplace_back(filename);

        // Set checkpoint file name(s) (if not already set by user).
        if (std::none_of(mdArgs_.cbegin(), mdArgs_.cend(), [](const std::string& arg) {
                return arg == "-cpi";
            }))
        {
            mdArgs_.emplace_back("-cpi");
            mdArgs_.emplace_back("state.cpt");
        }
        if (std::none_of(mdArgs_.cbegin(), mdArgs_.cend(), [](const std::string& arg) {
                return arg == "-cpo";
            }))
        {
            mdArgs_.emplace_back("-cpo");
            mdArgs_.emplace_back("state.cpt");
        }
        if (std::none_of(mdArgs_.cbegin(), mdArgs_.cend(), [](const std::string& arg) {
                return arg == "-o";
            }))
        {
            mdArgs_.emplace_back("-o");
            mdArgs_.emplace_back("traj.trr");
        }
        /* Note: we normalize the file names, but not the full paths.
         * A future gmxapi version should manage files (paths) in an absolute and
         * immutable way, with abstraction provided through the Context chain-of-responsibility.
         * TODO: API abstractions for initializing simulations that may be new or partially
         * complete. Reference gmxapi milestone 13 at https://gitlab.com/gromacs/gromacs/-/issues/2585
         */

        // Create a mock argv. Note that argv[0] is expected to hold the program name.
        const int  offset = 1;
        const auto argc   = static_cast<size_t>(mdArgs_.size() + offset);
        auto       argv   = std::vector<char*>(argc, nullptr);
        // argv[0] is ignored, but should be a valid string (e.g. null terminated array of char)
        argv[0]  = new char[1];
        *argv[0] = '\0';
        for (size_t argvIndex = offset; argvIndex < argc; ++argvIndex)
        {
            const auto&  mdArg   = mdArgs_[argvIndex - offset];
            const size_t argSize = mdArg.length() + 1;
            argv[argvIndex]      = new char[argSize];
            strncpy(argv[argvIndex], mdArg.c_str(), argSize);
        }

        auto mdModules = std::make_unique<MDModules>();

        const char* desc[] = { "gmxapi placeholder text" };

        // LegacyMdrunOptions needs to be kept alive for the life of ContextImpl,
        // so we use a data member for now.
        gmx::LegacyMdrunOptions& options = options_;
        if (options.updateFromCommandLine(argc, argv.data(), desc) == 0)
        {
            return nullptr;
        }

        ArrayRef<const std::string> multiSimDirectoryNames = opt2fnsIfOptionSet(
                "-multidir", gmx::ssize(options.filenames), options.filenames.data());


        // The SimulationContext is necessary with gmxapi so that
        // resources owned by the client code can have suitable
        // lifetime. The gmx wrapper binary uses the same infrastructure,
        // but the lifetime is now trivially that of the invocation of the
        // wrapper binary.
        //
        // For now, this should match the communicator used for hardware
        // detection. There's no way to assert this is true.
        auto communicator = mpi_.communicator();
        // Confirm the precondition for simulationContext().
        GMX_ASSERT(communicator == MPI_COMM_NULL ? !GMX_LIB_MPI : GMX_LIB_MPI,
                   "Context communicator does not have an appropriate value for the environment.");
        SimulationContext simulationContext(communicator, multiSimDirectoryNames);


        StartingBehavior startingBehavior = StartingBehavior::NewSimulation;
        LogFilePtr       logFileGuard     = nullptr;
        gmx_multisim_t*  ms               = simulationContext.multiSimulation_.get();
        std::tie(startingBehavior, logFileGuard) =
                handleRestart(findIsSimulationMainRank(ms, simulationContext.simulationCommunicator_),
                              simulationContext.simulationCommunicator_,
                              ms,
                              options.mdrunOptions.appendingBehavior,
                              gmx::ssize(options.filenames),
                              options.filenames.data());

        auto builder = MdrunnerBuilder(std::move(mdModules),
                                       compat::not_null<SimulationContext*>(&simulationContext));
        builder.addHardwareDetectionResult(hardwareInformation_.get());
        builder.addSimulationMethod(options.mdrunOptions, options.pforce, startingBehavior);
        builder.addDomainDecomposition(options.domdecOptions);
        // \todo pass by value
        builder.addNonBonded(options.nbpu_opt_choices[0]);
        // \todo pass by value
        builder.addElectrostatics(options.pme_opt_choices[0], options.pme_fft_opt_choices[0]);
        builder.addBondedTaskAssignment(options.bonded_opt_choices[0]);
        builder.addUpdateTaskAssignment(options.update_opt_choices[0]);
        builder.addNeighborList(options.nstlist_cmdline);
        builder.addReplicaExchange(options.replExParams);
        // Need to establish run-time values from various inputs to provide a resource handle to Mdrunner
        builder.addHardwareOptions(options.hw_opt);

        // \todo File names are parameters that should be managed modularly through further factoring.
        builder.addFilenames(options.filenames);
        // TODO: Remove `s` and `-cpi` from LegacyMdrunOptions before launch(). #3652
        auto simulationInput = makeSimulationInput(options);
        builder.addInput(simulationInput);

        // Note: The gmx_output_env_t life time is not managed after the call to parse_common_args.
        // \todo Implement lifetime management for gmx_output_env_t.
        // \todo Output environment should be configured outside of Mdrunner and provided as a resource.
        builder.addOutputEnvironment(options.oenv);
        builder.addLogFile(logFileGuard.get());

        // Note, creation is not mature enough to be exposed in the external API yet.
        launchedSession = createSession(
                shared_from_this(), std::move(builder), std::move(simulationContext), std::move(logFileGuard));

        // Clean up argv once builder is no longer in use
        // TODO: Remove long-lived references to argv so this is no longer necessary.
        //       Ref https://gitlab.com/gromacs/gromacs/-/issues/2877
        for (auto&& string : argv)
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

std::shared_ptr<Session> Context::launch(const Workflow& work)
{
    return impl_->launch(work);
}

Context::Context(std::shared_ptr<ContextImpl> impl) : impl_{ std::move(impl) }
{
    if (!impl_)
    {
        throw UsageError("Context requires a non-null implementation member.");
    }
}

void Context::setMDArgs(const MDArgs& mdArgs)
{
    impl_->mdArgs_ = mdArgs;
}

Context::~Context() = default;

} // end namespace gmxapi
