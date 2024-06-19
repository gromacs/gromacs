/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief Declares the routine running the inetgrators.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_RUNNER_H
#define GMX_MDRUN_RUNNER_H

#include <cstdio>

#include <array>
#include <memory>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/options.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/simulationinputhandle.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/unique_cptr.h"

#include "replicaexchange.h"

struct gmx_multisim_t;
struct gmx_output_env_t;
struct ReplicaExchangeParameters;
struct t_fileio;

namespace gmx
{

// Todo: move to forward declaration headers...
class MDModules;
class IRestraintPotential; // defined in restraint/restraintpotential.h
class RestraintManager;
class SimulationContext;
class StopHandlerBuilder;

/*! \libinternal \brief Runner object for supporting setup and execution of mdrun.
 *
 * This class has responsibility for the lifetime of data structures
 * that exist for the life of the simulation, e.g. for logging and
 * communication.
 *
 * It is also responsible for initializing data members that
 * e.g. correspond to values potentially set by commmand-line
 * options. Later these will be obtained directly from modules, and
 * the results of command-line option handling returned directly to
 * the modules, rather than propagated to them by data members of this
 * class.
 *
 * \todo Most of the attributes should be declared by specific modules
 * as command-line options. Accordingly, they do not conform to the
 * naming scheme, because that would make for a lot of noise in the
 * diff, only to have it change again when the options move to their
 * modules.
 *
 * \todo Preparing logging and MPI contexts could probably be a
 * higher-level responsibility, so that an Mdrunner would get made
 * without needing to re-initialize these components (as currently
 * happens always for the main rank, and differently for the spawned
 * ranks with thread-MPI).
 *
 * \ingroup module_mdrun
 */
class Mdrunner
{
public:
    /*! \brief Builder class to manage object creation.
     *
     * This class is a member of gmx::Mdrunner so that it can initialize
     * private members of gmx::Mdrunner.
     *
     * It is non-trivial to establish an initialized gmx::Mdrunner invariant,
     * so objects can be obtained by clients using a Builder or a move.
     * Clients cannot default initialize or copy gmx::Mdrunner.
     */
    class BuilderImplementation;

    ~Mdrunner();

    /*!
     * \brief Copy not allowed.
     *
     * An Mdrunner has unique resources and it is not clear whether any of
     * one of those resources should be duplicated or shared unless the
     * specific use case is known. Either build a fresh runner or use a
     * helper function for clearly indicated behavior. API clarification may
     * allow unambiguous initialization by copy in future versions.
     *
     * \{
     */
    Mdrunner(const Mdrunner&) = delete;
    Mdrunner& operator=(const Mdrunner&) = delete;
    /* \} */

    /*!
     * \brief Mdrunner objects can be passed by value via move semantics.
     *
     * \param handle runner instance to be moved from.
     * \{
     */
    Mdrunner(Mdrunner&& handle) noexcept;
    Mdrunner& operator=(Mdrunner&& handle) noexcept;
    /* \} */

    /*! \brief Driver routine, that calls the different simulation methods. */
    /*!
     * Currently, thread-MPI does not spawn threads until during mdrunner() and parallelism
     * is not initialized until some time during this call...
     */
    int mdrunner();

    /*!
     * \brief Add a potential to be evaluated during MD integration.
     *
     * \param restraint MD restraint potential to apply
     * \param name User-friendly plain-text name to uniquely identify the puller
     *
     * This implementation attaches an object providing the gmx::IRestraintPotential
     * interface.
     * \todo Mdrunner should fetch such resources from the SimulationContext
     * rather than offering this public interface.
     */
    void addPotential(std::shared_ptr<IRestraintPotential> restraint, const std::string& name);

    /*! \brief Prepare the thread-MPI communicator to have \c
     * numThreadsToLaunch ranks, by spawning new thread-MPI
     * threads.
     *
     * Called by mdrunner() to start a specific number of threads
     * (including the main thread) for thread-parallel runs. This
     * in turn calls mdrunner() for each thread. */
    void spawnThreads(int numThreadsToLaunch);

    /*! \brief Initializes a new Mdrunner from the main.
     *
     * Run this method in a new thread from a main runner to get additional
     * workers on spawned threads.
     *
     * \returns New Mdrunner instance suitable for thread-MPI work on new ranks.
     *
     * \internal
     * \todo clarify (multiple) invariants during MD runner start-up.
     * The runner state before and after launching threads is distinct enough that
     * it should be codified in the invariants of different classes. That would
     * mean that the object returned by this method would be of a different type
     * than the object held by the client up to the point of call, and its name
     * would be changed to "launchOnSpawnedThread" or something not including the
     * word "clone".
     */
    Mdrunner cloneOnSpawnedThread() const;

private:
    /*! \brief Constructor. */
    explicit Mdrunner(std::unique_ptr<MDModules> mdModules);

    //! Parallelism-related user options.
    gmx_hw_opt_t hw_opt;

    //! Filenames and properties from command-line argument values.
    ArrayRef<const t_filenm> filenames;

    /*! \brief Output context for writing text files
     *
     * \internal
     * \todo push this data member down when the information can be queried from an encapsulated resource.
     */
    gmx_output_env_t* oenv = nullptr;
    //! Ongoing collection of mdrun options
    MdrunOptions mdrunOptions;
    //! Options for the domain decomposition.
    DomdecOptions domdecOptions;

    /*! \brief Target short-range interactions for "cpu", "gpu", or "auto". Default is "auto".
     *
     * \internal
     * \todo replace with string or enum class and initialize with sensible value.
     */
    const char* nbpu_opt = nullptr;

    /*! \brief Target long-range interactions for "cpu", "gpu", or "auto". Default is "auto".
     *
     * \internal
     * \todo replace with string or enum class and initialize with sensible value.
     */
    const char* pme_opt = nullptr;

    /*! \brief Target long-range interactions FFT/solve stages for "cpu", "gpu", or "auto". Default is "auto".
     *
     * \internal
     * \todo replace with string or enum class and initialize with sensible value.
     */
    const char* pme_fft_opt = nullptr;

    /*! \brief Target bonded interactions for "cpu", "gpu", or "auto". Default is "auto".
     *
     * \internal
     * \todo replace with string or enum class and initialize with sensible value.
     */
    const char* bonded_opt = nullptr;

    /*! \brief Target update calculation for "cpu", "gpu", or "auto". Default is "auto".
     *
     * \internal
     * \todo replace with string or enum class and initialize with sensible value.
     */
    const char* update_opt = nullptr;

    //! Command-line override for the duration of a neighbor list with the Verlet scheme.
    int nstlist_cmdline = 0;
    //! Parameters for replica-exchange simulations.
    ReplicaExchangeParameters replExParams;
    //! Print a warning if any force is larger than this (in kJ/mol nm).
    real pforce = -1;

    //! Handle to file used for logging.
    LogFilePtr logFileGuard = nullptr;
    //! \brief Non-owning handle to file used for logging.
    t_fileio* logFileHandle = nullptr;

    /*! \brief Non-owning handle to world communication data structure for task assigment.
     *
     * With real MPI, gets a value from the SimulationContext
     * supplied to the MdrunnerBuilder. With thread-MPI gets a
     * value after threads have been spawned. */
    MPI_Comm libraryWorldCommunicator = MPI_COMM_NULL;

    /*! \brief Non-owning handle to communication data structure for the current simulation.
     *
     * With real MPI, gets a value from the SimulationContext
     * supplied to the MdrunnerBuilder. With thread-MPI gets a
     * value after threads have been spawned. */
    MPI_Comm simulationCommunicator = MPI_COMM_NULL;

    //! \brief Non-owning handle to multi-simulation handler.
    gmx_multisim_t* ms = nullptr;

    //! Whether the simulation will start afresh, or restart with/without appending.
    StartingBehavior startingBehavior = StartingBehavior::NewSimulation;

    /*!
     * \brief Handle to restraints manager for the current process.
     *
     * \internal
     * Use opaque pointer for this implementation detail.
     */
    std::unique_ptr<RestraintManager> restraintManager_;

    /*!
     * \brief Builder for stop signal handler
     *
     * Optionally provided through MdrunnerBuilder. Client may create a
     * StopHandlerBuilder and register any number of signal providers before
     * launching the Mdrunner.
     *
     * Default is an empty signal handler that will have local signal issuers
     * added after being passed into the integrator.
     *
     * \internal
     * We do not need a full type specification here, so we use an opaque pointer.
     */
    std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder_;
    //! The modules that comprise mdrun.
    std::unique_ptr<MDModules> mdModules_;

    //! Non-owning handle to the results of the hardware detection.
    const gmx_hw_info_t* hwinfo_ = nullptr;

    /*!
     * \brief Holds simulation input specification provided by client, if any.
     *
     * If present on any instance (rank) of a simulation runner, an identical
     * (or compatible) SimulationInput must be held on all cooperating instances.
     */
    SimulationInputHandle inputHolder_;
};

/*! \libinternal
 * \brief Build a gmx::Mdrunner.
 *
 * Client code (such as `gmx mdrun`) uses this builder to get an initialized Mdrunner.
 *
 * A builder allows the library to ensure that client code cannot obtain an
 * uninitialized or partially initialized runner by refusing to build() if the
 * client has not provided sufficient or self-consistent direction. Director
 * code can be implemented for different user interfaces, encapsulating any
 * run-time functionality that does not belong in the library MD code, such
 * as command-line option processing or interfacing to external libraries.
 *
 * \ingroup module_mdrun
 *
 * \internal
 *
 * The initial Builder implementation is neither extensible at run time nor
 * at compile time. Future implementations should evolve to compose the runner,
 * rather than just consolidating the parameters for initialization, but there
 * is not yet a firm design for how flexibly module code will be coupled to
 * the builder and how much of the client interface will be in this Builder
 * versus Builders provided by the various modules.
 *
 * The named components for the initial builder implementation are descriptive
 * of the state of mdrun at the time, and are not intended to be prescriptive of
 * future design.
 * The probable course of GROMACS development is for the modular components that
 * support MD simulation to independently express their input parameters (required
 * and optional) and to provide some sort of help to the UI for input preparation.
 * If each module provides or aids the instantiation of a Director
 * for the client code, the Directors could be constructed with a handle to this
 * Builder and it would not need a public interface.
 *
 * As the modules are more clearly encapsulated, each module can provide its own
 * builder, user interface helpers, and/or composable Director code.
 * The runner and client code will also have to be updated as appropriate
 * default behavior is clarified for
 * (a) default behavior of client when user does not provide input,
 * (b) default behavior of builder when client does not provide input, and
 * (c) default behavior of runner when builder does not provide input.
 */
class MdrunnerBuilder final
{
public:
    /*!
     * \brief Constructor requires a handle to a SimulationContext to share.
     *
     * \param mdModules  The handle to the set of modules active in mdrun
     * \param context    Required handle to simulation context
     *
     * The calling code must guarantee that the
     * pointer remains valid for the lifetime of the builder, and that the
     * resources retrieved from the context remain valid for the lifetime of
     * the runner produced.
     */
    explicit MdrunnerBuilder(std::unique_ptr<MDModules>           mdModules,
                             compat::not_null<SimulationContext*> context);

    //! \cond
    MdrunnerBuilder()                       = delete;
    MdrunnerBuilder(const MdrunnerBuilder&) = delete;
    MdrunnerBuilder& operator=(const MdrunnerBuilder&) = delete;
    //! \endcond

    /*! \brief Allow transfer of ownership with move semantics.
     *
     * \param builder source object to transfer.
     *
     * \{
     */
    MdrunnerBuilder(MdrunnerBuilder&& builder) noexcept;
    MdrunnerBuilder& operator=(MdrunnerBuilder&& builder) noexcept;
    //! \}

    /*!
     * \brief Get ownership of an initialized gmx::Mdrunner.
     *
     * After build() is called, the Builder object should not be used
     * again. It is an error to call build without first calling all builder
     * methods described as "required."
     *
     * \return A new Mdrunner.
     *
     * \throws APIError if a required component has not been added before calling build().
     */
    Mdrunner build();

    /*!
     * \brief Supply the result of hardware detection to the gmx::Mdrunner
     *
     * \param hwinfo  Non-owning not-null handle to result of hardware detection.
     *
     * \todo It would be better to express this as either a not-null const pointer or
     * a const reference, but neither of those is consistent with incremental
     * building of an object. This motivates future work to be able to make a deep copy
     * of the detection result. See https://gitlab.com/gromacs/gromacs/-/issues/3650 */
    MdrunnerBuilder& addHardwareDetectionResult(const gmx_hw_info_t* hwinfo);

    /*!
     * \brief Set up non-bonded short-range force calculations.
     *
     * Required. Director code must provide valid options for the non-bonded
     * interaction code. The builder does not apply any defaults.
     *
     * \param nbpu_opt Target short-range interactions for "cpu", "gpu", or "auto".
     *
     * Calling must guarantee that the pointed-to C string is valid through
     * simulation launch.
     *
     * \internal
     * \todo Replace with string or enum that we can have sensible defaults for.
     * \todo Either the Builder or modular Director code should provide sensible defaults.
     */
    MdrunnerBuilder& addNonBonded(const char* nbpu_opt);

    /*!
     * \brief Set up long-range electrostatics calculations.
     *
     * Required. Director code should provide valid options for PME electrostatics,
     * whether or not PME electrostatics are used. The builder does not apply
     * any defaults, so client code should be prepared to provide (e.g.) "auto"
     * in the event no user input or logic provides an alternative argument.
     *
     * \param pme_opt Target long-range interactions for "cpu", "gpu", or "auto".
     * \param pme_fft_opt Target long-range interactions FFT/solve stages for "cpu", "gpu", or "auto".
     *
     * Calling must guarantee that the pointed-to C strings are valid through
     * simulation launch.
     *
     * \internal
     * The arguments are passed as references to elements of arrays of C strings.
     * \todo Replace with modern strings or (better) enum classes.
     * \todo Make optional and/or encapsulate into electrostatics module.
     */
    MdrunnerBuilder& addElectrostatics(const char* pme_opt, const char* pme_fft_opt);

    /*!
     * \brief Assign responsibility for tasks for bonded interactions.
     *
     * Required. Director code should provide valid options for
     * bonded interaction task assignment, whether or not such
     * interactions are present. The builder does not apply any
     * defaults, so client code should be prepared to provide
     * (e.g.) "auto" in the event no user input or logic provides
     * an alternative argument.
     *
     * \param bonded_opt Target bonded interactions for "cpu", "gpu", or "auto".
     *
     * Calling must guarantee that the pointed-to C strings are valid through
     * simulation launch.
     *
     * \internal
     * The arguments are passed as references to elements of arrays of C strings.
     * \todo Replace with modern strings or (better) enum classes.
     * \todo Make optional and/or encapsulate into task assignment module.
     */
    MdrunnerBuilder& addBondedTaskAssignment(const char* bonded_opt);

    /*! \brief
     * Assign responsibility for tasks for update and constrain calculation.
     *
     * Required. Director code should provide valid options for
     * update and constraint task assignment. The builder does not apply any
     * defaults, so client code should be prepared to provide
     * (e.g.) "auto" in the event no user input or logic provides
     * an alternative argument.
     *
     * \param[in] update_opt Target update calculation for "cpu", "gpu", or "auto".
     *
     * Calling must guarantee that the pointed-to C strings are valid through
     * simulation launch.
     *
     * \internal
     * The arguments are passed as references to elements of arrays of C strings.
     * \todo Replace with modern strings or (better) enum classes.
     * \todo Make optional and/or encapsulate into task assignment module.
     */
    MdrunnerBuilder& addUpdateTaskAssignment(const char* update_opt);

    /*!
     * \brief Set MD options not owned by some other module.
     *
     * Optional. Override simulation parameters
     *
     * \param options structure to copy
     * \param forceWarningThreshold Print a warning if any force is larger than this (in kJ/mol nm) (default -1)
     * \param startingBehavior Whether the simulation will start afresh, or restart with/without appending.
     *
     * \internal
     * \todo Map these parameters to more appropriate encapsulating types.
     * Find a better way to indicate "unspecified" than a magic value of the parameter type.
     */
    MdrunnerBuilder& addSimulationMethod(const MdrunOptions& options,
                                         real                forceWarningThreshold,
                                         StartingBehavior    startingBehavior);

    /*!
     * \brief Set the domain decomposition module.
     *
     * Optional. Overrides default constructed DomdecOptions if provided.
     *
     * \param options options with which to construct domain decomposition.
     *
     * \internal
     * \todo revisit whether we should be passing this parameter struct or a higher-level handle of some sort.
     */
    MdrunnerBuilder& addDomainDecomposition(const DomdecOptions& options);

    /*!
     * \brief Set Verlet list manager.
     *
     * Optional. Neighbor list existence, type, and parameters are mostly determined
     * by the simulation parameters loaded elsewhere. This is just an override.
     *
     * \param rebuildInterval override for the duration of a neighbor list with the Verlet scheme.
     */
    MdrunnerBuilder& addNeighborList(int rebuildInterval);

    /*!
     * \brief Set replica exchange manager.
     *
     * Optional. For guidance on preparing a valid ReplicaExchangeParameters
     * value, refer to the details in mdrun.cpp, the `t_pargs pa[]` defined there,
     * and the action of parse_common_args() with regards to that structure.
     * If not provided by client, a default constructed ReplicaExchangeParameters
     * is used.
     *
     * \param params parameters with which to set up replica exchange.
     *
     * \internal
     * \todo revisit whether we should be passing this parameter struct or a higher-level handle of some sort.
     */
    MdrunnerBuilder& addReplicaExchange(const ReplicaExchangeParameters& params);

    /*!
     * \brief Specify parameters determining hardware resource allocation.
     *
     * Optional. If not provided, default-constructed gmx_hw_opt_t will be used.
     *
     * \param hardwareOptions Parallelism-related user options.
     */
    MdrunnerBuilder& addHardwareOptions(const gmx_hw_opt_t& hardwareOptions);

    /*!
     * \brief Provide the filenames options structure with option values chosen
     *
     * Required. The object is assumed to have been updated by
     * parse_common_args or equivalent.
     *
     * \param filenames Filenames and properties from command-line argument values or defaults.
     *
     * \internal
     * \todo Modules should manage their own filename options and defaults.
     */
    MdrunnerBuilder& addFilenames(ArrayRef<const t_filenm> filenames);

    /*!
     * \brief Provide parameters for setting up output environment.
     *
     * Required. Handle is assumed to have been produced by output_env_init
     * as in parse_common_args.
     *
     * \param outputEnvironment Output context for writing text files.
     *
     * \internal
     * \todo Allow client code to set up output environment and provide as a resource.
     * This parameter is used to set up resources that are dependent on the execution
     * environment and API context. Such resources should be retrieved by the simulator
     * from a client-provided resource, but currently the resources are only fully
     * initialized in Mdrunner.
     */
    MdrunnerBuilder& addOutputEnvironment(gmx_output_env_t* outputEnvironment);

    /*!
     * \brief Provide the filehandle pointer to be used for the MD log.
     *
     * Required. Either nullptr if no log should be written, or
     * valid and open reading for writing.
     *
     * \param logFileHandle Non-owning handle to file used for logging.
     * \internal
     */
    MdrunnerBuilder& addLogFile(t_fileio* logFileHandle);

    /*!
     * \brief Provide a StopHandlerBuilder for the MD stop signal handling.
     *
     * Optional. Defaults to empty.
     *
     * Client may provide additional (non-default) issuers of simulation stop
     * signals by preconfiguring the StopHandlerBuilder used later when the
     * simulation runs.
     *
     * \param builder
     */
    MdrunnerBuilder& addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder);

    /*!
     * \brief Acquire a handle to the SimulationInput.
     *
     * Required. SimulationInput will be taking responsibility for some of the
     * input provided through other methods, such as addFilenames.
     *
     * See also issue https://gitlab.com/gromacs/gromacs/-/issues/3374
     *
     * \param input Shared ownership of a SimulationInput.
     */
    MdrunnerBuilder& addInput(SimulationInputHandle input);

    ~MdrunnerBuilder();

private:
    std::unique_ptr<Mdrunner::BuilderImplementation> impl_;
};

} // namespace gmx

#endif // GMX_MDRUN_RUNNER_H
