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
#include "gromacs/compat/make_unique.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/logging.h"
#include "gromacs/mdrun/mdfilenames.h"
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

#include "context-impl.h"
#include "createsession.h"
#include "session-impl.h"
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

///////////////////////////////////////////////////////////
//
// START: front-end options processing copied from CLI code
//
// As this section is encapsulated, updates can be ported to
// the CLI front-end. This should happen as soon as possible,
// because this divergent code will be very hard to maintain.
// Note that modernization of the CLI code options handling
// will also allow updates in the code here.
//
// Reference gmxapi milestones 24 and 28, described at
// https://redmine.gromacs.org/issues/2585
//

        // Ongoing collection of mdrun options
        MdrunOptions                     mdrunOptions;
        // Options for the domain decomposition.
        DomdecOptions                    domdecOptions;
        // Parallelism-related user options.
        gmx_hw_opt_t                     hw_opt;
        // Command-line override for the duration of a neighbor list with the Verlet scheme.
        int                              nstlist_cmdline = 0;
        // Parameters for replica-exchange simulations.
        ReplicaExchangeParameters        replExParams;

        // Filenames and properties from command-line argument values.
        mdFilenames_ = gmx::MdFilenames();
        auto &filenames = mdFilenames_;

        // Print a warning if any force is larger than this (in kJ/mol nm).
        real                             pforce = -1;

        // Output context for writing text files
        // \todo Clarify initialization, ownership, and lifetime.
        gmx_output_env_t                *oenv = nullptr;

        /* Command line options */
        rvec              realddxyz                                               = {0, 0, 0};
        const char       *ddrank_opt_choices[static_cast<int>(DdRankOrder::nr)+1] =
        { nullptr, "interleave", "pp_pme", "cartesian", nullptr };
        const char       *dddlb_opt_choices[static_cast<int>(DlbOption::nr)+1] =
        { nullptr, "auto", "no", "yes", nullptr };
        const char       *thread_aff_opt_choices[threadaffNR+1] =
        { nullptr, "auto", "on", "off", nullptr };
        const char       *nbpu_opt_choices[] =
        { nullptr, "auto", "cpu", "gpu", nullptr };
        const char       *pme_opt_choices[] =
        { nullptr, "auto", "cpu", "gpu", nullptr };
        const char       *pme_fft_opt_choices[] =
        { nullptr, "auto", "cpu", "gpu", nullptr };
        gmx_bool          bTryToAppendFiles     = TRUE;
        const char       *gpuIdsAvailable       = "";
        const char       *userGpuTaskAssignment = "";

        ImdOptions       &imdOptions = mdrunOptions.imdOptions;

        t_pargs           pa[] = {

            { "-dd",      FALSE, etRVEC, {&realddxyz},
              "Domain decomposition grid, 0 is optimize" },
            { "-ddorder", FALSE, etENUM, {ddrank_opt_choices},
              "DD rank order" },
            { "-npme",    FALSE, etINT, {&domdecOptions.numPmeRanks},
              "Number of separate ranks to be used for PME, -1 is guess" },
            { "-nt",      FALSE, etINT, {&hw_opt.nthreads_tot},
              "Total number of threads to start (0 is guess)" },
            { "-ntmpi",   FALSE, etINT, {&hw_opt.nthreads_tmpi},
              "Number of thread-MPI ranks to start (0 is guess)" },
            { "-ntomp",   FALSE, etINT, {&hw_opt.nthreads_omp},
              "Number of OpenMP threads per MPI rank to start (0 is guess)" },
            { "-ntomp_pme", FALSE, etINT, {&hw_opt.nthreads_omp_pme},
              "Number of OpenMP threads per MPI rank to start (0 is -ntomp)" },
            { "-pin",     FALSE, etENUM, {thread_aff_opt_choices},
              "Whether mdrun should try to set thread affinities" },
            { "-pinoffset", FALSE, etINT, {&hw_opt.core_pinning_offset},
              "The lowest logical core number to which mdrun should pin the first thread" },
            { "-pinstride", FALSE, etINT, {&hw_opt.core_pinning_stride},
              "Pinning distance in logical cores for threads, use 0 to minimize the number of threads per physical core" },
            { "-gpu_id",  FALSE, etSTR, {&gpuIdsAvailable},
              "List of unique GPU device IDs available to use" },
            { "-gputasks",  FALSE, etSTR, {&userGpuTaskAssignment},
              "List of GPU device IDs, mapping each PP task on each node to a device" },
            { "-ddcheck", FALSE, etBOOL, {&domdecOptions.checkBondedInteractions},
              "Check for all bonded interactions with DD" },
            { "-ddbondcomm", FALSE, etBOOL, {&domdecOptions.useBondedCommunication},
              "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
            { "-rdd",     FALSE, etREAL, {&domdecOptions.minimumCommunicationRange},
              "The maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
            { "-rcon",    FALSE, etREAL, {&domdecOptions.constraintCommunicationRange},
              "Maximum distance for P-LINCS (nm), 0 is estimate" },
            { "-dlb",     FALSE, etENUM, {dddlb_opt_choices},
              "Dynamic load balancing (with DD)" },
            { "-dds",     FALSE, etREAL, {&domdecOptions.dlbScaling},
              "Fraction in (0,1) by whose reciprocal the initial DD cell size will be increased in order to "
              "provide a margin in which dynamic load balancing can act while preserving the minimum cell size." },
            { "-ddcsx",   FALSE, etSTR, {&domdecOptions.cellSizeX},
              "HIDDENA string containing a vector of the relative sizes in the x "
              "direction of the corresponding DD cells. Only effective with static "
              "load balancing." },
            { "-ddcsy",   FALSE, etSTR, {&domdecOptions.cellSizeY},
              "HIDDENA string containing a vector of the relative sizes in the y "
              "direction of the corresponding DD cells. Only effective with static "
              "load balancing." },
            { "-ddcsz",   FALSE, etSTR, {&domdecOptions.cellSizeZ},
              "HIDDENA string containing a vector of the relative sizes in the z "
              "direction of the corresponding DD cells. Only effective with static "
              "load balancing." },
            { "-gcom",    FALSE, etINT, {&mdrunOptions.globalCommunicationInterval},
              "Global communication frequency" },
            { "-nb",      FALSE, etENUM, {nbpu_opt_choices},
              "Calculate non-bonded interactions on" },
            { "-nstlist", FALSE, etINT, {&nstlist_cmdline},
              "Set nstlist when using a Verlet buffer tolerance (0 is guess)" },
            { "-tunepme", FALSE, etBOOL, {&mdrunOptions.tunePme},
              "Optimize PME load between PP/PME ranks or GPU/CPU (only with the Verlet cut-off scheme)" },
            { "-pme",     FALSE, etENUM, {pme_opt_choices},
              "Perform PME calculations on" },
            { "-pmefft", FALSE, etENUM, {pme_fft_opt_choices},
              "Perform PME FFT calculations on" },
            { "-v",       FALSE, etBOOL, {&mdrunOptions.verbose},
              "Be loud and noisy" },
            { "-pforce",  FALSE, etREAL, {&pforce},
              "Print all forces larger than this (kJ/mol nm)" },
            { "-reprod",  FALSE, etBOOL, {&mdrunOptions.reproducible},
              "Try to avoid optimizations that affect binary reproducibility" },
            { "-cpt",     FALSE, etREAL, {&mdrunOptions.checkpointOptions.period},
              "Checkpoint interval (minutes)" },
            { "-cpnum",   FALSE, etBOOL, {&mdrunOptions.checkpointOptions.keepAndNumberCheckpointFiles},
              "Keep and number checkpoint files" },
            { "-append",  FALSE, etBOOL, {&bTryToAppendFiles},
              "Append to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
            { "-nsteps",  FALSE, etINT64, {&mdrunOptions.numStepsCommandline},
              "Run this number of steps, overrides .mdp file option (-1 means infinite, -2 means use mdp option, smaller is invalid)" },
            { "-maxh",   FALSE, etREAL, {&mdrunOptions.maximumHoursToRun},
              "Terminate after 0.99 times this time (hours)" },
            { "-replex",  FALSE, etINT, {&replExParams.exchangeInterval},
              "Attempt replica exchange periodically with this period (steps)" },
            { "-nex",  FALSE, etINT, {&replExParams.numExchanges},
              "Number of random exchanges to carry out each exchange interval (N^3 is one suggestion).  -nex zero or not specified gives neighbor replica exchange." },
            { "-reseed",  FALSE, etINT, {&replExParams.randomSeed},
              "Seed for replica exchange, -1 is generate a seed" },
            { "-imdport",    FALSE, etINT, {&imdOptions.port},
              "HIDDENIMD listening port" },
            { "-imdwait",  FALSE, etBOOL, {&imdOptions.wait},
              "HIDDENPause the simulation while no IMD client is connected" },
            { "-imdterm",  FALSE, etBOOL, {&imdOptions.terminatable},
              "HIDDENAllow termination of the simulation from IMD client" },
            { "-imdpull",  FALSE, etBOOL, {&imdOptions.pull},
              "HIDDENAllow pulling in the simulation from IMD client" },
            { "-rerunvsite", FALSE, etBOOL, {&mdrunOptions.rerunConstructVsites},
              "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
            { "-confout", FALSE, etBOOL, {&mdrunOptions.writeConfout},
              "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
            { "-stepout", FALSE, etINT, {&mdrunOptions.verboseStepPrintInterval},
              "HIDDENFrequency of writing the remaining wall clock time for the run" },
            { "-resetstep", FALSE, etINT, {&mdrunOptions.timingOptions.resetStep},
              "HIDDENReset cycle counters after these many time steps" },
            { "-resethway", FALSE, etBOOL, {&mdrunOptions.timingOptions.resetHalfway},
              "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" }
        };

        // pointer-to-t_commrec is the de facto handle type for communications record.
        // Complete shared / borrowed ownership requires a reference to this stack variable
        // (or pointer-to-pointer-to-t_commrec) since borrowing code may update the pointer.
        // \todo Define the ownership and lifetime management semantics for a communication record, handle or value type.

        // Note: this communications record initialization acts directly on
        // MPI_COMM_WORLD and is incompatible with MPI environment sharing in
        // gmxapi through 0.0.7, at least.
        t_commrec       * cr = init_commrec();

        unsigned long     PCA_Flags = PCA_CAN_SET_DEFFNM;
        // With -multidir, the working directory still needs to be
        // changed, so we can't check for the existence of files during
        // parsing.  It isn't useful to do any completion based on file
        // system contents, either.
        bool is_multisim_option_set = false;
        for (auto i = 0; i < static_cast<decltype(i)>(argc); ++i)
        {
            if (strcmp(argv[i], "-multidir") == 0)
            {
                is_multisim_option_set = true;
            }
        }
        if (is_multisim_option_set)
        {
            PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
        }

        const char       *desc[]  = {"gmxapi placeholder text"};
        int               argcInt = static_cast<int>(argc);
        if (!parse_common_args(&argcInt, argv.data(), PCA_Flags,
                               static_cast<int>(filenames().size()), filenames().data(), asize(pa), pa,
                               asize(desc), desc, 0, nullptr, &oenv))
        {
            sfree(cr);
            return nullptr;
        }

        // Handle the options that permits the user to either declare
        // which compatible GPUs are availble for use, or to select a GPU
        // task assignment. Either could be in an environment variable (so
        // that there is a way to customize it, when using MPI in
        // heterogeneous contexts).
        {
            // TODO Argument parsing can't handle std::string. We should
            // fix that by changing the parsing, once more of the roles of
            // handling, validating and implementing defaults for user
            // command-line options have been seperated.
            hw_opt.gpuIdsAvailable       = gpuIdsAvailable;
            hw_opt.userGpuTaskAssignment = userGpuTaskAssignment;

            const char *env = getenv("GMX_GPU_ID");
            if (env != nullptr)
            {
                if (!hw_opt.gpuIdsAvailable.empty())
                {
                    gmx_fatal(FARGS, "GMX_GPU_ID and -gpu_id can not be used at the same time");
                }
                hw_opt.gpuIdsAvailable = env;
            }

            env = getenv("GMX_GPUTASKS");
            if (env != nullptr)
            {
                if (!hw_opt.userGpuTaskAssignment.empty())
                {
                    gmx_fatal(FARGS, "GMX_GPUTASKS and -gputasks can not be used at the same time");
                }
                hw_opt.userGpuTaskAssignment = env;
            }

            if (!hw_opt.gpuIdsAvailable.empty() && !hw_opt.userGpuTaskAssignment.empty())
            {
                gmx_fatal(FARGS, "-gpu_id and -gputasks cannot be used at the same time");
            }
        }

        hw_opt.thread_affinity = nenum(thread_aff_opt_choices);

        // now check for a multi-simulation
        gmx::ArrayRef<const std::string> multidir = opt2fnsIfOptionSet("-multidir",
                                                                       static_cast<int>(filenames().size()),
                                                                       filenames().data());

        if (replExParams.exchangeInterval != 0 && multidir.size() < 2)
        {
            gmx_fatal(FARGS, "Need at least two replicas for replica exchange (use option -multidir)");
        }

        if (replExParams.numExchanges < 0)
        {
            gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
        }

        gmx_multisim_t* ms = init_multisystem(MPI_COMM_WORLD, multidir);

        /* Prepare the intra-simulation communication */
        // TODO consolidate this with init_commrec, after changing the
        // relative ordering of init_commrec and init_multisystem
#if GMX_MPI
        if (ms != nullptr)
        {
            cr->nnodes = cr->nnodes / ms->nsim;
            MPI_Comm_split(MPI_COMM_WORLD, ms->sim, cr->sim_nodeid, &cr->mpi_comm_mysim);
            cr->mpi_comm_mygroup = cr->mpi_comm_mysim;
            MPI_Comm_rank(cr->mpi_comm_mysim, &cr->sim_nodeid);
            MPI_Comm_rank(cr->mpi_comm_mygroup, &cr->nodeid);
        }
#endif

        if (!opt2bSet("-cpi",
                      static_cast<int>(filenames().size()), filenames().data()))
        {
            // If we are not starting from a checkpoint we never allow files to be appended
            // to, since that has caused a ton of strange behaviour and bugs in the past.
            if (opt2parg_bSet("-append", asize(pa), pa))
            {
                // If the user explicitly used the -append option, explain that it is not possible.
                gmx_fatal(FARGS, "GROMACS can only append to files when restarting from a checkpoint.");
            }
            else
            {
                // If the user did not say anything explicit, just disable appending.
                bTryToAppendFiles = FALSE;
            }
        }

        ContinuationOptions &continuationOptions = mdrunOptions.continuationOptions;

        continuationOptions.appendFilesOptionSet = opt2parg_bSet("-append", asize(pa), pa);

        handleRestart(cr, ms, bTryToAppendFiles,
                      static_cast<int>(filenames().size()),
                      filenames().data(),
                      &continuationOptions.appendFiles,
                      &continuationOptions.startedFromCheckpoint);

        mdrunOptions.rerun            = opt2bSet("-rerun",
                                                 static_cast<int>(filenames().size()),
                                                 filenames().data());
        mdrunOptions.ntompOptionIsSet = opt2parg_bSet("-ntomp", asize(pa), pa);

        LogFilePtr logFileGuard(nullptr);
        if (MASTER(cr))
        {
            logFileGuard = openLogFile(ftp2fn(efLOG,
                                              filenames().size(),
                                              filenames().data()),
                                       continuationOptions.appendFiles,
                                       cr->nodeid,
                                       cr->nnodes);
        }

        domdecOptions.rankOrder    = static_cast<DdRankOrder>(nenum(ddrank_opt_choices));
        domdecOptions.dlbOption    = static_cast<DlbOption>(nenum(dddlb_opt_choices));
        domdecOptions.numCells[XX] = roundToInt(realddxyz[XX]);
        domdecOptions.numCells[YY] = roundToInt(realddxyz[YY]);
        domdecOptions.numCells[ZZ] = roundToInt(realddxyz[ZZ]);

//
// END: front-end options processing copied from CLI code
/////////////////////////////////////////////////////////
        auto simulationContext = createSimulationContext(cr);

        auto builder = MdrunnerBuilder(compat::not_null<decltype( &simulationContext)>(&simulationContext));
        builder.addSimulationMethod(mdrunOptions, pforce);
        builder.addDomainDecomposition(domdecOptions);
        // \todo pass by value
        builder.addNonBonded(nbpu_opt_choices[0]);
        // \todo pass by value
        builder.addElectrostatics(pme_opt_choices[0], pme_fft_opt_choices[0]);
        builder.addNeighborList(nstlist_cmdline);
        builder.addReplicaExchange(replExParams);
        // \todo take ownership of multisim resources (ms)
        builder.addMultiSim(ms);
        // \todo Provide parallelism resources through SimulationContext.
        // Need to establish run-time values from various inputs to provide a resource handle to Mdrunner
        builder.addHardwareOptions(hw_opt);
        // \todo File names are parameters that should be managed modularly through further factoring.
        builder.addFilenames(filenames());
        // Note: The gmx_output_env_t life time is not managed after the call to parse_common_args.
        // \todo Implement lifetime management for gmx_output_env_t.
        // \todo Output environment should be configured outside of Mdrunner and provided as a resource.
        builder.addOutputEnvironment(oenv);
        builder.addLogFile(logFileGuard.get());

        // Note, creation is not mature enough to be exposed in the external API yet.
        launchedSession = createSession(shared_from_this(),
                                        std::move(builder),
                                        simulationContext,
                                        std::move(logFileGuard),
                                        ms);

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
