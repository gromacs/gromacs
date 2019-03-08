/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Implements the MD runner routine calling all integrators.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "runner.h"

#include "config.h"

#include <cassert>
#include <cinttypes>
#include <csignal>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/ewald/ewald-utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-gpu-program.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/clfftinitializer.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/printhardware.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/gpubonded.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_tuning.h"
#include "gromacs/mdlib/ppforceworkload.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdrun/legacymdrunoptions.h"
#include "gromacs/mdrun/logging.h"
#include "gromacs/mdrun/multisim.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/restraint/manager.h"
#include "gromacs/restraint/restraintmdmodule.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/taskassignment/decidegpuusage.h"
#include "gromacs/taskassignment/resourcedivision.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "integrator.h"
#include "replicaexchange.h"

#if GMX_FAHCORE
#include "corewrap.h"
#endif

namespace gmx
{

/*! \brief Barrier for safe simultaneous thread access to mdrunner data
 *
 * Used to ensure that the master thread does not modify mdrunner during copy
 * on the spawned threads. */
static void threadMpiMdrunnerAccessBarrier()
{
#if GMX_THREAD_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

Mdrunner Mdrunner::cloneOnSpawnedThread() const
{
    auto newRunner = Mdrunner();

    // All runners in the same process share a restraint manager resource because it is
    // part of the interface to the client code, which is associated only with the
    // original thread. Handles to the same resources can be obtained by copy.
    {
        newRunner.restraintManager_ = compat::make_unique<RestraintManager>(*restraintManager_);
    }

    // Copy original cr pointer before master thread can pass the thread barrier
    newRunner.cr  = reinitialize_commrec_for_this_thread(cr);

    // Copy members of master runner.
    // \todo Replace with builder when Simulation context and/or runner phases are better defined.
    // Ref https://redmine.gromacs.org/issues/2587 and https://redmine.gromacs.org/issues/2375
    newRunner.hw_opt    = hw_opt;
    newRunner.filenames = filenames;

    newRunner.oenv                = oenv;
    newRunner.mdrunOptions        = mdrunOptions;
    newRunner.domdecOptions       = domdecOptions;
    newRunner.nbpu_opt            = nbpu_opt;
    newRunner.pme_opt             = pme_opt;
    newRunner.pme_fft_opt         = pme_fft_opt;
    newRunner.bonded_opt          = bonded_opt;
    newRunner.nstlist_cmdline     = nstlist_cmdline;
    newRunner.replExParams        = replExParams;
    newRunner.pforce              = pforce;
    newRunner.ms                  = ms;
    newRunner.stopHandlerBuilder_ = compat::make_unique<StopHandlerBuilder>(*stopHandlerBuilder_);

    threadMpiMdrunnerAccessBarrier();

    GMX_RELEASE_ASSERT(!MASTER(newRunner.cr), "cloneOnSpawnedThread should only be called on spawned threads");

    return newRunner;
}

/*! \brief The callback used for running on spawned threads.
 *
 * Obtains the pointer to the master mdrunner object from the one
 * argument permitted to the thread-launch API call, copies it to make
 * a new runner for this thread, reinitializes necessary data, and
 * proceeds to the simulation. */
static void mdrunner_start_fn(const void *arg)
{
    try
    {
        auto masterMdrunner = reinterpret_cast<const gmx::Mdrunner *>(arg);
        /* copy the arg list to make sure that it's thread-local. This
           doesn't copy pointed-to items, of course; fnm, cr and fplog
           are reset in the call below, all others should be const. */
        gmx::Mdrunner mdrunner = masterMdrunner->cloneOnSpawnedThread();
        mdrunner.mdrunner();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}


/*! \brief Start thread-MPI threads.
 *
 * Called by mdrunner() to start a specific number of threads
 * (including the main thread) for thread-parallel runs. This in turn
 * calls mdrunner() for each thread. All options are the same as for
 * mdrunner(). */
t_commrec *Mdrunner::spawnThreads(int numThreadsToLaunch) const
{

    /* first check whether we even need to start tMPI */
    if (numThreadsToLaunch < 2)
    {
        return cr;
    }

#if GMX_THREAD_MPI
    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns, we set thread affinity later */
    if (tMPI_Init_fn(TRUE, numThreadsToLaunch, TMPI_AFFINITY_NONE,
                     mdrunner_start_fn, static_cast<const void*>(this)) != TMPI_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to spawn thread-MPI threads"));
    }

    threadMpiMdrunnerAccessBarrier();
#else
    GMX_UNUSED_VALUE(mdrunner_start_fn);
#endif

    return reinitialize_commrec_for_this_thread(cr);
}

}  // namespace gmx

/*! \brief Initialize variables for Verlet scheme simulation */
static void prepare_verlet_scheme(FILE                           *fplog,
                                  t_commrec                      *cr,
                                  t_inputrec                     *ir,
                                  int                             nstlist_cmdline,
                                  const gmx_mtop_t               *mtop,
                                  const matrix                    box,
                                  bool                            makeGpuPairList,
                                  const gmx::CpuInfo             &cpuinfo)
{
    /* For NVE simulations, we will retain the initial list buffer */
    if (EI_DYNAMICS(ir->eI) &&
        ir->verletbuf_tol > 0 &&
        !(EI_MD(ir->eI) && ir->etc == etcNO))
    {
        /* Update the Verlet buffer size for the current run setup */

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        ListSetupType      listType  = (makeGpuPairList ? ListSetupType::Gpu : ListSetupType::CpuSimdWhenSupported);
        VerletbufListSetup listSetup = verletbufGetSafeListSetup(listType);

        real               rlist_new;
        calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1, -1, &listSetup, nullptr, &rlist_new);

        if (rlist_new != ir->rlist)
        {
            if (fplog != nullptr)
            {
                fprintf(fplog, "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist, rlist_new,
                        listSetup.cluster_size_i, listSetup.cluster_size_j);
            }
            ir->rlist     = rlist_new;
        }
    }

    if (nstlist_cmdline > 0 && (!EI_DYNAMICS(ir->eI) || ir->verletbuf_tol <= 0))
    {
        gmx_fatal(FARGS, "Can not set nstlist without %s",
                  !EI_DYNAMICS(ir->eI) ? "dynamics" : "verlet-buffer-tolerance");
    }

    if (EI_DYNAMICS(ir->eI))
    {
        /* Set or try nstlist values */
        increaseNstlist(fplog, cr, ir, nstlist_cmdline, mtop, box, makeGpuPairList, cpuinfo);
    }
}

/*! \brief Override the nslist value in inputrec
 *
 * with value passed on the command line (if any)
 */
static void override_nsteps_cmdline(const gmx::MDLogger &mdlog,
                                    int64_t              nsteps_cmdline,
                                    t_inputrec          *ir)
{
    assert(ir);

    /* override with anything else than the default -2 */
    if (nsteps_cmdline > -2)
    {
        char sbuf_steps[STEPSTRSIZE];
        char sbuf_msg[STRLEN];

        ir->nsteps = nsteps_cmdline;
        if (EI_DYNAMICS(ir->eI) && nsteps_cmdline != -1)
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps, %.3g ps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps),
                    fabs(nsteps_cmdline*ir->delta_t));
        }
        else
        {
            sprintf(sbuf_msg, "Overriding nsteps with value passed on the command line: %s steps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps));
        }

        GMX_LOG(mdlog.warning).asParagraph().appendText(sbuf_msg);
    }
    else if (nsteps_cmdline < -2)
    {
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %" PRId64,
                  nsteps_cmdline);
    }
    /* Do nothing if nsteps_cmdline == -2 */
}

namespace gmx
{

/*! \brief Return whether GPU acceleration of nonbondeds is supported with the given settings.
 *
 * If not, and if a warning may be issued, logs a warning about
 * falling back to CPU code. With thread-MPI, only the first
 * call to this function should have \c issueWarning true. */
static bool gpuAccelerationOfNonbondedIsUseful(const MDLogger   &mdlog,
                                               const t_inputrec *ir,
                                               bool              issueWarning)
{
    if (ir->opts.ngener - ir->nwall > 1)
    {
        /* The GPU code does not support more than one energy group.
         * If the user requested GPUs explicitly, a fatal error is given later.
         */
        if (issueWarning)
        {
            GMX_LOG(mdlog.warning).asParagraph()
                .appendText("Multiple energy groups is not implemented for GPUs, falling back to the CPU. "
                            "For better performance, run on the GPU without energy groups and then do "
                            "gmx mdrun -rerun option on the trajectory with an energy group .tpr file.");
        }
        return false;
    }
    return true;
}

//! Initializes the logger for mdrun.
static gmx::LoggerOwner buildLogger(FILE *fplog, const t_commrec *cr)
{
    gmx::LoggerBuilder builder;
    if (fplog != nullptr)
    {
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fplog);
    }
    if (cr == nullptr || SIMMASTER(cr))
    {
        builder.addTargetStream(gmx::MDLogger::LogLevel::Warning,
                                &gmx::TextOutputFile::standardError());
    }
    return builder.build();
}

//! Make a TaskTarget from an mdrun argument string.
static TaskTarget findTaskTarget(const char *optionString)
{
    TaskTarget returnValue = TaskTarget::Auto;

    if (strncmp(optionString, "auto", 3) == 0)
    {
        returnValue = TaskTarget::Auto;
    }
    else if (strncmp(optionString, "cpu", 3) == 0)
    {
        returnValue = TaskTarget::Cpu;
    }
    else if (strncmp(optionString, "gpu", 3) == 0)
    {
        returnValue = TaskTarget::Gpu;
    }
    else
    {
        GMX_ASSERT(false, "Option string should have been checked for sanity already");
    }

    return returnValue;
}

int Mdrunner::mdrunner()
{
    matrix                    box;
    t_nrnb                   *nrnb;
    t_forcerec               *fr               = nullptr;
    t_fcdata                 *fcd              = nullptr;
    real                      ewaldcoeff_q     = 0;
    real                      ewaldcoeff_lj    = 0;
    int                       nChargePerturbed = -1, nTypePerturbed = 0;
    gmx_wallcycle_t           wcycle;
    gmx_walltime_accounting_t walltime_accounting = nullptr;
    int                       rc;
    int64_t                   reset_counters;
    int                       nthreads_pme = 1;
    gmx_membed_t *            membed       = nullptr;
    gmx_hw_info_t            *hwinfo       = nullptr;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    std::unique_ptr<gmx::MDModules> mdModules(new gmx::MDModules);
    t_inputrec                      inputrecInstance;
    t_inputrec                     *inputrec = &inputrecInstance;
    gmx_mtop_t                      mtop;

    bool doMembed = opt2bSet("-membed", filenames.size(), filenames.data());
    bool doRerun  = mdrunOptions.rerun;

    // Handle task-assignment related user options.
    EmulateGpuNonbonded emulateGpuNonbonded = (getenv("GMX_EMULATE_GPU") != nullptr ?
                                               EmulateGpuNonbonded::Yes : EmulateGpuNonbonded::No);
    std::vector<int>    gpuIdsAvailable;
    try
    {
        gpuIdsAvailable = parseUserGpuIds(hw_opt.gpuIdsAvailable);
        // TODO We could put the GPU IDs into a std::map to find
        // duplicates, but for the small numbers of IDs involved, this
        // code is simple and fast.
        for (size_t i = 0; i != gpuIdsAvailable.size(); ++i)
        {
            for (size_t j = i+1; j != gpuIdsAvailable.size(); ++j)
            {
                if (gpuIdsAvailable[i] == gpuIdsAvailable[j])
                {
                    GMX_THROW(InvalidInputError(formatString("The string of available GPU device IDs '%s' may not contain duplicate device IDs", hw_opt.gpuIdsAvailable.c_str())));
                }
            }
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    std::vector<int> userGpuTaskAssignment;
    try
    {
        userGpuTaskAssignment = parseUserGpuIds(hw_opt.userGpuTaskAssignment);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    auto       nonbondedTarget = findTaskTarget(nbpu_opt);
    auto       pmeTarget       = findTaskTarget(pme_opt);
    auto       pmeFftTarget    = findTaskTarget(pme_fft_opt);
    auto       bondedTarget    = findTaskTarget(bonded_opt);
    PmeRunMode pmeRunMode      = PmeRunMode::None;

    // Here we assume that SIMMASTER(cr) does not change even after the
    // threads are started.

    FILE *fplog = nullptr;
    // If we are appending, we don't write log output because we need
    // to check that the old log file matches what the checkpoint file
    // expects. Otherwise, we should start to write log output now if
    // there is a file ready for it.
    if (logFileHandle != nullptr && !mdrunOptions.continuationOptions.appendFiles)
    {
        fplog = gmx_fio_getfp(logFileHandle);
    }
    gmx::LoggerOwner logOwner(buildLogger(fplog, cr));
    gmx::MDLogger    mdlog(logOwner.logger());

    // TODO The thread-MPI master rank makes a working
    // PhysicalNodeCommunicator here, but it gets rebuilt by all ranks
    // after the threads have been launched. This works because no use
    // is made of that communicator until after the execution paths
    // have rejoined. But it is likely that we can improve the way
    // this is expressed, e.g. by expressly running detection only the
    // master rank for thread-MPI, rather than relying on the mutex
    // and reference count.
    PhysicalNodeCommunicator physicalNodeComm(MPI_COMM_WORLD, gmx_physicalnode_id_hash());
    hwinfo = gmx_detect_hardware(mdlog, physicalNodeComm);

    gmx_print_detected_hardware(fplog, cr, ms, mdlog, hwinfo);

    std::vector<int> gpuIdsToUse;
    auto             compatibleGpus = getCompatibleGpus(hwinfo->gpu_info);
    if (gpuIdsAvailable.empty())
    {
        gpuIdsToUse = compatibleGpus;
    }
    else
    {
        for (const auto &availableGpuId : gpuIdsAvailable)
        {
            bool availableGpuIsCompatible = false;
            for (const auto &compatibleGpuId : compatibleGpus)
            {
                if (availableGpuId == compatibleGpuId)
                {
                    availableGpuIsCompatible = true;
                    break;
                }
            }
            if (!availableGpuIsCompatible)
            {
                gmx_fatal(FARGS, "You limited the set of compatible GPUs to a set that included ID #%d, but that ID is not for a compatible GPU. List only compatible GPUs.", availableGpuId);
            }
            gpuIdsToUse.push_back(availableGpuId);
        }
    }

    if (fplog != nullptr)
    {
        /* Print references after all software/hardware printing */
        please_cite(fplog, "Abraham2015");
        please_cite(fplog, "Pall2015");
        please_cite(fplog, "Pronk2013");
        please_cite(fplog, "Hess2008b");
        please_cite(fplog, "Spoel2005a");
        please_cite(fplog, "Lindahl2001a");
        please_cite(fplog, "Berendsen95a");
        writeSourceDoi(fplog);
    }

    std::unique_ptr<t_state> globalState;

    if (SIMMASTER(cr))
    {
        /* Only the master rank has the global state */
        globalState = compat::make_unique<t_state>();

        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPR, filenames.size(), filenames.data()), inputrec, globalState.get(), &mtop);

        /* In rerun, set velocities to zero if present */
        if (doRerun && ((globalState->flags & (1 << estV)) != 0))
        {
            // rerun does not use velocities
            GMX_LOG(mdlog.info).asParagraph().appendText(
                    "Rerun trajectory contains velocities. Rerun does only evaluate "
                    "potential energy and forces. The velocities will be ignored.");
            for (int i = 0; i < globalState->natoms; i++)
            {
                clear_rvec(globalState->v[i]);
            }
            globalState->flags &= ~(1 << estV);
        }

        if (inputrec->cutoff_scheme != ecutsVERLET)
        {
            if (nstlist_cmdline > 0)
            {
                gmx_fatal(FARGS, "Can not set nstlist with the group cut-off scheme");
            }

            if (!compatibleGpus.empty())
            {
                GMX_LOG(mdlog.warning).asParagraph().appendText(
                        "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                        "      To use a GPU, set the mdp option: cutoff-scheme = Verlet");
            }
        }
    }

    /* Check and update the hardware options for internal consistency */
    check_and_update_hw_opt_1(mdlog, &hw_opt, cr, domdecOptions.numPmeRanks);

    /* Early check for externally set process affinity. */
    gmx_check_thread_affinity_set(mdlog, cr,
                                  &hw_opt, hwinfo->nthreads_hw_avail, FALSE);

    if (GMX_THREAD_MPI && SIMMASTER(cr))
    {
        if (domdecOptions.numPmeRanks > 0 && hw_opt.nthreads_tmpi <= 0)
        {
            gmx_fatal(FARGS, "You need to explicitly specify the number of MPI threads (-ntmpi) when using separate PME ranks");
        }

        /* Since the master knows the cut-off scheme, update hw_opt for this.
         * This is done later for normal MPI and also once more with tMPI
         * for all tMPI ranks.
         */
        check_and_update_hw_opt_2(&hw_opt, inputrec->cutoff_scheme);

        bool useGpuForNonbonded = false;
        bool useGpuForPme       = false;
        try
        {
            // If the user specified the number of ranks, then we must
            // respect that, but in default mode, we need to allow for
            // the number of GPUs to choose the number of ranks.
            auto canUseGpuForNonbonded = buildSupportsNonbondedOnGpu(nullptr);
            useGpuForNonbonded = decideWhetherToUseGpusForNonbondedWithThreadMpi
                    (nonbondedTarget, gpuIdsToUse, userGpuTaskAssignment, emulateGpuNonbonded,
                    canUseGpuForNonbonded,
                    inputrec->cutoff_scheme == ecutsVERLET,
                    gpuAccelerationOfNonbondedIsUseful(mdlog, inputrec, GMX_THREAD_MPI),
                    hw_opt.nthreads_tmpi);
            useGpuForPme = decideWhetherToUseGpusForPmeWithThreadMpi
                    (useGpuForNonbonded, pmeTarget, gpuIdsToUse, userGpuTaskAssignment,
                    *hwinfo, *inputrec, mtop, hw_opt.nthreads_tmpi, domdecOptions.numPmeRanks);

        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

        /* Determine how many thread-MPI ranks to start.
         *
         * TODO Over-writing the user-supplied value here does
         * prevent any possible subsequent checks from working
         * correctly. */
        hw_opt.nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                &hw_opt,
                                                gpuIdsToUse,
                                                useGpuForNonbonded,
                                                useGpuForPme,
                                                inputrec, &mtop,
                                                mdlog,
                                                doMembed);

        // Now start the threads for thread MPI.
        cr = spawnThreads(hw_opt.nthreads_tmpi);
        /* The main thread continues here with a new cr. We don't deallocate
           the old cr because other threads may still be reading it. */
        // TODO Both master and spawned threads call dup_tfn and
        // reinitialize_commrec_for_this_thread. Find a way to express
        // this better.
        physicalNodeComm = PhysicalNodeCommunicator(MPI_COMM_WORLD, gmx_physicalnode_id_hash());
    }
    // END OF CAUTION: cr and physicalNodeComm are now reliable

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(cr, inputrec, &mtop);
    }

    // Now each rank knows the inputrec that SIMMASTER read and used,
    // and (if applicable) cr->nnodes has been assigned the number of
    // thread-MPI ranks that have been chosen. The ranks can now all
    // run the task-deciding functions and will agree on the result
    // without needing to communicate.
    //
    // TODO Should we do the communication in debug mode to support
    // having an assertion?
    //
    // Note that these variables describe only their own node.
    //
    // Note that when bonded interactions run on a GPU they always run
    // alongside a nonbonded task, so do not influence task assignment
    // even though they affect the force calculation workload.
    bool useGpuForNonbonded = false;
    bool useGpuForPme       = false;
    bool useGpuForBonded    = false;
    try
    {
        // It's possible that there are different numbers of GPUs on
        // different nodes, which is the user's responsibilty to
        // handle. If unsuitable, we will notice that during task
        // assignment.
        bool gpusWereDetected      = hwinfo->ngpu_compatible_tot > 0;
        bool usingVerletScheme     = inputrec->cutoff_scheme == ecutsVERLET;
        auto canUseGpuForNonbonded = buildSupportsNonbondedOnGpu(nullptr);
        useGpuForNonbonded = decideWhetherToUseGpusForNonbonded(nonbondedTarget, userGpuTaskAssignment,
                                                                emulateGpuNonbonded,
                                                                canUseGpuForNonbonded,
                                                                usingVerletScheme,
                                                                gpuAccelerationOfNonbondedIsUseful(mdlog, inputrec, !GMX_THREAD_MPI),
                                                                gpusWereDetected);
        useGpuForPme = decideWhetherToUseGpusForPme(useGpuForNonbonded, pmeTarget, userGpuTaskAssignment,
                                                    *hwinfo, *inputrec, mtop,
                                                    cr->nnodes, domdecOptions.numPmeRanks,
                                                    gpusWereDetected);
        auto canUseGpuForBonded = buildSupportsGpuBondeds(nullptr) && inputSupportsGpuBondeds(*inputrec, mtop, nullptr);
        useGpuForBonded =
            decideWhetherToUseGpusForBonded(useGpuForNonbonded, useGpuForPme, usingVerletScheme,
                                            bondedTarget, canUseGpuForBonded,
                                            EVDW_PME(inputrec->vdwtype),
                                            EEL_PME_EWALD(inputrec->coulombtype),
                                            domdecOptions.numPmeRanks, gpusWereDetected);

        pmeRunMode   = (useGpuForPme ? PmeRunMode::GPU : PmeRunMode::CPU);
        if (pmeRunMode == PmeRunMode::GPU)
        {
            if (pmeFftTarget == TaskTarget::Cpu)
            {
                pmeRunMode = PmeRunMode::Mixed;
            }
        }
        else if (pmeFftTarget == TaskTarget::Gpu)
        {
            gmx_fatal(FARGS, "Assigning FFTs to GPU requires PME to be assigned to GPU as well. With PME on CPU you should not be using -pmefft.");
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    // Build restraints.
    // TODO: hide restraint implementation details from Mdrunner.
    // There is nothing unique about restraints at this point as far as the
    // Mdrunner is concerned. The Mdrunner should just be getting a sequence of
    // factory functions from the SimulationContext on which to call mdModules->add().
    // TODO: capture all restraints into a single RestraintModule, passed to the runner builder.
    for (auto && restraint : restraintManager_->getRestraints())
    {
        auto module = RestraintMDModule::create(restraint,
                                                restraint->sites());
        mdModules->add(std::move(module));
    }

    // TODO: Error handling
    mdModules->assignOptionsToModules(*inputrec->params, nullptr);

    if (fplog != nullptr)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec, FALSE);
        fprintf(fplog, "\n");
    }

    if (SIMMASTER(cr))
    {
        /* now make sure the state is initialized and propagated */
        set_state_entries(globalState.get(), inputrec);
    }

    /* NM and TPI parallelize over force/energy calculations, not atoms,
     * so we need to initialize and broadcast the global state.
     */
    if (inputrec->eI == eiNM || inputrec->eI == eiTPI)
    {
        if (!MASTER(cr))
        {
            globalState = compat::make_unique<t_state>();
        }
        broadcastStateWithoutDynamics(cr, globalState.get());
    }

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(cr) && (domdecOptions.numCells[XX] > 1 ||
                     domdecOptions.numCells[YY] > 1 ||
                     domdecOptions.numCells[ZZ] > 1 ||
                     domdecOptions.numPmeRanks > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#if !GMX_MPI
                  "but %s was compiled without threads or MPI enabled", output_env_get_program_display_name(oenv));
#else
#if GMX_THREAD_MPI
                  "but the number of MPI-threads (option -ntmpi) is not set or is 1");
#else
                  "but %s was not started through mpirun/mpiexec or only one rank was requested through mpirun/mpiexec", output_env_get_program_display_name(oenv));
#endif
#endif
    }

    if (doRerun &&
        (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS, "The .mdp file specified an energy mininization or normal mode algorithm, and these are not compatible with mdrun -rerun");
    }

    if (can_use_allvsall(inputrec, TRUE, cr, fplog) && DOMAINDECOMP(cr))
    {
        gmx_fatal(FARGS, "All-vs-all loops do not work with domain decomposition, use a single MPI rank");
    }

    if (!(EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype)))
    {
        if (domdecOptions.numPmeRanks > 0)
        {
            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "PME-only ranks are requested, but the system does not use PME for electrostatics or LJ");
        }

        domdecOptions.numPmeRanks = 0;
    }

    if (useGpuForNonbonded && domdecOptions.numPmeRanks < 0)
    {
        /* With NB GPUs we don't automatically use PME-only CPU ranks. PME ranks can
         * improve performance with many threads per GPU, since our OpenMP
         * scaling is bad, but it's difficult to automate the setup.
         */
        domdecOptions.numPmeRanks = 0;
    }
    if (useGpuForPme)
    {
        if (domdecOptions.numPmeRanks < 0)
        {
            domdecOptions.numPmeRanks = 0;
            // TODO possibly print a note that one can opt-in for a separate PME GPU rank?
        }
        else
        {
            GMX_RELEASE_ASSERT(domdecOptions.numPmeRanks <= 1, "PME GPU decomposition is not supported");
        }
    }

#if GMX_FAHCORE
    if (MASTER(cr))
    {
        fcRegisterSteps(inputrec->nsteps, inputrec->init_step);
    }
#endif

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */
    snew(fcd, 1);

    /* This needs to be called before read_checkpoint to extend the state */
    init_disres(fplog, &mtop, inputrec, cr, ms, fcd, globalState.get(), replExParams.exchangeInterval > 0);

    init_orires(fplog, &mtop, inputrec, cr, ms, globalState.get(), &(fcd->orires));

    auto                 deform = prepareBoxDeformation(globalState->box, cr, *inputrec);

    ObservablesHistory   observablesHistory = {};

    ContinuationOptions &continuationOptions = mdrunOptions.continuationOptions;

    if (continuationOptions.startedFromCheckpoint)
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        gmx_bool bReadEkin;

        if (mdrunOptions.numStepsCommandline > -2)
        {
            /* Temporarily set the number of steps to unmlimited to avoid
             * triggering the nsteps check in load_checkpoint().
             * This hack will go away soon when the -nsteps option is removed.
             */
            inputrec->nsteps = -1;
        }

        load_checkpoint(opt2fn_master("-cpi", filenames.size(), filenames.data(), cr),
                        logFileHandle,
                        cr, domdecOptions.numCells,
                        inputrec, globalState.get(),
                        &bReadEkin, &observablesHistory,
                        continuationOptions.appendFiles,
                        continuationOptions.appendFilesOptionSet,
                        mdrunOptions.reproducible);

        if (bReadEkin)
        {
            continuationOptions.haveReadEkin = true;
        }

        if (continuationOptions.appendFiles && logFileHandle)
        {
            // Now we can start normal logging to the truncated log file.
            fplog    = gmx_fio_getfp(logFileHandle);
            prepareLogAppending(fplog);
            logOwner = buildLogger(fplog, cr);
            mdlog    = logOwner.logger();
        }
    }

    if (mdrunOptions.numStepsCommandline > -2)
    {
        GMX_LOG(mdlog.info).asParagraph().
            appendText("The -nsteps functionality is deprecated, and may be removed in a future version. "
                       "Consider using gmx convert-tpr -nsteps or changing the appropriate .mdp file field.");
    }
    /* override nsteps with value set on the commamdline */
    override_nsteps_cmdline(mdlog, mdrunOptions.numStepsCommandline, inputrec);

    if (SIMMASTER(cr))
    {
        copy_mat(globalState->box, box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box), box, cr);
    }

    /* Update rlist and nstlist. */
    if (inputrec->cutoff_scheme == ecutsVERLET)
    {
        prepare_verlet_scheme(fplog, cr, inputrec, nstlist_cmdline, &mtop, box,
                              useGpuForNonbonded || (emulateGpuNonbonded == EmulateGpuNonbonded::Yes), *hwinfo->cpuInfo);
    }

    LocalAtomSetManager atomSets;

    if (PAR(cr) && !(EI_TPI(inputrec->eI) ||
                     inputrec->eI == eiNM))
    {
        cr->dd = init_domain_decomposition(mdlog, cr, domdecOptions, mdrunOptions,
                                           &mtop, inputrec,
                                           box, positionsFromStatePointer(globalState.get()),
                                           &atomSets);
        // Note that local state still does not exist yet.
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty      = (DUTY_PP | DUTY_PME);

        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=screw is only implemented with domain decomposition");
        }
    }

    if (PAR(cr))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog, cr);
    }

#if GMX_MPI
    if (isMultiSim(ms))
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "This is simulation %d out of %d running as a composite GROMACS\n"
                "multi-simulation job. Setup for this simulation:\n",
                ms->sim, ms->nsim);
    }
    GMX_LOG(mdlog.warning).appendTextFormatted(
            "Using %d MPI %s\n",
            cr->nnodes,
#if GMX_THREAD_MPI
            cr->nnodes == 1 ? "thread" : "threads"
#else
            cr->nnodes == 1 ? "process" : "processes"
#endif
            );
    fflush(stderr);
#endif

    /* Check and update hw_opt for the cut-off scheme */
    check_and_update_hw_opt_2(&hw_opt, inputrec->cutoff_scheme);

    /* Check and update the number of OpenMP threads requested */
    checkAndUpdateRequestedNumOpenmpThreads(&hw_opt, *hwinfo, cr, ms, physicalNodeComm.size_,
                                            pmeRunMode, mtop);

    gmx_omp_nthreads_init(mdlog, cr,
                          hwinfo->nthreads_hw_avail,
                          physicalNodeComm.size_,
                          hw_opt.nthreads_omp,
                          hw_opt.nthreads_omp_pme,
                          !thisRankHasDuty(cr, DUTY_PP),
                          inputrec->cutoff_scheme == ecutsVERLET);

    // Enable FP exception detection for the Verlet scheme, but not in
    // Release mode and not for compilers with known buggy FP
    // exception support (clang with any optimization) or suspected
    // buggy FP exception support (gcc 7.* with optimization).
#if !defined NDEBUG && \
    !((defined __clang__ || (defined(__GNUC__) && !defined(__ICC) && __GNUC__ == 7)) \
    && defined __OPTIMIZE__)
    const bool bEnableFPE = inputrec->cutoff_scheme == ecutsVERLET;
#else
    const bool bEnableFPE = false;
#endif
    //FIXME - reconcile with gmx_feenableexcept() call from CommandLineModuleManager::run()
    if (bEnableFPE)
    {
        gmx_feenableexcept();
    }

    // Build a data structure that expresses which kinds of non-bonded
    // task are handled by this rank.
    //
    // TODO Later, this might become a loop over all registered modules
    // relevant to the mdp inputs, to find those that have such tasks.
    //
    // TODO This could move before init_domain_decomposition() as part
    // of refactoring that separates the responsibility for duty
    // assignment from setup for communication between tasks, and
    // setup for tasks handled with a domain (ie including short-ranged
    // tasks, bonded tasks, etc.).
    //
    // Note that in general useGpuForNonbonded, etc. can have a value
    // that is inconsistent with the presence of actual GPUs on any
    // rank, and that is not known to be a problem until the
    // duty of the ranks on a node become known.
    //
    // TODO Later we might need the concept of computeTasksOnThisRank,
    // from which we construct gpuTasksOnThisRank.
    //
    // Currently the DD code assigns duty to ranks that can
    // include PP work that currently can be executed on a single
    // GPU, if present and compatible.  This has to be coordinated
    // across PP ranks on a node, with possible multiple devices
    // or sharing devices on a node, either from the user
    // selection, or automatically.
    auto                 haveGpus = !gpuIdsToUse.empty();
    std::vector<GpuTask> gpuTasksOnThisRank;
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        if (useGpuForNonbonded)
        {
            // Note that any bonded tasks on a GPU always accompany a
            // non-bonded task.
            if (haveGpus)
            {
                gpuTasksOnThisRank.push_back(GpuTask::Nonbonded);
            }
            else if (nonbondedTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS, "Cannot run short-ranged nonbonded interactions on a GPU because there is none detected.");
            }
            else if (bondedTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS, "Cannot run bonded interactions on a GPU because there is none detected.");
            }
        }
    }
    // TODO cr->duty & DUTY_PME should imply that a PME algorithm is active, but currently does not.
    if (EEL_PME(inputrec->coulombtype) && (thisRankHasDuty(cr, DUTY_PME)))
    {
        if (useGpuForPme)
        {
            if (haveGpus)
            {
                gpuTasksOnThisRank.push_back(GpuTask::Pme);
            }
            else if (pmeTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS, "Cannot run PME on a GPU because there is none detected.");
            }
        }
    }

    GpuTaskAssignment gpuTaskAssignment;
    try
    {
        // Produce the task assignment for this rank.
        gpuTaskAssignment = runTaskAssignment(gpuIdsToUse, userGpuTaskAssignment, *hwinfo,
                                              mdlog, cr, ms, physicalNodeComm, gpuTasksOnThisRank,
                                              useGpuForBonded, pmeRunMode);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    /* Prevent other ranks from continuing after an issue was found
     * and reported as a fatal error.
     *
     * TODO This function implements a barrier so that MPI runtimes
     * can organize an orderly shutdown if one of the ranks has had to
     * issue a fatal error in various code already run. When we have
     * MPI-aware error handling and reporting, this should be
     * improved. */
#if GMX_MPI
    if (PAR(cr))
    {
        MPI_Barrier(cr->mpi_comm_mysim);
    }
    if (isMultiSim(ms))
    {
        if (SIMMASTER(cr))
        {
            MPI_Barrier(ms->mpi_comm_masters);
        }
        /* We need another barrier to prevent non-master ranks from contiuing
         * when an error occured in a different simulation.
         */
        MPI_Barrier(cr->mpi_comm_mysim);
    }
#endif

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, !gpuTaskAssignment.empty(), mdrunOptions.ntompOptionIsSet,
                                       cr, mdlog);

    gmx_device_info_t *nonbondedDeviceInfo = nullptr;

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        // This works because only one task of each type is currently permitted.
        auto nbGpuTaskMapping = std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(),
                                             hasTaskType<GpuTask::Nonbonded>);
        if (nbGpuTaskMapping != gpuTaskAssignment.end())
        {
            int nonbondedDeviceId = nbGpuTaskMapping->deviceId_;
            nonbondedDeviceInfo = getDeviceInfo(hwinfo->gpu_info, nonbondedDeviceId);
            init_gpu(nonbondedDeviceInfo);

            if (DOMAINDECOMP(cr))
            {
                /* When we share GPUs over ranks, we need to know this for the DLB */
                dd_setup_dlb_resource_sharing(cr, nonbondedDeviceId);
            }

        }
    }

    std::unique_ptr<ClfftInitializer> initializedClfftLibrary;

    gmx_device_info_t                *pmeDeviceInfo = nullptr;
    // Later, this program could contain kernels that might be later
    // re-used as auto-tuning progresses, or subsequent simulations
    // are invoked.
    PmeGpuProgramStorage pmeGpuProgram;
    // This works because only one task of each type is currently permitted.
    auto                 pmeGpuTaskMapping     = std::find_if(gpuTaskAssignment.begin(), gpuTaskAssignment.end(), hasTaskType<GpuTask::Pme>);
    const bool           thisRankHasPmeGpuTask = (pmeGpuTaskMapping != gpuTaskAssignment.end());
    if (thisRankHasPmeGpuTask)
    {
        pmeDeviceInfo = getDeviceInfo(hwinfo->gpu_info, pmeGpuTaskMapping->deviceId_);
        init_gpu(pmeDeviceInfo);
        pmeGpuProgram = buildPmeGpuProgram(pmeDeviceInfo);
        // TODO It would be nice to move this logic into the factory
        // function. See Redmine #2535.
        bool isMasterThread = !GMX_THREAD_MPI || MASTER(cr);
        if (pmeRunMode == PmeRunMode::GPU && !initializedClfftLibrary && isMasterThread)
        {
            initializedClfftLibrary = initializeClfftLibrary();
        }
    }

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

    int numThreadsOnThisRank;
    /* threads on this MPI process or TMPI thread */
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        numThreadsOnThisRank = gmx_omp_nthreads_get(emntNonbonded);
    }
    else
    {
        numThreadsOnThisRank = nthreads_pme;
    }

    checkHardwareOversubscription(numThreadsOnThisRank, cr->nodeid,
                                  *hwinfo->hardwareTopology,
                                  physicalNodeComm, mdlog);

    if (hw_opt.thread_affinity != threadaffOFF)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(mdlog, cr,
                                      &hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        int numThreadsOnThisNode, intraNodeThreadOffset;
        analyzeThreadsOnThisNode(physicalNodeComm, numThreadsOnThisRank, &numThreadsOnThisNode,
                                 &intraNodeThreadOffset);

        /* Set the CPU affinity */
        gmx_set_thread_affinity(mdlog, cr, &hw_opt, *hwinfo->hardwareTopology,
                                numThreadsOnThisRank, numThreadsOnThisNode,
                                intraNodeThreadOffset, nullptr);
    }

    if (mdrunOptions.timingOptions.resetStep > -1)
    {
        GMX_LOG(mdlog.info).asParagraph().
            appendText("The -resetstep functionality is deprecated, and may be removed in a future version.");
    }
    wcycle = wallcycle_init(fplog, mdrunOptions.timingOptions.resetStep, cr);

    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast_sim(sizeof(reset_counters), &reset_counters, cr);
        wcycle_set_reset_counters(wcycle, reset_counters);
    }

    // Membrane embedding must be initialized before we call init_forcerec()
    if (doMembed)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "Initializing membed");
        }
        /* Note that membed cannot work in parallel because mtop is
         * changed here. Fix this if we ever want to make it run with
         * multiple ranks. */
        membed = init_membed(fplog, filenames.size(), filenames.data(), &mtop, inputrec, globalState.get(), cr,
                             &mdrunOptions
                                 .checkpointOptions.period);
    }

    std::unique_ptr<MDAtoms>     mdAtoms;
    std::unique_ptr<gmx_vsite_t> vsite;

    snew(nrnb, 1);
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Initiate forcerecord */
        fr                 = mk_forcerec();
        fr->forceProviders = mdModules->initForceProviders();
        init_forcerec(fplog, mdlog, fr, fcd,
                      inputrec, &mtop, cr, box,
                      opt2fn("-table", filenames.size(), filenames.data()),
                      opt2fn("-tablep", filenames.size(), filenames.data()),
                      opt2fns("-tableb", filenames.size(), filenames.data()),
                      *hwinfo, nonbondedDeviceInfo,
                      useGpuForBonded,
                      FALSE,
                      pforce);

        /* Initialize the mdAtoms structure.
         * mdAtoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdAtoms = makeMDAtoms(fplog, mtop, *inputrec, thisRankHasPmeGpuTask);
        if (globalState && thisRankHasPmeGpuTask)
        {
            // The pinning of coordinates in the global state object works, because we only use
            // PME on GPU without DD or on a separate PME rank, and because the local state pointer
            // points to the global state object without DD.
            // FIXME: MD and EM separately set up the local state - this should happen in the same function,
            // which should also perform the pinning.
            changePinningPolicy(&globalState->x, pme_get_pinning_policy());
        }

        /* Initialize the virtual site communication */
        vsite = initVsite(mtop, cr);

        calc_shifts(box, fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr) &&
            !(inputrec->ePBC != epbcNONE && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->ePBC != epbcNONE)
            {
                do_pbc_first_mtop(fplog, inputrec->ePBC, box, &mtop, globalState->x.rvec_array());
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                constructVsitesGlobal(mtop, globalState->x);
            }
        }

        if (EEL_PME(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype))
        {
            ewaldcoeff_q  = fr->ic->ewaldcoeff_q;
            ewaldcoeff_lj = fr->ic->ewaldcoeff_lj;
        }
    }
    else
    {
        /* This is a PME only node */

        GMX_ASSERT(globalState == nullptr, "We don't need the state on a PME only rank and expect it to be unitialized");

        ewaldcoeff_q  = calc_ewaldcoeff_q(inputrec->rcoulomb, inputrec->ewald_rtol);
        ewaldcoeff_lj = calc_ewaldcoeff_lj(inputrec->rvdw, inputrec->ewald_rtol_lj);
    }

    gmx_pme_t *sepPmeData = nullptr;
    // This reference hides the fact that PME data is owned by runner on PME-only ranks and by forcerec on other ranks
    GMX_ASSERT(thisRankHasDuty(cr, DUTY_PP) == (fr != nullptr), "Double-checking that only PME-only ranks have no forcerec");
    gmx_pme_t * &pmedata = fr ? fr->pmedata : sepPmeData;

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype))
    {
        if (mdAtoms && mdAtoms->mdatoms())
        {
            nChargePerturbed = mdAtoms->mdatoms()->nChargePerturbed;
            if (EVDW_PME(inputrec->vdwtype))
            {
                nTypePerturbed   = mdAtoms->mdatoms()->nTypePerturbed;
            }
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed(FEP on Q) and nTypePerturbed(FEP on LJ)*/
            gmx_bcast_sim(sizeof(nChargePerturbed), &nChargePerturbed, cr);
            gmx_bcast_sim(sizeof(nTypePerturbed), &nTypePerturbed, cr);
        }

        if (thisRankHasDuty(cr, DUTY_PME))
        {
            try
            {
                pmedata = gmx_pme_init(cr,
                                       getNumPmeDomains(cr->dd),
                                       inputrec,
                                       mtop.natoms, nChargePerturbed != 0, nTypePerturbed != 0,
                                       mdrunOptions.reproducible,
                                       ewaldcoeff_q, ewaldcoeff_lj,
                                       nthreads_pme,
                                       pmeRunMode, nullptr,
                                       pmeDeviceInfo, pmeGpuProgram.get(), mdlog);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        }
    }


    if (EI_DYNAMICS(inputrec->eI))
    {
        /* Turn on signal handling on all nodes */
        /*
         * (A user signal from the PME nodes (if any)
         * is communicated to the PP nodes.
         */
        signal_handler_install();
    }

    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Assumes uniform use of the number of OpenMP threads */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntDefault));

        if (inputrec->bPull)
        {
            /* Initialize pull code */
            inputrec->pull_work =
                init_pull(fplog, inputrec->pull, inputrec,
                          &mtop, cr, &atomSets, inputrec->fepvals->init_lambda);
            if (inputrec->pull->bXOutAverage || inputrec->pull->bFOutAverage)
            {
                initPullHistory(inputrec->pull_work, &observablesHistory);
            }
            if (EI_DYNAMICS(inputrec->eI) && MASTER(cr))
            {
                init_pull_output_files(inputrec->pull_work,
                                       filenames.size(), filenames.data(), oenv,
                                       continuationOptions);
            }
        }

        std::unique_ptr<EnforcedRotation> enforcedRotation;
        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            enforcedRotation = init_rot(fplog,
                                        inputrec,
                                        filenames.size(),
                                        filenames.data(),
                                        cr,
                                        &atomSets,
                                        globalState.get(),
                                        &mtop,
                                        oenv,
                                        mdrunOptions);
        }

        if (inputrec->eSwapCoords != eswapNO)
        {
            /* Initialize ion swapping code */
            init_swapcoords(fplog, inputrec, opt2fn_master("-swap", filenames.size(), filenames.data(), cr),
                            &mtop, globalState.get(), &observablesHistory,
                            cr, &atomSets, oenv, mdrunOptions);
        }

        /* Let makeConstraints know whether we have essential dynamics constraints.
         * TODO: inputrec should tell us whether we use an algorithm, not a file option or the checkpoint
         */
        bool doEssentialDynamics = (opt2fn_null("-ei", filenames.size(), filenames.data()) != nullptr
                                    || observablesHistory.edsamHistory);
        auto constr              = makeConstraints(mtop, *inputrec, doEssentialDynamics,
                                                   fplog, *mdAtoms->mdatoms(),
                                                   cr, ms, nrnb, wcycle, fr->bMolPBC);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because fr->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog, cr->dd, &mtop, vsite.get(), inputrec,
                            domdecOptions.checkBondedInteractions,
                            fr->cginfo_mb);
        }

        // TODO This is not the right place to manage the lifetime of
        // this data structure, but currently it's the easiest way to
        // make it work. Later, it should probably be made/updated
        // after the workload for the lifetime of a PP domain is
        // understood.
        PpForceWorkload ppForceWorkload;

        GMX_ASSERT(stopHandlerBuilder_, "Runner must provide StopHandlerBuilder to integrator.");
        /* Now do whatever the user wants us to do (how flexible...) */
        Integrator integrator {
            fplog, cr, ms, mdlog, static_cast<int>(filenames.size()), filenames.data(),
            oenv,
            mdrunOptions,
            vsite.get(), constr.get(),
            enforcedRotation ? enforcedRotation->getLegacyEnfrot() : nullptr,
            deform.get(),
            mdModules->outputProvider(),
            inputrec, &mtop,
            fcd,
            globalState.get(),
            &observablesHistory,
            mdAtoms.get(), nrnb, wcycle, fr,
            &ppForceWorkload,
            replExParams,
            membed,
            walltime_accounting,
            std::move(stopHandlerBuilder_)
        };
        integrator.run(inputrec->eI, doRerun);

        if (inputrec->bPull)
        {
            finish_pull(inputrec->pull_work);
        }

    }
    else
    {
        GMX_RELEASE_ASSERT(pmedata, "pmedata was NULL while cr->duty was not DUTY_PP");
        /* do PME only */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntPME));
        gmx_pmeonly(pmedata, cr, nrnb, wcycle, walltime_accounting, inputrec, pmeRunMode);
    }

    wallcycle_stop(wcycle, ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog, mdlog, cr,
               inputrec, nrnb, wcycle, walltime_accounting,
               fr ? fr->nbv : nullptr,
               pmedata,
               EI_DYNAMICS(inputrec->eI) && !isMultiSim(ms));

    // Free PME data
    if (pmedata)
    {
        gmx_pme_destroy(pmedata);
        pmedata = nullptr;
    }

    // FIXME: this is only here to manually unpin mdAtoms->chargeA_ and state->x,
    // before we destroy the GPU context(s) in free_gpu_resources().
    // Pinned buffers are associated with contexts in CUDA.
    // As soon as we destroy GPU contexts after mdrunner() exits, these lines should go.
    mdAtoms.reset(nullptr);
    globalState.reset(nullptr);
    mdModules.reset(nullptr);   // destruct force providers here as they might also use the GPU

    /* Free GPU memory and set a physical node tMPI barrier (which should eventually go away) */
    free_gpu_resources(fr, physicalNodeComm);
    free_gpu(nonbondedDeviceInfo);
    free_gpu(pmeDeviceInfo);
    done_forcerec(fr, mtop.molblock.size(), mtop.groups.grps[egcENER].nr);
    sfree(fcd);

    if (doMembed)
    {
        free_membed(membed);
    }

    gmx_hardware_info_free();

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltime_accounting);
    sfree(nrnb);

    // Ensure log file content is written
    if (logFileHandle)
    {
        gmx_fio_flush(logFileHandle);
    }

    /* Reset FPEs (important for unit tests) by disabling them. Assumes no
     * exceptions were enabled before function was called. */
    if (bEnableFPE)
    {
        gmx_fedisableexcept();
    }

    rc = static_cast<int>(gmx_get_stop_condition());

#if GMX_THREAD_MPI
    /* we need to join all threads. The sub-threads join when they
       exit this function, but the master thread needs to be told to
       wait for that. */
    if (PAR(cr) && MASTER(cr))
    {
        done_commrec(cr);
        tMPI_Finalize();
    }
#endif

    return rc;
}

Mdrunner::~Mdrunner()
{
    // Clean up of the Manager.
    // This will end up getting called on every thread-MPI rank, which is unnecessary,
    // but okay as long as threads synchronize some time before adding or accessing
    // a new set of restraints.
    if (restraintManager_)
    {
        restraintManager_->clear();
        GMX_ASSERT(restraintManager_->countRestraints() == 0,
                   "restraints added during runner life time should be cleared at runner destruction.");
    }
};

void Mdrunner::addPotential(std::shared_ptr<gmx::IRestraintPotential> puller,
                            std::string                               name)
{
    GMX_ASSERT(restraintManager_, "Mdrunner must have a restraint manager.");
    // Not sure if this should be logged through the md logger or something else,
    // but it is helpful to have some sort of INFO level message sent somewhere.
    //    std::cout << "Registering restraint named " << name << std::endl;

    // When multiple restraints are used, it may be wasteful to register them separately.
    // Maybe instead register an entire Restraint Manager as a force provider.
    restraintManager_->addToSpec(std::move(puller),
                                 std::move(name));
}

Mdrunner::Mdrunner(Mdrunner &&) noexcept = default;

//NOLINTNEXTLINE(performance-noexcept-move-constructor) working around GCC bug 58265
Mdrunner &Mdrunner::operator=(Mdrunner && /*handle*/) noexcept(BUGFREE_NOEXCEPT_STRING) = default;

class Mdrunner::BuilderImplementation
{
    public:
        BuilderImplementation() = delete;
        explicit BuilderImplementation(SimulationContext* context);
        ~BuilderImplementation();

        BuilderImplementation &setExtraMdrunOptions(const MdrunOptions &options,
                                                    real                forceWarningThreshold);

        void addDomdec(const DomdecOptions &options);

        void addVerletList(int nstlist);

        void addReplicaExchange(const ReplicaExchangeParameters &params);

        void addMultiSim(gmx_multisim_t* multisim);

        void addNonBonded(const char* nbpu_opt);

        void addPME(const char* pme_opt_, const char* pme_fft_opt_);

        void addBondedTaskAssignment(const char* bonded_opt);

        void addHardwareOptions(const gmx_hw_opt_t &hardwareOptions);

        void addFilenames(ArrayRef <const t_filenm> filenames);

        void addOutputEnvironment(gmx_output_env_t* outputEnvironment);

        void addLogFile(t_fileio *logFileHandle);

        void addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder);

        Mdrunner build();

    private:
        // Default parameters copied from runner.h
        // \todo Clarify source(s) of default parameters.

        const char* nbpu_opt_    = nullptr;
        const char* pme_opt_     = nullptr;
        const char* pme_fft_opt_ = nullptr;
        const char *bonded_opt_  = nullptr;

        MdrunOptions                          mdrunOptions_;

        DomdecOptions                         domdecOptions_;

        ReplicaExchangeParameters             replicaExchangeParameters_;

        //! Command-line override for the duration of a neighbor list with the Verlet scheme.
        int         nstlist_ = 0;

        //! Non-owning multisim communicator handle.
        std::unique_ptr<gmx_multisim_t*>      multisim_ = nullptr;

        //! Print a warning if any force is larger than this (in kJ/mol nm).
        real forceWarningThreshold_ = -1;

        /*! \brief  Non-owning pointer to SimulationContext (owned and managed by client)
         *
         * \internal
         * \todo Establish robust protocol to make sure resources remain valid.
         * SimulationContext will likely be separated into multiple layers for
         * different levels of access and different phases of execution. Ref
         * https://redmine.gromacs.org/issues/2375
         * https://redmine.gromacs.org/issues/2587
         */
        SimulationContext* context_ = nullptr;

        //! \brief Parallelism information.
        gmx_hw_opt_t hardwareOptions_;

        //! filename options for simulation.
        ArrayRef<const t_filenm> filenames_;

        /*! \brief Handle to output environment.
         *
         * \todo gmx_output_env_t needs lifetime management.
         */
        gmx_output_env_t*    outputEnvironment_ = nullptr;

        /*! \brief Non-owning handle to MD log file.
         *
         * \todo Context should own output facilities for client.
         * \todo Improve log file handle management.
         * \internal
         * Code managing the FILE* relies on the ability to set it to
         * nullptr to check whether the filehandle is valid.
         */
        t_fileio* logFileHandle_ = nullptr;

        /*!
         * \brief Builder for simulation stop signal handler.
         */
        std::unique_ptr<StopHandlerBuilder> stopHandlerBuilder_ = nullptr;
};

Mdrunner::BuilderImplementation::BuilderImplementation(SimulationContext* context) :
    context_(context)
{
    GMX_ASSERT(context_, "Bug found. It should not be possible to construct builder without a valid context.");
}

Mdrunner::BuilderImplementation::~BuilderImplementation() = default;

Mdrunner::BuilderImplementation &
Mdrunner::BuilderImplementation::setExtraMdrunOptions(const MdrunOptions &options,
                                                      real                forceWarningThreshold)
{
    mdrunOptions_          = options;
    forceWarningThreshold_ = forceWarningThreshold;
    return *this;
}

void Mdrunner::BuilderImplementation::addDomdec(const DomdecOptions &options)
{
    domdecOptions_ = options;
}

void Mdrunner::BuilderImplementation::addVerletList(int nstlist)
{
    nstlist_ = nstlist;
}

void Mdrunner::BuilderImplementation::addReplicaExchange(const ReplicaExchangeParameters &params)
{
    replicaExchangeParameters_ = params;
}

void Mdrunner::BuilderImplementation::addMultiSim(gmx_multisim_t* multisim)
{
    multisim_ = compat::make_unique<gmx_multisim_t*>(multisim);
}

Mdrunner Mdrunner::BuilderImplementation::build()
{
    auto newRunner = Mdrunner();

    GMX_ASSERT(context_, "Bug found. It should not be possible to call build() without a valid context.");

    newRunner.mdrunOptions    = mdrunOptions_;
    newRunner.domdecOptions   = domdecOptions_;

    // \todo determine an invariant to check or confirm that all gmx_hw_opt_t objects are valid
    newRunner.hw_opt          = hardwareOptions_;

    // No invariant to check. This parameter exists to optionally override other behavior.
    newRunner.nstlist_cmdline = nstlist_;

    newRunner.replExParams    = replicaExchangeParameters_;

    newRunner.filenames = filenames_;

    GMX_ASSERT(context_->communicationRecord_, "SimulationContext communications not initialized.");
    newRunner.cr = context_->communicationRecord_;

    if (multisim_)
    {
        // nullptr is a valid value for the multisim handle, so we don't check the pointed-to pointer.
        newRunner.ms = *multisim_;
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addMultiSim() is required before build()"));
    }

    // \todo Clarify ownership and lifetime management for gmx_output_env_t
    // \todo Update sanity checking when output environment has clearly specified invariants.
    // Initialization and default values for oenv are not well specified in the current version.
    if (outputEnvironment_)
    {
        newRunner.oenv  = outputEnvironment_;
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addOutputEnvironment() is required before build()"));
    }

    newRunner.logFileHandle = logFileHandle_;

    if (nbpu_opt_)
    {
        newRunner.nbpu_opt    = nbpu_opt_;
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addNonBonded() is required before build()"));
    }

    if (pme_opt_ && pme_fft_opt_)
    {
        newRunner.pme_opt     = pme_opt_;
        newRunner.pme_fft_opt = pme_fft_opt_;
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addElectrostatics() is required before build()"));
    }

    if (bonded_opt_)
    {
        newRunner.bonded_opt = bonded_opt_;
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addBondedTaskAssignment() is required before build()"));
    }

    newRunner.restraintManager_ = compat::make_unique<gmx::RestraintManager>();

    if (stopHandlerBuilder_)
    {
        newRunner.stopHandlerBuilder_ = std::move(stopHandlerBuilder_);
    }
    else
    {
        newRunner.stopHandlerBuilder_ = compat::make_unique<StopHandlerBuilder>();
    }

    return newRunner;
}

void Mdrunner::BuilderImplementation::addNonBonded(const char* nbpu_opt)
{
    nbpu_opt_ = nbpu_opt;
}

void Mdrunner::BuilderImplementation::addPME(const char* pme_opt,
                                             const char* pme_fft_opt)
{
    pme_opt_     = pme_opt;
    pme_fft_opt_ = pme_fft_opt;
}

void Mdrunner::BuilderImplementation::addBondedTaskAssignment(const char* bonded_opt)
{
    bonded_opt_ = bonded_opt;
}

void Mdrunner::BuilderImplementation::addHardwareOptions(const gmx_hw_opt_t &hardwareOptions)
{
    hardwareOptions_ = hardwareOptions;
}

void Mdrunner::BuilderImplementation::addFilenames(ArrayRef<const t_filenm> filenames)
{
    filenames_ = filenames;
}

void Mdrunner::BuilderImplementation::addOutputEnvironment(gmx_output_env_t* outputEnvironment)
{
    outputEnvironment_ = outputEnvironment;
}

void Mdrunner::BuilderImplementation::addLogFile(t_fileio *logFileHandle)
{
    logFileHandle_ = logFileHandle;
}

void Mdrunner::BuilderImplementation::addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder)
{
    stopHandlerBuilder_ = std::move(builder);
}

MdrunnerBuilder::MdrunnerBuilder(compat::not_null<SimulationContext*> context) :
    impl_ {gmx::compat::make_unique<Mdrunner::BuilderImplementation>(context)}
{
}

MdrunnerBuilder::~MdrunnerBuilder() = default;

MdrunnerBuilder &MdrunnerBuilder::addSimulationMethod(const MdrunOptions &options,
                                                      real                forceWarningThreshold)
{
    impl_->setExtraMdrunOptions(options, forceWarningThreshold);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addDomainDecomposition(const DomdecOptions &options)
{
    impl_->addDomdec(options);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addNeighborList(int nstlist)
{
    impl_->addVerletList(nstlist);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addReplicaExchange(const ReplicaExchangeParameters &params)
{
    impl_->addReplicaExchange(params);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addMultiSim(gmx_multisim_t* multisim)
{
    impl_->addMultiSim(multisim);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addNonBonded(const char* nbpu_opt)
{
    impl_->addNonBonded(nbpu_opt);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addElectrostatics(const char* pme_opt,
                                                    const char* pme_fft_opt)
{
    // The builder method may become more general in the future, but in this version,
    // parameters for PME electrostatics are both required and the only parameters
    // available.
    if (pme_opt && pme_fft_opt)
    {
        impl_->addPME(pme_opt, pme_fft_opt);
    }
    else
    {
        GMX_THROW(gmx::InvalidInputError("addElectrostatics() arguments must be non-null pointers."));
    }
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addBondedTaskAssignment(const char* bonded_opt)
{
    impl_->addBondedTaskAssignment(bonded_opt);
    return *this;
}

Mdrunner MdrunnerBuilder::build()
{
    return impl_->build();
}

MdrunnerBuilder &MdrunnerBuilder::addHardwareOptions(const gmx_hw_opt_t &hardwareOptions)
{
    impl_->addHardwareOptions(hardwareOptions);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addFilenames(ArrayRef<const t_filenm> filenames)
{
    impl_->addFilenames(filenames);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addOutputEnvironment(gmx_output_env_t* outputEnvironment)
{
    impl_->addOutputEnvironment(outputEnvironment);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addLogFile(t_fileio *logFileHandle)
{
    impl_->addLogFile(logFileHandle);
    return *this;
}

MdrunnerBuilder &MdrunnerBuilder::addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder)
{
    impl_->addStopHandlerBuilder(std::move(builder));
    return *this;
}

MdrunnerBuilder::MdrunnerBuilder(MdrunnerBuilder &&) noexcept = default;

MdrunnerBuilder &MdrunnerBuilder::operator=(MdrunnerBuilder &&) noexcept = default;

} // namespace gmx
