/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011-2019,2020, by the GROMACS development team, led by
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
#include <memory>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/builder.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme_gpu_program.h"
#include "gromacs/ewald/pme_only.h"
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_context.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/printhardware.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/listed_forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist_tuning.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/restraint/manager.h"
#include "gromacs/restraint/restraintmdmodule.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/taskassignment/decidegpuusage.h"
#include "gromacs/taskassignment/decidesimulationworkload.h"
#include "gromacs/taskassignment/resourcedivision.h"
#include "gromacs/taskassignment/taskassignment.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/wallcyclereporting.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/mdmodulenotification.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "isimulator.h"
#include "replicaexchange.h"
#include "simulatorbuilder.h"

#if GMX_FAHCORE
#    include "corewrap.h"
#endif

namespace gmx
{


/*! \brief Manage any development feature flag variables encountered
 *
 * The use of dev features indicated by environment variables is
 * logged in order to ensure that runs with such features enabled can
 * be identified from their log and standard output. Any cross
 * dependencies are also checked, and if unsatisfied, a fatal error
 * issued.
 *
 * Note that some development features overrides are applied already here:
 * the GPU communication flags are set to false in non-tMPI and non-CUDA builds.
 *
 * \param[in]  mdlog                Logger object.
 * \param[in]  useGpuForNonbonded   True if the nonbonded task is offloaded in this run.
 * \param[in]  pmeRunMode           The PME run mode for this run
 * \returns                         The object populated with development feature flags.
 */
static DevelopmentFeatureFlags manageDevelopmentFeatures(const gmx::MDLogger& mdlog,
                                                         const bool           useGpuForNonbonded,
                                                         const PmeRunMode     pmeRunMode)
{
    DevelopmentFeatureFlags devFlags;

    // Some builds of GCC 5 give false positive warnings that these
    // getenv results are ignored when clearly they are used.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
    devFlags.enableGpuBufferOps = (getenv("GMX_USE_GPU_BUFFER_OPS") != nullptr)
                                  && (GMX_GPU == GMX_GPU_CUDA) && useGpuForNonbonded;
    devFlags.forceGpuUpdateDefault = (getenv("GMX_FORCE_UPDATE_DEFAULT_GPU") != nullptr);
    devFlags.enableGpuHaloExchange =
            (getenv("GMX_GPU_DD_COMMS") != nullptr && GMX_THREAD_MPI && (GMX_GPU == GMX_GPU_CUDA));
    devFlags.enableGpuPmePPComm =
            (getenv("GMX_GPU_PME_PP_COMMS") != nullptr && GMX_THREAD_MPI && (GMX_GPU == GMX_GPU_CUDA));
#pragma GCC diagnostic pop

    if (devFlags.enableGpuBufferOps)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "This run uses the 'GPU buffer ops' feature, enabled by the "
                        "GMX_USE_GPU_BUFFER_OPS environment variable.");
    }

    if (devFlags.forceGpuUpdateDefault)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "This run will default to '-update gpu' as requested by the "
                        "GMX_FORCE_UPDATE_DEFAULT_GPU environment variable. GPU update with domain "
                        "decomposition lacks substantial testing and should be used with caution.");
    }

    if (devFlags.enableGpuHaloExchange)
    {
        if (useGpuForNonbonded)
        {
            if (!devFlags.enableGpuBufferOps)
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendTextFormatted(
                                "Enabling GPU buffer operations required by GMX_GPU_DD_COMMS "
                                "(equivalent with GMX_USE_GPU_BUFFER_OPS=1).");
                devFlags.enableGpuBufferOps = true;
            }
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "This run has requested the 'GPU halo exchange' feature, enabled by "
                            "the "
                            "GMX_GPU_DD_COMMS environment variable.");
        }
        else
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "GMX_GPU_DD_COMMS environment variable detected, but the 'GPU "
                            "halo exchange' feature will not be enabled as nonbonded interactions "
                            "are not offloaded.");
            devFlags.enableGpuHaloExchange = false;
        }
    }

    if (devFlags.enableGpuPmePPComm)
    {
        if (pmeRunMode == PmeRunMode::GPU)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "This run uses the 'GPU PME-PP communications' feature, enabled "
                            "by the GMX_GPU_PME_PP_COMMS environment variable.");
        }
        else
        {
            std::string clarification;
            if (pmeRunMode == PmeRunMode::Mixed)
            {
                clarification =
                        "PME FFT and gather are not offloaded to the GPU (PME is running in mixed "
                        "mode).";
            }
            else
            {
                clarification = "PME is not offloaded to the GPU.";
            }
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GMX_GPU_PME_PP_COMMS environment variable detected, but the "
                            "'GPU PME-PP communications' feature was not enabled as "
                            + clarification);
            devFlags.enableGpuPmePPComm = false;
        }
    }

    return devFlags;
}

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
    auto newRunner = Mdrunner(std::make_unique<MDModules>());

    // All runners in the same process share a restraint manager resource because it is
    // part of the interface to the client code, which is associated only with the
    // original thread. Handles to the same resources can be obtained by copy.
    {
        newRunner.restraintManager_ = std::make_unique<RestraintManager>(*restraintManager_);
    }

    // Copy members of master runner.
    // \todo Replace with builder when Simulation context and/or runner phases are better defined.
    // Ref https://gitlab.com/gromacs/gromacs/-/issues/2587 and https://gitlab.com/gromacs/gromacs/-/issues/2375
    newRunner.hw_opt    = hw_opt;
    newRunner.filenames = filenames;

    newRunner.oenv            = oenv;
    newRunner.mdrunOptions    = mdrunOptions;
    newRunner.domdecOptions   = domdecOptions;
    newRunner.nbpu_opt        = nbpu_opt;
    newRunner.pme_opt         = pme_opt;
    newRunner.pme_fft_opt     = pme_fft_opt;
    newRunner.bonded_opt      = bonded_opt;
    newRunner.update_opt      = update_opt;
    newRunner.nstlist_cmdline = nstlist_cmdline;
    newRunner.replExParams    = replExParams;
    newRunner.pforce          = pforce;
    // Give the spawned thread the newly created valid communicator
    // for the simulation.
    newRunner.communicator        = MPI_COMM_WORLD;
    newRunner.ms                  = ms;
    newRunner.startingBehavior    = startingBehavior;
    newRunner.stopHandlerBuilder_ = std::make_unique<StopHandlerBuilder>(*stopHandlerBuilder_);

    threadMpiMdrunnerAccessBarrier();

    return newRunner;
}

/*! \brief The callback used for running on spawned threads.
 *
 * Obtains the pointer to the master mdrunner object from the one
 * argument permitted to the thread-launch API call, copies it to make
 * a new runner for this thread, reinitializes necessary data, and
 * proceeds to the simulation. */
static void mdrunner_start_fn(const void* arg)
{
    try
    {
        auto masterMdrunner = reinterpret_cast<const gmx::Mdrunner*>(arg);
        /* copy the arg list to make sure that it's thread-local. This
           doesn't copy pointed-to items, of course; fnm, cr and fplog
           are reset in the call below, all others should be const. */
        gmx::Mdrunner mdrunner = masterMdrunner->cloneOnSpawnedThread();
        mdrunner.mdrunner();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}


void Mdrunner::spawnThreads(int numThreadsToLaunch)
{
#if GMX_THREAD_MPI
    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns. Thread affinity is handled later. */
    if (tMPI_Init_fn(TRUE, numThreadsToLaunch, TMPI_AFFINITY_NONE, mdrunner_start_fn,
                     static_cast<const void*>(this))
        != TMPI_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to spawn thread-MPI threads"));
    }

    // Give the master thread the newly created valid communicator for
    // the simulation.
    communicator = MPI_COMM_WORLD;
    threadMpiMdrunnerAccessBarrier();
#else
    GMX_UNUSED_VALUE(numThreadsToLaunch);
    GMX_UNUSED_VALUE(mdrunner_start_fn);
#endif
}

} // namespace gmx

/*! \brief Initialize variables for Verlet scheme simulation */
static void prepare_verlet_scheme(FILE*               fplog,
                                  t_commrec*          cr,
                                  t_inputrec*         ir,
                                  int                 nstlist_cmdline,
                                  const gmx_mtop_t*   mtop,
                                  const matrix        box,
                                  bool                makeGpuPairList,
                                  const gmx::CpuInfo& cpuinfo)
{
    // We checked the cut-offs in grompp, but double-check here.
    // We have PME+LJcutoff kernels for rcoulomb>rvdw.
    if (EEL_PME_EWALD(ir->coulombtype) && ir->vdwtype == eelCUT)
    {
        GMX_RELEASE_ASSERT(ir->rcoulomb >= ir->rvdw,
                           "With Verlet lists and PME we should have rcoulomb>=rvdw");
    }
    else
    {
        GMX_RELEASE_ASSERT(ir->rcoulomb == ir->rvdw,
                           "With Verlet lists and no PME rcoulomb and rvdw should be identical");
    }
    /* For NVE simulations, we will retain the initial list buffer */
    if (EI_DYNAMICS(ir->eI) && ir->verletbuf_tol > 0 && !(EI_MD(ir->eI) && ir->etc == etcNO))
    {
        /* Update the Verlet buffer size for the current run setup */

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        ListSetupType listType =
                (makeGpuPairList ? ListSetupType::Gpu : ListSetupType::CpuSimdWhenSupported);
        VerletbufListSetup listSetup = verletbufGetSafeListSetup(listType);

        const real rlist_new =
                calcVerletBufferSize(*mtop, det(box), *ir, ir->nstlist, ir->nstlist - 1, -1, listSetup);

        if (rlist_new != ir->rlist)
        {
            if (fplog != nullptr)
            {
                fprintf(fplog,
                        "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist, rlist_new, listSetup.cluster_size_i, listSetup.cluster_size_j);
            }
            ir->rlist = rlist_new;
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
static void override_nsteps_cmdline(const gmx::MDLogger& mdlog, int64_t nsteps_cmdline, t_inputrec* ir)
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
            sprintf(sbuf_msg,
                    "Overriding nsteps with value passed on the command line: %s steps, %.3g ps",
                    gmx_step_str(nsteps_cmdline, sbuf_steps), fabs(nsteps_cmdline * ir->delta_t));
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
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %" PRId64, nsteps_cmdline);
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
static bool gpuAccelerationOfNonbondedIsUseful(const MDLogger& mdlog, const t_inputrec& ir, bool issueWarning)
{
    bool        gpuIsUseful = true;
    std::string warning;

    if (ir.opts.ngener - ir.nwall > 1)
    {
        /* The GPU code does not support more than one energy group.
         * If the user requested GPUs explicitly, a fatal error is given later.
         */
        gpuIsUseful = false;
        warning =
                "Multiple energy groups is not implemented for GPUs, falling back to the CPU. "
                "For better performance, run on the GPU without energy groups and then do "
                "gmx mdrun -rerun option on the trajectory with an energy group .tpr file.";
    }

    if (EI_TPI(ir.eI))
    {
        gpuIsUseful = false;
        warning     = "TPI is not implemented for GPUs.";
    }

    if (!gpuIsUseful && issueWarning)
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText(warning);
    }

    return gpuIsUseful;
}

//! Initializes the logger for mdrun.
static gmx::LoggerOwner buildLogger(FILE* fplog, const bool isSimulationMasterRank)
{
    gmx::LoggerBuilder builder;
    if (fplog != nullptr)
    {
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fplog);
    }
    if (isSimulationMasterRank)
    {
        builder.addTargetStream(gmx::MDLogger::LogLevel::Warning, &gmx::TextOutputFile::standardError());
    }
    return builder.build();
}

//! Make a TaskTarget from an mdrun argument string.
static TaskTarget findTaskTarget(const char* optionString)
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

//! Finish run, aggregate data to print performance info.
static void finish_run(FILE*                     fplog,
                       const gmx::MDLogger&      mdlog,
                       const t_commrec*          cr,
                       const t_inputrec*         inputrec,
                       t_nrnb                    nrnb[],
                       gmx_wallcycle_t           wcycle,
                       gmx_walltime_accounting_t walltime_accounting,
                       nonbonded_verlet_t*       nbv,
                       const gmx_pme_t*          pme,
                       gmx_bool                  bWriteStat)
{
    double delta_t = 0;
    double nbfs = 0, mflop = 0;
    double elapsed_time, elapsed_time_over_all_ranks, elapsed_time_over_all_threads,
            elapsed_time_over_all_threads_over_all_ranks;
    /* Control whether it is valid to print a report. Only the
       simulation master may print, but it should not do so if the run
       terminated e.g. before a scheduled reset step. This is
       complicated by the fact that PME ranks are unaware of the
       reason why they were sent a pmerecvqxFINISH. To avoid
       communication deadlocks, we always do the communication for the
       report, even if we've decided not to write the report, because
       how long it takes to finish the run is not important when we've
       decided not to report on the simulation performance.

       Further, we only report performance for dynamical integrators,
       because those are the only ones for which we plan to
       consider doing any optimizations. */
    bool printReport = EI_DYNAMICS(inputrec->eI) && SIMMASTER(cr);

    if (printReport && !walltime_accounting_get_valid_finish(walltime_accounting))
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("Simulation ended prematurely, no performance report will be written.");
        printReport = false;
    }

    t_nrnb*                 nrnb_tot;
    std::unique_ptr<t_nrnb> nrnbTotalStorage;
    if (cr->nnodes > 1)
    {
        nrnbTotalStorage = std::make_unique<t_nrnb>();
        nrnb_tot         = nrnbTotalStorage.get();
#if GMX_MPI
        MPI_Allreduce(nrnb->n, nrnb_tot->n, eNRNB, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mysim);
#endif
    }
    else
    {
        nrnb_tot = nrnb;
    }

    elapsed_time = walltime_accounting_get_time_since_reset(walltime_accounting);
    elapsed_time_over_all_threads =
            walltime_accounting_get_time_since_reset_over_all_threads(walltime_accounting);
    if (cr->nnodes > 1)
    {
#if GMX_MPI
        /* reduce elapsed_time over all MPI ranks in the current simulation */
        MPI_Allreduce(&elapsed_time, &elapsed_time_over_all_ranks, 1, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mysim);
        elapsed_time_over_all_ranks /= cr->nnodes;
        /* Reduce elapsed_time_over_all_threads over all MPI ranks in the
         * current simulation. */
        MPI_Allreduce(&elapsed_time_over_all_threads, &elapsed_time_over_all_threads_over_all_ranks,
                      1, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mysim);
#endif
    }
    else
    {
        elapsed_time_over_all_ranks                  = elapsed_time;
        elapsed_time_over_all_threads_over_all_ranks = elapsed_time_over_all_threads;
    }

    if (printReport)
    {
        print_flop(fplog, nrnb_tot, &nbfs, &mflop);
    }

    if (thisRankHasDuty(cr, DUTY_PP) && DOMAINDECOMP(cr))
    {
        print_dd_statistics(cr, inputrec, fplog);
    }

    /* TODO Move the responsibility for any scaling by thread counts
     * to the code that handled the thread region, so that there's a
     * mechanism to keep cycle counting working during the transition
     * to task parallelism. */
    int nthreads_pp  = gmx_omp_nthreads_get(emntNonbonded);
    int nthreads_pme = gmx_omp_nthreads_get(emntPME);
    wallcycle_scale_by_num_threads(wcycle, thisRankHasDuty(cr, DUTY_PME) && !thisRankHasDuty(cr, DUTY_PP),
                                   nthreads_pp, nthreads_pme);
    auto cycle_sum(wallcycle_sum(cr, wcycle));

    if (printReport)
    {
        auto nbnxn_gpu_timings =
                (nbv != nullptr && nbv->useGpu()) ? Nbnxm::gpu_get_timings(nbv->gpu_nbv) : nullptr;
        gmx_wallclock_gpu_pme_t pme_gpu_timings = {};

        if (pme_gpu_task_enabled(pme))
        {
            pme_gpu_get_timings(pme, &pme_gpu_timings);
        }
        wallcycle_print(fplog, mdlog, cr->nnodes, cr->npmenodes, nthreads_pp, nthreads_pme,
                        elapsed_time_over_all_ranks, wcycle, cycle_sum, nbnxn_gpu_timings,
                        &pme_gpu_timings);

        if (EI_DYNAMICS(inputrec->eI))
        {
            delta_t = inputrec->delta_t;
        }

        if (fplog)
        {
            print_perf(fplog, elapsed_time_over_all_threads_over_all_ranks, elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done_since_reset(walltime_accounting),
                       delta_t, nbfs, mflop);
        }
        if (bWriteStat)
        {
            print_perf(stderr, elapsed_time_over_all_threads_over_all_ranks, elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done_since_reset(walltime_accounting),
                       delta_t, nbfs, mflop);
        }
    }
}

int Mdrunner::mdrunner()
{
    matrix                    box;
    t_forcerec*               fr               = nullptr;
    t_fcdata*                 fcd              = nullptr;
    real                      ewaldcoeff_q     = 0;
    real                      ewaldcoeff_lj    = 0;
    int                       nChargePerturbed = -1, nTypePerturbed = 0;
    gmx_wallcycle_t           wcycle;
    gmx_walltime_accounting_t walltime_accounting = nullptr;
    gmx_membed_t*             membed              = nullptr;
    gmx_hw_info_t*            hwinfo              = nullptr;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    gmx_mtop_t mtop;

    /* TODO: inputrec should tell us whether we use an algorithm, not a file option */
    const bool doEssentialDynamics = opt2bSet("-ei", filenames.size(), filenames.data());
    const bool doMembed            = opt2bSet("-membed", filenames.size(), filenames.data());
    const bool doRerun             = mdrunOptions.rerun;

    // Handle task-assignment related user options.
    EmulateGpuNonbonded emulateGpuNonbonded =
            (getenv("GMX_EMULATE_GPU") != nullptr ? EmulateGpuNonbonded::Yes : EmulateGpuNonbonded::No);

    std::vector<int> userGpuTaskAssignment;
    try
    {
        userGpuTaskAssignment = parseUserTaskAssignmentString(hw_opt.userGpuTaskAssignment);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    auto nonbondedTarget = findTaskTarget(nbpu_opt);
    auto pmeTarget       = findTaskTarget(pme_opt);
    auto pmeFftTarget    = findTaskTarget(pme_fft_opt);
    auto bondedTarget    = findTaskTarget(bonded_opt);
    auto updateTarget    = findTaskTarget(update_opt);

    FILE* fplog = nullptr;
    // If we are appending, we don't write log output because we need
    // to check that the old log file matches what the checkpoint file
    // expects. Otherwise, we should start to write log output now if
    // there is a file ready for it.
    if (logFileHandle != nullptr && startingBehavior != StartingBehavior::RestartWithAppending)
    {
        fplog = gmx_fio_getfp(logFileHandle);
    }
    const bool       isSimulationMasterRank = findIsSimulationMasterRank(ms, communicator);
    gmx::LoggerOwner logOwner(buildLogger(fplog, isSimulationMasterRank));
    gmx::MDLogger    mdlog(logOwner.logger());

    // TODO The thread-MPI master rank makes a working
    // PhysicalNodeCommunicator here, but it gets rebuilt by all ranks
    // after the threads have been launched. This works because no use
    // is made of that communicator until after the execution paths
    // have rejoined. But it is likely that we can improve the way
    // this is expressed, e.g. by expressly running detection only the
    // master rank for thread-MPI, rather than relying on the mutex
    // and reference count.
    PhysicalNodeCommunicator physicalNodeComm(communicator, gmx_physicalnode_id_hash());
    hwinfo = gmx_detect_hardware(mdlog, physicalNodeComm);

    gmx_print_detected_hardware(fplog, isSimulationMasterRank && isMasterSim(ms), mdlog, hwinfo);

    std::vector<int> gpuIdsToUse = makeGpuIdsToUse(hwinfo->gpu_info, hw_opt.gpuIdsAvailable);

    // Print citation requests after all software/hardware printing
    pleaseCiteGromacs(fplog);

    // TODO Replace this by unique_ptr once t_inputrec is C++
    t_inputrec               inputrecInstance;
    t_inputrec*              inputrec = nullptr;
    std::unique_ptr<t_state> globalState;

    auto partialDeserializedTpr = std::make_unique<PartialDeserializedTprFile>();

    if (isSimulationMasterRank)
    {
        /* Only the master rank has the global state */
        globalState = std::make_unique<t_state>();

        /* Read (nearly) all data required for the simulation
         * and keep the partly serialized tpr contents to send to other ranks later
         */
        *partialDeserializedTpr = read_tpx_state(ftp2fn(efTPR, filenames.size(), filenames.data()),
                                                 &inputrecInstance, globalState.get(), &mtop);
        inputrec                = &inputrecInstance;
    }

    /* Check and update the hardware options for internal consistency */
    checkAndUpdateHardwareOptions(mdlog, &hw_opt, isSimulationMasterRank, domdecOptions.numPmeRanks,
                                  inputrec);

    if (GMX_THREAD_MPI && isSimulationMasterRank)
    {
        bool useGpuForNonbonded = false;
        bool useGpuForPme       = false;
        try
        {
            GMX_RELEASE_ASSERT(inputrec != nullptr, "Keep the compiler happy");

            // If the user specified the number of ranks, then we must
            // respect that, but in default mode, we need to allow for
            // the number of GPUs to choose the number of ranks.
            auto canUseGpuForNonbonded = buildSupportsNonbondedOnGpu(nullptr);
            useGpuForNonbonded         = decideWhetherToUseGpusForNonbondedWithThreadMpi(
                    nonbondedTarget, gpuIdsToUse, userGpuTaskAssignment, emulateGpuNonbonded,
                    canUseGpuForNonbonded,
                    gpuAccelerationOfNonbondedIsUseful(mdlog, *inputrec, GMX_THREAD_MPI),
                    hw_opt.nthreads_tmpi);
            useGpuForPme = decideWhetherToUseGpusForPmeWithThreadMpi(
                    useGpuForNonbonded, pmeTarget, gpuIdsToUse, userGpuTaskAssignment, *hwinfo,
                    *inputrec, mtop, hw_opt.nthreads_tmpi, domdecOptions.numPmeRanks);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

        /* Determine how many thread-MPI ranks to start.
         *
         * TODO Over-writing the user-supplied value here does
         * prevent any possible subsequent checks from working
         * correctly. */
        hw_opt.nthreads_tmpi = get_nthreads_mpi(hwinfo, &hw_opt, gpuIdsToUse, useGpuForNonbonded,
                                                useGpuForPme, inputrec, &mtop, mdlog, doMembed);

        // Now start the threads for thread MPI.
        spawnThreads(hw_opt.nthreads_tmpi);
        // The spawned threads enter mdrunner() and execution of
        // master and spawned threads joins at the end of this block.
        physicalNodeComm = PhysicalNodeCommunicator(communicator, gmx_physicalnode_id_hash());
    }

    GMX_RELEASE_ASSERT(communicator == MPI_COMM_WORLD, "Must have valid world communicator");
    CommrecHandle crHandle = init_commrec(communicator, ms);
    t_commrec*    cr       = crHandle.get();
    GMX_RELEASE_ASSERT(cr != nullptr, "Must have valid commrec");

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        if (!isSimulationMasterRank)
        {
            inputrec = &inputrecInstance;
        }
        init_parallel(cr->mpi_comm_mygroup, MASTER(cr), inputrec, &mtop, partialDeserializedTpr.get());
    }
    GMX_RELEASE_ASSERT(inputrec != nullptr, "All ranks should have a valid inputrec now");
    partialDeserializedTpr.reset(nullptr);

    // Now the number of ranks is known to all ranks, and each knows
    // the inputrec read by the master rank. The ranks can now all run
    // the task-deciding functions and will agree on the result
    // without needing to communicate.
    const bool useDomainDecomposition = (PAR(cr) && !(EI_TPI(inputrec->eI) || inputrec->eI == eiNM));

    // Note that these variables describe only their own node.
    //
    // Note that when bonded interactions run on a GPU they always run
    // alongside a nonbonded task, so do not influence task assignment
    // even though they affect the force calculation workload.
    bool useGpuForNonbonded = false;
    bool useGpuForPme       = false;
    bool useGpuForBonded    = false;
    bool useGpuForUpdate    = false;
    bool gpusWereDetected   = hwinfo->ngpu_compatible_tot > 0;
    try
    {
        // It's possible that there are different numbers of GPUs on
        // different nodes, which is the user's responsibility to
        // handle. If unsuitable, we will notice that during task
        // assignment.
        auto canUseGpuForNonbonded = buildSupportsNonbondedOnGpu(nullptr);
        useGpuForNonbonded         = decideWhetherToUseGpusForNonbonded(
                nonbondedTarget, userGpuTaskAssignment, emulateGpuNonbonded, canUseGpuForNonbonded,
                gpuAccelerationOfNonbondedIsUseful(mdlog, *inputrec, !GMX_THREAD_MPI), gpusWereDetected);
        useGpuForPme = decideWhetherToUseGpusForPme(
                useGpuForNonbonded, pmeTarget, userGpuTaskAssignment, *hwinfo, *inputrec, mtop,
                cr->nnodes, domdecOptions.numPmeRanks, gpusWereDetected);
        auto canUseGpuForBonded = buildSupportsGpuBondeds(nullptr)
                                  && inputSupportsGpuBondeds(*inputrec, mtop, nullptr);
        useGpuForBonded = decideWhetherToUseGpusForBonded(
                useGpuForNonbonded, useGpuForPme, bondedTarget, canUseGpuForBonded,
                EVDW_PME(inputrec->vdwtype), EEL_PME_EWALD(inputrec->coulombtype),
                domdecOptions.numPmeRanks, gpusWereDetected);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    const PmeRunMode pmeRunMode = determinePmeRunMode(useGpuForPme, pmeFftTarget, *inputrec);

    // Initialize development feature flags that enabled by environment variable
    // and report those features that are enabled.
    const DevelopmentFeatureFlags devFlags =
            manageDevelopmentFeatures(mdlog, useGpuForNonbonded, pmeRunMode);

    const bool useModularSimulator = checkUseModularSimulator(
            false, inputrec, doRerun, mtop, ms, replExParams, nullptr, doEssentialDynamics, doMembed);

    // Build restraints.
    // TODO: hide restraint implementation details from Mdrunner.
    // There is nothing unique about restraints at this point as far as the
    // Mdrunner is concerned. The Mdrunner should just be getting a sequence of
    // factory functions from the SimulationContext on which to call mdModules_->add().
    // TODO: capture all restraints into a single RestraintModule, passed to the runner builder.
    for (auto&& restraint : restraintManager_->getRestraints())
    {
        auto module = RestraintMDModule::create(restraint, restraint->sites());
        mdModules_->add(std::move(module));
    }

    // TODO: Error handling
    mdModules_->assignOptionsToModules(*inputrec->params, nullptr);
    // now that the MdModules know their options, they know which callbacks to sign up to
    mdModules_->subscribeToSimulationSetupNotifications();
    const auto& mdModulesNotifier = mdModules_->notifier().simulationSetupNotifications_;

    if (inputrec->internalParameters != nullptr)
    {
        mdModulesNotifier.notify(*inputrec->internalParameters);
    }

    if (fplog != nullptr)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec, FALSE);
        fprintf(fplog, "\n");
    }

    if (SIMMASTER(cr))
    {
        /* In rerun, set velocities to zero if present */
        if (doRerun && ((globalState->flags & (1 << estV)) != 0))
        {
            // rerun does not use velocities
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText(
                            "Rerun trajectory contains velocities. Rerun does only evaluate "
                            "potential energy and forces. The velocities will be ignored.");
            for (int i = 0; i < globalState->natoms; i++)
            {
                clear_rvec(globalState->v[i]);
            }
            globalState->flags &= ~(1 << estV);
        }

        /* now make sure the state is initialized and propagated */
        set_state_entries(globalState.get(), inputrec, useModularSimulator);
    }

    /* NM and TPI parallelize over force/energy calculations, not atoms,
     * so we need to initialize and broadcast the global state.
     */
    if (inputrec->eI == eiNM || inputrec->eI == eiTPI)
    {
        if (!MASTER(cr))
        {
            globalState = std::make_unique<t_state>();
        }
        broadcastStateWithoutDynamics(cr->mpi_comm_mygroup, DOMAINDECOMP(cr), PAR(cr), globalState.get());
    }

    /* A parallel command line option consistency check that we can
       only do after any threads have started. */
    if (!PAR(cr)
        && (domdecOptions.numCells[XX] > 1 || domdecOptions.numCells[YY] > 1
            || domdecOptions.numCells[ZZ] > 1 || domdecOptions.numPmeRanks > 0))
    {
        gmx_fatal(FARGS,
                  "The -dd or -npme option request a parallel simulation, "
#if !GMX_MPI
                  "but %s was compiled without threads or MPI enabled",
                  output_env_get_program_display_name(oenv));
#elif GMX_THREAD_MPI
                  "but the number of MPI-threads (option -ntmpi) is not set or is 1");
#else
                  "but %s was not started through mpirun/mpiexec or only one rank was requested "
                  "through mpirun/mpiexec",
                  output_env_get_program_display_name(oenv));
#endif
    }

    if (doRerun && (EI_ENERGY_MINIMIZATION(inputrec->eI) || eiNM == inputrec->eI))
    {
        gmx_fatal(FARGS,
                  "The .mdp file specified an energy mininization or normal mode algorithm, and "
                  "these are not compatible with mdrun -rerun");
    }

    if (!(EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype)))
    {
        if (domdecOptions.numPmeRanks > 0)
        {
            gmx_fatal_collective(FARGS, cr->mpi_comm_mysim, MASTER(cr),
                                 "PME-only ranks are requested, but the system does not use PME "
                                 "for electrostatics or LJ");
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
            GMX_RELEASE_ASSERT(domdecOptions.numPmeRanks <= 1,
                               "PME GPU decomposition is not supported");
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
    init_disres(fplog, &mtop, inputrec, DisResRunMode::MDRun, MASTER(cr) ? DDRole::Master : DDRole::Agent,
                PAR(cr) ? NumRanks::Multiple : NumRanks::Single, cr->mpi_comm_mysim, ms, fcd,
                globalState.get(), replExParams.exchangeInterval > 0);

    init_orires(fplog, &mtop, inputrec, cr, ms, globalState.get(), &(fcd->orires));

    auto deform = prepareBoxDeformation(globalState->box, MASTER(cr) ? DDRole::Master : DDRole::Agent,
                                        PAR(cr) ? NumRanks::Multiple : NumRanks::Single,
                                        cr->mpi_comm_mygroup, *inputrec);

    ObservablesHistory observablesHistory = {};

    if (startingBehavior != StartingBehavior::NewSimulation)
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        if (mdrunOptions.numStepsCommandline > -2)
        {
            /* Temporarily set the number of steps to unlimited to avoid
             * triggering the nsteps check in load_checkpoint().
             * This hack will go away soon when the -nsteps option is removed.
             */
            inputrec->nsteps = -1;
        }

        load_checkpoint(opt2fn_master("-cpi", filenames.size(), filenames.data(), cr),
                        logFileHandle, cr, domdecOptions.numCells, inputrec, globalState.get(),
                        &observablesHistory, mdrunOptions.reproducible, mdModules_->notifier());

        if (startingBehavior == StartingBehavior::RestartWithAppending && logFileHandle)
        {
            // Now we can start normal logging to the truncated log file.
            fplog = gmx_fio_getfp(logFileHandle);
            prepareLogAppending(fplog);
            logOwner = buildLogger(fplog, MASTER(cr));
            mdlog    = logOwner.logger();
        }
    }

    if (mdrunOptions.numStepsCommandline > -2)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -nsteps functionality is deprecated, and may be removed in a future "
                        "version. "
                        "Consider using gmx convert-tpr -nsteps or changing the appropriate .mdp "
                        "file field.");
    }
    /* override nsteps with value set on the commandline */
    override_nsteps_cmdline(mdlog, mdrunOptions.numStepsCommandline, inputrec);

    if (SIMMASTER(cr))
    {
        copy_mat(globalState->box, box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box), box, cr->mpi_comm_mygroup);
    }

    if (inputrec->cutoff_scheme != ecutsVERLET)
    {
        gmx_fatal(FARGS,
                  "This group-scheme .tpr file can no longer be run by mdrun. Please update to the "
                  "Verlet scheme, or use an earlier version of GROMACS if necessary.");
    }
    /* Update rlist and nstlist. */
    prepare_verlet_scheme(fplog, cr, inputrec, nstlist_cmdline, &mtop, box,
                          useGpuForNonbonded || (emulateGpuNonbonded == EmulateGpuNonbonded::Yes),
                          *hwinfo->cpuInfo);

    const bool prefer1DAnd1PulseDD = (devFlags.enableGpuHaloExchange && useGpuForNonbonded);
    // This builder is necessary while we have multi-part construction
    // of DD. Before DD is constructed, we use the existence of
    // the builder object to indicate that further construction of DD
    // is needed.
    std::unique_ptr<DomainDecompositionBuilder> ddBuilder;
    if (useDomainDecomposition)
    {
        ddBuilder = std::make_unique<DomainDecompositionBuilder>(
                mdlog, cr, domdecOptions, mdrunOptions, prefer1DAnd1PulseDD, mtop, *inputrec, box,
                positionsFromStatePointer(globalState.get()));
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty      = (DUTY_PP | DUTY_PME);

        if (inputrec->pbcType == PbcType::Screw)
        {
            gmx_fatal(FARGS, "pbc=screw is only implemented with domain decomposition");
        }
    }

    // Produce the task assignment for this rank - done after DD is constructed
    GpuTaskAssignments gpuTaskAssignments = GpuTaskAssignmentsBuilder::build(
            gpuIdsToUse, userGpuTaskAssignment, *hwinfo, communicator, physicalNodeComm,
            nonbondedTarget, pmeTarget, bondedTarget, updateTarget, useGpuForNonbonded,
            useGpuForPme, thisRankHasDuty(cr, DUTY_PP),
            // TODO cr->duty & DUTY_PME should imply that a PME
            // algorithm is active, but currently does not.
            EEL_PME(inputrec->coulombtype) && thisRankHasDuty(cr, DUTY_PME));

    // Get the device handles for the modules, nullptr when no task is assigned.
    int                deviceId   = -1;
    DeviceInformation* deviceInfo = gpuTaskAssignments.initDevice(&deviceId);

    // timing enabling - TODO put this in gpu_utils (even though generally this is just option handling?)
    bool useTiming = true;
    if (GMX_GPU == GMX_GPU_CUDA)
    {
        /* WARNING: CUDA timings are incorrect with multiple streams.
         *          This is the main reason why they are disabled by default.
         */
        // TODO: Consider turning on by default when we can detect nr of streams.
        useTiming = (getenv("GMX_ENABLE_GPU_TIMING") != nullptr);
    }
    else if (GMX_GPU == GMX_GPU_OPENCL)
    {
        useTiming = (getenv("GMX_DISABLE_GPU_TIMING") == nullptr);
    }

    // TODO Currently this is always built, yet DD partition code
    // checks if it is built before using it. Probably it should
    // become an MDModule that is made only when another module
    // requires it (e.g. pull, CompEl, density fitting), so that we
    // don't update the local atom sets unilaterally every step.
    LocalAtomSetManager atomSets;
    if (ddBuilder)
    {
        // TODO Pass the GPU streams to ddBuilder to use in buffer
        // transfers (e.g. halo exchange)
        cr->dd = ddBuilder->build(&atomSets);
        // The builder's job is done, so destruct it
        ddBuilder.reset(nullptr);
        // Note that local state still does not exist yet.
    }

    // The GPU update is decided here because we need to know whether the constraints or
    // SETTLEs can span accross the domain borders (i.e. whether or not update groups are
    // defined). This is only known after DD is initialized, hence decision on using GPU
    // update is done so late.
    try
    {
        const bool useUpdateGroups = cr->dd ? ddUsesUpdateGroups(*cr->dd) : false;

        useGpuForUpdate = decideWhetherToUseGpuForUpdate(
                useDomainDecomposition, useUpdateGroups, pmeRunMode, domdecOptions.numPmeRanks > 0,
                useGpuForNonbonded, updateTarget, gpusWereDetected, *inputrec, mtop,
                doEssentialDynamics, gmx_mtop_ftype_count(mtop, F_ORIRES) > 0,
                replExParams.exchangeInterval > 0, doRerun, devFlags, mdlog);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    const bool printHostName = (cr->nnodes > 1);
    gpuTaskAssignments.reportGpuUsage(mdlog, printHostName, useGpuForBonded, pmeRunMode, useGpuForUpdate);

    std::unique_ptr<DeviceStreamManager> deviceStreamManager = nullptr;

    if (deviceInfo != nullptr)
    {
        if (DOMAINDECOMP(cr) && thisRankHasDuty(cr, DUTY_PP))
        {
            dd_setup_dlb_resource_sharing(cr, deviceId);
        }
        deviceStreamManager = std::make_unique<DeviceStreamManager>(
                *deviceInfo, useGpuForPme, useGpuForNonbonded, havePPDomainDecomposition(cr),
                useGpuForUpdate, useTiming);
    }

    // If the user chose a task assignment, give them some hints
    // where appropriate.
    if (!userGpuTaskAssignment.empty())
    {
        gpuTaskAssignments.logPerformanceHints(mdlog, ssize(gpuIdsToUse));
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
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "This is simulation %d out of %d running as a composite GROMACS\n"
                        "multi-simulation job. Setup for this simulation:\n",
                        ms->sim, ms->nsim);
    }
    GMX_LOG(mdlog.warning)
            .appendTextFormatted("Using %d MPI %s\n", cr->nnodes,
#    if GMX_THREAD_MPI
                                 cr->nnodes == 1 ? "thread" : "threads"
#    else
                                 cr->nnodes == 1 ? "process" : "processes"
#    endif
            );
    fflush(stderr);
#endif

    // If mdrun -pin auto honors any affinity setting that already
    // exists. If so, it is nice to provide feedback about whether
    // that existing affinity setting was from OpenMP or something
    // else, so we run this code both before and after we initialize
    // the OpenMP support.
    gmx_check_thread_affinity_set(mdlog, &hw_opt, hwinfo->nthreads_hw_avail, FALSE);
    /* Check and update the number of OpenMP threads requested */
    checkAndUpdateRequestedNumOpenmpThreads(&hw_opt, *hwinfo, cr, ms, physicalNodeComm.size_,
                                            pmeRunMode, mtop, *inputrec);

    gmx_omp_nthreads_init(mdlog, cr, hwinfo->nthreads_hw_avail, physicalNodeComm.size_,
                          hw_opt.nthreads_omp, hw_opt.nthreads_omp_pme, !thisRankHasDuty(cr, DUTY_PP));

    // Enable FP exception detection, but not in
    // Release mode and not for compilers with known buggy FP
    // exception support (clang with any optimization) or suspected
    // buggy FP exception support (gcc 7.* with optimization).
#if !defined NDEBUG                                                                         \
        && !((defined __clang__ || (defined(__GNUC__) && !defined(__ICC) && __GNUC__ == 7)) \
             && defined __OPTIMIZE__)
    const bool bEnableFPE = true;
#else
    const bool bEnableFPE = false;
#endif
    // FIXME - reconcile with gmx_feenableexcept() call from CommandLineModuleManager::run()
    if (bEnableFPE)
    {
        gmx_feenableexcept();
    }

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, gpuTaskAssignments.thisRankHasAnyGpuTask(),
                                       mdrunOptions.ntompOptionIsSet, cr, mdlog);

    /* getting number of PP/PME threads on this MPI / tMPI rank.
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    const int numThreadsOnThisRank = thisRankHasDuty(cr, DUTY_PP) ? gmx_omp_nthreads_get(emntNonbonded)
                                                                  : gmx_omp_nthreads_get(emntPME);
    checkHardwareOversubscription(numThreadsOnThisRank, cr->nodeid, *hwinfo->hardwareTopology,
                                  physicalNodeComm, mdlog);

    // Enable Peer access between GPUs where available
    // Only for DD, only master PP rank needs to perform setup, and only if thread MPI plus
    // any of the GPU communication features are active.
    if (DOMAINDECOMP(cr) && MASTER(cr) && thisRankHasDuty(cr, DUTY_PP) && GMX_THREAD_MPI
        && (devFlags.enableGpuHaloExchange || devFlags.enableGpuPmePPComm))
    {
        setupGpuDevicePeerAccess(gpuIdsToUse, mdlog);
    }

    if (hw_opt.threadAffinity != ThreadAffinity::Off)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(mdlog, &hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        int numThreadsOnThisNode, intraNodeThreadOffset;
        analyzeThreadsOnThisNode(physicalNodeComm, numThreadsOnThisRank, &numThreadsOnThisNode,
                                 &intraNodeThreadOffset);

        /* Set the CPU affinity */
        gmx_set_thread_affinity(mdlog, cr, &hw_opt, *hwinfo->hardwareTopology, numThreadsOnThisRank,
                                numThreadsOnThisNode, intraNodeThreadOffset, nullptr);
    }

    if (mdrunOptions.timingOptions.resetStep > -1)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -resetstep functionality is deprecated, and may be removed in a "
                        "future version.");
    }
    wcycle = wallcycle_init(fplog, mdrunOptions.timingOptions.resetStep, cr);

    if (PAR(cr))
    {
        /* Master synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        int64_t reset_counters = wcycle_get_reset_counters(wcycle);
        gmx_bcast(sizeof(reset_counters), &reset_counters, cr->mpi_comm_mysim);
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
        membed = init_membed(fplog, filenames.size(), filenames.data(), &mtop, inputrec,
                             globalState.get(), cr, &mdrunOptions.checkpointOptions.period);
    }

    const bool                   thisRankHasPmeGpuTask = gpuTaskAssignments.thisRankHasPmeGpuTask();
    std::unique_ptr<MDAtoms>     mdAtoms;
    std::unique_ptr<gmx_vsite_t> vsite;
    std::unique_ptr<GpuBonded>   gpuBonded;

    t_nrnb nrnb;
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        mdModulesNotifier.notify(*cr);
        mdModulesNotifier.notify(&atomSets);
        mdModulesNotifier.notify(inputrec->pbcType);
        mdModulesNotifier.notify(SimulationTimeStep{ inputrec->delta_t });
        /* Initiate forcerecord */
        fr                 = new t_forcerec;
        fr->forceProviders = mdModules_->initForceProviders();
        init_forcerec(fplog, mdlog, fr, fcd, inputrec, &mtop, cr, box,
                      opt2fn("-table", filenames.size(), filenames.data()),
                      opt2fn("-tablep", filenames.size(), filenames.data()),
                      opt2fns("-tableb", filenames.size(), filenames.data()), pforce);

        // Save a handle to device stream manager to use elsewhere in the code
        // TODO: Forcerec is not a correct place to store it.
        fr->deviceStreamManager = deviceStreamManager.get();

        if (devFlags.enableGpuPmePPComm && !thisRankHasDuty(cr, DUTY_PME))
        {
            GMX_RELEASE_ASSERT(
                    deviceStreamManager != nullptr,
                    "GPU device stream manager should be valid in order to use PME-PP direct "
                    "communications.");
            GMX_RELEASE_ASSERT(
                    deviceStreamManager->streamIsValid(DeviceStreamType::PmePpTransfer),
                    "GPU PP-PME stream should be valid in order to use GPU PME-PP direct "
                    "communications.");
            fr->pmePpCommGpu = std::make_unique<gmx::PmePpCommGpu>(
                    cr->mpi_comm_mysim, cr->dd->pme_nodeid, deviceStreamManager->context(),
                    deviceStreamManager->stream(DeviceStreamType::PmePpTransfer));
        }

        fr->nbv = Nbnxm::init_nb_verlet(mdlog, inputrec, fr, cr, *hwinfo, useGpuForNonbonded,
                                        deviceStreamManager.get(), &mtop, box, wcycle);
        // TODO: Move the logic below to a GPU bonded builder
        if (useGpuForBonded)
        {
            GMX_RELEASE_ASSERT(deviceStreamManager != nullptr,
                               "GPU device stream manager should be valid in order to use GPU "
                               "version of bonded forces.");
            gpuBonded = std::make_unique<GpuBonded>(
                    mtop.ffparams, fr->ic->epsfac * fr->fudgeQQ, deviceStreamManager->context(),
                    deviceStreamManager->bondedStream(havePPDomainDecomposition(cr)), wcycle);
            fr->gpuBonded = gpuBonded.get();
        }

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
            // FIXME: MD and EM separately set up the local state - this should happen in the same
            // function, which should also perform the pinning.
            changePinningPolicy(&globalState->x, pme_get_pinning_policy());
        }

        /* Initialize the virtual site communication */
        vsite = initVsite(mtop, cr);

        calc_shifts(box, fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MASTER(cr)
            && !(inputrec->pbcType != PbcType::No && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->pbcType != PbcType::No)
            {
                do_pbc_first_mtop(fplog, inputrec->pbcType, box, &mtop, globalState->x.rvec_array());
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

        GMX_ASSERT(globalState == nullptr,
                   "We don't need the state on a PME only rank and expect it to be unitialized");

        ewaldcoeff_q  = calc_ewaldcoeff_q(inputrec->rcoulomb, inputrec->ewald_rtol);
        ewaldcoeff_lj = calc_ewaldcoeff_lj(inputrec->rvdw, inputrec->ewald_rtol_lj);
    }

    gmx_pme_t* sepPmeData = nullptr;
    // This reference hides the fact that PME data is owned by runner on PME-only ranks and by forcerec on other ranks
    GMX_ASSERT(thisRankHasDuty(cr, DUTY_PP) == (fr != nullptr),
               "Double-checking that only PME-only ranks have no forcerec");
    gmx_pme_t*& pmedata = fr ? fr->pmedata : sepPmeData;

    // TODO should live in ewald module once its testing is improved
    //
    // Later, this program could contain kernels that might be later
    // re-used as auto-tuning progresses, or subsequent simulations
    // are invoked.
    PmeGpuProgramStorage pmeGpuProgram;
    if (thisRankHasPmeGpuTask)
    {
        GMX_RELEASE_ASSERT(
                (deviceStreamManager != nullptr),
                "GPU device stream manager should be initialized in order to use GPU for PME.");
        GMX_RELEASE_ASSERT((deviceInfo != nullptr),
                           "GPU device should be initialized in order to use GPU for PME.");
        pmeGpuProgram = buildPmeGpuProgram(deviceStreamManager->context());
    }

    /* Initiate PME if necessary,
     * either on all nodes or on dedicated PME nodes only. */
    if (EEL_PME(inputrec->coulombtype) || EVDW_PME(inputrec->vdwtype))
    {
        if (mdAtoms && mdAtoms->mdatoms())
        {
            nChargePerturbed = mdAtoms->mdatoms()->nChargePerturbed;
            if (EVDW_PME(inputrec->vdwtype))
            {
                nTypePerturbed = mdAtoms->mdatoms()->nTypePerturbed;
            }
        }
        if (cr->npmenodes > 0)
        {
            /* The PME only nodes need to know nChargePerturbed(FEP on Q) and nTypePerturbed(FEP on LJ)*/
            gmx_bcast(sizeof(nChargePerturbed), &nChargePerturbed, cr->mpi_comm_mysim);
            gmx_bcast(sizeof(nTypePerturbed), &nTypePerturbed, cr->mpi_comm_mysim);
        }

        if (thisRankHasDuty(cr, DUTY_PME))
        {
            try
            {
                // TODO: This should be in the builder.
                GMX_RELEASE_ASSERT(!useGpuForPme || (deviceStreamManager != nullptr),
                                   "Device stream manager should be valid in order to use GPU "
                                   "version of PME.");
                GMX_RELEASE_ASSERT(
                        !useGpuForPme || deviceStreamManager->streamIsValid(DeviceStreamType::Pme),
                        "GPU PME stream should be valid in order to use GPU version of PME.");

                const DeviceContext* deviceContext =
                        useGpuForPme ? &deviceStreamManager->context() : nullptr;
                const DeviceStream* pmeStream =
                        useGpuForPme ? &deviceStreamManager->stream(DeviceStreamType::Pme) : nullptr;

                pmedata = gmx_pme_init(cr, getNumPmeDomains(cr->dd), inputrec, nChargePerturbed != 0,
                                       nTypePerturbed != 0, mdrunOptions.reproducible, ewaldcoeff_q,
                                       ewaldcoeff_lj, gmx_omp_nthreads_get(emntPME), pmeRunMode,
                                       nullptr, deviceContext, pmeStream, pmeGpuProgram.get(), mdlog);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
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

    pull_t* pull_work = nullptr;
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Assumes uniform use of the number of OpenMP threads */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntDefault));

        if (inputrec->bPull)
        {
            /* Initialize pull code */
            pull_work = init_pull(fplog, inputrec->pull, inputrec, &mtop, cr, &atomSets,
                                  inputrec->fepvals->init_lambda);
            if (inputrec->pull->bXOutAverage || inputrec->pull->bFOutAverage)
            {
                initPullHistory(pull_work, &observablesHistory);
            }
            if (EI_DYNAMICS(inputrec->eI) && MASTER(cr))
            {
                init_pull_output_files(pull_work, filenames.size(), filenames.data(), oenv, startingBehavior);
            }
        }

        std::unique_ptr<EnforcedRotation> enforcedRotation;
        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            enforcedRotation =
                    init_rot(fplog, inputrec, filenames.size(), filenames.data(), cr, &atomSets,
                             globalState.get(), &mtop, oenv, mdrunOptions, startingBehavior);
        }

        t_swap* swap = nullptr;
        if (inputrec->eSwapCoords != eswapNO)
        {
            /* Initialize ion swapping code */
            swap = init_swapcoords(fplog, inputrec,
                                   opt2fn_master("-swap", filenames.size(), filenames.data(), cr),
                                   &mtop, globalState.get(), &observablesHistory, cr, &atomSets,
                                   oenv, mdrunOptions, startingBehavior);
        }

        /* Let makeConstraints know whether we have essential dynamics constraints. */
        auto constr = makeConstraints(mtop, *inputrec, pull_work, doEssentialDynamics, fplog,
                                      *mdAtoms->mdatoms(), cr, ms, &nrnb, wcycle, fr->bMolPBC);

        /* Energy terms and groups */
        gmx_enerdata_t enerd(mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                             inputrec->fepvals->n_lambda);

        /* Kinetic energy data */
        gmx_ekindata_t ekind;
        init_ekindata(fplog, &mtop, &(inputrec->opts), &ekind, inputrec->cos_accel);

        /* Set up interactive MD (IMD) */
        auto imdSession =
                makeImdSession(inputrec, cr, wcycle, &enerd, ms, &mtop, mdlog,
                               MASTER(cr) ? globalState->x.rvec_array() : nullptr, filenames.size(),
                               filenames.data(), oenv, mdrunOptions.imdOptions, startingBehavior);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because fr->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog, cr->dd, mtop, vsite.get(), inputrec,
                            domdecOptions.checkBondedInteractions, fr->cginfo_mb);
        }

        // TODO This is not the right place to manage the lifetime of
        // this data structure, but currently it's the easiest way to
        // make it work.
        MdrunScheduleWorkload runScheduleWork;
        // Also populates the simulation constant workload description.
        runScheduleWork.simulationWork =
                createSimulationWorkload(*inputrec, useGpuForNonbonded, pmeRunMode, useGpuForBonded,
                                         useGpuForUpdate, devFlags.enableGpuBufferOps,
                                         devFlags.enableGpuHaloExchange, devFlags.enableGpuPmePPComm);

        std::unique_ptr<gmx::StatePropagatorDataGpu> stateGpu;
        if (gpusWereDetected
            && ((useGpuForPme && thisRankHasDuty(cr, DUTY_PME))
                || runScheduleWork.simulationWork.useGpuBufferOps))
        {
            GpuApiCallBehavior transferKind = (inputrec->eI == eiMD && !doRerun && !useModularSimulator)
                                                      ? GpuApiCallBehavior::Async
                                                      : GpuApiCallBehavior::Sync;
            GMX_RELEASE_ASSERT(deviceStreamManager != nullptr,
                               "GPU device stream manager should be initialized to use GPU.");
            stateGpu = std::make_unique<gmx::StatePropagatorDataGpu>(
                    *deviceStreamManager, transferKind, pme_gpu_get_block_size(fr->pmedata), wcycle);
            fr->stateGpu = stateGpu.get();
        }

        GMX_ASSERT(stopHandlerBuilder_, "Runner must provide StopHandlerBuilder to simulator.");
        SimulatorBuilder simulatorBuilder;

        // build and run simulator object based on user-input
        auto simulator = simulatorBuilder.build(
                useModularSimulator, fplog, cr, ms, mdlog, static_cast<int>(filenames.size()),
                filenames.data(), oenv, mdrunOptions, startingBehavior, vsite.get(), constr.get(),
                enforcedRotation ? enforcedRotation->getLegacyEnfrot() : nullptr, deform.get(),
                mdModules_->outputProvider(), mdModules_->notifier(), inputrec, imdSession.get(),
                pull_work, swap, &mtop, fcd, globalState.get(), &observablesHistory, mdAtoms.get(),
                &nrnb, wcycle, fr, &enerd, &ekind, &runScheduleWork, replExParams, membed,
                walltime_accounting, std::move(stopHandlerBuilder_), doRerun);
        simulator->run();

        if (fr->pmePpCommGpu)
        {
            // destroy object since it is no longer required. (This needs to be done while the GPU context still exists.)
            fr->pmePpCommGpu.reset();
        }

        if (inputrec->bPull)
        {
            finish_pull(pull_work);
        }
        finish_swapcoords(swap);
    }
    else
    {
        GMX_RELEASE_ASSERT(pmedata, "pmedata was NULL while cr->duty was not DUTY_PP");
        /* do PME only */
        walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(emntPME));
        gmx_pmeonly(pmedata, cr, &nrnb, wcycle, walltime_accounting, inputrec, pmeRunMode,
                    deviceStreamManager.get());
    }

    wallcycle_stop(wcycle, ewcRUN);

    /* Finish up, write some stuff
     * if rerunMD, don't write last frame again
     */
    finish_run(fplog, mdlog, cr, inputrec, &nrnb, wcycle, walltime_accounting,
               fr ? fr->nbv.get() : nullptr, pmedata, EI_DYNAMICS(inputrec->eI) && !isMultiSim(ms));

    // clean up cycle counter
    wallcycle_destroy(wcycle);

    deviceStreamManager.reset(nullptr);
    // Free PME data
    if (pmedata)
    {
        gmx_pme_destroy(pmedata);
        pmedata = nullptr;
    }

    // FIXME: this is only here to manually unpin mdAtoms->chargeA_ and state->x,
    // before we destroy the GPU context(s) in free_gpu().
    // Pinned buffers are associated with contexts in CUDA.
    // As soon as we destroy GPU contexts after mdrunner() exits, these lines should go.
    mdAtoms.reset(nullptr);
    globalState.reset(nullptr);
    mdModules_.reset(nullptr); // destruct force providers here as they might also use the GPU
    gpuBonded.reset(nullptr);
    /* Free pinned buffers in *fr */
    delete fr;
    fr = nullptr;

    if (hwinfo->gpu_info.n_dev > 0)
    {
        /* stop the GPU profiler (only CUDA) */
        stopGpuProfiler();
    }

    /* With tMPI we need to wait for all ranks to finish deallocation before
     * destroying the CUDA context in free_gpu() as some tMPI ranks may be sharing
     * GPU and context.
     *
     * This is not a concern in OpenCL where we use one context per rank which
     * is freed in nbnxn_gpu_free().
     *
     * Note: it is safe to not call the barrier on the ranks which do not use GPU,
     * but it is easier and more futureproof to call it on the whole node.
     *
     * Note that this function needs to be called even if GPUs are not used
     * in this run because the PME ranks have no knowledge of whether GPUs
     * are used or not, but all ranks need to enter the barrier below.
     * \todo Remove this physical node barrier after making sure
     * that it's not needed anymore (with a shared GPU run).
     */
    if (GMX_THREAD_MPI)
    {
        physicalNodeComm.barrier();
    }

    free_gpu(deviceInfo);
    sfree(fcd);

    if (doMembed)
    {
        free_membed(membed);
    }

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltime_accounting);

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

    auto rc = static_cast<int>(gmx_get_stop_condition());

#if GMX_THREAD_MPI
    /* we need to join all threads. The sub-threads join when they
       exit this function, but the master thread needs to be told to
       wait for that. */
    if (PAR(cr) && MASTER(cr))
    {
        tMPI_Finalize();
    }
#endif
    return rc;
} // namespace gmx

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
                   "restraints added during runner life time should be cleared at runner "
                   "destruction.");
    }
};

void Mdrunner::addPotential(std::shared_ptr<gmx::IRestraintPotential> puller, const std::string& name)
{
    GMX_ASSERT(restraintManager_, "Mdrunner must have a restraint manager.");
    // Not sure if this should be logged through the md logger or something else,
    // but it is helpful to have some sort of INFO level message sent somewhere.
    //    std::cout << "Registering restraint named " << name << std::endl;

    // When multiple restraints are used, it may be wasteful to register them separately.
    // Maybe instead register an entire Restraint Manager as a force provider.
    restraintManager_->addToSpec(std::move(puller), name);
}

Mdrunner::Mdrunner(std::unique_ptr<MDModules> mdModules) : mdModules_(std::move(mdModules)) {}

Mdrunner::Mdrunner(Mdrunner&&) noexcept = default;

//NOLINTNEXTLINE(performance-noexcept-move-constructor) working around GCC bug 58265
Mdrunner& Mdrunner::operator=(Mdrunner&& /*handle*/) noexcept(BUGFREE_NOEXCEPT_STRING) = default;

class Mdrunner::BuilderImplementation
{
public:
    BuilderImplementation() = delete;
    BuilderImplementation(std::unique_ptr<MDModules> mdModules, compat::not_null<SimulationContext*> context);
    ~BuilderImplementation();

    BuilderImplementation& setExtraMdrunOptions(const MdrunOptions& options,
                                                real                forceWarningThreshold,
                                                StartingBehavior    startingBehavior);

    void addDomdec(const DomdecOptions& options);

    void addVerletList(int nstlist);

    void addReplicaExchange(const ReplicaExchangeParameters& params);

    void addNonBonded(const char* nbpu_opt);

    void addPME(const char* pme_opt_, const char* pme_fft_opt_);

    void addBondedTaskAssignment(const char* bonded_opt);

    void addUpdateTaskAssignment(const char* update_opt);

    void addHardwareOptions(const gmx_hw_opt_t& hardwareOptions);

    void addFilenames(ArrayRef<const t_filenm> filenames);

    void addOutputEnvironment(gmx_output_env_t* outputEnvironment);

    void addLogFile(t_fileio* logFileHandle);

    void addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder);

    Mdrunner build();

private:
    // Default parameters copied from runner.h
    // \todo Clarify source(s) of default parameters.

    const char* nbpu_opt_    = nullptr;
    const char* pme_opt_     = nullptr;
    const char* pme_fft_opt_ = nullptr;
    const char* bonded_opt_  = nullptr;
    const char* update_opt_  = nullptr;

    MdrunOptions mdrunOptions_;

    DomdecOptions domdecOptions_;

    ReplicaExchangeParameters replicaExchangeParameters_;

    //! Command-line override for the duration of a neighbor list with the Verlet scheme.
    int nstlist_ = 0;

    //! Multisim communicator handle.
    gmx_multisim_t* multiSimulation_;

    //! mdrun communicator
    MPI_Comm communicator_ = MPI_COMM_NULL;

    //! Print a warning if any force is larger than this (in kJ/mol nm).
    real forceWarningThreshold_ = -1;

    //! Whether the simulation will start afresh, or restart with/without appending.
    StartingBehavior startingBehavior_ = StartingBehavior::NewSimulation;

    //! The modules that comprise the functionality of mdrun.
    std::unique_ptr<MDModules> mdModules_;

    //! \brief Parallelism information.
    gmx_hw_opt_t hardwareOptions_;

    //! filename options for simulation.
    ArrayRef<const t_filenm> filenames_;

    /*! \brief Handle to output environment.
     *
     * \todo gmx_output_env_t needs lifetime management.
     */
    gmx_output_env_t* outputEnvironment_ = nullptr;

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

Mdrunner::BuilderImplementation::BuilderImplementation(std::unique_ptr<MDModules> mdModules,
                                                       compat::not_null<SimulationContext*> context) :
    mdModules_(std::move(mdModules))
{
    communicator_    = context->communicator_;
    multiSimulation_ = context->multiSimulation_.get();
}

Mdrunner::BuilderImplementation::~BuilderImplementation() = default;

Mdrunner::BuilderImplementation&
Mdrunner::BuilderImplementation::setExtraMdrunOptions(const MdrunOptions&    options,
                                                      const real             forceWarningThreshold,
                                                      const StartingBehavior startingBehavior)
{
    mdrunOptions_          = options;
    forceWarningThreshold_ = forceWarningThreshold;
    startingBehavior_      = startingBehavior;
    return *this;
}

void Mdrunner::BuilderImplementation::addDomdec(const DomdecOptions& options)
{
    domdecOptions_ = options;
}

void Mdrunner::BuilderImplementation::addVerletList(int nstlist)
{
    nstlist_ = nstlist;
}

void Mdrunner::BuilderImplementation::addReplicaExchange(const ReplicaExchangeParameters& params)
{
    replicaExchangeParameters_ = params;
}

Mdrunner Mdrunner::BuilderImplementation::build()
{
    auto newRunner = Mdrunner(std::move(mdModules_));

    newRunner.mdrunOptions     = mdrunOptions_;
    newRunner.pforce           = forceWarningThreshold_;
    newRunner.startingBehavior = startingBehavior_;
    newRunner.domdecOptions    = domdecOptions_;

    // \todo determine an invariant to check or confirm that all gmx_hw_opt_t objects are valid
    newRunner.hw_opt = hardwareOptions_;

    // No invariant to check. This parameter exists to optionally override other behavior.
    newRunner.nstlist_cmdline = nstlist_;

    newRunner.replExParams = replicaExchangeParameters_;

    newRunner.filenames = filenames_;

    newRunner.communicator = communicator_;

    // nullptr is a valid value for the multisim handle
    newRunner.ms = multiSimulation_;

    // \todo Clarify ownership and lifetime management for gmx_output_env_t
    // \todo Update sanity checking when output environment has clearly specified invariants.
    // Initialization and default values for oenv are not well specified in the current version.
    if (outputEnvironment_)
    {
        newRunner.oenv = outputEnvironment_;
    }
    else
    {
        GMX_THROW(gmx::APIError(
                "MdrunnerBuilder::addOutputEnvironment() is required before build()"));
    }

    newRunner.logFileHandle = logFileHandle_;

    if (nbpu_opt_)
    {
        newRunner.nbpu_opt = nbpu_opt_;
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
        GMX_THROW(gmx::APIError(
                "MdrunnerBuilder::addBondedTaskAssignment() is required before build()"));
    }

    if (update_opt_)
    {
        newRunner.update_opt = update_opt_;
    }
    else
    {
        GMX_THROW(gmx::APIError(
                "MdrunnerBuilder::addUpdateTaskAssignment() is required before build()  "));
    }


    newRunner.restraintManager_ = std::make_unique<gmx::RestraintManager>();

    if (stopHandlerBuilder_)
    {
        newRunner.stopHandlerBuilder_ = std::move(stopHandlerBuilder_);
    }
    else
    {
        newRunner.stopHandlerBuilder_ = std::make_unique<StopHandlerBuilder>();
    }

    return newRunner;
}

void Mdrunner::BuilderImplementation::addNonBonded(const char* nbpu_opt)
{
    nbpu_opt_ = nbpu_opt;
}

void Mdrunner::BuilderImplementation::addPME(const char* pme_opt, const char* pme_fft_opt)
{
    pme_opt_     = pme_opt;
    pme_fft_opt_ = pme_fft_opt;
}

void Mdrunner::BuilderImplementation::addBondedTaskAssignment(const char* bonded_opt)
{
    bonded_opt_ = bonded_opt;
}

void Mdrunner::BuilderImplementation::addUpdateTaskAssignment(const char* update_opt)
{
    update_opt_ = update_opt;
}

void Mdrunner::BuilderImplementation::addHardwareOptions(const gmx_hw_opt_t& hardwareOptions)
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

void Mdrunner::BuilderImplementation::addLogFile(t_fileio* logFileHandle)
{
    logFileHandle_ = logFileHandle;
}

void Mdrunner::BuilderImplementation::addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder)
{
    stopHandlerBuilder_ = std::move(builder);
}

MdrunnerBuilder::MdrunnerBuilder(std::unique_ptr<MDModules>           mdModules,
                                 compat::not_null<SimulationContext*> context) :
    impl_{ std::make_unique<Mdrunner::BuilderImplementation>(std::move(mdModules), context) }
{
}

MdrunnerBuilder::~MdrunnerBuilder() = default;

MdrunnerBuilder& MdrunnerBuilder::addSimulationMethod(const MdrunOptions&    options,
                                                      real                   forceWarningThreshold,
                                                      const StartingBehavior startingBehavior)
{
    impl_->setExtraMdrunOptions(options, forceWarningThreshold, startingBehavior);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addDomainDecomposition(const DomdecOptions& options)
{
    impl_->addDomdec(options);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addNeighborList(int nstlist)
{
    impl_->addVerletList(nstlist);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addReplicaExchange(const ReplicaExchangeParameters& params)
{
    impl_->addReplicaExchange(params);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addNonBonded(const char* nbpu_opt)
{
    impl_->addNonBonded(nbpu_opt);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addElectrostatics(const char* pme_opt, const char* pme_fft_opt)
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
        GMX_THROW(
                gmx::InvalidInputError("addElectrostatics() arguments must be non-null pointers."));
    }
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addBondedTaskAssignment(const char* bonded_opt)
{
    impl_->addBondedTaskAssignment(bonded_opt);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addUpdateTaskAssignment(const char* update_opt)
{
    impl_->addUpdateTaskAssignment(update_opt);
    return *this;
}

Mdrunner MdrunnerBuilder::build()
{
    return impl_->build();
}

MdrunnerBuilder& MdrunnerBuilder::addHardwareOptions(const gmx_hw_opt_t& hardwareOptions)
{
    impl_->addHardwareOptions(hardwareOptions);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addFilenames(ArrayRef<const t_filenm> filenames)
{
    impl_->addFilenames(filenames);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addOutputEnvironment(gmx_output_env_t* outputEnvironment)
{
    impl_->addOutputEnvironment(outputEnvironment);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addLogFile(t_fileio* logFileHandle)
{
    impl_->addLogFile(logFileHandle);
    return *this;
}

MdrunnerBuilder& MdrunnerBuilder::addStopHandlerBuilder(std::unique_ptr<StopHandlerBuilder> builder)
{
    impl_->addStopHandlerBuilder(std::move(builder));
    return *this;
}

MdrunnerBuilder::MdrunnerBuilder(MdrunnerBuilder&&) noexcept = default;

MdrunnerBuilder& MdrunnerBuilder::operator=(MdrunnerBuilder&&) noexcept = default;

} // namespace gmx
