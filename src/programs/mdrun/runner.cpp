/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "runner.h"

#include "config.h"

#include <cassert>
#include <csignal>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/ewald-utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/printhardware.h"
#include "gromacs/listed-forces/disre.h"
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/integrator.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/minimize.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_tuning.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/tpi.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_rotation.h"
#include "gromacs/taskassignment/hardwareassign.h"
#include "gromacs/taskassignment/resourcedivision.h"
#include "gromacs/taskassignment/usergpuids.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/filestream.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/loggerbuilder.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

#include "deform.h"
#include "md.h"
#include "membed.h"
#include "repl_ex.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

//! First step used in pressure scaling
gmx_int64_t         deform_init_init_step_tpx;
//! Initial box for pressure scaling
matrix              deform_init_box_tpx;
//! MPI variable for use in pressure scaling
tMPI_Thread_mutex_t deform_init_box_mutex = TMPI_THREAD_MUTEX_INITIALIZER;

#if GMX_THREAD_MPI
/* The minimum number of atoms per tMPI thread. With fewer atoms than this,
 * the number of threads will get lowered.
 */
#define MIN_ATOMS_PER_MPI_THREAD    90
#define MIN_ATOMS_PER_GPU           900

namespace gmx
{

void Mdrunner::reinitializeOnSpawnedThread()
{
    // TODO This duplication is formally necessary if any thread might
    // modify any memory in fnm or the pointers it contains. If the
    // contents are ever provably const, then we can remove this
    // allocation (and memory leak).
    // TODO This should probably become part of a copy constructor for
    // Mdrunner.
    fnm = dup_tfn(nfile, fnm);

    cr  = reinitialize_commrec_for_this_thread(cr);

    if (!MASTER(cr))
    {
        // Only the master rank writes to the log files
        fplog = nullptr;
    }
}

/*! \brief The callback used for running on spawned threads.
 *
 * Obtains the pointer to the master mdrunner object from the one
 * argument permitted to the thread-launch API call, copies it to make
 * a new runner for this thread, reinitializes necessary data, and
 * proceeds to the simulation. */
static void mdrunner_start_fn(void *arg)
{
    try
    {
        auto masterMdrunner = reinterpret_cast<const gmx::Mdrunner *>(arg);
        /* copy the arg list to make sure that it's thread-local. This
           doesn't copy pointed-to items, of course, but those are all
           const. */
        gmx::Mdrunner mdrunner = *masterMdrunner;
        mdrunner.reinitializeOnSpawnedThread();
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
t_commrec *Mdrunner::spawnThreads(int numThreadsToLaunch)
{

    /* first check whether we even need to start tMPI */
    if (numThreadsToLaunch < 2)
    {
        return cr;
    }

    gmx::Mdrunner spawnedMdrunner = *this;
    // TODO This duplication is formally necessary if any thread might
    // modify any memory in fnm or the pointers it contains. If the
    // contents are ever provably const, then we can remove this
    // allocation (and memory leak).
    // TODO This should probably become part of a copy constructor for
    // Mdrunner.
    spawnedMdrunner.fnm = dup_tfn(this->nfile, fnm);

    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns, we set thread affinity later */
    if (tMPI_Init_fn(TRUE, numThreadsToLaunch, TMPI_AFFINITY_NONE,
                     mdrunner_start_fn, static_cast<void*>(&spawnedMdrunner)) != TMPI_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to spawn thread-MPI threads"));
    }

    return reinitialize_commrec_for_this_thread(cr);
}

}      // namespace

#endif /* GMX_THREAD_MPI */

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
        verletbuf_list_setup_t ls;
        real                   rlist_new;

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        verletbuf_get_list_setup(true, makeGpuPairList, &ls);

        calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1, -1, &ls, nullptr, &rlist_new);

        if (rlist_new != ir->rlist)
        {
            if (fplog != nullptr)
            {
                fprintf(fplog, "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist, rlist_new,
                        ls.cluster_size_i, ls.cluster_size_j);
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
                                    gmx_int64_t          nsteps_cmdline,
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
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %d",
                  nsteps_cmdline);
    }
    /* Do nothing if nsteps_cmdline == -2 */
}

namespace gmx
{

//! Halt the run if there are inconsistences between user choices to run with GPUs and/or hardware detection.
static void exitIfCannotForceGpuRun(bool                requirePhysicalGpu,
                                    EmulateGpuNonbonded emulateGpuNonbonded,
                                    bool                useVerletScheme,
                                    bool                compatibleGpusFound)
{
    /* Was GPU acceleration either explicitly (-nb gpu) or implicitly
     * (gpu ID passed) requested? */
    if (!requirePhysicalGpu)
    {
        return;
    }

    if (GMX_GPU == GMX_GPU_NONE)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but %s was compiled without GPU support!",
                  gmx::getProgramContext().displayName());
    }

    if (emulateGpuNonbonded == EmulateGpuNonbonded::Yes)
    {
        gmx_fatal(FARGS, "GPU emulation cannot be requested together with GPU acceleration!");
    }

    if (!useVerletScheme)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but can't be used without cutoff-scheme=Verlet");
    }

    if (!compatibleGpusFound)
    {
        gmx_fatal(FARGS, "GPU acceleration requested, but no compatible GPUs were detected.");
    }
}

/*! \brief Return whether GPU acceleration is useful with the given settings.
 *
 * If not, logs a message about falling back to CPU code. */
static bool gpuAccelerationIsUseful(const MDLogger   &mdlog,
                                    const t_inputrec *ir,
                                    bool              doRerun)
{
    if (doRerun && ir->opts.ngener > 1)
    {
        /* Rerun execution time is dominated by I/O and pair search,
         * so GPUs are not very useful, plus they do not support more
         * than one energy group. If the user requested GPUs
         * explicitly, a fatal error is given later.  With non-reruns,
         * we fall back to a single whole-of system energy group
         * (which runs much faster than a multiple-energy-groups
         * implementation would), and issue a note in the .log
         * file. Users can re-run if they want the information. */
        GMX_LOG(mdlog.warning).asParagraph().appendText("Multiple energy groups is not implemented for GPUs, so is not useful for this rerun, so falling back to the CPU");
        return false;
    }

    return true;
}

//! \brief Return the correct integrator function.
static integrator_t *my_integrator(unsigned int ei)
{
    switch (ei)
    {
        case eiMD:
        case eiBD:
        case eiSD1:
        case eiVV:
        case eiVVAK:
            if (!EI_DYNAMICS(ei))
            {
                GMX_THROW(APIError("do_md integrator would be called for a non-dynamical integrator"));
            }
            return do_md;
        case eiSteep:
            return do_steep;
        case eiCG:
            return do_cg;
        case eiNM:
            return do_nm;
        case eiLBFGS:
            return do_lbfgs;
        case eiTPI:
        case eiTPIC:
            if (!EI_TPI(ei))
            {
                GMX_THROW(APIError("do_tpi integrator would be called for a non-TPI integrator"));
            }
            return do_tpi;
        case eiSD2_REMOVED:
            GMX_THROW(NotImplementedError("SD2 integrator has been removed"));
        default:
            GMX_THROW(APIError("Non existing integrator selected"));
    }
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

int Mdrunner::mdrunner()
{
    matrix                    box;
    gmx_ddbox_t               ddbox = {0};
    int                       npme_major, npme_minor;
    t_nrnb                   *nrnb;
    gmx_mtop_t               *mtop          = nullptr;
    t_forcerec               *fr            = nullptr;
    t_fcdata                 *fcd           = nullptr;
    real                      ewaldcoeff_q  = 0;
    real                      ewaldcoeff_lj = 0;
    gmx_vsite_t              *vsite         = nullptr;
    gmx_constr_t              constr;
    int                       nChargePerturbed = -1, nTypePerturbed = 0;
    gmx_wallcycle_t           wcycle;
    gmx_walltime_accounting_t walltime_accounting = nullptr;
    int                       rc;
    gmx_int64_t               reset_counters;
    int                       nthreads_pme = 1;
    gmx_membed_t *            membed       = nullptr;
    gmx_hw_info_t            *hwinfo       = nullptr;

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    gmx::MDModules mdModules;
    t_inputrec     inputrecInstance;
    t_inputrec    *inputrec = &inputrecInstance;
    snew(mtop, 1);

    if (mdrunOptions.continuationOptions.appendFiles)
    {
        fplog = nullptr;
    }

    bool doMembed = opt2bSet("-membed", nfile, fnm);
    bool doRerun  = mdrunOptions.rerun;

    /* Handle GPU-related user options. Later, we check consistency
     * with things like whether support is compiled, or tMPI thread
     * count. */
    EmulateGpuNonbonded emulateGpuNonbonded = (getenv("GMX_EMULATE_GPU") != nullptr ?
                                               EmulateGpuNonbonded::Yes : EmulateGpuNonbonded::No);
    std::vector<int>    userGpuIds;
    try
    {
        userGpuIds = parseUserGpuIds(hw_opt.gpuIdTaskAssignment);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    bool                forceUseCpu           = (strncmp(nbpu_opt, "cpu", 3) == 0);
    if (!userGpuIds.empty() && forceUseCpu)
    {
        gmx_fatal(FARGS, "GPU IDs were specified, and short-ranged interactions were assigned to the CPU. Make no more than one of these choices.");
    }
    bool forceUsePhysicalGpu = (strncmp(nbpu_opt, "gpu", 3) == 0) || !userGpuIds.empty();
    bool tryUsePhysicalGpu   = (strncmp(nbpu_opt, "auto", 4) == 0) && userGpuIds.empty() && (emulateGpuNonbonded == EmulateGpuNonbonded::No);
    GMX_RELEASE_ASSERT(!(forceUsePhysicalGpu && tryUsePhysicalGpu), "Must either force use of "
                       "GPUs for short-ranged interactions, or try to use them, not both.");
    const PmeRunMode pmeRunMode = PmeRunMode::CPU;
    //TODO this is a placeholder as PME on GPU is not permitted yet
    //TODO should there exist a PmeRunMode::None value for consistency?

    // Here we assume that SIMMASTER(cr) does not change even after the
    // threads are started.
    gmx::LoggerOwner logOwner(buildLogger(fplog, cr));
    gmx::MDLogger    mdlog(logOwner.logger());

    hwinfo = gmx_detect_hardware(mdlog, cr);

    gmx_print_detected_hardware(fplog, cr, mdlog, hwinfo);

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
    }

    std::unique_ptr<t_state> globalState;

    if (SIMMASTER(cr))
    {
        /* Only the master rank has the global state */
        globalState = std::unique_ptr<t_state>(new t_state);

        /* Read (nearly) all data required for the simulation */
        read_tpx_state(ftp2fn(efTPR, nfile, fnm), inputrec, globalState.get(), mtop);

        exitIfCannotForceGpuRun(forceUsePhysicalGpu,
                                emulateGpuNonbonded,
                                inputrec->cutoff_scheme == ecutsVERLET,
                                compatibleGpusFound(hwinfo->gpu_info));

        if (inputrec->cutoff_scheme == ecutsVERLET)
        {
            /* TODO This logic could run later, e.g. before -npme -1
               is handled. If inputrec has already been communicated,
               then the resulting tryUsePhysicalGpu does not need to
               be communicated. */
            if ((tryUsePhysicalGpu || forceUsePhysicalGpu) &&
                !gpuAccelerationIsUseful(mdlog, inputrec, doRerun))
            {
                /* Fallback message printed by nbnxn_acceleration_supported */
                if (forceUsePhysicalGpu)
                {
                    gmx_fatal(FARGS, "GPU acceleration requested, but not supported with the given input settings");
                }
                tryUsePhysicalGpu = false;
            }
        }
        else
        {
            if (nstlist_cmdline > 0)
            {
                gmx_fatal(FARGS, "Can not set nstlist with the group cut-off scheme");
            }

            if (compatibleGpusFound(hwinfo->gpu_info))
            {
                GMX_LOG(mdlog.warning).asParagraph().appendText(
                        "NOTE: GPU(s) found, but the current simulation can not use GPUs\n"
                        "      To use a GPU, set the mdp option: cutoff-scheme = Verlet");
            }
            tryUsePhysicalGpu = false;
        }
    }
    bool nonbondedOnGpu = (tryUsePhysicalGpu || forceUsePhysicalGpu) && compatibleGpusFound(hwinfo->gpu_info);

    /* Check and update the hardware options for internal consistency */
    check_and_update_hw_opt_1(&hw_opt, cr, domdecOptions.numPmeRanks);

    /* Early check for externally set process affinity. */
    gmx_check_thread_affinity_set(mdlog, cr,
                                  &hw_opt, hwinfo->nthreads_hw_avail, FALSE);

#if GMX_THREAD_MPI
    if (SIMMASTER(cr))
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

        /* Determine how many thread-MPI ranks to start.
         *
         * TODO Over-writing the user-supplied value here does
         * prevent any possible subsequent checks from working
         * correctly. */
        hw_opt.nthreads_tmpi = get_nthreads_mpi(hwinfo,
                                                &hw_opt,
                                                userGpuIds,
                                                domdecOptions.numPmeRanks,
                                                nonbondedOnGpu,
                                                inputrec, mtop,
                                                mdlog,
                                                doMembed);

        // Now start the threads for thread MPI.
        cr = spawnThreads(hw_opt.nthreads_tmpi);
        /* The main thread continues here with a new cr. We don't deallocate
           the old cr because other threads may still be reading it. */
        // TODO Both master and spawned threads call dup_tfn and
        // reinitialize_commrec_for_this_thread. Find a way to express
        // this better.
    }
#endif
    /* END OF CAUTION: cr is now reliable */

    if (PAR(cr))
    {
        /* now broadcast everything to the non-master nodes/threads: */
        init_parallel(cr, inputrec, mtop);

        gmx_bcast_sim(sizeof(nonbondedOnGpu), &nonbondedOnGpu, cr);
    }
    // TODO: Error handling
    mdModules.assignOptionsToModules(*inputrec->params, nullptr);

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
            globalState = std::unique_ptr<t_state>(new t_state);
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
                  "but %s was compiled without threads or MPI enabled"
#else
#if GMX_THREAD_MPI
                  "but the number of MPI-threads (option -ntmpi) is not set or is 1"
#else
                  "but %s was not started through mpirun/mpiexec or only one rank was requested through mpirun/mpiexec"
#endif
#endif
                  , output_env_get_program_display_name(oenv)
                  );
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

    if (nonbondedOnGpu && domdecOptions.numPmeRanks < 0)
    {
        /* With GPUs we don't automatically use PME-only ranks. PME ranks can
         * improve performance with many threads per GPU, since our OpenMP
         * scaling is bad, but it's difficult to automate the setup.
         */
        domdecOptions.numPmeRanks = 0;
    }

#ifdef GMX_FAHCORE
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
    init_disres(fplog, mtop, inputrec, cr, fcd, globalState.get(), replExParams.exchangeInterval > 0);

    init_orires(fplog, mtop, inputrec, cr, globalState.get(), &(fcd->orires));

    if (inputrecDeform(inputrec))
    {
        /* Store the deform reference box before reading the checkpoint */
        if (SIMMASTER(cr))
        {
            copy_mat(globalState->box, box);
        }
        if (PAR(cr))
        {
            gmx_bcast(sizeof(box), box, cr);
        }
        /* Because we do not have the update struct available yet
         * in which the reference values should be stored,
         * we store them temporarily in static variables.
         * This should be thread safe, since they are only written once
         * and with identical values.
         */
        tMPI_Thread_mutex_lock(&deform_init_box_mutex);
        deform_init_init_step_tpx = inputrec->init_step;
        copy_mat(box, deform_init_box_tpx);
        tMPI_Thread_mutex_unlock(&deform_init_box_mutex);
    }

    ObservablesHistory   observablesHistory = {};

    ContinuationOptions &continuationOptions = mdrunOptions.continuationOptions;

    if (continuationOptions.startedFromCheckpoint)
    {
        /* Check if checkpoint file exists before doing continuation.
         * This way we can use identical input options for the first and subsequent runs...
         */
        gmx_bool bReadEkin;

        load_checkpoint(opt2fn_master("-cpi", nfile, fnm, cr), &fplog,
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
    }

    if (SIMMASTER(cr) && continuationOptions.appendFiles)
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr,
                     continuationOptions.appendFiles, &fplog);
        logOwner = buildLogger(fplog, nullptr);
        mdlog    = logOwner.logger();
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
        prepare_verlet_scheme(fplog, cr, inputrec, nstlist_cmdline, mtop, box,
                              nonbondedOnGpu || (emulateGpuNonbonded == EmulateGpuNonbonded::Yes), *hwinfo->cpuInfo);
    }

    if (PAR(cr) && !(EI_TPI(inputrec->eI) ||
                     inputrec->eI == eiNM))
    {
        const rvec *xOnMaster = (SIMMASTER(cr) ? as_rvec_array(globalState->x.data()) : nullptr);

        cr->dd = init_domain_decomposition(fplog, cr, domdecOptions, mdrunOptions,
                                           mtop, inputrec,
                                           box, xOnMaster,
                                           &ddbox, &npme_major, &npme_minor);
        // Note that local state still does not exist yet.
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->npmenodes = 0;
        cr->duty      = (DUTY_PP | DUTY_PME);
        npme_major    = 1;
        npme_minor    = 1;

        if (inputrec->ePBC == epbcSCREW)
        {
            gmx_fatal(FARGS,
                      "pbc=%s is only implemented with domain decomposition",
                      epbc_names[inputrec->ePBC]);
        }
    }

    if (PAR(cr))
    {
        /* After possible communicator splitting in make_dd_communicators.
         * we can set up the intra/inter node communication.
         */
        gmx_setup_nodecomm(fplog, cr);
    }

    /* Initialize per-physical-node MPI process/thread ID and counters. */
    gmx_init_intranode_counters(cr);
#if GMX_MPI
    if (MULTISIM(cr))
    {
        GMX_LOG(mdlog.warning).asParagraph().appendTextFormatted(
                "This is simulation %d out of %d running as a composite GROMACS\n"
                "multi-simulation job. Setup for this simulation:\n",
                cr->ms->sim, cr->ms->nsim);
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

    /* Check and update hw_opt for the number of MPI ranks */
    check_and_update_hw_opt_3(&hw_opt);

    gmx_omp_nthreads_init(mdlog, cr,
                          hwinfo->nthreads_hw_avail,
                          hw_opt.nthreads_omp,
                          hw_opt.nthreads_omp_pme,
                          !thisRankHasDuty(cr, DUTY_PP),
                          inputrec->cutoff_scheme == ecutsVERLET);

#ifndef NDEBUG
    if (EI_TPI(inputrec->eI) &&
        inputrec->cutoff_scheme == ecutsVERLET)
    {
        gmx_feenableexcept();
    }
#endif

    // Contains the ID of the GPU used by each PP rank on this node,
    // indexed by that rank. Empty if no GPUs are selected for use on
    // this node.
    std::vector<int> gpuTaskAssignment;
    if (nonbondedOnGpu)
    {
        // Currently the DD code assigns duty to ranks that can
        // include PP work that currently can be executed on a single
        // GPU, if present and compatible.  This has to be coordinated
        // across PP ranks on a node, with possible multiple devices
        // or sharing devices on a node, either from the user
        // selection, or automatically.
        //
        // GPU ID assignment strings, if provided, cover all the ranks on
        // a node. If nodes or the process placement on them are
        // heterogeneous, then the GMX_GPU_ID environment variable must be
        // set by a user who also wishes to direct GPU ID assignment.
        // Thus the implementation of task assignment can assume it has a
        // GPU ID assignment appropriate for the node upon which its
        // process is running.
        //
        // Valid GPU ID assignments are an ordered set of digits that
        // identify GPU device IDs (e.g. as understood by the GPU runtime,
        // and subject to environment modification such as with
        // CUDA_VISIBLE_DEVICES) that will be used for the GPU-suitable
        // tasks on all of the ranks of that node.
        bool rankCanUseGpu = thisRankHasDuty(cr, DUTY_PP);
        gpuTaskAssignment = mapPpRanksToGpus(rankCanUseGpu, cr, hwinfo->gpu_info, hwinfo->compatibleGpus, userGpuIds);
    }

    reportGpuUsage(mdlog, hwinfo->gpu_info, !userGpuIds.empty(),
                   gpuTaskAssignment, cr->nrank_pp_intranode, cr->nnodes > 1);

    if (!gpuTaskAssignment.empty())
    {
        GMX_RELEASE_ASSERT(cr->nrank_pp_intranode == static_cast<int>(gpuTaskAssignment.size()),
                           "The number of PP ranks on each node must equal the number of GPU tasks used on each node");
    }

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
#endif

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo, hw_opt.nthreads_tot, !gpuTaskAssignment.empty(), mdrunOptions.ntompOptionIsSet,
                                       cr, mdlog);

    gmx_device_info_t *shortRangedDeviceInfo = nullptr;
    int                shortRangedDeviceId   = -1;
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        if (!gpuTaskAssignment.empty())
        {
            shortRangedDeviceId   = gpuTaskAssignment[cr->rank_pp_intranode];
            shortRangedDeviceInfo = getDeviceInfo(hwinfo->gpu_info, shortRangedDeviceId);
        }
    }

    if (DOMAINDECOMP(cr))
    {
        /* When we share GPUs over ranks, we need to know this for the DLB */
        dd_setup_dlb_resource_sharing(cr, shortRangedDeviceId);
    }

    /* getting number of PP/PME threads
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    nthreads_pme = gmx_omp_nthreads_get(emntPME);

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
        membed = init_membed(fplog, nfile, fnm, mtop, inputrec, globalState.get(), cr, &mdrunOptions.checkpointOptions.period);
    }

    std::unique_ptr<MDAtoms> mdAtoms;

    snew(nrnb, 1);
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Initiate forcerecord */
        fr                 = mk_forcerec();
        fr->forceProviders = mdModules.initForceProviders();
        init_forcerec(fplog, mdlog, fr, fcd,
                      inputrec, mtop, cr, box,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      getFilenm("-tableb", nfile, fnm),
                      shortRangedDeviceInfo,
                      FALSE,
                      pforce);

        /* Initialize QM-MM */
        if (fr->bQMMM)
        {
            init_QMMMrec(cr, mtop, inputrec, fr);
        }

        /* Initialize the mdAtoms structure.
         * mdAtoms is not filled with atom data,
         * as this can not be done now with domain decomposition.
         */
        mdAtoms = makeMDAtoms(fplog, *mtop, *inputrec);

        /* Initialize the virtual site communication */
        vsite = initVsite(*mtop, cr);

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
                rvec *xGlobal = as_rvec_array(globalState->x.data());
                do_pbc_first_mtop(fplog, inputrec->ePBC, box, mtop, xGlobal);
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                constructVsitesGlobal(*mtop, globalState->x);
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

    if (hw_opt.thread_affinity != threadaffOFF)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(mdlog, cr,
                                      &hw_opt, hwinfo->nthreads_hw_avail, TRUE);

        int nthread_local;
        /* threads on this MPI process or TMPI thread */
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            nthread_local = gmx_omp_nthreads_get(emntNonbonded);
        }
        else
        {
            nthread_local = gmx_omp_nthreads_get(emntPME);
        }

        /* Set the CPU affinity */
        gmx_set_thread_affinity(mdlog, cr, &hw_opt, *hwinfo->hardwareTopology,
                                nthread_local, nullptr);
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
                gmx_device_info_t *pmeGpuInfo = nullptr;
                pmedata = gmx_pme_init(cr, npme_major, npme_minor, inputrec,
                                       mtop ? mtop->natoms : 0, nChargePerturbed, nTypePerturbed,
                                       mdrunOptions.reproducible,
                                       ewaldcoeff_q, ewaldcoeff_lj,
                                       nthreads_pme,
                                       pmeRunMode, nullptr, pmeGpuInfo, mdlog);
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
                init_pull(fplog, inputrec->pull, inputrec, nfile, fnm,
                          mtop, cr, oenv, inputrec->fepvals->init_lambda,
                          EI_DYNAMICS(inputrec->eI) && MASTER(cr),
                          continuationOptions);
        }

        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            init_rot(fplog, inputrec, nfile, fnm, cr, globalState.get(), mtop, oenv, mdrunOptions);
        }

        /* Let init_constraints know whether we have essential dynamics constraints.
         * TODO: inputrec should tell us whether we use an algorithm, not a file option or the checkpoint
         */
        bool doEdsam = (opt2fn_null("-ei", nfile, fnm) != nullptr || observablesHistory.edsamHistory);

        constr = init_constraints(fplog, mtop, inputrec, doEdsam, cr);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because fr->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog, cr->dd, mtop, vsite, inputrec,
                            domdecOptions.checkBondedInteractions,
                            fr->cginfo_mb);
        }

        /* Now do whatever the user wants us to do (how flexible...) */
        my_integrator(inputrec->eI) (fplog, cr, mdlog, nfile, fnm,
                                     oenv,
                                     mdrunOptions,
                                     vsite, constr,
                                     mdModules.outputProvider(),
                                     inputrec, mtop,
                                     fcd,
                                     globalState.get(),
                                     &observablesHistory,
                                     mdAtoms.get(), nrnb, wcycle, fr,
                                     replExParams,
                                     membed,
                                     walltime_accounting);

        if (inputrec->bRot)
        {
            finish_rot(inputrec->rot);
        }

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
               EI_DYNAMICS(inputrec->eI) && !MULTISIM(cr));

    // Free PME data
    if (pmedata)
    {
        gmx_pme_destroy(pmedata);
        pmedata = nullptr;
    }

    /* Free GPU memory and context */
    free_gpu_resources(fr, cr, shortRangedDeviceInfo);

    if (doMembed)
    {
        free_membed(membed);
    }

    gmx_hardware_info_free();

    /* Does what it says */
    print_date_and_time(fplog, cr->nodeid, "Finished mdrun", gmx_gettime());
    walltime_accounting_destroy(walltime_accounting);

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && continuationOptions.appendFiles)
    {
        gmx_log_close(fplog);
    }

    rc = (int)gmx_get_stop_condition();

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
}

} // namespace gmx
