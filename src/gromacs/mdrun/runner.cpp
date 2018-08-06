/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include <csignal>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include "gromacs/restraint/manager.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/mdrunutility/handlerestart.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/ewald/ewald-utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme-gpu-program.h"
#include "gromacs/fileio/checkpoint.h"
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
#include "gromacs/listed-forces/orires.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/main.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdlib/nbnxn_tuning.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/repl_ex.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/tpxstate.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/restraint/restraintmdmodule.h"
#include "gromacs/pulling/pull_rotation.h"
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
#include "context.h"

#ifdef GMX_FAHCORE
#include "corewrap.h"
#endif

template<typename T>
bool stringIsEmpty(T string)
{ return string == nullptr; };

template<>
bool stringIsEmpty<std::string>(std::string string)
{ return string.empty(); };


/*! \brief Return whether the command-line parameter that
 *  will trigger a multi-simulation is set */
static bool is_multisim_option_set(int argc, const char *const argv[])
{
    for (int i = 0; i < argc; ++i)
    {
        if (strcmp(argv[i], "-multidir") == 0)
        {
            return true;
        }
    }
    return false;
}

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

std::unique_ptr<Mdrunner> Mdrunner::cloneOnSpawnedThread() const
{
    auto newRunner = gmx::compat::make_unique<Mdrunner>();

    // Todo: how to handle the restraint manager or parameters not in inputrec?

    newRunner->hw_opt = hw_opt;
    // TODO This duplication is formally necessary if any thread might
    // modify any memory in fnm or the pointers it contains. If the
    // contents are ever provably const, then we can remove this
    // allocation (and memory leak).
//    newRunner->fnm = dup_tfn(nfile, fnm);
    newRunner->oenv = oenv;
//    newRunner->bVerbose = bVerbose;
//    newRunner->nstglobalcomm = nstglobalcomm;
//    newRunner->ddxyz[0] = ddxyz[0];
//    newRunner->ddxyz[1] = ddxyz[1];
//    newRunner->ddxyz[2] = ddxyz[2];
//    newRunner->dd_rank_order = dd_rank_order;
//    newRunner->npme = npme;
//    newRunner->rdd = rdd;
//    newRunner->rconstr = rconstr;
//    newRunner->dddlb_opt = dddlb_opt;
//    newRunner->dlb_scale = dlb_scale;
//    newRunner->ddcsx = ddcsx;
//    newRunner->ddcsy = ddcsy;
//    newRunner->ddcsz = ddcsz;
    newRunner->nbpu_opt        = nbpu_opt;
    newRunner->nstlist_cmdline = nstlist_cmdline;
//    newRunner->nsteps_cmdline = nsteps_cmdline;
//    newRunner->nstepout = nstepout;
//    newRunner->resetstep = resetstep;
//    newRunner->nmultisim = nmultisim;
    newRunner->replExParams = replExParams;
    newRunner->pforce       = pforce;
//    newRunner->cpt_period = cpt_period;
//    newRunner->max_hours = max_hours;
//    newRunner->imdport = imdport;
//    newRunner->Flags = Flags;
    newRunner->cr  = reinitialize_commrec_for_this_thread(cr);
    // Don't copy fplog file pointer.

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
           doesn't copy pointed-to items, of course, but those are all
           const. */
        std::unique_ptr<gmx::Mdrunner> mdrunner = masterMdrunner->cloneOnSpawnedThread();
        mdrunner->mdrunner();
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
#endif  /* GMX_THREAD_MPI */

    return reinitialize_commrec_for_this_thread(cr);
}

}      // namespace gmx

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
        gmx_fatal(FARGS, "Invalid nsteps value passed on the command line: %ld",
                  nsteps_cmdline);
    }
    /* Do nothing if nsteps_cmdline == -2 */
}

void gmx::Mdrunner::initFromAPI(const std::vector<std::string> &args)
{
    if ((tpxState_ == nullptr) || (!tpxState_->isInitialized()))
    {
        gmx_fatal(FARGS, "Need initialized input record to initialize runner.");
    }
    // Until the options processing gets picked apart more (at least the fnm handling)
    // we're just spoofing argv and wrapping initFromCLI. Note that a non-const argv is deeply
    // embedded in GROMACS.
    constexpr int       offset = 3; // need placeholders for argv[0] and tpr file
    int                 argc   = offset + static_cast<int>(args.size());
    std::vector<char *> argv {
        static_cast<size_t>(argc)
    };

    argv[0] = new char[1]; // Start with an empty string (doesn't really matter)
    strcpy(argv[0], "");
    argv[1] = new char[3];
    strcpy(argv[1], "-s");
    argv[2] = new char[strlen(tpxState_->filename()) + 1];
    strcpy(argv[2], tpxState_->filename());
    argv[2][strlen(tpxState_->filename()) + 1] = 0;

    for (size_t idx_args = 0; idx_args < args.size(); ++idx_args)
    {
        auto idx_argv = idx_args + offset;
        argv[idx_argv] = new char[args[idx_args].length() + 1];
        strcpy(argv[idx_argv], args[idx_args].c_str());
    }

    initFromCLI(argc, argv.data());

    for (auto && string : argv)
    {
        delete[] string;
        string = nullptr;
    }
}

void gmx::Mdrunner::initFromCLI(int argc, char *argv[])
{
    const char       *desc[] = {""};
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
    int               rc;

    cr = init_commrec();

    unsigned long PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multidir, the working directory still needs to be
    // changed, so we can't check for the existence of files during
    // parsing.  It isn't useful to do any completion based on file
    // system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    if (!parse_common_args(&argc, argv, PCA_Flags, nfile, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(cr);
        return;
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
    gmx::ArrayRef<const std::string> multidir = opt2fnsIfOptionSet("-multidir", nfile, fnm);

    if (replExParams.exchangeInterval != 0 && multidir.size() < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (use option -multidir)");
    }

    if (replExParams.numExchanges < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    ms = init_multisystem(MPI_COMM_WORLD, multidir);

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

    if (!opt2bSet("-cpi", nfile, fnm))
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

    handleRestart(cr, ms, bTryToAppendFiles, nfile, fnm, &continuationOptions.appendFiles, &continuationOptions.startedFromCheckpoint);

    mdrunOptions.rerun            = opt2bSet("-rerun", nfile, fnm);
    mdrunOptions.ntompOptionIsSet = opt2parg_bSet("-ntomp", asize(pa), pa);

    /* We postpone opening the log file if we are appending, so we can
       first truncate the old log file and append to the correct position
       there instead.  */
    if (MASTER(cr) && !continuationOptions.appendFiles)
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr,
                     continuationOptions.appendFiles, &fplog);
    }
    else
    {
        fplog = nullptr;
    }

}

namespace gmx
{

int Mdrunner::mainFunction(int argc, char *argv[])
{
    const char   *desc[] = {
        "[THISMODULE] is the main computational chemistry engine",
        "within GROMACS. Obviously, it performs Molecular Dynamics simulations,",
        "but it can also perform Stochastic Dynamics, Energy Minimization,",
        "test particle insertion or (re)calculation of energies.",
        "Normal mode analysis is another option. In this case [TT]mdrun[tt]",
        "builds a Hessian matrix from single conformation.",
        "For usual Normal Modes-like calculations, make sure that",
        "the structure provided is properly energy-minimized.",
        "The generated matrix can be diagonalized by [gmx-nmeig].[PAR]",
        "The [TT]mdrun[tt] program reads the run input file ([TT]-s[tt])",
        "and distributes the topology over ranks if needed.",
        "[TT]mdrun[tt] produces at least four output files.",
        "A single log file ([TT]-g[tt]) is written.",
        "The trajectory file ([TT]-o[tt]), contains coordinates, velocities and",
        "optionally forces.",
        "The structure file ([TT]-c[tt]) contains the coordinates and",
        "velocities of the last step.",
        "The energy file ([TT]-e[tt]) contains energies, the temperature,",
        "pressure, etc, a lot of these things are also printed in the log file.",
        "Optionally coordinates can be written to a compressed trajectory file",
        "([TT]-x[tt]).[PAR]",
        "The option [TT]-dhdl[tt] is only used when free energy calculation is",
        "turned on.[PAR]",
        "Running mdrun efficiently in parallel is a complex topic topic,",
        "many aspects of which are covered in the online User Guide. You",
        "should look there for practical advice on using many of the options",
        "available in mdrun.[PAR]",
        "ED (essential dynamics) sampling and/or additional flooding potentials",
        "are switched on by using the [TT]-ei[tt] flag followed by an [REF].edi[ref]",
        "file. The [REF].edi[ref] file can be produced with the [TT]make_edi[tt] tool",
        "or by using options in the essdyn menu of the WHAT IF program.",
        "[TT]mdrun[tt] produces a [REF].xvg[ref] output file that",
        "contains projections of positions, velocities and forces onto selected",
        "eigenvectors.[PAR]",
        "When user-defined potential functions have been selected in the",
        "[REF].mdp[ref] file the [TT]-table[tt] option is used to pass [TT]mdrun[tt]",
        "a formatted table with potential functions. The file is read from",
        "either the current directory or from the [TT]GMXLIB[tt] directory.",
        "A number of pre-formatted tables are presented in the [TT]GMXLIB[tt] dir,",
        "for 6-8, 6-9, 6-10, 6-11, 6-12 Lennard-Jones potentials with",
        "normal Coulomb.",
        "When pair interactions are present, a separate table for pair interaction",
        "functions is read using the [TT]-tablep[tt] option.[PAR]",
        "When tabulated bonded functions are present in the topology,",
        "interaction functions are read using the [TT]-tableb[tt] option.",
        "For each different tabulated interaction type used, a table file name must",
        "be given. For the topology to work, a file name given here must match a",
        "character sequence before the file extension. That sequence is: an underscore,",
        "then a 'b' for bonds, an 'a' for angles or a 'd' for dihedrals,",
        "and finally the matching table number index used in the topology.[PAR]",
        "The options [TT]-px[tt] and [TT]-pf[tt] are used for writing pull COM",
        "coordinates and forces when pulling is selected",
        "in the [REF].mdp[ref] file.[PAR]",
        "Finally some experimental algorithms can be tested when the",
        "appropriate options have been given. Currently under",
        "investigation are: polarizability.",
        "[PAR]",
        "The option [TT]-membed[tt] does what used to be g_membed, i.e. embed",
        "a protein into a membrane. This module requires a number of settings",
        "that are provided in a data file that is the argument of this option.",
        "For more details in membrane embedding, see the documentation in the",
        "user guide. The options [TT]-mn[tt] and [TT]-mp[tt] are used to provide",
        "the index and topology files used for the embedding.",
        "[PAR]",
        "The option [TT]-pforce[tt] is useful when you suspect a simulation",
        "crashes due to too large forces. With this option coordinates and",
        "forces of atoms with a force larger than a certain value will",
        "be printed to stderr. It will also terminate the run when non-finite",
        "forces are present.",
        "[PAR]",
        "Checkpoints containing the complete state of the system are written",
        "at regular intervals (option [TT]-cpt[tt]) to the file [TT]-cpo[tt],",
        "unless option [TT]-cpt[tt] is set to -1.",
        "The previous checkpoint is backed up to [TT]state_prev.cpt[tt] to",
        "make sure that a recent state of the system is always available,",
        "even when the simulation is terminated while writing a checkpoint.",
        "With [TT]-cpnum[tt] all checkpoint files are kept and appended",
        "with the step number.",
        "A simulation can be continued by reading the full state from file",
        "with option [TT]-cpi[tt]. This option is intelligent in the way that",
        "if no checkpoint file is found, GROMACS just assumes a normal run and",
        "starts from the first step of the [REF].tpr[ref] file. By default the output",
        "will be appending to the existing output files. The checkpoint file",
        "contains checksums of all output files, such that you will never",
        "loose data when some output files are modified, corrupt or removed.",
        "There are three scenarios with [TT]-cpi[tt]:[PAR]",
        "[TT]*[tt] no files with matching names are present: new output files are written[PAR]",
        "[TT]*[tt] all files are present with names and checksums matching those stored",
        "in the checkpoint file: files are appended[PAR]",
        "[TT]*[tt] otherwise no files are modified and a fatal error is generated[PAR]",
        "With [TT]-noappend[tt] new output files are opened and the simulation",
        "part number is added to all output file names.",
        "Note that in all cases the checkpoint file itself is not renamed",
        "and will be overwritten, unless its name does not match",
        "the [TT]-cpo[tt] option.",
        "[PAR]",
        "With checkpointing the output is appended to previously written",
        "output files, unless [TT]-noappend[tt] is used or none of the previous",
        "output files are present (except for the checkpoint file).",
        "The integrity of the files to be appended is verified using checksums",
        "which are stored in the checkpoint file. This ensures that output can",
        "not be mixed up or corrupted due to file appending. When only some",
        "of the previous output files are present, a fatal error is generated",
        "and no old output files are modified and no new output files are opened.",
        "The result with appending will be the same as from a single run.",
        "The contents will be binary identical, unless you use a different number",
        "of ranks or dynamic load balancing or the FFT library uses optimizations",
        "through timing.",
        "[PAR]",
        "With option [TT]-maxh[tt] a simulation is terminated and a checkpoint",
        "file is written at the first neighbor search step where the run time",
        "exceeds [TT]-maxh[tt]\\*0.99 hours. This option is particularly useful in",
        "combination with setting [TT]nsteps[tt] to -1 either in the mdp or using the",
        "similarly named command line option. This results in an infinite run,",
        "terminated only when the time limit set by [TT]-maxh[tt] is reached (if any)"
        "or upon receiving a signal."
        "[PAR]",
        "When [TT]mdrun[tt] receives a TERM or INT signal (e.g. when ctrl+C is",
        "pressed), it will stop at the next neighbor search step or at the",
        "second global communication step, whichever happens later.",
        "When [TT]mdrun[tt] receives a second TERM or INT signal and",
        "reproducibility is not requested, it will stop at the first global",
        "communication step.",
        "In both cases all the usual output will be written to file and",
        "a checkpoint file is written at the last step.",
        "When [TT]mdrun[tt] receives an ABRT signal or the third TERM or INT signal,",
        "it will abort directly without writing a new checkpoint file.",
        "When running with MPI, a signal to one of the [TT]mdrun[tt] ranks",
        "is sufficient, this signal should not be sent to mpirun or",
        "the [TT]mdrun[tt] process that is the parent of the others.",
        "[PAR]",
        "Interactive molecular dynamics (IMD) can be activated by using at least one",
        "of the three IMD switches: The [TT]-imdterm[tt] switch allows one to terminate",
        "the simulation from the molecular viewer (e.g. VMD). With [TT]-imdwait[tt],",
        "[TT]mdrun[tt] pauses whenever no IMD client is connected. Pulling from the",
        "IMD remote can be turned on by [TT]-imdpull[tt].",
        "The port [TT]mdrun[tt] listens to can be altered by [TT]-imdport[tt].The",
        "file pointed to by [TT]-if[tt] contains atom indices and forces if IMD",
        "pulling is used."
        "[PAR]",
        "When [TT]mdrun[tt] is started with MPI, it does not run niced by default."
    };

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
    int               rc;

    cr = init_commrec();

    unsigned long PCA_Flags = PCA_CAN_SET_DEFFNM;
    // With -multidir, the working directory still needs to be
    // changed, so we can't check for the existence of files during
    // parsing.  It isn't useful to do any completion based on file
    // system contents, either.
    if (is_multisim_option_set(argc, argv))
    {
        PCA_Flags |= PCA_DISABLE_INPUT_FILE_CHECKING;
    }

    if (!parse_common_args(&argc, argv, PCA_Flags, nfile, fnm, asize(pa), pa,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(cr);
        return 0;
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
    gmx::ArrayRef<const std::string> multidir = opt2fnsIfOptionSet("-multidir", nfile, fnm);

    if (replExParams.exchangeInterval != 0 && multidir.size() < 2)
    {
        gmx_fatal(FARGS, "Need at least two replicas for replica exchange (use option -multidir)");
    }

    if (replExParams.numExchanges < 0)
    {
        gmx_fatal(FARGS, "Replica exchange number of exchanges needs to be positive");
    }

    ms = init_multisystem(MPI_COMM_WORLD, multidir);

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

    if (!opt2bSet("-cpi", nfile, fnm))
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

    handleRestart(cr, ms, bTryToAppendFiles, nfile, fnm, &continuationOptions.appendFiles, &continuationOptions.startedFromCheckpoint);

    mdrunOptions.rerun            = opt2bSet("-rerun", nfile, fnm);
    mdrunOptions.ntompOptionIsSet = opt2parg_bSet("-ntomp", asize(pa), pa);

    /* We postpone opening the log file if we are appending, so we can
       first truncate the old log file and append to the correct position
       there instead.  */
    if (MASTER(cr) && !continuationOptions.appendFiles)
    {
        gmx_log_open(ftp2fn(efLOG, nfile, fnm), cr,
                     continuationOptions.appendFiles, &fplog);
    }
    else
    {
        fplog = nullptr;
    }

    //////////////////////////////////////////
    domdecOptions.rankOrder    = static_cast<DdRankOrder>(nenum(ddrank_opt_choices));
    domdecOptions.dlbOption    = static_cast<DlbOption>(nenum(dddlb_opt_choices));
    domdecOptions.numCells[XX] = (int)(realddxyz[XX] + 0.5);
    domdecOptions.numCells[YY] = (int)(realddxyz[YY] + 0.5);
    domdecOptions.numCells[ZZ] = (int)(realddxyz[ZZ] + 0.5);

    nbpu_opt    = nbpu_opt_choices[0];
    pme_opt     = pme_opt_choices[0];
    pme_fft_opt = pme_fft_opt_choices[0];


    rc = mdrunner();

    /* Log file has to be closed in mdrunner if we are appending to it
       (fplog not set here) */
    if (fplog != nullptr)
    {
        gmx_log_close(fplog);
    }

    if (GMX_LIB_MPI)
    {
        done_commrec(cr);
    }
    done_multisim(ms);
    return rc;
}

//! Halt the run if there are inconsistences between user choices to run with GPUs and/or hardware detection.
static void exitIfCannotForceGpuRun(bool requirePhysicalGpu,
                                    bool emulateGpu,
                                    bool useVerletScheme,
                                    bool compatibleGpusFound)
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

    if (emulateGpu)
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

////! \brief Return the correct integrator function.
//static integrator_t *my_integrator(unsigned int ei)
//{
//    switch (ei)
//    {
//        case eiMD:
//        case eiBD:
//        case eiSD1:
//        case eiVV:
//        case eiVVAK:
//            if (!EI_DYNAMICS(ei))
//            {
//                GMX_THROW(APIError("do_md integrator would be called for a non-dynamical integrator"));
//            }
//            return do_md;
//        case eiSteep:
//            return do_steep;
//        case eiCG:
//            return do_cg;
//        case eiNM:
//            return do_nm;
//        case eiLBFGS:
//            return do_lbfgs;
//        case eiTPI:
//        case eiTPIC:
//            if (!EI_TPI(ei))
//            {
//                GMX_THROW(APIError("do_tpi integrator would be called for a non-TPI integrator"));
//            }
//            return do_tpi;
//        case eiSD2_REMOVED:
//            GMX_THROW(NotImplementedError("SD2 integrator has been removed"));
//        default:
//            GMX_THROW(APIError("Non existing integrator selected"));
//    }
//}

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
    gmx_vsite_t              *vsite            = nullptr;
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
    std::unique_ptr<gmx::MDModules> mdModules(new gmx::MDModules);
    t_inputrec                      inputrecInstance;
    t_inputrec                     *inputrec = &inputrecInstance;
    gmx_mtop_t                      mtop;

    if (mdrunOptions.continuationOptions.appendFiles)
    {
        fplog = nullptr;
    }

    bool doMembed = opt2bSet("-membed", nfile, fnm);
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
    PmeRunMode pmeRunMode      = PmeRunMode::None;

    // Here we assume that SIMMASTER(cr) does not change even after the
    // threads are started.
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
    }

    if (tpxState_ == nullptr)
    {
        // Todo: move to Mdrunner constructor
        tpxState_ = std::make_shared<TpxState>();
    }
    if (SIMMASTER(cr))
    {
        /* Only the master rank has the global state */
        globalState = compat::make_unique<t_state>();

        /* Read (nearly) all data required for the simulation */
        const char* filename = ftp2fn(efTPR, nfile, fnm);
        if (!stringIsEmpty(filename) && !tpxState_->isInitialized())
        {
            // Todo: move out of mdrunner() to a setup routine.
            // e.g.
            //    tpxState_ = TpxState::fromFile(filename); // in Mdrunner::mainFunction()
            //    ...
            //    tpxState_.extractDataLocally(); // here
            //    // or just rely on communicator to let TpxState decide what to do when getters are called.
            tpxState_ = TpxState::initializeFromFile(filename);
        }

        t_inputrec              *inputrec            = tpxState_->getRawInputrec();
        gmx_mtop_t             * mtop                = tpxState_->getRawMtop();
        t_state *                globalState         = tpxState_->getRawState();

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

            useGpuForNonbonded = decideWhetherToUseGpusForNonbondedWithThreadMpi
                    (nonbondedTarget, gpuIdsToUse, userGpuTaskAssignment, emulateGpuNonbonded,
                    inputrec->cutoff_scheme == ecutsVERLET,
                    gpuAccelerationOfNonbondedIsUseful(mdlog, inputrec, GMX_THREAD_MPI),
                    hw_opt.nthreads_tmpi);
            auto canUseGpuForPme   = pme_gpu_supports_build(nullptr) && pme_gpu_supports_input(inputrec, nullptr);
            useGpuForPme = decideWhetherToUseGpusForPmeWithThreadMpi
                    (useGpuForNonbonded, pmeTarget, gpuIdsToUse, userGpuTaskAssignment,
                    canUseGpuForPme, hw_opt.nthreads_tmpi, domdecOptions.numPmeRanks);

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
    bool useGpuForNonbonded = false;
    bool useGpuForPme       = false;
    try
    {
        // It's possible that there are different numbers of GPUs on
        // different nodes, which is the user's responsibilty to
        // handle. If unsuitable, we will notice that during task
        // assignment.
        bool gpusWereDetected = hwinfo->ngpu_compatible_tot > 0;
        useGpuForNonbonded = decideWhetherToUseGpusForNonbonded(nonbondedTarget, userGpuTaskAssignment,
                                                                emulateGpuNonbonded, inputrec->cutoff_scheme == ecutsVERLET,
                                                                gpuAccelerationOfNonbondedIsUseful(mdlog, inputrec, !GMX_THREAD_MPI),
                                                                gpusWereDetected);
        auto canUseGpuForPme   = pme_gpu_supports_build(nullptr) && pme_gpu_supports_input(inputrec, nullptr);
        useGpuForPme = decideWhetherToUseGpusForPme(useGpuForNonbonded, pmeTarget, userGpuTaskAssignment,
                                                    canUseGpuForPme, cr->nnodes, domdecOptions.numPmeRanks,
                                                    gpusWereDetected);

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
    // Currently there is at most one restraint modules.
    auto pullers = restraintManager_->getSpec();
    if (!pullers.empty())
    {
        for (auto && puller : pullers)
        {
            auto module = ::gmx::RestraintMDModule::create(puller,
                                                           puller->sites());
            mdModules.add(std::move(module));
        }
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
        cr->dd = init_domain_decomposition(fplog, cr, domdecOptions, mdrunOptions,
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

    // Enable FP exception but not in Release mode and not for compilers
    // with known buggy FP exception support (clang with any optimization)
    // or suspected buggy FP exception support (gcc 7.* with optimization).
#if !defined NDEBUG && \
    !((defined __clang__ || (defined(__GNUC__) && !defined(__ICC) && __GNUC__ == 7)) \
    && defined __OPTIMIZE__)
    const bool bEnableFPE = !EI_TPI(inputrec->eI) &&
        inputrec->cutoff_scheme == ecutsVERLET;
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
    // duty of the ranks on a node become node.
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
            if (haveGpus)
            {
                gpuTasksOnThisRank.push_back(GpuTask::Nonbonded);
            }
            else if (nonbondedTarget == TaskTarget::Gpu)
            {
                gmx_fatal(FARGS, "Cannot run short-ranged nonbonded interactions on a GPU because there is none detected.");
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
                                              mdlog, cr, ms, physicalNodeComm, gpuTasksOnThisRank);
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
            init_gpu(mdlog, nonbondedDeviceInfo);

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
        init_gpu(mdlog, pmeDeviceInfo);
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
        membed = init_membed(fplog, nfile, fnm, &mtop, inputrec, globalState.get(), cr, &mdrunOptions.checkpointOptions.period);
    }

    std::unique_ptr<MDAtoms> mdAtoms;

    snew(nrnb, 1);
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        /* Initiate forcerecord */
        fr                 = mk_forcerec();
        fr->forceProviders = mdModules.initForceProviders();
        // Threads have been launched and DD initialized
        // Todo: restraintManager can provide a proper IMDModule interface later.
//        fr->forceProviders->addForceProvider(restraintManager_);
        init_forcerec(fplog, mdlog, fr, fcd,
                      inputrec, &mtop, cr, box,
                      opt2fn("-table", nfile, fnm),
                      opt2fn("-tablep", nfile, fnm),
                      opt2fns("-tableb", nfile, fnm),
                      *hwinfo, nonbondedDeviceInfo,
                      FALSE,
                      pforce);

        /* Initialize QM-MM */
        if (fr->bQMMM)
        {
            GMX_LOG(mdlog.info).asParagraph().
                appendText("Large parts of the QM/MM support is deprecated, and may be removed in a future "
                           "version. Please get in touch with the developers if you find the support useful, "
                           "as help is needed if the functionality is to continue to be available.");
            init_QMMMrec(cr, &mtop, inputrec, fr);
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
                rvec *xGlobal = as_rvec_array(globalState->x.data());
                do_pbc_first_mtop(fplog, inputrec->ePBC, box, &mtop, xGlobal);
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
                                       mtop.natoms, nChargePerturbed, nTypePerturbed,
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
            if (EI_DYNAMICS(inputrec->eI) && MASTER(cr))
            {
                init_pull_output_files(inputrec->pull_work,
                                       nfile, fnm, oenv,
                                       continuationOptions);
            }
        }
        // If old MDP traditional MDP pulling options were used, the pull code
        // wrapped up in gmx::LegacyPullPack can be used.
        if (inputrec->bPull && inputrec->pull != nullptr)
        {
            // TODO: move to constructor when initializing runner is decoupled from reading TPR.
            /* Initialize pull code structures */
            auto pull_work =
                init_pull(fplog, inputrec->pull, inputrec, nfile, fnm,
                          mtop, cr, oenv, real(inputrec->fepvals->init_lambda),
                          EI_DYNAMICS(inputrec->eI) && MASTER(cr), Flags.to_ulong());
            auto legacyPullers = gmx::compat::make_unique<gmx::LegacyPuller>(pull_work);
            auto restraints    = gmx::restraint::Manager::instance();
            // Maybe the error is here. If the results of init_pull are different on each thread,
            // then they probably get merged accidentally here.
            restraints->add(std::move(legacyPullers), std::string("old"));
        }
        // If we need an initialization hook, we can put it here.
        //pullers_->startRun();

        std::unique_ptr<EnforcedRotation> enforcedRotation;
        if (inputrec->bRot)
        {
            /* Initialize enforced rotation code */
            enforcedRotation = init_rot(fplog, inputrec, nfile, fnm, cr, globalState.get(), &mtop, oenv, mdrunOptions);
        }

        /* Let makeConstraints know whether we have essential dynamics constraints.
         * TODO: inputrec should tell us whether we use an algorithm, not a file option or the checkpoint
         */
        bool doEssentialDynamics = (opt2fn_null("-ei", nfile, fnm) != nullptr || observablesHistory.edsamHistory);
        auto constr              = makeConstraints(mtop, *inputrec, doEssentialDynamics,
                                                   fplog, *mdAtoms->mdatoms(),
                                                   cr, *ms, nrnb, wcycle, fr->bMolPBC);

        if (DOMAINDECOMP(cr))
        {
            GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
            /* This call is not included in init_domain_decomposition mainly
             * because fr->cginfo_mb is set later.
             */
            dd_init_bondeds(fplog, cr->dd, &mtop, vsite, inputrec,
                            domdecOptions.checkBondedInteractions,
                            fr->cginfo_mb);
        }

        auto       context = gmx::md::Context(*this);
        /* Now do whatever the user wants us to do (how flexible...) */
        Integrator integrator {
            fplog, cr, ms, mdlog, nfile, fnm,
            oenv,
            mdrunOptions,
            vsite, constr.get(),
            enforcedRotation ? enforcedRotation->getLegacyEnfrot() : nullptr,
            deform.get(),
            mdModules->outputProvider(),
            inputrec, &mtop,
            fcd,
            globalState.get(),
            &observablesHistory,
            mdAtoms.get(), nrnb, wcycle, fr,
            replExParams,
            membed,
            walltime_accounting
        };
        integrator.run(inputrec->eI);

        if (inputrec->bPull)
        {
            auto puller = gmx::restraint::Manager::instance();
            puller->finish();
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

    /* Close logfile already here if we were appending to it */
    if (MASTER(cr) && continuationOptions.appendFiles)
    {
        gmx_log_close(fplog);
        fplog = nullptr;
    }

    /* Reset FPEs (important for unit tests) by disabling them. Assumes no
     * exceptions were enabled before function was called. */
    if (bEnableFPE)
    {
        gmx_fedisableexcept();
    }

    rc = (int)gmx_get_stop_condition();

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

    // If log file is open, try to flush it before we return control to the API
    if (MASTER(cr) && fplog != nullptr)
    {
        // If fplog is already closed, but has not been set to nullptr, we expect errno to be set, but we don't care,
        // so we will make sure to leave it in the same state we found it.
        const auto tempErrno = errno;
        fflush(fplog);
        errno = tempErrno;
    }

    return rc;
}

Mdrunner::Mdrunner()
{
    // Assume ownership of the Manager singleton
    restraintManager_ = ::gmx::restraint::Manager::instance();

    cr = init_commrec();
    // oenv initialized by parse_commond_args

    // dd_rank_order set according to argument processing logic (e.g. int(1))
    dd_rank_order = 1;

    // handleRestart sets &bDoAppendFiles, &bStartFromCpt

    // Flags set with lots of processing
    // Note: We cannot extract e.g. opt2parg_bSet("-append", asize(pa), pa) from this block.
    Flags.set(rerun, false);
    Flags.set(ddBondCheck, true);
    Flags.set(ddBondComm, true);
    Flags.set(tunePME, true);
    Flags.set(confOut, true);
    Flags.set(rerunVSite, false);
    Flags.set(reproducible, false);
    Flags.set(appendFiles, false);
    //Flags = Flags | (opt2parg_bSet("-append", asize(pa), pa) ? MD_APPENDFILESSET : 0);
    Flags.set(appendFilesSet, false);
    Flags.set(keepAndNumCpt, false);
    Flags.set(startFromCpt, false);
    Flags.set(resetCountersHalfWay, false);
    //Flags = Flags | (opt2parg_bSet("-ntomp", asize(pa), pa) ? MD_NTOMPSET : 0);
    Flags.set(ntompSet, false);
    Flags.set(imdWait, false);
    Flags.set(imdTerm, false);
    Flags.set(imdPull, false);

    // log opened to fplog if MASTER(cr) && !bDoAppendFiles

    // dddlb_opt set from processed options (e.g. const char* "auto")
    dddlb_opt = "auto";
    // nbpu_opt set from processed options (e.g. const char* "auto")
    nbpu_opt = "auto";
};

Mdrunner::~Mdrunner()
{
    // Clean up of the Manager singleton.
    // This will end up getting called on every thread-MPI rank, which is okay, but unnecessary. There should probably
    // be a simulation shutdown hook and this manager probably shouldn't be a singleton.
    restraintManager_->clear();
    assert(restraintManager_->countRestraints() == 0);

    /* Log file has to be closed in mdrunner if we are appending to it
       (fplog not set here) */
    // assert(cr != nullptr); // Todo: can we just initialize the cr in the constructor and keep it initialized?
    if (cr != nullptr && MASTER(cr) && !(Flags.test(appendFiles)))
    {
        gmx_log_close(fplog);
        fplog = nullptr;
    }
    sfree(cr);
}

void Mdrunner::setTpx(std::shared_ptr<gmx::TpxState> newState)
{
    if (newState->isDirty())
    {
        GMX_THROW(gmx::InvalidInputError("Attempting to assign from a dirty state."));
    }
    // No good way to lock with default constructor and default moves for Mdrunner.
    // std::lock_guard<std::mutex> lock(stateAccess);
    // Todo: thread-safety
    // Locking the gmx::Mdrunner to serialize state updates would be nice, but
    // it would be sufficient to guarantee that the gmx::Mdrunner is thread-local
    //assert(tpxState_ != nullptr); // It is probably preferable to be able to assume tpxState_ is valid even if empty.
    if (tpxState_ != nullptr && tpxState_->isDirty())
    {
        // Calling code has a logic error: the old state is in use somewhere.
        GMX_THROW(gmx::APIError("Attempting to replace a state that may be in use (isDirty() == true)"));
    }
    tpxState_ = std::move(newState);
}

void Mdrunner::addPullPotential(std::shared_ptr<gmx::IRestraintPotential> puller,
                                std::string                               name)
{
    assert(restraintManager_ != nullptr);
    std::cout << "Registering restraint named " << name << std::endl;

    // When multiple restraints are used, it may be wasteful to register them separately.
    // Maybe instead register a Restraint Manager as a force provider.
    restraintManager_->addToSpec(std::move(puller),
                                 std::move(name));
}

void Mdrunner::declareFinalStep()
{
    simulationSignals_[eglsSTOPCOND].sig = true;
}

SimulationSignals *Mdrunner::signals() const
{
    return &simulationSignals_;
}

Mdrunner &Mdrunner::operator=(Mdrunner &&) = default;

Mdrunner::Mdrunner(Mdrunner &&) = default;

} // namespace gmx
