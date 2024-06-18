/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include <optional>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/builder.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/domdec/makebondedlinks.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/domdec/reversetopology.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_gpu_program.h"
#include "gromacs/ewald/pme_only.h"
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpueventsynchronizer_helpers.h"
#include "gromacs/gpu_utils/nvshmem_utils.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hardwaretopology.h"
#include "gromacs/hardware/printhardware.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/disre.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
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
#include "gromacs/mdlib/gpuforcereduction.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdgraph_gpu.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdrun/simulationcontext.h"
#include "gromacs/mdrun/simulationinput.h"
#include "gromacs/mdrun/simulationinputhandle.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/logging.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/mdtypes/checkpointdata.h"
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
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/modularsimulator/modularsimulator.h"
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
#include "gromacs/topology/topology.h"
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
#include "gromacs/utility/mpiinfo.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "isimulator.h"
#include "membedholder.h"
#include "replicaexchange.h"
#include "simulatorbuilder.h"

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
 * \param[in]  pmeRunMode   Run mode indicating what resource is PME executed on.
 * \param[in]  numRanksPerSimulation   The number of ranks in each simulation.
 * \param[in]  numPmeRanksPerSimulation   The number of PME ranks in each simulation, can be -1
 * \param[in]  gpuAwareMpiStatus  Minimum level of GPU-aware MPI support across all ranks
 * \returns                         The object populated with development feature flags.
 */
static DevelopmentFeatureFlags manageDevelopmentFeatures(const gmx::MDLogger& mdlog,
                                                         const bool           useGpuForNonbonded,
                                                         const PmeRunMode     pmeRunMode,
                                                         const int            numRanksPerSimulation,
                                                         const int numPmeRanksPerSimulation,
                                                         gmx::GpuAwareMpiStatus gpuAwareMpiStatus)
{
    DevelopmentFeatureFlags devFlags;

    devFlags.enableGpuBufferOps = (GMX_GPU_CUDA || GMX_GPU_SYCL) && useGpuForNonbonded
                                  && (getenv("GMX_USE_GPU_BUFFER_OPS") != nullptr);

    if (getenv("GMX_CUDA_GRAPH") != nullptr)
    {
        if (GMX_HAVE_GPU_GRAPH_SUPPORT)
        {
            devFlags.enableCudaGraphs = true;
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GMX_CUDA_GRAPH environment variable is detected. "
                            "The experimental CUDA Graphs feature will be used if run conditions "
                            "allow.");
        }
        else
        {
            devFlags.enableCudaGraphs = false;
            std::string errorReason;
            if (GMX_GPU_CUDA)
            {
                errorReason = "the CUDA version in use is below the minimum requirement (11.1)";
            }
            else
            {
                errorReason = "GROMACS is built without CUDA";
            }
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "GMX_CUDA_GRAPH environment variable is detected, but %s. GPU Graphs "
                            "will be disabled.",
                            errorReason.c_str());
        }
    }

    // Flag use to enable GPU-aware MPI depenendent features such PME GPU decomposition
    // GPU-aware MPI is marked available if it has been detected by GROMACS or detection fails but
    // user wants to force its use
    devFlags.canUseGpuAwareMpi = false;

    // Direct GPU comm path is being used with GPU-aware MPI
    // make sure underlying MPI implementation is GPU-aware

    if (GMX_LIB_MPI && (GMX_GPU_CUDA || GMX_GPU_SYCL))
    {
        // Allow overriding the detection for GPU-aware MPI
        if (getenv("GMX_FORCE_CUDA_AWARE_MPI") != nullptr)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GMX_FORCE_CUDA_AWARE_MPI environment variable is inactive. "
                            "Please use GMX_FORCE_GPU_AWARE_MPI instead.");
        }

        devFlags.canUseGpuAwareMpi = (gpuAwareMpiStatus == gmx::GpuAwareMpiStatus::Supported
                                      || gpuAwareMpiStatus == gmx::GpuAwareMpiStatus::Forced);
        if (getenv("GMX_ENABLE_DIRECT_GPU_COMM") != nullptr)
        {
            if (gpuAwareMpiStatus == gmx::GpuAwareMpiStatus::Forced)
            {
                // GPU-aware support not detected in MPI library but, user has forced its use
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "This run has forced use of 'GPU-aware MPI'. "
                                "However, GROMACS cannot determine if underlying MPI is GPU-aware. "
                                "Check the GROMACS install guide for recommendations for GPU-aware "
                                "support. If you observe failures at runtime, try unsetting the "
                                "GMX_FORCE_GPU_AWARE_MPI environment variable.");
            }

            if (devFlags.canUseGpuAwareMpi)
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "GMX_ENABLE_DIRECT_GPU_COMM environment variable detected, "
                                "enabling direct GPU communication using GPU-aware MPI.");
            }
            else
            {
                GMX_LOG(mdlog.warning)
                        .asParagraph()
                        .appendText(
                                "GPU-aware MPI was not detected, will not use direct GPU "
                                "communication. Check the GROMACS install guide for "
                                "recommendations "
                                "for GPU-aware support. If you are certain about GPU-aware support "
                                "in your MPI library, you can force its use by setting the "
                                "GMX_FORCE_GPU_AWARE_MPI environment variable.");
            }
        }
        else if (gpuAwareMpiStatus == gmx::GpuAwareMpiStatus::Supported)
        {
            // GPU-aware MPI was detected, let the user know that using it may improve performance
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GPU-aware MPI detected, but by default GROMACS will not "
                            "make use the direct GPU communication capabilities of MPI. "
                            "For improved performance try enabling the feature by setting "
                            "the GMX_ENABLE_DIRECT_GPU_COMM environment variable.");
        }
    }
    else
    {
        if (getenv("GMX_FORCE_GPU_AWARE_MPI") != nullptr)
        {
            // Cannot force use of GPU-aware MPI in this build configuration
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText(
                            "A CUDA or SYCL build with an external MPI library is required in "
                            "order to benefit from GMX_FORCE_GPU_AWARE_MPI. That environment "
                            "variable is being ignored because such a build is not in use.");
        }
    }

    if (getenv("GMX_ENABLE_NVSHMEM") != nullptr)
    {
        if (GMX_LIB_MPI && GMX_NVSHMEM)
        {
            devFlags.enableNvshmem = true;
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GMX_ENABLE_NVSHMEM environment variable is detected. "
                            "The experimental NVSHMEM feature will be used if run conditions "
                            "allow.");
        }
        else
        {
            devFlags.enableNvshmem = false;
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendText(
                            "GMX_ENABLE_NVSHMEM environment variable is detected, "
                            "but GROMACS was built without NVSHMEM support. "
                            "Direct use of NVSHMEM will be disabled. "
                            "NVSHMEM may still be used indirectly if cuFFTMp is enabled. ");
        }
    }

    if (devFlags.enableGpuBufferOps)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted(
                        "This run uses the 'GPU buffer ops' feature, enabled by the "
                        "GMX_USE_GPU_BUFFER_OPS environment variable.");
    }

    // PME decomposition is supported only with CUDA or SYCL and also
    // needs GPU-aware MPI support for it to work.
    const bool pmeGpuDecompositionRequested =
            (pmeRunMode == PmeRunMode::GPU || pmeRunMode == PmeRunMode::Mixed)
            && ((numRanksPerSimulation > 1 && numPmeRanksPerSimulation == 0)
                || numPmeRanksPerSimulation > 1);
    const bool pmeGpuDecompositionSupported =
            (devFlags.canUseGpuAwareMpi && (GMX_GPU_CUDA || GMX_GPU_SYCL)
             && ((pmeRunMode == PmeRunMode::GPU && (GMX_USE_Heffte || GMX_USE_cuFFTMp))
                 || pmeRunMode == PmeRunMode::Mixed));

    const bool forcePmeGpuDecomposition = getenv("GMX_GPU_PME_DECOMPOSITION") != nullptr;

    if (pmeGpuDecompositionSupported && pmeGpuDecompositionRequested)
    {
        // PME decomposition is supported only when it is forced using GMX_GPU_PME_DECOMPOSITION
        if (forcePmeGpuDecomposition)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted(
                            "This run has requested the 'GPU PME decomposition' feature, enabled "
                            "by the GMX_GPU_PME_DECOMPOSITION environment variable. "
                            "PME decomposition lacks substantial testing "
                            "and should be used with caution.");
        }
        else
        {
            gmx_fatal(FARGS,
                      "Multiple PME tasks were required to run on GPUs, "
                      "but that is not supported. "
                      "Use GMX_GPU_PME_DECOMPOSITION environment variable to enable it.");
        }
    }

    if (!pmeGpuDecompositionSupported && pmeGpuDecompositionRequested)
    {
        if (GMX_GPU_CUDA)
        {
            gmx_fatal(FARGS,
                      "PME tasks were required to run on more than one CUDA-devices. To enable "
                      "this feature, "
                      "use MPI with CUDA-aware support and build GROMACS with cuFFTMp support.");
        }
        else
        {
            gmx_fatal(
                    FARGS,
                    "PME tasks were required to run on GPUs, but that is not implemented with "
                    "more than one PME rank. Use a single rank simulation, or a separate PME rank, "
                    "or permit PME tasks to be assigned to the CPU.");
        }
    }

    devFlags.enableGpuPmeDecomposition =
            forcePmeGpuDecomposition && pmeGpuDecompositionRequested && pmeGpuDecompositionSupported;

    return devFlags;
}

/*! \brief Barrier for safe simultaneous thread access to mdrunner data
 *
 * Used to ensure that the main thread does not modify mdrunner during copy
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

    // Copy members of main runner.
    // \todo Replace with builder when Simulation context and/or runner phases are better defined.
    // Ref https://gitlab.com/gromacs/gromacs/-/issues/2587 and https://gitlab.com/gromacs/gromacs/-/issues/2375
    newRunner.hw_opt    = hw_opt;
    newRunner.filenames = filenames;

    newRunner.hwinfo_         = hwinfo_;
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
    newRunner.libraryWorldCommunicator = MPI_COMM_WORLD;
    newRunner.simulationCommunicator   = MPI_COMM_WORLD;
    newRunner.ms                       = ms;
    newRunner.startingBehavior         = startingBehavior;
    newRunner.stopHandlerBuilder_      = std::make_unique<StopHandlerBuilder>(*stopHandlerBuilder_);
    newRunner.inputHolder_             = inputHolder_;

    threadMpiMdrunnerAccessBarrier();

    return newRunner;
}

/*! \brief The callback used for running on spawned threads.
 *
 * Obtains the pointer to the main mdrunner object from the one
 * argument permitted to the thread-launch API call, copies it to make
 * a new runner for this thread, reinitializes necessary data, and
 * proceeds to the simulation. */
static void mdrunner_start_fn(const void* arg)
{
    try
    {
        const auto* mainMdrunner = reinterpret_cast<const gmx::Mdrunner*>(arg);
        /* copy the arg list to make sure that it's thread-local. This
           doesn't copy pointed-to items, of course; fnm, cr and fplog
           are reset in the call below, all others should be const. */
        gmx::Mdrunner mdrunner = mainMdrunner->cloneOnSpawnedThread();
        mdrunner.mdrunner();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}


void Mdrunner::spawnThreads(int numThreadsToLaunch)
{
#if GMX_THREAD_MPI
    /* now spawn new threads that start mdrunner_start_fn(), while
       the main thread returns. Thread affinity is handled later. */
    if (tMPI_Init_fn(TRUE, numThreadsToLaunch, TMPI_AFFINITY_NONE, mdrunner_start_fn, static_cast<const void*>(this))
        != TMPI_SUCCESS)
    {
        GMX_THROW(gmx::InternalError("Failed to spawn thread-MPI threads"));
    }

    // Give the main thread the newly created valid communicator for
    // the simulation.
    libraryWorldCommunicator = MPI_COMM_WORLD;
    simulationCommunicator   = MPI_COMM_WORLD;
    threadMpiMdrunnerAccessBarrier();
#else
    GMX_UNUSED_VALUE(numThreadsToLaunch);
    GMX_UNUSED_VALUE(mdrunner_start_fn);
#endif
}

} // namespace gmx

/*! \brief Initialize variables for Verlet scheme simulation */
static void prepare_verlet_scheme(FILE*                          fplog,
                                  t_commrec*                     cr,
                                  t_inputrec*                    ir,
                                  int                            nstlist_cmdline,
                                  const gmx_mtop_t&              mtop,
                                  gmx::ArrayRef<const gmx::RVec> coordinates,
                                  const matrix                   box,
                                  bool                           makeGpuPairList,
                                  const gmx::CpuInfo&            cpuinfo)
{
    // We checked the cut-offs in grompp, but double-check here.
    // We have PME+LJcutoff kernels for rcoulomb>rvdw.
    if (usingPmeOrEwald(ir->coulombtype) && ir->vdwtype == VanDerWaalsType::Cut)
    {
        GMX_RELEASE_ASSERT(ir->rcoulomb >= ir->rvdw,
                           "With Verlet lists and PME we should have rcoulomb>=rvdw");
    }
    else
    {
        GMX_RELEASE_ASSERT(ir->rcoulomb == ir->rvdw,
                           "With Verlet lists and no PME rcoulomb and rvdw should be identical");
    }

    std::optional<real> effectiveAtomDensity;
    if (EI_DYNAMICS(ir->eI))
    {
        effectiveAtomDensity = computeEffectiveAtomDensity(
                coordinates, box, std::max(ir->rcoulomb, ir->rvdw), cr->mpiDefaultCommunicator);
    }

    /* For NVE simulations, we will retain the initial list buffer */
    if (EI_DYNAMICS(ir->eI) && ir->verletbuf_tol > 0
        && !(EI_MD(ir->eI) && ir->etc == TemperatureCoupling::No))
    {
        /* Update the Verlet buffer size for the current run setup */

        /* Here we assume SIMD-enabled kernels are being used. But as currently
         * calc_verlet_buffer_size gives the same results for 4x8 and 4x4
         * and 4x2 gives a larger buffer than 4x4, this is ok.
         */
        ListSetupType listType =
                (makeGpuPairList ? ListSetupType::Gpu : ListSetupType::CpuSimdWhenSupported);
        VerletbufListSetup listSetup = verletbufGetSafeListSetup(listType);

        const real rlist_new = calcVerletBufferSize(mtop,
                                                    effectiveAtomDensity.value(),
                                                    *ir,
                                                    ir->verletBufferPressureTolerance,
                                                    ir->nstlist,
                                                    ir->nstlist - 1,
                                                    -1,
                                                    listSetup);

        if (rlist_new != ir->rlist)
        {
            if (fplog != nullptr)
            {
                fprintf(fplog,
                        "\nChanging rlist from %g to %g for non-bonded %dx%d atom kernels\n\n",
                        ir->rlist,
                        rlist_new,
                        listSetup.cluster_size_i,
                        listSetup.cluster_size_j);
            }
            ir->rlist = rlist_new;
        }
    }

    if (nstlist_cmdline > 0 && (!EI_DYNAMICS(ir->eI) || ir->verletbuf_tol <= 0))
    {
        gmx_fatal(FARGS,
                  "Can not set nstlist without %s",
                  !EI_DYNAMICS(ir->eI) ? "dynamics" : "verlet-buffer-tolerance");
    }

    if (EI_DYNAMICS(ir->eI))
    {
        /* Set or try nstlist values */
        increaseNstlist(
                fplog, cr, ir, nstlist_cmdline, &mtop, box, effectiveAtomDensity.value(), makeGpuPairList, cpuinfo);
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
                    gmx_step_str(nsteps_cmdline, sbuf_steps),
                    fabs(nsteps_cmdline * ir->delta_t));
        }
        else
        {
            sprintf(sbuf_msg,
                    "Overriding nsteps with value passed on the command line: %s steps",
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
static bool gpuAccelerationOfNonbondedIsUseful(const MDLogger&   mdlog,
                                               const t_inputrec& ir,
                                               const bool        issueWarning,
                                               const bool        doRerun)
{
    bool        gpuIsUseful = true;
    std::string warning;

    if (ir.opts.ngener - ir.nwall > 1)
    {
        /* The GPU code does not support more than one energy group.
         * If the user requested GPUs explicitly, a fatal error is given later.
         */
        gpuIsUseful = false;
        if (!doRerun)
        {
            warning =
                    "Multiple energy groups is not implemented for GPUs, falling back to the CPU. "
                    "For better performance, run on the GPU without energy groups and then do "
                    "gmx mdrun -rerun option on the trajectory with an energy group .tpr file.";
        }
    }

    /* There are resource handling issues in the GPU code paths with MTS on anything else than only
     * PME. Also those code paths need more testing.
     */
    MtsLevel mtsLevelOnlyPme;
    mtsLevelOnlyPme.forceGroups.set(static_cast<int>(MtsForceGroups::LongrangeNonbonded));
    if (ir.useMts && !(ir.mtsLevels.size() == 2 && ir.mtsLevels[1].forceGroups == mtsLevelOnlyPme.forceGroups))
    {
        gpuIsUseful = false;
        warning     = gmx::formatString(
                "Multiple time stepping is only supported with GPUs when MTS is only applied to %s "
                "forces.",
                mtsForceGroupNames[MtsForceGroups::LongrangeNonbonded].c_str());
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
static gmx::LoggerOwner buildLogger(FILE* fplog, const bool isSimulationMainRank)
{
    gmx::LoggerBuilder builder;
    if (fplog != nullptr)
    {
        builder.addTargetFile(gmx::MDLogger::LogLevel::Info, fplog);
    }
    if (isSimulationMainRank)
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
                       const t_inputrec&         inputrec,
                       t_nrnb                    nrnb[],
                       gmx_wallcycle*            wcycle,
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
       simulation main may print, but it should not do so if the run
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
    bool printReport = EI_DYNAMICS(inputrec.eI) && SIMMAIN(cr);

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
        MPI_Allreduce(nrnb->n.data(), nrnb_tot->n.data(), eNRNB, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mysim);
#endif
    }
    else
    {
        nrnb_tot = nrnb;
    }

    elapsed_time = walltime_accounting_get_time_since_reset(walltime_accounting);
    elapsed_time_over_all_threads =
            walltime_accounting_get_time_since_reset_over_all_threads(walltime_accounting);
    if (GMX_MPI && cr->nnodes > 1)
    {
#if GMX_MPI
        /* reduce elapsed_time over all MPI ranks in the current simulation */
        MPI_Allreduce(&elapsed_time, &elapsed_time_over_all_ranks, 1, MPI_DOUBLE, MPI_SUM, cr->mpi_comm_mysim);
        elapsed_time_over_all_ranks /= cr->nnodes;
        /* Reduce elapsed_time_over_all_threads over all MPI ranks in the
         * current simulation. */
        MPI_Allreduce(&elapsed_time_over_all_threads,
                      &elapsed_time_over_all_threads_over_all_ranks,
                      1,
                      MPI_DOUBLE,
                      MPI_SUM,
                      cr->mpi_comm_mysim);
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

    if (thisRankHasDuty(cr, DUTY_PP) && haveDDAtomOrdering(*cr))
    {
        print_dd_statistics(cr, inputrec, fplog);
    }

    /* TODO Move the responsibility for any scaling by thread counts
     * to the code that handled the thread region, so that there's a
     * mechanism to keep cycle counting working during the transition
     * to task parallelism. */
    int nthreads_pp  = gmx_omp_nthreads_get(ModuleMultiThread::Nonbonded);
    int nthreads_pme = gmx_omp_nthreads_get(ModuleMultiThread::Pme);
    wallcycle_scale_by_num_threads(
            wcycle, thisRankHasDuty(cr, DUTY_PME) && !thisRankHasDuty(cr, DUTY_PP), nthreads_pp, nthreads_pme);
    auto cycle_sum(wallcycle_sum(cr, wcycle));

    if (printReport)
    {
        auto* nbnxn_gpu_timings =
                (nbv != nullptr && nbv->useGpu()) ? Nbnxm::gpu_get_timings(nbv->gpuNbv()) : nullptr;
        gmx_wallclock_gpu_pme_t pme_gpu_timings = {};

        if (pme_gpu_task_enabled(pme))
        {
            pme_gpu_get_timings(pme, &pme_gpu_timings);
        }
        wallcycle_print(fplog,
                        mdlog,
                        cr->nnodes,
                        cr->npmenodes,
                        nthreads_pp,
                        nthreads_pme,
                        elapsed_time_over_all_ranks,
                        wcycle,
                        cycle_sum,
                        nbnxn_gpu_timings,
                        &pme_gpu_timings);

        if (EI_DYNAMICS(inputrec.eI))
        {
            delta_t = inputrec.delta_t;
        }

        if (fplog)
        {
            print_perf(fplog,
                       elapsed_time_over_all_threads_over_all_ranks,
                       elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done_since_reset(walltime_accounting),
                       delta_t,
                       nbfs,
                       mflop);
        }
        if (bWriteStat)
        {
            print_perf(stderr,
                       elapsed_time_over_all_threads_over_all_ranks,
                       elapsed_time_over_all_ranks,
                       walltime_accounting_get_nsteps_done_since_reset(walltime_accounting),
                       delta_t,
                       nbfs,
                       mflop);
        }
    }
}

int Mdrunner::mdrunner()
{
    std::unique_ptr<t_forcerec> fr;
    real                        ewaldcoeff_q     = 0;
    real                        ewaldcoeff_lj    = 0;
    int                         nChargePerturbed = -1, nTypePerturbed = 0;
    gmx_walltime_accounting_t   walltime_accounting = nullptr;
    MembedHolder                membedHolder(filenames.size(), filenames.data());

    /* CAUTION: threads may be started later on in this function, so
       cr doesn't reflect the final parallel state right now */
    gmx_mtop_t mtop;

    /* TODO: inputrec should tell us whether we use an algorithm, not a file option */
    const bool doEssentialDynamics = opt2bSet("-ei", filenames.size(), filenames.data());
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
    const bool       isSimulationMainRank = findIsSimulationMainRank(ms, simulationCommunicator);
    gmx::LoggerOwner logOwner(buildLogger(fplog, isSimulationMainRank));
    gmx::MDLogger    mdlog(logOwner.logger());

    gmx_print_detected_hardware(fplog, isSimulationMainRank && isMainSim(ms), mdlog, hwinfo_);

    std::vector<int> availableDevices =
            makeListOfAvailableDevices(hwinfo_->deviceInfoList, hw_opt.devicesSelectedByUser);
    const int numAvailableDevices = gmx::ssize(availableDevices);

    // Print citation requests after all software/hardware printing
    pleaseCiteGromacs(fplog);

    // Note: legacy program logic relies on checking whether these pointers are assigned.
    // Objects may or may not be allocated later.
    std::unique_ptr<t_inputrec> inputrec;
    std::unique_ptr<t_state>    globalState;

    auto partialDeserializedTpr = std::make_unique<PartialDeserializedTprFile>();

    if (isSimulationMainRank)
    {
        // Allocate objects to be initialized by later function calls.
        /* Only the main rank has the global state */
        globalState = std::make_unique<t_state>();
        inputrec    = std::make_unique<t_inputrec>();

        /* Read (nearly) all data required for the simulation
         * and keep the partly serialized tpr contents to send to other ranks later
         */
        applyGlobalSimulationState(
                *inputHolder_.get(), partialDeserializedTpr.get(), globalState.get(), inputrec.get(), &mtop);

        static_assert(sc_trrMaxAtomCount == sc_checkpointMaxAtomCount);
        if (mtop.natoms > sc_checkpointMaxAtomCount)
        {
            gmx_fatal(FARGS,
                      "System has %d atoms, which is more than can be stored in checkpoint and trr "
                      "files (max %" PRId64 ")",
                      mtop.natoms,
                      sc_checkpointMaxAtomCount);
        }

        // The XTC format has been updated to support up to 2^31-1 atoms, which is anyway the
        // largest supported by GROMACS, so no need for any particular check here.
    }

    /* Check and update the hardware options for internal consistency */
    checkAndUpdateHardwareOptions(
            mdlog, &hw_opt, isSimulationMainRank, domdecOptions.numPmeRanks, inputrec.get());

    if (GMX_THREAD_MPI && isSimulationMainRank)
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
                    nonbondedTarget,
                    numAvailableDevices > 0,
                    userGpuTaskAssignment,
                    emulateGpuNonbonded,
                    canUseGpuForNonbonded,
                    gpuAccelerationOfNonbondedIsUseful(mdlog, *inputrec, GMX_THREAD_MPI, doRerun),
                    mdrunOptions.reproducible,
                    hw_opt.nthreads_tmpi);
            useGpuForPme = decideWhetherToUseGpusForPmeWithThreadMpi(useGpuForNonbonded,
                                                                     pmeTarget,
                                                                     pmeFftTarget,
                                                                     numAvailableDevices,
                                                                     userGpuTaskAssignment,
                                                                     *inputrec,
                                                                     hw_opt.nthreads_tmpi,
                                                                     domdecOptions.numPmeRanks);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

        /* Determine how many thread-MPI ranks to start.
         *
         * TODO Over-writing the user-supplied value here does
         * prevent any possible subsequent checks from working
         * correctly. */
        hw_opt.nthreads_tmpi = get_nthreads_mpi(hwinfo_,
                                                &hw_opt,
                                                numAvailableDevices,
                                                useGpuForNonbonded,
                                                useGpuForPme,
                                                inputrec.get(),
                                                mtop,
                                                mdlog,
                                                membedHolder.doMembed());

        // Now start the threads for thread MPI.
        spawnThreads(hw_opt.nthreads_tmpi);
        // The spawned threads enter mdrunner() and execution of
        // main and spawned threads joins at the end of this block.
    }

    GMX_RELEASE_ASSERT(!GMX_MPI || ms || simulationCommunicator != MPI_COMM_NULL,
                       "Must have valid communicator unless running a multi-simulation");
    std::unique_ptr<t_commrec> crHandle(init_commrec(simulationCommunicator));
    t_commrec*                 cr = crHandle.get();
    GMX_RELEASE_ASSERT(cr != nullptr, "Must have valid commrec");

    PhysicalNodeCommunicator physicalNodeComm(libraryWorldCommunicator, gmx_physicalnode_id_hash());

    if (PAR(cr))
    {
        /* now broadcast everything to the non-main nodes/threads: */
        if (!isSimulationMainRank)
        {
            // Until now, only the main rank has a non-null pointer.
            // On non-main ranks, allocate the object that will receive data in the following call.
            inputrec = std::make_unique<t_inputrec>();
        }
        init_parallel(
                cr->mpiDefaultCommunicator, MAIN(cr), inputrec.get(), &mtop, partialDeserializedTpr.get());
    }
    GMX_RELEASE_ASSERT(inputrec != nullptr, "All ranks should have a valid inputrec now");
    partialDeserializedTpr.reset(nullptr);

    // Note that these variables describe only their own node.
    //
    // Note that when bonded interactions run on a GPU they always run
    // alongside a nonbonded task, so do not influence task assignment
    // even though they affect the force calculation workload.
    bool useGpuForNonbonded = false;
    bool useGpuForPme       = false;
    bool useGpuForBonded    = false;
    bool useGpuForUpdate    = false;
    bool gpusWereDetected   = hwinfo_->ngpu_compatible_tot > 0;
    try
    {
        // It's possible that there are different numbers of GPUs on
        // different nodes, which is the user's responsibility to
        // handle. If unsuitable, we will notice that during task
        // assignment.
        auto canUseGpuForNonbonded = buildSupportsNonbondedOnGpu(nullptr);
        useGpuForNonbonded         = decideWhetherToUseGpusForNonbonded(
                nonbondedTarget,
                userGpuTaskAssignment,
                emulateGpuNonbonded,
                canUseGpuForNonbonded,
                gpuAccelerationOfNonbondedIsUseful(mdlog, *inputrec, !GMX_THREAD_MPI, doRerun),
                mdrunOptions.reproducible,
                gpusWereDetected);
        useGpuForPme    = decideWhetherToUseGpusForPme(useGpuForNonbonded,
                                                    pmeTarget,
                                                    pmeFftTarget,
                                                    userGpuTaskAssignment,
                                                    *inputrec,
                                                    cr->sizeOfDefaultCommunicator,
                                                    domdecOptions.numPmeRanks,
                                                    gpusWereDetected);
        useGpuForBonded = decideWhetherToUseGpusForBonded(
                useGpuForNonbonded, useGpuForPme, bondedTarget, *inputrec, mtop, domdecOptions.numPmeRanks, gpusWereDetected);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    const PmeRunMode pmeRunMode = determinePmeRunMode(useGpuForPme, pmeFftTarget, *inputrec);

    // Initialize development feature flags that enabled by environment variable
    // and report those features that are enabled.
    // We are using the minimal supported level of GPU-aware MPI
    // support as an approximation of whether such communication
    // will work. It likely would not work in cases where ranks
    // have heterogeneous device types or vendors unless the MPI
    // library supported that.
    const DevelopmentFeatureFlags devFlags = manageDevelopmentFeatures(mdlog,
                                                                       useGpuForNonbonded,
                                                                       pmeRunMode,
                                                                       cr->sizeOfDefaultCommunicator,
                                                                       domdecOptions.numPmeRanks,
                                                                       hwinfo_->minGpuAwareMpiStatus);

    const bool useModularSimulator = checkUseModularSimulator(false,
                                                              inputrec.get(),
                                                              doRerun,
                                                              mtop,
                                                              ms,
                                                              replExParams,
                                                              nullptr,
                                                              doEssentialDynamics,
                                                              membedHolder.doMembed(),
                                                              updateTarget == TaskTarget::Gpu);

    // Now the number of ranks is known to all ranks, and each knows
    // the inputrec read by the main rank. The ranks can now all run
    // the task-deciding functions and will agree on the result
    // without needing to communicate.
    // The LBFGS minimizer, test-particle insertion, normal modes and shell dynamics don't support DD
    const bool hasCustomParallelization =
            (EI_TPI(inputrec->eI) || inputrec->eI == IntegrationAlgorithm::NM);
    const bool canUseDomainDecomposition =
            (inputrec->eI != IntegrationAlgorithm::LBFGS && !hasCustomParallelization
             && gmx_mtop_particletype_count(mtop)[ParticleType::Shell] == 0);
    GMX_RELEASE_ASSERT(!PAR(cr) || hasCustomParallelization || canUseDomainDecomposition,
                       "A parallel run should not arrive here without DD support");

    int useDDWithSingleRank = -1;
    if (const char* ddSingleRankEnv = getenv("GMX_DD_SINGLE_RANK"))
    {
        useDDWithSingleRank = std::strtol(ddSingleRankEnv, nullptr, 10);
    }

    // The overhead of DD partitioning is only compensated when we have both non-bondeds and PME on the CPU
    const bool useDomainDecomposition =
            canUseDomainDecomposition
            && (PAR(cr)
                || (!useGpuForNonbonded && usingFullElectrostatics(inputrec->coulombtype)
                    && useDDWithSingleRank != 0)
                || useDDWithSingleRank == 1);

    ObservablesReducerBuilder observablesReducerBuilder;

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
    // now that the MDModules know their options, they know which callbacks to sign up to
    mdModules_->subscribeToSimulationSetupNotifications();
    const auto& setupNotifier = mdModules_->notifiers().simulationSetupNotifier_;

    // Notify MdModules of existing logger
    setupNotifier.notify(mdlog);

    // Notify MdModules of internal parameters, saved into KVT
    if (inputrec->internalParameters != nullptr)
    {
        setupNotifier.notify(*inputrec->internalParameters);
    }

    // Let MdModules know the .tpr input and .edr output filenames
    {
        gmx::MdRunInputFilename mdRunInputFilename = { ftp2fn(efTPR, filenames.size(), filenames.data()) };
        setupNotifier.notify(mdRunInputFilename);
        gmx::EdrOutputFilename edrOutputFilename = { ftp2fn(efEDR, filenames.size(), filenames.data()) };
        setupNotifier.notify(edrOutputFilename);
    }

    if (fplog != nullptr)
    {
        pr_inputrec(fplog, 0, "Input Parameters", inputrec.get(), FALSE);
        fprintf(fplog, "\n");
    }

    if (SIMMAIN(cr))
    {
        /* In rerun, set velocities to zero if present */
        if (doRerun && globalState->hasEntry(StateEntry::V))
        {
            // rerun does not use velocities
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText(
                            "Rerun trajectory contains velocities. Rerun does only evaluate "
                            "potential energy and forces. The velocities will be ignored.");
            for (int i = 0; i < globalState->numAtoms(); i++)
            {
                clear_rvec(globalState->v[i]);
            }
            globalState->setFlags(globalState->flags() & ~enumValueToBitMask(StateEntry::V));
        }

        /* now make sure the state is initialized and propagated */
        set_state_entries(globalState.get(), inputrec.get(), useModularSimulator);
    }

    /* NM and TPI parallelize over force/energy calculations, not atoms,
     * so we need to initialize and broadcast the global state.
     */
    if (inputrec->eI == IntegrationAlgorithm::NM || inputrec->eI == IntegrationAlgorithm::TPI)
    {
        if (!MAIN(cr))
        {
            globalState = std::make_unique<t_state>();
        }
        broadcastStateWithoutDynamics(
                cr->mpiDefaultCommunicator, haveDDAtomOrdering(*cr), PAR(cr), globalState.get());
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

    if (doRerun && (EI_ENERGY_MINIMIZATION(inputrec->eI) || IntegrationAlgorithm::NM == inputrec->eI))
    {
        gmx_fatal(FARGS,
                  "The .mdp file specified an energy mininization or normal mode algorithm, and "
                  "these are not compatible with mdrun -rerun");
    }

    /* NMR restraints must be initialized before load_checkpoint,
     * since with time averaging the history is added to t_state.
     * For proper consistency check we therefore need to extend
     * t_state here.
     * So the PME-only nodes (if present) will also initialize
     * the distance restraints.
     */

    /* This needs to be called before read_checkpoint to extend the state */
    t_disresdata* disresdata;
    snew(disresdata, 1);
    init_disres(fplog,
                mtop,
                inputrec.get(),
                DisResRunMode::MDRun,
                MAIN(cr) ? DDRole::Main : DDRole::Agent,
                PAR(cr) ? NumRanks::Multiple : NumRanks::Single,
                cr->mpi_comm_mysim,
                ms,
                disresdata,
                globalState.get(),
                replExParams.exchangeInterval > 0);

    if (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0 && isSimulationMainRank)
    {
        extendStateWithOriresHistory(mtop, *inputrec, globalState.get());
    }

#if GMX_FAHCORE
    /* We have to remember the generation's first step before reading checkpoint.
       This way, we can report to the F@H core both the generation's first step
       and the restored first step, thus making it able to distinguish between
       an interruption/resume and start of the n-th generation simulation.
       Having this information, the F@H core can correctly calculate and report
       the progress.
     */
    int gen_first_step = 0;
    if (MAIN(cr))
    {
        gen_first_step = inputrec->init_step;
    }
#endif

    ObservablesHistory observablesHistory = {};

    auto modularSimulatorCheckpointData = std::make_unique<ReadCheckpointDataHolder>();
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

        // Finish applying initial simulation state information from external sources on all ranks.
        // Reconcile checkpoint file data with Mdrunner state established up to this point.
        applyLocalState(*inputHolder_.get(),
                        logFileHandle,
                        cr,
                        domdecOptions.numCells,
                        inputrec.get(),
                        globalState.get(),
                        &observablesHistory,
                        mdrunOptions.reproducible,
                        mdModules_->notifiers(),
                        modularSimulatorCheckpointData.get(),
                        useModularSimulator);
        // TODO: (#3652) Synchronize filesystem state, SimulationInput contents, and program
        //  invariants
        //  on all code paths.
        // Write checkpoint or provide hook to update SimulationInput.
        // If there was a checkpoint file, SimulationInput contains more information
        // than if there wasn't. At this point, we have synchronized the in-memory
        // state with the filesystem state only for restarted simulations. We should
        // be calling applyLocalState unconditionally and expect that the completeness
        // of SimulationInput is not dependent on its creation method.

        if (startingBehavior == StartingBehavior::RestartWithAppending && logFileHandle)
        {
            // Now we can start normal logging to the truncated log file.
            fplog = gmx_fio_getfp(logFileHandle);
            prepareLogAppending(fplog);
            logOwner = buildLogger(fplog, MAIN(cr));
            mdlog    = logOwner.logger();
        }
    }

#if GMX_FAHCORE
    if (MAIN(cr))
    {
        fcRegisterSteps(inputrec->nsteps + inputrec->init_step, gen_first_step);
    }
#endif

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
    override_nsteps_cmdline(mdlog, mdrunOptions.numStepsCommandline, inputrec.get());

    matrix box;
    if (isSimulationMainRank)
    {
        copy_mat(globalState->box, box);
    }

    if (PAR(cr))
    {
        gmx_bcast(sizeof(box), box, cr->mpiDefaultCommunicator);
    }

    if (inputrec->cutoff_scheme != CutoffScheme::Verlet)
    {
        gmx_fatal(FARGS,
                  "This group-scheme .tpr file can no longer be run by mdrun. Please update to the "
                  "Verlet scheme, or use an earlier version of GROMACS if necessary.");
    }
    /* Update rlist and nstlist. */
    /* Note: prepare_verlet_scheme is calling increaseNstlist(...), which (while attempting to
     * increase rlist) tries to check if the newly chosen value fits with the DD scheme. As this is
     * run before any DD scheme is set up, this check is never executed. See #3334 for more details.
     */
    prepare_verlet_scheme(fplog,
                          cr,
                          inputrec.get(),
                          nstlist_cmdline,
                          mtop,
                          MAIN(cr) ? globalState->x : gmx::ArrayRef<const gmx::RVec>(),
                          box,
                          useGpuForNonbonded || (emulateGpuNonbonded == EmulateGpuNonbonded::Yes),
                          *hwinfo_->cpuInfo);

    // We need to decide on update groups early, as this affects
    // inter-domain communication distances.
    auto         updateGroupingsPerMoleculeTypeResult = makeUpdateGroupingsPerMoleculeType(mtop);
    UpdateGroups updateGroups;
    if (std::holds_alternative<std::string>(updateGroupingsPerMoleculeTypeResult))
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendTextFormatted("Update groups can not be used for this system because %s",
                                     std::get<std::string>(updateGroupingsPerMoleculeTypeResult).c_str());
    }
    else
    {
        auto updateGroupingsPerMoleculeType =
                std::get<std::vector<RangePartitioning>>(updateGroupingsPerMoleculeTypeResult);
        const real maxUpdateGroupRadius = computeMaxUpdateGroupRadius(
                mtop, updateGroupingsPerMoleculeType, maxReferenceTemperature(*inputrec));
        const real cutoffMargin = std::sqrt(max_cutoff2(inputrec->pbcType, box)) - inputrec->rlist;
        updateGroups            = makeUpdateGroups(mdlog,
                                        std::move(updateGroupingsPerMoleculeType),
                                        maxUpdateGroupRadius,
                                        doRerun,
                                        useDomainDecomposition,
                                        systemHasConstraintsOrVsites(mtop),
                                        cutoffMargin);
    }

    try
    {
        const bool haveFrozenAtoms = inputrecFrozenAtoms(inputrec.get());

        useGpuForUpdate = decideWhetherToUseGpuForUpdate(useDomainDecomposition,
                                                         updateGroups.useUpdateGroups(),
                                                         pmeRunMode,
                                                         domdecOptions.numPmeRanks > 0,
                                                         useGpuForNonbonded,
                                                         updateTarget,
                                                         gpusWereDetected,
                                                         *inputrec,
                                                         mtop,
                                                         doEssentialDynamics,
                                                         gmx_mtop_ftype_count(mtop, F_ORIRES) > 0,
                                                         haveFrozenAtoms,
                                                         useModularSimulator,
                                                         doRerun,
                                                         mdlog);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    const bool canUseDirectGpuComm = decideWhetherDirectGpuCommunicationCanBeUsed(
            devFlags, inputrec->useMts, (inputrec->eSwapCoords != SwapType::No), mdlog);

    bool useGpuDirectHalo = false;

    if (useGpuForNonbonded)
    {
        // cr->npmenodes is not yet initialized.
        // domdecOptions.numPmeRanks == -1 results in 0 separate PME ranks when useGpuForNonbonded is true.
        // Todo: remove this assumption later once auto mode has support for separate PME rank
        const int numPmeRanks = domdecOptions.numPmeRanks > 0 ? domdecOptions.numPmeRanks : 0;
        bool      havePPDomainDecomposition = (cr->sizeOfDefaultCommunicator - numPmeRanks) > 1;
        useGpuDirectHalo = decideWhetherToUseGpuForHalo(havePPDomainDecomposition,
                                                        useGpuForNonbonded,
                                                        canUseDirectGpuComm,
                                                        useModularSimulator,
                                                        doRerun,
                                                        EI_ENERGY_MINIMIZATION(inputrec->eI),
                                                        mdlog);
    }

    // This builder is necessary while we have multi-part construction
    // of DD. Before DD is constructed, we use the existence of
    // the builder object to indicate that further construction of DD
    // is needed.
    std::unique_ptr<DomainDecompositionBuilder> ddBuilder;
    if (useDomainDecomposition)
    {
        // The DD builder will disable useGpuDirectHalo if the Y or Z component of any domain is
        // smaller than twice the communication distance, since GPU-direct communication presently
        // only works with a single pulse in these dimensions, and we want to avoid box scaling
        // resulting in fatal errors far into the simulation. Such small systems will not
        // perform well on multiple GPUs in any case, but it is important that our core functionality
        // (in particular for testing) does not break depending on GPU direct communication being enabled.
        ddBuilder = std::make_unique<DomainDecompositionBuilder>(
                mdlog,
                cr,
                domdecOptions,
                mdrunOptions,
                mtop,
                *inputrec,
                mdModules_->notifiers(),
                box,
                updateGroups.updateGroupingPerMoleculeType(),
                updateGroups.useUpdateGroups(),
                updateGroups.maxUpdateGroupRadius(),
                positionsFromStatePointer(globalState.get()),
                useGpuForNonbonded,
                useGpuForPme,
                useGpuForUpdate,
                useGpuDirectHalo,
                devFlags.enableGpuPmeDecomposition);
    }
    else
    {
        /* PME, if used, is done on all nodes with 1D decomposition */
        cr->mpi_comm_mygroup = cr->mpiDefaultCommunicator;
        cr->mpi_comm_mysim   = cr->mpiDefaultCommunicator;

        cr->nnodes     = cr->sizeOfDefaultCommunicator;
        cr->sim_nodeid = cr->rankInDefaultCommunicator;
        cr->nodeid     = cr->rankInDefaultCommunicator;
        cr->npmenodes  = 0;
        cr->duty       = (DUTY_PP | DUTY_PME);

        if (inputrec->pbcType == PbcType::Screw)
        {
            gmx_fatal(FARGS, "pbc=screw is only implemented with domain decomposition");
        }
    }

    // Produce the task assignment for all ranks on this node - done after DD is constructed
    GpuTaskAssignments gpuTaskAssignments = GpuTaskAssignmentsBuilder::build(
            availableDevices,
            userGpuTaskAssignment,
            *hwinfo_,
            simulationCommunicator,
            physicalNodeComm,
            nonbondedTarget,
            pmeTarget,
            bondedTarget,
            updateTarget,
            useGpuForNonbonded,
            useGpuForPme,
            thisRankHasDuty(cr, DUTY_PP),
            // TODO cr->duty & DUTY_PME should imply that a PME
            // algorithm is active, but currently does not.
            usingPme(inputrec->coulombtype) && thisRankHasDuty(cr, DUTY_PME));

    // Get the device handle for the modules on this rank, nullptr
    // when no task is assigned.
    int                deviceId   = -1;
    DeviceInformation* deviceInfo = gpuTaskAssignments.initDevice(&deviceId);

    // TODO Currently this is always built, yet DD partition code
    // checks if it is built before using it. Probably it should
    // become an MDModule that is made only when another module
    // requires it (e.g. pull, CompEl, density fitting), so that we
    // don't update the local atom sets unilaterally every step.
    LocalAtomSetManager atomSets;

    // Local state and topology are declared (and perhaps constructed)
    // now, because DD needs them for the LocalTopologyChecker, but
    // they do not contain valid data until after the first DD
    // partition.
    std::unique_ptr<t_state> localStateInstance;
    t_state*                 localState;
    gmx_localtop_t           localTopology(mtop.ffparams);

    if (ddBuilder)
    {
        localStateInstance = std::make_unique<t_state>();
        localState         = localStateInstance.get();
        // TODO Pass the GPU streams to ddBuilder to use in buffer
        // transfers (e.g. halo exchange)
        cr->setDD(ddBuilder->build(&atomSets, localTopology, *localState, &observablesReducerBuilder));
        // The builder's job is done, so destruct it
        ddBuilder.reset(nullptr);
        // Note that local state still does not exist yet.
    }
    else
    {
        // Without DD, the local state is merely an alias to the global state,
        // so we don't need to allocate anything.
        localState = globalState.get();
    }

    // Ensure that all atoms within the same update group are in the
    // same periodic image. Otherwise, a simulation that did not use
    // update groups (e.g. a single-rank simulation) cannot always be
    // correctly restarted in a way that does use update groups
    // (e.g. a multi-rank simulation).
    if (isSimulationMainRank)
    {
        const bool useUpdateGroups = cr->dd ? ddUsesUpdateGroups(*cr->dd) : false;
        if (useUpdateGroups)
        {
            putUpdateGroupAtomsInSamePeriodicImage(*cr->dd, mtop, globalState->box, globalState->x);
        }
    }

    const bool disableNonbondedCalculation = (getenv("GMX_NO_NONBONDED") != nullptr);
    if (disableNonbondedCalculation)
    {
        /* turn off non-bonded calculations */
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText(
                        "Found environment variable GMX_NO_NONBONDED.\n"
                        "Disabling nonbonded calculations.");
    }

    const NumPmeDomains numPmeDomains = getNumPmeDomains(cr->dd);
    const bool useGpuPmeDecomposition = numPmeDomains.x * numPmeDomains.y > 1 && useGpuForPme;
    GMX_RELEASE_ASSERT(!useGpuPmeDecomposition || devFlags.enableGpuPmeDecomposition,
                       "GPU PME decomposition works only in the cases where it is supported");

    MdrunScheduleWorkload runScheduleWork;

    // Also populates the simulation constant workload description.
    // Note: currently the default duty is DUTY_PP | DUTY_PME for all simulations, including those without PME,
    // so this boolean is sufficient on all ranks to determine whether separate PME ranks are used,
    // but this will no longer be the case if cr->duty is changed for !usingPme(fr->ic->eeltype).
    const bool haveSeparatePmeRank = (!thisRankHasDuty(cr, DUTY_PP) || !thisRankHasDuty(cr, DUTY_PME));
    runScheduleWork.simulationWork = createSimulationWorkload(*inputrec,
                                                              disableNonbondedCalculation,
                                                              devFlags,
                                                              havePPDomainDecomposition(cr),
                                                              haveSeparatePmeRank,
                                                              useGpuForNonbonded,
                                                              pmeRunMode,
                                                              useGpuForBonded,
                                                              useGpuForUpdate,
                                                              useGpuDirectHalo,
                                                              canUseDirectGpuComm,
                                                              useGpuPmeDecomposition);

    if (runScheduleWork.simulationWork.useGpuDirectCommunication && GMX_GPU_CUDA)
    {
        // Don't enable event counting with GPU Direct comm, see #3988.
        gmx::internal::disableCudaEventConsumptionCounting();
    }

    if (isSimulationMainRank && GMX_GPU_SYCL)
    {
        const SimulationWorkload& simWorkload    = runScheduleWork.simulationWork;
        bool                      haveAnyGpuWork = simWorkload.useGpuPme || simWorkload.useGpuBonded
                              || simWorkload.useGpuNonbonded || simWorkload.useGpuUpdate;
        if (haveAnyGpuWork)
        {
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText(
                            "\nNOTE: SYCL GPU support in GROMACS, and the compilers, libraries,\n"
                            "and drivers that it depends on are fairly new.\n"
                            "Please, pay extra attention to the correctness of your results,\n"
                            "and update to the latest GROMACS patch version if warranted.");
        }
    }

    const bool printHostName = (cr->nnodes > 1);
    gpuTaskAssignments.reportGpuUsage(mdlog, printHostName, pmeRunMode, runScheduleWork.simulationWork);

    std::unique_ptr<DeviceStreamManager> deviceStreamManager = nullptr;

    if (deviceInfo != nullptr)
    {
        if (runScheduleWork.simulationWork.havePpDomainDecomposition && thisRankHasDuty(cr, DUTY_PP))
        {
            dd_setup_dlb_resource_sharing(cr, deviceId);
        }
        const bool useGpuTiming = decideGpuTimingsUsage();
        deviceStreamManager     = std::make_unique<DeviceStreamManager>(
                *deviceInfo, runScheduleWork.simulationWork, useGpuTiming);
    }

    // If the user chose a task assignment, give them some hints
    // where appropriate.
    if (!userGpuTaskAssignment.empty())
    {
        gpuTaskAssignments.logPerformanceHints(mdlog, numAvailableDevices);
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
                        ms->simulationIndex_,
                        ms->numSimulations_);
    }
    GMX_LOG(mdlog.warning)
            .appendTextFormatted("Using %d MPI %s\n",
                                 cr->nnodes,
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
    gmx_check_thread_affinity_set(
            mdlog, &hw_opt, hwinfo_->hardwareTopology->maxThreads(), FALSE, libraryWorldCommunicator);
    /* Check and update the number of OpenMP threads requested */
    checkAndUpdateRequestedNumOpenmpThreads(
            &hw_opt, *hwinfo_, cr, ms, physicalNodeComm.size_, pmeRunMode, mtop, *inputrec);

    gmx_omp_nthreads_init(mdlog,
                          cr,
                          hwinfo_->hardwareTopology->maxThreads(),
                          physicalNodeComm.size_,
                          hw_opt.nthreads_omp,
                          hw_opt.nthreads_omp_pme,
                          !thisRankHasDuty(cr, DUTY_PP));

    const bool bEnableFPE = gmxShouldEnableFPExceptions();
    // FIXME - reconcile with gmx_feenableexcept() call from CommandLineModuleManager::run()
    if (bEnableFPE)
    {
        gmx_feenableexcept();
    }

    /* Now that we know the setup is consistent, check for efficiency */
    check_resource_division_efficiency(hwinfo_, gpuTaskAssignments.thisRankHasAnyGpuTask(), cr, mdlog);

    /* getting number of PP/PME threads on this MPI / tMPI rank.
       PME: env variable should be read only on one node to make sure it is
       identical everywhere;
     */
    const int numThreadsOnThisRank = thisRankHasDuty(cr, DUTY_PP)
                                             ? gmx_omp_nthreads_get(ModuleMultiThread::Nonbonded)
                                             : gmx_omp_nthreads_get(ModuleMultiThread::Pme);
    checkHardwareOversubscription(
            numThreadsOnThisRank, cr->nodeid, *hwinfo_->hardwareTopology, physicalNodeComm, mdlog);

    // Enable Peer access between GPUs where available
    // Only for DD, only main PP rank needs to perform setup, and only if thread MPI plus
    // any of the GPU communication features are active.
    if (haveDDAtomOrdering(*cr) && MAIN(cr) && thisRankHasDuty(cr, DUTY_PP) && GMX_THREAD_MPI
        && (runScheduleWork.simulationWork.useGpuHaloExchange
            || runScheduleWork.simulationWork.useGpuPmePpCommunication))
    {
        setupGpuDevicePeerAccess(gpuTaskAssignments.deviceIdsAssigned(), mdlog);
    }

    if (mdrunOptions.timingOptions.resetStep > -1)
    {
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -resetstep functionality is deprecated, and may be removed in a "
                        "future version.");
    }
    std::unique_ptr<gmx_wallcycle> wcycle =
            wallcycle_init(fplog, mdrunOptions.timingOptions.resetStep, cr);

    if (PAR(cr))
    {
        /* Main synchronizes its value of reset_counters with all nodes
         * including PME only nodes */
        int64_t reset_counters = wcycle_get_reset_counters(wcycle.get());
        gmx_bcast(sizeof(reset_counters), &reset_counters, cr->mpi_comm_mysim);
        wcycle_set_reset_counters(wcycle.get(), reset_counters);
    }

    // Membrane embedding must be initialized before we call init_forcerec()
    membedHolder.initializeMembed(fplog,
                                  filenames.size(),
                                  filenames.data(),
                                  &mtop,
                                  inputrec.get(),
                                  globalState.get(),
                                  cr,
                                  &mdrunOptions.checkpointOptions.period);

    const bool thisRankHasPmeGpuTask = gpuTaskAssignments.thisRankHasPmeGpuTask();
    std::unique_ptr<BoxDeformation>      deform;
    std::unique_ptr<MDAtoms>             mdAtoms;
    std::unique_ptr<VirtualSitesHandler> vsite;

    t_nrnb nrnb;
    if (thisRankHasDuty(cr, DUTY_PP))
    {
        setupNotifier.notify(*cr);
        setupNotifier.notify(&atomSets);
        setupNotifier.notify(mtop);
        setupNotifier.notify(inputrec->pbcType);
        setupNotifier.notify(SimulationTimeStep{ inputrec->delta_t });

        /* Initiate forcerecord */
        fr                 = std::make_unique<t_forcerec>();
        fr->forceProviders = mdModules_->initForceProviders();
        init_forcerec(fplog,
                      mdlog,
                      runScheduleWork.simulationWork,
                      fr.get(),
                      *inputrec,
                      mtop,
                      cr,
                      box,
                      opt2fn("-table", filenames.size(), filenames.data()),
                      opt2fn("-tablep", filenames.size(), filenames.data()),
                      opt2fns("-tableb", filenames.size(), filenames.data()),
                      pforce);
        // Dirty hack, for fixing disres and orires should be made mdmodules
        fr->fcdata->disres = disresdata;
        if (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0)
        {
            fr->fcdata->orires = std::make_unique<t_oriresdata>(
                    fplog, mtop, *inputrec, ms, globalState.get(), &atomSets);
        }

        deform = buildBoxDeformation(globalState != nullptr
                                             ? createMatrix3x3FromLegacyMatrix(globalState->box)
                                             : diagonalMatrix<real, 3, 3>(0.0),
                                     MAIN(cr) ? DDRole::Main : DDRole::Agent,
                                     PAR(cr) ? NumRanks::Multiple : NumRanks::Single,
                                     cr->mpi_comm_mygroup,
                                     *inputrec);

        // Save a handle to device stream manager to use elsewhere in the code
        // TODO: Forcerec is not a correct place to store it.
        fr->deviceStreamManager = deviceStreamManager.get();

        if (runScheduleWork.simulationWork.useGpuPmePpCommunication && !thisRankHasDuty(cr, DUTY_PME))
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
                    cr->mpi_comm_mysim,
                    cr->dd->pme_nodeid,
                    &cr->dd->pmeForceReceiveBuffer,
                    deviceStreamManager->context(),
                    deviceStreamManager->stream(DeviceStreamType::PmePpTransfer),
                    runScheduleWork.simulationWork.useNvshmem);
        }

        fr->nbv = Nbnxm::init_nb_verlet(
                mdlog,
                *inputrec,
                *fr,
                cr,
                *hwinfo_,
                runScheduleWork.simulationWork.useGpuNonbonded,
                deviceStreamManager.get(),
                mtop,
                PAR(cr) ? &observablesReducerBuilder : nullptr,
                isSimulationMainRank ? globalState->x : gmx::ArrayRef<const gmx::RVec>(),
                box,
                wcycle.get());
        // TODO: Move the logic below to a GPU bonded builder
        if (runScheduleWork.simulationWork.useGpuBonded)
        {
            GMX_RELEASE_ASSERT(deviceStreamManager != nullptr,
                               "GPU device stream manager should be valid in order to use GPU "
                               "version of bonded forces.");
            fr->listedForcesGpu =
                    std::make_unique<ListedForcesGpu>(mtop.ffparams,
                                                      fr->ic->epsfac * fr->fudgeQQ,
                                                      inputrec->opts.ngener - inputrec->nwall,
                                                      *deviceInfo,
                                                      deviceStreamManager->context(),
                                                      deviceStreamManager->bondedStream(),
                                                      wcycle.get());
        }
        fr->longRangeNonbondeds = std::make_unique<CpuPpLongRangeNonbondeds>(fr->n_tpi,
                                                                             fr->ic->ewaldcoeff_q,
                                                                             fr->ic->epsilon_r,
                                                                             fr->qsum,
                                                                             fr->ic->eeltype,
                                                                             fr->ic->vdwtype,
                                                                             *inputrec,
                                                                             &nrnb,
                                                                             wcycle.get(),
                                                                             fplog);

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
        vsite = makeVirtualSitesHandler(
                mtop, cr, fr->pbcType, updateGroups.updateGroupingPerMoleculeType());

        calc_shifts(box, fr->shift_vec);

        /* With periodic molecules the charge groups should be whole at start up
         * and the virtual sites should not be far from their proper positions.
         */
        if (!inputrec->bContinuation && MAIN(cr)
            && !(inputrec->pbcType != PbcType::No && inputrec->bPeriodicMols))
        {
            /* Make molecules whole at start of run */
            if (fr->pbcType != PbcType::No)
            {
                do_pbc_first_mtop(fplog,
                                  inputrec->pbcType,
                                  ir_haveBoxDeformation(*inputrec),
                                  inputrec->deform,
                                  box,
                                  &mtop,
                                  globalState->x,
                                  globalState->v);
            }
            if (vsite)
            {
                /* Correct initial vsite positions are required
                 * for the initial distribution in the domain decomposition
                 * and for the initial shell prediction.
                 */
                constructVirtualSitesGlobal(mtop, globalState->x);
            }
        }
        // Make the DD reverse topology, now that any vsites that are present are available
        if (haveDDAtomOrdering(*cr))
        {
            dd_make_reverse_top(fplog, cr->dd, mtop, vsite.get(), *inputrec, domdecOptions.ddBondedChecking);
        }

        if (usingPme(fr->ic->eeltype) || usingLJPme(fr->ic->vdwtype))
        {
            ewaldcoeff_q  = fr->ic->ewaldcoeff_q;
            ewaldcoeff_lj = fr->ic->ewaldcoeff_lj;
        }
    }
    else
    {
        /* This is a PME only node */

        GMX_ASSERT(globalState == nullptr,
                   "We don't need the state on a PME only rank and expect it to be uninitialized");

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
    if (usingPme(inputrec->coulombtype) || usingLJPme(inputrec->vdwtype))
    {
        if (mdAtoms && mdAtoms->mdatoms())
        {
            nChargePerturbed = mdAtoms->mdatoms()->nChargePerturbed;
            if (usingLJPme(inputrec->vdwtype))
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
                GMX_RELEASE_ASSERT(!runScheduleWork.simulationWork.useGpuPme
                                           || (deviceStreamManager != nullptr),
                                   "Device stream manager should be valid in order to use GPU "
                                   "version of PME.");
                GMX_RELEASE_ASSERT(
                        !runScheduleWork.simulationWork.useGpuPme
                                || deviceStreamManager->streamIsValid(DeviceStreamType::Pme),
                        "GPU PME stream should be valid in order to use GPU version of PME.");

                const DeviceContext* deviceContext = runScheduleWork.simulationWork.useGpuPme
                                                             ? &deviceStreamManager->context()
                                                             : nullptr;
                const DeviceStream*  pmeStream =
                        runScheduleWork.simulationWork.useGpuPme
                                 ? &deviceStreamManager->stream(DeviceStreamType::Pme)
                                 : nullptr;

                const t_inputrec* ir = inputrec.get();
                /* For each atom we allow a relative (cut-off) error of up to ewald_rtol.
                 * Thus we can also tolerate an error of an order of magnitude less due to
                 * atoms being slightly outside the halo extent. This will only cause a fraction
                 * of the charge to be missing on the grid. So we pass ewald_rtol as the allowed
                 * chance per atom (ChanceTarget::Atom) to be outside the halo extent.
                 */
                const real haloExtentForAtomDisplacement =
                        updateGroups.maxUpdateGroupRadius()
                        + minCellSizeForAtomDisplacement(mtop,
                                                         *ir,
                                                         updateGroups.updateGroupingPerMoleculeType(),
                                                         ir->ewald_rtol,
                                                         ChanceTarget::Atom);
                pmedata = gmx_pme_init(cr,
                                       getNumPmeDomains(cr->dd),
                                       ir,
                                       box,
                                       haloExtentForAtomDisplacement,
                                       nChargePerturbed != 0,
                                       nTypePerturbed != 0,
                                       mdrunOptions.reproducible,
                                       ewaldcoeff_q,
                                       ewaldcoeff_lj,
                                       gmx_omp_nthreads_get(ModuleMultiThread::Pme),
                                       pmeRunMode,
                                       nullptr,
                                       deviceContext,
                                       pmeStream,
                                       pmeGpuProgram.get(),
                                       mdlog);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
    }

    std::unique_ptr<gmxNvshmemHandle> nvshmemHandlePtr;
    if (runScheduleWork.simulationWork.useNvshmem)
    {
        nvshmemHandlePtr = std::make_unique<gmxNvshmemHandle>(cr->mpiDefaultCommunicator);
    }

    /* Set thread affinity after gmx_pme_init(), otherwise with cuFFTMp the NVSHMEM helper thread
     * can be pinned to the same core as the PME thread, causing performance degradation.
     */
    if (hw_opt.threadAffinity != ThreadAffinity::Off)
    {
        /* Before setting affinity, check whether the affinity has changed
         * - which indicates that probably the OpenMP library has changed it
         * since we first checked).
         */
        gmx_check_thread_affinity_set(
                mdlog, &hw_opt, hwinfo_->hardwareTopology->maxThreads(), TRUE, libraryWorldCommunicator);

        int numThreadsOnThisNode, intraNodeThreadOffset;
        analyzeThreadsOnThisNode(
                physicalNodeComm, numThreadsOnThisRank, &numThreadsOnThisNode, &intraNodeThreadOffset);

        /* Set the CPU affinity */
        gmx_set_thread_affinity(mdlog,
                                cr,
                                &hw_opt,
                                *hwinfo_->hardwareTopology,
                                numThreadsOnThisRank,
                                numThreadsOnThisNode,
                                intraNodeThreadOffset,
                                nullptr);
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

    try
    {
        pull_t* pull_work = nullptr;
        if (thisRankHasDuty(cr, DUTY_PP))
        {
            /* Assumes uniform use of the number of OpenMP threads */
            walltime_accounting =
                    walltime_accounting_init(gmx_omp_nthreads_get(ModuleMultiThread::Default));

            if (inputrec->bPull)
            {
                /* Initialize pull code */
                pull_work = init_pull(fplog,
                                      inputrec->pull.get(),
                                      inputrec.get(),
                                      mtop,
                                      cr,
                                      &atomSets,
                                      inputrec->fepvals->init_lambda);
                if (inputrec->pull->bXOutAverage || inputrec->pull->bFOutAverage)
                {
                    initPullHistory(pull_work, &observablesHistory);
                }
                if (EI_DYNAMICS(inputrec->eI) && MAIN(cr))
                {
                    init_pull_output_files(
                            pull_work, filenames.size(), filenames.data(), oenv, startingBehavior);
                }
            }

            std::unique_ptr<EnforcedRotation> enforcedRotation;
            if (inputrec->bRot)
            {
                /* Initialize enforced rotation code */
                enforcedRotation = init_rot(fplog,
                                            inputrec.get(),
                                            filenames.size(),
                                            filenames.data(),
                                            cr,
                                            &atomSets,
                                            globalState.get(),
                                            mtop,
                                            oenv,
                                            mdrunOptions,
                                            startingBehavior);
            }

            t_swap* swap = nullptr;
            if (inputrec->eSwapCoords != SwapType::No)
            {
                /* Initialize ion swapping code */
                swap = init_swapcoords(fplog,
                                       inputrec.get(),
                                       opt2fn_main("-swap", filenames.size(), filenames.data(), cr),
                                       mtop,
                                       globalState.get(),
                                       &observablesHistory,
                                       cr,
                                       &atomSets,
                                       oenv,
                                       mdrunOptions,
                                       startingBehavior);
            }

            /* Let makeConstraints know whether we have essential dynamics constraints. */
            auto constr = makeConstraints(mtop,
                                          *inputrec,
                                          pull_work,
                                          pull_work != nullptr ? pull_have_constraint(*pull_work) : false,
                                          doEssentialDynamics,
                                          fplog,
                                          cr,
                                          updateGroups.useUpdateGroups(),
                                          ms,
                                          &nrnb,
                                          wcycle.get(),
                                          fr->bMolPBC,
                                          PAR(cr) ? &observablesReducerBuilder : nullptr);

            /* Energy terms and groups */
            gmx_enerdata_t enerd(mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                                 &inputrec->fepvals->all_lambda);

            // cos acceleration is only supported by md, but older tpr
            // files might still combine it with other integrators
            GMX_RELEASE_ASSERT(inputrec->cos_accel == 0.0 || inputrec->eI == IntegrationAlgorithm::MD,
                               "cos_acceleration is only supported by integrator=md");

            /* Kinetic energy data */
            gmx_ekindata_t ekind(gmx::constArrayRefFromArray(inputrec->opts.ref_t, inputrec->opts.ngtc),
                                 inputrec->ensembleTemperatureSetting,
                                 inputrec->ensembleTemperature,
                                 fr->haveBoxDeformation,
                                 inputrec->cos_accel,
                                 gmx_omp_nthreads_get(ModuleMultiThread::Update));

            /* Set up interactive MD (IMD) */
            auto imdSession = makeImdSession(inputrec.get(),
                                             cr,
                                             wcycle.get(),
                                             &enerd,
                                             ms,
                                             mtop,
                                             mdlog,
                                             MAIN(cr) ? globalState->x : gmx::ArrayRef<gmx::RVec>(),
                                             filenames.size(),
                                             filenames.data(),
                                             oenv,
                                             mdrunOptions.imdOptions,
                                             startingBehavior);

            if (haveDDAtomOrdering(*cr))
            {
                GMX_RELEASE_ASSERT(fr, "fr was NULL while cr->duty was DUTY_PP");
                /* This call is not included in init_domain_decomposition
                 * because fr->atomInfoForEachMoleculeBlock is set later.
                 */
                makeBondedLinks(cr->dd, mtop, fr->atomInfoForEachMoleculeBlock);
            }

            if (runScheduleWork.simulationWork.useGpuFBufferOpsWhenAllowed)
            {
                fr->gpuForceReduction[gmx::AtomLocality::Local] = std::make_unique<gmx::GpuForceReduction>(
                        deviceStreamManager->context(),
                        deviceStreamManager->stream(gmx::DeviceStreamType::NonBondedLocal),
                        wcycle.get());

                if (runScheduleWork.simulationWork.havePpDomainDecomposition)
                {
                    fr->gpuForceReduction[gmx::AtomLocality::NonLocal] =
                            std::make_unique<gmx::GpuForceReduction>(
                                    deviceStreamManager->context(),
                                    deviceStreamManager->stream(gmx::DeviceStreamType::NonBondedNonLocal),
                                    wcycle.get());
                }

                if (runScheduleWork.simulationWork.useMdGpuGraph)
                {
                    fr->mdGraph[MdGraphEvenOrOddStep::EvenStep] =
                            std::make_unique<gmx::MdGpuGraph>(*fr->deviceStreamManager,
                                                              runScheduleWork.simulationWork,
                                                              cr->mpi_comm_mygroup,
                                                              MdGraphEvenOrOddStep::EvenStep,
                                                              wcycle.get());

                    fr->mdGraph[MdGraphEvenOrOddStep::OddStep] =
                            std::make_unique<gmx::MdGpuGraph>(*fr->deviceStreamManager,
                                                              runScheduleWork.simulationWork,
                                                              cr->mpi_comm_mygroup,
                                                              MdGraphEvenOrOddStep::OddStep,
                                                              wcycle.get());

                    fr->mdGraph[MdGraphEvenOrOddStep::EvenStep]->setAlternateStepPpTaskCompletionEvent(
                            fr->mdGraph[MdGraphEvenOrOddStep::OddStep]->getPpTaskCompletionEvent());
                    fr->mdGraph[MdGraphEvenOrOddStep::OddStep]->setAlternateStepPpTaskCompletionEvent(
                            fr->mdGraph[MdGraphEvenOrOddStep::EvenStep]->getPpTaskCompletionEvent());
                }
            }

            std::unique_ptr<gmx::StatePropagatorDataGpu> stateGpu;
            if (gpusWereDetected && gmx::needStateGpu(runScheduleWork.simulationWork))
            {
                GpuApiCallBehavior transferKind =
                        (inputrec->eI == IntegrationAlgorithm::MD && !doRerun && !useModularSimulator)
                                ? GpuApiCallBehavior::Async
                                : GpuApiCallBehavior::Sync;
                GMX_RELEASE_ASSERT(deviceStreamManager != nullptr,
                                   "GPU device stream manager should be initialized to use GPU.");
                stateGpu = std::make_unique<gmx::StatePropagatorDataGpu>(
                        *deviceStreamManager, transferKind, pme_gpu_get_block_size(fr->pmedata), wcycle.get());
                fr->stateGpu = stateGpu.get();
            }

            GMX_ASSERT(stopHandlerBuilder_, "Runner must provide StopHandlerBuilder to simulator.");
            SimulatorBuilder simulatorBuilder;

            simulatorBuilder.add(SimulatorStateData(
                    globalState.get(), localState, &observablesHistory, &enerd, &ekind));
            simulatorBuilder.add(std::move(membedHolder));
            simulatorBuilder.add(std::move(stopHandlerBuilder_));
            simulatorBuilder.add(SimulatorConfig(mdrunOptions, startingBehavior, &runScheduleWork));


            simulatorBuilder.add(SimulatorEnv(fplog, cr, ms, mdlog, oenv, &observablesReducerBuilder));
            simulatorBuilder.add(Profiling(&nrnb, walltime_accounting, wcycle.get()));
            simulatorBuilder.add(ConstraintsParam(
                    constr.get(),
                    enforcedRotation ? enforcedRotation->getLegacyEnfrot() : nullptr,
                    vsite.get()));
            // TODO: Separate `fr` to a separate add, and make the `build` handle the coupling sensibly.
            simulatorBuilder.add(LegacyInput(
                    static_cast<int>(filenames.size()), filenames.data(), inputrec.get(), fr.get()));
            simulatorBuilder.add(ReplicaExchangeParameters(replExParams));
            simulatorBuilder.add(InteractiveMD(imdSession.get()));
            simulatorBuilder.add(SimulatorModules(mdModules_->outputProvider(), mdModules_->notifiers()));
            simulatorBuilder.add(CenterOfMassPulling(pull_work));
            // Todo move to an MDModule
            simulatorBuilder.add(IonSwapping(swap));
            simulatorBuilder.add(TopologyData(mtop, &localTopology, mdAtoms.get()));
            simulatorBuilder.add(BoxDeformationHandle(deform.get()));
            simulatorBuilder.add(std::move(modularSimulatorCheckpointData));

            // build and run simulator object based on user-input
            auto simulator = simulatorBuilder.build(useModularSimulator);
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
            walltime_accounting = walltime_accounting_init(gmx_omp_nthreads_get(ModuleMultiThread::Pme));
            gmx_pmeonly(&pmedata,
                        cr,
                        &nrnb,
                        wcycle.get(),
                        walltime_accounting,
                        inputrec.get(),
                        pmeRunMode,
                        runScheduleWork.simulationWork.useGpuPmePpCommunication,
                        runScheduleWork.simulationWork.useNvshmem,
                        deviceStreamManager.get());
        }

        if (!hwinfo_->deviceInfoList.empty())
        {
            /* stop the GPU profiler (only CUDA);
             * Doing so here avoids including a lot of cleanup/freeing API calls in the trace. */
            stopGpuProfiler();
        }

        wallcycle_stop(wcycle.get(), WallCycleCounter::Run);

        /* Finish up, write some stuff
         * if rerunMD, don't write last frame again
         */
        finish_run(fplog,
                   mdlog,
                   cr,
                   *inputrec,
                   &nrnb,
                   wcycle.get(),
                   walltime_accounting,
                   fr ? fr->nbv.get() : nullptr,
                   pmedata,
                   EI_DYNAMICS(inputrec->eI) && !isMultiSim(ms));
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    try
    {
        // Free PME data
        if (pmedata)
        {
            gmx_pme_destroy(pmedata);
            pmedata = nullptr;
        }

        // FIXME: this is only here to manually unpin mdAtoms->chargeA_ and state->x,
        // before we destroy the GPU context(s)
        // Pinned buffers are associated with contexts in CUDA.
        // As soon as we destroy GPU contexts after mdrunner() exits, these lines should go.
        cr->destroyDD();
        mdAtoms.reset(nullptr);
        globalState.reset(nullptr);
        localStateInstance.reset(nullptr);
        mdModules_.reset(nullptr); // destruct force providers here as they might also use the GPU
        fr.reset(nullptr);         // destruct forcerec before gpu
        // TODO convert to C++ so we can get rid of these frees
        sfree(disresdata);

        // Destroy streams after all the structures using them
        deviceStreamManager.reset(nullptr);

        /* With tMPI we need to wait for all ranks to finish deallocation before
         * destroying the CUDA context as some tMPI ranks may be sharing
         * GPU and context.
         *
         * This is not a concern in OpenCL where we use one context per rank.
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

        const bool haveDetectedOrForcedCudaAwareMpi =
                (gmx::checkMpiCudaAwareSupport() == gmx::GpuAwareMpiStatus::Supported
                 || gmx::checkMpiCudaAwareSupport() == gmx::GpuAwareMpiStatus::Forced);
        if (!haveDetectedOrForcedCudaAwareMpi)
        {
            // Don't reset GPU in case of GPU-AWARE MPI
            // UCX creates GPU buffers which are cleaned-up as part of MPI_Finalize()
            // resetting the device before MPI_Finalize() results in crashes inside UCX
            // This can also cause issues in tests that invoke mdrunner() multiple
            // times in the same process; ref #3952.
            releaseDevice(deviceInfo);
        }
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

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
       exit this function, but the main thread needs to be told to
       wait for that. */
    if (MAIN(cr))
    {
        tMPI_Finalize();
    }
#endif

    return rc;
} // Mdrunner::mdrunner

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
}

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

Mdrunner& Mdrunner::operator=(Mdrunner&& /*handle*/) noexcept = default;

class Mdrunner::BuilderImplementation
{
public:
    BuilderImplementation() = delete;
    BuilderImplementation(std::unique_ptr<MDModules> mdModules, compat::not_null<SimulationContext*> context);
    ~BuilderImplementation();

    BuilderImplementation& setExtraMdrunOptions(const MdrunOptions& options,
                                                real                forceWarningThreshold,
                                                StartingBehavior    startingBehavior);

    void addHardwareDetectionResult(const gmx_hw_info_t* hwinfo);

    void addDomdec(const DomdecOptions& options);

    void addInput(SimulationInputHandle inputHolder);

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

    //! World communicator, used for hardware detection and task assignment
    MPI_Comm libraryWorldCommunicator_ = MPI_COMM_NULL;

    //! Multisim communicator handle.
    gmx_multisim_t* multiSimulation_;

    //! mdrun communicator
    MPI_Comm simulationCommunicator_ = MPI_COMM_NULL;

    //! Print a warning if any force is larger than this (in kJ/mol nm).
    real forceWarningThreshold_ = -1;

    //! Whether the simulation will start afresh, or restart with/without appending.
    StartingBehavior startingBehavior_ = StartingBehavior::NewSimulation;

    //! The modules that comprise the functionality of mdrun.
    std::unique_ptr<MDModules> mdModules_;

    //! Detected hardware.
    const gmx_hw_info_t* hwinfo_ = nullptr;

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

    /*!
     * \brief Sources for initial simulation state.
     *
     * See issue #3652 for near-term refinements to the SimulationInput interface.
     *
     * See issue #3379 for broader discussion on API aspects of simulation inputs and outputs.
     */
    SimulationInputHandle inputHolder_;
};

Mdrunner::BuilderImplementation::BuilderImplementation(std::unique_ptr<MDModules> mdModules,
                                                       compat::not_null<SimulationContext*> context) :
    mdModules_(std::move(mdModules))
{
    libraryWorldCommunicator_ = context->libraryWorldCommunicator_;
    simulationCommunicator_   = context->simulationCommunicator_;
    multiSimulation_          = context->multiSimulation_.get();
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

    newRunner.libraryWorldCommunicator = libraryWorldCommunicator_;

    newRunner.simulationCommunicator = simulationCommunicator_;

    // nullptr is a valid value for the multisim handle
    newRunner.ms = multiSimulation_;

    if (hwinfo_)
    {
        newRunner.hwinfo_ = hwinfo_;
    }
    else
    {
        GMX_THROW(gmx::APIError(
                "MdrunnerBuilder::addHardwareDetectionResult() is required before build()"));
    }

    if (inputHolder_)
    {
        newRunner.inputHolder_ = std::move(inputHolder_);
    }
    else
    {
        GMX_THROW(gmx::APIError("MdrunnerBuilder::addInput() is required before build()."));
    }

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

void Mdrunner::BuilderImplementation::addHardwareDetectionResult(const gmx_hw_info_t* hwinfo)
{
    hwinfo_ = hwinfo;
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

void Mdrunner::BuilderImplementation::addInput(SimulationInputHandle inputHolder)
{
    inputHolder_ = std::move(inputHolder);
}

MdrunnerBuilder::MdrunnerBuilder(std::unique_ptr<MDModules>           mdModules,
                                 compat::not_null<SimulationContext*> context) :
    impl_{ std::make_unique<Mdrunner::BuilderImplementation>(std::move(mdModules), context) }
{
}

MdrunnerBuilder::~MdrunnerBuilder() = default;

MdrunnerBuilder& MdrunnerBuilder::addHardwareDetectionResult(const gmx_hw_info_t* hwinfo)
{
    impl_->addHardwareDetectionResult(hwinfo);
    return *this;
}

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

MdrunnerBuilder& MdrunnerBuilder::addInput(SimulationInputHandle input)
{
    impl_->addInput(std::move(input));
    return *this;
}

MdrunnerBuilder::MdrunnerBuilder(MdrunnerBuilder&&) noexcept = default;

MdrunnerBuilder& MdrunnerBuilder::operator=(MdrunnerBuilder&&) noexcept = default;

} // namespace gmx
