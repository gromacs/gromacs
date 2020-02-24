/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Common functions for the different NBNXN GPU implementations.
 *
 * \author Szilard Pall <pall.szilard@gmail.com>
 *
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_GPU_COMMON_H
#define GMX_NBNXM_GPU_COMMON_H

#include "config.h"

#include <string>

#if GMX_GPU == GMX_GPU_CUDA
#    include "cuda/nbnxm_cuda_types.h"
#endif

#if GMX_GPU == GMX_GPU_OPENCL
#    include "opencl/nbnxm_ocl_types.h"
#endif

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/listed_forces/gpubonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "gpu_common_utils.h"
#include "nbnxm_gpu.h"

namespace gmx
{
class GpuBonded;
}

namespace Nbnxm
{

/*! \brief Check that atom locality values are valid for the GPU module.
 *
 *  In the GPU module atom locality "all" is not supported, the local and
 *  non-local ranges are treated separately.
 *
 *  \param[in] atomLocality atom locality specifier
 */
static inline void validateGpuAtomLocality(const AtomLocality atomLocality)
{
    std::string str = gmx::formatString(
            "Invalid atom locality passed (%d); valid here is only "
            "local (%d) or nonlocal (%d)",
            static_cast<int>(atomLocality), static_cast<int>(AtomLocality::Local),
            static_cast<int>(AtomLocality::NonLocal));

    GMX_ASSERT(atomLocality == AtomLocality::Local || atomLocality == AtomLocality::NonLocal, str.c_str());
}

/*! \brief Convert atom locality to interaction locality.
 *
 *  In the current implementation the this is straightforward conversion:
 *  local to local, non-local to non-local.
 *
 *  \param[in] atomLocality Atom locality specifier
 *  \returns                Interaction locality corresponding to the atom locality passed.
 */
static inline InteractionLocality gpuAtomToInteractionLocality(const AtomLocality atomLocality)
{
    validateGpuAtomLocality(atomLocality);

    /* determine interaction locality from atom locality */
    if (atomLocality == AtomLocality::Local)
    {
        return InteractionLocality::Local;
    }
    else if (atomLocality == AtomLocality::NonLocal)
    {
        return InteractionLocality::NonLocal;
    }
    else
    {
        gmx_incons("Wrong locality");
    }
}


//NOLINTNEXTLINE(misc-definitions-in-headers)
void setupGpuShortRangeWork(NbnxmGpu* nb, const gmx::GpuBonded* gpuBonded, const gmx::InteractionLocality iLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    // There is short-range work if the pair list for the provided
    // interaction locality contains entries or if there is any
    // bonded work (as this is not split into local/nonlocal).
    nb->haveWork[iLocality] = ((nb->plist[iLocality]->nsci != 0)
                               || (gpuBonded != nullptr && gpuBonded->haveInteractions()));
}

/*! \brief Returns true if there is GPU short-range work for the given interaction locality.
 *
 * Note that as, unlike nonbonded tasks, bonded tasks are not split into local/nonlocal,
 * and therefore if there are GPU offloaded bonded interactions, this function will return
 * true for all interaction localities.
 *
 * \param[inout]  nb        Pointer to the nonbonded GPU data structure
 * \param[in]     iLocality Interaction locality identifier
 */
static bool haveGpuShortRangeWork(const NbnxmGpu& nb, const gmx::InteractionLocality iLocality)
{
    return nb.haveWork[iLocality];
}

//NOLINTNEXTLINE(misc-definitions-in-headers)
bool haveGpuShortRangeWork(const NbnxmGpu* nb, const gmx::AtomLocality aLocality)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    return haveGpuShortRangeWork(*nb, gpuAtomToInteractionLocality(aLocality));
}


/*! \brief Calculate atom range and return start index and length.
 *
 * \param[in] atomData Atom descriptor data structure
 * \param[in] atomLocality Atom locality specifier
 * \param[out] atomRangeBegin Starting index of the atom range in the atom data array.
 * \param[out] atomRangeLen Atom range length in the atom data array.
 */
template<typename AtomDataT>
static inline void getGpuAtomRange(const AtomDataT*   atomData,
                                   const AtomLocality atomLocality,
                                   int*               atomRangeBegin,
                                   int*               atomRangeLen)
{
    assert(atomData);
    validateGpuAtomLocality(atomLocality);

    /* calculate the atom data index range based on locality */
    if (atomLocality == AtomLocality::Local)
    {
        *atomRangeBegin = 0;
        *atomRangeLen   = atomData->natoms_local;
    }
    else
    {
        *atomRangeBegin = atomData->natoms_local;
        *atomRangeLen   = atomData->natoms - atomData->natoms_local;
    }
}


/*! \brief Count pruning kernel time if either kernel has been triggered
 *
 *  We do the accounting for either of the two pruning kernel flavors:
 *   - 1st pass prune: ran during the current step (prior to the force kernel);
 *   - rolling prune:  ran at the end of the previous step (prior to the current step H2D xq);
 *
 * Note that the resetting of cu_timers_t::didPrune and cu_timers_t::didRollingPrune should happen
 * after calling this function.
 *
 * \param[in] timers   structs with GPU timer objects
 * \param[inout] timings  GPU task timing data
 * \param[in] iloc        interaction locality
 */
template<typename GpuTimers>
static void countPruneKernelTime(GpuTimers*                 timers,
                                 gmx_wallclock_gpu_nbnxn_t* timings,
                                 const InteractionLocality  iloc)
{
    gpu_timers_t::Interaction& iTimers = timers->interaction[iloc];

    // We might have not done any pruning (e.g. if we skipped with empty domains).
    if (!iTimers.didPrune && !iTimers.didRollingPrune)
    {
        return;
    }

    if (iTimers.didPrune)
    {
        timings->pruneTime.c++;
        timings->pruneTime.t += iTimers.prune_k.getLastRangeTime();
    }

    if (iTimers.didRollingPrune)
    {
        timings->dynamicPruneTime.c++;
        timings->dynamicPruneTime.t += iTimers.rollingPrune_k.getLastRangeTime();
    }
}

/*! \brief Reduce data staged internally in the nbnxn module.
 *
 * Shift forces and electrostatic/LJ energies copied from the GPU into
 * a module-internal staging area are immediately reduced (CPU-side buffers passed)
 * after having waited for the transfers' completion.
 *
 * Note that this function should always be called after the transfers into the
 * staging buffers has completed.
 *
 * \tparam     StagingData    Type of staging data
 * \param[in]  nbst           Nonbonded staging data
 * \param[in]  iLocality      Interaction locality specifier
 * \param[in]  reduceEnergies True if energy reduction should be done
 * \param[in]  reduceFshift   True if shift force reduction should be done
 * \param[out] e_lj           Variable to accumulate LJ energy into
 * \param[out] e_el           Variable to accumulate electrostatic energy into
 * \param[out] fshift         Pointer to the array of shift forces to accumulate into
 */
template<typename StagingData>
static inline void gpu_reduce_staged_outputs(const StagingData&        nbst,
                                             const InteractionLocality iLocality,
                                             const bool                reduceEnergies,
                                             const bool                reduceFshift,
                                             real*                     e_lj,
                                             real*                     e_el,
                                             rvec*                     fshift)
{
    /* add up energies and shift forces (only once at local F wait) */
    if (iLocality == InteractionLocality::Local)
    {
        if (reduceEnergies)
        {
            *e_lj += *nbst.e_lj;
            *e_el += *nbst.e_el;
        }

        if (reduceFshift)
        {
            for (int i = 0; i < SHIFTS; i++)
            {
                rvec_inc(fshift[i], nbst.fshift[i]);
            }
        }
    }
}

/*! \brief Do the per-step timing accounting of the nonbonded tasks.
 *
 *  Does timing accumulation and call-count increments for the nonbonded kernels.
 *  Note that this function should be called after the current step's nonbonded
 *  nonbonded tasks have completed with the exception of the rolling pruning kernels
 *  that are accounted for during the following step.
 *
 * NOTE: if timing with multiple GPUs (streams) becomes possible, the
 *      counters could end up being inconsistent due to not being incremented
 *      on some of the node when this is skipped on empty local domains!
 *
 * \tparam     GpuTimers         GPU timers type
 * \tparam     GpuPairlist       Pair list type
 * \param[out] timings           Pointer to the NB GPU timings data
 * \param[in]  timers            Pointer to GPU timers data
 * \param[in]  plist             Pointer to the pair list data
 * \param[in]  atomLocality      Atom locality specifier
 * \param[in]  stepWork          Force schedule flags
 * \param[in]  doTiming          True if timing is enabled.
 *
 */
template<typename GpuTimers, typename GpuPairlist>
static inline void gpu_accumulate_timings(gmx_wallclock_gpu_nbnxn_t* timings,
                                          GpuTimers*                 timers,
                                          const GpuPairlist*         plist,
                                          AtomLocality               atomLocality,
                                          const gmx::StepWorkload&   stepWork,
                                          bool                       doTiming)
{
    /* timing data accumulation */
    if (!doTiming)
    {
        return;
    }

    /* determine interaction locality from atom locality */
    const InteractionLocality iLocality        = gpuAtomToInteractionLocality(atomLocality);
    const bool                didEnergyKernels = stepWork.computeEnergy;

    /* only increase counter once (at local F wait) */
    if (iLocality == InteractionLocality::Local)
    {
        timings->nb_c++;
        timings->ktime[plist->haveFreshList ? 1 : 0][didEnergyKernels ? 1 : 0].c += 1;
    }

    /* kernel timings */
    timings->ktime[plist->haveFreshList ? 1 : 0][didEnergyKernels ? 1 : 0].t +=
            timers->interaction[iLocality].nb_k.getLastRangeTime();

    /* X/q H2D and F D2H timings */
    timings->nb_h2d_t += timers->xf[atomLocality].nb_h2d.getLastRangeTime();
    timings->nb_d2h_t += timers->xf[atomLocality].nb_d2h.getLastRangeTime();

    /* Count the pruning kernel times for both cases:1st pass (at search step)
       and rolling pruning (if called at the previous step).
       We do the accounting here as this is the only sync point where we
       know (without checking or additional sync-ing) that prune tasks in
       in the current stream have completed (having just blocking-waited
       for the force D2H). */
    countPruneKernelTime(timers, timings, iLocality);

    /* only count atdat at pair-search steps (add only once, at local F wait) */
    if (stepWork.doNeighborSearch && atomLocality == AtomLocality::Local)
    {
        /* atdat transfer timing */
        timings->pl_h2d_c++;
        timings->pl_h2d_t += timers->atdat.getLastRangeTime();
    }

    /* only count pair-list H2D when actually performed */
    if (timers->interaction[iLocality].didPairlistH2D)
    {
        timings->pl_h2d_t += timers->interaction[iLocality].pl_h2d.getLastRangeTime();

        /* Clear the timing flag for the next step */
        timers->interaction[iLocality].didPairlistH2D = false;
    }
}

/*! \brief Attempts to complete nonbonded GPU task.
 *
 * See documentation in nbnxm_gpu.h for details.
 *
 * \todo Move into shared source file with gmx_compile_cpp_as_cuda
 */
//NOLINTNEXTLINE(misc-definitions-in-headers)
bool gpu_try_finish_task(NbnxmGpu*                nb,
                         const gmx::StepWorkload& stepWork,
                         const AtomLocality       aloc,
                         real*                    e_lj,
                         real*                    e_el,
                         gmx::ArrayRef<gmx::RVec> shiftForces,
                         GpuTaskCompletion        completionKind,
                         gmx_wallcycle*           wcycle)
{
    GMX_ASSERT(nb, "Need a valid nbnxn_gpu object");

    /* determine interaction locality from atom locality */
    const InteractionLocality iLocality = gpuAtomToInteractionLocality(aloc);


    // Transfers are launched and therefore need to be waited on if:
    // - buffer ops is not offloaded
    // - energies or virials are needed (on the local stream)
    //
    // (Note that useGpuFBufferOps and computeVirial are mutually exclusive
    // in current code as virial steps do CPU reduction.)
    const bool haveResultToWaitFor =
            (!stepWork.useGpuFBufferOps
             || (aloc == AtomLocality::Local && (stepWork.computeEnergy || stepWork.computeVirial)));

    //  We skip when during the non-local phase there was actually no work to do.
    //  This is consistent with nbnxn_gpu_launch_kernel but it also considers possible
    //  bonded GPU work.
    if ((iLocality == InteractionLocality::Local) || haveGpuShortRangeWork(*nb, iLocality))
    {
        // Query the state of the GPU stream and return early if we're not done
        if (completionKind == GpuTaskCompletion::Check)
        {
            // To get the wcycle call count right, when in GpuTaskCompletion::Check mode,
            // we start without counting and only when the task finished we issue a
            // start/stop to increment.
            // GpuTaskCompletion::Wait mode the timing is expected to be done in the caller.
            wallcycle_start_nocount(wcycle, ewcWAIT_GPU_NB_L);

            if (!haveStreamTasksCompleted(*nb->deviceStreams[iLocality]))
            {
                wallcycle_stop(wcycle, ewcWAIT_GPU_NB_L);

                // Early return to skip the steps below that we have to do only
                // after the NB task completed
                return false;
            }

            wallcycle_increment_event_count(wcycle, ewcWAIT_GPU_NB_L);
        }
        else if (haveResultToWaitFor)
        {
            nb->deviceStreams[iLocality]->synchronize();
        }

        // TODO: this needs to be moved later because conditional wait could brake timing
        // with a future OpenCL implementation, but with CUDA timing is anyway disabled
        // in all cases where we skip the wait.
        gpu_accumulate_timings(nb->timings, nb->timers, nb->plist[iLocality], aloc, stepWork, nb->bDoTime);

        if (stepWork.computeEnergy || stepWork.computeVirial)
        {
            gpu_reduce_staged_outputs(nb->nbst, iLocality, stepWork.computeEnergy, stepWork.computeVirial,
                                      e_lj, e_el, as_rvec_array(shiftForces.data()));
        }
    }

    /* Always reset both pruning flags (doesn't hurt doing it even when timing is off). */
    nb->timers->interaction[iLocality].didPrune =
            nb->timers->interaction[iLocality].didRollingPrune = false;

    /* Turn off initial list pruning (doesn't hurt if this is not pair-search step). */
    nb->plist[iLocality]->haveFreshList = false;

    return true;
}

/*! \brief
 * Wait for the asynchronously launched nonbonded tasks and data
 * transfers to finish.
 *
 * Also does timing accounting and reduction of the internal staging buffers.
 * As this is called at the end of the step, it also resets the pair list and
 * pruning flags.
 *
 * \param[in] nb The nonbonded data GPU structure
 * \param[in]  stepWork     Force schedule flags
 * \param[in] aloc Atom locality identifier
 * \param[out] e_lj Pointer to the LJ energy output to accumulate into
 * \param[out] e_el Pointer to the electrostatics energy output to accumulate into
 * \param[out] shiftForces Shift forces buffer to accumulate into
 * \param[out] wcycle Pointer to wallcycle data structure
 * \return            The number of cycles the gpu wait took
 */
//NOLINTNEXTLINE(misc-definitions-in-headers) TODO: move into source file
float gpu_wait_finish_task(NbnxmGpu*                nb,
                           const gmx::StepWorkload& stepWork,
                           AtomLocality             aloc,
                           real*                    e_lj,
                           real*                    e_el,
                           gmx::ArrayRef<gmx::RVec> shiftForces,
                           gmx_wallcycle*           wcycle)
{
    auto cycleCounter = (gpuAtomToInteractionLocality(aloc) == InteractionLocality::Local)
                                ? ewcWAIT_GPU_NB_L
                                : ewcWAIT_GPU_NB_NL;

    wallcycle_start(wcycle, cycleCounter);
    gpu_try_finish_task(nb, stepWork, aloc, e_lj, e_el, shiftForces, GpuTaskCompletion::Wait, wcycle);
    float waitTime = wallcycle_stop(wcycle, cycleCounter);

    return waitTime;
}

} // namespace Nbnxm

#endif
