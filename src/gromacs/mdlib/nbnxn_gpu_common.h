/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Implements common routines for different NBNXN GPU implementations
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_NBNXN_GPU_COMMON_H
#define GMX_MDLIB_NBNXN_GPU_COMMON_H

#include "config.h"

#include <string>

#if GMX_GPU == GMX_GPU_CUDA
#include "nbnxn_cuda/nbnxn_cuda_types.h"
#endif

#if GMX_GPU == GMX_GPU_OPENCL
#include "nbnxn_ocl/nbnxn_ocl_types.h"
#endif

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/stringutil.h"

#include "nbnxn_gpu.h"


/*! \brief An early return condition for empty NB GPU workloads
 *
 * This is currently used for non-local kernels/transfers only.
 * Skipping the local kernel is more complicated, since the
 * local part of the force array also depends on the non-local kernel.
 * The skip of the local kernel is taken care of separately.
 */
static inline bool canSkipWork(const gmx_nbnxn_gpu_t *nb, int iloc)
{
    assert(nb && nb->plist[iloc]);
    return (iloc == eintNonlocal) && (nb->plist[iloc]->nsci == 0);
}

/*! \brief Check that atom locality values are valid for the GPU module.
 *
 *  In the GPU module atom locality "all" is not supported, the local and
 *  non-local ranges are treated separately.
 *
 *  \param[in] atomLocality atom locality specifier
 */
static inline void validateGpuAtomLocality(int atomLocality)
{
    std::string str = gmx::formatString("Invalid atom locality passed (%d); valid here is only "
                                        "local (%d) or nonlocal (%d)", atomLocality, eatLocal, eatNonlocal);

    GMX_ASSERT(LOCAL_OR_NONLOCAL_A(atomLocality), str.c_str());
}

/*! \brief Convert atom locality to interaction locality.
 *
 *  In the current implementation the this is straightforward conversion:
 *  local to local, non-local to non-local.
 *
 *  \param[in] atomLocality Atom locality specifier
 *  \returns                Interaction locality corresponding to the atom locality passed.
 */
static inline int gpuAtomToInteractionLocality(int atomLocality)
{
    validateGpuAtomLocality(atomLocality);

    /* determine interaction locality from atom locality */
    if (LOCAL_A(atomLocality))
    {
        return eintLocal;
    }
    else if (NONLOCAL_A(atomLocality))
    {
        return eintNonlocal;
    }
    else
    {
        // can't be reached
        assert(false);
        return -1;
    }
}

/*! \brief Calculate atom range and return start index and length.
 *
 * \param[in] atomData Atom descriptor data structure
 * \param[in] atomLocality Atom locality specifier
 * \param[out] atomRangeBegin Starting index of the atom range in the atom data array.
 * \param[out] atomRangeLen Atom range length in the atom data array.
 */
template <typename AtomDataT>
static inline void getGpuAtomRange(const AtomDataT *atomData,
                                   int              atomLocality,
                                   int             &atomRangeBegin,
                                   int             &atomRangeLen)
{
    assert(atomData);
    validateGpuAtomLocality(atomLocality);

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(atomLocality))
    {
        atomRangeBegin  = 0;
        atomRangeLen    = atomData->natoms_local;
    }
    else
    {
        atomRangeBegin  = atomData->natoms_local;
        atomRangeLen    = atomData->natoms - atomData->natoms_local;
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
template <typename GpuTimers>
static void countPruneKernelTime(GpuTimers                 *timers,
                                 gmx_wallclock_gpu_nbnxn_t *timings,
                                 const int                  iloc)
{
    // We might have not done any pruning (e.g. if we skipped with empty domains).
    if (!timers->didPrune[iloc] && !timers->didRollingPrune[iloc])
    {
        return;
    }

    if (timers->didPrune[iloc])
    {
        timings->pruneTime.c++;
        timings->pruneTime.t += timers->prune_k[iloc].getLastRangeTime();
    }

    if (timers->didRollingPrune[iloc])
    {
        timings->dynamicPruneTime.c++;
        timings->dynamicPruneTime.t += timers->rollingPrune_k[iloc].getLastRangeTime();
    }
}

template <typename GpuTimers, typename GpuPairlist>
static void nbnxn_gpu_accumulate_timings(struct gmx_wallclock_gpu_nbnxn_t  *timings,
                                         GpuTimers                         *timers,
                                         const GpuPairlist                 *plist,
                                         int                                atomLocality,
                                         bool                               didEnergyKernels,
                                         bool                               doTiming)
{
    /* timing data accumulation */
    if (!doTiming)
    {
        return;
    }

    /* determine interaction locality from atom locality */
    int iLocality = gpuAtomToInteractionLocality(atomLocality);

    /* only increase counter once (at local F wait) */
    if (LOCAL_I(iLocality))
    {
        timings->nb_c++;
        timings->ktime[plist->haveFreshList ? 1 : 0][didEnergyKernels ? 1 : 0].c += 1;
    }

    /* kernel timings */
    timings->ktime[plist->haveFreshList ? 1 : 0][didEnergyKernels ? 1 : 0].t +=
        timers->nb_k[iLocality].getLastRangeTime();

    /* X/q H2D and F D2H timings */
    timings->nb_h2d_t += timers->nb_h2d[iLocality].getLastRangeTime();
    timings->nb_d2h_t += timers->nb_d2h[iLocality].getLastRangeTime();

    /* Count the pruning kernel times for both cases:1st pass (at search step)
       and rolling pruning (if called at the previous step).
       We do the accounting here as this is the only sync point where we
       know (without checking or additional sync-ing) that prune tasks in
       in the current stream have completed (having just blocking-waited
       for the force D2H). */
    countPruneKernelTime(timers, timings, iLocality);

    /* only count atdat and pair-list H2D at pair-search step */
    if (timers->didPairlistH2D[iLocality])
    {
        /* atdat transfer timing (add only once, at local F wait) */
        if (LOCAL_A(atomLocality))
        {
            timings->pl_h2d_c++;
            timings->pl_h2d_t += timers->atdat.getLastRangeTime();
        }

        timings->pl_h2d_t += timers->pl_h2d[iLocality].getLastRangeTime();

        /* Clear the timing flag for the next step */
        timers->didPairlistH2D[iLocality] = false;
    }
}

template <typename StagingData>
static void nbnxn_gpu_reduce_staged_outputs(const StagingData     &nbst,
                                            int                    iLocality,
                                            bool                   reduceEnergies,
                                            bool                   reduceFshift,
                                            real                  &e_lj,
                                            real                  &e_el,
                                            rvec                  *fshift)
{
    /* add up energies and shift forces (only once at local F wait) */
    if (LOCAL_I(iLocality))
    {
        if (reduceEnergies)
        {
            e_lj += *nbst.e_lj;
            e_el += *nbst.e_el;
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

template <typename NbnxnGpu>
void nbnxn_gpu_wait_for_gpu(NbnxnGpu *nb,
                            int       flags,
                            int       aloc,
                            real     *e_lj,
                            real     *e_el,
                            rvec     *fshift)
{
    /* determine interaction locality from atom locality */
    int iLocality = gpuAtomToInteractionLocality(aloc);

    /* Launch wait/update timers & counters and do reduction into staging buffers
       BUT skip it when during the non-local phase there was actually no work to do.
       This is consistent with nbnxn_gpu_launch_kernel.

       NOTE: if timing with multiple GPUs (streams) becomes possible, the
       counters could end up being inconsistent due to not being incremented
       on some of the nodes! */
    if (!canSkipWork(nb, iLocality))
    {
        gpuStreamSynchronize(nb->stream[iLocality]);

        bool calcEner   = flags & GMX_FORCE_ENERGY;
        bool calcFshift = flags & GMX_FORCE_VIRIAL;

        /* turn energy calculation always on/off (for debugging/testing only) */
        //calcEner = (calcEner || always_ener) && !never_ener;

        nbnxn_gpu_accumulate_timings(nb->timings, nb->timers, nb->plist[iLocality], aloc, calcEner, nb->bDoTime);

        nbnxn_gpu_reduce_staged_outputs(nb->nbst, iLocality, calcEner, calcFshift, *e_lj, *e_el, fshift);

    }

    /* Always reset both pruning flags (doesn't hurt doing it even when timing is off). */
    nb->timers->didPrune[iLocality] = nb->timers->didRollingPrune[iLocality] = false;

    /* Turn off initial list pruning (doesn't hurt if this is not pair-search step). */
    nb->plist[iLocality]->haveFreshList = false;
}

#if GMX_GPU == GMX_GPU_CUDA
template void nbnxn_gpu_wait_for_gpu<gmx_nbnxn_cuda_t>(gmx_nbnxn_cuda_t *, int, int, real *, real *, rvec *);
#endif
#if GMX_GPU == GMX_GPU_OPENCL
template void nbnxn_gpu_wait_for_gpu<gmx_nbnxn_ocl_t>(gmx_nbnxn_ocl_t *, int, int, real *, real *, rvec *);
#endif

#endif
