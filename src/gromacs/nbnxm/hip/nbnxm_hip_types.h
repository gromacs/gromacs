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
 *  \brief
 *  Data types used internally in the nbnxn_hip module.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \author Paul Bauer <paul.bauer.q@gmail.com>
 *  \ingroup module_nbnxm
 */

#ifndef NBNXM_HIP_TYPES_H
#define NBNXM_HIP_TYPES_H

#include <memory>
#include <variant>

#include "gromacs/gpu_utils/gpueventsynchronizer.h"
#include "gromacs/nbnxm/gpu_types_common.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/enumerationhelpers.h"

namespace gmx
{

/*! \internal
 * \brief Main data structure for HIP nonbonded force calculations.
 */
struct NbnxmGpu
{
    /*! \brief GPU device context.
     *
     * \todo Make it constant reference, once NbnxmGpu is a proper class.
     */
    const DeviceContext* deviceContext_;
    /*! \brief true if doing both local/non-local NB work on GPU */
    bool bUseTwoStreams = false;
    //! true indicates that the nonlocal_done event was marked
    bool bNonLocalStreamDoneMarked = false;

    /*! \brief atom data */
    NBAtomDataGpu* atdat = nullptr;
    /*! \brief array of atom indices */
    int* atomIndices = nullptr;
    /*! \brief size of atom indices */
    int atomIndicesSize = 0;
    /*! \brief size of atom indices allocated in device buffer */
    int atomIndicesSize_alloc = 0;
    /*! \brief x buf ops num of atoms */
    int* cxy_na = nullptr;
    /*! \brief number of elements in cxy_na */
    int ncxy_na = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_na_alloc = 0;
    /*! \brief x buf ops cell index mapping */
    int* cxy_ind = nullptr;
    /*! \brief number of elements in cxy_ind */
    int ncxy_ind = 0;
    /*! \brief number of elements allocated allocated in device buffer */
    int ncxy_ind_alloc = 0;
    /*! \brief parameters required for the non-bonded calc. */
    NBParamGpu* nbparam = nullptr;
    /*! \brief pair-list data structures (local and non-local) */
    std::variant<EnumerationArray<InteractionLocality, std::unique_ptr<GpuPairlist<PairlistType::Hierarchical8x8x8>>>> plist;
    /*! \brief staging area where fshift/energies get downloaded */
    NBStagingData nbst;
    /*! \brief local and non-local GPU streams */
    EnumerationArray<InteractionLocality, const DeviceStream*> deviceStreams;

    //! true when a pair-list transfer has been done at this step
    EnumerationArray<InteractionLocality, bool> didPairlistH2D = { { false } };
    //! true when we we did pruning on this step
    EnumerationArray<InteractionLocality, bool> didPrune = { { false } };
    //! true when we did rolling pruning (at the previous step)
    EnumerationArray<InteractionLocality, bool> didRollingPrune = { { false } };

    /*! \brief Event triggered when the non-local non-bonded
     * kernel is done (and the local transfer can proceed) */
    GpuEventSynchronizer nonlocal_done;
    /*! \brief Event triggered when the tasks issued in the local
     * stream that need to precede the non-local force or buffer
     * operation calculations are done (e.g. f buffer 0-ing, local
     * x/q H2D, buffer op initialization in local stream that is
     * required also by nonlocal stream ) */
    GpuEventSynchronizer misc_ops_and_local_H2D_done;

    /*! \brief True if there is work for the current domain in the
     * respective locality.
     *
     * This includes local/nonlocal GPU work, either bonded or
     * nonbonded, scheduled to be executed in the current
     * domain. As long as bonded work is not split up into
     * local/nonlocal, if there is bonded GPU work, both flags
     * will be true. */
    EnumerationArray<InteractionLocality, bool> haveWork = { { false } };

    /*! \brief True if event-based timing is enabled. */
    bool bDoTime = false;
    /*! \brief HIP event-based timers. */
    GpuTimers* timers = nullptr;
    /*! \brief Timing data. TODO: deprecate this and query timers for accumulated data instead */
    gmx_wallclock_gpu_nbnxn_t* timings = nullptr;
};

} // namespace gmx

#endif /* NBNXN_HIP_TYPES_H */
