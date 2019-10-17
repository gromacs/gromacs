/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019, by the GROMACS development team, led by
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
 * \brief Implements common internal types for different NBNXN GPU implementations
 *
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \ingroup module_nbnxm
 */

#ifndef GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H
#define GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H

#include "config.h"

#include "gromacs/mdtypes/locality.h"
#include "gromacs/utility/enumerationhelpers.h"

#include "pairlist.h"

#if GMX_GPU == GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gpuregiontimer_ocl.h"
#endif

#if GMX_GPU == GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gpuregiontimer.cuh"
#endif

namespace Nbnxm
{

using gmx::AtomLocality;
using gmx::InteractionLocality;

/*! \internal
 * \brief GPU region timers used for timing GPU kernels and H2D/D2H transfers.
 *
 * The two-sized arrays hold the local and non-local values and should always
 * be indexed with eintLocal/eintNonlocal.
 */
struct gpu_timers_t
{
    /*! \internal
     * \brief Timers for local or non-local coordinate/force transfers
     */
    struct XFTransfers
    {
        GpuRegionTimer nb_h2d; /**< timer for x/q H2D transfers (l/nl, every step) */
        GpuRegionTimer nb_d2h; /**< timer for f D2H transfer (l/nl, every step) */
    };

    /*! \internal
     * \brief Timers for local or non-local interaction related operations
     */
    struct Interaction
    {
        GpuRegionTimer pl_h2d;       /**< timer for pair-list H2D transfers (l/nl, every PS step) */
        bool didPairlistH2D = false; /**< true when a pair-list transfer has been done at this step */
        GpuRegionTimer nb_k;         /**< timer for non-bonded kernels (l/nl, every step)         */
        GpuRegionTimer prune_k; /**< timer for the 1st pass list pruning kernel (l/nl, every PS step) */
        bool didPrune = false; /**< true when we timed pruning and the timings need to be accounted for */
        GpuRegionTimer rollingPrune_k; /**< timer for rolling pruning kernels (l/nl, frequency depends on chunk size)  */
        bool           didRollingPrune =
                false; /**< true when we timed rolling pruning (at the previous step) and the timings need to be accounted for */
    };

    //! timer for atom data transfer (every PS step)
    GpuRegionTimer atdat;
    //! timers for coordinate/force transfers (every step)
    gmx::EnumerationArray<AtomLocality, XFTransfers> xf;
    //! timers for interaction related transfers
    gmx::EnumerationArray<InteractionLocality, Nbnxm::gpu_timers_t::Interaction> interaction;
};

struct gpu_plist
{
    int na_c; /**< number of atoms per cluster                  */

    int                       nsci;       /**< size of sci, # of i clusters in the list     */
    int                       sci_nalloc; /**< allocation size of sci                       */
    DeviceBuffer<nbnxn_sci_t> sci;        /**< list of i-cluster ("super-clusters")         */

    int                       ncj4;          /**< total # of 4*j clusters                      */
    int                       cj4_nalloc;    /**< allocation size of cj4                       */
    DeviceBuffer<nbnxn_cj4_t> cj4;           /**< 4*j cluster list, contains j cluster number
                                                and index into the i cluster list            */
    int                        nimask;       /**< # of 4*j clusters * # of warps               */
    int                        imask_nalloc; /**< allocation size of imask                     */
    DeviceBuffer<unsigned int> imask;        /**< imask for 2 warps for each 4*j cluster group */
    DeviceBuffer<nbnxn_excl_t> excl;         /**< atom interaction bits                        */
    int                        nexcl;        /**< count for excl                               */
    int                        excl_nalloc;  /**< allocation size of excl                      */

    /* parameter+variables for normal and rolling pruning */
    bool haveFreshList; /**< true after search, indictes that initial pruning with outer prunning is needed */
    int  rollingPruningNumParts; /**< the number of parts/steps over which one cyle of roling pruning takes places */
    int  rollingPruningPart; /**< the next part to which the roling pruning needs to be applied */
};

} // namespace Nbnxm

#endif
