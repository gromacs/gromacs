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
 * \brief Implements common internal types for different NBNXN GPU implementations
 *
 * \author Szilárd Páll <pall.szilard@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H
#define GMX_MDLIB_NBNXN_GPU_COMMON_TYPES_H

#include "config.h"

#if GMX_GPU == GMX_GPU_OPENCL
#include "gromacs/gpu_utils/gpuregiontimer_ocl.h"
#endif

#if GMX_GPU == GMX_GPU_CUDA
#include "gromacs/gpu_utils/gpuregiontimer.cuh"
#endif

/*! \internal
 * \brief GPU region timers used for timing GPU kernels and H2D/D2H transfers.
 *
 * The two-sized arrays hold the local and non-local values and should always
 * be indexed with eintLocal/eintNonlocal.
 */
struct nbnxn_gpu_timers_t
{
    GpuRegionTimer atdat;              /**< timer for atom data transfer (every PS step)            */
    GpuRegionTimer nb_h2d[2];          /**< timer for x/q H2D transfers (l/nl, every step)          */
    GpuRegionTimer nb_d2h[2];          /**< timer for f D2H transfer (l/nl, every step)             */
    GpuRegionTimer pl_h2d[2];          /**< timer for pair-list H2D transfers (l/nl, every PS step) */
    bool           didPairlistH2D[2];  /**< true when a pair-list transfer has been done at this step */
    GpuRegionTimer nb_k[2];            /**< timer for non-bonded kernels (l/nl, every step)         */
    GpuRegionTimer prune_k[2];         /**< timer for the 1st pass list pruning kernel (l/nl, every PS step)   */
    bool           didPrune[2];        /**< true when we timed pruning and the timings need to be accounted for */
    GpuRegionTimer rollingPrune_k[2];  /**< timer for rolling pruning kernels (l/nl, frequency depends on chunk size)  */
    bool           didRollingPrune[2]; /**< true when we timed rolling pruning (at the previous step) and the timings need to be accounted for */
};

#endif
