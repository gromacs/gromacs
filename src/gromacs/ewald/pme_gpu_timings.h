/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018,2019, by the GROMACS development team, led by
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
 *  \brief Defines PME GPU timing functions.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_TIMINGS_H
#define GMX_EWALD_PME_GPU_TIMINGS_H

#include "config.h"

#if GMX_GPU == GMX_GPU_CUDA
#    include "gromacs/gpu_utils/gputraits.cuh"
#elif GMX_GPU == GMX_GPU_OPENCL
#    include "gromacs/gpu_utils/gputraits_ocl.h"
#endif

struct PmeGpu;

/*! \libinternal \brief
 * Starts timing the certain PME GPU stage during a single computation (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] PMEStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
void pme_gpu_start_timing(const PmeGpu* pmeGpu, size_t PMEStageId);

/*! \libinternal \brief
 * Returns raw timing event from the corresponding GpuRegionTimer (if timings are enabled).
 * In CUDA result can be nullptr stub, per GpuRegionTimer implementation.
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] PMEStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
CommandEvent* pme_gpu_fetch_timing_event(const PmeGpu* pmeGpu, size_t PMEStageId);

/*! \libinternal \brief
 * Stops timing the certain PME GPU stage during a single computation (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] PMEStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
void pme_gpu_stop_timing(const PmeGpu* pmeGpu, size_t PMEStageId);

#endif
