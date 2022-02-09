/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 *  \brief Defines PME GPU timing functions.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_GPU_TIMINGS_H
#define GMX_EWALD_PME_GPU_TIMINGS_H

#include <cstddef>

struct gmx_wallclock_gpu_pme_t;
struct PmeGpu;

enum class PmeStage : int;

/*! \libinternal \brief
 * Starts timing the certain PME GPU stage during a single computation (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] pmeStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
void pme_gpu_start_timing(const PmeGpu* pmeGpu, PmeStage pmeStageId);

/*! \libinternal \brief
 * Stops timing the certain PME GPU stage during a single computation (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \param[in] pmeStageId     The PME GPU stage gtPME_ index from the enum in src/gromacs/timing/gpu_timing.h
 */
void pme_gpu_stop_timing(const PmeGpu* pmeGpu, PmeStage pmeStageId);

/*! \brief
 * Tells if CUDA-based performance tracking is enabled for PME.
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 * \returns                  True if timings are enabled, false otherwise.
 */
bool pme_gpu_timings_enabled(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Finalizes all the active PME GPU stage timings for the current computation. Should be called at the end of every computation.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
void pme_gpu_update_timings(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Updates the internal list of active PME GPU stages (if timings are enabled).
 *
 * \param[in] pmeGpu         The PME GPU data structure.
 */
void pme_gpu_reinit_timings(const PmeGpu* pmeGpu);

/*! \brief
 * Resets the PME GPU timings. To be called at the reset MD step.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 */
void pme_gpu_reset_timings(const PmeGpu* pmeGpu);

/*! \libinternal \brief
 * Copies the PME GPU timings to the gmx_wallclock_gpu_t structure (for log output). To be called at the run end.
 *
 * \param[in] pmeGpu         The PME GPU structure.
 * \param[in] timings        The gmx_wallclock_gpu_pme_t structure.
 */
void pme_gpu_get_timings(const PmeGpu* pmeGpu, gmx_wallclock_gpu_pme_t* timings);

#endif
