/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
 *  \brief Implements PME GPU timing events in CUDA.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

#include "gmxpre.h"

#include "pme-timings.cuh"

#include "gromacs/utility/gmxassert.h"

#include "pme.cuh"

/*! \brief \internal
 * Tells if CUDA-based performance tracking is enabled for PME.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  True if timings are enabled, false otherwise.
 */
gmx_inline bool pme_gpu_timings_enabled(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->archSpecific->useTiming;
}

void pme_gpu_start_timing(const pme_gpu_t *pmeGPU, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_ASSERT(PMEStageId < pmeGPU->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGPU->archSpecific->timingEvents[PMEStageId].openTimingRegion(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_stop_timing(const pme_gpu_t *pmeGPU, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_ASSERT(PMEStageId < pmeGPU->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGPU->archSpecific->timingEvents[PMEStageId].closeTimingRegion(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_get_timings(const pme_gpu_t *pmeGPU, gmx_wallclock_gpu_pme_t *timings)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_RELEASE_ASSERT(timings, "Null GPU timing pointer");
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            timings->timing[i].t = pmeGPU->archSpecific->timingEvents[i].getTotalTime();
            timings->timing[i].c = pmeGPU->archSpecific->timingEvents[i].getCallCount();
        }
    }
}

void pme_gpu_update_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        pme_gpu_synchronize(pmeGPU);

        for (const size_t &activeTimer : pmeGPU->archSpecific->activeTimers)
        {
            pmeGPU->archSpecific->timingEvents[activeTimer].getLastRangeTime();
        }
    }
}

void pme_gpu_reinit_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        pmeGPU->archSpecific->activeTimers.clear();
        pmeGPU->archSpecific->activeTimers.insert(gtPME_SPLINEANDSPREAD);
        // TODO: no separate gtPME_SPLINE and gtPME_SPREAD as they are not used currently
        if (pme_gpu_performs_FFT(pmeGPU))
        {
            pmeGPU->archSpecific->activeTimers.insert(gtPME_FFT_C2R);
            pmeGPU->archSpecific->activeTimers.insert(gtPME_FFT_R2C);
        }
        if (pme_gpu_performs_solve(pmeGPU))
        {
            pmeGPU->archSpecific->activeTimers.insert(gtPME_SOLVE);
        }
        if (pme_gpu_performs_gather(pmeGPU))
        {
            pmeGPU->archSpecific->activeTimers.insert(gtPME_GATHER);
        }
    }
}

void pme_gpu_reset_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            pmeGPU->archSpecific->timingEvents[i].reset();
        }
    }
}
