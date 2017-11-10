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
gmx_inline bool pme_gpu_timings_enabled(const PmeGpu *pmeGpu)
{
    return pmeGpu->archSpecific->useTiming;
}

void pme_gpu_start_timing(const PmeGpu *pmeGpu, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        GMX_ASSERT(PMEStageId < pmeGpu->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGpu->archSpecific->timingEvents[PMEStageId].openTimingRegion(pmeGpu->archSpecific->pmeStream);
    }
}

void pme_gpu_stop_timing(const PmeGpu *pmeGpu, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        GMX_ASSERT(PMEStageId < pmeGpu->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGpu->archSpecific->timingEvents[PMEStageId].closeTimingRegion(pmeGpu->archSpecific->pmeStream);
    }
}

void pme_gpu_get_timings(const PmeGpu *pmeGpu, gmx_wallclock_gpu_pme_t *timings)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        GMX_RELEASE_ASSERT(timings, "Null GPU timing pointer");
        for (size_t i = 0; i < pmeGpu->archSpecific->timingEvents.size(); i++)
        {
            timings->timing[i].t = pmeGpu->archSpecific->timingEvents[i].getTotalTime();
            timings->timing[i].c = pmeGpu->archSpecific->timingEvents[i].getCallCount();
        }
    }
}

void pme_gpu_update_timings(const PmeGpu *pmeGpu)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        for (const size_t &activeTimer : pmeGpu->archSpecific->activeTimers)
        {
            pmeGpu->archSpecific->timingEvents[activeTimer].getLastRangeTime();
        }
    }
}

void pme_gpu_reinit_timings(const PmeGpu *pmeGpu)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        pmeGpu->archSpecific->activeTimers.clear();
        pmeGpu->archSpecific->activeTimers.insert(gtPME_SPLINEANDSPREAD);
        // TODO: no separate gtPME_SPLINE and gtPME_SPREAD as they are not used currently
        if (pme_gpu_performs_FFT(pmeGpu))
        {
            pmeGpu->archSpecific->activeTimers.insert(gtPME_FFT_C2R);
            pmeGpu->archSpecific->activeTimers.insert(gtPME_FFT_R2C);
        }
        if (pme_gpu_performs_solve(pmeGpu))
        {
            pmeGpu->archSpecific->activeTimers.insert(gtPME_SOLVE);
        }
        if (pme_gpu_performs_gather(pmeGpu))
        {
            pmeGpu->archSpecific->activeTimers.insert(gtPME_GATHER);
        }
    }
}

void pme_gpu_reset_timings(const PmeGpu *pmeGpu)
{
    if (pme_gpu_timings_enabled(pmeGpu))
    {
        for (size_t i = 0; i < pmeGpu->archSpecific->timingEvents.size(); i++)
        {
            pmeGpu->archSpecific->timingEvents[i].reset();
        }
    }
}
