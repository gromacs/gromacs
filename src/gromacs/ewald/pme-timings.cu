/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include <cuda.h>

#include "gromacs/ewald/pme.h"
#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"

/* The pme_gpu_timing class implementation */

pme_gpu_timing::pme_gpu_timing()
{
    initialized = false;
    reset();
}

pme_gpu_timing::~pme_gpu_timing()
{
    if (initialized)
    {
        cudaError_t stat;
        stat = cudaEventDestroy(event_start);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        stat = cudaEventDestroy(event_stop);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        initialized = false;
    }
}

void pme_gpu_timing::enable()
{
    if (!initialized)
    {
        cudaError_t stat;
        stat = cudaEventCreate(&event_start, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        stat = cudaEventCreate(&event_stop, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        initialized = true;
    }
}

void pme_gpu_timing::start_recording(cudaStream_t s)
{
    if (initialized)
    {
        cudaError_t stat = cudaEventRecord(event_start, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
    }
}

void pme_gpu_timing::stop_recording(cudaStream_t s)
{
    if (initialized)
    {
        cudaError_t stat = cudaEventRecord(event_stop, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
        call_count++;
    }
}

void pme_gpu_timing::reset()
{
    total_milliseconds = 0.0;
    call_count         = 0;
}

void pme_gpu_timing::update()
{
    if (initialized && (call_count > 0)) /* Only the touched events needed */
    {
        real        milliseconds = 0.0;
        cudaError_t stat         = cudaEventElapsedTime(&milliseconds, event_start, event_stop);
        CU_RET_ERR(stat, "PME timing cudaEventElapsedTime fail");
        total_milliseconds += milliseconds;
    }
}

real pme_gpu_timing::get_total_time_milliseconds()
{
    return total_milliseconds;
}

unsigned int pme_gpu_timing::get_call_count()
{
    return call_count;
}

/* The general PME GPU timing functions */

/*! \brief \internal
 * Tells if CUDA-based performance tracking is enabled for PME.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if timings are enabled, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_timings_enabled(const gmx_pme_t *pme)
{
    return pme_gpu_enabled(pme) && pme->gpu->archSpecific->bTiming;
}

void pme_gpu_start_timing(const gmx_pme_t *pme, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pme))
    {
        GMX_ASSERT(PMEStageId < gtPME_EVENT_COUNT, "Wrong PME GPU timing event index");
        pme->gpu->archSpecific->timingEvents[PMEStageId]->start_recording(pme->gpu->archSpecific->pmeStream);
    }
}

void pme_gpu_stop_timing(const gmx_pme_t *pme, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pme))
    {
        GMX_ASSERT(PMEStageId < gtPME_EVENT_COUNT, "Wrong PME GPU timing event index");
        pme->gpu->archSpecific->timingEvents[PMEStageId]->stop_recording(pme->gpu->archSpecific->pmeStream);
    }
}

void pme_gpu_get_timings(const gmx_pme_t *pme, gmx_wallclock_gpu_t **timings)
{
    if (pme_gpu_timings_enabled(pme))
    {
        GMX_RELEASE_ASSERT(timings, "Null GPU timing pointer");
        if (!*timings)
        {
            // alloc for PME-only run
            snew(*timings, 1);
            // init_timings(*timings);
            // frankly, it's just memset..
        }
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            (*timings)->pme.timing[i].t = pme->gpu->archSpecific->timingEvents[i]->get_total_time_milliseconds();
            (*timings)->pme.timing[i].c = pme->gpu->archSpecific->timingEvents[i]->get_call_count();
        }
    }
}

void pme_gpu_update_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_timings_enabled(pme))
    {
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            pme->gpu->archSpecific->timingEvents[i]->update();
        }
    }
}

void pme_gpu_init_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_timings_enabled(pme))
    {
        cudaStreamSynchronize(pme->gpu->archSpecific->pmeStream);
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            pme->gpu->archSpecific->timingEvents[i] = new pme_gpu_timing();
            pme->gpu->archSpecific->timingEvents[i]->enable();
        }
    }
}

void pme_gpu_destroy_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_timings_enabled(pme))
    {
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            delete pme->gpu->archSpecific->timingEvents[i];
        }
        memset(pme->gpu->archSpecific->timingEvents, 0, sizeof(pme->gpu->archSpecific->timingEvents));
    }
}

void pme_gpu_reset_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_timings_enabled(pme))
    {
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            pme->gpu->archSpecific->timingEvents[i]->reset();
        }
    }
}
