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

#include "pme-timings.cuh"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"

/* The pme_gpu_timing class implementation */

pme_gpu_timing::pme_gpu_timing()
{
    _initialized = FALSE;
    reset();
}

pme_gpu_timing::~pme_gpu_timing()
{
    if (_initialized)
    {
        cudaError_t stat;
        stat = cudaEventDestroy(_eventStart);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        stat = cudaEventDestroy(_eventStop);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        _initialized = FALSE;
    }
}

void pme_gpu_timing::enable()
{
    if (!_initialized)
    {
        cudaError_t stat;
        stat = cudaEventCreate(&_eventStart, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        stat = cudaEventCreate(&_eventStop, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        _initialized = TRUE;
    }
}

void pme_gpu_timing::start_recording(cudaStream_t s)
{
    if (_initialized)
    {
        cudaError_t stat = cudaEventRecord(_eventStart, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
    }
}

void pme_gpu_timing::stop_recording(cudaStream_t s)
{
    if (_initialized)
    {
        cudaError_t stat = cudaEventRecord(_eventStop, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
        _callCount++;
    }
}

void pme_gpu_timing::reset()
{
    _totalMilliseconds = 0.0;
    _callCount         = 0;
}

void pme_gpu_timing::update()
{
    if (_initialized && (_callCount > 0)) /* Only the touched events needed */
    {
        float        milliseconds = 0.0;
        cudaError_t  stat         = cudaEventElapsedTime(&milliseconds, _eventStart, _eventStop);
        CU_RET_ERR(stat, "PME timing cudaEventElapsedTime fail");
        _totalMilliseconds += milliseconds;
    }
}

float pme_gpu_timing::get_total_time_milliseconds()
{
    return _totalMilliseconds;
}

unsigned int pme_gpu_timing::get_call_count()
{
    return _callCount;
}

/* The general PME GPU timing functions */

/*! \brief \internal
 * Tells if CUDA-based performance tracking is enabled for PME.
 *
 * \param[in] pme            The PME data structure.
 * \returns                  TRUE if timings are enabled, FALSE otherwise.
 */
gmx_inline gmx_bool pme_gpu_timings_enabled(const pme_gpu_t *pmeGPU)
{
    return pmeGPU->archSpecific->bTiming;
}

void pme_gpu_start_timing(const pme_gpu_t *pmeGPU, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_ASSERT(PMEStageId < pmeGPU->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGPU->archSpecific->timingEvents[PMEStageId]->start_recording(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_stop_timing(const pme_gpu_t *pmeGPU, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_ASSERT(PMEStageId < pmeGPU->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGPU->archSpecific->timingEvents[PMEStageId]->stop_recording(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_get_timings(const pme_gpu_t *pmeGPU, gmx_wallclock_gpu_t **timings)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_RELEASE_ASSERT(timings, "Null GPU timing pointer");
        if (!*timings)
        {
            // alloc for PME-only run
            snew(*timings, 1);
            /* FIXME: this is not freed, should be shared */
            // init_timings(*timings);
        }
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            (*timings)->pme.timing[i].t = pmeGPU->archSpecific->timingEvents[i]->get_total_time_milliseconds();
            (*timings)->pme.timing[i].c = pmeGPU->archSpecific->timingEvents[i]->get_call_count();
        }
    }
}

void pme_gpu_update_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            pmeGPU->archSpecific->timingEvents[i]->update();
        }
    }
}

void pme_gpu_init_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        pme_gpu_synchronize(pmeGPU);
        for (size_t i = 0; i < gtPME_EVENT_COUNT; i++)
        {
            pmeGPU->archSpecific->timingEvents.push_back(std::unique_ptr<pme_gpu_timing>(new pme_gpu_timing()));
            pmeGPU->archSpecific->timingEvents[i]->enable();
        }
    }
}

void pme_gpu_destroy_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        pmeGPU->archSpecific->timingEvents.resize(0);
    }
}

void pme_gpu_reset_timings(const pme_gpu_t *pmeGPU)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            pmeGPU->archSpecific->timingEvents[i]->reset();
        }
    }
}
