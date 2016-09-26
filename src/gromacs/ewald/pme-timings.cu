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

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "pme.cuh"

/* The pme_gpu_timing class implementation */

pme_gpu_timing::pme_gpu_timing()
{
    initialized_ = false;
    reset();
}

pme_gpu_timing::~pme_gpu_timing()
{
    if (initialized_)
    {
        cudaError_t stat;
        stat = cudaEventDestroy(eventStart_);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        stat = cudaEventDestroy(eventStop_);
        CU_RET_ERR(stat, "PME timing cudaEventDestroy fail");
        initialized_ = false;
    }
}

void pme_gpu_timing::enable()
{
    if (!initialized_)
    {
        cudaError_t stat;
        stat = cudaEventCreate(&eventStart_, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        stat = cudaEventCreate(&eventStop_, cudaEventDefault);
        CU_RET_ERR(stat, "PME timing cudaEventCreate fail");
        initialized_ = true;
    }
}

void pme_gpu_timing::startRecording(cudaStream_t s)
{
    if (initialized_)
    {
        cudaError_t stat = cudaEventRecord(eventStart_, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
    }
}

void pme_gpu_timing::stopRecording(cudaStream_t s)
{
    if (initialized_)
    {
        cudaError_t stat = cudaEventRecord(eventStop_, s);
        CU_RET_ERR(stat, "PME timing cudaEventRecord fail");
        callCount_++;
    }
}

void pme_gpu_timing::reset()
{
    totalMilliseconds_ = 0.0;
    callCount_         = 0;
}

void pme_gpu_timing::update()
{
    if (initialized_ && (callCount_ > 0)) /* Only the touched events needed */
    {
        float        milliseconds = 0.0;
        cudaError_t  stat         = cudaEventElapsedTime(&milliseconds, eventStart_, eventStop_);
        CU_RET_ERR(stat, "PME timing cudaEventElapsedTime fail");
        totalMilliseconds_ += milliseconds;
    }
}

float pme_gpu_timing::getTotalTimeMilliseconds()
{
    return totalMilliseconds_;
}

unsigned int pme_gpu_timing::getCallCount()
{
    return callCount_;
}

/* The general PME GPU timing functions */

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
        pmeGPU->archSpecific->timingEvents[PMEStageId]->startRecording(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_stop_timing(const pme_gpu_t *pmeGPU, size_t PMEStageId)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_ASSERT(PMEStageId < pmeGPU->archSpecific->timingEvents.size(), "Wrong PME GPU timing event index");
        pmeGPU->archSpecific->timingEvents[PMEStageId]->stopRecording(pmeGPU->archSpecific->pmeStream);
    }
}

void pme_gpu_get_timings(const pme_gpu_t *pmeGPU, gmx_wallclock_gpu_pme_t *timings)
{
    if (pme_gpu_timings_enabled(pmeGPU))
    {
        GMX_RELEASE_ASSERT(timings, "Null GPU timing pointer");
        for (size_t i = 0; i < pmeGPU->archSpecific->timingEvents.size(); i++)
        {
            timings->timing[i].t = pmeGPU->archSpecific->timingEvents[i]->getTotalTimeMilliseconds();
            timings->timing[i].c = pmeGPU->archSpecific->timingEvents[i]->getCallCount();
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
