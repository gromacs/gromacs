/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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

#include "gmxpre.h"

#include "dlbtiming.h"

#include <cstring>

#include <array>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "domdec_internal.h"

/*! \brief Struct for timing the region for dynamic load balancing */
class BalanceRegion::Impl
{
public:
    /*! \brief Constructor */
    Impl() :
        isOpen(false), isOpenOnCpu(false), isOpenOnGpu(false), cyclesOpenCpu(0), cyclesLastCpu(0)
    {
    }

    bool         isOpen;        /**< Are we in an open balancing region? */
    bool         isOpenOnCpu;   /**< Is the, currently open, region still open on the CPU side? */
    bool         isOpenOnGpu;   /**< Is the, currently open, region open on the GPU side? */
    gmx_cycles_t cyclesOpenCpu; /**< Cycle count when opening the CPU region */
    gmx_cycles_t cyclesLastCpu; /**< Cycle count at the last call to \p ddCloseBalanceRegionCpu() */
};

BalanceRegion::BalanceRegion()
{
    impl_ = std::make_unique<Impl>();
}

BalanceRegion::~BalanceRegion() = default;

/*! \brief Returns the pointer to the balance region.
 *
 * This should be replaced by a properly managed BalanceRegion class,
 * but that requires a lot of refactoring in domdec.cpp.
 */
static BalanceRegion::Impl* getBalanceRegion(const gmx_domdec_t* dd)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion::Impl* region = dd->comm->balanceRegion.impl_.get();
    GMX_ASSERT(region != nullptr, "Balance region should be initialized before use");
    return region;
}

void DDBalanceRegionHandler::openRegionCpuImpl(DdAllowBalanceRegionReopen gmx_unused allowReopen) const
{
    BalanceRegion::Impl* reg = getBalanceRegion(dd_);
    if (dd_->comm->ddSettings.recordLoad)
    {
        GMX_ASSERT(allowReopen == DdAllowBalanceRegionReopen::yes || !reg->isOpen,
                   "Should not open an already opened region");

        reg->cyclesOpenCpu = gmx_cycles_read();
        reg->isOpen        = true;
        reg->isOpenOnCpu   = true;
        reg->isOpenOnGpu   = false;
    }
}

void DDBalanceRegionHandler::openRegionGpuImpl() const
{
    BalanceRegion::Impl* reg = getBalanceRegion(dd_);
    GMX_ASSERT(reg->isOpen, "Can only open a GPU region inside an open CPU region");
    GMX_ASSERT(!reg->isOpenOnGpu, "Can not re-open a GPU balance region");
    reg->isOpenOnGpu = true;
}

void ddReopenBalanceRegionCpu(const gmx_domdec_t* dd)
{
    BalanceRegion::Impl* reg = getBalanceRegion(dd);
    /* If the GPU is busy, don't reopen as we are overlapping with work */
    if (reg->isOpen && !reg->isOpenOnGpu)
    {
        reg->cyclesOpenCpu = gmx_cycles_read();
    }
}

void DDBalanceRegionHandler::closeRegionCpuImpl() const
{
    BalanceRegion::Impl* reg = getBalanceRegion(dd_);
    if (reg->isOpen && reg->isOpenOnCpu)
    {
        GMX_ASSERT(reg->isOpenOnCpu, "Can only close an open region");
        gmx_cycles_t cycles = gmx_cycles_read();
        reg->isOpenOnCpu    = false;

        if (reg->isOpenOnGpu)
        {
            /* Store the cycles for estimating the GPU/CPU overlap time */
            reg->cyclesLastCpu = cycles;
        }
        else
        {
            /* We can close the region */
            float cyclesCpu = cycles - reg->cyclesOpenCpu;
            dd_cycles_add(dd_, cyclesCpu, ddCyclF);
            reg->isOpen = false;
        }
    }
}

void DDBalanceRegionHandler::closeRegionGpuImpl(float waitGpuCyclesInCpuRegion,
                                                DdBalanceRegionWaitedForGpu waitedForGpu) const
{
    BalanceRegion::Impl* reg = getBalanceRegion(dd_);
    if (reg->isOpen)
    {
        GMX_ASSERT(reg->isOpenOnGpu, "Can not close a non-open GPU balance region");
        GMX_ASSERT(!reg->isOpenOnCpu,
                   "The GPU region should be closed after closing the CPU region");

        float waitGpuCyclesEstimate = gmx_cycles_read() - reg->cyclesLastCpu;
        if (waitedForGpu == DdBalanceRegionWaitedForGpu::no)
        {
            /* The actual time could be anywhere between 0 and
             * waitCyclesEstimate. Using half is the best we can do.
             */
            const float unknownWaitEstimateFactor = 0.5F;
            waitGpuCyclesEstimate *= unknownWaitEstimateFactor;
        }

        float cyclesCpu = reg->cyclesLastCpu - reg->cyclesOpenCpu;
        dd_cycles_add(dd_, cyclesCpu + waitGpuCyclesEstimate, ddCyclF);

        /* Register the total GPU wait time, to redistribute with GPU sharing */
        dd_cycles_add(dd_, waitGpuCyclesInCpuRegion + waitGpuCyclesEstimate, ddCyclWaitGPU);

        /* Close the region */
        reg->isOpenOnGpu = false;
        reg->isOpen      = false;
    }
}

//! Accumulates flop counts for force calculations.
static double force_flop_count(const t_nrnb* nrnb)
{
    double sum = 0;
    for (int i = 0; i < eNR_NBKERNEL_FREE_ENERGY; i++)
    {
        /* To get closer to the real timings, we half the count
         * for the normal loops and again half it for water loops.
         */
        const char* name = nrnb_str(i);
        if (strstr(name, "W3") != nullptr || strstr(name, "W4") != nullptr)
        {
            sum += nrnb->n[i] * 0.25 * cost_nrnb(i);
        }
        else
        {
            sum += nrnb->n[i] * 0.50 * cost_nrnb(i);
        }
    }
    for (int i = eNR_NBKERNEL_FREE_ENERGY; i <= eNR_NB14; i++)
    {
        const char* name = nrnb_str(i);
        if (strstr(name, "W3") != nullptr || strstr(name, "W4") != nullptr)
        {
            sum += nrnb->n[i] * cost_nrnb(i);
        }
    }
    for (int i = eNR_BONDS; i <= eNR_WALLS; i++)
    {
        sum += nrnb->n[i] * cost_nrnb(i);
    }

    return sum;
}

void dd_force_flop_start(gmx_domdec_t* dd, t_nrnb* nrnb)
{
    if (dd->comm->ddSettings.eFlop)
    {
        dd->comm->flop -= force_flop_count(nrnb);
    }
}

void dd_force_flop_stop(gmx_domdec_t* dd, t_nrnb* nrnb)
{
    if (dd->comm->ddSettings.eFlop)
    {
        dd->comm->flop += force_flop_count(nrnb);
        dd->comm->flop_n++;
    }
}

void clear_dd_cycle_counts(gmx_domdec_t* dd)
{
    for (int i = 0; i < ddCyclNr; i++)
    {
        dd->comm->cycl[i]     = 0;
        dd->comm->cycl_n[i]   = 0;
        dd->comm->cycl_max[i] = 0;
    }
    dd->comm->flop   = 0;
    dd->comm->flop_n = 0;
}
