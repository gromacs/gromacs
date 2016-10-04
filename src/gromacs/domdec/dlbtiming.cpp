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

#include "gmxpre.h"

#include "dlbtiming.h"

#include "gromacs/domdec/domdec.h"
#include "gromacs/utility/gmxassert.h"

#include "domdec_internal.h"

void ddOpenBalanceRegionCpu(const gmx_domdec_t *dd,
                            bool gmx_unused     allowReopen)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    GMX_ASSERT(allowReopen || !reg.isOpen, "Should not open an already opened region");
    if (dd->comm->bRecordLoad)
    {
        reg.cyclesOpenCpu = gmx_cycles_read();
        reg.isOpen        = true;
        reg.isOpenOnCpu   = true;
        reg.isOpenOnGpu   = false;
    }
}

void ddOpenBalanceRegionGpu(const gmx_domdec_t *dd)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    if (reg.isOpen)
    {
        GMX_ASSERT(!reg.isOpenOnGpu, "Can not re-open a GPU balance region");
        reg.isOpenOnGpu = true;
    }
}

void ddReopenBalanceRegionCpu(const gmx_domdec_t *dd)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    /* If the GPU is busy, don't reopen as we are overlapping with work */
    if (reg.isOpen && !reg.isOpenOnGpu)
    {
        reg.cyclesOpenCpu = gmx_cycles_read();
    }
}

void ddCloseBalanceRegionCpu(const gmx_domdec_t *dd)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    if (reg.isOpen)
    {
        gmx_cycles_t cycles   = gmx_cycles_read();

        if (reg.isOpenOnCpu)
        {
            reg.isOpenOnCpu   = false;
            float cyclesCpu   = cycles - reg.cyclesOpenCpu;
            dd_cycles_add(dd, cyclesCpu, ddCyclF);
        }

        if (reg.isOpenOnGpu)
        {
            /* Store the cycles for estimating the GPU/CPU overlap time */
            reg.cyclesLastCpu = cycles;
        }
        else
        {
            /* We can close the region */
            reg.isOpen        = false;
        }
    }
}

void ddCloseBalanceRegionGpu(const gmx_domdec_t          *dd,
                             float                        waitCyclesToAdd,
                             DdBalanceRegionWaitedForGpu  waitedForGpu)
{
    GMX_ASSERT(dd != nullptr && dd->comm != nullptr, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    if (reg.isOpen)
    {
        GMX_ASSERT(reg.isOpenOnGpu, "Can not close a non-open GPU balance region");
        GMX_ASSERT(!reg.isOpenOnCpu, "The GPU region should be closed after closing the CPU region");

        float waitCyclesEstimate = gmx_cycles_read() - reg.cyclesLastCpu;
        if (waitedForGpu == DdBalanceRegionWaitedForGpu::no)
        {
            /* The actual time could be anywhere between 0 and
             * waitCyclesEstimate. Using half is the best we can do.
             */
            waitCyclesEstimate *= 0.5f;
        }

        float cyclesCpu  = reg.cyclesLastCpu - reg.cyclesOpenCpu;
        dd_cycles_add(dd, cyclesCpu + waitCyclesEstimate, ddCyclF);

        /* Register the total GPU wait time, to rebalance with GPU sharing */
        waitCyclesToAdd += waitCyclesEstimate;
        dd_cycles_add(dd, waitCyclesToAdd, ddCyclWaitGPU);

        /* Close the region */
        reg.isOpenOnGpu = false;
        reg.isOpen      = false;
    }
}
