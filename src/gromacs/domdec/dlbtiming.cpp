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

void dd_openBalanceRegion(const gmx_domdec_t      *dd,
                          DdBalanceRegionUsingGpu  usingGpu,
                          bool gmx_unused          allowReopen)
{
    GMX_ASSERT(dd != NULL, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    GMX_ASSERT(allowReopen || !reg.isOpen, "Should not open an already opened region");
    if (dd->comm->bRecordLoad)
    {
        reg.cyclesOpen  = gmx_cycles_read();
        reg.usingGpu    = usingGpu;
        reg.isOpen      = true;
        reg.isOpenOnCpu = true;
    }
}

void dd_reopenBalanceRegion(const gmx_domdec_t *dd)
{
    GMX_ASSERT(dd != NULL, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    if (reg.isOpen)
    {
        reg.cyclesOpen = gmx_cycles_read();
    }
}

void dd_closeBalanceRegionCpu(const gmx_domdec_t *dd)
{
    GMX_ASSERT(dd != NULL, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    if (reg.isOpen)
    {
        gmx_cycles_t cycles   = gmx_cycles_read();

        if (reg.isOpenOnCpu)
        {
            reg.isOpenOnCpu   = false;
            float cyclesCpu   = cycles - reg.cyclesOpen;
            dd_cycles_add(dd, cyclesCpu, ddCyclF);
        }

        if (dd->comm->balanceRegion.usingGpu == DdBalanceRegionUsingGpu::yes)
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

void dd_closeBalanceRegionGpu(const gmx_domdec_t          *dd,
                              float                        waitCyclesToAdd,
                              DdBalanceRegionWaitedForGpu  waitedForGpu)
{
    GMX_ASSERT(dd != NULL, "Balance regions should only be used with DD");
    BalanceRegion &reg = dd->comm->balanceRegion;
    GMX_ASSERT(reg.usingGpu == DdBalanceRegionUsingGpu::yes, "The close GPU region function should only be called when we registered the use of a GPU");
    if (reg.isOpen)
    {
        GMX_ASSERT(!reg.isOpenOnCpu, "The GPU region should be closed after closing the CPU region");
        float waitCyclesEstimate = gmx_cycles_read() - reg.cyclesLastCpu;
        if (waitedForGpu == DdBalanceRegionWaitedForGpu::no)
        {
            /* The actual time could be anywhere between 0 and
             * waitCyclesEstimate. Using half is the best we can do.
             */
            waitCyclesEstimate *= 0.5f;
        }
        waitCyclesToAdd += waitCyclesEstimate;
        /* Register the GPU wait time, need to rebalance with GPU sharing */
        dd_cycles_add(dd, waitCyclesToAdd, ddCyclWaitGPU);
        /* Close the region */
        reg.isOpen = false;
    }
}
