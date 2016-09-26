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
 * \brief Implements high-level PME GPU functions which do not require GPU framework-specific code.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "gromacs/ewald/pme.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "pme-gpu-internal.h"
#include "pme-grid.h"
#include "pme-internal.h"
#include "pme-solve.h"

bool pme_gpu_task_enabled(const gmx_pme_t *pme)
{
    return (pme != nullptr) && (pme->runMode != PmeRunMode::CPU);
}

void pme_gpu_reset_timings(const gmx_pme_t *pme)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_reset_timings(pme->gpu);
    }
}

void pme_gpu_get_timings(const gmx_pme_t *pme, gmx_wallclock_gpu_pme_t *timings)
{
    if (pme_gpu_active(pme))
    {
        pme_gpu_get_timings(pme->gpu, timings);
    }
}

void pme_gpu_get_results(const gmx_pme_t *pme,
                         gmx_wallcycle_t  wcycle,
                         matrix           vir_q,
                         real            *energy_q)
{
    GMX_ASSERT(pme_gpu_active(pme), "This should be a GPU run of PME but it is not enabled.");

    const bool haveComputedEnergyAndVirial = pme->gpu->settings.stepFlags & GMX_PME_CALC_ENER_VIR;
    const bool haveComputedForces          = pme->gpu->settings.stepFlags & GMX_PME_CALC_F;

    // The wallcycle hierarchy is messy - we pause ewcFORCE just to exclude the PME GPU sync time
    wallcycle_stop(wcycle, ewcFORCE);

    wallcycle_start(wcycle, ewcWAIT_GPU_PME_GATHER);
    pme_gpu_finish_step(pme->gpu, haveComputedForces, haveComputedEnergyAndVirial);
    wallcycle_stop(wcycle, ewcWAIT_GPU_PME_GATHER);

    wallcycle_start_nocount(wcycle, ewcFORCE);

    if (haveComputedEnergyAndVirial)
    {
        if (pme->doCoulomb)
        {
            if (pme_gpu_performs_solve(pme->gpu))
            {
                pme_gpu_get_energy_virial(pme->gpu, energy_q, vir_q);
            }
            else
            {
                get_pme_ener_vir_q(pme->solve_work, pme->nthread, energy_q, vir_q);
            }
        }
        else
        {
            *energy_q = 0;
        }
    }
    /* No additional haveComputedForces code since forces are copied to the output host buffer with no transformation. */
}
