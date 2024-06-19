/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *
 * \brief This file implements functions to interact with the dynamic load
 * balancing machinery.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "dlb.h"

#include <array>
#include <memory>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "domdec_internal.h"
#include "utility.h"

float dd_pme_f_ratio(const gmx_domdec_t* dd)
{
    GMX_ASSERT(DDMAIN(dd), "This function should only be called on the main rank");

    if (dd->comm->load[0].mdf > 0 && dd->comm->cycl_n[ddCyclPME] > 0)
    {
        return dd->comm->load[0].pme / dd->comm->load[0].mdf;
    }
    else
    {
        return -1.0;
    }
}

void set_dlb_limits(gmx_domdec_t* dd)

{
    for (int d = 0; d < dd->ndim; d++)
    {
        /* Set the number of pulses to the value for DLB */
        dd->comm->cd[d].ind.resize(dd->comm->cd[d].np_dlb);

        dd->comm->cellsize_min[dd->dim[d]] = dd->comm->cellsize_min_dlb[dd->dim[d]];
    }
}

void dd_dlb_set_should_check_whether_to_turn_dlb_on(gmx_domdec_t* dd, gmx_bool bValue)
{
    if (dd->comm->dlbState == DlbState::offCanTurnOn)
    {
        dd->comm->bCheckWhetherToTurnDlbOn = bValue;

        if (bValue)
        {
            /* Store the DD partitioning count, so we can ignore cycle counts
             * over the next nstlist steps, which are often slower.
             */
            dd->comm->ddPartioningCountFirstDlbOff = dd->ddp_count;
        }
    }
}

gmx_bool dd_dlb_get_should_check_whether_to_turn_dlb_on(gmx_domdec_t* dd)
{
    if (dd->comm->dlbState != DlbState::offCanTurnOn)
    {
        return FALSE;
    }

    if (dd->ddp_count <= dd->comm->ddPartioningCountFirstDlbOff)
    {
        /* We ignore the first nstlist steps at the start of the run
         * or after PME load balancing or after turning DLB off, since
         * these often have extra allocation or cache miss overhead.
         */
        return FALSE;
    }

    if (dd->comm->cycl_n[ddCyclStep] == 0)
    {
        /* We can have zero timed steps when dd_partition_system is called
         * more than once at the same step, e.g. with replica exchange.
         * Turning on DLB would trigger an assertion failure later, but is
         * also useless right after exchanging replicas.
         */
        return FALSE;
    }

    /* We should check whether we should use DLB directly after
     * unlocking DLB. */
    if (dd->comm->bCheckWhetherToTurnDlbOn)
    {
        /* This flag was set when the PME load-balancing routines
           unlocked DLB, and should now be cleared. */
        dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, FALSE);
        return TRUE;
    }
    /* We check whether we should use DLB every c_checkTurnDlbOnInterval
     * partitionings (we do not do this every partioning, so that we
     * avoid excessive communication). */
    return dd->comm->n_load_have % c_checkTurnDlbOnInterval == c_checkTurnDlbOnInterval - 1;
}

gmx_bool dd_dlb_is_on(const gmx_domdec_t* dd)
{
    return isDlbOn(dd->comm->dlbState);
}

gmx_bool dd_dlb_is_locked(const gmx_domdec_t* dd)
{
    return (dd->comm->dlbState == DlbState::offTemporarilyLocked);
}

void dd_dlb_lock(gmx_domdec_t* dd)
{
    /* We can only lock the DLB when it is set to auto, otherwise don't do anything */
    if (dd->comm->dlbState == DlbState::offCanTurnOn)
    {
        dd->comm->dlbState = DlbState::offTemporarilyLocked;
    }
}

void dd_dlb_unlock(gmx_domdec_t* dd)
{
    /* We can only lock the DLB when it is set to auto, otherwise don't do anything */
    if (dd->comm->dlbState == DlbState::offTemporarilyLocked)
    {
        dd->comm->dlbState = DlbState::offCanTurnOn;
        dd_dlb_set_should_check_whether_to_turn_dlb_on(dd, TRUE);
    }
}
