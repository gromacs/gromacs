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
/*! \libinternal \file
 *
 * \brief This file declares functions to interact with the dynamic load
 * balancing machinery.
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_DLB_H
#define GMX_DOMDEC_DLB_H

#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct t_commrec;


/*! \brief We check if to turn on DLB at the first and every 100 DD partitionings.
 * With large imbalance DLB will turn on at the first step, so we can
 * make the interval so large that the MPI overhead of the check is negligible.
 */
constexpr int c_checkTurnDlbOnInterval = 100;
/*! \brief We need to check if DLB results in worse performance and then turn it off.
 * We check this more often then for turning DLB on, because the DLB can scale
 * the domains very rapidly, so if unlucky the load imbalance can go up quickly
 * and furthermore, we are already synchronizing often with DLB, so
 * the overhead of the MPI Bcast is not that high.
 */
constexpr int c_checkTurnDlbOffInterval = 20;


/*! \brief Return the PME/PP force load ratio, or -1 if nothing was measured.
 *
 * Should only be called on the DD master node.
 */
float dd_pme_f_ratio(const gmx_domdec_t* dd);

//! Sets the cell size limits for DD to suit dynamic load balancing.
void set_dlb_limits(gmx_domdec_t* dd);

/*! \brief Limit DLB to preserve the option of returning to the current cut-off.
 *
 * Domain boundary changes due to the DD dynamic load balancing can limit
 * the cut-off distance that can be set in change_dd_cutoff. This function
 * sets/changes the DLB limit such that using the passed (pair-list) cut-off
 * should still be possible after subsequently setting a shorter cut-off
 * with change_dd_cutoff.
 */
void set_dd_dlb_max_cutoff(struct t_commrec* cr, real cutoff);

/*! \brief Sets whether we should later check the load imbalance data, so that
 * we can trigger dynamic load balancing if enough imbalance has
 * arisen.
 *
 * Used after PME load balancing unlocks DLB, so that the check
 * whether DLB will be useful can happen immediately.
 */
void dd_dlb_set_should_check_whether_to_turn_dlb_on(gmx_domdec_t* dd, bool bValue);

/*! \brief Returns if we should check whether there has been enough
 * load imbalance to trigger dynamic load balancing.
 *
 * We need to check whether we check because it might be always off.
 */
bool dd_dlb_get_should_check_whether_to_turn_dlb_on(gmx_domdec_t* dd);

/*! \brief Return if we are currently using dynamic load balancing */
bool dd_dlb_is_on(const gmx_domdec_t* dd);

/*! \brief Return if the DLB lock is set */
bool dd_dlb_is_locked(const gmx_domdec_t* dd);

/*! \brief Set a lock such that with DLB=auto DLB cannot get turned on */
void dd_dlb_lock(struct gmx_domdec_t* dd);

/*! \brief Clear a lock such that with DLB=auto DLB may get turned on later */
void dd_dlb_unlock(struct gmx_domdec_t* dd);

#endif
