/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief
 * Implements classes from iforceprovider.h.
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \ingroup module_mdtypes
 */

#include "gmxpre.h"
#include "iforceschedule.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdtypes/state.h"

using namespace gmx;

void gmx::DefaultVerletSchedule::computeForces(int flags)
{
    /* modify force flag if not doing nonbonded */
    if (!fr_->bNonbonded)
    {
        flags &= ~GMX_FORCE_NONBONDED;
    }

    do_force_cutsVERLET(log_, cr_, ms_, inputrec_,
                        awh_, enforcedRotation_, *(sp_->step_),
                        nrnb_, wcycle_, top_, groups_,
                        state_->box, state_->x.arrayRefWithPadding(), &(state_->hist),
                        force_->arrayRefWithPadding(), *vir_force_,
                        mdatoms_,
                        enerd_, fcd_,
                        state_->lambda.data(), graph_,
                        fr_, ppForceWorkload_, fr_->ic,
                        vsite_, *mu_tot_,
                        *(sp_->t_), ed_,
                        flags,
                        *(sp_->ddOpenBalanceRegion_),
                        *(sp_->ddCloseBalanceRegion_));

    /* In case we don't have constraints and are using GPUs, the next balancing
     * region starts here.
     * Some "special" work at the end of do_force_cuts?, such as vsite spread,
     * virial calculation and COM pulling, is not thus not included in
     * the balance timing, which is ok as most tasks do communication.
     */
    if (*(sp_->ddOpenBalanceRegion_) == DdOpenBalanceRegionBeforeForceComputation::yes)
    {
        ddOpenBalanceRegionCpu(cr_->dd, DdAllowBalanceRegionReopen::no);
    }
}

void gmx::DefaultGroupSchedule::computeForces(int flags)
{
    /* modify force flag if not doing nonbonded */
    if (!fr_->bNonbonded)
    {
        flags &= ~GMX_FORCE_NONBONDED;
    }


    do_force_cutsGROUP(log_, cr_, ms_, inputrec_,
                       awh_, enforcedRotation_, *(sp_->step_),
                       nrnb_, wcycle_, top_, groups_,
                       state_->box, state_->x.arrayRefWithPadding(), &(state_->hist),
                       force_->arrayRefWithPadding(), *vir_force_,
                       mdatoms_,
                       enerd_, fcd_,
                       state_->lambda.data(), graph_,
                       fr_, vsite_, *mu_tot_,
                       *(sp_->t_), ed_,
                       flags,
                       *(sp_->ddOpenBalanceRegion_),
                       *(sp_->ddCloseBalanceRegion_));

    /* In case we don't have constraints and are using GPUs, the next balancing
     * region starts here.
     * Some "special" work at the end of do_force_cuts?, such as vsite spread,
     * virial calculation and COM pulling, is not thus not included in
     * the balance timing, which is ok as most tasks do communication.
     */
    if (*(sp_->ddOpenBalanceRegion_) == DdOpenBalanceRegionBeforeForceComputation::yes)
    {
        ddOpenBalanceRegionCpu(cr_->dd, DdAllowBalanceRegionReopen::no);
    }
}
