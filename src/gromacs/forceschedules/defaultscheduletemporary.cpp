/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements default force-calculation schedule
 *
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 *
 * \todo Later we will do a simple move of sim_util.cpp to become
 * defaultschedule.cpp so that git rebase will automatically be able
 * to move changes in work not yet in master. Then we will have a
 * further commit that changes the function declaration to match what
 * is currently in defaultscheduletemporary.cpp.
 *
 * \ingroup module_forceschedules
 */

#include "gmxpre.h"

#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/state.h"

#include "defaultschedule.h"

namespace gmx
{

// TODO Later the contents of do_force_cutsVERLET will implement this function.
void gmx::DefaultSchedule::computeForces(int flags)
{
    do_force(log_, cr_, ms_, inputrec_,
             awh_, enforcedRotation_, *(sp_->step_),
             nrnb_, wcycle_, top_, groups_,
             state_->box, state_->x.arrayRefWithPadding(), &(state_->hist),
             force_->arrayRefWithPadding(), *vir_force_,
             mdatoms_,
             enerd_, fcd_,
             state_->lambda, graph_,
             fr_, ppForceWorkload_,
             vsite_, *mu_tot_,
             *(sp_->t_), ed_,
             flags,
             *ddBalanceRegionHandler_);
}

}   // namespace gmx
