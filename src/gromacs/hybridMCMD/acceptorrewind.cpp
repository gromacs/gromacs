/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 * Implements the AcceptOrRewind class and an interface to a Metropolis step.
 *
 * \author Sebastian Wingbermuehle
 * \ingroup module_hybridMCMD
 */

#include "gmxpre.h"

#include "acceptorrewind.h"

#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/topology.h"

namespace gmx
{

AcceptOrRewind::AcceptOrRewind(IMetropolisStep *metropolisStep, t_state *globalState, t_commrec *cr) : metropolisStep_(metropolisStep),
                                                                                                       accepted_(false), // for testing
                                                                                                       totalProposedConfigurations_(0),
                                                                                                       acceptedProposedConfigurations_(0)
{
    // to make sure globalState_ is initialised properly
    if (MASTER(cr))
    {
        globalState_ = *globalState;
    }
}

/*! \brief Updates the class t_state according to the result of the Metropolis step.
 *         If the configuration is accepted, it makes a back-up; if the configuration is rejected, it rewinds to the last accepted configuration (and restarts domain decomposition).
 */
bool AcceptOrRewind::updateTState(const int64_t         step,                                  // arguments for core functionality
                                  t_state              *localState,
                                  gmx_enerdata_t       *enerd,
                                  FILE                 *fplog,                                 // arguments that are forwarded to domain decomposition functions
                                  const t_commrec      *cr,
                                  const gmx_mtop_t     *top_global,
                                  const t_inputrec     *ir,
                                  gmx_enfrot           *enforcedRotation,
                                  PaddedRVecVector     *f,
                                  MDAtoms              *mdAtoms,
                                  gmx_localtop_t       *top,
                                  t_forcerec           *fr,
                                  gmx_vsite_t          *vsite,
                                  Constraints          *constr,
                                  t_nrnb               *nrnb,
                                  gmx_wallcycle        *wcycle,
                                  gmx_update_t         *upd,
                                  bool                 *shouldCheckNumberOfBondedInteractions) // bool required for signalling in md.cpp
{
    /* compute_globals uses MPI_Allreduce
     * => all ranks/nodes/processors are on the same page as far as enerd is concerned
     * => can provide the same input data on all (PP) nodes
     * step is incremented on all (PP) nodes in do_md()
     * => random number reproducible
     * => evaluation of bool accepted_ can be done on all nodes (no global communication required)
     */
    accepted_ = metropolisStep_->accept(step, enerd);

    // back-up local state to avoid communication if configuration is accepted
    if (accepted_)
    {
        /*  deep copy of class is required
         *  TODO find solution for running averages in history_t
         */
        localStateBackUp_ = *localState;
        acceptedProposedConfigurations_++;

        /*  To collect the local states in a global state, the domain decomposition structure containing the mapping must be up-to-date
         *  -> a back-up of this structure is very complicated (too many important quantities are handled by pointers)
         *  => back-up global state here (although this requires an extra communication step)
         */
        if (DOMAINDECOMP(cr))
        {
            dd_collect_state(cr->dd, &localStateBackUp_, &globalState_);
        }
    }
    /* TODO back-up state of domain decomposition, dynamic loadbalancing, pme-tuning (and neighbour searching) to allow for local rewinding
     *      -> dynamic loadbalancing and pme-tuning are currently disabled if hybrid MC/MD is used
     * in case of rejection, construct global state from localStateBackUp, rewind and start all over with dd
     */
    else
    {
        if (DOMAINDECOMP(cr))
        {
            /* explanations of parameter values:
               - bMasterState = TRUE: we want to start over new with domain decomposition after rewinding
               - nstglobalcomm = 1: only seems to ensure that dynamic load balancing is done => do we want this?
               - bVerbose = FALSE (obvious)
             */
            dd_partition_system(fplog, step, cr, TRUE, 1,
                                &globalState_, top_global, ir,
                                enforcedRotation,
                                localState, f, mdAtoms, top, fr,
                                vsite, constr,
                                nrnb, wcycle, FALSE);
            /* In the current implementation, I go back to the beginning of the main MD loop in md.cpp in case of rejection
             * and the check should be done together with the central compute_globals call for md-vv. So, signalling is enough here.
             */
            *shouldCheckNumberOfBondedInteractions = true;
            update_realloc(upd, localState->natoms);
        }
        else
        {
            *localState = localStateBackUp_;
        }
    }
    totalProposedConfigurations_++;

    // In the first call, the configuration always has to be accepted and used for back-up
#ifdef DEBUG
    if (totalProposedConfigurations_ == 1)
    {
        GMX_ASSERT(acceptedProposedConfigurations_ == 1, "The first back-up of the state has not been made. Logic is broken.");
    }
#endif

    return accepted_;
}

//! \brief Prints the result of the Metropolis step(s) to the log file
void AcceptOrRewind::printOutput(FILE *fplog, const int nstlog, const int64_t step)
{
    // do not let the initial back-up enter statistics
    if (do_per_step(step, nstlog) && fplog && totalProposedConfigurations_ > 1)
    {
        fprintf(fplog,
                "   %19s\n"
                "   %34s %9s\n"
                "   %16s %7.5f\n\n",
                "Hybrid Monte Carlo:",
                "The current configuration has been", (accepted_ ? "accepted." : "rejected."),
                "Empirical P_acc:", (acceptedProposedConfigurations_-1)/(totalProposedConfigurations_-1));
    }
}

} // namespace
