/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief Defines the domain decomposition helper for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "domdechelper.h"

#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"

#include "statepropagatordata.h"
#include "topologyholder.h"

namespace gmx
{
DomDecHelper::DomDecHelper(bool                               isVerbose,
                           int                                verbosePrintInterval,
                           StatePropagatorData*               statePropagatorData,
                           TopologyHolder*                    topologyHolder,
                           CheckBondedInteractionsCallbackPtr checkBondedInteractionsCallback,
                           int                                nstglobalcomm,
                           FILE*                              fplog,
                           t_commrec*                         cr,
                           const MDLogger&                    mdlog,
                           Constraints*                       constr,
                           t_inputrec*                        inputrec,
                           MDAtoms*                           mdAtoms,
                           t_nrnb*                            nrnb,
                           gmx_wallcycle*                     wcycle,
                           t_forcerec*                        fr,
                           gmx_vsite_t*                       vsite,
                           ImdSession*                        imdSession,
                           pull_t*                            pull_work) :
    nextNSStep_(-1),
    isVerbose_(isVerbose),
    verbosePrintInterval_(verbosePrintInterval),
    nstglobalcomm_(nstglobalcomm),
    statePropagatorData_(statePropagatorData),
    topologyHolder_(topologyHolder),
    checkBondedInteractionsCallback_(std::move(checkBondedInteractionsCallback)),
    fplog_(fplog),
    cr_(cr),
    mdlog_(mdlog),
    constr_(constr),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms),
    nrnb_(nrnb),
    wcycle_(wcycle),
    fr_(fr),
    vsite_(vsite),
    imdSession_(imdSession),
    pull_work_(pull_work)
{
    GMX_ASSERT(DOMAINDECOMP(cr), "Domain decomposition Helper constructed in non-DD simulation");
}

void DomDecHelper::setup()
{
    std::unique_ptr<t_state> localState   = statePropagatorData_->localState();
    t_state*                 globalState  = statePropagatorData_->globalState();
    PaddedHostVector<RVec>*  forcePointer = statePropagatorData_->forcePointer();

    // constant choices for this call to dd_partition_system
    const bool     verbose       = false;
    const bool     isMasterState = true;
    const int      nstglobalcomm = 1;
    gmx_wallcycle* wcycle        = nullptr;

    // Distribute the charge groups over the nodes from the master node
    dd_partition_system(fplog_, mdlog_, inputrec_->init_step, cr_, isMasterState, nstglobalcomm,
                        globalState, topologyHolder_->globalTopology(), inputrec_, imdSession_,
                        pull_work_, localState.get(), forcePointer, mdAtoms_,
                        topologyHolder_->localTopology_.get(), fr_, vsite_, constr_, nrnb_, wcycle,
                        verbose);
    topologyHolder_->updateLocalTopology();
    (*checkBondedInteractionsCallback_)();
    statePropagatorData_->setLocalState(std::move(localState));
}

void DomDecHelper::run(Step step, Time gmx_unused time)
{
    if (step != nextNSStep_ || (step == inputrec_->init_step && inputrec_->bContinuation))
    {
        return;
    }
    std::unique_ptr<t_state> localState   = statePropagatorData_->localState();
    t_state*                 globalState  = statePropagatorData_->globalState();
    PaddedHostVector<RVec>*  forcePointer = statePropagatorData_->forcePointer();

    // constant choices for this call to dd_partition_system
    const bool verbose = isVerbose_ && (step % verbosePrintInterval_ == 0 || step == inputrec_->init_step);
    bool isMasterState = false;

    // Correct the new box if it is too skewed
    if (inputrecDynamicBox(inputrec_))
    {
        // TODO: Correcting the box is done here (if using DD) or in ForceElement (non-DD simulations).
        //       Think about unifying this responsibility, could this be done in one place?
        t_graph* graph = nullptr;
        if (correct_box(fplog_, step, localState->box, graph))
        {
            isMasterState = true;
        }
    }
    if (isMasterState)
    {
        dd_collect_state(cr_->dd, localState.get(), globalState);
    }

    // Distribute the charge groups over the nodes from the master node
    dd_partition_system(fplog_, mdlog_, step, cr_, isMasterState, nstglobalcomm_, globalState,
                        topologyHolder_->globalTopology(), inputrec_, imdSession_, pull_work_,
                        localState.get(), forcePointer, mdAtoms_, topologyHolder_->localTopology_.get(),
                        fr_, vsite_, constr_, nrnb_, wcycle_, verbose);
    topologyHolder_->updateLocalTopology();
    (*checkBondedInteractionsCallback_)();
    statePropagatorData_->setLocalState(std::move(localState));
}

SignallerCallbackPtr DomDecHelper::registerNSCallback()
{
    return std::make_unique<SignallerCallback>(
            [this](Step step, Time gmx_unused time) { this->nextNSStep_ = step; });
}

} // namespace gmx
