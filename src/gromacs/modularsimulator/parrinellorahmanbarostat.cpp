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
 * \brief Defines the Parrinello-Rahman barostat for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "parrinellorahmanbarostat.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"

#include "energyelement.h"
#include "statepropagatordata.h"

namespace gmx
{

ParrinelloRahmanBarostat::ParrinelloRahmanBarostat(int                   nstpcouple,
                                                   int                   offset,
                                                   real                  couplingTimeStep,
                                                   Step                  initStep,
                                                   ArrayRef<rvec>        scalingTensor,
                                                   PropagatorCallbackPtr propagatorCallback,
                                                   StatePropagatorData*  statePropagatorData,
                                                   EnergyElement*        energyElement,
                                                   FILE*                 fplog,
                                                   const t_inputrec*     inputrec,
                                                   const MDAtoms*        mdAtoms,
                                                   const t_state*        globalState,
                                                   t_commrec*            cr,
                                                   bool                  isRestart) :
    nstpcouple_(nstpcouple),
    offset_(offset),
    couplingTimeStep_(couplingTimeStep),
    initStep_(initStep),
    scalingTensor_(scalingTensor),
    propagatorCallback_(std::move(propagatorCallback)),
    statePropagatorData_(statePropagatorData),
    energyElement_(energyElement),
    fplog_(fplog),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms)
{
    clear_mat(mu_);
    clear_mat(boxRel_);
    clear_mat(boxVelocity_);

    // TODO: This is only needed to restore the thermostatIntegral_ from cpt. Remove this when
    //       switching to purely client-based checkpointing.
    if (isRestart)
    {
        if (MASTER(cr))
        {
            copy_mat(globalState->boxv, boxVelocity_);
            copy_mat(globalState->box_rel, boxRel_);
        }
        if (DOMAINDECOMP(cr))
        {
            dd_bcast(cr->dd, sizeof(boxVelocity_), boxVelocity_);
            dd_bcast(cr->dd, sizeof(boxRel_), boxRel_);
        }
    }
}

void ParrinelloRahmanBarostat::scheduleTask(gmx::Step step,
                                            gmx::Time gmx_unused               time,
                                            const gmx::RegisterRunFunctionPtr& registerRunFunction)
{
    const bool scaleOnNextStep = do_per_step(step + nstpcouple_ + offset_ + 1, nstpcouple_);
    const bool scaleOnThisStep = do_per_step(step + nstpcouple_ + offset_, nstpcouple_);

    if (scaleOnThisStep)
    {
        (*registerRunFunction)(
                std::make_unique<SimulatorRunFunction>([this]() { scaleBoxAndPositions(); }));
    }
    if (scaleOnNextStep)
    {
        (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                [this, step]() { integrateBoxVelocityEquations(step); }));
        // let propagator know that it will have to scale on next step
        (*propagatorCallback_)(step + 1);
    }
}

void ParrinelloRahmanBarostat::integrateBoxVelocityEquations(Step step)
{
    auto box = statePropagatorData_->constBox();
    parrinellorahman_pcoupl(fplog_, step, inputrec_, couplingTimeStep_, energyElement_->pressure(step),
                            box, boxRel_, boxVelocity_, scalingTensor_.data(), mu_, false);
    // multiply matrix by the coupling time step to avoid having the propagator needing to know about that
    msmul(scalingTensor_.data(), couplingTimeStep_, scalingTensor_.data());
}

void ParrinelloRahmanBarostat::scaleBoxAndPositions()
{
    // Propagate the box by the box velocities
    auto box = statePropagatorData_->box();
    for (int i = 0; i < DIM; i++)
    {
        for (int m = 0; m <= i; m++)
        {
            box[i][m] += couplingTimeStep_ * boxVelocity_[i][m];
        }
    }
    preserve_box_shape(inputrec_, boxRel_, box);

    // Scale the coordinates
    const int start  = 0;
    const int homenr = mdAtoms_->mdatoms()->homenr;
    auto      x      = as_rvec_array(statePropagatorData_->positionsView().paddedArrayRef().data());
    for (int n = start; n < start + homenr; n++)
    {
        tmvmul_ur0(mu_, x[n], x[n]);
    }
}

void ParrinelloRahmanBarostat::elementSetup()
{
    if (inputrecPreserveShape(inputrec_))
    {
        auto      box  = statePropagatorData_->box();
        const int ndim = inputrec_->epct == epctSEMIISOTROPIC ? 2 : 3;
        do_box_rel(ndim, inputrec_->deform, boxRel_, box, true);
    }

    const bool scaleOnInitStep = do_per_step(initStep_ + nstpcouple_ + offset_, nstpcouple_);
    if (scaleOnInitStep)
    {
        // If we need to scale on the first step, we need to set the scaling matrix using the current
        // box velocity. If this is a fresh start, we will hence not move the box (this does currently
        // never happen as the offset is set to -1 in all cases). If this is a restart, we will use
        // the saved box velocity which we would have updated right before checkpointing.
        // Setting bFirstStep = true in parrinellorahman_pcoupl (last argument) makes sure that only
        // the scaling matrix is calculated, without updating the box velocities.
        // The call to parrinellorahman_pcoupl is using nullptr for fplog (since we don't expect any
        // output here) and for the pressure (since it might not be calculated yet, and we don't need it).
        auto box = statePropagatorData_->constBox();
        parrinellorahman_pcoupl(nullptr, initStep_, inputrec_, couplingTimeStep_, nullptr, box,
                                boxRel_, boxVelocity_, scalingTensor_.data(), mu_, true);
        // multiply matrix by the coupling time step to avoid having the propagator needing to know about that
        msmul(scalingTensor_.data(), couplingTimeStep_, scalingTensor_.data());

        (*propagatorCallback_)(initStep_);
    }
}

const rvec* ParrinelloRahmanBarostat::boxVelocities() const
{
    return boxVelocity_;
}

void ParrinelloRahmanBarostat::writeCheckpoint(t_state* localState, t_state gmx_unused* globalState)
{
    copy_mat(boxVelocity_, localState->boxv);
    copy_mat(boxRel_, localState->box_rel);
    localState->flags |= (1U << estBOXV) | (1U << estBOX_REL);
}


} // namespace gmx
