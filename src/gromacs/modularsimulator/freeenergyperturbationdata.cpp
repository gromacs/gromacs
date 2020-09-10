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
 * \brief Defines the free energy perturbation element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "freeenergyperturbationdata.h"

#include "gromacs/mdlib/freeenergyparameters.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"

#include "modularsimulator.h"
#include "simulatoralgorithm.h"

namespace gmx
{

FreeEnergyPerturbationData::FreeEnergyPerturbationData(FILE* fplog, const t_inputrec* inputrec, MDAtoms* mdAtoms) :
    element_(std::make_unique<Element>(this, inputrec->fepvals->delta_lambda)),
    lambda_(),
    currentFEPState_(0),
    fplog_(fplog),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms)
{
    lambda_.fill(0);
    initialize_lambdas(fplog_, *inputrec_, true, &currentFEPState_, lambda_);
    update_mdatoms(mdAtoms_->mdatoms(), lambda_[efptMASS]);
}

void FreeEnergyPerturbationData::Element::scheduleTask(Step step,
                                                       Time gmx_unused            time,
                                                       const RegisterRunFunction& registerRunFunction)
{
    if (lambdasChange_)
    {
        registerRunFunction([this, step]() { freeEnergyPerturbationData_->updateLambdas(step); });
    }
}

void FreeEnergyPerturbationData::updateLambdas(Step step)
{
    // at beginning of step (if lambdas change...)
    lambda_ = currentLambdas(step, *(inputrec_->fepvals), currentFEPState_);
    update_mdatoms(mdAtoms_->mdatoms(), lambda_[efptMASS]);
}

ArrayRef<real> FreeEnergyPerturbationData::lambdaView()
{
    return lambda_;
}

ArrayRef<const real> FreeEnergyPerturbationData::constLambdaView()
{
    return lambda_;
}

int FreeEnergyPerturbationData::currentFEPState()
{
    return currentFEPState_;
}

void FreeEnergyPerturbationData::Element::writeCheckpoint(t_state* localState, t_state gmx_unused* globalState)
{
    localState->fep_state = freeEnergyPerturbationData_->currentFEPState_;
    localState->lambda    = freeEnergyPerturbationData_->lambda_;
    localState->flags |= (1U << estLAMBDA) | (1U << estFEPSTATE);
}

FreeEnergyPerturbationData::Element::Element(FreeEnergyPerturbationData* freeEnergyPerturbationElement,
                                             double                      deltaLambda) :
    freeEnergyPerturbationData_(freeEnergyPerturbationElement),
    lambdasChange_(deltaLambda != 0)
{
}

FreeEnergyPerturbationData::Element* FreeEnergyPerturbationData::element()
{
    return element_.get();
}

ISimulatorElement* FreeEnergyPerturbationData::Element::getElementPointerImpl(
        LegacySimulatorData gmx_unused*        legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper gmx_unused* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData gmx_unused*      energyData,
        FreeEnergyPerturbationData* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper)
{
    return freeEnergyPerturbationData->element();
}

} // namespace gmx
