/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Defines the free energy perturbation element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "freeenergyperturbationdata.h"

#include <algorithm>
#include <functional>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/mdlib/freeenergyparameters.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

#include "modularsimulator.h"
#include "simulatoralgorithm.h"

namespace gmx
{
template<CheckpointDataOperation operation>
class CheckpointData;

FreeEnergyPerturbationData::FreeEnergyPerturbationData(FILE*             fplog,
                                                       const t_inputrec& inputrec,
                                                       MDAtoms*          mdAtoms,
                                                       gmx_ekindata_t*   ekindata) :
    element_(std::make_unique<Element>(this, inputrec.fepvals->delta_lambda)),
    lambda_(),
    currentFEPState_(0),
    fplog_(fplog),
    inputrec_(inputrec),
    mdAtoms_(mdAtoms)
{
    std::fill(lambda_.begin(), lambda_.end(), 0);
    // The legacy implementation only filled the lambda vector in state_global, which is only
    // available on main. We have the lambda vector available everywhere, so we pass a `true`
    // for isMain on all ranks. See #3647.
    initialize_lambdas(fplog_,
                       inputrec_.efep,
                       inputrec_.bSimTemp,
                       *inputrec_.fepvals,
                       inputrec_.simtempvals->temperatures,
                       ekindata,
                       true,
                       &currentFEPState_,
                       lambda_);
}

void FreeEnergyPerturbationData::Element::scheduleTask(Step                       step,
                                                       Time gmx_unused            time,
                                                       const RegisterRunFunction& registerRunFunction)
{
    // If we do slow growth, we update lambda every step
    // If it's set externally, we get notified, so we only update when necessary (at nextLambdaSettingStep_)
    // However, if we reload from checkpoint, it might be that checkpointing happened right between the
    // external caller setting the state and us applying it, so we also check newFepStateStep_.
    const bool needToSetExternalState = externalFepStateSetting_
                                        && ((step == externalFepStateSetting_->nextFepStateSettingStep)
                                            || (step == externalFepStateSetting_->newFepStateStep));
    if (doSlowGrowth_)
    {
        registerRunFunction([this, step]() { freeEnergyPerturbationData_->updateLambdas(step); });
    }
    else if (needToSetExternalState)
    {
        registerRunFunction([this, step]() {
            GMX_ASSERT(step == externalFepStateSetting_->newFepStateStep,
                       "FEP state setting step mismatch");
            freeEnergyPerturbationData_->setLambdaState(step, externalFepStateSetting_->newFepState);
        });
    }
}

void FreeEnergyPerturbationData::updateLambdas(Step step)
{
    // at beginning of step (if lambdas change...)
    lambda_ = currentLambdas(step, *(inputrec_.fepvals), currentFEPState_);
    updateMDAtoms();
}

void FreeEnergyPerturbationData::setLambdaState(Step step, int newState)
{
    currentFEPState_ = newState;
    updateLambdas(step);
}

ArrayRef<real> FreeEnergyPerturbationData::lambdaView()
{
    return lambda_;
}

ArrayRef<const real> FreeEnergyPerturbationData::constLambdaView() const
{
    return lambda_;
}

int FreeEnergyPerturbationData::currentFEPState() const
{
    return currentFEPState_;
}

void FreeEnergyPerturbationData::updateMDAtoms()
{
    update_mdatoms(mdAtoms_->mdatoms(), lambda_[FreeEnergyPerturbationCouplingType::Mass]);
}

FepStateSetting* FreeEnergyPerturbationData::enableExternalFepStateSetting() const
{
    return element_->enableExternalFepStateSetting();
}

namespace
{
/*!
 * \brief Enum describing the contents FreeEnergyPerturbationData::Element writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class CheckpointVersion
{
    Base,                       //!< First version of modular checkpointing
    AddedExternalLambdaSetting, //!< Additional values to ensure no state setting info is lost
    Count                       //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersion = CheckpointVersion(int(CheckpointVersion::Count) - 1);
} // namespace

template<CheckpointDataOperation operation>
void FreeEnergyPerturbationData::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "FreeEnergyPerturbationData version", c_currentVersion);

    checkpointData->scalar("current FEP state", &currentFEPState_);
    checkpointData->arrayRef("lambda vector", makeCheckpointArrayRef<operation>(lambda_));
}

template<CheckpointDataOperation operation>
void FreeEnergyPerturbationData::Element::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    CheckpointVersion fileVersion = c_currentVersion;
    if (operation == CheckpointDataOperation::Read)
    {
        // We can read the same key as above - we can only write it once, though!
        fileVersion = checkpointVersion(
                checkpointData, "FreeEnergyPerturbationData version", c_currentVersion);
    }

    if (fileVersion >= CheckpointVersion::AddedExternalLambdaSetting)
    {
        // If checkpointing happens between receiving the request and actually setting the new
        // lambda state, we need to preserve this information.
        bool externalLambdaSetting = externalFepStateSetting_.has_value();
        checkpointData->scalar("External lambda setting", &externalLambdaSetting);
        if constexpr (operation == CheckpointDataOperation::Read)
        {
            if (externalLambdaSetting)
            {
                externalFepStateSetting_ = FepStateSetting();
            }
        }
        if (externalFepStateSetting_.has_value()) // NOLINT(readability-misleading-indentation)
        {
            checkpointData->scalar("Requested new FEP state", &externalFepStateSetting_->newFepState);
            checkpointData->scalar("Step at which new FEP state is applied",
                                   &externalFepStateSetting_->newFepStateStep);
        }
    }
}

void FreeEnergyPerturbationData::Element::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                                              const t_commrec*                   cr)
{
    if (MAIN(cr))
    {
        freeEnergyPerturbationData_->doCheckpointData<CheckpointDataOperation::Write>(
                &checkpointData.value());
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void FreeEnergyPerturbationData::Element::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                                 const t_commrec* cr)
{
    if (MAIN(cr))
    {
        freeEnergyPerturbationData_->doCheckpointData<CheckpointDataOperation::Read>(
                &checkpointData.value());
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_bcast(cr->dd, sizeof(int), &freeEnergyPerturbationData_->currentFEPState_);
        dd_bcast(cr->dd,
                 ssize(freeEnergyPerturbationData_->lambda_) * int(sizeof(real)),
                 freeEnergyPerturbationData_->lambda_.data());
        auto externalLambdaSetting = int(externalFepStateSetting_.has_value());
        dd_bcast(cr->dd, sizeof(int), &externalLambdaSetting);
        if (!MAIN(cr) && externalLambdaSetting)
        {
            // Main rank constructed this while reading the
            // checkpoint, but other ranks have to do this now.
            externalFepStateSetting_ = FepStateSetting();
        }
        if (externalFepStateSetting_.has_value())
        {
            dd_bcast(cr->dd,
                     sizeof(externalFepStateSetting_->newFepState),
                     &externalFepStateSetting_->newFepState);
            dd_bcast(cr->dd,
                     sizeof(externalFepStateSetting_->newFepStateStep),
                     &externalFepStateSetting_->newFepStateStep);
        }
    }
    restoredFromCheckpoint_ = true;
}

const std::string& FreeEnergyPerturbationData::Element::clientID()
{
    return FreeEnergyPerturbationData::checkpointID();
}

DomDecCallback FreeEnergyPerturbationData::Element::registerDomDecCallback()
{
    return [this]() { freeEnergyPerturbationData_->updateMDAtoms(); };
}

FreeEnergyPerturbationData::Element::Element(FreeEnergyPerturbationData* freeEnergyPerturbationElement,
                                             double                      deltaLambda) :
    freeEnergyPerturbationData_(freeEnergyPerturbationElement),
    doSlowGrowth_(deltaLambda != 0),
    numExternalFepStateSettingClients_(0),
    restoredFromCheckpoint_(false)
{
}

void FreeEnergyPerturbationData::Element::elementSetup()
{
    // Sanity check ensuring that checkpointed data matches the algorithms that were set up
    GMX_RELEASE_ASSERT(!(externalFepStateSetting_.has_value() && numExternalFepStateSettingClients_ == 0),
                       "Checkpoint mismatch: Checkpointed simulation used external lambda setting, "
                       "while the current simulation does not.");
    freeEnergyPerturbationData_->updateMDAtoms();
}

FreeEnergyPerturbationData::Element* FreeEnergyPerturbationData::element()
{
    return element_.get();
}

FepStateSetting* FreeEnergyPerturbationData::Element::enableExternalFepStateSetting()
{
    GMX_RELEASE_ASSERT(!doSlowGrowth_,
                       "External FEP state setting is incompatible with slow growth.");
    // This could be implemented with some sanity checks, but there's no use case right now
    GMX_RELEASE_ASSERT(numExternalFepStateSettingClients_ == 0,
                       "External FEP state setting by more than one element not supported.");
    numExternalFepStateSettingClients_++;
    // externalFepStateSetting_ may already have been constructed from checkpoint
    if (!externalFepStateSetting_.has_value())
    {
        GMX_RELEASE_ASSERT(!restoredFromCheckpoint_,
                           "Checkpoint mismatch: Checkpointed simulation did not use external "
                           "lambda setting, while the current simulation does.");
        externalFepStateSetting_ = FepStateSetting();
    }
    return &externalFepStateSetting_.value();
}

void FepStateSetting::signalSettingStep(Step step)
{
    nextFepStateSettingStep = step;
}

void FepStateSetting::setNewState(int state, Step step)
{
    newFepState     = state;
    newFepStateStep = step;
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

void FreeEnergyPerturbationData::readCheckpointToTrxFrame(t_trxframe* trxFrame,
                                                          std::optional<ReadCheckpointData> readCheckpointData)
{
    if (readCheckpointData)
    {
        FreeEnergyPerturbationData freeEnergyPerturbationData(nullptr, t_inputrec(), nullptr, nullptr);
        freeEnergyPerturbationData.doCheckpointData(&readCheckpointData.value());
        trxFrame->lambda = freeEnergyPerturbationData.lambda_[FreeEnergyPerturbationCouplingType::Fep];
        trxFrame->fep_state = freeEnergyPerturbationData.currentFEPState_;
    }
    else
    {
        trxFrame->lambda    = 0;
        trxFrame->fep_state = 0;
    }
    trxFrame->bLambda = true;
}

const std::string& FreeEnergyPerturbationData::checkpointID()
{
    static const std::string identifier = "FreeEnergyPerturbationData";
    return identifier;
}

} // namespace gmx
