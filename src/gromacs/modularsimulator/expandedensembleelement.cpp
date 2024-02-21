/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief Defines the expanded ensemble element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "expandedensembleelement.h"

#include "gromacs/domdec/distribute.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"

#include "energydata.h"
#include "simulatoralgorithm.h"

namespace gmx
{

void ExpandedEnsembleElement::apply(Step step, bool doLambdaStep, bool doLog)
{
    if (doLambdaStep)
    {
        const int newFepState =
                expandedEnsembleUpdateLambdaState(fplog_,
                                                  inputrec_,
                                                  energyData_->enerdata(),
                                                  freeEnergyPerturbationData_->currentFEPState(),
                                                  dfhist_.get(),
                                                  step);
        // Set new state at next step
        fepStateSetting_->setNewState(newFepState, step + 1);
    }
    if (doLog)
    {
        /* only needed if doing expanded ensemble */
        PrintFreeEnergyInfoToFile(fplog_,
                                  inputrec_->fepvals.get(),
                                  inputrec_->expandedvals.get(),
                                  inputrec_->bSimTemp ? inputrec_->simtempvals.get() : nullptr,
                                  dfhist_.get(),
                                  freeEnergyPerturbationData_->currentFEPState(),
                                  inputrec_->nstlog,
                                  step);
    }
}

void ExpandedEnsembleElement::elementSetup()
{
    // Check nstexpanded here, because the grompp check was broken (#2714)
    if (inputrec_->expandedvals->nstexpanded % inputrec_->nstcalcenergy != 0)
    {
        gmx_fatal(FARGS,
                  "With expanded ensemble, nstexpanded should be a multiple of nstcalcenergy");
    }
    init_expanded_ensemble(restoredFromCheckpoint_, inputrec_, dfhist_.get());
}

void ExpandedEnsembleElement::scheduleTask(Step step, Time /*unused*/, const RegisterRunFunction& registerRunFunction)
{
    const bool isFirstStep  = (step == initialStep_ && !restoredFromCheckpoint_);
    const bool doLambdaStep = (do_per_step(step, frequency_) && !isFirstStep);
    const bool doLog        = (isMainRank_ && step == nextLogWritingStep_ && (fplog_ != nullptr));

    if (doLambdaStep || doLog)
    {
        registerRunFunction([this, step, doLambdaStep, doLog]() { apply(step, doLambdaStep, doLog); });
    }
    if (doLambdaStep)
    {
        // We'll compute a new lambda state and want it applied for next step
        fepStateSetting_->signalSettingStep(step + 1);
    }
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
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_currentVersion = CheckpointVersion(int(CheckpointVersion::Count) - 1);
} // namespace

template<CheckpointDataOperation operation>
void ExpandedEnsembleElement::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "ExpandedEnsembleElement version", c_currentVersion);

    dfhist_->doCheckpoint<operation>(checkpointData->subCheckpointData("dfhist"),
                                     inputrec_->expandedvals->elamstats);
}

void ExpandedEnsembleElement::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                                  const t_commrec*                   cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void ExpandedEnsembleElement::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                     const t_commrec*                  cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_distribute_dfhist(cr->dd, dfhist_.get());
    }
    restoredFromCheckpoint_ = true;
}

const std::string& ExpandedEnsembleElement::clientID()
{
    return identifier_;
}

std::optional<SignallerCallback> ExpandedEnsembleElement::registerLoggingCallback()
{
    if (isMainRank_)
    {
        return [this](Step step, Time /*unused*/) { nextLogWritingStep_ = step; };
    }
    else
    {
        return std::nullopt;
    }
}

ExpandedEnsembleElement::ExpandedEnsembleElement(bool                              isMainRank,
                                                 Step                              initialStep,
                                                 int                               frequency,
                                                 const EnergyData*                 energyData,
                                                 const FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                 FILE*                             fplog,
                                                 const t_inputrec*                 inputrec) :
    fepStateSetting_(freeEnergyPerturbationData->enableExternalFepStateSetting()),
    isMainRank_(isMainRank),
    initialStep_(initialStep),
    frequency_(frequency),
    nextLogWritingStep_(-1),
    dfhist_(std::make_unique<df_history_t>()),
    restoredFromCheckpoint_(false),
    energyData_(energyData),
    freeEnergyPerturbationData_(freeEnergyPerturbationData),
    fplog_(fplog),
    inputrec_(inputrec)
{
    init_df_history(dfhist_.get(), inputrec_->fepvals->n_lambda);
}

ExpandedEnsembleElement::~ExpandedEnsembleElement() = default;

ISimulatorElement* ExpandedEnsembleElement::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData*     freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer* /*observablesReducer*/)
{
    return builderHelper->storeElement(std::make_unique<ExpandedEnsembleElement>(
            MAIN(legacySimulatorData->cr_),
            legacySimulatorData->inputRec_->init_step,
            legacySimulatorData->inputRec_->expandedvals->nstexpanded,
            energyData,
            freeEnergyPerturbationData,
            legacySimulatorData->fpLog_,
            legacySimulatorData->inputRec_));
}

} // namespace gmx
