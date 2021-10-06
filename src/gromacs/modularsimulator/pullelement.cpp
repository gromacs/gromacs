/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * \brief Defines the pull element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
#include "gromacs/pulling/pull.h"

#include "pullelement.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"

namespace gmx
{

PullElement::PullElement(bool                 setPbcRefToPrevStepCOM,
                         PbcType              pbcType,
                         StatePropagatorData* statePropagatorData,
                         pull_t*              pullWork,
                         const t_commrec*     commrec,
                         const MDAtoms*       mdAtoms) :
    setPbcRefToPrevStepCOM_(setPbcRefToPrevStepCOM),
    pbcType_(pbcType),
    restoredFromCheckpoint_(false),
    statePropagatorData_(statePropagatorData),
    pullWork_(pullWork),
    commrec_(commrec),
    mdAtoms_(mdAtoms)
{
}

void PullElement::elementSetup()
{
    if (setPbcRefToPrevStepCOM_ && !restoredFromCheckpoint_)
    {
        preparePrevStepPullComNewSimulation(
                commrec_,
                pullWork_,
                arrayRefFromArray(mdAtoms_->mdatoms()->massT, mdAtoms_->mdatoms()->nr),
                statePropagatorData_->constPositionsView().unpaddedArrayRef(),
                statePropagatorData_->constBox(),
                pbcType_,
                std::nullopt);
    }
}

void PullElement::scheduleTask(Step /*unused*/, Time /*unused*/, const RegisterRunFunction& registerRunFunction)
{
    if (setPbcRefToPrevStepCOM_)
    {
        registerRunFunction([this]() { updatePrevStepPullCom(pullWork_, std::nullopt); });
    }
}

void PullElement::schedulePostStep(Step step, Time time, const RegisterRunFunction& registerRunFunction)
{
    // Printing output must happen after all external pull potentials
    // (currently only AWH) were applied, so execute this after step
    if (MASTER(commrec_))
    {
        registerRunFunction([this, step, time]() { pull_print_output(pullWork_, step, time); });
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
static void doCheckpointData(CheckpointData<operation>* checkpointData, ArrayRef<double> previousStepCom)
{
    checkpointVersion(checkpointData, "PullElement version", c_currentVersion);
    checkpointData->arrayRef("Previous step COM positions",
                             makeCheckpointArrayRef<operation>(previousStepCom));
}

void PullElement::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr)
{
    if (MASTER(cr))
    {
        auto previousStepCom = prevStepPullCom(pullWork_);
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value(), previousStepCom);
    }
}

void PullElement::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                         const t_commrec*                  cr)
{
    auto previousStepCom = prevStepPullCom(pullWork_);
    if (MASTER(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value(), previousStepCom);
    }
    if (haveDDAtomOrdering(*cr))
    {
        gmx_bcast(sizeof(double) * previousStepCom.size(), previousStepCom.data(), cr->mpi_comm_mygroup);
    }
    setPrevStepPullCom(pullWork_, previousStepCom);
    restoredFromCheckpoint_ = true;
}

const std::string& PullElement::clientID()
{
    return identifier_;
}

ISimulatorElement* PullElement::getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                      ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                      StatePropagatorData* statePropagatorData,
                                                      EnergyData* /*energyData*/,
                                                      FreeEnergyPerturbationData* /*freeEnergyPerturbationData*/,
                                                      GlobalCommunicationHelper* /*globalCommunicationHelper*/,
                                                      ObservablesReducer* /*observablesReducer*/)
{
    auto* pullElement = builderHelper->storeElement(std::make_unique<PullElement>(
            legacySimulatorData->inputrec->pull->bSetPbcRefToPrevStepCOM,
            legacySimulatorData->inputrec->pbcType,
            statePropagatorData,
            legacySimulatorData->pull_work,
            legacySimulatorData->cr,
            legacySimulatorData->mdAtoms));
    // Printing output is scheduled after the step
    builderHelper->registerPostStepScheduling(
            [pullElement](Step step, Time time, const RegisterRunFunction& registerRunFunction) {
                pullElement->schedulePostStep(step, time, registerRunFunction);
            });
    return pullElement;
}

} // namespace gmx
