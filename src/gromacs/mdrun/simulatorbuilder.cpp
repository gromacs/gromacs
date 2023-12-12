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
 * \brief Defines the simulator builder for mdrun
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun
 */

#include "gmxpre.h"

#include "simulatorbuilder.h"

#include <memory>

#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/modularsimulator/modularsimulator.h"
#include "gromacs/topology/topology.h"

#include "legacysimulator.h"
#include "membedholder.h"
#include "replicaexchange.h"


namespace gmx
{

//! \brief Build a Simulator object
std::unique_ptr<ISimulator> SimulatorBuilder::build(bool useModularSimulator)
{
    // TODO: Reduce protocol complexity.
    //     Investigate individual parameters. Identify default-constructable parameters and clarify
    //     usage requirements.
    if (!stopHandlerBuilder_)
    {
        throw APIError("You must add a StopHandlerBuilder before calling build().");
    }
    if (!membedHolder_)
    {
        throw APIError("You must add a MembedHolder before calling build().");
    }
    if (!simulatorStateData_)
    {
        throw APIError("Simulator State Data has not been added to the builder");
    }
    if (!simulatorConfig_)
    {
        throw APIError("Simulator config should be set before building the simulator");
    }
    if (!simulatorEnv_)
    {
        throw APIError("You must add a SimulatorEnv before calling build().");
    }
    if (!profiling_)
    {
        throw APIError("You must add a Profiling before calling build().");
    }
    if (!constraintsParam_)
    {
        throw APIError("You must add a ConstraintsParam before calling build().");
    }
    if (!legacyInput_)
    {
        throw APIError("You must add a LegacyInput before calling build().");
    }
    if (!replicaExchangeParameters_)
    {
        throw APIError("You must add a ReplicaExchangeParameters before calling build().");
    }
    if (!interactiveMD_)
    {
        throw APIError("You must add a InteractiveMD before calling build().");
    }
    if (!simulatorModules_)
    {
        throw APIError("You must add a SimulatorModules before calling build().");
    }
    if (!centerOfMassPulling_)
    {
        throw APIError("You must add a CenterOfMassPulling before calling build().");
    }
    if (!ionSwapping_)
    {
        throw APIError("You must add a IonSwapping before calling build().");
    }
    if (!topologyData_)
    {
        throw APIError("You must add a TopologyData before calling build().");
    }

    if (useModularSimulator)
    {
        // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
        return std::unique_ptr<ModularSimulator>(new ModularSimulator(
                std::make_unique<LegacySimulatorData>(simulatorEnv_->fplog_,
                                                      simulatorEnv_->commRec_,
                                                      simulatorEnv_->multisimCommRec_,
                                                      simulatorEnv_->logger_,
                                                      legacyInput_->numFile,
                                                      legacyInput_->filenames,
                                                      simulatorEnv_->outputEnv_,
                                                      simulatorConfig_->mdrunOptions_,
                                                      simulatorConfig_->startingBehavior_,
                                                      constraintsParam_->vsite,
                                                      constraintsParam_->constr,
                                                      constraintsParam_->enforcedRotation_,
                                                      boxDeformation_->deform,
                                                      simulatorModules_->outputProvider,
                                                      simulatorModules_->mdModulesNotifiers,
                                                      legacyInput_->inputrec,
                                                      interactiveMD_->imdSession_,
                                                      centerOfMassPulling_->pull_work,
                                                      ionSwapping_->ionSwap_,
                                                      topologyData_->globalTopology_,
                                                      topologyData_->localTopology_,
                                                      simulatorStateData_->globalState_p,
                                                      simulatorStateData_->localState_p,
                                                      simulatorStateData_->observablesHistory_p,
                                                      topologyData_->mdAtoms_,
                                                      profiling_->nrnb_,
                                                      profiling_->wallCycle_,
                                                      legacyInput_->forceRec_,
                                                      simulatorStateData_->enerdata_p,
                                                      simulatorEnv_->observablesReducerBuilder_,
                                                      simulatorStateData_->ekindata_p,
                                                      simulatorConfig_->runScheduleWork_,
                                                      *replicaExchangeParameters_,
                                                      membedHolder_->membed(),
                                                      profiling_->wallTimeAccounting_,
                                                      std::move(stopHandlerBuilder_),
                                                      simulatorConfig_->mdrunOptions_.rerun),
                std::move(modularSimulatorCheckpointData_)));
    }
    // NOLINTNEXTLINE(modernize-make-unique): make_unique does not work with private constructor
    return std::unique_ptr<LegacySimulator>(new LegacySimulator(simulatorEnv_->fplog_,
                                                                simulatorEnv_->commRec_,
                                                                simulatorEnv_->multisimCommRec_,
                                                                simulatorEnv_->logger_,
                                                                legacyInput_->numFile,
                                                                legacyInput_->filenames,
                                                                simulatorEnv_->outputEnv_,
                                                                simulatorConfig_->mdrunOptions_,
                                                                simulatorConfig_->startingBehavior_,
                                                                constraintsParam_->vsite,
                                                                constraintsParam_->constr,
                                                                constraintsParam_->enforcedRotation_,
                                                                boxDeformation_->deform,
                                                                simulatorModules_->outputProvider,
                                                                simulatorModules_->mdModulesNotifiers,
                                                                legacyInput_->inputrec,
                                                                interactiveMD_->imdSession_,
                                                                centerOfMassPulling_->pull_work,
                                                                ionSwapping_->ionSwap_,
                                                                topologyData_->globalTopology_,
                                                                topologyData_->localTopology_,
                                                                simulatorStateData_->globalState_p,
                                                                simulatorStateData_->localState_p,
                                                                simulatorStateData_->observablesHistory_p,
                                                                topologyData_->mdAtoms_,
                                                                profiling_->nrnb_,
                                                                profiling_->wallCycle_,
                                                                legacyInput_->forceRec_,
                                                                simulatorStateData_->enerdata_p,
                                                                simulatorEnv_->observablesReducerBuilder_,
                                                                simulatorStateData_->ekindata_p,
                                                                simulatorConfig_->runScheduleWork_,
                                                                *replicaExchangeParameters_,
                                                                membedHolder_->membed(),
                                                                profiling_->wallTimeAccounting_,
                                                                std::move(stopHandlerBuilder_),
                                                                simulatorConfig_->mdrunOptions_.rerun));
}

void SimulatorBuilder::add(MembedHolder&& membedHolder)
{
    membedHolder_ = std::make_unique<MembedHolder>(std::move(membedHolder));
}

void SimulatorBuilder::add(ReplicaExchangeParameters&& replicaExchangeParameters)
{
    replicaExchangeParameters_ = std::make_unique<ReplicaExchangeParameters>(replicaExchangeParameters);
}

void SimulatorBuilder::add(std::unique_ptr<ReadCheckpointDataHolder> modularSimulatorCheckpointData)
{
    modularSimulatorCheckpointData_ = std::move(modularSimulatorCheckpointData);
}


} // namespace gmx
