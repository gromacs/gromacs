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
 * \brief Defines the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "modularsimulator.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdrun/replicaexchange.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"

#include "computeglobalselement.h"
#include "constraintelement.h"
#include "forceelement.h"
#include "parrinellorahmanbarostat.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"
#include "vrescalethermostat.h"

namespace gmx
{
void ModularSimulator::run()
{
    GMX_LOG(legacySimulatorData_->mdlog.info)
            .asParagraph()
            .appendText("Using the modular simulator.");

    ModularSimulatorAlgorithmBuilder algorithmBuilder(compat::make_not_null(legacySimulatorData_));
    addIntegrationElements(&algorithmBuilder);
    auto algorithm = algorithmBuilder.build();

    while (const auto* task = algorithm.getNextTask())
    {
        // execute task
        (*task)();
    }
}

void ModularSimulator::addIntegrationElements(ModularSimulatorAlgorithmBuilder* builder)
{
    if (legacySimulatorData_->inputrec->eI == eiMD)
    {
        // The leap frog integration algorithm
        builder->add<ForceElement>();
        builder->add<StatePropagatorData::Element>();
        if (legacySimulatorData_->inputrec->etc == etcVRESCALE)
        {
            builder->add<VRescaleThermostat>(-1, UseFullStepKE::No, ReportPreviousStepConservedEnergy::No);
        }
        builder->add<Propagator<IntegrationStep::LeapFrog>>(legacySimulatorData_->inputrec->delta_t,
                                                            RegisterWithThermostat::True,
                                                            RegisterWithBarostat::True);
        if (legacySimulatorData_->constr)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Positions>>();
        }
        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>>();
        builder->add<EnergyData::Element>();
        if (legacySimulatorData_->inputrec->epc == epcPARRINELLORAHMAN)
        {
            builder->add<ParrinelloRahmanBarostat>(-1);
        }
    }
    else if (legacySimulatorData_->inputrec->eI == eiVV)
    {
        // The velocity verlet integration algorithm
        builder->add<ForceElement>();
        builder->add<Propagator<IntegrationStep::VelocitiesOnly>>(
                0.5 * legacySimulatorData_->inputrec->delta_t, RegisterWithThermostat::False,
                RegisterWithBarostat::True);
        if (legacySimulatorData_->constr)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Velocities>>();
        }
        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();
        builder->add<StatePropagatorData::Element>();
        if (legacySimulatorData_->inputrec->etc == etcVRESCALE)
        {
            builder->add<VRescaleThermostat>(0, UseFullStepKE::Yes, ReportPreviousStepConservedEnergy::Yes);
        }
        builder->add<Propagator<IntegrationStep::VelocityVerletPositionsAndVelocities>>(
                legacySimulatorData_->inputrec->delta_t, RegisterWithThermostat::True,
                RegisterWithBarostat::False);
        if (legacySimulatorData_->constr)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Positions>>();
        }
        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();
        builder->add<EnergyData::Element>();
        if (legacySimulatorData_->inputrec->epc == epcPARRINELLORAHMAN)
        {
            builder->add<ParrinelloRahmanBarostat>(-1);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Integrator not implemented for the modular simulator.");
    }
}

bool ModularSimulator::isInputCompatible(bool                             exitOnFailure,
                                         const t_inputrec*                inputrec,
                                         bool                             doRerun,
                                         const gmx_mtop_t&                globalTopology,
                                         const gmx_multisim_t*            ms,
                                         const ReplicaExchangeParameters& replExParams,
                                         const t_fcdata*                  fcd,
                                         bool                             doEssentialDynamics,
                                         bool                             doMembed)
{
    auto conditionalAssert = [exitOnFailure](bool condition, const char* message) {
        if (exitOnFailure)
        {
            GMX_RELEASE_ASSERT(condition, message);
        }
        return condition;
    };

    bool isInputCompatible = true;

    // GMX_USE_MODULAR_SIMULATOR allows to use modular simulator also for non-standard uses,
    // such as the leap-frog integrator
    const auto modularSimulatorExplicitlyTurnedOn = (getenv("GMX_USE_MODULAR_SIMULATOR") != nullptr);
    // GMX_USE_MODULAR_SIMULATOR allows to use disable modular simulator for all uses,
    // including the velocity-verlet integrator used by default
    const auto modularSimulatorExplicitlyTurnedOff = (getenv("GMX_DISABLE_MODULAR_SIMULATOR") != nullptr);

    GMX_RELEASE_ASSERT(
            !(modularSimulatorExplicitlyTurnedOn && modularSimulatorExplicitlyTurnedOff),
            "Cannot have both GMX_USE_MODULAR_SIMULATOR=ON and GMX_DISABLE_MODULAR_SIMULATOR=ON. "
            "Unset one of the two environment variables to explicitly chose which simulator to "
            "use, "
            "or unset both to recover default behavior.");

    GMX_RELEASE_ASSERT(
            !(modularSimulatorExplicitlyTurnedOff && inputrec->eI == eiVV
              && inputrec->epc == epcPARRINELLORAHMAN),
            "Cannot use a Parrinello-Rahman barostat with md-vv and "
            "GMX_DISABLE_MODULAR_SIMULATOR=ON, "
            "as the Parrinello-Rahman barostat is not implemented in the legacy simulator. Unset "
            "GMX_DISABLE_MODULAR_SIMULATOR or use a different pressure control algorithm.");

    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       inputrec->eI == eiMD || inputrec->eI == eiVV,
                       "Only integrators md and md-vv are supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->eI != eiMD || modularSimulatorExplicitlyTurnedOn,
                                             "Set GMX_USE_MODULAR_SIMULATOR=ON to use the modular "
                                             "simulator with integrator md.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doRerun, "Rerun is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       inputrec->etc == etcNO || inputrec->etc == etcVRESCALE,
                       "Only v-rescale thermostat is supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       inputrec->epc == epcNO || inputrec->epc == epcPARRINELLORAHMAN,
                       "Only Parrinello-Rahman barostat is supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       !(inputrecNptTrotter(inputrec) || inputrecNphTrotter(inputrec)
                         || inputrecNvtTrotter(inputrec)),
                       "Legacy Trotter decomposition is not supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->efep == efepNO || inputrec->efep == efepYES
                                                     || inputrec->efep == efepSLOWGROWTH,
                                             "Expanded ensemble free energy calculation is not "
                                             "supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(!inputrec->bPull,
                                             "Pulling is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(inputrec->opts.ngacc == 1 && inputrec->opts.acc[0][XX] == 0.0
                                         && inputrec->opts.acc[0][YY] == 0.0
                                         && inputrec->opts.acc[0][ZZ] == 0.0 && inputrec->cos_accel == 0.0,
                                 "Acceleration is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(inputrec->opts.ngfrz == 1 && inputrec->opts.nFreeze[0][XX] == 0
                                         && inputrec->opts.nFreeze[0][YY] == 0
                                         && inputrec->opts.nFreeze[0][ZZ] == 0,
                                 "Freeze groups are not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       inputrec->deform[XX][XX] == 0.0 && inputrec->deform[XX][YY] == 0.0
                               && inputrec->deform[XX][ZZ] == 0.0 && inputrec->deform[YY][XX] == 0.0
                               && inputrec->deform[YY][YY] == 0.0 && inputrec->deform[YY][ZZ] == 0.0
                               && inputrec->deform[ZZ][XX] == 0.0 && inputrec->deform[ZZ][YY] == 0.0
                               && inputrec->deform[ZZ][ZZ] == 0.0,
                       "Deformation is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(gmx_mtop_interaction_count(globalTopology, IF_VSITE) == 0,
                                 "Virtual sites are not supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(!inputrec->bDoAwh,
                                             "AWH is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(gmx_mtop_ftype_count(globalTopology, F_DISRES) == 0,
                                 "Distance restraints are not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       gmx_mtop_ftype_count(globalTopology, F_ORIRES) == 0,
                       "Orientation restraints are not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(ms == nullptr,
                                 "Multi-sim are not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(replExParams.exchangeInterval == 0,
                                 "Replica exchange is not supported by the modular simulator.");

    int numEnsembleRestraintSystems;
    if (fcd)
    {
        numEnsembleRestraintSystems = fcd->disres->nsystems;
    }
    else
    {
        auto distantRestraintEnsembleEnvVar = getenv("GMX_DISRE_ENSEMBLE_SIZE");
        numEnsembleRestraintSystems =
                (ms != nullptr && distantRestraintEnsembleEnvVar != nullptr)
                        ? static_cast<int>(strtol(distantRestraintEnsembleEnvVar, nullptr, 10))
                        : 0;
    }
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(numEnsembleRestraintSystems <= 1,
                                 "Ensemble restraints are not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doSimulatedAnnealing(inputrec),
                                 "Simulated annealing is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!inputrec->bSimTemp,
                                 "Simulated tempering is not supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(!inputrec->bExpanded,
                                             "Expanded ensemble simulations are not supported by "
                                             "the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doEssentialDynamics,
                                 "Essential dynamics is not supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->eSwapCoords == eswapNO,
                                             "Ion / water position swapping is not supported by "
                                             "the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!inputrec->bIMD,
                                 "Interactive MD is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doMembed,
                                 "Membrane embedding is not supported by the modular simulator.");
    // TODO: Change this to the boolean passed when we merge the user interface change for the GPU update.
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                       getenv("GMX_FORCE_UPDATE_DEFAULT_GPU") == nullptr,
                       "Integration on the GPU is not supported by the modular simulator.");
    // Modular simulator is centered around NS updates
    // TODO: think how to handle nstlist == 0
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->nstlist != 0,
                                             "Simulations without neighbor list update are not "
                                             "supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(!GMX_FAHCORE,
                                             "GMX_FAHCORE not supported by the modular simulator.");

    return isInputCompatible;
}

ModularSimulator::ModularSimulator(std::unique_ptr<LegacySimulatorData>      legacySimulatorData,
                                   std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder) :
    legacySimulatorData_(std::move(legacySimulatorData)),
    checkpointDataHolder_(std::move(checkpointDataHolder))
{
    checkInputForDisabledFunctionality();
}

void ModularSimulator::checkInputForDisabledFunctionality()
{
    isInputCompatible(true, legacySimulatorData_->inputrec, legacySimulatorData_->mdrunOptions.rerun,
                      *legacySimulatorData_->top_global, legacySimulatorData_->ms,
                      legacySimulatorData_->replExParams, legacySimulatorData_->fr->fcdata.get(),
                      opt2bSet("-ei", legacySimulatorData_->nfile, legacySimulatorData_->fnm),
                      legacySimulatorData_->membed != nullptr);
    if (legacySimulatorData_->observablesHistory->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }
}

} // namespace gmx
