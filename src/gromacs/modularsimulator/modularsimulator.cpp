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
#include "gromacs/gmxlib/network.h"
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

#include "compositesimulatorelement.h"
#include "computeglobalselement.h"
#include "constraintelement.h"
#include "energyelement.h"
#include "forceelement.h"
#include "freeenergyperturbationelement.h"
#include "parrinellorahmanbarostat.h"
#include "propagator.h"
#include "signallers.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"
#include "trajectoryelement.h"
#include "vrescalethermostat.h"

namespace gmx
{
void ModularSimulator::run()
{
    GMX_LOG(mdlog.info).asParagraph().appendText("Using the modular simulator.");

    ModularSimulatorAlgorithmBuilder algorithmBuilder(compat::make_not_null(this));
    auto                             algorithm = algorithmBuilder.build();

    while (const auto* task = algorithm.getNextTask())
    {
        // execute task
        (*task)();
    }
}

std::unique_ptr<ISimulatorElement> ModularSimulatorAlgorithmBuilder::buildForces(
        SignallerBuilder<NeighborSearchSignaller>* neighborSearchSignallerBuilder,
        SignallerBuilder<EnergySignaller>*         energySignallerBuilder,
        StatePropagatorData*                       statePropagatorDataPtr,
        EnergyElement*                             energyElementPtr,
        FreeEnergyPerturbationElement*             freeEnergyPerturbationElement,
        TopologyHolder*                            topologyHolder)
{
    const bool isVerbose    = mdrunOptions.verbose;
    const bool isDynamicBox = inputrecDynamicBox(inputrec);

    auto forceElement = std::make_unique<ForceElement>(
            statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElement, isVerbose,
            isDynamicBox, fplog, cr, inputrec, mdAtoms, nrnb, fr, wcycle, runScheduleWork, vsite,
            imdSession, pull_work, constr, top_global, enforcedRotation);
    topologyHolder->registerClient(forceElement.get());
    neighborSearchSignallerBuilder->registerSignallerClient(compat::make_not_null(forceElement.get()));
    energySignallerBuilder->registerSignallerClient(compat::make_not_null(forceElement.get()));

    // std::move *should* not be needed with c++-14, but clang-3.6 still requires it
    return std::move(forceElement);
}

std::unique_ptr<ISimulatorElement> ModularSimulatorAlgorithmBuilder::buildIntegrator(
        SignallerBuilder<NeighborSearchSignaller>* neighborSearchSignallerBuilder,
        SignallerBuilder<EnergySignaller>*         energySignallerBuilder,
        SignallerBuilder<LoggingSignaller>*        loggingSignallerBuilder,
        SignallerBuilder<TrajectorySignaller>*     trajectorySignallerBuilder,
        std::vector<ICheckpointHelperClient*>*     checkpointClients,
        CheckBondedInteractionsCallbackPtr*        checkBondedInteractionsCallback,
        compat::not_null<StatePropagatorData*>     statePropagatorDataPtr,
        compat::not_null<EnergyElement*>           energyElementPtr,
        FreeEnergyPerturbationElement*             freeEnergyPerturbationElementPtr,
        bool                                       hasReadEkinState,
        TopologyHolder*                            topologyHolder,
        SimulationSignals*                         signals)
{
    auto forceElement = buildForces(neighborSearchSignallerBuilder, energySignallerBuilder,
                                    statePropagatorDataPtr, energyElementPtr,
                                    freeEnergyPerturbationElementPtr, topologyHolder);

    // list of elements owned by the simulator composite object
    std::vector<std::unique_ptr<ISimulatorElement>> elementsOwnershipList;
    // call list of the simulator composite object
    std::vector<compat::not_null<ISimulatorElement*>> elementCallList;

    std::function<void()> needToCheckNumberOfBondedInteractions;
    if (inputrec->eI == eiMD)
    {
        auto computeGlobalsElement =
                std::make_unique<ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>>(
                        statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
                        signals, nstglobalcomm_, fplog, mdlog, cr, inputrec, mdAtoms, nrnb, wcycle,
                        fr, top_global, constr, hasReadEkinState);
        topologyHolder->registerClient(computeGlobalsElement.get());
        energySignallerBuilder->registerSignallerClient(compat::make_not_null(computeGlobalsElement.get()));
        trajectorySignallerBuilder->registerSignallerClient(
                compat::make_not_null(computeGlobalsElement.get()));

        *checkBondedInteractionsCallback =
                computeGlobalsElement->getCheckNumberOfBondedInteractionsCallback();

        auto propagator = std::make_unique<Propagator<IntegrationStep::LeapFrog>>(
                inputrec->delta_t, statePropagatorDataPtr, mdAtoms, wcycle);

        addToCallListAndMove(std::move(forceElement), elementCallList, elementsOwnershipList);
        addToCallList(statePropagatorDataPtr, elementCallList); // we have a full microstate at time t here!
        if (inputrec->etc == etcVRESCALE)
        {
            // TODO: With increased complexity of the propagator, this will need further development,
            //       e.g. using propagators templated for velocity propagation policies and a builder
            propagator->setNumVelocityScalingVariables(inputrec->opts.ngtc);
            auto thermostat = std::make_unique<VRescaleThermostat>(
                    inputrec->nsttcouple, -1, false, inputrec->ld_seed, inputrec->opts.ngtc,
                    inputrec->delta_t * inputrec->nsttcouple, inputrec->opts.ref_t, inputrec->opts.tau_t,
                    inputrec->opts.nrdf, energyElementPtr, propagator->viewOnVelocityScaling(),
                    propagator->velocityScalingCallback(), state_global, cr, inputrec->bContinuation);
            checkpointClients->emplace_back(thermostat.get());
            energyElementPtr->setVRescaleThermostat(thermostat.get());
            addToCallListAndMove(std::move(thermostat), elementCallList, elementsOwnershipList);
        }

        std::unique_ptr<ParrinelloRahmanBarostat> prBarostat = nullptr;
        if (inputrec->epc == epcPARRINELLORAHMAN)
        {
            // Building the PR barostat here since it needs access to the propagator
            // and we want to be able to move the propagator object
            prBarostat = std::make_unique<ParrinelloRahmanBarostat>(
                    inputrec->nstpcouple, -1, inputrec->delta_t * inputrec->nstpcouple,
                    inputrec->init_step, propagator->viewOnPRScalingMatrix(),
                    propagator->prScalingCallback(), statePropagatorDataPtr, energyElementPtr,
                    fplog, inputrec, mdAtoms, state_global, cr, inputrec->bContinuation);
            energyElementPtr->setParrinelloRahamnBarostat(prBarostat.get());
            checkpointClients->emplace_back(prBarostat.get());
        }
        addToCallListAndMove(std::move(propagator), elementCallList, elementsOwnershipList);
        if (constr)
        {
            auto constraintElement = std::make_unique<ConstraintsElement<ConstraintVariable::Positions>>(
                    constr, statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
                    MASTER(cr), fplog, inputrec, mdAtoms->mdatoms());
            auto constraintElementPtr = compat::make_not_null(constraintElement.get());
            energySignallerBuilder->registerSignallerClient(constraintElementPtr);
            trajectorySignallerBuilder->registerSignallerClient(constraintElementPtr);
            loggingSignallerBuilder->registerSignallerClient(constraintElementPtr);

            addToCallListAndMove(std::move(constraintElement), elementCallList, elementsOwnershipList);
        }

        addToCallListAndMove(std::move(computeGlobalsElement), elementCallList, elementsOwnershipList);
        addToCallList(energyElementPtr, elementCallList); // we have the energies at time t here!
        if (prBarostat)
        {
            addToCallListAndMove(std::move(prBarostat), elementCallList, elementsOwnershipList);
        }
    }
    else if (inputrec->eI == eiVV)
    {
        auto computeGlobalsElement =
                std::make_unique<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>(
                        statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
                        signals, nstglobalcomm_, fplog, mdlog, cr, inputrec, mdAtoms, nrnb, wcycle,
                        fr, &topologyHolder->globalTopology(), constr, hasReadEkinState);
        topologyHolder->registerClient(computeGlobalsElement.get());
        energySignallerBuilder->registerSignallerClient(compat::make_not_null(computeGlobalsElement.get()));
        trajectorySignallerBuilder->registerSignallerClient(
                compat::make_not_null(computeGlobalsElement.get()));

        *checkBondedInteractionsCallback =
                computeGlobalsElement->getCheckNumberOfBondedInteractionsCallback();

        auto propagatorVelocities = std::make_unique<Propagator<IntegrationStep::VelocitiesOnly>>(
                inputrec->delta_t * 0.5, statePropagatorDataPtr, mdAtoms, wcycle);
        auto propagatorVelocitiesAndPositions =
                std::make_unique<Propagator<IntegrationStep::VelocityVerletPositionsAndVelocities>>(
                        inputrec->delta_t, statePropagatorDataPtr, mdAtoms, wcycle);

        addToCallListAndMove(std::move(forceElement), elementCallList, elementsOwnershipList);

        std::unique_ptr<ParrinelloRahmanBarostat> prBarostat = nullptr;
        if (inputrec->epc == epcPARRINELLORAHMAN)
        {
            // Building the PR barostat here since it needs access to the propagator
            // and we want to be able to move the propagator object
            prBarostat = std::make_unique<ParrinelloRahmanBarostat>(
                    inputrec->nstpcouple, -1, inputrec->delta_t * inputrec->nstpcouple,
                    inputrec->init_step, propagatorVelocities->viewOnPRScalingMatrix(),
                    propagatorVelocities->prScalingCallback(), statePropagatorDataPtr, energyElementPtr,
                    fplog, inputrec, mdAtoms, state_global, cr, inputrec->bContinuation);
            energyElementPtr->setParrinelloRahamnBarostat(prBarostat.get());
            checkpointClients->emplace_back(prBarostat.get());
        }
        addToCallListAndMove(std::move(propagatorVelocities), elementCallList, elementsOwnershipList);
        if (constr)
        {
            auto constraintElement = std::make_unique<ConstraintsElement<ConstraintVariable::Velocities>>(
                    constr, statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
                    MASTER(cr), fplog, inputrec, mdAtoms->mdatoms());
            energySignallerBuilder->registerSignallerClient(compat::make_not_null(constraintElement.get()));
            trajectorySignallerBuilder->registerSignallerClient(
                    compat::make_not_null(constraintElement.get()));
            loggingSignallerBuilder->registerSignallerClient(
                    compat::make_not_null(constraintElement.get()));

            addToCallListAndMove(std::move(constraintElement), elementCallList, elementsOwnershipList);
        }
        addToCallList(compat::make_not_null(computeGlobalsElement.get()), elementCallList);
        addToCallList(statePropagatorDataPtr, elementCallList); // we have a full microstate at time t here!
        if (inputrec->etc == etcVRESCALE)
        {
            // TODO: With increased complexity of the propagator, this will need further development,
            //       e.g. using propagators templated for velocity propagation policies and a builder
            propagatorVelocitiesAndPositions->setNumVelocityScalingVariables(inputrec->opts.ngtc);
            auto thermostat = std::make_unique<VRescaleThermostat>(
                    inputrec->nsttcouple, 0, true, inputrec->ld_seed, inputrec->opts.ngtc,
                    inputrec->delta_t * inputrec->nsttcouple, inputrec->opts.ref_t,
                    inputrec->opts.tau_t, inputrec->opts.nrdf, energyElementPtr,
                    propagatorVelocitiesAndPositions->viewOnVelocityScaling(),
                    propagatorVelocitiesAndPositions->velocityScalingCallback(), state_global, cr,
                    inputrec->bContinuation);
            checkpointClients->emplace_back(thermostat.get());
            energyElementPtr->setVRescaleThermostat(thermostat.get());
            addToCallListAndMove(std::move(thermostat), elementCallList, elementsOwnershipList);
        }
        addToCallListAndMove(std::move(propagatorVelocitiesAndPositions), elementCallList,
                             elementsOwnershipList);
        if (constr)
        {
            auto constraintElement = std::make_unique<ConstraintsElement<ConstraintVariable::Positions>>(
                    constr, statePropagatorDataPtr, energyElementPtr, freeEnergyPerturbationElementPtr,
                    MASTER(cr), fplog, inputrec, mdAtoms->mdatoms());
            energySignallerBuilder->registerSignallerClient(compat::make_not_null(constraintElement.get()));
            trajectorySignallerBuilder->registerSignallerClient(
                    compat::make_not_null(constraintElement.get()));
            loggingSignallerBuilder->registerSignallerClient(
                    compat::make_not_null(constraintElement.get()));

            addToCallListAndMove(std::move(constraintElement), elementCallList, elementsOwnershipList);
        }
        addToCallListAndMove(std::move(computeGlobalsElement), elementCallList, elementsOwnershipList);
        addToCallList(energyElementPtr, elementCallList); // we have the energies at time t here!
        if (prBarostat)
        {
            addToCallListAndMove(std::move(prBarostat), elementCallList, elementsOwnershipList);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Integrator not implemented for the modular simulator.");
    }

    auto integrator = std::make_unique<CompositeSimulatorElement>(std::move(elementCallList),
                                                                  std::move(elementsOwnershipList));
    // std::move *should* not be needed with c++-14, but clang-3.6 still requires it
    return std::move(integrator);
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

void ModularSimulator::checkInputForDisabledFunctionality()
{
    isInputCompatible(true, inputrec, mdrunOptions.rerun, *top_global, ms, replExParams,
                      &fr->listedForces->fcdata(), opt2bSet("-ei", nfile, fnm), membed != nullptr);
    if (observablesHistory->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }
}

} // namespace gmx
