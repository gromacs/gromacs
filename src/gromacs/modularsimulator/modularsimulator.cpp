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
 * \brief Defines the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "modularsimulator.h"

#include <filesystem>
#include <optional>
#include <string_view>

#include "gromacs/commandline/filenm.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdrun/replicaexchange.h"
#include "gromacs/mdrun/shellfc.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/modularsimulator/energydata.h"
#include "gromacs/modularsimulator/freeenergyperturbationdata.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/modularsimulator/propagator.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/int64_to_int.h"
#include "gromacs/utility/logger.h"

#include "andersentemperaturecoupling.h"
#include "computeglobalselement.h"
#include "constraintelement.h"
#include "expandedensembleelement.h"
#include "firstorderpressurecoupling.h"
#include "forceelement.h"
#include "mttk.h"
#include "nosehooverchains.h"
#include "parrinellorahmanbarostat.h"
#include "pullelement.h"
#include "simulatoralgorithm.h"
#include "statepropagatordata.h"
#include "velocityscalingtemperaturecoupling.h"

struct gmx_multisim_t;

namespace gmx
{
void ModularSimulator::run()
{
    GMX_LOG(legacySimulatorData_->mdLog_.info)
            .asParagraph()
            .appendText("Using the modular simulator.");

    ModularSimulatorAlgorithmBuilder algorithmBuilder(compat::make_not_null(legacySimulatorData_),
                                                      std::move(checkpointDataHolder_));
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
    const bool isTrotter = inputrecNvtTrotter(legacySimulatorData_->inputRec_)
                           || inputrecNptTrotter(legacySimulatorData_->inputRec_)
                           || inputrecNphTrotter(legacySimulatorData_->inputRec_);
    if (legacySimulatorData_->inputRec_->eI == IntegrationAlgorithm::MD)
    {
        // The leap frog integration algorithm
        builder->add<ForceElement>();
        builder->add<StatePropagatorData::Element>();
        if (legacySimulatorData_->inputRec_->etc == TemperatureCoupling::VRescale
            || legacySimulatorData_->inputRec_->etc == TemperatureCoupling::Berendsen
            || legacySimulatorData_->inputRec_->etc == TemperatureCoupling::NoseHoover)
        {
            builder->add<VelocityScalingTemperatureCoupling>(Offset(-1),
                                                             UseFullStepKE::No,
                                                             ReportPreviousStepConservedEnergy::No,
                                                             PropagatorTag("LeapFrogPropagator"));
        }
        builder->add<Propagator<IntegrationStage::LeapFrog>>(
                PropagatorTag("LeapFrogPropagator"), TimeStep(legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->constr_)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Positions>>();
        }

        if (legacySimulatorData_->inputRec_->bPull)
        {
            builder->add<PullElement>();
        }

        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>>();
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman)
        {
            builder->add<ParrinelloRahmanBarostat>(Offset(-1), PropagatorTag("LeapFrogPropagator"));
        }
        else if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Berendsen
                 || legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::CRescale)
        {
            builder->add<FirstOrderPressureCoupling>(0, ReportPreviousStepConservedEnergy::No);
        }
    }
    else if (legacySimulatorData_->inputRec_->eI == IntegrationAlgorithm::VV && !isTrotter)
    {
        // The velocity verlet integration algorithm
        builder->add<ForceElement>();
        builder->add<Propagator<IntegrationStage::VelocitiesOnly>>(
                PropagatorTag("VelocityHalfStep"),
                TimeStep(0.5 * legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->constr_)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Velocities>>();
        }
        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();
        // Here, we have x / v / f at the full time step
        builder->add<StatePropagatorData::Element>();
        if (legacySimulatorData_->inputRec_->bExpanded)
        {
            builder->add<ExpandedEnsembleElement>();
        }
        if (legacySimulatorData_->inputRec_->etc == TemperatureCoupling::VRescale
            || legacySimulatorData_->inputRec_->etc == TemperatureCoupling::Berendsen)
        {
            builder->add<VelocityScalingTemperatureCoupling>(
                    Offset(0),
                    UseFullStepKE::Yes,
                    ReportPreviousStepConservedEnergy::Yes,
                    PropagatorTag("VelocityHalfAndPositionFullStep"));
        }
        else if (ETC_ANDERSEN(legacySimulatorData_->inputRec_->etc))
        {
            builder->add<AndersenTemperatureCoupling>();
        }
        builder->add<Propagator<IntegrationStage::VelocityVerletPositionsAndVelocities>>(
                PropagatorTag("VelocityHalfAndPositionFullStep"),
                TimeStep(legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->constr_)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Positions>>();
        }

        if (legacySimulatorData_->inputRec_->bPull)
        {
            builder->add<PullElement>();
        }

        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman)
        {
            builder->add<ParrinelloRahmanBarostat>(Offset(-1), PropagatorTag("VelocityHalfStep"));
        }
        else if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Berendsen
                 || legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::CRescale)
        {
            builder->add<FirstOrderPressureCoupling>(0, ReportPreviousStepConservedEnergy::Yes);
        }
    }
    else if (legacySimulatorData_->inputRec_->eI == IntegrationAlgorithm::VV && isTrotter)
    {
        // For a new simulation, avoid the first Trotter half step
        const auto scheduleTrotterFirstHalfOnInitStep =
                ((legacySimulatorData_->startingBehavior_ == StartingBehavior::NewSimulation)
                         ? ScheduleOnInitStep::No
                         : ScheduleOnInitStep::Yes);
        // Define the tags and offsets for MTTK pressure scaling
        const MttkPropagatorConnectionDetails mttkPropagatorConnectionDetails = {
            PropagatorTag("ScaleMTTKXPre"),  PropagatorTag("ScaleMTTKXPost"),  Offset(0),
            PropagatorTag("ScaleMTTKVPre1"), PropagatorTag("ScaleMTTKVPost1"), Offset(1),
            PropagatorTag("ScaleMTTKVPre2"), PropagatorTag("ScaleMTTKVPost2"), Offset(0)
        };

        builder->add<ForceElement>();
        // Propagate velocities from t-dt/2 to t
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<Propagator<IntegrationStage::ScaleVelocities>>(
                    PropagatorTag("ScaleMTTKVPre1"));
        }
        builder->add<Propagator<IntegrationStage::VelocitiesOnly>>(
                PropagatorTag("VelocityHalfStep1"),
                TimeStep(0.5 * legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<Propagator<IntegrationStage::ScaleVelocities>>(
                    PropagatorTag("ScaleMTTKVPost1"));
        }
        if (legacySimulatorData_->constr_)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Velocities>>();
        }
        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();

        // Propagate extended system variables from t-dt/2 to t
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<MttkElement>(
                    Offset(-1), scheduleTrotterFirstHalfOnInitStep, mttkPropagatorConnectionDetails);
        }
        if (legacySimulatorData_->inputRec_->etc == TemperatureCoupling::NoseHoover)
        {
            builder->add<NoseHooverChainsElement>(NhcUsage::System,
                                                  Offset(-1),
                                                  UseFullStepKE::Yes,
                                                  scheduleTrotterFirstHalfOnInitStep,
                                                  PropagatorTag("ScaleNHC"));
            builder->add<Propagator<IntegrationStage::ScaleVelocities>>(PropagatorTag("ScaleNHC"));
        }
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<NoseHooverChainsElement>(NhcUsage::Barostat,
                                                  Offset(-1),
                                                  UseFullStepKE::Yes,
                                                  scheduleTrotterFirstHalfOnInitStep,
                                                  mttkPropagatorConnectionDetails);
        }
        // We have a full state at time t here
        builder->add<StatePropagatorData::Element>();
        if (legacySimulatorData_->inputRec_->bExpanded)
        {
            builder->add<ExpandedEnsembleElement>();
        }

        // Propagate extended system variables from t to t+dt/2
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<NoseHooverChainsElement>(NhcUsage::Barostat,
                                                  Offset(0),
                                                  UseFullStepKE::Yes,
                                                  ScheduleOnInitStep::Yes,
                                                  mttkPropagatorConnectionDetails);
        }
        if (legacySimulatorData_->inputRec_->etc == TemperatureCoupling::NoseHoover)
        {
            builder->add<NoseHooverChainsElement>(NhcUsage::System,
                                                  Offset(0),
                                                  UseFullStepKE::Yes,
                                                  ScheduleOnInitStep::Yes,
                                                  PropagatorTag("VelocityHalfStep2"));
        }
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<MttkElement>(Offset(0), ScheduleOnInitStep::Yes, mttkPropagatorConnectionDetails);
            builder->add<Propagator<IntegrationStage::ScaleVelocities>>(
                    PropagatorTag("ScaleMTTKVPre2"));
        }

        // Propagate velocities from t to t+dt/2
        builder->add<Propagator<IntegrationStage::VelocitiesOnly>>(
                PropagatorTag("VelocityHalfStep2"),
                TimeStep(0.5 * legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<Propagator<IntegrationStage::ScaleVelocities>>(
                    PropagatorTag("ScaleMTTKVPost2"));
            builder->add<Propagator<IntegrationStage::ScalePositions>>(
                    PropagatorTag("ScaleMTTKXPre"));
        }
        // Propagate positions from t to t+dt
        builder->add<Propagator<IntegrationStage::PositionsOnly>>(
                PropagatorTag("PositionFullStep"), TimeStep(legacySimulatorData_->inputRec_->delta_t));
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<Propagator<IntegrationStage::ScalePositions>>(
                    PropagatorTag("ScaleMTTKXPost"));
        }
        if (legacySimulatorData_->constr_)
        {
            builder->add<ConstraintsElement<ConstraintVariable::Positions>>();
        }

        if (legacySimulatorData_->inputRec_->bPull)
        {
            builder->add<PullElement>();
        }

        builder->add<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>();

        // Propagate box from t to t+dt
        if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        {
            builder->add<MttkBoxScaling>(mttkPropagatorConnectionDetails);
        }
        else if (legacySimulatorData_->inputRec_->pressureCouplingOptions.epc == PressureCoupling::CRescale)
        {
            // Legacy implementation allows combination of C-Rescale with Trotter Nose-Hoover
            builder->add<FirstOrderPressureCoupling>(0, ReportPreviousStepConservedEnergy::Yes);
        }
    }
    else
    {
        gmx_fatal(FARGS, "Integrator not implemented for the modular simulator.");
    }
    builder->add<EnergyData::Element>();
}

bool ModularSimulator::isInputCompatible(bool                             exitOnFailure,
                                         const t_inputrec*                inputrec,
                                         bool                             doRerun,
                                         const gmx_mtop_t&                globalTopology,
                                         const gmx_multisim_t*            ms,
                                         const ReplicaExchangeParameters& replExParams,
                                         const t_fcdata*                  fcd,
                                         bool                             doEssentialDynamics,
                                         bool                             doMembed,
                                         bool                             useGpuForUpdate)
{
    auto conditionalAssert = [exitOnFailure](bool condition, const char* message) {
        if (exitOnFailure)
        {
            GMX_RELEASE_ASSERT(condition, message);
        }
        return condition;
    };

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
            !(modularSimulatorExplicitlyTurnedOff && inputrec->eI == IntegrationAlgorithm::VV
              && inputrec->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman),
            "Cannot use a Parrinello-Rahman barostat with md-vv and "
            "GMX_DISABLE_MODULAR_SIMULATOR=ON, "
            "as the Parrinello-Rahman barostat is not implemented in the legacy simulator. Unset "
            "GMX_DISABLE_MODULAR_SIMULATOR or use a different pressure control algorithm.");

    bool isInputCompatible = conditionalAssert(
            inputrec->eI == IntegrationAlgorithm::MD || inputrec->eI == IntegrationAlgorithm::VV,
            "Only integrators md and md-vv are supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->eI != IntegrationAlgorithm::MD
                                                     || modularSimulatorExplicitlyTurnedOn,
                                             "Set GMX_USE_MODULAR_SIMULATOR=ON to use the modular "
                                             "simulator with integrator md.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                    !inputrec->useMts,
                    "Multiple time stepping is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doRerun, "Rerun is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!inputrec->useConstantAcceleration && inputrec->cos_accel == 0.0,
                                 "Acceleration is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!inputrecFrozenAtoms(inputrec),
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
        auto* distantRestraintEnsembleEnvVar = getenv("GMX_DISRE_ENSEMBLE_SIZE");
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
            && conditionalAssert(!doSimulatedAnnealing(*inputrec),
                                 "Simulated annealing is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!inputrec->bSimTemp,
                                 "Simulated tempering is not supported by the modular simulator.");
    isInputCompatible =
            isInputCompatible
            && conditionalAssert(!doEssentialDynamics,
                                 "Essential dynamics is not supported by the modular simulator.");
    isInputCompatible = isInputCompatible
                        && conditionalAssert(inputrec->eSwapCoords == SwapType::No,
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

    isInputCompatible =
            isInputCompatible
            && conditionalAssert(
                    !useGpuForUpdate,
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
    if (!isInputCompatible
        && (inputrec->eI == IntegrationAlgorithm::VV
            && inputrec->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman))
    {
        gmx_fatal(FARGS,
                  "Requested Parrinello-Rahman barostat with md-vv. This combination is only "
                  "available in the modular simulator. Some other selected options are, however, "
                  "only available in the legacy simulator. Use a different pressure control "
                  "algorithm.");
    }
    return isInputCompatible;
}

ModularSimulator::ModularSimulator(std::unique_ptr<LegacySimulatorData>      legacySimulatorData,
                                   std::unique_ptr<ReadCheckpointDataHolder> checkpointDataHolder) :
    legacySimulatorData_(std::move(legacySimulatorData)),
    checkpointDataHolder_(std::move(checkpointDataHolder))
{
    checkInputForDisabledFunctionality();
}

ModularSimulator::~ModularSimulator() = default;

void ModularSimulator::checkInputForDisabledFunctionality()
{
    isInputCompatible(true,
                      legacySimulatorData_->inputRec_,
                      legacySimulatorData_->mdrunOptions_.rerun,
                      legacySimulatorData_->topGlobal_,
                      legacySimulatorData_->ms_,
                      legacySimulatorData_->replExParams_,
                      legacySimulatorData_->fr_->fcdata.get(),
                      opt2bSet("-ei", legacySimulatorData_->nFile_, legacySimulatorData_->fnm_),
                      legacySimulatorData_->membed_ != nullptr,
                      false);
    if (legacySimulatorData_->observablesHistory_->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }
}

void ModularSimulator::readCheckpointToTrxFrame(t_trxframe*               fr,
                                                ReadCheckpointDataHolder* readCheckpointDataHolder,
                                                const CheckpointHeaderContents& checkpointHeaderContents)
{
    GMX_RELEASE_ASSERT(checkpointHeaderContents.isModularSimulatorCheckpoint,
                       "ModularSimulator::readCheckpointToTrxFrame can only read checkpoints "
                       "written by modular simulator.");
    fr->bStep = true;
    fr->step = int64_to_int(checkpointHeaderContents.step, "conversion of checkpoint to trajectory");
    fr->bTime = true;
    fr->time  = checkpointHeaderContents.t;

    fr->bAtoms = false;

    StatePropagatorData::readCheckpointToTrxFrame(
            fr, readCheckpointDataHolder->checkpointData(StatePropagatorData::checkpointID()));
    if (readCheckpointDataHolder->keyExists(FreeEnergyPerturbationData::checkpointID()))
    {
        FreeEnergyPerturbationData::readCheckpointToTrxFrame(
                fr, readCheckpointDataHolder->checkpointData(FreeEnergyPerturbationData::checkpointID()));
    }
    else
    {
        FreeEnergyPerturbationData::readCheckpointToTrxFrame(fr, std::nullopt);
    }
}

} // namespace gmx
