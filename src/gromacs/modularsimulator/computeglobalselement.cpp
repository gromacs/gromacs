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
 * \brief Defines the global reduction element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "computeglobalselement.h"

#include <any>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gromacs/domdec/domdec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdrun/isimulator.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/modularsimulator/energydata.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/modularsimulator/statepropagatordata.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

#include "freeenergyperturbationdata.h"
#include "modularsimulator.h"
#include "simulatoralgorithm.h"

struct t_forcerec;

namespace gmx
{
template<ComputeGlobalsAlgorithm algorithm>
ComputeGlobalsElement<algorithm>::ComputeGlobalsElement(StatePropagatorData* statePropagatorData,
                                                        EnergyData*          energyData,
                                                        FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                        SimulationSignals*  signals,
                                                        int                 nstglobalcomm,
                                                        FILE*               fplog,
                                                        const MDLogger&     mdlog,
                                                        t_commrec*          cr,
                                                        const t_inputrec*   inputrec,
                                                        const MDAtoms*      mdAtoms,
                                                        t_nrnb*             nrnb,
                                                        gmx_wallcycle*      wcycle,
                                                        t_forcerec*         fr,
                                                        const gmx_mtop_t&   global_top,
                                                        Constraints*        constr,
                                                        ObservablesReducer* observablesReducer) :
    energyReductionStep_(-1),
    virialReductionStep_(-1),
    vvSchedulingStep_(-1),
    doStopCM_(inputrec->comm_mode != ComRemovalAlgorithm::No),
    nstcomm_(inputrec->nstcomm),
    nstglobalcomm_(nstglobalcomm),
    lastStep_(inputrec->nsteps + inputrec->init_step),
    initStep_(inputrec->init_step),
    nullSignaller_(std::make_unique<SimulationSignaller>(nullptr, nullptr, nullptr, false, false)),
    statePropagatorData_(statePropagatorData),
    energyData_(energyData),
    freeEnergyPerturbationData_(freeEnergyPerturbationData),
    vcm_(global_top.groups, *inputrec),
    signals_(signals),
    fplog_(fplog),
    mdlog_(mdlog),
    cr_(cr),
    inputrec_(inputrec),
    top_global_(global_top),
    mdAtoms_(mdAtoms),
    constr_(constr),
    nrnb_(nrnb),
    wcycle_(wcycle),
    fr_(fr),
    observablesReducer_(observablesReducer)
{
    reportComRemovalInfo(fplog, vcm_);
    gstat_ = global_stat_init(inputrec_);
}

template<ComputeGlobalsAlgorithm algorithm>
ComputeGlobalsElement<algorithm>::~ComputeGlobalsElement()
{
    global_stat_destroy(gstat_);
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::elementSetup()
{
    if (doStopCM_ && !inputrec_->bContinuation)
    {
        // To minimize communication, compute_globals computes the COM velocity
        // and the kinetic energy for the velocities without COM motion removed.
        // Thus to get the kinetic energy without the COM contribution, we need
        // to call compute_globals twice.

        compute(-1, CGLO_GSTAT | CGLO_STOPCM, nullSignaller_.get(), false, true);
        // Clean up after pre-step use of compute()
        observablesReducer_->markAsReadyToReduce();

        auto v = statePropagatorData_->velocitiesView();
        // At initialization, do not pass x with acceleration-correction mode
        // to avoid (incorrect) correction of the initial coordinates.
        auto x = vcm_.mode == ComRemovalAlgorithm::LinearAccelerationCorrection
                         ? ArrayRefWithPadding<RVec>()
                         : statePropagatorData_->positionsView();
        process_and_stopcm_grp(
                fplog_, &vcm_, *mdAtoms_->mdatoms(), x.unpaddedArrayRef(), v.unpaddedArrayRef());
        inc_nrnb(nrnb_, eNR_STOPCM, mdAtoms_->mdatoms()->homenr);
    }

    unsigned int cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT
                               | (energyData_->hasReadEkinFromCheckpoint() ? CGLO_READEKIN : 0));

    if (algorithm == ComputeGlobalsAlgorithm::VelocityVerlet)
    {
        cglo_flags |= CGLO_PRESSURE | CGLO_CONSTRAINT;
    }

    compute(-1, cglo_flags, nullSignaller_.get(), false, true);

    // Calculate the initial half step temperature, and save the ekinh_old
    for (int i = 0; (i < inputrec_->opts.ngtc); i++)
    {
        copy_mat(energyData_->ekindata()->tcstat[i].ekinh, energyData_->ekindata()->tcstat[i].ekinh_old);
    }

    // Clean up after pre-step use of compute()
    observablesReducer_->markAsReadyToReduce();
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::scheduleTask(Step                       step,
                                                    Time gmx_unused            time,
                                                    const RegisterRunFunction& registerRunFunction)
{
    const bool needComReduction    = doStopCM_ && do_per_step(step, nstcomm_);
    const bool needGlobalReduction = step == energyReductionStep_ || step == virialReductionStep_
                                     || needComReduction || do_per_step(step, nstglobalcomm_)
                                     || (EI_VV(inputrec_->eI) && inputrecNvtTrotter(inputrec_)
                                         && do_per_step(step - 1, nstglobalcomm_));

    // TODO: CGLO_GSTAT is only used for needToSumEkinhOld_, i.e. to signal that we do or do not
    //       sum the previous kinetic energy. We should simplify / clarify this.

    if (algorithm == ComputeGlobalsAlgorithm::LeapFrog)
    {
        // With Leap-Frog we can skip compute_globals at
        // non-communication steps, but we need to calculate
        // the kinetic energy one step before communication.

        // With leap-frog we also need to compute the half-step
        // kinetic energy at the step before we need to write
        // the full-step kinetic energy
        const bool needEkinAtNextStep = (do_per_step(step + 1, nstglobalcomm_) || step + 1 == lastStep_);

        if (!needGlobalReduction && !needEkinAtNextStep)
        {
            return;
        }

        const bool doEnergy = step == energyReductionStep_;
        int        flags    = (needGlobalReduction ? CGLO_GSTAT : 0) | (doEnergy ? CGLO_ENERGY : 0)
                    | (needComReduction ? CGLO_STOPCM : 0) | CGLO_TEMPERATURE | CGLO_PRESSURE
                    | CGLO_CONSTRAINT;

        // Since we're already communicating at this step, we
        // can propagate intra-simulation signals. Note that
        // check_nstglobalcomm has the responsibility for
        // choosing the value of nstglobalcomm which satisfies
        // the need of the different signallers.
        const bool doIntraSimSignal = true;
        // Disable functionality
        const bool doInterSimSignal = false;

        // Make signaller to signal stop / reset / checkpointing signals
        auto signaller = std::make_shared<SimulationSignaller>(
                signals_, cr_, nullptr, doInterSimSignal, doIntraSimSignal);

        registerRunFunction([this, step, flags, signaller = std::move(signaller)]() {
            compute(step, flags, signaller.get(), true);
        });
    }
    else if (algorithm == ComputeGlobalsAlgorithm::VelocityVerlet)
    {
        // For VV, we schedule two calls to compute globals per step.
        if (step != vvSchedulingStep_)
        {
            // This is the first scheduling call for this step (positions & velocities at full time
            // step) Set this as the current scheduling step
            vvSchedulingStep_ = step;

            // For vv, the state at the beginning of the step is positions at time t, velocities at time t - dt/2
            // The first velocity propagation (+dt/2) therefore actually corresponds to the previous step.
            // So we need information from the last step in the first half of the integration
            if (!needGlobalReduction && !do_per_step(step - 1, nstglobalcomm_))
            {
                return;
            }

            const bool doTemperature = step != initStep_ || inputrec_->bContinuation;
            const bool doEnergy      = step == energyReductionStep_;

            int flags = (needGlobalReduction ? CGLO_GSTAT : 0) | (doEnergy ? CGLO_ENERGY : 0)
                        | (doTemperature ? CGLO_TEMPERATURE : 0) | CGLO_PRESSURE | CGLO_CONSTRAINT
                        | (needComReduction ? CGLO_STOPCM : 0) | CGLO_SCALEEKIN;

            registerRunFunction(
                    [this, step, flags]() { compute(step, flags, nullSignaller_.get(), false); });
        }
        else
        {
            // second call to compute_globals for this step
            // Reset the scheduling step to avoid confusion if scheduling needs
            // to be repeated (in case of unexpected simulation termination)
            vvSchedulingStep_ = -1;

            if (!needGlobalReduction)
            {
                return;
            }
            int flags = CGLO_GSTAT | CGLO_CONSTRAINT;

            // Since we're already communicating at this step, we
            // can propagate intra-simulation signals. Note that
            // check_nstglobalcomm has the responsibility for
            // choosing the value of nstglobalcomm which satisfies
            // the need of the different signallers.
            const bool doIntraSimSignal = true;
            // Disable functionality
            const bool doInterSimSignal = false;

            auto signaller = std::make_shared<SimulationSignaller>(
                    signals_, cr_, nullptr, doInterSimSignal, doIntraSimSignal);

            registerRunFunction([this, step, flags, signaller = std::move(signaller)]() {
                compute(step, flags, signaller.get(), true);
            });
        }
    }
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::compute(gmx::Step            step,
                                               unsigned int         flags,
                                               SimulationSignaller* signaller,
                                               bool                 useLastBox,
                                               bool                 isInit)
{
    auto        x       = statePropagatorData_->positionsView().unpaddedArrayRef();
    auto        v       = statePropagatorData_->velocitiesView().unpaddedArrayRef();
    const auto* box     = statePropagatorData_->constBox();
    const auto* lastbox = useLastBox ? statePropagatorData_->constPreviousBox()
                                     : statePropagatorData_->constBox();

    compute_globals(gstat_,
                    cr_,
                    inputrec_,
                    fr_,
                    energyData_->ekindata(),
                    x,
                    v,
                    box,
                    mdAtoms_->mdatoms(),
                    nrnb_,
                    &vcm_,
                    step != -1 ? wcycle_ : nullptr,
                    energyData_->enerdata(),
                    energyData_->forceVirial(step),
                    energyData_->constraintVirial(step),
                    energyData_->totalVirial(step),
                    energyData_->pressure(step),
                    signaller,
                    lastbox,
                    energyData_->needToSumEkinhOld(),
                    flags,
                    step,
                    observablesReducer_);
    if (flags & CGLO_STOPCM && !isInit)
    {
        process_and_stopcm_grp(fplog_, &vcm_, *mdAtoms_->mdatoms(), x, v);
        inc_nrnb(nrnb_, eNR_STOPCM, mdAtoms_->mdatoms()->homenr);
    }
}

template<ComputeGlobalsAlgorithm algorithm>
std::optional<SignallerCallback> ComputeGlobalsElement<algorithm>::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { energyReductionStep_ = step; };
    }
    if (event == EnergySignallerEvent::VirialCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { virialReductionStep_ = step; };
    }
    return std::nullopt;
}

template<ComputeGlobalsAlgorithm algorithm>
std::optional<SignallerCallback>
ComputeGlobalsElement<algorithm>::registerTrajectorySignallerCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep)
    {
        return [this](Step step, Time /*unused*/) { energyReductionStep_ = step; };
    }
    return std::nullopt;
}

namespace
{

/*! \brief Schedule a function for actions that must happen at the end of each step
 *
 * After reduction, an ObservablesReducer is marked as unavailable for
 * further reduction this step. This needs to be reset in order to be
 * used on the next step.
 *
 * \param[in]  observablesReducer The ObservablesReducer to mark as ready for use
 */
SchedulingFunction registerPostStepSchedulingFunction(ObservablesReducer* observablesReducer)
{
    SchedulingFunction postStepSchedulingFunction =
            [observablesReducer](
                    Step /*step*/, Time /*time*/, const RegisterRunFunction& registerRunFunction) {
                SimulatorRunFunction completeObservablesReducerStep = [&observablesReducer]() {
                    observablesReducer->markAsReadyToReduce();
                };
                registerRunFunction(completeObservablesReducerStep);
            };
    return postStepSchedulingFunction;
}

} // namespace

//! Explicit template instantiation
//! \{
template class ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>;
template class ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>;
//! \}

template<>
ISimulatorElement* ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData*             freeEnergyPerturbationData,
        GlobalCommunicationHelper*              globalCommunicationHelper,
        ObservablesReducer*                     observablesReducer)
{
    ComputeGlobalsElement* element = builderHelper->storeElement(
            std::make_unique<ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>>(
                    statePropagatorData,
                    energyData,
                    freeEnergyPerturbationData,
                    globalCommunicationHelper->simulationSignals(),
                    globalCommunicationHelper->nstglobalcomm(),
                    legacySimulatorData->fpLog_,
                    legacySimulatorData->mdLog_,
                    legacySimulatorData->cr_,
                    legacySimulatorData->inputRec_,
                    legacySimulatorData->mdAtoms_,
                    legacySimulatorData->nrnb_,
                    legacySimulatorData->wallCycleCounters_,
                    legacySimulatorData->fr_,
                    legacySimulatorData->topGlobal_,
                    legacySimulatorData->constr_,
                    observablesReducer));
    builderHelper->registerPostStepScheduling(
            registerPostStepSchedulingFunction(element->observablesReducer_));

    return element;
}

template<>
ISimulatorElement* ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>::getElementPointerImpl(
        LegacySimulatorData*                    simulator,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData*             freeEnergyPerturbationData,
        GlobalCommunicationHelper*              globalCommunicationHelper,
        ObservablesReducer*                     observablesReducer)
{
    // We allow this element to be added multiple times to the call list, but we only want one
    // actual element built
    static const std::string key("vvComputeGlobalsElement");

    const std::optional<std::any> cachedValue = builderHelper->builderData(key);

    if (cachedValue)
    {
        return std::any_cast<ComputeGlobalsElement*>(cachedValue.value());
    }
    else
    {
        ComputeGlobalsElement* vvComputeGlobalsElement = builderHelper->storeElement(
                std::make_unique<ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>>(
                        statePropagatorData,
                        energyData,
                        freeEnergyPerturbationData,
                        globalCommunicationHelper->simulationSignals(),
                        globalCommunicationHelper->nstglobalcomm(),
                        simulator->fpLog_,
                        simulator->mdLog_,
                        simulator->cr_,
                        simulator->inputRec_,
                        simulator->mdAtoms_,
                        simulator->nrnb_,
                        simulator->wallCycleCounters_,
                        simulator->fr_,
                        simulator->topGlobal_,
                        simulator->constr_,
                        observablesReducer));
        builderHelper->storeBuilderData(key, vvComputeGlobalsElement);
        builderHelper->registerPostStepScheduling(
                registerPostStepSchedulingFunction(vvComputeGlobalsElement->observablesReducer_));
        return vvComputeGlobalsElement;
    }
}
} // namespace gmx
