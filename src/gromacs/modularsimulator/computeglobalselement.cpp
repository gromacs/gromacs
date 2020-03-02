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
 * \brief Defines the global reduction element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "computeglobalselement.h"

#include "gromacs/domdec/partition.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/topology.h"

#include "freeenergyperturbationelement.h"

namespace gmx
{
template<ComputeGlobalsAlgorithm algorithm>
ComputeGlobalsElement<algorithm>::ComputeGlobalsElement(StatePropagatorData* statePropagatorData,
                                                        EnergyElement*       energyElement,
                                                        FreeEnergyPerturbationElement* freeEnergyPerturbationElement,
                                                        SimulationSignals* signals,
                                                        int                nstglobalcomm,
                                                        FILE*              fplog,
                                                        const MDLogger&    mdlog,
                                                        t_commrec*         cr,
                                                        const t_inputrec*  inputrec,
                                                        const MDAtoms*     mdAtoms,
                                                        t_nrnb*            nrnb,
                                                        gmx_wallcycle*     wcycle,
                                                        t_forcerec*        fr,
                                                        const gmx_mtop_t*  global_top,
                                                        Constraints*       constr,
                                                        bool               hasReadEkinState) :
    energyReductionStep_(-1),
    virialReductionStep_(-1),
    vvSchedulingStep_(-1),
    doStopCM_(inputrec->comm_mode != ecmNO),
    nstcomm_(inputrec->nstcomm),
    nstglobalcomm_(nstglobalcomm),
    lastStep_(inputrec->nsteps + inputrec->init_step),
    initStep_(inputrec->init_step),
    nullSignaller_(std::make_unique<SimulationSignaller>(nullptr, nullptr, nullptr, false, false)),
    hasReadEkinState_(hasReadEkinState),
    totalNumberOfBondedInteractions_(0),
    shouldCheckNumberOfBondedInteractions_(false),
    statePropagatorData_(statePropagatorData),
    energyElement_(energyElement),
    localTopology_(nullptr),
    freeEnergyPerturbationElement_(freeEnergyPerturbationElement),
    vcm_(global_top->groups, *inputrec),
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
    fr_(fr)
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
    GMX_ASSERT(localTopology_, "Setup called before local topology was set.");

    if (doStopCM_ && !inputrec_->bContinuation)
    {
        // To minimize communication, compute_globals computes the COM velocity
        // and the kinetic energy for the velocities without COM motion removed.
        // Thus to get the kinetic energy without the COM contribution, we need
        // to call compute_globals twice.

        compute(-1, CGLO_GSTAT | CGLO_STOPCM, nullSignaller_.get(), false, true);

        auto v = statePropagatorData_->velocitiesView();
        // At initialization, do not pass x with acceleration-correction mode
        // to avoid (incorrect) correction of the initial coordinates.
        auto x = vcm_.mode == ecmLINEAR_ACCELERATION_CORRECTION ? ArrayRefWithPadding<RVec>()
                                                                : statePropagatorData_->positionsView();
        process_and_stopcm_grp(fplog_, &vcm_, *mdAtoms_->mdatoms(), x.unpaddedArrayRef(),
                               v.unpaddedArrayRef());
        inc_nrnb(nrnb_, eNR_STOPCM, mdAtoms_->mdatoms()->homenr);
    }

    unsigned int cglo_flags = (CGLO_TEMPERATURE | CGLO_GSTAT | (hasReadEkinState_ ? CGLO_READEKIN : 0));

    if (algorithm == ComputeGlobalsAlgorithm::VelocityVerlet)
    {
        cglo_flags |= CGLO_PRESSURE | CGLO_CONSTRAINT;
    }

    compute(-1, cglo_flags, nullSignaller_.get(), false, true);

    // Calculate the initial half step temperature, and save the ekinh_old
    for (int i = 0; (i < inputrec_->opts.ngtc); i++)
    {
        copy_mat(energyElement_->ekindata()->tcstat[i].ekinh,
                 energyElement_->ekindata()->tcstat[i].ekinh_old);
    }
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::scheduleTask(Step step,
                                                    Time gmx_unused               time,
                                                    const RegisterRunFunctionPtr& registerRunFunction)
{
    const bool needComReduction    = doStopCM_ && do_per_step(step, nstcomm_);
    const bool needGlobalReduction = step == energyReductionStep_ || step == virialReductionStep_
                                     || needComReduction || do_per_step(step, nstglobalcomm_);

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
        auto signaller = std::make_shared<SimulationSignaller>(signals_, cr_, nullptr,
                                                               doInterSimSignal, doIntraSimSignal);

        (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                [this, step, flags, signaller = std::move(signaller)]() {
                    compute(step, flags, signaller.get(), true);
                }));
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

            (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                    [this, step, flags]() { compute(step, flags, nullSignaller_.get(), false); }));
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

            (*registerRunFunction)(std::make_unique<SimulatorRunFunction>(
                    [this, step, flags, signaller = std::move(signaller)]() {
                        compute(step, flags, signaller.get(), true);
                    }));
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
    auto x       = statePropagatorData_->positionsView().unpaddedArrayRef();
    auto v       = statePropagatorData_->velocitiesView().unpaddedArrayRef();
    auto box     = statePropagatorData_->constBox();
    auto lastbox = useLastBox ? statePropagatorData_->constPreviousBox()
                              : statePropagatorData_->constBox();

    const real vdwLambda = freeEnergyPerturbationElement_
                                   ? freeEnergyPerturbationElement_->constLambdaView()[efptVDW]
                                   : 0;

    compute_globals(
            gstat_, cr_, inputrec_, fr_, energyElement_->ekindata(), x, v, box, vdwLambda,
            mdAtoms_->mdatoms(), nrnb_, &vcm_, step != -1 ? wcycle_ : nullptr, energyElement_->enerdata(),
            energyElement_->forceVirial(step), energyElement_->constraintVirial(step),
            energyElement_->totalVirial(step), energyElement_->pressure(step), constr_, signaller,
            lastbox, &totalNumberOfBondedInteractions_, energyElement_->needToSumEkinhOld(),
            flags | (shouldCheckNumberOfBondedInteractions_ ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS : 0));
    checkNumberOfBondedInteractions(mdlog_, cr_, totalNumberOfBondedInteractions_, top_global_,
                                    localTopology_, x, box, &shouldCheckNumberOfBondedInteractions_);
    if (flags & CGLO_STOPCM && !isInit)
    {
        process_and_stopcm_grp(fplog_, &vcm_, *mdAtoms_->mdatoms(), x, v);
        inc_nrnb(nrnb_, eNR_STOPCM, mdAtoms_->mdatoms()->homenr);
    }
}

template<ComputeGlobalsAlgorithm algorithm>
CheckBondedInteractionsCallbackPtr ComputeGlobalsElement<algorithm>::getCheckNumberOfBondedInteractionsCallback()
{
    return std::make_unique<CheckBondedInteractionsCallback>(
            [this]() { needToCheckNumberOfBondedInteractions(); });
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::needToCheckNumberOfBondedInteractions()
{
    shouldCheckNumberOfBondedInteractions_ = true;
}

template<ComputeGlobalsAlgorithm algorithm>
void ComputeGlobalsElement<algorithm>::setTopology(const gmx_localtop_t* top)
{
    localTopology_ = top;
}

template<ComputeGlobalsAlgorithm algorithm>
SignallerCallbackPtr ComputeGlobalsElement<algorithm>::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { energyReductionStep_ = step; });
    }
    if (event == EnergySignallerEvent::VirialCalculationStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { virialReductionStep_ = step; });
    }
    return nullptr;
}

template<ComputeGlobalsAlgorithm algorithm>
SignallerCallbackPtr ComputeGlobalsElement<algorithm>::registerTrajectorySignallerCallback(TrajectoryEvent event)
{
    if (event == TrajectoryEvent::EnergyWritingStep)
    {
        return std::make_unique<SignallerCallback>(
                [this](Step step, Time /*unused*/) { energyReductionStep_ = step; });
    }
    return nullptr;
}

//! Explicit template instantiation
//! @{
template class ComputeGlobalsElement<ComputeGlobalsAlgorithm::LeapFrog>;
template class ComputeGlobalsElement<ComputeGlobalsAlgorithm::VelocityVerlet>;
//! @}
} // namespace gmx
