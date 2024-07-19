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
 * \brief Defines classes related to MTTK pressure coupling
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "mttk.h"

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/ifunc.h"

#include "energydata.h"
#include "nosehooverchains.h"
#include "simulatoralgorithm.h"
#include "trotterhelperfunctions.h"
#include "velocityscalingtemperaturecoupling.h"

namespace gmx
{

void MttkData::build(LegacySimulatorData*                    legacySimulatorData,
                     ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                     StatePropagatorData*                    statePropagatorData,
                     EnergyData*                             energyData,
                     const MttkPropagatorConnectionDetails&  mttkPropagatorConnectionDetails)
{
    // Uses reference temperature of first T-group
    const real ensembleTemperature = constantEnsembleTemperature(*legacySimulatorData->inputRec_);
    const real referencePressure =
            ::trace(legacySimulatorData->inputRec_->pressureCouplingOptions.ref_p) / DIM;
    // Weights are set based on initial volume
    real initialVolume = det(statePropagatorData->constBox());

    // When using domain decomposition, statePropagatorData might not have the initial
    // box yet, so we get it from the legacy state_global instead.
    // TODO: Make sure we have a valid state in statePropagatorData at all times (#3421)
    if (haveDDAtomOrdering(*legacySimulatorData->cr_))
    {
        if (MAIN(legacySimulatorData->cr_))
        {
            initialVolume = det(legacySimulatorData->stateGlobal_->box);
        }
        dd_bcast(legacySimulatorData->cr_->dd, int(sizeof(real)), &initialVolume);
    }

    GMX_RELEASE_ASSERT(
            !builderHelper->simulationData<MttkPropagatorConnection>(MttkPropagatorConnection::dataID()),
            "Attempted to build MttkPropagatorConnection more than once.");
    MttkPropagatorConnection::build(builderHelper,
                                    mttkPropagatorConnectionDetails.propagatorTagPrePosition,
                                    mttkPropagatorConnectionDetails.propagatorTagPostPosition,
                                    mttkPropagatorConnectionDetails.positionOffset,
                                    mttkPropagatorConnectionDetails.propagatorTagPreVelocity1,
                                    mttkPropagatorConnectionDetails.propagatorTagPostVelocity1,
                                    mttkPropagatorConnectionDetails.velocityOffset1,
                                    mttkPropagatorConnectionDetails.propagatorTagPreVelocity2,
                                    mttkPropagatorConnectionDetails.propagatorTagPostVelocity2,
                                    mttkPropagatorConnectionDetails.velocityOffset2);
    auto* mttkPropagatorConnection =
            builderHelper
                    ->simulationData<MttkPropagatorConnection>(MttkPropagatorConnection::dataID())
                    .value();

    builderHelper->storeSimulationData(
            MttkData::dataID(),
            MttkData(ensembleTemperature,
                     referencePressure,
                     legacySimulatorData->inputRec_->pressureCouplingOptions.nstpcouple
                             * legacySimulatorData->inputRec_->delta_t,
                     legacySimulatorData->inputRec_->pressureCouplingOptions.tau_p,
                     initialVolume,
                     legacySimulatorData->inputRec_->opts.nrdf[0],
                     legacySimulatorData->inputRec_->delta_t,
                     legacySimulatorData->inputRec_->pressureCouplingOptions.compress,
                     statePropagatorData,
                     mttkPropagatorConnection));
    auto* ptrToDataObject = builderHelper->simulationData<MttkData>(MttkData::dataID()).value();

    energyData->addConservedEnergyContribution([ptrToDataObject](Step /*unused*/, Time time) {
        return ptrToDataObject->temperatureCouplingIntegral(time);
    });
    energyData->setParrinelloRahmanBoxVelocities(
            [ptrToDataObject]() { return ptrToDataObject->boxVelocity_; });
    builderHelper->registerReferenceTemperatureUpdate(
            [ptrToDataObject](ArrayRef<const real> temperatures, ReferenceTemperatureChangeAlgorithm algorithm) {
                ptrToDataObject->updateReferenceTemperature(temperatures[0], algorithm);
            });
}

std::string MttkData::dataID()
{
    return "MttkData";
}

MttkData::MttkData(real                       ensembleTemperature,
                   real                       referencePressure,
                   real                       couplingTimeStep,
                   real                       couplingTime,
                   real                       initialVolume,
                   real                       numDegreesOfFreedom,
                   real                       simulationTimeStep,
                   const tensor               compressibility,
                   const StatePropagatorData* statePropagatorData,
                   MttkPropagatorConnection*  mttkPropagatorConnection) :
    couplingTimeStep_(couplingTimeStep),
    etaVelocity_(0.0),
    invMass_((c_presfac * ::trace(compressibility) * c_boltz * ensembleTemperature)
             / (DIM * initialVolume * gmx::square(couplingTime / M_2PI))),
    etaVelocityTime_(0.0),
    temperatureCouplingIntegral_(0.0),
    integralTime_(0.0),
    referencePressure_(referencePressure),
    boxVelocity_{ { 0 } },
    numDegreesOfFreedom_(numDegreesOfFreedom),
    simulationTimeStep_(simulationTimeStep),
    ensembleTemperature_(ensembleTemperature),
    statePropagatorData_(statePropagatorData),
    mttkPropagatorConnection_(mttkPropagatorConnection)
{
    // Set integral based on initial volume
    calculateIntegral(initialVolume);
}

MttkData::MttkData(const MttkData& other) :
    couplingTimeStep_(other.couplingTimeStep_),
    etaVelocity_(other.etaVelocity_),
    invMass_(other.invMass_),
    etaVelocityTime_(other.etaVelocityTime_),
    temperatureCouplingIntegral_(other.temperatureCouplingIntegral_),
    integralTime_(other.integralTime_),
    referencePressure_(other.referencePressure_),
    numDegreesOfFreedom_(other.numDegreesOfFreedom_),
    simulationTimeStep_(other.simulationTimeStep_),
    statePropagatorData_(other.statePropagatorData_),
    mttkPropagatorConnection_(other.mttkPropagatorConnection_)
{
    copy_mat(other.boxVelocity_, boxVelocity_);
}

void MttkData::calculateIntegralIfNeeded()
{
    // Check whether coordinate time divided by the time step is close to integer
    const bool calculationNeeded = timesClose(
            std::lround(etaVelocityTime_ / couplingTimeStep_) * couplingTimeStep_, etaVelocityTime_);

    if (calculationNeeded)
    {
        const real volume = det(statePropagatorData_->constBox());
        // Calculate current value of barostat integral
        calculateIntegral(volume);
    }
}

void MttkData::calculateIntegral(real volume)
{
    temperatureCouplingIntegral_ = kineticEnergy() + volume * referencePressure_ / c_presfac;
    integralTime_                = etaVelocityTime_;
}

real MttkData::kineticEnergy() const
{
    return 0.5 * etaVelocity_ * etaVelocity_ / invMass_;
}

void MttkData::scale(real scalingFactor, bool scalingAtFullCouplingTimeStep)
{
    etaVelocity_ *= scalingFactor;
    if (scalingAtFullCouplingTimeStep)
    {
        calculateIntegralIfNeeded();
    }
    updateScalingFactors();
}

real MttkData::etaVelocity() const
{
    return etaVelocity_;
}

real MttkData::invEtaMass() const
{
    return invMass_;
}

void MttkData::setEtaVelocity(real etaVelocity, real etaVelocityTimeIncrement)
{
    etaVelocity_ = etaVelocity;
    etaVelocityTime_ += etaVelocityTimeIncrement;
    calculateIntegralIfNeeded();
    updateScalingFactors();
}

double MttkData::temperatureCouplingIntegral(Time gmx_used_in_debug time) const
{
    /* When using nstpcouple >= nstcalcenergy, we accept that the coupling
     * integral might be ahead of the current energy calculation step. The
     * extended system degrees of freedom are either in sync or ahead of the
     * rest of the system.
     */
    GMX_ASSERT(time <= integralTime_ || timesClose(integralTime_, time),
               "MttkData conserved energy time mismatch.");
    return temperatureCouplingIntegral_;
}

real MttkData::referencePressure() const
{
    return referencePressure_;
}

rvec* MttkData::boxVelocities()
{
    return boxVelocity_;
}

void MttkData::updateReferenceTemperature(real temperature,
                                          ReferenceTemperatureChangeAlgorithm gmx_unused algorithm)
{
    // Currently, we don't know about any temperature change algorithms, so we assert this never gets called
    GMX_ASSERT(false, "MttkData: Unknown ReferenceTemperatureChangeAlgorithm.");
    invMass_ *= temperature / ensembleTemperature_;
    ensembleTemperature_ = temperature;
}

namespace
{
/*!
 * \brief Enum describing the contents MttkData writes to modular checkpoint
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
void MttkData::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "MttkData version", c_currentVersion);
    checkpointData->scalar("veta", &etaVelocity_);
    // Mass is calculated from initial volume, so need to save it for exact continuation
    checkpointData->scalar("mass", &invMass_);
    checkpointData->scalar("time", &etaVelocityTime_);
    checkpointData->scalar("integral", &temperatureCouplingIntegral_);
    checkpointData->scalar("integralTime", &integralTime_);
}

void MttkData::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void MttkData::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_bcast(cr->dd, int(sizeof(real)), &etaVelocity_);
        dd_bcast(cr->dd, int(sizeof(real)), &invMass_);
        dd_bcast(cr->dd, int(sizeof(Time)), &etaVelocityTime_);
        dd_bcast(cr->dd, int(sizeof(double)), &temperatureCouplingIntegral_);
        dd_bcast(cr->dd, int(sizeof(Time)), &integralTime_);
    }
}

const std::string& MttkData::clientID()
{
    return identifier_;
}

void MttkData::propagatorCallback(Step step) const
{
    mttkPropagatorConnection_->propagatorCallback(step);
}

void MttkPropagatorConnection::build(ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                     const PropagatorTag& propagatorTagPrePosition,
                                     const PropagatorTag& propagatorTagPostPosition,
                                     int                  positionOffset,
                                     const PropagatorTag& propagatorTagPreVelocity1,
                                     const PropagatorTag& propagatorTagPostVelocity1,
                                     int                  velocityOffset1,
                                     const PropagatorTag& propagatorTagPreVelocity2,
                                     const PropagatorTag& propagatorTagPostVelocity2,
                                     int                  velocityOffset2)
{
    GMX_RELEASE_ASSERT(!(propagatorTagPrePosition == propagatorTagPostPosition
                         && propagatorTagPrePosition != PropagatorTag("")),
                       "Pre- and post-step position scaling in same element is not supported.");
    GMX_RELEASE_ASSERT(!((propagatorTagPreVelocity1 == propagatorTagPostVelocity1
                          && propagatorTagPreVelocity1 != PropagatorTag(""))
                         || (propagatorTagPreVelocity2 == propagatorTagPostVelocity2
                             && propagatorTagPreVelocity2 != PropagatorTag(""))),
                       "Pre- and post-step velocity scaling in same element is not implemented.");

    // Store object with simulation algorithm for safe pointer capturing
    builderHelper->storeSimulationData(MttkPropagatorConnection::dataID(), MttkPropagatorConnection());
    auto* object = builderHelper
                           ->simulationData<MttkPropagatorConnection>(MttkPropagatorConnection::dataID())
                           .value();

    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPrePosition, positionOffset](const PropagatorConnection& connection) {
                object->connectWithPropagatorPositionPreStepScaling(
                        connection, propagatorTagPrePosition, positionOffset);
            });
    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPostPosition, positionOffset](const PropagatorConnection& connection) {
                object->connectWithPropagatorPositionPostStepScaling(
                        connection, propagatorTagPostPosition, positionOffset);
            });
    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPreVelocity1, velocityOffset1](const PropagatorConnection& connection) {
                object->connectWithPropagatorVelocityPreStepScaling(
                        connection, propagatorTagPreVelocity1, velocityOffset1);
            });
    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPostVelocity1, velocityOffset1](const PropagatorConnection& connection) {
                object->connectWithPropagatorVelocityPostStepScaling(
                        connection, propagatorTagPostVelocity1, velocityOffset1);
            });
    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPreVelocity2, velocityOffset2](const PropagatorConnection& connection) {
                object->connectWithPropagatorVelocityPreStepScaling(
                        connection, propagatorTagPreVelocity2, velocityOffset2);
            });
    builderHelper->registerTemperaturePressureControl(
            [object, propagatorTagPostVelocity2, velocityOffset2](const PropagatorConnection& connection) {
                object->connectWithPropagatorVelocityPostStepScaling(
                        connection, propagatorTagPostVelocity2, velocityOffset2);
            });
}

void MttkPropagatorConnection::propagatorCallback(Step step) const
{
    for (const auto& callback : propagatorCallbacks_)
    {
        std::get<0>(callback)(step + std::get<1>(callback));
    }
}

void MttkPropagatorConnection::setPositionScaling(real preStepScaling, real postStepScaling)
{
    for (const auto& scalingFactor : startPositionScalingFactors_)
    {
        std::fill(scalingFactor.begin(), scalingFactor.end(), preStepScaling);
    }
    for (const auto& scalingFactor : endPositionScalingFactors_)
    {
        std::fill(scalingFactor.begin(), scalingFactor.end(), postStepScaling);
    }
}

void MttkPropagatorConnection::setVelocityScaling(real preStepScaling, real postStepScaling)
{
    for (const auto& scalingFactor : startVelocityScalingFactors_)
    {
        std::fill(scalingFactor.begin(), scalingFactor.end(), preStepScaling);
    }
    for (const auto& scalingFactor : endVelocityScalingFactors_)
    {
        std::fill(scalingFactor.begin(), scalingFactor.end(), postStepScaling);
    }
}

std::string MttkPropagatorConnection::dataID()
{
    return "MttkPropagatorConnection";
}

void MttkPropagatorConnection::connectWithPropagatorVelocityPreStepScaling(const PropagatorConnection& connectionData,
                                                                           const PropagatorTag& propagatorTag,
                                                                           int offset)
{
    if (connectionData.tag == propagatorTag && connectionData.hasStartVelocityScaling())
    {
        connectionData.setNumVelocityScalingVariables(1, ScaleVelocities::PreStepOnly);
        startVelocityScalingFactors_.emplace_back(connectionData.getViewOnStartVelocityScaling());
        propagatorCallbacks_.emplace_back(connectionData.getVelocityScalingCallback(), offset);
    }
}

void MttkPropagatorConnection::connectWithPropagatorVelocityPostStepScaling(const PropagatorConnection& connectionData,
                                                                            const PropagatorTag& propagatorTag,
                                                                            int offset)
{
    if (connectionData.tag == propagatorTag && connectionData.hasStartVelocityScaling())
    {
        // Although we're using this propagator for scaling after the update, we're using
        // getViewOnStartVelocityScaling() - getViewOnEndVelocityScaling() is only
        // used for propagators doing BOTH start and end scaling
        connectionData.setNumVelocityScalingVariables(1, ScaleVelocities::PreStepOnly);
        endVelocityScalingFactors_.emplace_back(connectionData.getViewOnStartVelocityScaling());
        propagatorCallbacks_.emplace_back(connectionData.getVelocityScalingCallback(), offset);
    }
}

void MttkPropagatorConnection::connectWithPropagatorPositionPreStepScaling(const PropagatorConnection& connectionData,
                                                                           const PropagatorTag& propagatorTag,
                                                                           int offset)
{
    if (connectionData.tag == propagatorTag && connectionData.hasPositionScaling())
    {
        connectionData.setNumPositionScalingVariables(1);
        startPositionScalingFactors_.emplace_back(connectionData.getViewOnPositionScaling());
        propagatorCallbacks_.emplace_back(connectionData.getPositionScalingCallback(), offset);
    }
}

void MttkPropagatorConnection::connectWithPropagatorPositionPostStepScaling(const PropagatorConnection& connectionData,
                                                                            const PropagatorTag& propagatorTag,
                                                                            int offset)
{
    if (connectionData.tag == propagatorTag && connectionData.hasPositionScaling())
    {
        connectionData.setNumPositionScalingVariables(1);
        endPositionScalingFactors_.emplace_back(connectionData.getViewOnPositionScaling());
        propagatorCallbacks_.emplace_back(connectionData.getPositionScalingCallback(), offset);
    }
}

void MttkData::updateScalingFactors()
{
    // Tuckerman et al. 2006, Eq 5.8
    // Note that we're using the dof of the first temperature group only
    const real alpha = 1.0 + DIM / (numDegreesOfFreedom_);
    /* Tuckerman et al. 2006, eqs 5.11 and 5.13:
     *
     * r(t+dt)   = r(t)*exp(v_eta*dt) + dt*v*exp(v_eta*dt/2) * [sinh(v_eta*dt/2) / (v_eta*dt/2)]
     * v(t+dt/2) = v(t)*exp(-a*v_eta*dt/2) +
     *             dt/2*f/m*exp(-a*v_eta*dt/4) * [sinh(a*v_eta*dt/4) / (a*v_eta*dt/4)]
     * with a = 1 + 1/Natoms
     *
     * For r, let
     *   s1 = exp(v_eta*dt/2)
     *   s2 = [sinh(v_eta*dt/2) / (v_eta*dt/2)]
     * so we can use
     *   r(t) *= s1/s2
     *   r(t+dt) = r(t) + dt*v
     *   r(t+dt) *= s1*s2  <=>  r(t+dt) = s1*s2 * (r(t)*s1/s2 + dt*v) = s1^2*r(t) + dt*v*s1*s2
     *
     * For v, let
     *   s1 = exp(-a*v_eta*dt/4)
     *   s2 = [sinh(a*v_eta*dt/4) / (a*v_eta*dt/4)]
     * so we can use
     *   v(t) *= s1/s2
     *   v(t+dt/2) = v(t) + dt/2*f/m
     *   v(t+dt/2) *= s1*s2  <=>  v(t+dt/2) = s1^2*v(t) + dt/2*f/m*s1*s2
     *
     * In legacy simulator, this scaling is applied every step, even if the barostat is updated
     * less frequently, so we are mirroring this by using the simulation time step for dt and
     * requesting scaling every step. This could likely be applied impulse-style by using the
     * coupling time step for dt and only applying it when the barostat gets updated.
     */
    const real scalingPos1 = std::exp(0.5 * simulationTimeStep_ * etaVelocity_);
    const real scalingPos2 = gmx::series_sinhx(0.5 * simulationTimeStep_ * etaVelocity_);
    const real scalingVel1 = std::exp(-alpha * 0.25 * simulationTimeStep_ * etaVelocity_);
    const real scalingVel2 = gmx::series_sinhx(alpha * 0.25 * simulationTimeStep_ * etaVelocity_);

    mttkPropagatorConnection_->setPositionScaling(scalingPos1 / scalingPos2, scalingPos1 * scalingPos2);
    mttkPropagatorConnection_->setVelocityScaling(scalingVel1 / scalingVel2, scalingVel1 * scalingVel2);
}

void MttkElement::propagateEtaVelocity(Step step)
{
    const auto* ekind         = energyData_->ekindata();
    const auto* virial        = energyData_->totalVirial(step);
    const real  currentVolume = det(statePropagatorData_->constBox());
    // Tuckerman et al. 2006, Eq 5.8
    // Note that we're using the dof of the first temperature group only
    const real alpha = 1.0 + DIM / (numDegreesOfFreedom_);
    // Also here, using first group only
    const real kineticEnergyFactor = alpha * ekind->tcstat[0].ekinscalef_nhc;
    // Now, we're using full system kinetic energy!
    tensor modifiedKineticEnergy;
    msmul(ekind->ekin, kineticEnergyFactor, modifiedKineticEnergy);

    tensor currentPressureTensor;

    const real currentPressure = calc_pres(
            pbcType_, numWalls_, statePropagatorData_->constBox(), modifiedKineticEnergy, virial, currentPressureTensor);

    const real etaAcceleration = DIM * currentVolume * (mttkData_->invEtaMass() / c_presfac)
                                 * (currentPressure - mttkData_->referencePressure());

    mttkData_->setEtaVelocity(mttkData_->etaVelocity() + propagationTimeStep_ * etaAcceleration,
                              propagationTimeStep_);
}

MttkElement::MttkElement(int                        nstcouple,
                         int                        offset,
                         real                       propagationTimeStep,
                         ScheduleOnInitStep         scheduleOnInitStep,
                         Step                       initStep,
                         const StatePropagatorData* statePropagatorData,
                         EnergyData*                energyData,
                         MttkData*                  mttkData,
                         PbcType                    pbcType,
                         int                        numWalls,
                         real                       numDegreesOfFreedom) :
    pbcType_(pbcType),
    numWalls_(numWalls),
    numDegreesOfFreedom_(numDegreesOfFreedom),
    nstcouple_(nstcouple),
    offset_(offset),
    propagationTimeStep_(propagationTimeStep),
    scheduleOnInitStep_(scheduleOnInitStep),
    initialStep_(initStep),
    statePropagatorData_(statePropagatorData),
    energyData_(energyData),
    mttkData_(mttkData)
{
}

void MttkElement::scheduleTask(Step step, Time /*unused*/, const RegisterRunFunction& registerRunFunction)
{
    if (step == initialStep_ && scheduleOnInitStep_ == ScheduleOnInitStep::No)
    {
        return;
    }
    if (do_per_step(step + nstcouple_ + offset_, nstcouple_))
    {
        // do T-coupling this step
        registerRunFunction([this, step]() { propagateEtaVelocity(step); });
    }

    // Let propagators know that we want to scale
    // (we're scaling every step - see comment in MttkData::updateScalingFactors())
    mttkData_->propagatorCallback(step);
}

ISimulatorElement* MttkElement::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer gmx_unused*         observablesReducer,
        Offset                                 offset,
        ScheduleOnInitStep                     scheduleOnInitStep,
        const MttkPropagatorConnectionDetails& mttkPropagatorConnectionDetails)
{
    // Data is now owned by the caller of this method, who will handle lifetime (see ModularSimulatorAlgorithm)
    if (!builderHelper->simulationData<MttkData>(MttkData::dataID()))
    {
        MttkData::build(legacySimulatorData, builderHelper, statePropagatorData, energyData, mttkPropagatorConnectionDetails);
    }
    auto* mttkData = builderHelper->simulationData<MttkData>(MttkData::dataID()).value();

    // Element is now owned by the caller of this method, who will handle lifetime (see ModularSimulatorAlgorithm)
    auto* element = static_cast<MttkElement*>(builderHelper->storeElement(std::make_unique<MttkElement>(
            legacySimulatorData->inputRec_->nsttcouple,
            offset,
            legacySimulatorData->inputRec_->delta_t
                    * legacySimulatorData->inputRec_->pressureCouplingOptions.nstpcouple / 2,
            scheduleOnInitStep,
            legacySimulatorData->inputRec_->init_step,
            statePropagatorData,
            energyData,
            mttkData,
            legacySimulatorData->inputRec_->pbcType,
            legacySimulatorData->inputRec_->nwall,
            legacySimulatorData->inputRec_->opts.nrdf[0])));

    return element;
}

MttkBoxScaling::MttkBoxScaling(real                 simulationTimeStep,
                               StatePropagatorData* statePropagatorData,
                               MttkData*            mttkData) :
    simulationTimeStep_(simulationTimeStep), statePropagatorData_(statePropagatorData), mttkData_(mttkData)
{
}

void MttkBoxScaling::scheduleTask(Step gmx_unused            step,
                                  gmx_unused Time            time,
                                  const RegisterRunFunction& registerRunFunction)
{
    registerRunFunction([this]() { scaleBox(); });
}

void MttkBoxScaling::scaleBox()
{
    auto* box = statePropagatorData_->box();

    /* DIM * eta = ln V.  so DIM*eta_new = DIM*eta_old + DIM*dt*veta =>
       ln V_new = ln V_old + 3*dt*veta => V_new = V_old*exp(3*dt*veta) =>
       Side length scales as exp(veta*dt) */
    msmul(box, std::exp(mttkData_->etaVelocity() * simulationTimeStep_), box);

    /* Relate veta to boxv.  veta = d(eta)/dT = (1/DIM)*1/V dV/dT.
       o               If we assume isotropic scaling, and box length scaling
       factor L, then V = L^DIM (det(M)).  So dV/dt = DIM
       L^(DIM-1) dL/dt det(M), and veta = (1/L) dL/dt.  The
       determinant of B is L^DIM det(M), and the determinant
       of dB/dt is (dL/dT)^DIM det (M).  veta will be
       (det(dB/dT)/det(B))^(1/3).  Then since M =
       B_new*(vol_new)^(1/3), dB/dT_new = (veta_new)*B(new). */
    msmul(box, mttkData_->etaVelocity(), mttkData_->boxVelocities());

    mttkData_->calculateIntegralIfNeeded();
}

ISimulatorElement* MttkBoxScaling::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData*                    statePropagatorData,
        EnergyData*                             energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer gmx_unused*         observablesReducer,
        const MttkPropagatorConnectionDetails& mttkPropagatorConnectionDetails)
{
    // Data is now owned by the caller of this method, who will handle lifetime (see ModularSimulatorAlgorithm)
    if (!builderHelper->simulationData<MttkData>(MttkData::dataID()))
    {
        MttkData::build(legacySimulatorData, builderHelper, statePropagatorData, energyData, mttkPropagatorConnectionDetails);
    }

    return builderHelper->storeElement(std::make_unique<MttkBoxScaling>(
            legacySimulatorData->inputRec_->delta_t,
            statePropagatorData,
            builderHelper->simulationData<MttkData>(MttkData::dataID()).value()));
}

} // namespace gmx
