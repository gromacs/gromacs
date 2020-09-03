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
 * \brief Defines a velocity-scaling temperature coupling element for
 * the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "velocityscalingtemperaturecoupling.h"

#include <numeric>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"

#include "modularsimulator.h"
#include "simulatoralgorithm.h"

namespace gmx
{

VelocityScalingTemperatureCoupling::VelocityScalingTemperatureCoupling(
        int                               nstcouple,
        int                               offset,
        UseFullStepKE                     useFullStepKE,
        ReportPreviousStepConservedEnergy reportPreviousConservedEnergy,
        int64_t                           seed,
        int                               numTemperatureGroups,
        double                            couplingTimeStep,
        const real*                       referenceTemperature,
        const real*                       couplingTime,
        const real*                       numDegreesOfFreedom,
        EnergyData*                       energyData) :
    nstcouple_(nstcouple),
    offset_(offset),
    useFullStepKE_(useFullStepKE),
    reportPreviousConservedEnergy_(reportPreviousConservedEnergy),
    seed_(seed),
    numTemperatureGroups_(numTemperatureGroups),
    couplingTimeStep_(couplingTimeStep),
    referenceTemperature_(referenceTemperature, referenceTemperature + numTemperatureGroups),
    couplingTime_(couplingTime, couplingTime + numTemperatureGroups),
    numDegreesOfFreedom_(numDegreesOfFreedom, numDegreesOfFreedom + numTemperatureGroups),
    thermostatIntegral_(numTemperatureGroups, 0.0),
    energyData_(energyData)
{
    if (reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
    {
        thermostatIntegralPreviousStep_ = thermostatIntegral_;
    }
    energyData->setVelocityScalingTemperatureCoupling(this);
}

void VelocityScalingTemperatureCoupling::connectWithPropagator(const PropagatorThermostatConnection& connectionData)
{
    connectionData.setNumVelocityScalingVariables(numTemperatureGroups_);
    lambda_             = connectionData.getViewOnVelocityScaling();
    propagatorCallback_ = connectionData.getVelocityScalingCallback();
}

void VelocityScalingTemperatureCoupling::elementSetup()
{
    if (!propagatorCallback_ || lambda_.empty())
    {
        throw MissingElementConnectionError(
                "Velocity scaling temperature coupling was not connected to a propagator.\n"
                "Connection to a propagator element is needed to scale the velocities.\n"
                "Use connectWithPropagator(...) before building the ModularSimulatorAlgorithm "
                "object.");
    }
}

void VelocityScalingTemperatureCoupling::scheduleTask(Step step,
                                                      Time gmx_unused            time,
                                                      const RegisterRunFunction& registerRunFunction)
{
    /* The thermostat will need a valid kinetic energy when it is running.
     * Currently, computeGlobalCommunicationPeriod() is making sure this
     * happens on time.
     * TODO: Once we're switching to a new global communication scheme, we
     *       will want the thermostat to signal that global reduction
     *       of the kinetic energy is needed.
     *
     */
    if (do_per_step(step + nstcouple_ + offset_, nstcouple_))
    {
        // do T-coupling this step
        registerRunFunction([this, step]() { setLambda(step); });

        // Let propagator know that we want to do T-coupling
        propagatorCallback_(step);
    }
}

void VelocityScalingTemperatureCoupling::setLambda(Step step)
{
    // if we report the previous energy, calculate before the step
    if (reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
    {
        thermostatIntegralPreviousStep_ = thermostatIntegral_;
    }

    real currentKineticEnergy, referenceKineticEnergy, newKineticEnergy;

    const auto* ekind = energyData_->ekindata();

    for (int i = 0; (i < numTemperatureGroups_); i++)
    {
        if (useFullStepKE_ == UseFullStepKE::Yes)
        {
            currentKineticEnergy = trace(ekind->tcstat[i].ekinf);
        }
        else
        {
            currentKineticEnergy = trace(ekind->tcstat[i].ekinh);
        }

        if (couplingTime_[i] >= 0 && numDegreesOfFreedom_[i] > 0 && currentKineticEnergy > 0)
        {
            referenceKineticEnergy = 0.5 * referenceTemperature_[i] * BOLTZ * numDegreesOfFreedom_[i];

            newKineticEnergy = vrescale_resamplekin(currentKineticEnergy, referenceKineticEnergy,
                                                    numDegreesOfFreedom_[i],
                                                    couplingTime_[i] / couplingTimeStep_, step, seed_);

            // Analytically newKineticEnergy >= 0, but we check for rounding errors
            if (newKineticEnergy <= 0)
            {
                lambda_[i] = 0.0;
            }
            else
            {
                lambda_[i] = std::sqrt(newKineticEnergy / currentKineticEnergy);
            }

            thermostatIntegral_[i] -= newKineticEnergy - currentKineticEnergy;

            if (debug)
            {
                fprintf(debug, "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n", i,
                        referenceKineticEnergy, currentKineticEnergy, newKineticEnergy, lambda_[i]);
            }
        }
        else
        {
            lambda_[i] = 1.0;
        }
    }
}

namespace
{
/*!
 * \brief Enum describing the contents VelocityScalingTemperatureCoupling writes to modular checkpoint
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
void VelocityScalingTemperatureCoupling::doCheckpointData(CheckpointData<operation>* checkpointData,
                                                          const t_commrec*           cr)
{
    if (MASTER(cr))
    {
        checkpointVersion(checkpointData, "VRescaleThermostat version", c_currentVersion);

        checkpointData->arrayRef("thermostat integral",
                                 makeCheckpointArrayRef<operation>(thermostatIntegral_));
    }
    if (operation == CheckpointDataOperation::Read && DOMAINDECOMP(cr))
    {
        dd_bcast(cr->dd, thermostatIntegral_.size() * sizeof(double), thermostatIntegral_.data());
    }
}

void VelocityScalingTemperatureCoupling::writeCheckpoint(WriteCheckpointData checkpointData,
                                                         const t_commrec*    cr)
{
    doCheckpointData<CheckpointDataOperation::Write>(&checkpointData, cr);
}

void VelocityScalingTemperatureCoupling::readCheckpoint(ReadCheckpointData checkpointData,
                                                        const t_commrec*   cr)
{
    doCheckpointData<CheckpointDataOperation::Read>(&checkpointData, cr);
}

const std::string& VelocityScalingTemperatureCoupling::clientID()
{
    return identifier_;
}

real VelocityScalingTemperatureCoupling::conservedEnergyContribution() const
{
    return (reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
                   ? std::accumulate(thermostatIntegralPreviousStep_.begin(),
                                     thermostatIntegralPreviousStep_.end(), 0.0)
                   : std::accumulate(thermostatIntegral_.begin(), thermostatIntegral_.end(), 0.0);
}

ISimulatorElement* VelocityScalingTemperatureCoupling::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData gmx_unused*     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        int                                   offset,
        UseFullStepKE                         useFullStepKE,
        ReportPreviousStepConservedEnergy     reportPreviousStepConservedEnergy)
{
    auto* element = builderHelper->storeElement(std::make_unique<VelocityScalingTemperatureCoupling>(
            legacySimulatorData->inputrec->nsttcouple, offset, useFullStepKE, reportPreviousStepConservedEnergy,
            legacySimulatorData->inputrec->ld_seed, legacySimulatorData->inputrec->opts.ngtc,
            legacySimulatorData->inputrec->delta_t * legacySimulatorData->inputrec->nsttcouple,
            legacySimulatorData->inputrec->opts.ref_t, legacySimulatorData->inputrec->opts.tau_t,
            legacySimulatorData->inputrec->opts.nrdf, energyData));
    auto* thermostat = static_cast<VelocityScalingTemperatureCoupling*>(element);
    builderHelper->registerThermostat([thermostat](const PropagatorThermostatConnection& connection) {
        thermostat->connectWithPropagator(connection);
    });
    return element;
}

} // namespace gmx
