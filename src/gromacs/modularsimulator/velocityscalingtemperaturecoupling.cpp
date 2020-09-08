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

/*! \internal
 * \brief Data used by the concrete temperature coupling implementations
 */
struct TemperatureCouplingData
{
    //! The coupling time step - simulation time step x nstcouple_
    const double couplingTimeStep;
    //! Coupling temperature per group
    ArrayRef<const real> referenceTemperature;
    //! Coupling time per group
    ArrayRef<const real> couplingTime;
    //! Number of degrees of freedom per group
    ArrayRef<const real> numDegreesOfFreedom;
    //! Work exerted by thermostat per group
    ArrayRef<const double> temperatureCouplingIntegral;
};

/*! \internal
 * \brief Interface for temperature coupling implementations
 */
class ITemperatureCouplingImpl
{
public:
    //! Allow access to the scaling vectors
    virtual void connectWithPropagator(const PropagatorThermostatConnection& connectionData,
                                       int numTemperatureGroups) = 0;

    /*! \brief Make a temperature control step
     *
     * \param step                     The current step
     * \param temperatureGroup         The current temperature group
     * \param currentKineticEnergy     The kinetic energy of the temperature group
     * \param currentTemperature       The temperature of the temperature group
     * \param temperatureCouplingData  Access to general temperature coupling data
     *
     * \return  The temperature coupling integral for the current temperature group
     */
    [[nodiscard]] virtual real apply(Step                           step,
                                     int                            temperatureGroup,
                                     real                           currentKineticEnergy,
                                     real                           currentTemperature,
                                     const TemperatureCouplingData& temperatureCouplingData) = 0;

    //! Write private data to checkpoint
    virtual void writeCheckpoint(WriteCheckpointData checkpointData, const t_commrec* cr) = 0;
    //! Read private data from checkpoint
    virtual void readCheckpoint(ReadCheckpointData checkpointData, const t_commrec* cr) = 0;

    //! Standard virtual destructor
    virtual ~ITemperatureCouplingImpl() = default;
};

/*! \internal
 * \brief Implements v-rescale temperature coupling
 */
class VRescaleTemperatureCoupling final : public ITemperatureCouplingImpl
{
public:
    //! Apply the v-rescale temperature control
    real apply(Step step,
               int  temperatureGroup,
               real currentKineticEnergy,
               real gmx_unused                currentTemperature,
               const TemperatureCouplingData& temperatureCouplingData) override
    {
        if (!(temperatureCouplingData.couplingTime[temperatureGroup] >= 0
              && temperatureCouplingData.numDegreesOfFreedom[temperatureGroup] > 0
              && currentKineticEnergy > 0))
        {
            lambdaStartVelocities_[temperatureGroup] = 1.0;
            return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup];
        }

        const real referenceKineticEnergy =
                0.5 * temperatureCouplingData.referenceTemperature[temperatureGroup] * BOLTZ
                * temperatureCouplingData.numDegreesOfFreedom[temperatureGroup];

        const real newKineticEnergy =
                vrescale_resamplekin(currentKineticEnergy, referenceKineticEnergy,
                                     temperatureCouplingData.numDegreesOfFreedom[temperatureGroup],
                                     temperatureCouplingData.couplingTime[temperatureGroup]
                                             / temperatureCouplingData.couplingTimeStep,
                                     step, seed_);

        // Analytically newKineticEnergy >= 0, but we check for rounding errors
        if (newKineticEnergy <= 0)
        {
            lambdaStartVelocities_[temperatureGroup] = 0.0;
        }
        else
        {
            lambdaStartVelocities_[temperatureGroup] = std::sqrt(newKineticEnergy / currentKineticEnergy);
        }

        if (debug)
        {
            fprintf(debug, "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n", temperatureGroup,
                    referenceKineticEnergy, currentKineticEnergy, newKineticEnergy,
                    lambdaStartVelocities_[temperatureGroup]);
        }

        return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup]
               - (newKineticEnergy - currentKineticEnergy);
    }

    //! Connect with propagator - v-rescale only scales start step velocities
    void connectWithPropagator(const PropagatorThermostatConnection& connectionData,
                               int                                   numTemperatureGroups) override
    {
        connectionData.setNumVelocityScalingVariables(numTemperatureGroups);
        lambdaStartVelocities_ = connectionData.getViewOnVelocityScaling();
    }

    //! No data to write to checkpoint
    void writeCheckpoint(WriteCheckpointData gmx_unused checkpointData, const t_commrec gmx_unused* cr) override
    {
    }
    //! No data to read from checkpoints
    void readCheckpoint(ReadCheckpointData gmx_unused checkpointData, const t_commrec gmx_unused* cr) override
    {
    }

    //! Constructor
    VRescaleTemperatureCoupling(int64_t seed) : seed_(seed) {}

private:
    //! The random seed
    const int64_t seed_;

    //! View on the scaling factor of the propagator (pre-step velocities)
    ArrayRef<real> lambdaStartVelocities_;
};

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
        EnergyData*                       energyData,
        TemperatureCouplingType           couplingType) :
    nstcouple_(nstcouple),
    offset_(offset),
    useFullStepKE_(useFullStepKE),
    reportPreviousConservedEnergy_(reportPreviousConservedEnergy),
    numTemperatureGroups_(numTemperatureGroups),
    couplingTimeStep_(couplingTimeStep),
    referenceTemperature_(referenceTemperature, referenceTemperature + numTemperatureGroups),
    couplingTime_(couplingTime, couplingTime + numTemperatureGroups),
    numDegreesOfFreedom_(numDegreesOfFreedom, numDegreesOfFreedom + numTemperatureGroups),
    temperatureCouplingIntegral_(numTemperatureGroups, 0.0),
    energyData_(energyData)
{
    if (reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
    {
        temperatureCouplingIntegralPreviousStep_ = temperatureCouplingIntegral_;
    }
    energyData->setVelocityScalingTemperatureCoupling(this);
    if (couplingType == etcVRESCALE)
    {
        temperatureCouplingImpl_ = std::make_unique<VRescaleTemperatureCoupling>(seed);
    }
    else
    {
        throw NotImplementedError("Temperature coupling " + std::string(ETCOUPLTYPE(couplingType))
                                  + " is not implemented for modular simulator.");
    }
}

void VelocityScalingTemperatureCoupling::connectWithPropagator(const PropagatorThermostatConnection& connectionData)
{
    temperatureCouplingImpl_->connectWithPropagator(connectionData, numTemperatureGroups_);
    propagatorCallback_ = connectionData.getVelocityScalingCallback();
}

void VelocityScalingTemperatureCoupling::elementSetup()
{
    if (!propagatorCallback_)
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
        temperatureCouplingIntegralPreviousStep_ = temperatureCouplingIntegral_;
    }

    const auto*             ekind          = energyData_->ekindata();
    TemperatureCouplingData thermostatData = { couplingTimeStep_, referenceTemperature_, couplingTime_,
                                               numDegreesOfFreedom_, temperatureCouplingIntegral_ };

    for (int temperatureGroup = 0; (temperatureGroup < numTemperatureGroups_); temperatureGroup++)
    {
        const real currentKineticEnergy = useFullStepKE_ == UseFullStepKE::Yes
                                                  ? trace(ekind->tcstat[temperatureGroup].ekinf)
                                                  : trace(ekind->tcstat[temperatureGroup].ekinh);
        const real currentTemperature = useFullStepKE_ == UseFullStepKE::Yes
                                                ? ekind->tcstat[temperatureGroup].T
                                                : ekind->tcstat[temperatureGroup].Th;

        temperatureCouplingIntegral_[temperatureGroup] = temperatureCouplingImpl_->apply(
                step, temperatureGroup, currentKineticEnergy, currentTemperature, thermostatData);
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
                                 makeCheckpointArrayRef<operation>(temperatureCouplingIntegral_));
    }
    if (operation == CheckpointDataOperation::Read && DOMAINDECOMP(cr))
    {
        dd_bcast(cr->dd, temperatureCouplingIntegral_.size() * sizeof(double),
                 temperatureCouplingIntegral_.data());
    }
}

void VelocityScalingTemperatureCoupling::writeCheckpoint(WriteCheckpointData checkpointData,
                                                         const t_commrec*    cr)
{
    doCheckpointData<CheckpointDataOperation::Write>(&checkpointData, cr);
    temperatureCouplingImpl_->writeCheckpoint(checkpointData.subCheckpointData("thermostat impl"), cr);
}

void VelocityScalingTemperatureCoupling::readCheckpoint(ReadCheckpointData checkpointData,
                                                        const t_commrec*   cr)
{
    doCheckpointData<CheckpointDataOperation::Read>(&checkpointData, cr);
    temperatureCouplingImpl_->readCheckpoint(checkpointData.subCheckpointData("thermostat impl"), cr);
}

const std::string& VelocityScalingTemperatureCoupling::clientID()
{
    return identifier_;
}

real VelocityScalingTemperatureCoupling::conservedEnergyContribution() const
{
    return (reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
                   ? std::accumulate(temperatureCouplingIntegralPreviousStep_.begin(),
                                     temperatureCouplingIntegralPreviousStep_.end(), 0.0)
                   : std::accumulate(temperatureCouplingIntegral_.begin(),
                                     temperatureCouplingIntegral_.end(), 0.0);
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
            legacySimulatorData->inputrec->opts.nrdf, energyData, legacySimulatorData->inputrec->etc));
    auto* thermostat = static_cast<VelocityScalingTemperatureCoupling*>(element);
    builderHelper->registerThermostat([thermostat](const PropagatorThermostatConnection& connection) {
        thermostat->connectWithPropagator(connection);
    });
    return element;
}

} // namespace gmx
