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
 * \brief Defines a velocity-scaling temperature coupling element for
 * the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "velocityscalingtemperaturecoupling.h"

#include <cmath>
#include <cstdio>

#include <algorithm>
#include <functional>
#include <numeric>
#include <type_traits>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdrun/isimulator.h"
#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/observablesreducer.h"
#include "gromacs/modularsimulator/energydata.h"
#include "gromacs/modularsimulator/modularsimulatorinterfaces.h"
#include "gromacs/modularsimulator/propagator.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strconvert.h"

#include "modularsimulator.h"
#include "simulatoralgorithm.h"

namespace gmx
{
class FreeEnergyPerturbationData;
class StatePropagatorData;
enum class ReferenceTemperatureChangeAlgorithm;
template<CheckpointDataOperation operation>
class CheckpointData;

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
    virtual void connectWithPropagator(const PropagatorConnection& connectionData,
                                       int                         numTemperatureGroups) = 0;

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
    virtual void writeCheckpoint(std::optional<WriteCheckpointData> checkpointData,
                                 const t_commrec*                   cr) = 0;
    //! Read private data from checkpoint
    virtual void readCheckpoint(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) = 0;

    //! Update the reference temperature and update and return the temperature coupling integral
    virtual real updateReferenceTemperatureAndIntegral(int  temperatureGroup,
                                                       real newTemperature,
                                                       ReferenceTemperatureChangeAlgorithm algorithm,
                                                       const TemperatureCouplingData& temperatureCouplingData) = 0;

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
    real apply(Step                           step,
               int                            temperatureGroup,
               real                           currentKineticEnergy,
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
                0.5 * temperatureCouplingData.referenceTemperature[temperatureGroup] * gmx::c_boltz
                * temperatureCouplingData.numDegreesOfFreedom[temperatureGroup];

        const real newKineticEnergy =
                vrescale_resamplekin(currentKineticEnergy,
                                     referenceKineticEnergy,
                                     temperatureCouplingData.numDegreesOfFreedom[temperatureGroup],
                                     temperatureCouplingData.couplingTime[temperatureGroup]
                                             / temperatureCouplingData.couplingTimeStep,
                                     step,
                                     seed_);

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
            fprintf(debug,
                    "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n",
                    temperatureGroup,
                    referenceKineticEnergy,
                    currentKineticEnergy,
                    newKineticEnergy,
                    lambdaStartVelocities_[temperatureGroup]);
        }

        return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup]
               - (newKineticEnergy - currentKineticEnergy);
    }

    //! Connect with propagator - v-rescale only scales start step velocities
    void connectWithPropagator(const PropagatorConnection& connectionData, int numTemperatureGroups) override
    {
        GMX_RELEASE_ASSERT(connectionData.hasStartVelocityScaling(),
                           "V-Rescale requires start velocity scaling.");
        connectionData.setNumVelocityScalingVariables(numTemperatureGroups, ScaleVelocities::PreStepOnly);
        lambdaStartVelocities_ = connectionData.getViewOnStartVelocityScaling();
    }

    //! No data to write to checkpoint
    void writeCheckpoint(std::optional<WriteCheckpointData> gmx_unused checkpointData,
                         const t_commrec gmx_unused* cr) override
    {
    }
    //! No data to read from checkpoints
    void readCheckpoint(std::optional<ReadCheckpointData> gmx_unused checkpointData,
                        const t_commrec gmx_unused* cr) override
    {
    }

    //! No changes needed
    real updateReferenceTemperatureAndIntegral(int             temperatureGroup,
                                               real gmx_unused newTemperature,
                                               ReferenceTemperatureChangeAlgorithm gmx_unused algorithm,
                                               const TemperatureCouplingData& temperatureCouplingData) override
    {
        return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup];
    }

    //! Constructor
    VRescaleTemperatureCoupling(int64_t seed) : seed_(seed) {}

private:
    //! The random seed
    const int64_t seed_;

    //! View on the scaling factor of the propagator (pre-step velocities)
    ArrayRef<real> lambdaStartVelocities_;
};

/*! \internal
 * \brief Implements Berendsen temperature coupling
 */
class BerendsenTemperatureCoupling final : public ITemperatureCouplingImpl
{
public:
    //! Apply the v-rescale temperature control
    real apply(Step gmx_unused                step,
               int                            temperatureGroup,
               real                           currentKineticEnergy,
               real                           currentTemperature,
               const TemperatureCouplingData& temperatureCouplingData) override
    {
        if (!(temperatureCouplingData.couplingTime[temperatureGroup] >= 0
              && temperatureCouplingData.numDegreesOfFreedom[temperatureGroup] > 0
              && currentKineticEnergy > 0))
        {
            lambdaStartVelocities_[temperatureGroup] = 1.0;
            return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup];
        }

        real lambda =
                std::sqrt(1.0
                          + (temperatureCouplingData.couplingTimeStep
                             / temperatureCouplingData.couplingTime[temperatureGroup])
                                    * (temperatureCouplingData.referenceTemperature[temperatureGroup] / currentTemperature
                                       - 1.0));
        lambdaStartVelocities_[temperatureGroup] =
                std::max<real>(std::min<real>(lambda, 1.25_real), 0.8_real);
        if (debug)
        {
            fprintf(debug,
                    "TC: group %d: T: %g, Lambda: %g\n",
                    temperatureGroup,
                    currentTemperature,
                    lambdaStartVelocities_[temperatureGroup]);
        }
        return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup]
               - (lambdaStartVelocities_[temperatureGroup] * lambdaStartVelocities_[temperatureGroup]
                  - 1) * currentKineticEnergy;
    }

    //! Connect with propagator - Berendsen only scales start step velocities
    void connectWithPropagator(const PropagatorConnection& connectionData, int numTemperatureGroups) override
    {
        GMX_RELEASE_ASSERT(connectionData.hasStartVelocityScaling(),
                           "Berendsen T-coupling requires start velocity scaling.");
        connectionData.setNumVelocityScalingVariables(numTemperatureGroups, ScaleVelocities::PreStepOnly);
        lambdaStartVelocities_ = connectionData.getViewOnStartVelocityScaling();
    }

    //! No data to write to checkpoint
    void writeCheckpoint(std::optional<WriteCheckpointData> gmx_unused checkpointData,
                         const t_commrec gmx_unused* cr) override
    {
    }
    //! No data to read from checkpoints
    void readCheckpoint(std::optional<ReadCheckpointData> gmx_unused checkpointData,
                        const t_commrec gmx_unused* cr) override
    {
    }

    //! No changes needed
    real updateReferenceTemperatureAndIntegral(int             temperatureGroup,
                                               real gmx_unused newTemperature,
                                               ReferenceTemperatureChangeAlgorithm gmx_unused algorithm,
                                               const TemperatureCouplingData& temperatureCouplingData) override
    {
        return temperatureCouplingData.temperatureCouplingIntegral[temperatureGroup];
    }

private:
    //! View on the scaling factor of the propagator (pre-step velocities)
    ArrayRef<real> lambdaStartVelocities_;
};

// Prepare NoseHooverTemperatureCoupling checkpoint data
namespace
{
/*!
 * \brief Enum describing the contents NoseHoover writes to modular checkpoint
 *
 * When changing the checkpoint content, add a new element just above Count, and adjust the
 * checkpoint functionality.
 */
enum class NHCheckpointVersion
{
    Base, //!< First version of modular checkpointing
    Count //!< Number of entries. Add new versions right above this!
};
constexpr auto c_nhCurrentVersion = NHCheckpointVersion(int(NHCheckpointVersion::Count) - 1);
} // namespace

/*! \internal
 * \brief Implements the Nose-Hoover temperature coupling
 */
class NoseHooverTemperatureCoupling final : public ITemperatureCouplingImpl
{
public:
    //! Calculate the current value of the temperature coupling integral
    real integral(int temperatureGroup, real numDegreesOfFreedom, real referenceTemperature)
    {
        return 0.5 * c_boltz * numDegreesOfFreedom
                       * (xiVelocities_[temperatureGroup] * xiVelocities_[temperatureGroup])
                       / invXiMass_[temperatureGroup]
               + numDegreesOfFreedom * xi_[temperatureGroup] * c_boltz * referenceTemperature;
    }

    //! Apply the Nose-Hoover temperature control
    real apply(Step gmx_unused                step,
               int                            temperatureGroup,
               real                           currentKineticEnergy,
               real                           currentTemperature,
               const TemperatureCouplingData& thermostatData) override
    {
        return applyLeapFrog(
                step, temperatureGroup, currentKineticEnergy, currentTemperature, thermostatData);
    }

    /*! \brief Apply for leap-frog
     *
     * This is called after the force calculation, before coordinate update
     *
     * We expect system to be at x(t), v(t-dt/2), f(t), T(t-dt/2)
     * Internal variables are at xi(t-dt), v_xi(t-dt)
     * Force on xi is calculated at time of system temperature
     * After calling this, we will have xi(t), v_xi(t)
     * The thermostat integral returned is a function of xi and v_xi,
     * and hence at time t.
     *
     * This performs an update of the thermostat variables calculated as
     *     a_xi(t-dt/2) = (T_sys(t-dt/2) - T_ref) / mass_xi;
     *     v_xi(t) = v_xi(t-dt) + dt_xi * a_xi(t-dt/2);
     *     xi(t) = xi(t-dt) + dt_xi * (v_xi(t-dt) + v_xi(t))/2;
     *
     * This will be followed by leap-frog integration of coordinates, calculated as
     *     v(t-dt/2) *= - 0.5 * dt * v_xi(t);  // scale previous velocities
     *     v(t+dt/2) = update_leapfrog_v(v(t-dt/2), f(t));  // do whatever LF does
     *     v(t+dt/2) *= 1 / (1 + 0.5 * dt * v_xi(t))  // scale new velocities
     *     x(t+dt) = update_leapfrog_x(x(t), v(t+dt/2));  // do whatever LF does
     */
    real applyLeapFrog(Step gmx_unused                step,
                       int                            temperatureGroup,
                       real                           currentKineticEnergy,
                       real                           currentTemperature,
                       const TemperatureCouplingData& thermostatData)
    {
        if (!(thermostatData.couplingTime[temperatureGroup] >= 0
              && thermostatData.numDegreesOfFreedom[temperatureGroup] > 0 && currentKineticEnergy > 0))
        {
            lambdaStartVelocities_[temperatureGroup] = 1.0;
            lambdaEndVelocities_[temperatureGroup]   = 1.0;
            return thermostatData.temperatureCouplingIntegral[temperatureGroup];
        }

        const auto oldXiVelocity = xiVelocities_[temperatureGroup];
        const auto xiAcceleration =
                invXiMass_[temperatureGroup]
                * (currentTemperature - thermostatData.referenceTemperature[temperatureGroup]);
        xiVelocities_[temperatureGroup] += thermostatData.couplingTimeStep * xiAcceleration;
        xi_[temperatureGroup] += thermostatData.couplingTimeStep
                                 * (oldXiVelocity + xiVelocities_[temperatureGroup]) * 0.5;
        lambdaStartVelocities_[temperatureGroup] =
                (1 - 0.5 * thermostatData.couplingTimeStep * xiVelocities_[temperatureGroup]);
        lambdaEndVelocities_[temperatureGroup] =
                1. / (1 + 0.5 * thermostatData.couplingTimeStep * xiVelocities_[temperatureGroup]);

        // Current value of the thermostat integral
        return integral(temperatureGroup,
                        thermostatData.numDegreesOfFreedom[temperatureGroup],
                        thermostatData.referenceTemperature[temperatureGroup]);
    }

    //! Connect with propagator - Nose-Hoover scales start and end step velocities
    void connectWithPropagator(const PropagatorConnection& connectionData, int numTemperatureGroups) override
    {
        GMX_RELEASE_ASSERT(
                connectionData.hasStartVelocityScaling() && connectionData.hasEndVelocityScaling(),
                "Nose-Hoover T-coupling requires both start and end velocity scaling.");
        connectionData.setNumVelocityScalingVariables(numTemperatureGroups,
                                                      ScaleVelocities::PreStepAndPostStep);
        lambdaStartVelocities_ = connectionData.getViewOnStartVelocityScaling();
        lambdaEndVelocities_   = connectionData.getViewOnEndVelocityScaling();
    }

    //! Constructor
    NoseHooverTemperatureCoupling(int                  numTemperatureGroups,
                                  ArrayRef<const real> referenceTemperature,
                                  ArrayRef<const real> couplingTime)
    {
        xi_.resize(numTemperatureGroups, 0.0);
        xiVelocities_.resize(numTemperatureGroups, 0.0);
        invXiMass_.resize(numTemperatureGroups, 0.0);
        for (auto temperatureGroup = 0; temperatureGroup < numTemperatureGroups; ++temperatureGroup)
        {
            if (referenceTemperature[temperatureGroup] > 0 && couplingTime[temperatureGroup] > 0)
            {
                // Note: This mass definition is equal to legacy md
                //       legacy md-vv divides the mass by ndof * kB
                invXiMass_[temperatureGroup] = 1.0
                                               / (gmx::square(couplingTime[temperatureGroup] / M_2PI)
                                                  * referenceTemperature[temperatureGroup]);
            }
        }
    }

    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData)
    {
        checkpointVersion(checkpointData, "Nose-Hoover version", c_nhCurrentVersion);
        checkpointData->arrayRef("xi", makeCheckpointArrayRef<operation>(xi_));
        checkpointData->arrayRef("xi velocities", makeCheckpointArrayRef<operation>(xiVelocities_));
    }

    //! Write thermostat dof to checkpoint
    void writeCheckpoint(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr) override
    {
        if (MAIN(cr))
        {
            doCheckpointData(&checkpointData.value());
        }
    }
    //! Read thermostat dof from checkpoint
    void readCheckpoint(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) override
    {
        if (MAIN(cr))
        {
            doCheckpointData(&checkpointData.value());
        }
        if (haveDDAtomOrdering(*cr))
        {
            dd_bcast(cr->dd, xi_.size() * sizeof(real), xi_.data());
            dd_bcast(cr->dd, xiVelocities_.size() * sizeof(real), xiVelocities_.data());
        }
    }

    //! Adapt masses
    real updateReferenceTemperatureAndIntegral(int             temperatureGroup,
                                               real gmx_unused newTemperature,
                                               ReferenceTemperatureChangeAlgorithm gmx_unused algorithm,
                                               const TemperatureCouplingData& temperatureCouplingData) override
    {
        // Currently, we don't know about any temperature change algorithms, so we assert this never gets called
        GMX_ASSERT(false,
                   "NoseHooverTemperatureCoupling: Unknown ReferenceTemperatureChangeAlgorithm.");
        const bool newTemperatureIsValid =
                (newTemperature > 0 && temperatureCouplingData.couplingTime[temperatureGroup] > 0
                 && temperatureCouplingData.numDegreesOfFreedom[temperatureGroup] > 0);
        const bool oldTemperatureIsValid =
                (temperatureCouplingData.referenceTemperature[temperatureGroup] > 0
                 && temperatureCouplingData.couplingTime[temperatureGroup] > 0
                 && temperatureCouplingData.numDegreesOfFreedom[temperatureGroup] > 0);
        GMX_RELEASE_ASSERT(newTemperatureIsValid == oldTemperatureIsValid,
                           "Cannot turn temperature coupling on / off during simulation run.");
        if (oldTemperatureIsValid && newTemperatureIsValid)
        {
            invXiMass_[temperatureGroup] *=
                    (temperatureCouplingData.referenceTemperature[temperatureGroup] / newTemperature);
            xiVelocities_[temperatureGroup] *= std::sqrt(
                    newTemperature / temperatureCouplingData.referenceTemperature[temperatureGroup]);
        }
        return integral(temperatureGroup,
                        temperatureCouplingData.numDegreesOfFreedom[temperatureGroup],
                        newTemperature);
    }

private:
    //! The thermostat degree of freedom
    std::vector<real> xi_;
    //! Velocity of the thermostat dof
    std::vector<real> xiVelocities_;
    //! Inverse mass of the thermostat dof
    std::vector<real> invXiMass_;

    //! View on the scaling factor of the propagator (pre-step velocities)
    ArrayRef<real> lambdaStartVelocities_;
    //! View on the scaling factor of the propagator (post-step velocities)
    ArrayRef<real> lambdaEndVelocities_;
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
        TemperatureCoupling               couplingType) :
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
    energyData_(energyData),
    nextEnergyCalculationStep_(-1)
{
    if (couplingType == TemperatureCoupling::VRescale)
    {
        temperatureCouplingImpl_ = std::make_unique<VRescaleTemperatureCoupling>(seed);
    }
    else if (couplingType == TemperatureCoupling::Berendsen)
    {
        temperatureCouplingImpl_ = std::make_unique<BerendsenTemperatureCoupling>();
    }
    else if (couplingType == TemperatureCoupling::NoseHoover)
    {
        temperatureCouplingImpl_ = std::make_unique<NoseHooverTemperatureCoupling>(
                numTemperatureGroups_, referenceTemperature_, couplingTime_);
    }
    else
    {
        throw NotImplementedError("Temperature coupling " + std::string(enumValueToString(couplingType))
                                  + " is not implemented for modular simulator.");
    }
    energyData->addConservedEnergyContribution([this](Step gmx_used_in_debug step, Time /*unused*/) {
        GMX_ASSERT(conservedEnergyContributionStep_ == step,
                   "VelocityScalingTemperatureCoupling conserved energy step mismatch.");
        return conservedEnergyContribution_;
    });
}

void VelocityScalingTemperatureCoupling::connectWithMatchingPropagator(const PropagatorConnection& connectionData,
                                                                       const PropagatorTag& propagatorTag)
{
    if (connectionData.tag == propagatorTag)
    {
        temperatureCouplingImpl_->connectWithPropagator(connectionData, numTemperatureGroups_);
        propagatorCallback_ = connectionData.getVelocityScalingCallback();
    }
}

void VelocityScalingTemperatureCoupling::elementSetup()
{
    if (!propagatorCallback_)
    {
        throw MissingElementConnectionError(
                "Velocity scaling temperature coupling was not connected to a propagator.\n"
                "Connection to a propagator element is needed to scale the velocities.\n"
                "Use connectWithMatchingPropagator(...) before building the "
                "ModularSimulatorAlgorithm "
                "object.");
    }
}

void VelocityScalingTemperatureCoupling::scheduleTask(Step                       step,
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
    if (step == nextEnergyCalculationStep_
        && reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::Yes)
    {
        // add conserved energy before we do T-coupling
        registerRunFunction([this, step]() {
            conservedEnergyContribution_     = conservedEnergyContribution();
            conservedEnergyContributionStep_ = step;
        });
    }
    if (do_per_step(step + nstcouple_ + offset_, nstcouple_))
    {
        // do T-coupling this step
        registerRunFunction([this, step]() { setLambda(step); });

        // Let propagator know that we want to do T-coupling
        propagatorCallback_(step);
    }
    if (step == nextEnergyCalculationStep_
        && reportPreviousConservedEnergy_ == ReportPreviousStepConservedEnergy::No)
    {
        // add conserved energy after we did T-coupling
        registerRunFunction([this, step]() {
            conservedEnergyContribution_     = conservedEnergyContribution();
            conservedEnergyContributionStep_ = step;
        });
    }
}

void VelocityScalingTemperatureCoupling::setLambda(Step step)
{
    const auto*             ekind          = energyData_->ekindata();
    TemperatureCouplingData thermostatData = {
        couplingTimeStep_, referenceTemperature_, couplingTime_, numDegreesOfFreedom_, temperatureCouplingIntegral_
    };

    for (int temperatureGroup = 0; (temperatureGroup < numTemperatureGroups_); temperatureGroup++)
    {
        const real currentKineticEnergy = useFullStepKE_ == UseFullStepKE::Yes
                                                  ? ::trace(ekind->tcstat[temperatureGroup].ekinf)
                                                  : ::trace(ekind->tcstat[temperatureGroup].ekinh);
        const real currentTemperature   = useFullStepKE_ == UseFullStepKE::Yes
                                                  ? ekind->tcstat[temperatureGroup].T
                                                  : ekind->tcstat[temperatureGroup].Th;

        temperatureCouplingIntegral_[temperatureGroup] = temperatureCouplingImpl_->apply(
                step, temperatureGroup, currentKineticEnergy, currentTemperature, thermostatData);
    }
}

void VelocityScalingTemperatureCoupling::updateReferenceTemperature(ArrayRef<const real> temperatures,
                                                                    ReferenceTemperatureChangeAlgorithm algorithm)
{
    TemperatureCouplingData thermostatData = {
        couplingTimeStep_, referenceTemperature_, couplingTime_, numDegreesOfFreedom_, temperatureCouplingIntegral_
    };
    for (int temperatureGroup = 0; (temperatureGroup < numTemperatureGroups_); temperatureGroup++)
    {
        temperatureCouplingIntegral_[temperatureGroup] =
                temperatureCouplingImpl_->updateReferenceTemperatureAndIntegral(
                        temperatureGroup, temperatures[temperatureGroup], algorithm, thermostatData);
    }
    // Currently, we don't know about any temperature change algorithms, so we assert this never gets called
    GMX_ASSERT(false,
               "VelocityScalingTemperatureCoupling: Unknown ReferenceTemperatureChangeAlgorithm.");
    std::copy(temperatures.begin(), temperatures.end(), referenceTemperature_.begin());
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
void VelocityScalingTemperatureCoupling::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "VRescaleThermostat version", c_currentVersion);

    checkpointData->arrayRef("thermostat integral",
                             makeCheckpointArrayRef<operation>(temperatureCouplingIntegral_));
}

void VelocityScalingTemperatureCoupling::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                                             const t_commrec*                   cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
    temperatureCouplingImpl_->writeCheckpoint(
            checkpointData
                    ? std::make_optional(checkpointData->subCheckpointData("thermostat impl"))
                    : std::nullopt,
            cr);
}

void VelocityScalingTemperatureCoupling::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                                const t_commrec* cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    if (haveDDAtomOrdering(*cr))
    {
        dd_bcast(cr->dd,
                 gmx::ssize(temperatureCouplingIntegral_) * int(sizeof(double)),
                 temperatureCouplingIntegral_.data());
    }
    temperatureCouplingImpl_->readCheckpoint(
            checkpointData
                    ? std::make_optional(checkpointData->subCheckpointData("thermostat impl"))
                    : std::nullopt,
            cr);
}

const std::string& VelocityScalingTemperatureCoupling::clientID()
{
    return identifier_;
}

real VelocityScalingTemperatureCoupling::conservedEnergyContribution() const
{
    return std::accumulate(temperatureCouplingIntegral_.begin(), temperatureCouplingIntegral_.end(), 0.0);
}

std::optional<SignallerCallback> VelocityScalingTemperatureCoupling::registerEnergyCallback(EnergySignallerEvent event)
{
    if (event == EnergySignallerEvent::EnergyCalculationStep)
    {
        return [this](Step step, Time /*unused*/) { nextEnergyCalculationStep_ = step; };
    }
    return std::nullopt;
}

ISimulatorElement* VelocityScalingTemperatureCoupling::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer* /*observablesReducer*/,
        Offset                            offset,
        UseFullStepKE                     useFullStepKE,
        ReportPreviousStepConservedEnergy reportPreviousStepConservedEnergy,
        const PropagatorTag&              propagatorTag)
{
    // Element is now owned by the caller of this method, who will handle lifetime (see ModularSimulatorAlgorithm)
    auto* element = builderHelper->storeElement(std::make_unique<VelocityScalingTemperatureCoupling>(
            legacySimulatorData->inputRec_->nsttcouple,
            offset,
            useFullStepKE,
            reportPreviousStepConservedEnergy,
            legacySimulatorData->inputRec_->ld_seed,
            legacySimulatorData->inputRec_->opts.ngtc,
            legacySimulatorData->inputRec_->delta_t * legacySimulatorData->inputRec_->nsttcouple,
            legacySimulatorData->inputRec_->opts.ref_t,
            legacySimulatorData->inputRec_->opts.tau_t,
            legacySimulatorData->inputRec_->opts.nrdf,
            energyData,
            legacySimulatorData->inputRec_->etc));
    auto* thermostat = static_cast<VelocityScalingTemperatureCoupling*>(element);
    // Capturing pointer is safe because lifetime is handled by caller
    builderHelper->registerTemperaturePressureControl(
            [thermostat, propagatorTag](const PropagatorConnection& connection) {
                thermostat->connectWithMatchingPropagator(connection, propagatorTag);
            });
    builderHelper->registerReferenceTemperatureUpdate(
            [thermostat](ArrayRef<const real> temperatures, ReferenceTemperatureChangeAlgorithm algorithm) {
                thermostat->updateReferenceTemperature(temperatures, algorithm);
            });
    return element;
}

} // namespace gmx
