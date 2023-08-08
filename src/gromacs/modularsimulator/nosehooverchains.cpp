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
 * \brief Defines classes related to Nose-Hoover chains for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#include "gmxpre.h"

#include "nosehooverchains.h"

#include <numeric>

#include "gromacs/domdec/domdec_network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/strconvert.h"

#include "energydata.h"
#include "mttk.h"
#include "simulatoralgorithm.h"
#include "trotterhelperfunctions.h"
#include "velocityscalingtemperaturecoupling.h"

namespace gmx
{
// Names of the NHC usage options
static constexpr EnumerationArray<NhcUsage, const char*> nhcUsageNames = { "System", "Barostat" };

//! The current state of the Nose-Hoover chain degree of freedom for a temperature group
class NoseHooverGroup final
{
public:
    //! Constructor
    NoseHooverGroup(int      chainLength,
                    real     referenceTemperature,
                    real     numDegreesOfFreedom,
                    real     couplingTime,
                    real     couplingTimeStep,
                    NhcUsage nhcUsage);

    //! Trotter operator for the NHC degrees of freedom
    real applyNhc(real currentKineticEnergy, real couplingTimeStep);

    //! Save to or restore from a CheckpointData object
    template<CheckpointDataOperation operation>
    void doCheckpoint(CheckpointData<operation>* checkpointData);
    //! Broadcast values read from checkpoint over DD ranks
    void broadcastCheckpointValues(const gmx_domdec_t* dd);

    //! Whether the coordinate time is at a full coupling time step
    bool isAtFullCouplingTimeStep() const;
    //! Update the value of the NHC integral with the current coordinates
    void calculateIntegral();
    //! Return the current NHC integral for the group
    double integral() const;
    //! Return the current time of the NHC integral for the group
    real integralTime() const;

    //! Set the reference temperature
    void updateReferenceTemperature(real temperature);

private:
    //! Increment coordinate time and update integral if applicable
    void finalizeUpdate(real couplingTimeStep);

    //! The reference temperature of this group
    real referenceTemperature_;
    //! The coupling time of this group
    const real couplingTime_;
    //! The number of degrees of freedom in this group
    const real numDegreesOfFreedom_;
    //! The chain length of this group
    const int chainLength_;
    //! The coupling time step, indicates when the coordinates are at a full step
    const real couplingTimeStep_;
    //! The thermostat degree of freedom
    std::vector<real> xi_;
    //! Velocity of the thermostat dof
    std::vector<real> xiVelocities_;
    //! Work exerted by thermostat per group
    double temperatureCouplingIntegral_;
    //! Inverse mass of the thermostat dof
    std::vector<real> invXiMass_;
    //! The current time of xi and xiVelocities_
    real coordinateTime_;
    //! The current time of the temperature integral
    real integralTime_;
};

NoseHooverGroup::NoseHooverGroup(int      chainLength,
                                 real     referenceTemperature,
                                 real     numDegreesOfFreedom,
                                 real     couplingTime,
                                 real     couplingTimeStep,
                                 NhcUsage nhcUsage) :
    referenceTemperature_(referenceTemperature),
    couplingTime_(couplingTime),
    numDegreesOfFreedom_(numDegreesOfFreedom),
    chainLength_(chainLength),
    couplingTimeStep_(couplingTimeStep),
    xi_(chainLength, 0),
    xiVelocities_(chainLength, 0),
    temperatureCouplingIntegral_(0),
    invXiMass_(chainLength, 0),
    coordinateTime_(0),
    integralTime_(0)
{
    if (referenceTemperature > 0 && couplingTime > 0 && numDegreesOfFreedom > 0)
    {
        for (auto chainPosition = 0; chainPosition < chainLength; ++chainPosition)
        {
            const real numDof = ((chainPosition == 0) ? numDegreesOfFreedom : 1);
            invXiMass_[chainPosition] =
                    1.0 / (gmx::square(couplingTime / M_2PI) * referenceTemperature * numDof * c_boltz);
            if (nhcUsage == NhcUsage::Barostat && chainPosition == 0)
            {
                invXiMass_[chainPosition] /= DIM * DIM;
            }
        }
    }
}

NoseHooverChainsData::NoseHooverChainsData(int                  numTemperatureGroups,
                                           real                 couplingTimeStep,
                                           int                  chainLength,
                                           ArrayRef<const real> referenceTemperature,
                                           ArrayRef<const real> couplingTime,
                                           ArrayRef<const real> numDegreesOfFreedom,
                                           NhcUsage             nhcUsage) :
    identifier_(formatString("NoseHooverChainsData-%s", nhcUsageNames[nhcUsage])),
    numTemperatureGroups_(numTemperatureGroups)
{
    if (nhcUsage == NhcUsage::System)
    {
        for (auto temperatureGroup = 0; temperatureGroup < numTemperatureGroups; ++temperatureGroup)
        {
            noseHooverGroups_.emplace_back(chainLength,
                                           referenceTemperature[temperatureGroup],
                                           numDegreesOfFreedom[temperatureGroup],
                                           couplingTime[temperatureGroup],
                                           couplingTimeStep,
                                           nhcUsage);
        }
    }
    else if (nhcUsage == NhcUsage::Barostat)
    {
        GMX_RELEASE_ASSERT(numTemperatureGroups == 1,
                           "There can only be one barostat for the system");
        // Barostat has a single degree of freedom
        const int degreesOfFreedom = 1;
        noseHooverGroups_.emplace_back(
                chainLength, referenceTemperature[0], degreesOfFreedom, couplingTime[0], couplingTimeStep, nhcUsage);
    }
}

NoseHooverChainsData::NoseHooverChainsData(const NoseHooverChainsData& other) :
    identifier_(other.identifier_),
    noseHooverGroups_(other.noseHooverGroups_),
    numTemperatureGroups_(other.numTemperatureGroups_)
{
}

void NoseHooverChainsData::build(NhcUsage                                nhcUsage,
                                 LegacySimulatorData*                    legacySimulatorData,
                                 ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                 EnergyData*                             energyData)
{
    // The caller of this method now owns the data and will handle its lifetime (see ModularSimulatorAlgorithm)
    if (nhcUsage == NhcUsage::System)
    {
        builderHelper->storeSimulationData(
                NoseHooverChainsData::dataID(nhcUsage),
                NoseHooverChainsData(
                        legacySimulatorData->inputRec_->opts.ngtc,
                        legacySimulatorData->inputRec_->delta_t * legacySimulatorData->inputRec_->nsttcouple,
                        legacySimulatorData->inputRec_->opts.nhchainlength,
                        constArrayRefFromArray(legacySimulatorData->inputRec_->opts.ref_t,
                                               legacySimulatorData->inputRec_->opts.ngtc),
                        constArrayRefFromArray(legacySimulatorData->inputRec_->opts.tau_t,
                                               legacySimulatorData->inputRec_->opts.ngtc),
                        constArrayRefFromArray(legacySimulatorData->inputRec_->opts.nrdf,
                                               legacySimulatorData->inputRec_->opts.ngtc),
                        nhcUsage));
    }
    else
    {
        const int numTemperatureGroups = 1;
        builderHelper->storeSimulationData(
                NoseHooverChainsData::dataID(nhcUsage),
                NoseHooverChainsData(
                        numTemperatureGroups,
                        legacySimulatorData->inputRec_->delta_t
                                * legacySimulatorData->inputRec_->pressureCouplingOptions.nstpcouple,
                        legacySimulatorData->inputRec_->opts.nhchainlength,
                        constArrayRefFromArray(legacySimulatorData->inputRec_->opts.ref_t, 1),
                        constArrayRefFromArray(legacySimulatorData->inputRec_->opts.tau_t, 1),
                        ArrayRef<real>(),
                        nhcUsage));
    }
    auto* nhcDataPtr =
            builderHelper
                    ->simulationData<NoseHooverChainsData>(NoseHooverChainsData::dataID(nhcUsage))
                    .value();
    builderHelper->registerReferenceTemperatureUpdate(
            [nhcDataPtr](ArrayRef<const real> temperatures, ReferenceTemperatureChangeAlgorithm algorithm) {
                nhcDataPtr->updateReferenceTemperature(temperatures, algorithm);
            });

    const auto* ptrToDataObject =
            builderHelper
                    ->simulationData<NoseHooverChainsData>(NoseHooverChainsData::dataID(nhcUsage))
                    .value();
    energyData->addConservedEnergyContribution([ptrToDataObject](Step /*unused*/, Time time) {
        return ptrToDataObject->temperatureCouplingIntegral(time);
    });
}

NoseHooverChainsData::~NoseHooverChainsData() = default;

void NoseHooverGroup::finalizeUpdate(real couplingTimeStep)
{
    coordinateTime_ += couplingTimeStep;
    if (isAtFullCouplingTimeStep())
    {
        calculateIntegral();
    }
}

inline int NoseHooverChainsData::numTemperatureGroups() const
{
    return numTemperatureGroups_;
}

inline bool NoseHooverChainsData::isAtFullCouplingTimeStep() const
{
    return std::all_of(noseHooverGroups_.begin(), noseHooverGroups_.end(), [](const auto& group) {
        return group.isAtFullCouplingTimeStep();
    });
}

void NoseHooverGroup::calculateIntegral()
{
    // Calculate current value of thermostat integral
    temperatureCouplingIntegral_ = 0.0;
    for (auto chainPosition = 0; chainPosition < chainLength_; ++chainPosition)
    {
        // Chain thermostats have only one degree of freedom
        const real numDegreesOfFreedomThisPosition = (chainPosition == 0) ? numDegreesOfFreedom_ : 1;
        temperatureCouplingIntegral_ +=
                0.5 * gmx::square(xiVelocities_[chainPosition]) / invXiMass_[chainPosition]
                + numDegreesOfFreedomThisPosition * xi_[chainPosition] * c_boltz * referenceTemperature_;
    }
    integralTime_ = coordinateTime_;
}

inline double NoseHooverGroup::integral() const
{
    return temperatureCouplingIntegral_;
}

inline real NoseHooverGroup::integralTime() const
{
    return integralTime_;
}

double NoseHooverChainsData::temperatureCouplingIntegral(Time gmx_used_in_debug time) const
{
    /* When using nsttcouple >= nstcalcenergy, we accept that the coupling
     * integral might be ahead of the current energy calculation step. The
     * extended system degrees of freedom are either in sync or ahead of the
     * rest of the system.
     */
    GMX_ASSERT(!std::any_of(noseHooverGroups_.begin(),
                            noseHooverGroups_.end(),
                            [time](const auto& group) {
                                return !(time <= group.integralTime()
                                         || timesClose(group.integralTime(), time));
                            }),
               "NoseHooverChainsData conserved energy time mismatch.");
    double result = 0;
    std::for_each(noseHooverGroups_.begin(), noseHooverGroups_.end(), [&result](const auto& group) {
        result += group.integral();
    });
    return result;
}

inline bool NoseHooverGroup::isAtFullCouplingTimeStep() const
{
    // Check whether coordinate time divided by the time step is close to integer
    return timesClose(std::lround(coordinateTime_ / couplingTimeStep_) * couplingTimeStep_, coordinateTime_);
}

void NoseHooverChainsData::updateReferenceTemperature(ArrayRef<const real> temperatures,
                                                      ReferenceTemperatureChangeAlgorithm gmx_unused algorithm)
{
    // Currently, we don't know about any temperature change algorithms, so we assert this never gets called
    GMX_ASSERT(false, "NoseHooverChainsData: Unknown ReferenceTemperatureChangeAlgorithm.");
    for (auto temperatureGroup = 0; temperatureGroup < numTemperatureGroups_; ++temperatureGroup)
    {
        noseHooverGroups_[temperatureGroup].updateReferenceTemperature(temperatures[temperatureGroup]);
        if (noseHooverGroups_[temperatureGroup].isAtFullCouplingTimeStep())
        {
            noseHooverGroups_[temperatureGroup].calculateIntegral();
        }
    }
}

void NoseHooverGroup::updateReferenceTemperature(real temperature)
{
    const bool newTemperatureIsValid = (temperature > 0 && couplingTime_ > 0 && numDegreesOfFreedom_ > 0);
    const bool oldTemperatureIsValid =
            (referenceTemperature_ > 0 && couplingTime_ > 0 && numDegreesOfFreedom_ > 0);
    GMX_RELEASE_ASSERT(newTemperatureIsValid == oldTemperatureIsValid,
                       "Cannot turn temperature coupling on / off during simulation run.");
    if (oldTemperatureIsValid && newTemperatureIsValid)
    {
        const real velocityFactor = std::sqrt(temperature / referenceTemperature_);
        for (auto chainPosition = 0; chainPosition < chainLength_; ++chainPosition)
        {
            invXiMass_[chainPosition] *= (referenceTemperature_ / temperature);
            xiVelocities_[chainPosition] *= velocityFactor;
        }
    }
    referenceTemperature_ = temperature;
    if (isAtFullCouplingTimeStep())
    {
        calculateIntegral();
    }
}

namespace
{
/*!
 * \brief Enum describing the contents NoseHooverChainsData writes to modular checkpoint
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
void NoseHooverChainsData::doCheckpointData(CheckpointData<operation>* checkpointData)
{
    checkpointVersion(checkpointData, "NoseHooverChainsData version", c_currentVersion);

    for (int temperatureGroup = 0; temperatureGroup < numTemperatureGroups_; ++temperatureGroup)
    {
        const auto temperatureGroupStr = "T-group #" + toString(temperatureGroup);
        auto       groupCheckpointData = checkpointData->subCheckpointData(temperatureGroupStr);
        noseHooverGroups_[temperatureGroup].doCheckpoint(&groupCheckpointData);
    }
}

template<CheckpointDataOperation operation>
void NoseHooverGroup::doCheckpoint(CheckpointData<operation>* checkpointData)
{
    checkpointData->arrayRef("xi", makeCheckpointArrayRef<operation>(xi_));
    checkpointData->arrayRef("xi velocities", makeCheckpointArrayRef<operation>(xiVelocities_));
    checkpointData->scalar("Coordinate time", &coordinateTime_);
}

//! Broadcast values read from checkpoint over DD ranks
void NoseHooverGroup::broadcastCheckpointValues(const gmx_domdec_t* dd)
{
    dd_bcast(dd, gmx::ssize(xi_) * int(sizeof(real)), xi_.data());
    dd_bcast(dd, gmx::ssize(xiVelocities_) * int(sizeof(real)), xiVelocities_.data());
    dd_bcast(dd, int(sizeof(real)), &coordinateTime_);
}

void NoseHooverChainsData::saveCheckpointState(std::optional<WriteCheckpointData> checkpointData,
                                               const t_commrec*                   cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Write>(&checkpointData.value());
    }
}

void NoseHooverChainsData::restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData,
                                                  const t_commrec*                  cr)
{
    if (MAIN(cr))
    {
        doCheckpointData<CheckpointDataOperation::Read>(&checkpointData.value());
    }
    for (auto& group : noseHooverGroups_)
    {
        if (haveDDAtomOrdering(*cr))
        {
            group.broadcastCheckpointValues(cr->dd);
        }
        group.calculateIntegral();
    }
}

const std::string& NoseHooverChainsData::clientID()
{
    return identifier_;
}

std::string NoseHooverChainsData::dataID(NhcUsage nhcUsage)
{
    return formatString("NoseHooverChainsData%s", nhcUsageNames[nhcUsage]);
}

/* This follows Tuckerman et al. 2006
 *
 * In NVT, the Trotter decomposition reads
 *   exp[iL dt] = exp[iLT dt/2] exp[iLv dt/2] exp[iLx dt] exp[iLv dt/2] exp[iLT dt/2]
 * iLv denotes the velocity propagation, iLx the position propagation
 * iLT denotes the thermostat propagation implemented here:
 *     v_xi[i](t-dt/2) = v_xi[i](t-dt) + dt_xi/2 * a_xi[i](t-dt);
 *     xi[i](t) = xi[i](t-dt) + dt_xi * v_xi[i](t-dt/2);
 *     v_sys *= exp(-dt/2 * v_xi[1](t-dt/2))
 *     v_xi[i](t) = v_xi[i](t-dt/2) + dt_xi/2 * a_xi[i](t);
 * where i = 1 ... N_chain, and
 *     a[i](t) = (M_xi * v_xi[i+1](t)^2 - 2*K_ref) / M_xi , i = 2 ... N_chain
 *     a[1](t) = (K_sys - K_ref) / M_xi
 * Note, iLT contains a term scaling the system velocities!
 *
 * In the legacy GROMACS simulator, the top of the loop marks the simulation
 * state at x(t), v(t-dt/2), f(t-1), mirroring the leap-frog implementation.
 * The loop then proceeds to calculate the forces at time t, followed by a
 * velocity half step (corresponding to the second exp[iLv dt/2] above).
 * For Tuckerman NHC NVT, this is followed by a thermostat propagation to reach
 * the full timestep t state. This is the state which is printed to file, so
 * we need to scale the velocities.
 * After writing to file, the next step effectively starts, by moving the thermostat
 * variables (half step), the velocities (half step) and the positions (full step),
 * which is equivalent to the first three terms of the Trotter decomposition above.
 * Currently, modular simulator is replicating the division of the simulator loop
 * used by the legacy simulator. The implementation here is independent of these
 * assumptions, but the builder of the simulator must be careful to ensure that
 * velocity scaling is applied before re-using the velocities after the thermostat.
 *
 * The time-scale separation between the particles and the thermostat requires the
 * NHC operator to have a higher-order factorization. The method used is the
 * Suzuki-Yoshida scheme which uses weighted time steps chosen to cancel out
 * lower-order error terms. Here, the fifth order SY scheme is used.
 */
real NoseHooverGroup::applyNhc(real currentKineticEnergy, const real couplingTimeStep)
{
    if (currentKineticEnergy < 0)
    {
        finalizeUpdate(couplingTimeStep);
        return 1.0;
    }

    constexpr unsigned int c_suzukiYoshidaOrder                         = 5;
    constexpr double       c_suzukiYoshidaWeights[c_suzukiYoshidaOrder] = {
        0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065, 0.2967324292201065
    };

    real velocityScalingFactor = 1.0;

    // Apply Suzuki-Yoshida scheme
    for (unsigned int syOuterLoop = 0; syOuterLoop < c_suzukiYoshidaOrder; ++syOuterLoop)
    {
        for (unsigned int syInnerLoop = 0; syInnerLoop < c_suzukiYoshidaOrder; ++syInnerLoop)
        {
            const real timeStep =
                    couplingTimeStep * c_suzukiYoshidaWeights[syInnerLoop] / c_suzukiYoshidaOrder;

            // Reverse loop - start from last thermostat in chain to update velocities,
            // because we need the new velocity to scale the next thermostat in the chain
            for (auto chainPosition = chainLength_ - 1; chainPosition >= 0; --chainPosition)
            {
                const real kineticEnergy2 =
                        ((chainPosition == 0) ? 2 * currentKineticEnergy
                                              : gmx::square(xiVelocities_[chainPosition - 1])
                                                        / invXiMass_[chainPosition - 1]);
                // DOF of temperature group or chain member
                const real numDof         = ((chainPosition == 0) ? numDegreesOfFreedom_ : 1);
                const real xiAcceleration = invXiMass_[chainPosition]
                                            * (kineticEnergy2 - numDof * c_boltz * referenceTemperature_);

                // We scale based on the next thermostat in chain.
                // Last thermostat in chain doesn't get scaled.
                const real localScalingFactor =
                        (chainPosition < chainLength_ - 1)
                                ? std::exp(-0.25 * timeStep * xiVelocities_[chainPosition + 1])
                                : 1.0;
                xiVelocities_[chainPosition] = localScalingFactor
                                               * (xiVelocities_[chainPosition] * localScalingFactor
                                                  + 0.5 * timeStep * xiAcceleration);
            }

            // Calculate the new system scaling factor
            const real systemScalingFactor = std::exp(-timeStep * xiVelocities_[0]);
            velocityScalingFactor *= systemScalingFactor;
            currentKineticEnergy *= systemScalingFactor * systemScalingFactor;

            // Forward loop - start from the system thermostat
            for (auto chainPosition = 0; chainPosition < chainLength_; ++chainPosition)
            {
                // Update thermostat positions
                xi_[chainPosition] += timeStep * xiVelocities_[chainPosition];

                // Kinetic energy of system or previous chain member
                const real kineticEnergy2 =
                        ((chainPosition == 0) ? 2 * currentKineticEnergy
                                              : gmx::square(xiVelocities_[chainPosition - 1])
                                                        / invXiMass_[chainPosition - 1]);
                // DOF of temperature group or chain member
                const real numDof         = ((chainPosition == 0) ? numDegreesOfFreedom_ : 1);
                const real xiAcceleration = invXiMass_[chainPosition]
                                            * (kineticEnergy2 - numDof * c_boltz * referenceTemperature_);

                // We scale based on the next thermostat in chain.
                // Last thermostat in chain doesn't get scaled.
                const real localScalingFactor =
                        (chainPosition < chainLength_ - 1)
                                ? std::exp(-0.25 * timeStep * xiVelocities_[chainPosition + 1])
                                : 1.0;
                xiVelocities_[chainPosition] = localScalingFactor
                                               * (xiVelocities_[chainPosition] * localScalingFactor
                                                  + 0.5 * timeStep * xiAcceleration);
            }
        }
    }
    finalizeUpdate(couplingTimeStep);
    return velocityScalingFactor;
}

real NoseHooverChainsData::applyNhc(int temperatureGroup, double propagationTimeStep, real currentKineticEnergy)
{
    return noseHooverGroups_[temperatureGroup].applyNhc(currentKineticEnergy, propagationTimeStep);
}

/*!
 * \brief Calculate the current kinetic energy
 *
 * \param tcstat  The group's kinetic energy structure
 * \return real   The current kinetic energy
 */
inline real NoseHooverChainsElement::currentKineticEnergy(const t_grp_tcstat& tcstat)
{
    if (nhcUsage_ == NhcUsage::System)
    {
        if (useFullStepKE_ == UseFullStepKE::Yes)
        {
            return ::trace(tcstat.ekinf) * tcstat.ekinscalef_nhc;
        }
        else
        {
            return ::trace(tcstat.ekinh) * tcstat.ekinscaleh_nhc;
        }
    }
    else if (nhcUsage_ == NhcUsage::Barostat)
    {
        GMX_RELEASE_ASSERT(useFullStepKE_ == UseFullStepKE::Yes,
                           "Barostat NHC only works with full step KE.");
        return mttkData_->kineticEnergy();
    }
    else
    {
        gmx_fatal(FARGS, "Unknown NhcUsage.");
    }
}

void NoseHooverChainsElement::propagateNhc()
{
    auto* ekind = energyData_->ekindata();

    for (int temperatureGroup = 0; (temperatureGroup < noseHooverChainData_->numTemperatureGroups());
         temperatureGroup++)
    {
        auto scalingFactor =
                noseHooverChainData_->applyNhc(temperatureGroup,
                                               propagationTimeStep_,
                                               currentKineticEnergy(ekind->tcstat[temperatureGroup]));

        if (nhcUsage_ == NhcUsage::System)
        {
            // Scale system velocities by scalingFactor
            lambdaStartVelocities_[temperatureGroup] = scalingFactor;
            // Scale kinetic energy by scalingFactor^2
            ekind->tcstat[temperatureGroup].ekinscaleh_nhc *= scalingFactor * scalingFactor;
            ekind->tcstat[temperatureGroup].ekinscalef_nhc *= scalingFactor * scalingFactor;
        }
        else if (nhcUsage_ == NhcUsage::Barostat)
        {
            // Scale eta velocities by scalingFactor
            mttkData_->scale(scalingFactor, noseHooverChainData_->isAtFullCouplingTimeStep());
        }
    }

    if (nhcUsage_ == NhcUsage::System && noseHooverChainData_->isAtFullCouplingTimeStep())
    {
        // We've set the scaling factors for the full time step, so scale
        // kinetic energy accordingly before it gets printed
        energyData_->updateKineticEnergy();
    }
}

NoseHooverChainsElement::NoseHooverChainsElement(int                   nstcouple,
                                                 int                   offset,
                                                 NhcUsage              nhcUsage,
                                                 UseFullStepKE         useFullStepKE,
                                                 double                propagationTimeStep,
                                                 ScheduleOnInitStep    scheduleOnInitStep,
                                                 Step                  initStep,
                                                 EnergyData*           energyData,
                                                 NoseHooverChainsData* noseHooverChainData,
                                                 MttkData*             mttkData) :
    nsttcouple_(nstcouple),
    offset_(offset),
    propagationTimeStep_(propagationTimeStep),
    nhcUsage_(nhcUsage),
    useFullStepKE_(useFullStepKE),
    scheduleOnInitStep_(scheduleOnInitStep),
    initialStep_(initStep),
    energyData_(energyData),
    noseHooverChainData_(noseHooverChainData),
    mttkData_(mttkData)
{
}

void NoseHooverChainsElement::elementSetup()
{
    GMX_RELEASE_ASSERT(
            !(nhcUsage_ == NhcUsage::System && !propagatorCallback_),
            "Nose-Hoover chain element was not connected to a propagator.\n"
            "Connection to a propagator element is needed to scale the velocities.\n"
            "Use connectWithPropagator(...) before building the ModularSimulatorAlgorithm "
            "object.");
}

void NoseHooverChainsElement::scheduleTask(Step step, Time /*unused*/, const RegisterRunFunction& registerRunFunction)
{
    if (step == initialStep_ && scheduleOnInitStep_ == ScheduleOnInitStep::No)
    {
        return;
    }
    if (do_per_step(step + nsttcouple_ + offset_, nsttcouple_))
    {
        // do T-coupling this step
        registerRunFunction([this]() { propagateNhc(); });

        if (propagatorCallback_)
        {
            // Let propagator know that we want to do T-coupling
            propagatorCallback_(step);
        }
    }
}

void NoseHooverChainsElement::connectWithPropagator(const PropagatorConnection& connectionData,
                                                    const PropagatorTag&        propagatorTag)
{
    if (connectionData.tag == propagatorTag)
    {
        GMX_RELEASE_ASSERT(connectionData.hasStartVelocityScaling(),
                           "Trotter NHC needs start velocity scaling.");
        connectionData.setNumVelocityScalingVariables(noseHooverChainData_->numTemperatureGroups(),
                                                      ScaleVelocities::PreStepOnly);
        lambdaStartVelocities_ = connectionData.getViewOnStartVelocityScaling();
        propagatorCallback_    = connectionData.getVelocityScalingCallback();
    }
}

//! \cond
// Doxygen gets confused by the overload
ISimulatorElement* NoseHooverChainsElement::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused*  globalCommunicationHelper,
        ObservablesReducer*                    observablesReducer,
        NhcUsage                               nhcUsage,
        Offset                                 offset,
        UseFullStepKE                          useFullStepKE,
        ScheduleOnInitStep                     scheduleOnInitStep,
        const MttkPropagatorConnectionDetails& mttkPropagatorConnectionDetails)
{
    GMX_RELEASE_ASSERT(nhcUsage == NhcUsage::Barostat, "System NHC element needs a propagator tag.");
    if (!builderHelper->simulationData<MttkData>(MttkData::dataID()))
    {
        MttkData::build(legacySimulatorData, builderHelper, statePropagatorData, energyData, mttkPropagatorConnectionDetails);
    }
    return getElementPointerImpl(legacySimulatorData,
                                 builderHelper,
                                 statePropagatorData,
                                 energyData,
                                 freeEnergyPerturbationData,
                                 globalCommunicationHelper,
                                 observablesReducer,
                                 nhcUsage,
                                 offset,
                                 useFullStepKE,
                                 scheduleOnInitStep,
                                 PropagatorTag(""));
}

ISimulatorElement* NoseHooverChainsElement::getElementPointerImpl(
        LegacySimulatorData*                    legacySimulatorData,
        ModularSimulatorAlgorithmBuilderHelper* builderHelper,
        StatePropagatorData gmx_unused* statePropagatorData,
        EnergyData*                     energyData,
        FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
        GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
        ObservablesReducer gmx_unused* observablesReducer,
        NhcUsage                       nhcUsage,
        Offset                         offset,
        UseFullStepKE                  useFullStepKE,
        ScheduleOnInitStep             scheduleOnInitStep,
        const PropagatorTag&           propagatorTag)
{
    if (!builderHelper->simulationData<NoseHooverChainsData>(NoseHooverChainsData::dataID(nhcUsage)))
    {
        NoseHooverChainsData::build(nhcUsage, legacySimulatorData, builderHelper, energyData);
    }
    auto* nhcData = builderHelper
                            ->simulationData<NoseHooverChainsData>(NoseHooverChainsData::dataID(nhcUsage))
                            .value();

    // MTTK data is only needed when connecting to a barostat
    MttkData* mttkData = nullptr;
    if (nhcUsage == NhcUsage::Barostat)
    {
        mttkData = builderHelper->simulationData<MttkData>(MttkData::dataID()).value();
    }

    // Element is now owned by the caller of this method, who will handle lifetime (see ModularSimulatorAlgorithm)
    auto* element = builderHelper->storeElement(std::make_unique<NoseHooverChainsElement>(
            legacySimulatorData->inputRec_->nsttcouple,
            offset,
            nhcUsage,
            useFullStepKE,
            legacySimulatorData->inputRec_->delta_t * legacySimulatorData->inputRec_->nsttcouple / 2,
            scheduleOnInitStep,
            legacySimulatorData->inputRec_->init_step,
            energyData,
            nhcData,
            mttkData));
    if (nhcUsage == NhcUsage::System)
    {
        auto* thermostat = static_cast<NoseHooverChainsElement*>(element);
        // Capturing pointer is safe because caller handles lifetime
        builderHelper->registerTemperaturePressureControl(
                [thermostat, propagatorTag](const PropagatorConnection& connection) {
                    thermostat->connectWithPropagator(connection, propagatorTag);
                });
    }
    else
    {
        GMX_RELEASE_ASSERT(propagatorTag == PropagatorTag(""),
                           "Propagator tag is unused for Barostat NHC element.");
    }
    return element;
}
//! \endcond

} // namespace gmx
