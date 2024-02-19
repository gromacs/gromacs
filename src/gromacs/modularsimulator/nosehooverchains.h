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
 * \brief Declares classes related to Nose-Hoover chains for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_NOSEHOOVERCHAINS_H
#define GMX_MODULARSIMULATOR_NOSEHOOVERCHAINS_H

#include <queue>

#include "modularsimulatorinterfaces.h"

struct t_grp_tcstat;

namespace gmx
{
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class ModularSimulatorAlgorithmBuilderHelper;
struct MttkPropagatorConnectionDetails;
class MttkData;
class NoseHooverGroup;
class ObservablesReducer;
class StatePropagatorData;
enum class UseFullStepKE;

//! Whether the element does schedule on the initial step
enum class ScheduleOnInitStep
{
    Yes,  //!< Schedule on first step
    No,   //!< Do not schedule on first step
    Count //!< Number of enum entries
};

//! The usages of Nose-Hoover chains
enum class NhcUsage
{
    System,   //!< Couple system to temperature bath
    Barostat, //!< Couple barostat to temperature bath
    Count     //!< Number of enum entries
};

/*! \internal
 * \brief Class holding data used by the Nose-Hoover chains
 *
 * As the Trotter update is split in several sub-steps (i.e. is updated
 * by several element instances), the NHC degrees of freedom must be
 * stored centrally rather than by the single elements.
 *
 * This class manages these extra degrees of freedom. It controls access
 * (making sure that only one element has write access at a time), keeps
 * track of the current time stamp of the dofs, calculates the energy
 * related to the dof at the requested times, and writes the data needed
 * for restarts to checkpoint. As this is not implementing the
 * ISimulatorElement interface, it is not part of the simulator loop, but
 * relies on callbacks to perform its duties.
 */
class NoseHooverChainsData final : public ICheckpointHelperClient
{
public:
    //! Constructor
    NoseHooverChainsData(int                  numTemperatureGroups,
                         real                 couplingTimeStep,
                         int                  chainLength,
                         ArrayRef<const real> referenceTemperature,
                         ArrayRef<const real> couplingTime,
                         ArrayRef<const real> numDegreesOfFreedom,
                         NhcUsage             nhcUsage);

    //! Explicit default destructor
    ~NoseHooverChainsData();

    //! Explicit copy constructor (interface has a standard destructor)
    NoseHooverChainsData(const NoseHooverChainsData& other);

    /*! \brief Propagate the NHC degrees of freedom for a temperature group and
     *        return its current velocity scaling value
     *
     * \param temperatureGroup  The temperature group to be propagated
     * \param propagationTimeStep  The time step by which the DOF are propagated
     * \param currentKineticEnergy  The current kinetic energy of the temperature group
     * \return  The current velocity scaling value for the temperature group
     */
    real applyNhc(int temperatureGroup, double propagationTimeStep, real currentKineticEnergy);

    //! The number of temperature groups
    int numTemperatureGroups() const;
    //! Whether the NHC degrees of freedom are at a full coupling time step
    bool isAtFullCouplingTimeStep() const;

    //! ICheckpointHelperClient write checkpoint implementation
    void saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient read checkpoint implementation
    void restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient key implementation
    const std::string& clientID() override;

    //! Build object and store in builder helper object
    static void build(NhcUsage                                nhcUsage,
                      LegacySimulatorData*                    legacySimulatorData,
                      ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                      EnergyData*                             energyData);

    //! Identifier used to store objects
    static std::string dataID(NhcUsage nhcUsage);

private:
    //! Return the value of the coupling integral at a specific time
    double temperatureCouplingIntegral(Time time) const;
    //! Update the reference temperature
    void updateReferenceTemperature(ArrayRef<const real>                temperatures,
                                    ReferenceTemperatureChangeAlgorithm algorithm);

    //! CheckpointHelper identifier
    const std::string identifier_;
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);

    //! List of temperature groups with Nose-Hoover chains
    std::vector<NoseHooverGroup> noseHooverGroups_;

    //! The number of temperature groups
    const int numTemperatureGroups_;
};

/*! \internal
 * \brief Element propagating the Nose-Hoover chains
 *
 * This propagates the Nose-Hoover chain degrees of freedom, and
 * transmits the scaling factor to a connected propagator.
 */
class NoseHooverChainsElement final : public ISimulatorElement
{
public:
    //! Constructor
    NoseHooverChainsElement(int                   nstcouple,
                            int                   offset,
                            NhcUsage              nhcUsage,
                            UseFullStepKE         useFullStepKE,
                            double                propagationTimeStep,
                            ScheduleOnInitStep    scheduleOnInitStep,
                            Step                  initStep,
                            EnergyData*           energyData,
                            NoseHooverChainsData* noseHooverChainData,
                            MttkData*             mttkData);

    /*! \brief Register run function for step / time
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! Sanity check at setup time
    void elementSetup() override;
    //! No element teardown needed
    void elementTeardown() override {}

    //! Connect this to propagator
    void connectWithPropagator(const PropagatorConnection& connectionData,
                               const PropagatorTag&        propagatorTag);

    /*! \brief Factory method implementation (no propagator connection)
     *
     * This signature is used to connect a Nose-Hoover chain to a barostat
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param nhcUsage  What the NHC is connected to - system or barostat
     * \param offset  The step offset at which the thermostat is applied
     * \param useFullStepKE  Whether full step or half step KE is used
     * \param scheduleOnInitStep  Whether the element is scheduled on the initial step
     * \param mttkPropagatorConnectionDetails  Connection information for the MTTK barostat
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement*
    getElementPointerImpl(LegacySimulatorData*                    legacySimulatorData,
                          ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                          StatePropagatorData*                    statePropagatorData,
                          EnergyData*                             energyData,
                          FreeEnergyPerturbationData*             freeEnergyPerturbationData,
                          GlobalCommunicationHelper*              globalCommunicationHelper,
                          ObservablesReducer*                     observablesReducer,
                          NhcUsage                                nhcUsage,
                          Offset                                  offset,
                          UseFullStepKE                           useFullStepKE,
                          ScheduleOnInitStep                      scheduleOnInitStep,
                          const MttkPropagatorConnectionDetails&  mttkPropagatorConnectionDetails);

    /*! \brief Factory method implementation (including propagator connection)
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param nhcUsage  What the NHC is connected to - system or barostat
     * \param offset  The step offset at which the thermostat is applied
     * \param useFullStepKE  Whether full step or half step KE is used
     * \param scheduleOnInitStep  Whether the element is scheduled on the initial step
     * \param propagatorTag  Tag of the propagator to connect to
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer,
                                                    NhcUsage                   nhcUsage,
                                                    Offset                     offset,
                                                    UseFullStepKE              useFullStepKE,
                                                    ScheduleOnInitStep         scheduleOnInitStep,
                                                    const PropagatorTag&       propagatorTag);

private:
    //! Propagate the NHC degrees of freedom
    void propagateNhc();
    //! Helper function returning the appropriate kinetic energy
    real currentKineticEnergy(const t_grp_tcstat& tcstat);

    //! View on the scaling factor of the propagator (pre-step velocities)
    ArrayRef<real> lambdaStartVelocities_;
    //! Callback to let propagator know that we will update lambda
    PropagatorCallback propagatorCallback_;

    //! The period at which the thermostat is applied
    const int nsttcouple_;
    //! If != 0, offset the step at which the thermostat is applied
    const int offset_;
    //! The propagation time step - by how much we propagate the NHC dof
    const double propagationTimeStep_;
    //! Whether this NHC is acting on the system or a barostat
    const NhcUsage nhcUsage_;
    //! Whether we're using full step kinetic energy
    const UseFullStepKE useFullStepKE_;
    //! Whether we're scheduling on the first step
    const ScheduleOnInitStep scheduleOnInitStep_;
    //! The initial step number
    const Step initialStep_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the energy data (for ekindata)
    EnergyData* energyData_;
    //! Pointer to the NHC data
    NoseHooverChainsData* noseHooverChainData_;
    //! Pointer to the MTTK data (nullptr if this is not connected to barostat)
    MttkData* mttkData_;
};


} // namespace gmx

#endif // GMX_MODULARSIMULATOR_NOSEHOOVERCHAINS_H
