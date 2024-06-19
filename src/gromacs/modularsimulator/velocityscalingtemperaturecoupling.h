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
 * \brief Declares a velocity-scaling temperature coupling element for
 * the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_VELOCITYSCALINGTEMPERATURECOUPLING_H
#define GMX_MODULARSIMULATOR_VELOCITYSCALINGTEMPERATURECOUPLING_H

#include <cstdint>

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "energydata.h"
#include "modularsimulatorinterfaces.h"
#include "propagator.h"

struct t_commrec;
enum class TemperatureCoupling : int;

namespace gmx
{
class ITemperatureCouplingImpl;
class LegacySimulatorData;
class ObservablesReducer;
struct TemperatureCouplingData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class ModularSimulatorAlgorithmBuilderHelper;
class StatePropagatorData;
enum class ReferenceTemperatureChangeAlgorithm;
template<CheckpointDataOperation operation>
class CheckpointData;

//! Enum describing whether the thermostat is using full or half step kinetic energy
enum class UseFullStepKE
{
    Yes,
    No,
    Count
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element implementing the a velocity-scaling thermostat
 *
 * This element takes a callback to the propagator and updates the velocity
 * scaling factor according to the internal temperature coupling implementation.
 *
 * Note that the concrete implementation is handled by the concrete
 * implementations of the ITemperatureCouplingImpl interface, while the element
 * handles the scheduling and interfacing with other elements.
 */
class VelocityScalingTemperatureCoupling final :
    public ISimulatorElement,
    public ICheckpointHelperClient,
    public IEnergySignallerClient
{
public:
    //! Constructor
    VelocityScalingTemperatureCoupling(int                               nstcouple,
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
                                       TemperatureCoupling               couplingType);

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
    void connectWithMatchingPropagator(const PropagatorConnection& connectionData,
                                       const PropagatorTag&        propagatorTag);

    //! ICheckpointHelperClient write checkpoint implementation
    void saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient read checkpoint implementation
    void restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient key implementation
    const std::string& clientID() override;

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param propagatorTag  Tag of the propagator to connect to
     * \param offset  The step offset at which the thermostat is applied
     * \param useFullStepKE  Whether full step or half step KE is used
     * \param reportPreviousStepConservedEnergy  Report the previous or the current step conserved energy
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
                          Offset                                  offset,
                          UseFullStepKE                           useFullStepKE,
                          ReportPreviousStepConservedEnergy       reportPreviousStepConservedEnergy,
                          const PropagatorTag&                    propagatorTag);

private:
    //! Update the reference temperature
    void updateReferenceTemperature(ArrayRef<const real>                temperatures,
                                    ReferenceTemperatureChangeAlgorithm algorithm);

    //! The frequency at which the thermostat is applied
    const int nstcouple_;
    //! If != 0, offset the step at which the thermostat is applied
    const int offset_;
    //! Whether we're using full step kinetic energy
    const UseFullStepKE useFullStepKE_;
    //! Whether we are reporting the conserved energy from the previous step
    const ReportPreviousStepConservedEnergy reportPreviousConservedEnergy_;

    //! The number of temperature groups
    const int numTemperatureGroups_;
    //! The coupling time step - simulation time step x nstcouple_
    const double couplingTimeStep_;
    //! Coupling temperature per group
    std::vector<real> referenceTemperature_;
    //! Coupling time per group
    const std::vector<real> couplingTime_;
    //! Number of degrees of freedom per group
    const std::vector<real> numDegreesOfFreedom_;
    //! Work exerted by thermostat per group
    std::vector<double> temperatureCouplingIntegral_;

    //! Current conserved energy contribution
    real conservedEnergyContribution_;
    //! Step of current conserved energy contribution
    Step conservedEnergyContributionStep_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the energy data (for ekindata)
    EnergyData* energyData_;

    //! Callback to let propagator know that we updated lambda
    PropagatorCallback propagatorCallback_;

    //! Set new lambda value (at T-coupling steps)
    void setLambda(Step step);
    //! Contribution to the conserved energy
    [[nodiscard]] real conservedEnergyContribution() const;

    //! The temperature coupling implementation
    std::unique_ptr<ITemperatureCouplingImpl> temperatureCouplingImpl_;

    //! CheckpointHelper identifier
    const std::string identifier_ = "VelocityScalingTemperatureCoupling";
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);

    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;
    //! The next communicated energy calculation step
    Step nextEnergyCalculationStep_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_VELOCITYSCALINGTEMPERATURECOUPLING_H
