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
 * \brief Declares the v-rescale thermostat for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_VRESCALETHERMOSTAT_H
#define GMX_MODULARSIMULATOR_VRESCALETHERMOSTAT_H

#include "gromacs/utility/arrayref.h"

#include "energydata.h"
#include "modularsimulatorinterfaces.h"
#include "propagator.h"

struct t_commrec;

namespace gmx
{
class LegacySimulatorData;

//! Enum describing whether the thermostat is using full or half step kinetic energy
enum class VRescaleThermostatUseFullStepKE
{
    Yes,
    No
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element implementing the v-rescale thermostat
 *
 * This element takes a callback to the propagator and updates the velocity
 * scaling factor according to the v-rescale thermostat.
 */
class VRescaleThermostat final : public ISimulatorElement, public ICheckpointHelperClient
{
public:
    //! Constructor
    VRescaleThermostat(int            nstcouple,
                       int            offset,
                       bool           useFullStepKE,
                       int64_t        seed,
                       int            numTemperatureGroups,
                       double         couplingTimeStep,
                       const real*    referenceTemperature,
                       const real*    couplingTime,
                       const real*    numDegreesOfFreedom,
                       EnergyData*    energyData,
                       const t_state* globalState,
                       t_commrec*     cr,
                       bool           isRestart);

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

    //! Getter for the thermostatIntegral
    const std::vector<double>& thermostatIntegral() const;

    //! Connect this to propagator
    void connectWithPropagator(const PropagatorThermostatConnection& connectionData);

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
     * \param offset  The step offset at which the thermostat is applied
     * \param useFullStepKE  Whether full step or half step KE is used
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    int                        offset,
                                                    VRescaleThermostatUseFullStepKE useFullStepKE);

private:
    //! The frequency at which the thermostat is applied
    const int nstcouple_;
    //! If != 0, offset the step at which the thermostat is applied
    const int offset_;
    //! Whether we're using full step kinetic energy
    const bool useFullStepKE_;
    //! The random seed
    const int64_t seed_;

    //! The number of temperature groups
    const int numTemperatureGroups_;
    //! The coupling time step - simulation time step x nstcouple_
    const double couplingTimeStep_;
    //! Coupling temperature per group
    const std::vector<real> referenceTemperature_;
    //! Coupling time per group
    const std::vector<real> couplingTime_;
    //! Number of degrees of freedom per group
    const std::vector<real> numDegreesOfFreedom_;
    //! Work exerted by thermostat
    std::vector<double> thermostatIntegral_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the energy data (for ekindata)
    EnergyData* energyData_;

    //! View on the scaling factor of the propagator
    ArrayRef<real> lambda_;
    //! Callback to let propagator know that we updated lambda
    PropagatorCallback propagatorCallback_;

    //! Set new lambda value (at T-coupling steps)
    void setLambda(Step step);

    //! ICheckpointHelperClient implementation
    void writeCheckpoint(t_state* localState, t_state* globalState) override;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_VRESCALETHERMOSTAT_H
