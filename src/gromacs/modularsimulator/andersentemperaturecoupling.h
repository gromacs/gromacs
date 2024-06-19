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
 * \brief Declares Andersen temperature coupling for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 */

#ifndef GMX_MODULARSIMULATOR_ANDERSENTHERMOSTAT_H
#define GMX_MODULARSIMULATOR_ANDERSENTHERMOSTAT_H

#include <cstdint>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "energydata.h"
#include "modularsimulatorinterfaces.h"
#include "propagator.h"

struct t_commrec;
struct t_mdatoms;

namespace gmx
{
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;
enum class ReferenceTemperatureChangeAlgorithm;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Element implementing the Andersen thermostat
 *
 */
class AndersenTemperatureCoupling final : public ISimulatorElement
{
public:
    //! Constructor
    AndersenTemperatureCoupling(double               simulationTimestep,
                                bool                 doMassive,
                                int64_t              seed,
                                ArrayRef<const real> referenceTemperature,
                                ArrayRef<const real> couplingTime,
                                StatePropagatorData* statePropagatorData,
                                const MDAtoms*       mdAtoms,
                                const t_commrec*     cr);

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

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer);

    //! Returns the frequency at which temperature coupling is performed
    [[nodiscard]] int frequency() const;

private:
    //! Update the reference temperature
    static void updateReferenceTemperature(ArrayRef<const real>                temperatures,
                                           ReferenceTemperatureChangeAlgorithm algorithm);

    //! Whether we're doing massive Andersen thermostatting
    const bool doMassive_;
    //! The rate
    const real randomizationRate_;
    //! The frequency at which the thermostat is applied
    const int couplingFrequency_;
    //! The random seed
    const int64_t seed_;
    //! Coupling temperature per group
    ArrayRef<const real> referenceTemperature_;
    //! Coupling time per group
    ArrayRef<const real> couplingTime_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Atom parameters for this domain.
    const t_mdatoms* mdAtoms_;
    //! Handles communication.
    const t_commrec* cr_;

    //! Apply the thermostat at step
    void apply(Step step);
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_ANDERSENTHERMOSTAT_H
