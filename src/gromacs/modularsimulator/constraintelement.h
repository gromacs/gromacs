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
 * \brief Declares the constraint element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_CONSTRAINTELEMENT_H
#define GMX_MODULARSIMULATOR_CONSTRAINTELEMENT_H

#include "gromacs/mdlib/constr.h"

#include "modularsimulatorinterfaces.h"

namespace gmx
{
class Constraints;
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Constraints element
 *
 * The ConstraintsElement is implemented for the position-and-velocity and the
 * velocity-only case. It does not change the constraint implementation itself,
 * but uses the current constraints implementation and the data management
 * introduced with the modular simulator.
 *
 * \tparam variable  The constraining variable
 */
template<ConstraintVariable variable>
class ConstraintsElement final :
    public ISimulatorElement,
    public IEnergySignallerClient,
    public ITrajectorySignallerClient,
    public ILoggingSignallerClient
{
public:
    //! Constructor
    ConstraintsElement(Constraints*                constr,
                       StatePropagatorData*        statePropagatorData,
                       EnergyData*                 energyData,
                       FreeEnergyPerturbationData* freeEnergyPerturbationData,
                       bool                        isMain,
                       FILE*                       fplog,
                       const t_inputrec*           inputrec,
                       const t_mdatoms*            mdAtoms);

    /*! \brief Register constraining function for step / time
     *
     * Under LF, this is expected to be run once, constraining positions and velocities
     * Under VV, this is expected to be run twice, once contraining velocities only,
     * a second time constraining positions and velocities.
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    /*! \brief Performs inital constraining
     *  \todo Should this rather happen at grompp time? Right position of this operation is currently
     *        depending on the integrator algorithm (after domdec, before compute globals...),
     *        so doing this earlier would be much more stable!
     */
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

private:
    //! The actual constraining computation
    void apply(Step step, bool calculateVirial, bool writeLog, bool writeEnergy);

    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;
    //! ITrajectorySignallerClient implementation
    std::optional<SignallerCallback> registerTrajectorySignallerCallback(TrajectoryEvent event) override;
    //! ILoggingSignallerClient implementation
    std::optional<SignallerCallback> registerLoggingCallback() override;

    //! The next energy calculation step
    Step nextVirialCalculationStep_;
    //! The next energy writing step
    Step nextEnergyWritingStep_;
    //! The next logging step
    Step nextLogWritingStep_;

    //! Whether we're main rank
    const bool isMainRank_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy data
    EnergyData* energyData_;
    //! Pointer to the free energy perturbation data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    // Access to ISimulator data
    //! Handles constraints.
    Constraints* constr_;
    //! Handles logging.
    FILE* fplog_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Atom parameters for this domain.
    const t_mdatoms* mdAtoms_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_CONSTRAINTELEMENT_H
