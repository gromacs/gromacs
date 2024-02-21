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
 * \brief Declares the expanded ensemble element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_EXPANDEDENSEMBLEELEMENT_H
#define GMX_MODULARSIMULATOR_EXPANDEDENSEMBLEELEMENT_H

#include "freeenergyperturbationdata.h"

struct df_history_t;
struct t_inputrec;

namespace gmx
{
enum class CheckpointDataOperation;
class EnergyData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief The expanded ensemble element
 *
 * This element periodically attempts Monte Carlo moves in lambda
 * space and sets the new lambda state in FreeEnergyPerturbationData::Element.
 */
class ExpandedEnsembleElement final : public ISimulatorElement, public ICheckpointHelperClient, public ILoggingSignallerClient
{
public:
    //! Constructor
    explicit ExpandedEnsembleElement(bool                              isMainRank,
                                     Step                              initialStep,
                                     int                               frequency,
                                     const EnergyData*                 energyData,
                                     const FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                     FILE*                             fplog,
                                     const t_inputrec*                 inputrec);
    //! Destructor
    ~ExpandedEnsembleElement();

    //! Attempt lambda MC step and write log
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;
    //! Set up FEP history object
    void elementSetup() override;
    //! No teardown needed
    void elementTeardown() override{};

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
     * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
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
    //! Use expanded ensemble to determine new FEP state or write log
    void apply(Step step, bool doLambdaStep, bool doLog);

    //! Helper object allowing to set new FEP state
    FepStateSetting* fepStateSetting_;

    //! Whether this runs on main
    const bool isMainRank_;
    //! The initial Step
    const Step initialStep_;
    //! The frequency of lambda MC steps
    const int frequency_;

    //! ILoggingSignallerClient implementation
    std::optional<SignallerCallback> registerLoggingCallback() override;
    //! The next logging step
    Step nextLogWritingStep_;

    //! The free energy sampling history
    std::unique_ptr<df_history_t> dfhist_;

    //! CheckpointHelper identifier
    const std::string identifier_ = "ExpandedEnsembleElement";
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);
    //! Whether this object was restored from checkpoint
    bool restoredFromCheckpoint_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the energy data
    const EnergyData* energyData_;
    //! Pointer to the free energy perturbation data
    const FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_EXPANDEDENSEMBLEELEMENT_H
