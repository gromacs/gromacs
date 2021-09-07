/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief Declares the element performing first-order pressure coupling for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_FIRSTORDERPRESSURECOUPLING_H
#define GMX_MODULARSIMULATOR_FIRSTORDERPRESSURECOUPLING_H

#include "modularsimulatorinterfaces.h"

struct t_inputrec;
struct t_nrnb;

enum class PressureCoupling;

namespace gmx
{
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;

/*! \internal
 * \brief Element performing first-order pressure coupling
 *
 * This element implements the first-order pressure coupling algorithms
 * (Berendesen, c-rescale).
 */
class FirstOrderPressureCoupling final : public ISimulatorElement, public ICheckpointHelperClient
{
public:
    //! Constructor
    FirstOrderPressureCoupling(int                               couplingFrequency,
                               int                               couplingOffset,
                               real                              couplingTimeStep,
                               StatePropagatorData*              statePropagatorData,
                               EnergyData*                       energyData,
                               FILE*                             fplog,
                               const t_inputrec*                 inputrec,
                               const MDAtoms*                    mdAtoms,
                               t_nrnb*                           nrnb,
                               ReportPreviousStepConservedEnergy reportPreviousStepConservedEnergy);

    void scheduleTask(Step step, Time time, const RegisterRunFunction& function) override;
    //! Setup - initialize relative box matrix
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
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     * \param offset  The step offset at which the barostat is applied
     * \param reportPreviousStepConservedEnergy  Report the previous or the current step conserved energy
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement*
    getElementPointerImpl(LegacySimulatorData*                    legacySimulatorData,
                          ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                          StatePropagatorData*                    statePropagatorData,
                          EnergyData*                             energyData,
                          FreeEnergyPerturbationData gmx_unused* freeEnergyPerturbationData,
                          GlobalCommunicationHelper gmx_unused* globalCommunicationHelper,
                          ObservablesReducer*                   observablesReducer,
                          int                                   offset,
                          ReportPreviousStepConservedEnergy     reportPreviousStepConservedEnergy);

private:
    //! Calculate the scaling matrix
    template<PressureCoupling pressureCouplingType>
    void calculateScalingMatrix(Step step);
    //! Scale the box and coordinates according to the current scaling matrix
    template<PressureCoupling pressureCouplingType>
    void scaleBoxAndCoordinates();
    //! Helper function returning the conserved energy contribution
    real conservedEnergyContribution(Step step);

    //! The pressure coupling type required
    const PressureCoupling pressureCouplingType_;
    //! The coupling time step (simulation time step x coupling frequency)
    const real couplingTimeStep_;
    //! The frequency at which pressure coupling is applied
    const int couplingFrequency_;
    //! The offset at which pressure coupling is applied
    const int couplingOffset_;

    //! The current box scaling matrix
    tensor boxScalingMatrix_;
    //! Relative box shape
    tensor boxRel_;
    //! Contribution to the conserved energy
    double conservedEnergyContribution_;
    //! Contribution to the conserved energy
    double previousStepConservedEnergyContribution_;
    //! Current step of the conserved energy contribution
    Step conservedEnergyContributionStep_;
    //! Whether we're reporting current step or previous step conserved energy
    const ReportPreviousStepConservedEnergy reportPreviousStepConservedEnergy_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy data
    EnergyData* energyData_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;

    //! CheckpointHelper identifier
    const std::string identifier_;
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_FIRSTORDERPRESSURECOUPLING_H
