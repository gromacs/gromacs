/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
 * \brief Declares the pull element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_PULLELEMENT_H
#define GMX_MODULARSIMULATOR_PULLELEMENT_H

#include "modularsimulatorinterfaces.h"

struct gmx_mtop_t;
struct pull_t;
struct t_inputrec;

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
 * \brief Element calling pull functionality
 */
class PullElement : public ISimulatorElement, public ICheckpointHelperClient
{
public:
    //! Constructor
    PullElement(bool                 setPbcRefToPrevStepCOM,
                PbcType              pbcType,
                StatePropagatorData* statePropagatorData,
                pull_t*              pullWork,
                const t_commrec*     commrec,
                const MDAtoms*       mdAtoms);
    //! Update annealing temperature
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;
    //! Set initial annealing temperature
    void elementSetup() override;
    //! No teardown needed
    void elementTeardown() override {}

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
    //! Schedule post step functionality
    void schedulePostStep(Step step, Time time, const RegisterRunFunction& registerRunFunction);

    //! Whether to use the COM of each group from the previous step as reference
    const bool setPbcRefToPrevStepCOM_;
    //! The PBC type
    const PbcType pbcType_;

    //! CheckpointHelper identifier
    const std::string identifier_ = "PullElement";
    //! Whether this object was restored from checkpoint
    bool restoredFromCheckpoint_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;

    // Access to LegacySimulatorData
    //! The pull work object.
    pull_t* pullWork_;
    //! Handles communication.
    const t_commrec* commrec_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
};
} // namespace gmx


#endif // GMX_MODULARSIMULATOR_PULLELEMENT_H
