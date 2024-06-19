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
 * \brief Declares the free energy perturbation element for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_FREEENERGYPERTURBATIONELEMENT_H
#define GMX_MODULARSIMULATOR_FREEENERGYPERTURBATIONELEMENT_H

#include <cstdio>

#include <memory>
#include <optional>
#include <string>

#include "gromacs/mdtypes/checkpointdata.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"

class gmx_ekindata_t;
struct t_inputrec;
struct t_trxframe;
struct t_commrec;

namespace gmx
{
enum class CheckpointDataOperation;
class EnergyData;
class FepStateSetting;
class GlobalCommunicationHelper;
class LegacySimulatorData;
class MDAtoms;
class ModularSimulatorAlgorithmBuilderHelper;
class StatePropagatorData;
template<CheckpointDataOperation operation>
class CheckpointData;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief The free energy perturbation data
 *
 * The lambda vector and the current FEP state are held by the
 * FreeEnergyPerturbationData, offering access to its values via getter
 * functions. The FreeEnergyPerturbationData::Element is responsible for
 * lambda update (if applicable) and checkpointing.
 */
class FreeEnergyPerturbationData final
{
public:
    //! Constructor
    FreeEnergyPerturbationData(FILE* fplog, const t_inputrec& inputrec, MDAtoms* mdAtoms, gmx_ekindata_t* ekindata);

    //! Get a view of the current lambda vector
    ArrayRef<real> lambdaView();
    //! Get a const view of the current lambda vector
    [[nodiscard]] ArrayRef<const real> constLambdaView() const;
    //! Get the current FEP state
    [[nodiscard]] int currentFEPState() const;

    /*! \brief Enable setting of the FEP state by an external object
     *
     * Currently, this can only be called once, usually during setup time.
     * Having more than one object setting the FEP state would require additional bookkeeping.
     *
     * \return Pointer to an object allowing to set new FEP state
     */
    [[nodiscard]] FepStateSetting* enableExternalFepStateSetting() const;

    //! The element taking part in the simulator loop
    class Element;
    //! Get pointer to element (whose lifetime is managed by this)
    Element* element();

    //! Read everything that can be stored in t_trxframe from a checkpoint file
    static void readCheckpointToTrxFrame(t_trxframe*                       trxFrame,
                                         std::optional<ReadCheckpointData> readCheckpointData);
    //! CheckpointHelper identifier
    static const std::string& checkpointID();

private:
    //! Update the lambda values
    void updateLambdas(Step step);
    //! Update the lambda values
    void setLambdaState(Step step, int newState);
    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);
    //! Update MDAtoms
    void updateMDAtoms();

    //! The element
    std::unique_ptr<Element> element_;

    //! The lambda vector
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real> lambda_;
    //! The current free energy state
    int currentFEPState_;

    //! Handles logging.
    FILE* fplog_;
    //! Contains user input mdp options.
    const t_inputrec& inputrec_;
    //! Atom parameters for this domain.
    MDAtoms* mdAtoms_;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Allows external clients to specify how to change the FEP state
 */
class FepStateSetting
{
public:
    //! Signal (during task scheduling) that a signal stepping step will happen
    void signalSettingStep(Step step);
    //! Set new state at specific step (called during simulation run)
    void setNewState(int state, Step step);

    // Allow private member access
    friend class FreeEnergyPerturbationData::Element;

private:
    //! The next external lambda setting step
    Step nextFepStateSettingStep = -1;
    //! The new FEP state set externally
    int newFepState = -1;
    //! The step at which the new FEP state gets used
    Step newFepStateStep = -1;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief The free energy perturbation data element
 *
 * The FreeEnergyPerturbationData::Element does update the lambda
 * values during the simulation run if lambda is non-static. It does
 * implement the checkpointing client interface to save its current
 * state for restart.
 */
class FreeEnergyPerturbationData::Element final :
    public ISimulatorElement,
    public ICheckpointHelperClient,
    public IDomDecHelperClient
{
public:
    //! Constructor
    explicit Element(FreeEnergyPerturbationData* freeEnergyPerturbationElement, double deltaLambda);

    //! Update lambda and mdatoms
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! Update the MdAtoms object
    void elementSetup() override;

    //! No teardown needed
    void elementTeardown() override{};

    //! ICheckpointHelperClient write checkpoint implementation
    void saveCheckpointState(std::optional<WriteCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient read checkpoint implementation
    void restoreCheckpointState(std::optional<ReadCheckpointData> checkpointData, const t_commrec* cr) override;
    //! ICheckpointHelperClient key implementation
    const std::string& clientID() override;

    //! Callback on domain decomposition repartitioning
    DomDecCallback registerDomDecCallback() override;

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper  Pointer to the \c GlobalCommunicationHelper object
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper);

    //! Enable setting of the FEP state by an external object
    FepStateSetting* enableExternalFepStateSetting();

private:
    //! The free energy data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;
    //! Whether lambda values change continuously
    const bool doSlowGrowth_;

    //! Information about external lambda setting, set only if external lambda setting is enabled
    std::optional<FepStateSetting> externalFepStateSetting_;
    //! The number of external lambda setting clients
    int numExternalFepStateSettingClients_;

    //! Helper function to read from / write to CheckpointData
    template<CheckpointDataOperation operation>
    void doCheckpointData(CheckpointData<operation>* checkpointData);
    //! Whether the element was restored from checkpoint
    bool restoredFromCheckpoint_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_FREEENERGYPERTURBATIONELEMENT_H
