/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief
 * Implements NNPot MDModule class
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "nnpot.h"

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "nnpotforceprovider.h"
#include "nnpotoptions.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief NNPot Module
 *
 * Class that implements the NNPot MDModule.
 */
class NNPotMDModule final : public IMDModule
{
public:
    explicit NNPotMDModule() = default;

    /*! \brief Requests to be notified during preprocessing.
     *
     * NNPot module subscribes to the following notifications:
     *   - the atom groups and their names specified in the index file (to specify NN input)
     *     by taking a const IndexGroupsAndNames& as parameter
     *   - the topology of the system, which has to be modified (to remove classical interactions)
     *     by taking a gmx_mtop_t* as parameter
     *   - writing the module parameters to the KVT for storage in the .tpr file
     *     by taking a KeyValueTreeObjectBuilder as parameter
     *   - access the MDLogger to log messages
     *     by taking a const MDLogger& as parameter
     *   - access the WarningHandler to output warnings
     *     by taking a WarningHandler* as parameter
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifier) override
    {
        if (!nnpotOptions_.isActive())
        {
            return;
        }

        const auto setInputGroupIndicesFunction = [this](const IndexGroupsAndNames& indexGroupsAndNames)
        { nnpotOptions_.setInputGroupIndices(indexGroupsAndNames); };
        notifier->preProcessingNotifier_.subscribe(setInputGroupIndicesFunction);

        const auto modifyTopologyFunction = [this](gmx_mtop_t* top)
        { nnpotOptions_.modifyTopology(top); };
        notifier->preProcessingNotifier_.subscribe(modifyTopologyFunction);

        const auto writeParamsToKvtFunction = [this](KeyValueTreeObjectBuilder kvt)
        { nnpotOptions_.writeParamsToKvt(kvt); };
        notifier->preProcessingNotifier_.subscribe(writeParamsToKvtFunction);

        // Set Logger during pre-processing
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { nnpotOptions_.setLogger(logger); };
        notifier->preProcessingNotifier_.subscribe(setLoggerFunction);

        // Set warning output during pre-processing
        const auto setWarninpFunction = [this](WarningHandler* wi) { nnpotOptions_.setWarninp(wi); };
        notifier->preProcessingNotifier_.subscribe(setWarninpFunction);
    }

    /*! \brief Requests to be notified during simulation setup.
     *
     * NNPot module subscribes to the following notifications:
     *   - the topology of the system
     *     by taking a const gmx_mtop_t& as parameter
     *   - the local atom set manager to construct local atom sets for NN input and its complement
     *     by taking a LocalAtomSetManager* as parameter
     *   - reading the module parameters from the KVT
     *     by taking a const KeyValueTreeObject& as parameter
     *   - the PBC type
     *     by taking a PbcType as parameter
     *   - access the MDLogger to log messages
     *     by taking a const MDLogger& as parameter
     *   - access the communication record
     *     by taking a const t_commrec& as parameter
     *   - notify when atoms are redistributed
     *     by taking a const MDModulesAtomsRedistributedSignal as parameter
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifier) override
    {
        if (!nnpotOptions_.isActive())
        {
            return;
        }

        const auto setTopologyFunction = [this](const gmx_mtop_t& top)
        { nnpotOptions_.setTopology(top); };
        notifier->simulationSetupNotifier_.subscribe(setTopologyFunction);

        // constructing local atom sets during simulation setup
        const auto setLocalAtomSetFunction = [this](LocalAtomSetManager* localAtomSetManager)
        {
            LocalAtomSet atomSet1 = localAtomSetManager->add(nnpotOptions_.parameters().inpIndices_);
            nnpotOptions_.setLocalInputAtomSet(atomSet1);
            LocalAtomSet atomSet2 = localAtomSetManager->add(nnpotOptions_.parameters().mmIndices_);
            nnpotOptions_.setLocalMMAtomSet(atomSet2);
        };
        notifier->simulationSetupNotifier_.subscribe(setLocalAtomSetFunction);

        const auto readParamsFromKvtFunction = [this](const KeyValueTreeObject& kvt)
        { nnpotOptions_.readParamsFromKvt(kvt); };
        notifier->simulationSetupNotifier_.subscribe(readParamsFromKvtFunction);

        const auto setPBCTypeFunction = [this](const PbcType& pbc) { nnpotOptions_.setPbcType(pbc); };
        notifier->simulationSetupNotifier_.subscribe(setPBCTypeFunction);

        // Add NN output to energy file
        const auto requestEnergyOutput = [](MDModulesEnergyOutputToNNPotRequestChecker* energyOutputRequest)
        { energyOutputRequest->energyOutputToNNPot_ = true; };
        notifier->simulationSetupNotifier_.subscribe(requestEnergyOutput);

        // Set Logger during simulation setup
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { nnpotOptions_.setLogger(logger); };
        notifier->simulationSetupNotifier_.subscribe(setLoggerFunction);

        // set communication record during simulation setup
        const auto setCommRecFunction = [this](const t_commrec& cr) { nnpotOptions_.setCommRec(cr); };
        notifier->simulationSetupNotifier_.subscribe(setCommRecFunction);

        // subscribe to DD notification to trigger atom number and index gathering
        const auto notifyDDFunction = [this](const MDModulesAtomsRedistributedSignal& /*signal*/)
        { nnpotForceProvider_->gatherAtomNumbersIndices(); };
        notifier->simulationSetupNotifier_.subscribe(notifyDDFunction);
    }

    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (!nnpotOptions_.isActive())
        {
            return;
        }

        nnpotForceProvider_ = std::make_unique<NNPotForceProvider>(nnpotOptions_.parameters(),
                                                                   nnpotOptions_.logger());
        forceProviders->addForceProvider(nnpotForceProvider_.get(), "NN potential");
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &nnpotOptions_; }
    //! From IMDModule
    //! NNPot doesn't need extra output
    IMDOutputProvider* outputProvider() override { return nullptr; }

private:
    NNPotOptions nnpotOptions_;

    std::unique_ptr<NNPotForceProvider> nnpotForceProvider_;
};

} // end namespace

std::unique_ptr<IMDModule> NNPotModuleInfo::create()
{
    return std::make_unique<NNPotMDModule>();
}

const std::string NNPotModuleInfo::name_ = "nnpot";

} // end namespace gmx
