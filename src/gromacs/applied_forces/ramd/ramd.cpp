/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Implements RAMD class that implements IMDModule interface
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "ramd.h"

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/logger.h"

#include "ramdforceprovider.h"
#include "ramdoptions.h"
#include "ramdoutputprovider.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Helper class that holds simulation data and
 * callback functions for simulation setup time notifications
 */
class RAMDSimulationParameterSetup
{
public:
    RAMDSimulationParameterSetup() = default;

    /*! \brief Set the periodic boundary condition via MdModuleNotifier.
     *
     * The pbc type is wrapped in PeriodicBoundaryConditionType to
     * allow the MdModuleNotifier to statically distinguish the callback
     * function type from other 'int' function callbacks.
     *
     * \param[in] pbc MdModuleNotification class that contains a variable
     *                that enumerates the periodic boundary condition.
     */
    void setPeriodicBoundaryConditionType(const PbcType& pbc)
    {
        pbcType_ = std::make_unique<PbcType>(pbc);
    }

    //! Get the periodic boundary conditions
    PbcType periodicBoundaryConditionType()
    {
        if (pbcType_ == nullptr)
        {
            GMX_THROW(InternalError("Periodic boundary condition enum not set for RAMD."));
        }
        return *pbcType_;
    }

    /*! \brief Set the logger for RAMD during mdrun
     * \param[in] logger Logger instance to be used for output
     */
    void setLogger(const MDLogger& logger) { logger_ = &logger; }

    //! Get the logger instance
    const MDLogger& logger() const
    {
        GMX_RELEASE_ASSERT(logger_, "Logger not set for RAMD.");
        return *logger_;
    }

    //! Set the local atom sets
    void setLocalAtomSets(const LocalAtomSet& localAtomSet)
    {
        localAtomSets_.emplace_back(std::make_unique<LocalAtomSet>(localAtomSet));
    }

    /*! \brief Return local atom sets
     * \throws InternalError if local atom set is not set
     */
    const std::vector<std::unique_ptr<LocalAtomSet>>& localAtomSets() const
    {
        // if (localAtomSets_.empty())
        // {
        //     GMX_THROW(InternalError("Local atom sets are not set for RAMD."));
        // }
        return localAtomSets_;
    }

    /*! \brief Set the topology for RAMD during mdrun
     * \param[in] top The topology of the system
     */
    void setTopology(const gmx_mtop_t& top) { top_ = &top; }

    /*! \brief Get the topology of the system
     * \throws InternalError if topology is not set
     */
    const gmx_mtop_t& topology() const
    {
        if (top_ == nullptr)
        {
            GMX_THROW(InternalError("Topology not set for RAMD."));
        }
        return *top_;
    }

private:
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;

    /*! \brief MDLogger during mdrun
     *
     * This is a pointer only because we need an "optional reference"
     * to a const MDLogger before the notification always provides the
     * actual reference. */
    const MDLogger* logger_ = nullptr;

    //! The local atom sets to act on
    std::vector<std::unique_ptr<LocalAtomSet>> localAtomSets_;

    //! The topology of the system
    const gmx_mtop_t* top_ = nullptr;

    GMX_DISALLOW_COPY_AND_ASSIGN(RAMDSimulationParameterSetup);
};


/*! \internal
 * \brief RAMD module
 */
class RAMD final : public IMDModule
{
public:
    //! \brief Construct the RAMD module.
    explicit RAMD() = default;

    //! Returns an interface for handling mdp input (and tpr I/O).
    IMdpOptionProvider* mdpOptionProvider() override { return &ramdOptions_; }

    //! Returns an interface for handling output files during simulation.
    IMDOutputProvider* outputProvider() override { return &ramdOutputProvider_; }

    //! Initializes force providers from this module.
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (ramdOptions_.active())
        {
            const auto& parameters = ramdOptions_.parameters();
            forceProvider_         = std::make_unique<RAMDForceProvider>(
                    parameters,
                    ramdSimulationParameters_.localAtomSets(),
                    ramdSimulationParameters_.topology(),
                    ramdSimulationParameters_.periodicBoundaryConditionType(),
                    ramdSimulationParameters_.logger(),
                    ramdOutputProvider_);
            forceProviders->addForceProvider(forceProvider_.get(), "RAMD");
        }
    }

    //! Subscribe to pre processing notifications
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override
    {
        if (!ramdOptions_.active())
        {
            return;
        }

        // Set input group indices during pre-processing
        const auto setInputGroupIndicesFunction = [this](const IndexGroupsAndNames& indexGroupsAndNames)
        { ramdOptions_.setInputGroupIndices(indexGroupsAndNames); };
        notifiers->preProcessingNotifier_.subscribe(setInputGroupIndicesFunction);

        // Writing internal parameters during pre-processing
        const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder)
        { ramdOptions_.writeInternalParametersToKvt(treeBuilder); };
        notifiers->preProcessingNotifier_.subscribe(writeInternalParametersFunction);

        // Set Logger during pre-processing
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { ramdOptions_.setLogger(logger); };
        notifiers->preProcessingNotifier_.subscribe(setLoggerFunction);
    }

    //! Subscribe to simulation setup notifications
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
    {
        if (!ramdOptions_.active())
        {
            return;
        }

        // Reading internal parameters during simulation setup
        const auto readInternalParametersFunction = [this](const KeyValueTreeObject& tree)
        { ramdOptions_.readInternalParametersFromKvt(tree); };
        notifiers->simulationSetupNotifier_.subscribe(readInternalParametersFunction);

        // Constructing local atom sets during simulation setup
        const auto setLocalAtomSetFunction = [this](LocalAtomSetManager* localAtomSetManager)
        {
            for (int g = 0; g < ramdOptions_.parameters().ngroups_; ++g)
            {
                LocalAtomSet atomSet =
                        localAtomSetManager->add(ramdOptions_.parameters().groups_[g].ligand_indices_);
                this->ramdSimulationParameters_.setLocalAtomSets(atomSet);
            }
        };
        notifiers->simulationSetupNotifier_.subscribe(setLocalAtomSetFunction);

        // Reading topology during simulation setup
        const auto setTopologyFunction = [this](const gmx_mtop_t& top)
        { this->ramdSimulationParameters_.setTopology(top); };
        notifiers->simulationSetupNotifier_.subscribe(setTopologyFunction);

        // Reading PBC parameters during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc)
        { this->ramdSimulationParameters_.setPeriodicBoundaryConditionType(pbc); };
        notifiers->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

        // Saving MDLogger during simulation setup
        const auto setLoggerFunction = [this](const MDLogger& logger)
        { this->ramdSimulationParameters_.setLogger(logger); };
        notifiers->simulationSetupNotifier_.subscribe(setLoggerFunction);
    }

    //! Subscribe to simulation run notifications
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

private:
    //! The output provider
    RAMDOutputProvider ramdOutputProvider_;

    //! The options provided for RAMD
    RAMDOptions ramdOptions_;

    //! Object that evaluates the forces
    std::unique_ptr<RAMDForceProvider> forceProvider_;

    //! RAMD Parameters that become available at simulation setup time.
    RAMDSimulationParameterSetup ramdSimulationParameters_;

    GMX_DISALLOW_COPY_AND_ASSIGN(RAMD);
};

} // namespace

std::unique_ptr<IMDModule> RAMDModuleInfo::create()
{
    return std::make_unique<RAMD>();
}

} // namespace gmx
