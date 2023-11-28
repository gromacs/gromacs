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
 * \brief
 * Implements QMMM class that implements IMDModule interface
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "qmmm.h"

#include <memory>
#include <numeric>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/logger.h"

#include "qmmmforceprovider.h"
#include "qmmmoptions.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Helper class that holds simulation data and
 * callback functions for simulation setup time notifications
 */
class QMMMSimulationParameterSetup
{
public:
    QMMMSimulationParameterSetup() = default;

    /*! \brief Set the local atom set for the QM atoms.
     * \param[in] localAtomSet of QM atoms
     */
    void setLocalQMAtomSet(const LocalAtomSet& localAtomSet)
    {
        localQMAtomSet_ = std::make_unique<LocalAtomSet>(localAtomSet);
    }

    /*! \brief Set the local atom set for the MM atoms.
     * \param[in] localAtomSet of MM atoms
     */
    void setLocalMMAtomSet(const LocalAtomSet& localAtomSet)
    {
        localMMAtomSet_ = std::make_unique<LocalAtomSet>(localAtomSet);
    }

    /*! \brief Return local atom set for QM atoms.
     * \throws InternalError if local atom set is not set
     * \returns local atom set for QM atoms
     */
    const LocalAtomSet& localQMAtomSet() const
    {
        if (localQMAtomSet_ == nullptr)
        {
            GMX_THROW(InternalError("Local atom set is not set for QM atoms."));
        }
        return *localQMAtomSet_;
    }

    /*! \brief Return local atom set for MM atoms.
     * \throws InternalError if local atom set is not set
     * \returns local atom set for MM atoms
     */
    const LocalAtomSet& localMMAtomSet() const
    {
        if (localMMAtomSet_ == nullptr)
        {
            GMX_THROW(InternalError("Local atom set is not set for MM atoms."));
        }
        return *localMMAtomSet_;
    }

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
            GMX_THROW(
                    InternalError("Periodic boundary condition enum not set for QMMM simulation."));
        }
        return *pbcType_;
    }

    /*! \brief Set the logger for QMMM during mdrun
     * \param[in] logger Logger instance to be used for output
     */
    void setLogger(const MDLogger& logger) { logger_ = &logger; }

    //! Get the logger instance
    const MDLogger& logger() const
    {
        if (logger_ == nullptr)
        {
            GMX_THROW(InternalError("Logger not set for QMMM simulation."));
        }
        return *logger_;
    }

private:
    //! The local QM atom set to act on
    std::unique_ptr<LocalAtomSet> localQMAtomSet_;
    //! The local MM atom set to act on
    std::unique_ptr<LocalAtomSet> localMMAtomSet_;
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;
    //! MDLogger for notifications during mdrun
    const MDLogger* logger_ = nullptr;

    GMX_DISALLOW_COPY_AND_ASSIGN(QMMMSimulationParameterSetup);
};


/*! \internal
 * \brief Handle file output for QMMM simulations.
 * empty implementation as QMMM does not use that
 */
class QMMMOutputProvider final : public IMDOutputProvider
{
public:
    //! Initialize output
    void initOutput(FILE* /*fplog*/,
                    int /*nfile*/,
                    const t_filenm /*fnm*/[],
                    bool /*bAppendFiles*/,
                    const gmx_output_env_t* /*oenv*/) override
    {
    }
    //! Finalizes output from a simulation run.
    void finishOutput() override {}
};


/*! \internal
 * \brief QMMM module
 *
 * Class that implements the QMMM forces and potential
 * \note the virial calculation is not yet implemented
 */
class QMMM final : public IMDModule
{
public:
    //! \brief Construct the QMMM module.
    explicit QMMM() = default;

    // Now callbacks for several kinds of MdModuleNotification are created
    // and subscribed, and will be dispatched correctly at run time
    // based on the type of the parameter required by the lambda.

    /*! \brief Requests to be notified during pre-processing.
     *
     * \param[in] notifier allows the module to subscribe to notifications from MdModules.
     *
     * The QMMM code subscribes to these notifications:
     *   - setting atom group indices in the qmmmOptions_ from an
     *     index group string by taking a parmeter const IndexGroupsAndNames &
     *   - storing its internal parameters in a tpr file by writing to a
     *     key-value-tree during pre-processing by a function taking a
     *     KeyValueTreeObjectBuilder as parameter
     *   - Modify topology according to QMMM rules using gmx_mtop_t notification
     *     and utilizing QMMMTopologyPreprocessor class
     *   - Access MDLogger for notifications output
     *   - Access warninp for for grompp warnings output
     *   - Coordinates, PBC and box for CP2K input generation
     *   - QM Input file provided with -qmi option of grompp
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifier) override
    {
        if (!qmmmOptions_.active())
        {
            return;
        }

        // Writing internal parameters during pre-processing
        const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder) {
            qmmmOptions_.writeInternalParametersToKvt(treeBuilder);
        };
        notifier->preProcessingNotifier_.subscribe(writeInternalParametersFunction);

        // Setting atom group indices
        const auto setQMMMGroupIndicesFunction = [this](const IndexGroupsAndNames& indexGroupsAndNames) {
            qmmmOptions_.setQMMMGroupIndices(indexGroupsAndNames);
        };
        notifier->preProcessingNotifier_.subscribe(setQMMMGroupIndicesFunction);

        // Set Logger during pre-processing
        const auto setLoggerFunction = [this](const MDLogger& logger) {
            qmmmOptions_.setLogger(logger);
        };
        notifier->preProcessingNotifier_.subscribe(setLoggerFunction);

        // Set warning output during pre-processing
        const auto setWarninpFunction = [this](WarningHandler* wi) { qmmmOptions_.setWarninp(wi); };
        notifier->preProcessingNotifier_.subscribe(setWarninpFunction);

        // Notification of the Coordinates, box and pbc during pre-processing
        const auto processCoordinatesFunction = [this](const CoordinatesAndBoxPreprocessed& coord) {
            qmmmOptions_.processCoordinates(coord);
        };
        notifier->preProcessingNotifier_.subscribe(processCoordinatesFunction);

        // Modification of the topology during pre-processing
        const auto modifyQMMMTopologyFunction = [this](gmx_mtop_t* mtop) {
            qmmmOptions_.modifyQMMMTopology(mtop);
        };
        notifier->preProcessingNotifier_.subscribe(modifyQMMMTopologyFunction);

        // Notification of the QM input file provided via -qmi option of grompp
        const auto setQMExternalInputFileNameFunction = [this](const QMInputFileName& qmInputFileName) {
            qmmmOptions_.setQMExternalInputFile(qmInputFileName);
        };
        notifier->preProcessingNotifier_.subscribe(setQMExternalInputFileNameFunction);
    }

    /*! \brief Requests to be notified during simulation setup.
     * The QMMM code subscribes to these notifications:
     *   - reading its internal parameters from a key-value-tree during
     *     simulation setup by taking a const KeyValueTreeObject & parameter
     *   - *.tpr filename for CP2K input generation
     *   - constructing local atom sets in the simulation parameter setup
     *     by taking a LocalAtomSetManager * as parameter
     *   - the type of periodic boundary conditions that are used
     *     by taking a PeriodicBoundaryConditionType as parameter
     *   - Access MDLogger for notifications output
     *   - Disable PME-only ranks for QMMM runs
     *   - Request QM energy output to md.log
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifier) override
    {
        if (!qmmmOptions_.active())
        {
            return;
        }

        // Reading internal parameters during simulation setup
        const auto readInternalParametersFunction = [this](const KeyValueTreeObject& tree) {
            qmmmOptions_.readInternalParametersFromKvt(tree);
        };
        notifier->simulationSetupNotifier_.subscribe(readInternalParametersFunction);

        // Process tpr filename
        const auto setTprFileNameFunction = [this](const MdRunInputFilename& tprName) {
            qmmmOptions_.processTprFilename(tprName);
        };
        notifier->simulationSetupNotifier_.subscribe(setTprFileNameFunction);

        // constructing local atom sets during simulation setup
        const auto setLocalAtomSetFunction = [this](LocalAtomSetManager* localAtomSetManager) {
            LocalAtomSet atomSet1 = localAtomSetManager->add(qmmmOptions_.parameters().qmIndices_);
            this->qmmmSimulationParameters_.setLocalQMAtomSet(atomSet1);
            LocalAtomSet atomSet2 = localAtomSetManager->add(qmmmOptions_.parameters().mmIndices_);
            this->qmmmSimulationParameters_.setLocalMMAtomSet(atomSet2);
        };
        notifier->simulationSetupNotifier_.subscribe(setLocalAtomSetFunction);

        // Reading PBC parameters during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc) {
            this->qmmmSimulationParameters_.setPeriodicBoundaryConditionType(pbc);
        };
        notifier->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

        // Saving MDLogger during simulation setup
        const auto setLoggerFunction = [this](const MDLogger& logger) {
            this->qmmmSimulationParameters_.setLogger(logger);
        };
        notifier->simulationSetupNotifier_.subscribe(setLoggerFunction);

        // Adding output to energy file
        const auto requestEnergyOutput = [](MDModulesEnergyOutputToQMMMRequestChecker* energyOutputRequest) {
            energyOutputRequest->energyOutputToQMMM_ = true;
        };
        notifier->simulationSetupNotifier_.subscribe(requestEnergyOutput);

        // Request to disable PME-only ranks, which are not compatible with CP2K
        const auto requestPmeRanks = [](SeparatePmeRanksPermitted* pmeRanksPermitted) {
            pmeRanksPermitted->disablePmeRanks(
                    "Separate PME-only ranks are not compatible with QMMM MdModule");
        };
        notifier->simulationSetupNotifier_.subscribe(requestPmeRanks);
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &qmmmOptions_; }

    //! Add this module to the force providers if active
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (!qmmmOptions_.active())
        {
            return;
        }

        const auto& parameters = qmmmOptions_.parameters();
        forceProvider_         = std::make_unique<QMMMForceProvider>(
                parameters,
                qmmmSimulationParameters_.localQMAtomSet(),
                qmmmSimulationParameters_.localMMAtomSet(),
                qmmmSimulationParameters_.periodicBoundaryConditionType(),
                qmmmSimulationParameters_.logger());
        forceProviders->addForceProvider(forceProvider_.get());
    }

    //! QMMM Module should not use OutputProvider as it will be removed in the future
    IMDOutputProvider* outputProvider() override { return &qmmmOutputProvider_; }

private:
    //! The output provider
    QMMMOutputProvider qmmmOutputProvider_;
    //! The options provided for QMMM
    QMMMOptions qmmmOptions_;
    //! Object that evaluates the forces
    std::unique_ptr<QMMMForceProvider> forceProvider_;
    /*! \brief Parameters for QMMM that become available at
     * simulation setup time.
     */
    QMMMSimulationParameterSetup qmmmSimulationParameters_;

    GMX_DISALLOW_COPY_AND_ASSIGN(QMMM);
};

} // namespace

std::unique_ptr<IMDModule> QMMMModuleInfo::create()
{
    return std::make_unique<QMMM>();
}

const std::string QMMMModuleInfo::name_ = c_qmmmCP2KModuleName;

} // namespace gmx
