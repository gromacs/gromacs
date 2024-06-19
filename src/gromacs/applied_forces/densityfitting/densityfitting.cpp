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
 * \brief
 * Declares data structure and utilities for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfitting.h"

#include <filesystem>
#include <functional>
#include <memory>
#include <numeric>
#include <vector>

#include "gromacs/applied_forces/densityfitting/densityfittingparameters.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/fileio/mrcdensitymap.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "densityfittingforceprovider.h"
#include "densityfittingoptions.h"
#include "densityfittingoutputprovider.h"

enum class PbcType : int;

namespace gmx
{

class IMdpOptionProvider;
class DensityFittingForceProvider;
class IMDOutputProvider;
class KeyValueTreeObject;

namespace
{

/*! \internal
 * \brief Collect density fitting parameters only available during simulation setup.
 *
 * \todo Implement builder pattern that will not use unqiue_ptr to check if
 *       parameters have been set or not.
 *
 * To build the density fitting force provider during simulation setup,
 * the DensityFitting class needs access to parameters that become available
 * only during simulation setup.
 *
 * This class collects these parameters via MDModulesNotifiers in the
 * simulation setup phase and provides a check if all necessary parameters have
 * been provided.
 */
class DensityFittingSimulationParameterSetup
{
public:
    DensityFittingSimulationParameterSetup() = default;

    /*! \brief Set the local atom set for the density fitting.
     * \param[in] localAtomSet of atoms to be fitted
     */
    void setLocalAtomSet(const LocalAtomSet& localAtomSet)
    {
        localAtomSet_ = std::make_unique<LocalAtomSet>(localAtomSet);
    }

    /*! \brief Return local atom set for density fitting.
     * \throws InternalError if local atom set is not set
     * \returns local atom set for density fitting
     */
    const LocalAtomSet& localAtomSet() const
    {
        if (localAtomSet_ == nullptr)
        {
            GMX_THROW(
                    InternalError("Local atom set is not set for density "
                                  "guided simulation."));
        }
        return *localAtomSet_;
    }

    /*! \brief Return transformation into density lattice.
     * \throws InternalError if transformation into density lattice is not set
     * \returns transformation into density lattice
     */
    const TranslateAndScale& transformationToDensityLattice() const
    {
        if (transformationToDensityLattice_ == nullptr)
        {
            GMX_THROW(InternalError(
                    "Transformation to reference density not set for density guided simulation."));
        }
        return *transformationToDensityLattice_;
    }
    /*! \brief Return reference density
     * \throws InternalError if reference density is not set
     * \returns the reference density
     */
    basic_mdspan<const float, dynamicExtents3D> referenceDensity() const
    {
        if (referenceDensity_ == nullptr)
        {
            GMX_THROW(InternalError("Reference density not set for density guided simulation."));
        }
        return referenceDensity_->asConstView();
    }

    /*! \brief Reads the reference density from file.
     *
     * Reads and check file, then set and communicate the internal
     * parameters related to the reference density with the file data.
     *
     * \throws FileIOError if reading from file was not successful
     */
    void readReferenceDensityFromFile(const std::string& referenceDensityFileName)
    {
        MrcDensityMapOfFloatFromFileReader reader(referenceDensityFileName);
        referenceDensity_ = std::make_unique<MultiDimArray<std::vector<float>, dynamicExtents3D>>(
                reader.densityDataCopy());
        transformationToDensityLattice_ =
                std::make_unique<TranslateAndScale>(reader.transformationToDensityLattice());
    }

    //! Normalize the reference density so that the sum over all voxels is unity
    void normalizeReferenceDensity()
    {
        if (referenceDensity_ == nullptr)
        {
            GMX_THROW(InternalError("Need to set reference density before normalizing it."));
        }

        normalizeSumPositiveValuesToUnity(referenceDensity_->toArrayRef());
    }
    /*! \brief Set the periodic boundary condition via MDModuleNotifier.
     *
     * The pbc type is wrapped in PeriodicBoundaryConditionType to
     * allow the MDModuleNotifier to statically distinguish the callback
     * function type from other 'int' function callbacks.
     *
     * \param[in] pbcType enumerates the periodic boundary condition.
     */
    void setPeriodicBoundaryConditionType(const PbcType& pbcType)
    {
        pbcType_ = std::make_unique<PbcType>(pbcType);
    }

    //! Get the periodic boundary conditions
    PbcType periodicBoundaryConditionType()
    {
        if (pbcType_ == nullptr)
        {
            GMX_THROW(InternalError(
                    "Periodic boundary condition enum not set for density guided simulation."));
        }
        return *pbcType_;
    }

    //! Set the simulation time step
    void setSimulationTimeStep(double timeStep) { simulationTimeStep_ = timeStep; }

    //! Return the simulation time step
    double simulationTimeStep() const { return simulationTimeStep_; }

private:
    //! The reference density to fit to
    std::unique_ptr<MultiDimArray<std::vector<float>, dynamicExtents3D>> referenceDensity_;
    //! The coordinate transformation into the reference density
    std::unique_ptr<TranslateAndScale> transformationToDensityLattice_;
    //! The local atom set to act on
    std::unique_ptr<LocalAtomSet> localAtomSet_;
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;
    //! The simulation time step
    double simulationTimeStep_ = 1;

    GMX_DISALLOW_COPY_AND_ASSIGN(DensityFittingSimulationParameterSetup);
};

/*! \internal
 * \brief Density fitting
 *
 * Class that implements the density fitting forces and potential
 * \note the virial calculation is not yet implemented
 */
class DensityFitting final : public IMDModule
{

public:
    //! Construct the density fitting module.
    explicit DensityFitting() = default;

    /*! \brief Request to be notified during pre-processing.
     *
     * \param[in] notifiers allows the module to subscribe to notifications from MDModules.
     *
     * The density fitting code subscribes to these notifications:
     *   - setting atom group indices in the densityFittingOptions_ from an
     *     index group string by taking a parmeter const IndexGroupsAndNames &
     *   - storing its internal parameters in a tpr file by writing to a
     *     key-value-tree during pre-processing by a function taking a
     *     KeyValueTreeObjectBuilder as parameter
     */
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) override
    {
        // Callbacks for several kinds of MDModulesNotifier are created
        // and subscribed, and will be dispatched correctly at run time
        // based on the type of the parameter required by the lambda.

        if (!densityFittingOptions_.active())
        {
            return;
        }

        // Setting the atom group indices from index group string
        const auto setFitGroupIndicesFunction = [this](const IndexGroupsAndNames& indexGroupsAndNames) {
            densityFittingOptions_.setFitGroupIndices(indexGroupsAndNames);
        };
        notifiers->preProcessingNotifier_.subscribe(setFitGroupIndicesFunction);

        // Writing internal parameters during pre-processing
        const auto writeInternalParametersFunction = [this](KeyValueTreeObjectBuilder treeBuilder) {
            densityFittingOptions_.writeInternalParametersToKvt(treeBuilder);
        };
        notifiers->preProcessingNotifier_.subscribe(writeInternalParametersFunction);

        // Checking for consistency with all .mdp options
        const auto checkEnergyCaluclationFrequencyFunction =
                [this](EnergyCalculationFrequencyErrors* energyCalculationFrequencyErrors) {
                    densityFittingOptions_.checkEnergyCaluclationFrequency(energyCalculationFrequencyErrors);
                };
        notifiers->preProcessingNotifier_.subscribe(checkEnergyCaluclationFrequencyFunction);
    }

    /*! \brief Request to be notified.
     * The density fitting code subscribes to these notifications:
     *   - reading its internal parameters from a key-value-tree during
     *     simulation setup by taking a const KeyValueTreeObject & parameter
     *   - constructing local atom sets in the simulation parameter setup
     *     by taking a LocalAtomSetManager * as parameter
     *   - the type of periodic boundary conditions that are used
     *     by taking a PeriodicBoundaryConditionType as parameter
     *   - the writing of checkpoint data
     *     by taking a MDModulesWriteCheckpointData as parameter
     *   - the reading of checkpoint data
     *     by taking a MDModulesCheckpointReadingDataOnMain as parameter
     *   - the broadcasting of checkpoint data
     *     by taking MDModulesCheckpointReadingBroadcast as parameter
     */
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
    {
        if (!densityFittingOptions_.active())
        {
            return;
        }

        // Reading internal parameters during simulation setup
        const auto readInternalParametersFunction = [this](const KeyValueTreeObject& tree) {
            densityFittingOptions_.readInternalParametersFromKvt(tree);
        };
        notifiers->simulationSetupNotifier_.subscribe(readInternalParametersFunction);

        // constructing local atom sets during simulation setup
        const auto setLocalAtomSetFunction = [this](LocalAtomSetManager* localAtomSetManager) {
            this->constructLocalAtomSet(localAtomSetManager);
        };
        notifiers->simulationSetupNotifier_.subscribe(setLocalAtomSetFunction);

        // constructing local atom sets during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc) {
            this->densityFittingSimulationParameters_.setPeriodicBoundaryConditionType(pbc);
        };
        notifiers->simulationSetupNotifier_.subscribe(setPeriodicBoundaryContionsFunction);

        // setting the simulation time step
        const auto setSimulationTimeStepFunction = [this](const SimulationTimeStep& simulationTimeStep) {
            this->densityFittingSimulationParameters_.setSimulationTimeStep(simulationTimeStep.delta_t);
        };
        notifiers->simulationSetupNotifier_.subscribe(setSimulationTimeStepFunction);

        // adding output to energy file
        const auto requestEnergyOutput =
                [this](MDModulesEnergyOutputToDensityFittingRequestChecker* energyOutputRequest) {
                    this->setEnergyOutputRequest(energyOutputRequest);
                };
        notifiers->simulationSetupNotifier_.subscribe(requestEnergyOutput);

        // writing checkpoint data
        const auto checkpointDataWriting = [this](MDModulesWriteCheckpointData checkpointData) {
            forceProvider_->writeCheckpointData(checkpointData, DensityFittingModuleInfo::name_);
        };
        notifiers->checkpointingNotifier_.subscribe(checkpointDataWriting);

        // reading checkpoint data
        const auto checkpointDataReading = [this](MDModulesCheckpointReadingDataOnMain checkpointData) {
            densityFittingState_.readState(checkpointData.checkpointedData_,
                                           DensityFittingModuleInfo::name_);
        };
        notifiers->checkpointingNotifier_.subscribe(checkpointDataReading);

        // broadcasting checkpoint data
        const auto checkpointDataBroadcast = [this](MDModulesCheckpointReadingBroadcast checkpointData) {
            densityFittingState_.broadcastState(checkpointData.communicator_, checkpointData.isParallelRun_);
        };
        notifiers->checkpointingNotifier_.subscribe(checkpointDataBroadcast);
    }

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &densityFittingOptions_; }

    //! Add this module to the force providers if active
    void initForceProviders(ForceProviders* forceProviders) override
    {
        if (densityFittingOptions_.active())
        {
            const auto& parameters = densityFittingOptions_.buildParameters();
            densityFittingSimulationParameters_.readReferenceDensityFromFile(
                    densityFittingOptions_.referenceDensityFileName());
            if (parameters.normalizeDensities_)
            {
                densityFittingSimulationParameters_.normalizeReferenceDensity();
            }
            forceProvider_ = std::make_unique<DensityFittingForceProvider>(
                    parameters,
                    densityFittingSimulationParameters_.referenceDensity(),
                    densityFittingSimulationParameters_.transformationToDensityLattice(),
                    densityFittingSimulationParameters_.localAtomSet(),
                    densityFittingSimulationParameters_.periodicBoundaryConditionType(),
                    densityFittingSimulationParameters_.simulationTimeStep(),
                    densityFittingState_);
            forceProviders->addForceProvider(forceProvider_.get());
        }
    }

    //! This MDModule provides its own output
    IMDOutputProvider* outputProvider() override { return &densityFittingOutputProvider_; }

    /*! \brief Set up the local atom sets that are used by this module.
     *
     * \note When density fitting is set up with MDModulesNotifier in
     *       the constructor, this function is called back.
     *
     * \param[in] localAtomSetManager the manager to add local atom sets.
     */
    void constructLocalAtomSet(LocalAtomSetManager* localAtomSetManager)
    {
        LocalAtomSet atomSet = localAtomSetManager->add(densityFittingOptions_.buildParameters().indices_);
        densityFittingSimulationParameters_.setLocalAtomSet(atomSet);
    }

    /*! \brief Request energy output to energy file during simulation.
     */
    void setEnergyOutputRequest(MDModulesEnergyOutputToDensityFittingRequestChecker* energyOutputRequest)
    {
        energyOutputRequest->energyOutputToDensityFitting_ = densityFittingOptions_.active();
    }

private:
    //! The output provider
    DensityFittingOutputProvider densityFittingOutputProvider_;
    //! The options provided for density fitting
    DensityFittingOptions densityFittingOptions_;
    //! Object that evaluates the forces
    std::unique_ptr<DensityFittingForceProvider> forceProvider_;
    /*! \brief Parameters for density fitting that become available at
     * simulation setup time.
     */
    DensityFittingSimulationParameterSetup densityFittingSimulationParameters_;
    //! The internal parameters of density fitting force provider
    DensityFittingForceProviderState densityFittingState_;

    GMX_DISALLOW_COPY_AND_ASSIGN(DensityFitting);
};

} // namespace

std::unique_ptr<IMDModule> DensityFittingModuleInfo::create()
{
    return std::make_unique<DensityFitting>();
}

const std::string DensityFittingModuleInfo::name_ = "density-guided-simulation";

} // namespace gmx
