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
 * Declares force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H

#include <memory>

#include "gromacs/fileio/checkpoint.h"
#include "gromacs/math/exponentialmovingaverage.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/iforceprovider.h"

enum class PbcType : int;

namespace gmx
{

class LocalAtomSet;
class TranslateAndScale;
struct DensityFittingParameters;

/*! \internal
 * \brief Parameters defining the internal density fitting force provider state.
 */
struct DensityFittingForceProviderState
{
    /*! \brief The steps since the last force calculation.
     *  Used if density fitting is to be calculated every N steps.
     */
    std::int64_t stepsSinceLastCalculation_ = 0;
    /*! \brief String naming variable holding the steps since last calculation.
     * \note Changing this name will break backwards compatibility for checkpoint file writing.
     */
    static const std::string stepsSinceLastCalculationName_;
    //! The state of the exponential moving average of the similarity measure
    ExponentialMovingAverageState exponentialMovingAverageState_ = {};
    /*! \brief String naming variable holding the exponential moving average.
     * \note Changing this name will break backwards compatibility for checkpoint file writing.
     */
    static const std::string exponentialMovingAverageStateName_;

    //! An additional factor scaling the force for adaptive force scaling
    real adaptiveForceConstantScale_ = 1.0_real;
    /*! \brief String naming variable holding the adaptive force constant scale.
     * \note Changing this name will break backwards compatibility for checkpoint file writing.
     */
    static const std::string adaptiveForceConstantScaleName_;

    /*! \brief Write internal density fitting data into a key value tree.
     * The entries to the kvt are identified with identifier, so that a variable
     * is indentified with the key "identifier-variablename"
     *
     * \param[in] kvtBuilder enables writing to the Key-Value-Tree
     *                              the state is written to
     *
     * \param[in] identifier denotes the module that is checkpointing the data
     */
    void writeState(KeyValueTreeObjectBuilder kvtBuilder, const std::string& identifier) const;

    /*! \brief Read the internal parameters from the checkpoint file on main
     * \param[in] kvtData holding the checkpoint information
     * \param[in] identifier identifies the data in a key-value-tree
     */
    void readState(const KeyValueTreeObject& kvtData, const std::string& identifier);

    /*! \brief Broadcast the internal parameters.
     *
     * \param[in] communicator to broadcast the state information
     * \param[in] isParallelRun to determine if anything has to be broadcast at all
     *
     */
    void broadcastState(MPI_Comm communicator, bool isParallelRun);
};

/*! \internal \brief
 * Implements IForceProvider for density-fitting forces.
 */
class DensityFittingForceProvider final : public IForceProvider
{
public:
    //! Construct force provider for density fitting from its parameters
    DensityFittingForceProvider(const DensityFittingParameters&             parameters,
                                basic_mdspan<const float, dynamicExtents3D> referenceDensity,
                                const TranslateAndScale& transformationToDensityLattice,
                                const LocalAtomSet&      localAtomSet,
                                PbcType                  pbcType,
                                double                   simulationTimeStep,
                                const DensityFittingForceProviderState& state);
    ~DensityFittingForceProvider();
    /*!\brief Calculate forces that maximise goodness-of-fit with a reference density map.
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    /*! \brief Write internal density fitting data to checkpoint file.
     * \param[in] checkpointWriting enables writing to the Key-Value-Tree
     *                              that is used for storing the checkpoint
     *                              information
     * \param[in] moduleName names the module that is checkpointing this force-provider
     *
     * \note The provided state to checkpoint has to change if checkpointing
     *       is moved before the force provider call in the MD-loop.
     */
    void writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting, const std::string& moduleName);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H
