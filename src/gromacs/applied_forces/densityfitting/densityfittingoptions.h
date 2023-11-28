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
 * Declares options for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_DENSITYFITTINGOPTIONS_H
#define GMX_APPLIED_FORCES_DENSITYFITTINGOPTIONS_H

#include <string>

#include "gromacs/mdtypes/imdpoptionprovider.h"

#include "densityfittingparameters.h"

namespace gmx
{

class EnergyCalculationFrequencyErrors;
class IndexGroupsAndNames;
class KeyValueTreeObject;
class KeyValueTreeBuilder;

/*! \internal
 * \brief Input data storage for density fitting
 */
class DensityFittingOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for density fitting to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     * \note This should be symmetrical to option initialization without
     *       employing manual prefixing with the section name string once
     *       the legacy code blocking this design is removed.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Report if this set of options is active
    bool active() const;

    //! Process input options to parameters, including input file reading.
    const DensityFittingParameters& buildParameters();

    /*! \brief Evaluate and store atom indices.
     *
     * During pre-processing, use the group string from the options to
     * evaluate the indices of the atoms to be subject to forces from this
     * module.
     */
    void setFitGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames);

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    //! Return the file name of the reference density
    const std::string& referenceDensityFileName() const;

    //! Check if input parameters are consistent with other simulation parameters
    void checkEnergyCaluclationFrequency(EnergyCalculationFrequencyErrors* energyCalculationFrequencyErrors) const;

private:
    const std::string c_activeTag_ = "active";

    /*! \brief Denote the .mdp option that defines the group of fit atoms.
     * \note Changing this string will break .tpr backwards compatibility
     */
    const std::string c_groupTag_  = "group";
    std::string       groupString_ = "protein";

    const std::string c_similarityMeasureTag_ = "similarity-measure";

    const std::string c_amplitudeMethodTag_ = "atom-spreading-weight";

    const std::string c_forceConstantTag_ = "force-constant";

    const std::string c_gaussianTransformSpreadingWidthTag_ = "gaussian-transform-spreading-width";
    const std::string c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_ =
            "gaussian-transform-spreading-range-in-multiples-of-width";

    const std::string c_referenceDensityFileNameTag_ = "reference-density-filename";
    std::string       referenceDensityFileName_      = "reference.mrc";

    const std::string c_everyNStepsTag_ = "nst";

    const std::string c_normalizeDensitiesTag_ = "normalize-densities";

    const std::string c_adaptiveForceScalingTag_ = "adaptive-force-scaling";

    const std::string c_adaptiveForceScalingTimeConstantTag_ =
            "adaptive-force-scaling-time-constant";

    const std::string c_translationTag_ = "shift-vector";

    const std::string c_transformationMatrixTag_ = "transformation-matrix";

    DensityFittingParameters parameters_;
};

} // namespace gmx

#endif
