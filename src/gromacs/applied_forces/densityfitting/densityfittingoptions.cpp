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
 * Implements force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingoptions.h"

#include "gromacs/applied_forces/densityfitting/densityfitting.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"

#include "densityfittingamplitudelookup.h"

namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(DensityFittingModuleInfo::sc_name);
}


const std::string c_activeTag = "active";

/*! \brief Denote the .mdp option that defines the group of fit atoms.
 * \note Changing this string will break .tpr backwards compatibility
 */
const std::string c_groupTag = "group";

const std::string c_similarityMeasureTag = "similarity-measure";

const std::string c_amplitudeMethodTag = "atom-spreading-weight";

const std::string c_forceConstantTag = "force-constant";

const std::string c_gaussianTransformSpreadingWidthTag = "gaussian-transform-spreading-width";
const std::string c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag =
        "gaussian-transform-spreading-range-in-multiples-of-width";

const std::string c_referenceDensityFileNameTag = "reference-density-filename";

const std::string c_everyNStepsTag = "nst";

const std::string c_normalizeDensitiesTag = "normalize-densities";

const std::string c_adaptiveForceScalingTag = "adaptive-force-scaling";

const std::string c_adaptiveForceScalingTimeConstantTag = "adaptive-force-scaling-time-constant";

const std::string c_translationTag = "shift-vector";

const std::string c_transformationMatrixTag = "transformation-matrix";

} // namespace

void DensityFittingOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(
            rules, &fromStdString<bool>, DensityFittingModuleInfo::sc_name, c_activeTag);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, DensityFittingModuleInfo::sc_name, c_groupTag);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, DensityFittingModuleInfo::sc_name, c_similarityMeasureTag);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, DensityFittingModuleInfo::sc_name, c_amplitudeMethodTag);
    addMdpTransformFromString<real>(
            rules, &fromStdString<real>, DensityFittingModuleInfo::sc_name, c_forceConstantTag);
    addMdpTransformFromString<real>(
            rules, &fromStdString<real>, DensityFittingModuleInfo::sc_name, c_gaussianTransformSpreadingWidthTag);
    addMdpTransformFromString<real>(rules,
                                    &fromStdString<real>,
                                    DensityFittingModuleInfo::sc_name,
                                    c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, DensityFittingModuleInfo::sc_name, c_referenceDensityFileNameTag);
    addMdpTransformFromString<std::int64_t>(
            rules, &fromStdString<std::int64_t>, DensityFittingModuleInfo::sc_name, c_everyNStepsTag);
    addMdpTransformFromString<bool>(
            rules, &fromStdString<bool>, DensityFittingModuleInfo::sc_name, c_normalizeDensitiesTag);
    addMdpTransformFromString<bool>(
            rules, &fromStdString<bool>, DensityFittingModuleInfo::sc_name, c_adaptiveForceScalingTag);
    addMdpTransformFromString<real>(rules,
                                    &fromStdString<real>,
                                    DensityFittingModuleInfo::sc_name,
                                    c_adaptiveForceScalingTimeConstantTag);

    const auto& stringRVecToStringRVecWithCheck = [](const std::string& str)
    {
        return stringIdentityTransformWithArrayCheck<real, 3>(
                str,
                "Reading three real values as vector while parsing the .mdp input failed in "
                        + moduleName() + ".");
    };
    addMdpTransformFromString<std::string>(
            rules, stringRVecToStringRVecWithCheck, DensityFittingModuleInfo::sc_name, c_translationTag);

    const auto& stringMatrixToStringMatrixWithCheck = [](const std::string& str)
    {
        return stringIdentityTransformWithArrayCheck<real, 9>(
                str,
                "Reading nine real values as vector while parsing the .mdp input failed in "
                        + moduleName() + ".");
    };
    addMdpTransformFromString<std::string>(rules,
                                           stringMatrixToStringMatrixWithCheck,
                                           DensityFittingModuleInfo::sc_name,
                                           c_transformationMatrixTag);
}

//! Name the methods that may be used to evaluate similarity between densities
static const EnumerationArray<DensitySimilarityMeasureMethod, const char*> c_densitySimilarityMeasureMethodNames = {
    { "inner-product", "relative-entropy", "cross-correlation" }
};
//! The names of the methods to determine the amplitude of the atoms to be spread on a grid
static const EnumerationArray<DensityFittingAmplitudeMethod, const char*> c_densityFittingAmplitudeMethodNames = {
    { "unity", "mass", "charge" }
};

void DensityFittingOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, DensityFittingModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(
            builder, DensityFittingModuleInfo::sc_name, "module", "; Density guided simulation");
    addMdpOutputValue(builder, DensityFittingModuleInfo::sc_name, c_activeTag, parameters_.active_);

    if (parameters_.active_)
    {
        addMdpOutputValue(builder, DensityFittingModuleInfo::sc_name, c_groupTag, groupString_);

        addMdpOutputComment(
                builder,
                DensityFittingModuleInfo::sc_name,
                c_similarityMeasureTag,
                "; Similarity measure between densities: inner-product, relative-entropy, or "
                "cross-correlation");
        addMdpOutputValue<std::string>(
                builder,
                DensityFittingModuleInfo::sc_name,
                c_similarityMeasureTag,
                c_densitySimilarityMeasureMethodNames[parameters_.similarityMeasureMethod_]);

        addMdpOutputComment(builder,
                            DensityFittingModuleInfo::sc_name,
                            c_amplitudeMethodTag,
                            "; Atom amplitude for spreading onto grid: unity, mass, or charge");
        addMdpOutputValue<std::string>(
                builder,
                DensityFittingModuleInfo::sc_name,
                c_amplitudeMethodTag,
                c_densityFittingAmplitudeMethodNames[parameters_.amplitudeLookupMethod_]);
        addMdpOutputValue(
                builder, DensityFittingModuleInfo::sc_name, c_forceConstantTag, parameters_.forceConstant_);
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_gaussianTransformSpreadingWidthTag,
                          parameters_.gaussianTransformSpreadingWidth_);
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag,
                          parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_);

        addMdpOutputComment(builder,
                            DensityFittingModuleInfo::sc_name,
                            c_referenceDensityFileNameTag,
                            "; Reference density file location as absolute path "
                            "or relative to the gmx mdrun calling location");
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_referenceDensityFileNameTag,
                          referenceDensityFileName_);
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_everyNStepsTag,
                          parameters_.calculationIntervalInSteps_);

        addMdpOutputComment(builder,
                            DensityFittingModuleInfo::sc_name,
                            c_normalizeDensitiesTag,
                            "; Normalize the sum of density voxel values to one");
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_normalizeDensitiesTag,
                          parameters_.normalizeDensities_);

        addMdpOutputComment(builder,
                            DensityFittingModuleInfo::sc_name,
                            c_adaptiveForceScalingTag,
                            "; Apply adaptive force scaling");
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_adaptiveForceScalingTag,
                          parameters_.adaptiveForceScaling_);

        addMdpOutputComment(builder,
                            DensityFittingModuleInfo::sc_name,
                            c_adaptiveForceScalingTimeConstantTag,
                            "; Time constant for adaptive force scaling in ps");
        addMdpOutputValue(builder,
                          DensityFittingModuleInfo::sc_name,
                          c_adaptiveForceScalingTimeConstantTag,
                          parameters_.adaptiveForceScalingTimeConstant_);
    }
}

void DensityFittingOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));

    section.addOption(BooleanOption(c_activeTag.c_str()).store(&parameters_.active_));
    section.addOption(StringOption(c_groupTag.c_str()).store(&groupString_));

    section.addOption(EnumOption<DensitySimilarityMeasureMethod>(c_similarityMeasureTag.c_str())
                              .enumValue(c_densitySimilarityMeasureMethodNames)
                              .store(&parameters_.similarityMeasureMethod_));

    section.addOption(EnumOption<DensityFittingAmplitudeMethod>(c_amplitudeMethodTag.c_str())
                              .enumValue(c_densityFittingAmplitudeMethodNames)
                              .store(&parameters_.amplitudeLookupMethod_));

    section.addOption(RealOption(c_forceConstantTag.c_str()).store(&parameters_.forceConstant_));
    section.addOption(RealOption(c_gaussianTransformSpreadingWidthTag.c_str())
                              .store(&parameters_.gaussianTransformSpreadingWidth_));
    section.addOption(RealOption(c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag.c_str())
                              .store(&parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_));
    section.addOption(StringOption(c_referenceDensityFileNameTag.c_str()).store(&referenceDensityFileName_));
    section.addOption(Int64Option(c_everyNStepsTag.c_str()).store(&parameters_.calculationIntervalInSteps_));
    section.addOption(BooleanOption(c_normalizeDensitiesTag.c_str()).store(&parameters_.normalizeDensities_));
    section.addOption(
            BooleanOption(c_adaptiveForceScalingTag.c_str()).store(&parameters_.adaptiveForceScaling_));
    section.addOption(RealOption(c_adaptiveForceScalingTimeConstantTag.c_str())
                              .store(&parameters_.adaptiveForceScalingTimeConstant_));
    section.addOption(StringOption(c_translationTag.c_str()).store(&parameters_.translationString_));
    section.addOption(
            StringOption(c_transformationMatrixTag.c_str()).store(&parameters_.transformationMatrixString_));
}

bool DensityFittingOptions::active() const
{
    return parameters_.active_;
}

const DensityFittingParameters& DensityFittingOptions::buildParameters()
{
    // the options modules does not know unsigned integers so any input of this
    // kind is rectified here
    if (parameters_.calculationIntervalInSteps_ < 1)
    {
        parameters_.calculationIntervalInSteps_ = 1;
    }
    return parameters_;
}

void DensityFittingOptions::setFitGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames)
{
    if (!parameters_.active_)
    {
        return;
    }
    parameters_.indices_ = indexGroupsAndNames.indices(groupString_);
}

void DensityFittingOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    auto groupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_groupTag);
    for (const auto& indexValue : parameters_.indices_)
    {
        groupIndexAdder.addValue(indexValue);
    }
}

void DensityFittingOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{
    if (!parameters_.active_)
    {
        return;
    }

    if (!tree.keyExists(moduleName() + "-" + c_groupTag))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find atom index vector required for density guided simulation."));
    }
    auto kvtIndexArray = tree[moduleName() + "-" + c_groupTag].asArray().values();
    parameters_.indices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(parameters_.indices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });
}

void DensityFittingOptions::checkEnergyCaluclationFrequency(
        EnergyCalculationFrequencyErrors* energyCalculationFrequencyErrors) const
{
    if (energyCalculationFrequencyErrors->energyCalculationIntervalInSteps() % parameters_.calculationIntervalInSteps_
        != 0)
    {
        energyCalculationFrequencyErrors->addError(
                "nstcalcenergy ("
                + toString(energyCalculationFrequencyErrors->energyCalculationIntervalInSteps())
                + ") is not a multiple of " + moduleName() + "-" + c_everyNStepsTag + " ("
                + toString(parameters_.calculationIntervalInSteps_) + ") .");
    }
}

const std::string& DensityFittingOptions::referenceDensityFileName() const
{
    return referenceDensityFileName_;
}
} // namespace gmx
