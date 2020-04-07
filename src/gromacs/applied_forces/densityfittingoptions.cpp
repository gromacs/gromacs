/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
 * \brief
 * Implements force provider for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "densityfittingoptions.h"

#include "gromacs/applied_forces/densityfitting.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/mdmodulenotification.h"
#include "gromacs/utility/strconvert.h"

#include "densityfittingamplitudelookup.h"

namespace gmx
{

namespace
{

/*! \brief Helper to declare mdp transform rules.
 *
 * Enforces uniform mdp options that are always prepended with the correct
 * string for the densityfitting mdp options.
 *
 * \tparam ToType type to be transformed to
 * \tparam TransformWithFunctionType type of transformation function to be used
 *
 * \param[in] rules KVT transformation rules
 * \param[in] transformationFunction the function to transform the flat kvt tree
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the density guided simulation
 */
template<class ToType, class TransformWithFunctionType>
void densityfittingMdpTransformFromString(IKeyValueTreeTransformRules* rules,
                                          TransformWithFunctionType    transformationFunction,
                                          const std::string&           optionTag)
{
    rules->addRule()
            .from<std::string>("/" + DensityFittingModuleInfo::name_ + "-" + optionTag)
            .to<ToType>("/" + DensityFittingModuleInfo::name_ + "/" + optionTag)
            .transformWith(transformationFunction);
}
/*! \brief Helper to declare mdp output.
 *
 * Enforces uniform mdp options output strings that are always prepended with the
 * correct string for the densityfitting mdp options and are consistent with the
 * options name and transformation type.
 *
 * \tparam OptionType the type of the mdp option
 * \param[in] builder the KVT builder to generate the output
 * \param[in] option the mdp option
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the density guided simulation
 */
template<class OptionType>
void addDensityFittingMdpOutputValue(KeyValueTreeObjectBuilder* builder,
                                     const OptionType&          option,
                                     const std::string&         optionTag)
{
    builder->addValue<OptionType>(DensityFittingModuleInfo::name_ + "-" + optionTag, option);
}

/*! \brief Helper to declare mdp output comments.
 *
 * Enforces uniform mdp options comment output strings that are always prepended
 * with the correct string for the densityfitting mdp options and are consistent
 * with the options name and transformation type.
 *
 * \param[in] builder the KVT builder to generate the output
 * \param[in] comment on the mdp option
 * \param[in] optionTag string tag that describes the mdp option
 */
void addDensityFittingMdpOutputValueComment(KeyValueTreeObjectBuilder* builder,
                                            const std::string&         comment,
                                            const std::string&         optionTag)
{
    builder->addValue<std::string>("comment-" + DensityFittingModuleInfo::name_ + "-" + optionTag, comment);
}

} // namespace

void DensityFittingOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    densityfittingMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_groupTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform,
                                                      c_similarityMeasureTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_amplitudeMethodTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>, c_forceConstantTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>,
                                               c_gaussianTransformSpreadingWidthTag_);
    densityfittingMdpTransformFromString<real>(
            rules, &fromStdString<real>, c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform,
                                                      c_referenceDensityFileNameTag_);
    densityfittingMdpTransformFromString<std::int64_t>(rules, &fromStdString<std::int64_t>,
                                                       c_everyNStepsTag_);
    densityfittingMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_normalizeDensitiesTag_);
    densityfittingMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_adaptiveForceScalingTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>,
                                               c_adaptiveForceScalingTimeConstantTag_);
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

    addDensityFittingMdpOutputValueComment(builder, "", "empty-line");
    addDensityFittingMdpOutputValueComment(builder, "; Density guided simulation", "module");

    addDensityFittingMdpOutputValue(builder, parameters_.active_, c_activeTag_);
    if (parameters_.active_)
    {
        addDensityFittingMdpOutputValue(builder, groupString_, c_groupTag_);

        addDensityFittingMdpOutputValueComment(
                builder,
                "; Similarity measure between densities: inner-product, relative-entropy, or "
                "cross-correlation",
                c_similarityMeasureTag_);
        addDensityFittingMdpOutputValue<std::string>(
                builder, c_densitySimilarityMeasureMethodNames[parameters_.similarityMeasureMethod_],
                c_similarityMeasureTag_);

        addDensityFittingMdpOutputValueComment(
                builder, "; Atom amplitude for spreading onto grid: unity, mass, or charge",
                c_amplitudeMethodTag_);
        addDensityFittingMdpOutputValue<std::string>(
                builder, c_densityFittingAmplitudeMethodNames[parameters_.amplitudeLookupMethod_],
                c_amplitudeMethodTag_);

        addDensityFittingMdpOutputValue(builder, parameters_.forceConstant_, c_forceConstantTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.gaussianTransformSpreadingWidth_,
                                        c_gaussianTransformSpreadingWidthTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_,
                                        c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_);
        addDensityFittingMdpOutputValueComment(builder,
                                               "; Reference density file location as absolute path "
                                               "or relative to the gmx mdrun calling location",
                                               c_referenceDensityFileNameTag_);
        addDensityFittingMdpOutputValue(builder, referenceDensityFileName_, c_referenceDensityFileNameTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.calculationIntervalInSteps_, c_everyNStepsTag_);
        addDensityFittingMdpOutputValueComment(
                builder, "; Normalize the sum of density voxel values to one", c_normalizeDensitiesTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.normalizeDensities_, c_normalizeDensitiesTag_);
        addDensityFittingMdpOutputValueComment(builder, "; Apply adaptive force scaling",
                                               c_adaptiveForceScalingTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.adaptiveForceScaling_,
                                        c_adaptiveForceScalingTag_);
        addDensityFittingMdpOutputValueComment(builder,
                                               "; Time constant for adaptive force scaling in ps",
                                               c_adaptiveForceScalingTimeConstantTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.adaptiveForceScalingTimeConstant_,
                                        c_adaptiveForceScalingTimeConstantTag_);
    }
}

void DensityFittingOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(DensityFittingModuleInfo::name_.c_str()));

    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&parameters_.active_));
    section.addOption(StringOption(c_groupTag_.c_str()).store(&groupString_));

    section.addOption(EnumOption<DensitySimilarityMeasureMethod>(c_similarityMeasureTag_.c_str())
                              .enumValue(c_densitySimilarityMeasureMethodNames)
                              .store(&parameters_.similarityMeasureMethod_));

    section.addOption(EnumOption<DensityFittingAmplitudeMethod>(c_amplitudeMethodTag_.c_str())
                              .enumValue(c_densityFittingAmplitudeMethodNames)
                              .store(&parameters_.amplitudeLookupMethod_));

    section.addOption(RealOption(c_forceConstantTag_.c_str()).store(&parameters_.forceConstant_));
    section.addOption(RealOption(c_gaussianTransformSpreadingWidthTag_.c_str())
                              .store(&parameters_.gaussianTransformSpreadingWidth_));
    section.addOption(RealOption(c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_.c_str())
                              .store(&parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_));
    section.addOption(StringOption(c_referenceDensityFileNameTag_.c_str()).store(&referenceDensityFileName_));
    section.addOption(Int64Option(c_everyNStepsTag_.c_str()).store(&parameters_.calculationIntervalInSteps_));
    section.addOption(BooleanOption(c_normalizeDensitiesTag_.c_str()).store(&parameters_.normalizeDensities_));
    section.addOption(
            BooleanOption(c_adaptiveForceScalingTag_.c_str()).store(&parameters_.adaptiveForceScaling_));
    section.addOption(RealOption(c_adaptiveForceScalingTimeConstantTag_.c_str())
                              .store(&parameters_.adaptiveForceScalingTimeConstant_));
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
    auto groupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(DensityFittingModuleInfo::name_
                                                                     + "-" + c_groupTag_);
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

    if (!tree.keyExists(DensityFittingModuleInfo::name_ + "-" + c_groupTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find atom index vector required for density guided simulation."));
    }
    auto kvtIndexArray = tree[DensityFittingModuleInfo::name_ + "-" + c_groupTag_].asArray().values();
    parameters_.indices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray), std::end(kvtIndexArray), std::begin(parameters_.indices_),
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
                + ") is not a multiple of " + DensityFittingModuleInfo::name_ + "-" + c_everyNStepsTag_
                + " (" + toString(parameters_.calculationIntervalInSteps_) + ") .");
    }
}

const std::string& DensityFittingOptions::referenceDensityFileName() const
{
    return referenceDensityFileName_;
}
} // namespace gmx
