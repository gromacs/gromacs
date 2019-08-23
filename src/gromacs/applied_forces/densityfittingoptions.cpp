/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
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
template <class ToType, class TransformWithFunctionType>
void densityfittingMdpTransformFromString(IKeyValueTreeTransformRules * rules,
                                          TransformWithFunctionType     transformationFunction,
                                          const std::string            &optionTag)
{
    rules->addRule()
        .from<std::string>("/" + DensityFittingModuleInfo::name_ + "-" + optionTag)
        .to<ToType>("/" + DensityFittingModuleInfo::name_ +"/" + optionTag)
        .transformWith(transformationFunction);
}
/*! \brief Helper to declare mdp output.
 *
 * Enforces uniform mdp options output sting that are always prepended with the
 * correct string for the densityfitting mdp options and is consistent with the
 * options name and transformation.
 *
 * \tparam OptionType the type of the mdp option
 * \param[in] builder the KVT builder to generate the output
 * \param[in] option the mdp option
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the density guided simulation
 */
template <class OptionType>
void addDensityFittingMdpOutputValue(KeyValueTreeObjectBuilder *builder,
                                     const OptionType          &option,
                                     const std::string         &optionTag)
{
    builder->addValue<OptionType>(DensityFittingModuleInfo::name_ + "-" + optionTag,
                                  option);
}

}   // namespace

void DensityFittingOptions::initMdpTransform(IKeyValueTreeTransformRules * rules)
{
    const auto &stringIdentityTransform = [](std::string s){
            return s;
        };
    densityfittingMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_groupTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_similarityMeasureTag_);
    densityfittingMdpTransformFromString<std::string>(rules, stringIdentityTransform, c_amplitudeMethodTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>, c_forceConstantTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>, c_gaussianTransformSpreadingWidthTag_);
    densityfittingMdpTransformFromString<real>(rules, &fromStdString<real>, c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_);
}

void DensityFittingOptions::buildMdpOutput(KeyValueTreeObjectBuilder *builder) const
{
    addDensityFittingMdpOutputValue(builder, parameters_.active_, c_activeTag_);
    if (parameters_.active_)
    {
        addDensityFittingMdpOutputValue(builder, groupString_, c_groupTag_);
        addDensityFittingMdpOutputValue<std::string>(builder,
                                                     c_densitySimilarityMeasureMethodNames[parameters_.similarityMeasureMethod_],
                                                     c_similarityMeasureTag_);

        addDensityFittingMdpOutputValue<std::string>(builder,
                                                     c_densityFittingAmplitudeMethodNames[parameters_.amplitudeLookupMethod_],
                                                     c_amplitudeMethodTag_);

        addDensityFittingMdpOutputValue(builder, parameters_.forceConstant_, c_forceConstantTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.gaussianTransformSpreadingWidth_, c_gaussianTransformSpreadingWidthTag_);
        addDensityFittingMdpOutputValue(builder, parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_, c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_);
    }
}

void DensityFittingOptions::initMdpOptions(IOptionsContainerWithSections *options)
{
    auto section = options->addSection(OptionSection(DensityFittingModuleInfo::name_.c_str()));

    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&parameters_.active_));
    section.addOption(StringOption(c_groupTag_.c_str()).store(&groupString_));

    section.addOption(EnumOption<DensitySimilarityMeasureMethod>(c_similarityMeasureTag_.c_str())
                          .enumValue(c_densitySimilarityMeasureMethodNames.m_elements)
                          .store(&parameters_.similarityMeasureMethod_));

    section.addOption(EnumOption<DensityFittingAmplitudeMethod>(c_amplitudeMethodTag_.c_str())
                          .enumValue(c_densityFittingAmplitudeMethodNames.m_elements)
                          .store(&parameters_.amplitudeLookupMethod_));

    section.addOption(RealOption(c_forceConstantTag_.c_str()).store(&parameters_.forceConstant_));
    section.addOption(RealOption(c_gaussianTransformSpreadingWidthTag_.c_str()).store(&parameters_.gaussianTransformSpreadingWidth_));
    section.addOption(RealOption(c_gaussianTransformSpreadingRangeInMultiplesOfWidthTag_.c_str()).store(&parameters_.gaussianTransformSpreadingRangeInMultiplesOfWidth_));
}

bool DensityFittingOptions::active() const
{
    return parameters_.active_;
}

const DensityFittingParameters &DensityFittingOptions::buildParameters()
{
    return parameters_;
}

void DensityFittingOptions::setFitGroupIndices(const IndexGroupsAndNames &indexGroupsAndNames)
{
    if (!parameters_.active_)
    {
        return;
    }
    parameters_.indices_ = indexGroupsAndNames.indices(groupString_);
}

void DensityFittingOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    auto groupIndexAdder = treeBuilder.addUniformArray<index>(DensityFittingModuleInfo::name_ + "-" + c_groupTag_);
    for (const auto &indexValue : parameters_.indices_)
    {
        groupIndexAdder.addValue(indexValue);
    }
}

void DensityFittingOptions::readInternalParametersFromKvt(const KeyValueTreeObject &tree)
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
                   [](const KeyValueTreeValue &val) { return val.cast<index>(); });
}

} // namespace gmx
