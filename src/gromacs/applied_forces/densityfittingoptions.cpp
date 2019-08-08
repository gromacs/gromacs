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
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

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
    densityfittingMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
}

void DensityFittingOptions::buildMdpOutput(KeyValueTreeObjectBuilder *builder) const
{
    addDensityFittingMdpOutputValue(builder, parameters_.active_, c_activeTag_);
}

void DensityFittingOptions::initMdpOptions(IOptionsContainerWithSections *options)
{
    auto section = options->addSection(OptionSection(DensityFittingModuleInfo::name_.c_str()));

    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&parameters_.active_));
}

bool DensityFittingOptions::active() const
{
    return parameters_.active_;
}

const DensityFittingParameters &DensityFittingOptions::buildParameters()
{
    return parameters_;
}

} // namespace gmx
