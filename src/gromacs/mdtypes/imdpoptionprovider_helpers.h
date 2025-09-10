/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 *  \brief Defines helper functions to support mdp option transformation and output.
 *
 *  \ingroup module_options
 */

#ifndef GMX_IMDPOPTIONPROVIDER_HELPERS_H
#define GMX_IMDPOPTIONPROVIDER_HELPERS_H

#include <string>
#include <string_view>

#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{
/*! \brief Helper to declare mdp transform rules.
 *
 * Enforces uniform mdp option transformation for a given module.
 *
 * \tparam ToType type to be transformed to
 * \tparam TransformWithFunctionType type of transformation function from std::string to ToType
 *
 * \param[in] rules KVT transformation rules
 * \param[in] transformationFunction function to transform string to ToType
 * \param[in] moduleName module name prepended to the option name
 * \param[in] optionName suffix of the mdp option name
 */
template<class ToType, class TransformWithFunctionType>
void addMdpTransformFromString(IKeyValueTreeTransformRules* rules,
                               TransformWithFunctionType    transformationFunction,
                               std::string_view             moduleName,
                               std::string_view             optionName)
{
    std::string prefix("/");
    prefix.append(moduleName);
    std::string fromKey(prefix);
    fromKey.append("-").append(optionName);
    std::string toKey(prefix);
    toKey.append("/").append(optionName);

    rules->addRule().from<std::string>(fromKey).to<ToType>(toKey).transformWith(transformationFunction);
}

/*! \brief Add an mdp value to output KVT with a module-prefixed name
 *
 * \tparam ValueType            Type of the value described by the option
 *
 * \param[in] builder            Builder for KVT output
 * \param[in] moduleName         Module name prefix
 * \param[in] optionName         Suffix of the mdp option name
 * \param[in] optionValue        Value to write
 */
template<class ValueType>
void addMdpOutputValue(KeyValueTreeObjectBuilder* builder,
                       std::string_view           moduleName,
                       std::string_view           optionName,
                       const ValueType&           optionValue)
{
    std::string key(moduleName);
    key.append("-").append(optionName);
    builder->addValue<ValueType>(key, optionValue);
}

/*! \brief Add a comment for an mdp value with a module-prefixed name
 *
 * \param[in] builder            Builder for KVT output
 * \param[in] moduleName         Module name prefix
 * \param[in] optionName         Suffix of the mdp option name
 * \param[in] comment            Comment text to write
 */
inline void addMdpOutputComment(KeyValueTreeObjectBuilder* builder,
                                std::string_view           moduleName,
                                std::string_view           optionName,
                                const std::string&         comment)
{
    std::string key("comment-");
    key.append(moduleName).append("-").append(optionName);
    builder->addValue<std::string>(key, comment);
}

} // namespace gmx

#endif // GMX_IMDPOPTIONPROVIDER_HELPERS_H
