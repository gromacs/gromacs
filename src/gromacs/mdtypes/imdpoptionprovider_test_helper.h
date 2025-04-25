/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Defines a helper template for use in tests of IMdpOptionsProviders
 *
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IMDPOPTIONPROVIDER_TEST_HELPER_H
#define GMX_MDTYPES_IMDPOPTIONPROVIDER_TEST_HELPER_H

#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringcompare.h"

namespace gmx
{

class IKeyValueTreeTransformRules;
class IOptionsContainerWithSections;
class KeyValueTreeObjectBuilder;

namespace test
{

/*! \brief Convenience method to convert mdp options from a key-value tree
 * to a module's IMdpOptionsProvider.
 *
 * Converts the set of .mdp option values for a module like grompp
 * does with the text-file input and returns them in the module's
 * IMdpOptionsProvider.
 */
template<typename MdpOptionsProvider>
MdpOptionsProvider fillOptionsFromMdpValuesTemplate(const KeyValueTreeObject& moduleMdpValues)
{
    MdpOptionsProvider moduleMdpOptionsProvider;
    Options            moduleOptions;
    // Fill the Options object with the module options
    moduleMdpOptionsProvider.initMdpOptions(&moduleOptions);

    // Add rules to transform mdp inputs to module data
    KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
    moduleMdpOptionsProvider.initMdpTransform(transform.rules());

    // Execute the transform on the moduleMdpValues
    auto transformedMdpValues = transform.transform(moduleMdpValues, nullptr);
    // Fill moduleOptions with values from the tree. Values for mdp
    // options are then stored in variables set up in
    // moduleMdpOptionsProvider by its initMdpOptions() method.
    assignOptionsFromKeyValueTree(&moduleOptions, transformedMdpValues.object(), nullptr);

    // Return the mdp options from the KVT now stored in the provider
    return moduleMdpOptionsProvider;
}

} // namespace test

} // namespace gmx

#endif
