/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for functionality of the NNPot MDModule
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpot.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringcompare.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

class NNPotTest : public ::testing::Test
{
public:
    void addMdpOptionNNPotActive()
    {
        mdpValueBuilder_.rootObject().addValue("nnpot-active", std::string("true"));
    }

    //! build an mdp options tree that sets the options for the NNPot module
    void makeNNPotModuleWithSetOptions()
    {
        KeyValueTreeObject mdpOptionsTree = mdpValueBuilder_.build();

        NNPotModule_ = NNPotModuleInfo::create();

        // set up options
        Options NNPotModuleOptions;
        NNPotModule_->mdpOptionProvider()->initMdpOptions(&NNPotModuleOptions);

        // Add rules to transform mdp inputs to NNPotModule data
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
        NNPotModule_->mdpOptionProvider()->initMdpTransform(transform.rules());

        // Execute the transform on the mdpValues
        auto transformedMdpValues = transform.transform(mdpOptionsTree, nullptr);
        assignOptionsFromKeyValueTree(&NNPotModuleOptions, transformedMdpValues.object(), nullptr);
    }

    void initializeForceProviders() { NNPotModule_->initForceProviders(&NNPotForces_); }

protected:
    KeyValueTreeBuilder        mdpValueBuilder_;
    ForceProviders             NNPotForces_;
    std::unique_ptr<IMDModule> NNPotModule_;
};

TEST_F(NNPotTest, ForceProviderLackingInputThrows)
{
    // Prepare MDP inputs
    addMdpOptionNNPotActive();

    // Initialize with default options
    makeNNPotModuleWithSetOptions();

    // Build the force provider
    EXPECT_ANY_THROW(initializeForceProviders());
}

} // namespace test

} // namespace gmx
