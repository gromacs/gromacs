/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for functionality of the Colvars module.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

class ColvarsTest : public ::testing::Test
{
public:
    void addMdpOptionColvarsActive()
    {
        mdpValueBuilder_.rootObject().addValue("colvars-active", std::string("yes"));
    }

    //! build an mdp options tree that sets the options for the Colvars module
    void makeColvarsModuleWithSetOptions()
    {
        KeyValueTreeObject mdpOptionsTree = mdpValueBuilder_.build();

        ColvarsModule_ = ColvarsModuleInfo::create();

        // set up options
        Options ColvarsModuleOptions;
        ColvarsModule_->mdpOptionProvider()->initMdpOptions(&ColvarsModuleOptions);

        // Add rules to transform mdp inputs to ColvarsModule data
        KeyValueTreeTransformer transform;
        transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
        ColvarsModule_->mdpOptionProvider()->initMdpTransform(transform.rules());

        // Execute the transform on the mdpValues
        auto transformedMdpValues = transform.transform(mdpOptionsTree, nullptr);
        assignOptionsFromKeyValueTree(&ColvarsModuleOptions, transformedMdpValues.object(), nullptr);
    }

    void intializeForceProviders() { ColvarsModule_->initForceProviders(&ColvarsForces_); }

protected:
    KeyValueTreeBuilder        mdpValueBuilder_;
    ForceProviders             ColvarsForces_;
    std::unique_ptr<IMDModule> ColvarsModule_;
};

TEST_F(ColvarsTest, ForceProviderLackingInputThrows)
{
    // Prepare MDP inputs
    addMdpOptionColvarsActive();

    // Initialize with default options
    makeColvarsModuleWithSetOptions();

    // Build the force provider with missing input data
    EXPECT_ANY_THROW(intializeForceProviders());
}

} // namespace

} // namespace gmx
