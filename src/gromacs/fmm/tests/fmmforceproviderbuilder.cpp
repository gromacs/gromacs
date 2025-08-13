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
 * \brief Tests for FmmForceProviderBuilder class.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#include "gmxpre.h"

#include "gromacs/fmm/fmmforceproviderbuilder.h"

#include <gtest/gtest.h>

#include "gromacs/fmm/fmm_mdmodule.h"
#include "gromacs/fmm/fmm_mdpoptions.h"
#include "gromacs/fmm/fmmforceprovider.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/logger.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

class FmmForceProviderBuilderTest : public ::testing::Test
{
public:
    void notifyFmmModuleSimulationSetup()
    {
        fmmModule_->subscribeToSimulationSetupNotifications(&notifiers_);
        notifiers_.simulationSetupNotifier_.notify(pbcType_);
        notifiers_.simulationSetupNotifier_.notify(mtop_);
        notifiers_.simulationSetupNotifier_.notify(logger_);
    }

    void initializeForceProviders() { fmmModule_->initForceProviders(&forceProviders_); }

protected:
    ForceProviders             forceProviders_;
    std::unique_ptr<IMDModule> fmmModule_;
    MDModulesNotifiers         notifiers_;

    // default settings
    MDLogger   logger_;
    PbcType    pbcType_ = PbcType::Xyz;
    gmx_mtop_t mtop_;
};

TEST_F(FmmForceProviderBuilderTest, ThrowsIfBuilderNotReady)
{
    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);
    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());

    // Throws (internal error) because FmmForceProviderBuilder is not ready
    EXPECT_THROW(initializeForceProviders(), InternalError);
}

TEST_F(FmmForceProviderBuilderTest, BuildsWithMinimalSetup)
{
    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);
    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModuleSimulationSetup();

    if (GMX_USE_EXT_FMM)
    {
        EXPECT_NO_THROW(initializeForceProviders());
        EXPECT_TRUE(forceProviders_.hasForceProvider());
    }
    else
    {
        // Throws (internal error) from FmmForceProvider stub due to FMM selected without -DGMX_USE_EXT_FMM=ON
        EXPECT_THROW(initializeForceProviders(), InternalError);
    }
}

} // namespace test

} // namespace gmx
